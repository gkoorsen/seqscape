#!/usr/bin/env python3
"""Select focal IDs that protect regional diversity during genome collapse.

The recombination panel is built from whole genomes, but candidate resistance
breaking may be carried by a short origin/IR tract. This helper chooses one
well-supported representative for each distinct regional sequence cluster, so
`seqscape collapse-panel` can collapse whole-genome near-duplicates without
discarding the regional haplotypes the screen is looking for.
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    current_id: str | None = None
    chunks: list[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, clean_seq("".join(chunks))))
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current_id is not None:
        records.append((current_id, clean_seq("".join(chunks))))
    return records


def clean_seq(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq).upper().replace("U", "T")


def load_support(path: str) -> dict[str, int]:
    if not path:
        return {}
    with Path(path).open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return {}
    id_col = "unique_id" if "unique_id" in rows[0] else "id"
    support_col = "support" if "support" in rows[0] else ""
    support: dict[str, int] = {}
    for row in rows:
        sid = row.get(id_col, "").strip()
        if not sid:
            continue
        raw = row.get(support_col, "1").strip() if support_col else "1"
        try:
            support[sid] = int(raw) if raw else 1
        except ValueError:
            support[sid] = 1
    return support


def load_ids(path: str) -> list[str]:
    if not path:
        return []
    return [
        line.strip()
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def identity_pct(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    matches = sum(1 for i in range(n) if a[i] == b[i])
    return 100.0 * matches / max(len(a), len(b))


def region_name(path: Path) -> str:
    name = path.name
    if name.startswith("region_"):
        name = name[len("region_") :]
    if name.endswith(".fasta"):
        name = name[: -len(".fasta")]
    return name


def select_for_region(
    records: list[tuple[str, str]],
    support: dict[str, int],
    threshold: float,
    min_support: int,
) -> tuple[list[str], dict[str, list[str]]]:
    filtered = [
        (sid, seq)
        for sid, seq in records
        if seq and support.get(sid, 1) >= min_support
    ]
    filtered.sort(key=lambda item: (-support.get(item[0], 1), -len(item[1]), item[0]))

    kept: list[tuple[str, str]] = []
    absorbed: dict[str, list[str]] = defaultdict(list)
    for sid, seq in filtered:
        representative = ""
        for kept_id, kept_seq in kept:
            if identity_pct(seq, kept_seq) >= threshold:
                representative = kept_id
                break
        if representative:
            absorbed[representative].append(sid)
        else:
            kept.append((sid, seq))
    return [sid for sid, _ in kept], absorbed


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Choose region-diversity IDs to protect in seqscape collapse-panel."
    )
    parser.add_argument(
        "--region-fasta",
        action="append",
        required=True,
        help="Regional FASTA. May be supplied more than once.",
    )
    parser.add_argument("--support-tsv", default="", help="Optional unique_id/support TSV")
    parser.add_argument(
        "--identity-threshold",
        type=float,
        default=98.1,
        help="Regional identity at or above which records are represented by one ID",
    )
    parser.add_argument("--min-support", type=int, default=1)
    parser.add_argument("--always-keep", default="", help="Optional IDs that must be included")
    parser.add_argument("--out-ids", required=True, help="One protected ID per line")
    parser.add_argument("--out-groups-tsv", required=True)
    args = parser.parse_args()

    support = load_support(args.support_tsv)
    protected: list[str] = []
    protected_seen: set[str] = set()
    rows: list[dict[str, str | int]] = []

    for fasta_arg in args.region_fasta:
        fasta = Path(fasta_arg)
        region = region_name(fasta)
        records = read_fasta(fasta)
        kept, absorbed = select_for_region(
            records,
            support=support,
            threshold=args.identity_threshold,
            min_support=args.min_support,
        )
        for sid in kept:
            if sid not in protected_seen:
                protected.append(sid)
                protected_seen.add(sid)
            rows.append(
                {
                    "region": region,
                    "representative_id": sid,
                    "representative_support": support.get(sid, 1),
                    "n_absorbed": len(absorbed.get(sid, [])),
                    "absorbed_ids": ";".join(sorted(absorbed.get(sid, []))),
                }
            )
        print(
            f"{region}: records={len(records)} protected={len(kept)} "
            f"threshold={args.identity_threshold}"
        )

    for sid in load_ids(args.always_keep):
        if sid not in protected_seen:
            protected.append(sid)
            protected_seen.add(sid)

    out_ids = Path(args.out_ids)
    out_groups = Path(args.out_groups_tsv)
    out_ids.parent.mkdir(parents=True, exist_ok=True)
    out_groups.parent.mkdir(parents=True, exist_ok=True)
    out_ids.write_text("\n".join(protected) + ("\n" if protected else ""))

    columns = [
        "region",
        "representative_id",
        "representative_support",
        "n_absorbed",
        "absorbed_ids",
    ]
    with out_groups.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"protected IDs: {len(protected)}")
    print(f"outputs: {out_ids}  {out_groups}")


if __name__ == "__main__":
    main()
