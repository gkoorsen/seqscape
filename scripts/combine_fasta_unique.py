#!/usr/bin/env python3
"""Combine FASTA files while keeping one record per ID."""
from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path):
    current_id = None
    current_header = ""
    chunks: list[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    yield current_id, current_header, "".join(chunks)
                current_header = line[1:]
                current_id = current_header.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if current_id is not None:
        yield current_id, current_header, "".join(chunks)


def load_exclude_ids(paths: list[str]) -> set[str]:
    ids: set[str] = set()
    for raw_path in paths:
        if not raw_path:
            continue
        for sid, _, _ in read_fasta(Path(raw_path)):
            ids.add(sid)
    return ids


def write_wrapped(handle, seq: str, width: int = 60) -> None:
    for i in range(0, len(seq), width):
        handle.write(seq[i : i + width] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Combine FASTA files with duplicate-ID removal.")
    parser.add_argument("--fasta", action="append", required=True, help="Input FASTA; may repeat")
    parser.add_argument(
        "--exclude-fasta",
        action="append",
        default=[],
        help="IDs present here are excluded from the combined output",
    )
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--duplicates-tsv", required=True)
    args = parser.parse_args()

    exclude = load_exclude_ids(args.exclude_fasta)
    seen: set[str] = set()
    written = 0
    excluded = 0
    duplicates: list[tuple[str, str, str]] = []

    out = Path(args.out_fasta)
    dup_path = Path(args.duplicates_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    dup_path.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w") as fh:
        for source in args.fasta:
            for sid, header, seq in read_fasta(Path(source)):
                if sid in exclude:
                    excluded += 1
                    duplicates.append((sid, source, "excluded_by_anchor_reference"))
                    continue
                if sid in seen:
                    duplicates.append((sid, source, "duplicate_id_skipped"))
                    continue
                seen.add(sid)
                fh.write(f">{header}\n")
                write_wrapped(fh, seq)
                written += 1

    with dup_path.open("w") as fh:
        fh.write("id\tsource\treason\n")
        for row in duplicates:
            fh.write("\t".join(row) + "\n")

    print(f"records written: {written}")
    print(f"records excluded by anchor references: {excluded}")
    print(f"duplicate/excluded rows: {len(duplicates)}")
    print(f"outputs: {out}  {dup_path}")


if __name__ == "__main__":
    main()
