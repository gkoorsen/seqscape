#!/usr/bin/env python3
"""Extract origin-centered circular-genome windows from a FASTA file.

This is intended for the ToCSV/TYLCV resistance-breaking screen, where the
published lesion is a short IR/origin-region recombinant tract. It deliberately
uses the conserved begomovirus nonanucleotide rather than GenBank feature
coordinates so the same extraction can be applied to local phased genomes and
public comparator FASTAs.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


DEFAULT_MOTIF = "TAATATTAC"


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
                    records.append((current_id, normalise_sequence("".join(chunks))))
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current_id is not None:
        records.append((current_id, normalise_sequence("".join(chunks))))
    return records


def normalise_sequence(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq).upper().replace("U", "T")


def load_ids(path: str) -> set[str] | None:
    if not path:
        return None
    ids = {
        line.strip()
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    return ids


def hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def best_motif_match(seq: str, motif: str) -> tuple[int, int, str]:
    """Return 0-based motif start, mismatch count, and observed motif string."""
    if not seq:
        return -1, len(motif), ""
    n = len(seq)
    wrap = seq + seq[: max(0, len(motif) - 1)]
    best_pos = 0
    best_mismatches = len(motif) + 1
    best_match = ""
    for pos in range(n):
        candidate = wrap[pos : pos + len(motif)]
        mismatches = hamming(candidate, motif)
        if mismatches < best_mismatches:
            best_pos = pos
            best_mismatches = mismatches
            best_match = candidate
            if mismatches == 0:
                break
    return best_pos, best_mismatches, best_match


def circular_slice(seq: str, start: int, length: int) -> str:
    if not seq or length <= 0:
        return ""
    n = len(seq)
    return "".join(seq[(start + i) % n] for i in range(length))


def write_wrapped(handle, seq: str, width: int = 60) -> None:
    for i in range(0, len(seq), width):
        handle.write(seq[i : i + width] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract fixed-length windows centered on the begomovirus origin motif."
    )
    parser.add_argument("--fasta", required=True, help="Input circular-genome FASTA")
    parser.add_argument("--ids", default="", help="Optional one-ID-per-line subset")
    parser.add_argument("--window-size", type=int, default=400, help="Total window length in nt")
    parser.add_argument("--motif", default=DEFAULT_MOTIF, help="Origin motif to centre on")
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=-1,
        help="Skip records whose best motif match has more mismatches; -1 keeps all records",
    )
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--manifest-tsv", required=True)
    args = parser.parse_args()

    motif = normalise_sequence(args.motif)
    if not motif:
        sys.exit("--motif must contain at least one nucleotide")
    if args.window_size < len(motif):
        sys.exit("--window-size must be at least as long as --motif")

    keep = load_ids(args.ids)
    records = read_fasta(Path(args.fasta))
    out_fasta = Path(args.out_fasta)
    manifest = Path(args.manifest_tsv)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    skipped_id = 0
    skipped_empty = 0
    skipped_motif = 0
    half = args.window_size // 2

    with out_fasta.open("w") as fa, manifest.open("w") as mf:
        mf.write(
            "\t".join(
                [
                    "id",
                    "source_length",
                    "window_size",
                    "window_start_1based",
                    "motif_start_1based",
                    "motif_start_in_window",
                    "motif_end_in_window",
                    "motif_mismatches",
                    "motif_match",
                ]
            )
            + "\n"
        )
        for record_id, seq in records:
            if keep is not None and record_id not in keep:
                skipped_id += 1
                continue
            if not seq:
                skipped_empty += 1
                continue
            motif_start, mismatches, match = best_motif_match(seq, motif)
            if args.max_mismatches >= 0 and mismatches > args.max_mismatches:
                skipped_motif += 1
                continue
            motif_center = motif_start + len(motif) // 2
            window_start = motif_center - half
            window = circular_slice(seq, window_start, args.window_size)
            motif_start_in_window = ((motif_start - window_start) % len(seq)) + 1
            motif_end_in_window = motif_start_in_window + len(motif) - 1
            fa.write(f">{record_id}\n")
            write_wrapped(fa, window)
            mf.write(
                "\t".join(
                    [
                        record_id,
                        str(len(seq)),
                        str(args.window_size),
                        str((window_start % len(seq)) + 1),
                        str(motif_start + 1),
                        str(motif_start_in_window),
                        str(motif_end_in_window),
                        str(mismatches),
                        match,
                    ]
                )
                + "\n"
            )
            written += 1

    print(f"records read: {len(records)}")
    print(f"records written: {written}")
    print(f"skipped by id filter: {skipped_id}")
    print(f"skipped empty: {skipped_empty}")
    print(f"skipped by motif mismatch filter: {skipped_motif}")
    print(f"outputs: {out_fasta}  {manifest}")


if __name__ == "__main__":
    main()
