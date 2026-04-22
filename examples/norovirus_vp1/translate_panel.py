"""Translate the nucleotide VP1 panel FASTA to amino acid sequences.

Norovirus VP1 reference sequences are not always deposited in-frame, and partial
sequences may be missing the initial ATG. To be robust we try all three forward
reading frames and keep whichever gives the longest ORF without a premature
internal stop codon. Terminal stop codons are stripped.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def translate_best_frame(seq: Seq) -> tuple[int, str, int] | None:
    """Return (frame, aa_string, stops_in_middle) for the best forward frame."""
    candidates: list[tuple[int, int, str, int]] = []
    for frame in range(3):
        s = seq[frame:]
        trimmed = s[: len(s) - (len(s) % 3)]
        if len(trimmed) == 0:
            continue
        aa = str(trimmed.translate())
        if aa.endswith("*"):
            aa_no_term = aa[:-1]
        else:
            aa_no_term = aa
        stops_middle = aa_no_term.count("*")
        candidates.append((stops_middle, -len(aa_no_term), aa_no_term, frame))
    if not candidates:
        return None
    candidates.sort()
    stops_middle, neg_len, aa_no_term, frame = candidates[0]
    return frame, aa_no_term, stops_middle


def main() -> int:
    ap = argparse.ArgumentParser(description="Translate nt VP1 panel to aa.")
    ap.add_argument("--nt-fasta", type=Path, required=True)
    ap.add_argument("--aa-fasta", type=Path, required=True)
    ap.add_argument("--report-tsv", type=Path, required=True)
    args = ap.parse_args()

    n_ok = 0
    n_with_internal_stop = 0
    n_skipped = 0
    with args.aa_fasta.open("w") as fa, args.report_tsv.open("w") as rep:
        rep.write("id\tnt_len\tframe\taa_len\tinternal_stops\n")
        for rec in SeqIO.parse(str(args.nt_fasta), "fasta"):
            result = translate_best_frame(rec.seq)
            if result is None:
                n_skipped += 1
                continue
            frame, aa, stops = result
            if stops > 0:
                n_with_internal_stop += 1
            n_ok += 1
            fa.write(f">{rec.id}\n{aa}\n")
            rep.write(f"{rec.id}\t{len(rec.seq)}\t{frame}\t{len(aa)}\t{stops}\n")

    print(
        f"translated={n_ok} skipped={n_skipped} "
        f"internal_stops>0={n_with_internal_stop}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
