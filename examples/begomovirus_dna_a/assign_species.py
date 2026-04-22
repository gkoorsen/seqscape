"""Assign Brown 2015 begomovirus species labels via BLASTn best-hit.

ICTV uses 91% pairwise nucleotide identity as the species demarcation threshold
(Brown et al., 2015). Strain (sub-species) threshold is 94%.

Usage:
  python assign_species.py \\
      --query runs/panel/representative_panel.fasta \\
      --ref-fasta data/brown_references_labeled.fasta \\
      --outdir runs/species_assignments
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


MIN_PIDENT = 75.0   # Lowered to catch below-species-threshold hits for "novel" reporting
MIN_QCOVS = 40.0


def run(cmd: list[str]) -> None:
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    if result.stderr:
        print(result.stderr.strip(), file=sys.stderr)


def makeblastdb(fasta: Path, dbdir: Path) -> Path:
    dbdir.mkdir(parents=True, exist_ok=True)
    db = dbdir / "refs"
    run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(db)])
    return db


def blastn(query: Path, db: Path, outfile: Path, threads: int = 8) -> None:
    outfile.parent.mkdir(parents=True, exist_ok=True)
    run([
        "blastn", "-query", str(query), "-db", str(db),
        "-out", str(outfile),
        "-outfmt", "6 qseqid sseqid pident qcovs bitscore",
        "-max_target_seqs", "5",
        "-num_threads", str(threads),
        "-dust", "no",
    ])


def parse_blast(blast_tsv: Path) -> dict[str, tuple[str, float, float]]:
    best: dict[str, tuple[float, str, float, float]] = {}
    with blast_tsv.open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            qid, sid, pident, qcovs, bitscore = (
                row[0], row[1], float(row[2]), float(row[3]), float(row[4])
            )
            if pident < MIN_PIDENT or qcovs < MIN_QCOVS:
                continue
            species = sid.split("|")[0] if "|" in sid else sid
            if qid not in best or bitscore > best[qid][0]:
                best[qid] = (bitscore, species, pident, qcovs)
    return {qid: (sp, pi, qc) for qid, (_, sp, pi, qc) in best.items()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", type=Path, required=True)
    ap.add_argument("--ref-fasta", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    db = makeblastdb(args.ref_fasta, args.outdir / "blastdb")
    print("BLAST database built.", file=sys.stderr)

    blast_out = args.outdir / "blast_hits.tsv"
    blastn(args.query, db, blast_out, threads=args.threads)
    print(f"BLAST done. Hits: {blast_out}", file=sys.stderr)

    assignments = parse_blast(blast_out)

    from Bio import SeqIO
    all_ids = [r.id for r in SeqIO.parse(str(args.query), "fasta")]

    out_tsv = args.outdir / "species_assignments.tsv"
    n_species_level = 0  # pident >= 91
    n_genus_level = 0    # 75 <= pident < 91
    n_unresolved = 0
    with out_tsv.open("w") as fh:
        fh.write("id\tspecies\tpident\tqcovs\tcategory\n")
        for qid in all_ids:
            if qid in assignments:
                sp, pi, qc = assignments[qid]
                if pi >= 91.0:
                    category = "species"
                    n_species_level += 1
                else:
                    category = "below_species"
                    n_genus_level += 1
                fh.write(f"{qid}\t{sp}\t{pi:.2f}\t{qc:.1f}\t{category}\n")
            else:
                fh.write(f"{qid}\tunassigned\t\t\tunassigned\n")
                n_unresolved += 1

    print(f"Species-level (pident≥91): {n_species_level}", file=sys.stderr)
    print(f"Below species (75-91%): {n_genus_level}", file=sys.stderr)
    print(f"Unassigned: {n_unresolved}", file=sys.stderr)
    print(f"→ {out_tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
