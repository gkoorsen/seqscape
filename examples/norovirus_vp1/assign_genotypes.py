"""Assign Chhabra 2019 genotype labels to VP1 sequences via BLASTn best-hit.

Steps:
  1. makeblastdb on CDC reference FASTA (produced by fetch_cdc_references.py)
  2. blastn query FASTA against database
  3. Best hit per query: highest bitscore with pident >= MIN_PIDENT and qcovs >= MIN_QCOVS
  4. Write TSV: id, genotype, pident, qcovs, sseqid

Usage:
  python assign_genotypes.py \\
      --query   runs/panel/representative_panel_aa_clean.fasta  \\  # or nt fasta
      --ref-fasta data/cdc_references.fasta \\
      --outdir  runs/genotype_assignments \\
      [--query-type nt|aa]
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


MIN_PIDENT = 60.0   # % identity floor (nt); lowered to catch distant genotypes
MIN_QCOVS  = 40.0   # query coverage floor


def run(cmd: list[str]) -> None:
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    if result.stderr:
        print(result.stderr.strip(), file=sys.stderr)


def makeblastdb(fasta: Path, dbdir: Path, dbtype: str = "nucl") -> Path:
    dbdir.mkdir(parents=True, exist_ok=True)
    db = dbdir / "refs"
    run(["makeblastdb", "-in", str(fasta), "-dbtype", dbtype, "-out", str(db)])
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


def tblastn(query: Path, db: Path, outfile: Path, threads: int = 8) -> None:
    """For aa query against nt db."""
    outfile.parent.mkdir(parents=True, exist_ok=True)
    run([
        "tblastn", "-query", str(query), "-db", str(db),
        "-out", str(outfile),
        "-outfmt", "6 qseqid sseqid pident qcovs bitscore",
        "-max_target_seqs", "5",
        "-num_threads", str(threads),
    ])


def parse_blast(blast_tsv: Path) -> dict[str, tuple[str, float, float]]:
    """Return {query_id: (genotype, pident, qcovs)} for best hit per query."""
    best: dict[str, tuple[float, str, float, float]] = {}  # id -> (bitscore, genotype, pident, qcovs)
    with blast_tsv.open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            qid, sid, pident, qcovs, bitscore = row[0], row[1], float(row[2]), float(row[3]), float(row[4])
            if pident < MIN_PIDENT or qcovs < MIN_QCOVS:
                continue
            # sseqid format: GII.4_Sydney|JX459908.1  (from makeblastdb parse_seqids, sid may be lcl|...)
            # strip any lcl| prefix
            sid_clean = sid.split("|")[-2] if "|" in sid else sid
            # The FASTA header was >GII.4_Sydney|JX459908.1 so after makeblastdb the id is the full header
            # Try to extract genotype from sseqid - the genotype is before the first |
            if "|" in sid:
                genotype = sid.split("|")[0]
            else:
                genotype = sid
            if qid not in best or bitscore > best[qid][0]:
                best[qid] = (bitscore, genotype, pident, qcovs)
    return {qid: (gt, pi, qc) for qid, (_, gt, pi, qc) in best.items()}


def main() -> int:
    ap = argparse.ArgumentParser(description="Assign genotypes via BLASTn best-hit.")
    ap.add_argument("--query", type=Path, required=True, help="Query nucleotide FASTA")
    ap.add_argument("--ref-fasta", type=Path, required=True, help="Labeled reference FASTA from fetch_cdc_references.py")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--query-type", choices=["nt", "aa"], default="nt")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Build BLAST db (nucl always — reference is nt)
    db = makeblastdb(args.ref_fasta, args.outdir / "blastdb", dbtype="nucl")
    print("BLAST database built.", file=sys.stderr)

    blast_out = args.outdir / "blast_hits.tsv"
    if args.query_type == "aa":
        tblastn(args.query, db, blast_out, threads=args.threads)
    else:
        blastn(args.query, db, blast_out, threads=args.threads)
    print(f"BLAST done. Hits: {blast_out}", file=sys.stderr)

    assignments = parse_blast(blast_out)

    # Count all query sequences
    from Bio import SeqIO
    all_ids = [r.id for r in SeqIO.parse(str(args.query), "fasta")]

    out_tsv = args.outdir / "genotype_assignments.tsv"
    n_typed = 0
    n_untyped = 0
    with out_tsv.open("w") as fh:
        fh.write("id\tgenotype\tpident\tqcovs\n")
        for qid in all_ids:
            if qid in assignments:
                gt, pi, qc = assignments[qid]
                fh.write(f"{qid}\t{gt}\t{pi:.1f}\t{qc:.1f}\n")
                n_typed += 1
            else:
                fh.write(f"{qid}\tuntyped\t\t\n")
                n_untyped += 1

    print(f"Typed={n_typed} untyped={n_untyped} → {out_tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
