"""Fetch Brown 2015 exemplar sequences missing from the main dataset.

Outputs data/brown_references.fasta (all 252 exemplars).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO

EMAIL = "g.koorsen@googlemail.com"
BATCH = 50
SLEEP = 0.4

HERE = Path(__file__).parent
EXEMPLAR_FILE = HERE / "brown_2015_group_exemplar_accessions.txt"
DATASET_FASTA = HERE / "data" / "begomovirus_dna_a.fasta"
PARTIAL_FASTA = HERE / "data" / "brown_references_partial.fasta"
OUT_FASTA = HERE / "data" / "brown_references.fasta"


def main() -> int:
    Entrez.email = EMAIL

    with open(EXEMPLAR_FILE) as f:
        brown_accs = [line.strip() for line in f if line.strip()]

    # Which are already in our dataset?
    in_dataset = {r.id for r in SeqIO.parse(str(PARTIAL_FASTA), "fasta")}
    to_fetch = [a for a in brown_accs if a not in in_dataset]
    print(f"Brown exemplars: {len(brown_accs)}")
    print(f"In dataset: {len(in_dataset)}")
    print(f"To fetch: {len(to_fetch)}")

    fetched: list = []
    for i in range(0, len(to_fetch), BATCH):
        batch = to_fetch[i:i + BATCH]
        handle = Entrez.efetch(
            db="nucleotide", id=",".join(batch), rettype="fasta", retmode="text"
        )
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        fetched.extend(records)
        print(f"  {min(i + BATCH, len(to_fetch))}/{len(to_fetch)} fetched")
        time.sleep(SLEEP)

    # Write combined output: dataset exemplars + fetched
    dataset_recs = list(SeqIO.parse(str(PARTIAL_FASTA), "fasta"))
    combined = dataset_recs + fetched
    SeqIO.write(combined, str(OUT_FASTA), "fasta")
    print(f"Written {len(combined)} sequences → {OUT_FASTA}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
