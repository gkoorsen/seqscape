"""Fetch ICTV species names for the 252 Brown 2015 exemplar accessions via Entrez.

Output: data/brown_reference_species.tsv (accession\tspecies) with real organism names.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO

EMAIL = "g.koorsen@googlemail.com"
BATCH = 100
SLEEP = 0.4

HERE = Path(__file__).parent
ACC_FILE = HERE / "brown_2015_group_exemplar_accessions.txt"
OUT_TSV = HERE / "data" / "brown_reference_species.tsv"


def main() -> int:
    Entrez.email = EMAIL
    with open(ACC_FILE) as f:
        accs = [line.strip() for line in f if line.strip()]
    print(f"Fetching organism for {len(accs)} accessions...", file=sys.stderr)

    acc_to_species: dict[str, str] = {}
    for i in range(0, len(accs), BATCH):
        batch = accs[i:i + BATCH]
        handle = Entrez.efetch(
            db="nucleotide", id=",".join(batch), rettype="gb", retmode="text"
        )
        for rec in SeqIO.parse(handle, "genbank"):
            acc_to_species[rec.id] = rec.annotations.get("organism", "unknown")
        handle.close()
        print(f"  {min(i + BATCH, len(accs))}/{len(accs)}", file=sys.stderr)
        time.sleep(SLEEP)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w") as fh:
        fh.write("accession\tspecies\n")
        for acc in accs:
            sp = acc_to_species.get(acc, "unknown")
            fh.write(f"{acc}\t{sp}\n")

    n_unknown = sum(1 for a in accs if acc_to_species.get(a, "unknown") == "unknown")
    print(f"Written {len(accs)} entries → {OUT_TSV} ({n_unknown} unknown)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
