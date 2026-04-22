"""Fetch Chhabra 302 sequences in GenBank format and extract VP1 CDS.

For sequences already in our 4700-seq dataset, copies those directly.
For the rest (complete genomes), fetches from NCBI and extracts VP1.
"""
from __future__ import annotations
import csv, sys, time
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = "g.koorsen@googlemail.com"

VP1_KEYWORDS = {"vp1", "capsid protein vp1", "major capsid protein", "orf2",
                "vp60", "coat protein", "structural protein"}

def is_vp1(feat) -> bool:
    for qual in ("product", "gene", "note"):
        v = " ".join(feat.qualifiers.get(qual, [])).lower()
        if any(kw in v for kw in VP1_KEYWORDS):
            return True
    return False

def extract_vp1_nt(gb_rec) -> str | None:
    for feat in gb_rec.features:
        if feat.type == "CDS" and is_vp1(feat):
            seq = feat.extract(gb_rec.seq)
            if 1400 <= len(seq) <= 2000:
                return str(seq)
    # fallback: full sequence if already VP1-sized
    if 1400 <= len(gb_rec.seq) <= 2000:
        return str(gb_rec.seq)
    return None

def main() -> None:
    # Load Chhabra accession -> genotype mapping
    acc_to_gt: dict[str, str] = {}
    with open("examples/norovirus_vp1/data/chhabra_305_accessions.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            acc_to_gt[row["accession"].split(".")[0]] = row["genotype"]

    # Load our existing 4700-seq dataset
    existing: dict[str, SeqRecord] = {}
    for rec in SeqIO.parse("examples/norovirus_vp1/data/norovirus_vp1.fasta", "fasta"):
        existing[rec.id.split(".")[0]] = rec

    need_fetch = [acc for acc in acc_to_gt if acc not in existing]
    print(f"In existing dataset: {len(acc_to_gt) - len(need_fetch)}", file=sys.stderr)
    print(f"Need to fetch: {len(need_fetch)}", file=sys.stderr)

    out = Path("examples/norovirus_vp1/data/chhabra_302_vp1.fasta")
    n_ok = n_fail = 0
    with out.open("w") as fa:
        # Write existing ones first
        for acc_base, gt in acc_to_gt.items():
            if acc_base in existing:
                rec = existing[acc_base]
                fa.write(f">{gt}|{rec.id}\n{rec.seq}\n")
                n_ok += 1

        # Fetch remaining in batches
        for i in range(0, len(need_fetch), 50):
            chunk = need_fetch[i:i+50]
            try:
                handle = Entrez.efetch(db="nucleotide", id=",".join(chunk),
                                       rettype="gb", retmode="text")
                recs = list(SeqIO.parse(handle, "genbank"))
                handle.close()
            except Exception as e:
                print(f"  batch {i} error: {e}", file=sys.stderr)
                n_fail += len(chunk)
                time.sleep(2)
                continue
            fetched = {r.name: r for r in recs}
            for acc in chunk:
                gb = fetched.get(acc)
                if gb is None:
                    n_fail += 1
                    continue
                vp1 = extract_vp1_nt(gb)
                if vp1 is None:
                    print(f"  no VP1 found: {acc} ({len(gb.seq)} nt)", file=sys.stderr)
                    n_fail += 1
                    continue
                gt = acc_to_gt[acc]
                fa.write(f">{gt}|{acc}\n{vp1}\n")
                n_ok += 1
            time.sleep(0.5)
            print(f"  {min(i+50, len(need_fetch))}/{len(need_fetch)} fetched", file=sys.stderr)

    print(f"Written={n_ok} failed={n_fail} → {out}", file=sys.stderr)

if __name__ == "__main__":
    main()
