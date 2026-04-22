"""Fetch a minimal norovirus VP1 reference panel of well-known genotype prototypes.

This is a stopgap reference panel for bootstrapping the SeqScape run. It is
superseded by the full Chhabra 2019 Table S1 reference set (305 sequences)
once that xlsx is downloaded manually and parsed via `parse_chhabra_s1.py`.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

PROTOTYPES: dict[str, str] = {
    "M87661.2":   "GI.1 Norwalk",
    "L07418.1":   "GI.2 Southampton",
    "U04469.1":   "GI.3 DSV395",
    "AB042808.1": "GI.4 Chiba",
    "AJ277614.1": "GI.5 Musgrove",
    "AF093797.1": "GI.6 Hesse",
    "AJ277609.1": "GI.7 Winchester",
    "U07611.3":   "GII.1 Hawaii",
    "AY134748.2": "GII.2 Snow Mountain",
    "U02030.1":   "GII.3 Toronto",
    "X76716.1":   "GII.4 Bristol",
    "JX459908.1": "GII.4 Sydney 2012",
    "AJ277620.1": "GII.6 Seacroft",
    "AJ277608.1": "GII.7 Leeds",
    "LC037415.1": "GII.17 Kawasaki",
}

VP1_KEYWORDS = ("vp1", "capsid", "orf2", "major capsid protein")
MIN_LEN = 1400
MAX_LEN = 1900


def extract_vp1(record: SeqRecord) -> SeqRecord | None:
    for feat in record.features:
        if feat.type != "CDS":
            continue
        product = " ".join(feat.qualifiers.get("product", [])).lower()
        gene = " ".join(feat.qualifiers.get("gene", [])).lower()
        note = " ".join(feat.qualifiers.get("note", [])).lower()
        hay = f"{product} {gene} {note}"
        if not any(kw in hay for kw in VP1_KEYWORDS):
            continue
        sub = feat.extract(record.seq)
        if MIN_LEN <= len(sub) <= MAX_LEN:
            return SeqRecord(sub, id=record.id, description=product.strip() or "VP1")
    if MIN_LEN <= len(record.seq) <= MAX_LEN:
        return SeqRecord(record.seq, id=record.id, description="VP1 (whole record)")
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description="Fetch norovirus VP1 prototype panel.")
    ap.add_argument("--email", required=True)
    ap.add_argument("--api-key", default=None)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    args.outdir.mkdir(parents=True, exist_ok=True)
    fasta_path = args.outdir / "norovirus_vp1_references.fasta"
    label_path = args.outdir / "norovirus_vp1_reference_labels.tsv"

    kept: list[tuple[str, str, str]] = []
    for acc, label in PROTOTYPES.items():
        try:
            h = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
            rec = SeqIO.read(h, "genbank")
            h.close()
        except Exception as exc:
            print(f"  [err]  {acc} ({label}): {exc}", file=sys.stderr)
            continue
        vp1 = extract_vp1(rec)
        if vp1 is None:
            print(f"  [skip] {acc} ({label}): no VP1 CDS", file=sys.stderr)
            continue
        kept.append((vp1.id, str(vp1.seq).upper(), label))
        print(f"  [ok]   {vp1.id} ({label}) {len(vp1.seq)} nt", file=sys.stderr)
        time.sleep(0.34)

    with fasta_path.open("w") as fa, label_path.open("w") as lb:
        lb.write("id\tlabel\n")
        for rid, seq, lbl in kept:
            fa.write(f">{rid} {lbl}\n{seq}\n")
            lb.write(f"{rid}\t{lbl}\n")

    print(f"wrote {len(kept)} prototypes to {fasta_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
