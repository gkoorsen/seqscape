"""Fetch norovirus complete VP1 sequences from NCBI for the SeqScape Case 1 benchmark.

Pulls Norovirus nuccore records of length 1500-2500 nt (the VP1/ORF2 range),
extracts the VP1 CDS from GenBank feature tables, filters for length and
ambiguity, deduplicates by sequence hash, and writes a FASTA plus a manifest.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import time
from pathlib import Path
from typing import Iterator

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

ORGANISM_QUERY = "Norovirus[Organism] AND 1500:2500[SLEN]"
VP1_PRODUCT_KEYWORDS = ("vp1", "capsid", "orf2", "major capsid protein")
MIN_VP1_LEN = 1500
MAX_VP1_LEN = 1800
MAX_AMBIG_FRAC = 0.01
BATCH_SIZE = 200


def esearch_ids(query: str) -> list[str]:
    handle = Entrez.esearch(db="nuccore", term=query, retmax=0, usehistory="y")
    initial = Entrez.read(handle)
    handle.close()
    total = int(initial["Count"])
    webenv = initial["WebEnv"]
    query_key = initial["QueryKey"]
    ids: list[str] = []
    for start in range(0, total, 10000):
        h = Entrez.esearch(
            db="nuccore",
            term=query,
            retstart=start,
            retmax=10000,
            webenv=webenv,
            query_key=query_key,
        )
        batch = Entrez.read(h)
        h.close()
        ids.extend(batch["IdList"])
        time.sleep(0.34)
    return ids


def efetch_records(ids: list[str]) -> Iterator[SeqRecord]:
    for start in range(0, len(ids), BATCH_SIZE):
        chunk = ids[start : start + BATCH_SIZE]
        attempt = 0
        while True:
            try:
                h = Entrez.efetch(
                    db="nuccore",
                    id=",".join(chunk),
                    rettype="gb",
                    retmode="text",
                )
                for rec in SeqIO.parse(h, "genbank"):
                    yield rec
                h.close()
                break
            except Exception as exc:
                attempt += 1
                if attempt >= 3:
                    raise
                time.sleep(2 * attempt)
        time.sleep(0.34)


def extract_vp1(record: SeqRecord) -> SeqRecord | None:
    for feat in record.features:
        if feat.type != "CDS":
            continue
        product = " ".join(feat.qualifiers.get("product", [])).lower()
        gene = " ".join(feat.qualifiers.get("gene", [])).lower()
        note = " ".join(feat.qualifiers.get("note", [])).lower()
        hay = f"{product} {gene} {note}"
        if not any(kw in hay for kw in VP1_PRODUCT_KEYWORDS):
            continue
        sub = feat.extract(record.seq)
        if not MIN_VP1_LEN <= len(sub) <= MAX_VP1_LEN:
            continue
        return SeqRecord(sub, id=record.id, description=product.strip() or "VP1")
    if MIN_VP1_LEN <= len(record.seq) <= MAX_VP1_LEN:
        return SeqRecord(record.seq, id=record.id, description="VP1 (whole record)")
    return None


def ambiguity_fraction(seq: str) -> float:
    seq = seq.upper()
    acgt = sum(seq.count(b) for b in "ACGT")
    return 1.0 - acgt / len(seq) if seq else 1.0


def main() -> int:
    ap = argparse.ArgumentParser(description="Fetch norovirus complete VP1 dataset.")
    ap.add_argument("--email", required=True, help="Contact email for NCBI Entrez.")
    ap.add_argument("--api-key", default=None, help="Optional NCBI API key.")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--max-records", type=int, default=None, help="Smoke-test cap.")
    args = ap.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    args.outdir.mkdir(parents=True, exist_ok=True)
    fasta_path = args.outdir / "norovirus_vp1.fasta"
    manifest_path = args.outdir / "norovirus_vp1_manifest.tsv"

    print(f"esearch: {ORGANISM_QUERY}", file=sys.stderr)
    ids = esearch_ids(ORGANISM_QUERY)
    print(f"  -> {len(ids)} UIDs", file=sys.stderr)
    if args.max_records:
        ids = ids[: args.max_records]

    seen_hashes: set[str] = set()
    kept = 0
    skipped_no_vp1 = 0
    skipped_ambig = 0
    skipped_dup = 0

    with fasta_path.open("w") as fa, manifest_path.open("w") as mf:
        mf.write("accession\tlength\tdescription\n")
        for rec in efetch_records(ids):
            vp1 = extract_vp1(rec)
            if vp1 is None:
                skipped_no_vp1 += 1
                continue
            seq_str = str(vp1.seq).upper()
            if ambiguity_fraction(seq_str) > MAX_AMBIG_FRAC:
                skipped_ambig += 1
                continue
            h = hashlib.sha1(seq_str.encode()).hexdigest()
            if h in seen_hashes:
                skipped_dup += 1
                continue
            seen_hashes.add(h)
            fa.write(f">{vp1.id}\n{seq_str}\n")
            mf.write(f"{vp1.id}\t{len(seq_str)}\t{vp1.description}\n")
            kept += 1
            if kept % 500 == 0:
                print(f"  kept {kept}", file=sys.stderr)

    print(
        f"done. kept={kept} no_vp1={skipped_no_vp1} "
        f"ambig={skipped_ambig} dup={skipped_dup}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
