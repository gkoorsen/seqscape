from __future__ import annotations

import csv
from pathlib import Path
from typing import Sequence

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_fasta_records(path: Path) -> dict[str, SeqRecord]:
    out: dict[str, SeqRecord] = {}
    with path.open() as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            seq = Seq(str(rec.seq).upper().replace("-", ""))
            out[str(rec.id)] = SeqRecord(seq, id=str(rec.id), name="", description="")
    if not out:
        raise RuntimeError(f"No sequences loaded from {path}")
    return out


def load_references(path: Path) -> list[dict]:
    refs = []
    seen = set()
    with path.open() as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            rid = str(rec.id)
            if rid in seen:
                continue
            seen.add(rid)
            label = rec.description.strip() if rec.description else rid
            refs.append(
                {
                    "ref_id": rid,
                    "ref_label": label,
                    "sequence": str(rec.seq).upper(),
                    "len": len(rec.seq),
                }
            )
    return refs


def read_assignments(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    needed = {"id", "label", "item_class", "cluster"}
    missing = needed - set(rows[0].keys())
    if missing:
        raise RuntimeError(f"Missing columns in {path}: {sorted(missing)}")
    return rows


def read_coords(path: Path) -> dict[str, tuple[float, float]]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    needed = {"ID", "AF_UMAP1", "AF_UMAP2"}
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    missing = needed - set(rows[0].keys())
    if missing:
        raise RuntimeError(f"Missing columns in {path}: {sorted(missing)}")
    out: dict[str, tuple[float, float]] = {}
    for row in rows:
        out[row["ID"]] = (float(row["AF_UMAP1"]), float(row["AF_UMAP2"]))
    return out


def write_fasta(records: list[SeqRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        SeqIO.write(records, fh, "fasta")


def write_tsv(path: Path, rows: list[dict], headers: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(headers), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in headers})
