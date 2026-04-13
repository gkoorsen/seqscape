from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

from Bio.SeqRecord import SeqRecord

from .io_utils import load_fasta_records, write_fasta, write_tsv


def safe_slug(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", text.strip())
    return slug.strip("_") or "cluster"


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    return rows


def load_manifest_map(path: Path | None) -> dict[str, dict]:
    if path is None:
        return {}
    rows = read_rows(path)
    if "id" not in rows[0]:
        raise RuntimeError(f"Manifest missing 'id' column: {path}")
    return {row["id"]: row for row in rows}


def load_records(sequence_fasta: Path, reference_fasta: Path | None) -> dict[str, SeqRecord]:
    records = load_fasta_records(sequence_fasta)
    if reference_fasta is not None:
        ref_records = load_fasta_records(reference_fasta)
        overlap = set(records) & set(ref_records)
        if overlap:
            raise RuntimeError(f"Overlapping ids between sequence and reference FASTA: {sorted(list(overlap))[:5]}")
        records.update(ref_records)
    return records


def select_cluster_rows(rows: list[dict], id_column: str, cluster_column: str, filter_column: str, filter_value: str) -> list[dict]:
    required = {id_column, cluster_column}
    missing = required - set(rows[0].keys())
    if missing:
        raise RuntimeError(f"Missing columns in cluster TSV: {sorted(missing)}")
    if filter_column:
        if filter_column not in rows[0]:
            raise RuntimeError(f"Missing filter column in cluster TSV: {filter_column}")
        rows = [row for row in rows if row.get(filter_column, "") == filter_value]
    if not rows:
        detail = f" after filtering {filter_column}={filter_value}" if filter_column else ""
        raise RuntimeError(f"No cluster rows found{detail}")
    return rows


def run(args) -> None:
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    sequence_fasta = Path(args.sequence_fasta).resolve()
    reference_fasta = Path(args.reference_fasta).resolve() if args.reference_fasta else None
    cluster_tsv = Path(args.cluster_tsv).resolve()
    manifest_tsv = Path(args.manifest_tsv).resolve() if args.manifest_tsv else None

    records = load_records(sequence_fasta, reference_fasta)
    manifest = load_manifest_map(manifest_tsv)
    cluster_rows = select_cluster_rows(
        read_rows(cluster_tsv),
        args.id_column,
        args.cluster_column,
        args.filter_column,
        args.filter_value,
    )

    wanted_clusters = {
        c.strip() for c in args.cluster_list.split(",") if c.strip()
    } if args.cluster_list else None

    cluster_to_ids: dict[str, list[str]] = defaultdict(list)
    for row in cluster_rows:
        cluster_id = row[args.cluster_column]
        if wanted_clusters is not None and cluster_id not in wanted_clusters:
            continue
        cluster_to_ids[cluster_id].append(row[args.id_column])

    if not cluster_to_ids:
        raise RuntimeError("No clusters selected for export")

    summary_rows: list[dict] = []
    export_prefix = safe_slug(args.filename_prefix) if args.filename_prefix else ""
    for idx, cluster_id in enumerate(sorted(cluster_to_ids), start=1):
        seq_ids = sorted(cluster_to_ids[cluster_id])
        out_records: list[SeqRecord] = []
        genome_count = 0
        reference_count = 0
        for seq_id in seq_ids:
            rec = records.get(seq_id)
            if rec is None:
                raise RuntimeError(f"Sequence id missing from provided FASTA inputs: {seq_id}")
            meta = manifest.get(seq_id, {})
            label = meta.get("label", seq_id)
            item_class = meta.get("item_class", "")
            source_cluster = meta.get("source_cluster", "")
            desc_bits = [f"id={seq_id}", f"cluster={cluster_id}"]
            if item_class:
                desc_bits.append(f"item_class={item_class}")
            if source_cluster:
                desc_bits.append(f"source_cluster={source_cluster}")
            out_records.append(SeqRecord(rec.seq, id=label, name="", description=" ".join(desc_bits)))
            if item_class == "reference":
                reference_count += 1
            else:
                genome_count += 1

        stem = f"{export_prefix + '_' if export_prefix else ''}{safe_slug(cluster_id)}"
        out_path = outdir / f"{stem}.fasta"
        write_fasta(out_records, out_path)
        summary_rows.append(
            {
                "cluster": cluster_id,
                "file_name": out_path.name,
                "members": len(out_records),
                "genomes": genome_count,
                "references": reference_count,
            }
        )

    write_tsv(outdir / "export_summary.tsv", summary_rows, ["cluster", "file_name", "members", "genomes", "references"])
    with (outdir / "summary.txt").open("w") as fh:
        fh.write(f"sequence_fasta\t{sequence_fasta}\n")
        fh.write(f"reference_fasta\t{reference_fasta if reference_fasta else ''}\n")
        fh.write(f"cluster_tsv\t{cluster_tsv}\n")
        fh.write(f"manifest_tsv\t{manifest_tsv if manifest_tsv else ''}\n")
        fh.write(f"id_column\t{args.id_column}\n")
        fh.write(f"cluster_column\t{args.cluster_column}\n")
        fh.write(f"filter_column\t{args.filter_column}\n")
        fh.write(f"filter_value\t{args.filter_value}\n")
        fh.write(f"cluster_list\t{args.cluster_list}\n")
        fh.write(f"exports\t{len(summary_rows)}\n")
