from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .umap_explorer import format_resolution
from .io_utils import load_fasta_records, read_assignments, read_coords, write_fasta, write_tsv
from .validation import compression_ratio, full_pairwise_count, panel_fraction_pct, pearson_cluster_proportion_preservation


def largest_remainder_alloc(weights: dict[str, int], total: int) -> dict[str, int]:
    grand = sum(weights.values())
    if grand <= 0:
        raise RuntimeError("Cannot allocate from zero total weight")
    exact = {k: (weights[k] * total / grand) for k in weights}
    base = {k: int(math.floor(v)) for k, v in exact.items()}
    remaining = total - sum(base.values())
    order = sorted(weights, key=lambda k: (-(exact[k] - base[k]), -weights[k], k))
    for k in order[:remaining]:
        base[k] += 1
    return base


def farthest_point_order(ids: list[str], coords_map: dict[str, tuple[float, float]]) -> tuple[list[str], dict[str, float]]:
    pts = np.array([coords_map[seq_id] for seq_id in ids], dtype=float)
    centroid = pts.mean(axis=0)
    d_centroid = np.linalg.norm(pts - centroid, axis=1)
    first_idx = int(np.argmin(d_centroid))
    selected = [first_idx]
    selected_set = {first_idx}
    min_dist = np.linalg.norm(pts - pts[first_idx], axis=1)
    while len(selected) < len(ids):
        best_idx = None
        best_score = None
        for idx in range(len(ids)):
            if idx in selected_set:
                continue
            score = (float(min_dist[idx]), float(d_centroid[idx]))
            if best_score is None or score > best_score:
                best_score = score
                best_idx = idx
        assert best_idx is not None
        selected.append(best_idx)
        selected_set.add(best_idx)
        min_dist = np.minimum(min_dist, np.linalg.norm(pts - pts[best_idx], axis=1))
    order_ids = [ids[idx] for idx in selected]
    centroid_map = {seq_id: float(d_centroid[idx]) for idx, seq_id in enumerate(ids)}
    return order_ids, centroid_map


def select_genome_panel(assignments: list[dict], coords_map: dict[str, tuple[float, float]], genome_budget: int) -> tuple[list[dict], list[dict]]:
    genome_rows = [row for row in assignments if row["item_class"] == "genome"]
    cluster_to_ids: dict[str, list[str]] = defaultdict(list)
    id_to_row = {row["id"]: row for row in genome_rows}
    for row in genome_rows:
        if row["id"] not in coords_map:
            raise RuntimeError(f"Missing AF coords for genome {row['id']}")
        cluster_to_ids[row["cluster"]].append(row["id"])
    quotas = largest_remainder_alloc({k: len(v) for k, v in cluster_to_ids.items()}, genome_budget)
    selected_rows = []
    quota_rows = []
    for cluster_id in sorted(cluster_to_ids):
        ids = cluster_to_ids[cluster_id]
        order_ids, centroid_map = farthest_point_order(ids, coords_map)
        quota = min(quotas.get(cluster_id, 0), len(order_ids))
        quota_rows.append({"cluster": cluster_id, "genomes_in_cluster": len(ids), "allocated_genome_samples": quota})
        for rank, seq_id in enumerate(order_ids[:quota], start=1):
            row = id_to_row[seq_id]
            x, y = coords_map[seq_id]
            selected_rows.append(
                {
                    "id": seq_id,
                    "label": row["label"],
                    "item_class": "genome",
                    "source_cluster": cluster_id,
                    "selection_rank": rank,
                    "selection_mode": "centroid" if rank == 1 else "coverage",
                    "af_umap1": f"{x:.6f}",
                    "af_umap2": f"{y:.6f}",
                    "centroid_distance": f"{centroid_map[seq_id]:.6f}",
                }
            )
    return selected_rows, quota_rows


def select_reference_panel(assignments: list[dict], coords_map: dict[str, tuple[float, float]]) -> list[dict]:
    refs = [row for row in assignments if row["item_class"] == "reference"]
    selected = []
    for row in refs:
        if row["id"] not in coords_map:
            raise RuntimeError(f"Missing AF coords for reference {row['id']}")
        x, y = coords_map[row["id"]]
        selected.append(
            {
                "id": row["id"],
                "label": row["label"],
                "item_class": "reference",
                "source_cluster": row["cluster"],
                "selection_rank": 0,
                "selection_mode": "reference",
                "af_umap1": f"{x:.6f}",
                "af_umap2": f"{y:.6f}",
                "centroid_distance": "",
            }
        )
    return selected


def run(args) -> None:
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    chosen_slug = (
        f"k{args.chosen_kmer}_n{args.chosen_neighbors}_d{str(args.chosen_min_dist).replace('.', 'p')}"
        f"_r{format_resolution(args.chosen_leiden_resolution)}"
    )
    state_dir = Path(args.umap_explorer_dir).resolve() / "states" / chosen_slug
    assignments_tsv = state_dir / "assignments.tsv"
    coords_csv = state_dir / "coords.csv"
    if not assignments_tsv.is_file() or not coords_csv.is_file():
        legacy_slug = f"k{args.chosen_kmer}_n{args.chosen_neighbors}_d{str(args.chosen_min_dist).replace('.', 'p')}"
        legacy_state_dir = Path(args.umap_explorer_dir).resolve() / "states" / legacy_slug
        legacy_assignments = legacy_state_dir / "assignments.tsv"
        legacy_coords = legacy_state_dir / "coords.csv"
        if legacy_assignments.is_file() and legacy_coords.is_file():
            chosen_slug = legacy_slug
            state_dir = legacy_state_dir
            assignments_tsv = legacy_assignments
            coords_csv = legacy_coords
        else:
            raise RuntimeError(f"Chosen AF state not found: {state_dir}")

    genome_records = load_fasta_records(Path(args.input_fasta).resolve())
    ref_records = load_fasta_records(Path(args.reference_fasta).resolve())
    assignments = read_assignments(assignments_tsv)
    coords_map = read_coords(coords_csv)
    ref_count = sum(1 for row in assignments if row["item_class"] == "reference")
    if args.sample_size <= ref_count:
        raise RuntimeError(f"Sample size {args.sample_size} must exceed reference count {ref_count}")
    genome_budget = args.sample_size - ref_count
    selected_genomes, quota_rows = select_genome_panel(assignments, coords_map, genome_budget)
    selected_refs = select_reference_panel(assignments, coords_map)
    selected_rows = selected_genomes + selected_refs

    panel_items = []
    panel_records: list[SeqRecord] = []
    for row in selected_rows:
        rec = ref_records.get(row["id"]) if row["item_class"] == "reference" else genome_records.get(row["id"])
        if rec is None:
            raise RuntimeError(f"Selected id {row['id']} not found in input FASTA")
        seq = str(rec.seq).upper()
        item = dict(row)
        item["length_bp"] = len(seq)
        panel_items.append(item)
        panel_records.append(SeqRecord(Seq(seq), id=row["id"], name="", description=""))

    manifest_tsv = outdir / "representative_panel_manifest.tsv"
    quota_tsv = outdir / "source_cluster_quotas.tsv"
    panel_fasta = outdir / "representative_panel.fasta"
    summary_txt = outdir / "summary.txt"
    write_fasta(panel_records, panel_fasta)
    write_tsv(manifest_tsv, sorted(panel_items, key=lambda row: (row["item_class"] != "reference", row["source_cluster"], int(row["selection_rank"]))), ["id", "label", "item_class", "source_cluster", "selection_rank", "selection_mode", "length_bp", "af_umap1", "af_umap2", "centroid_distance"])
    write_tsv(quota_tsv, quota_rows, ["cluster", "genomes_in_cluster", "allocated_genome_samples"])
    full_cluster_counts = {row["cluster"]: int(row["genomes_in_cluster"]) for row in quota_rows}
    panel_cluster_counts = {row["cluster"]: int(row["allocated_genome_samples"]) for row in quota_rows}
    pearson_r, pearson_p = pearson_cluster_proportion_preservation(full_cluster_counts, panel_cluster_counts)
    validation_metrics = {
        "N": len(genome_records),
        "sample_size_total": args.sample_size,
        "sample_size_genomes": len(selected_genomes),
        "sample_size_references": len(selected_refs),
        "panel_fraction_pct": panel_fraction_pct(args.sample_size, len(genome_records)),
        "full_pairwise_count": full_pairwise_count(len(genome_records)),
        "compression_ratio_vs_panel": compression_ratio(len(genome_records), args.sample_size),
        "pearson_r_full_vs_panel_cluster_proportions": pearson_r,
        "pearson_pvalue": pearson_p,
    }
    (outdir / "validation_metrics.json").write_text(json.dumps(validation_metrics, indent=2))
    with summary_txt.open("w") as fh:
        fh.write(f"input_fasta\t{Path(args.input_fasta).resolve()}\n")
        fh.write(f"reference_fasta\t{Path(args.reference_fasta).resolve()}\n")
        fh.write(f"umap_explorer_dir\t{Path(args.umap_explorer_dir).resolve()}\n")
        fh.write(f"chosen_state\t{chosen_slug}\n")
        fh.write(f"sample_size\t{args.sample_size}\n")
        fh.write(f"references_retained\t{len(selected_refs)}\n")
        fh.write(f"genome_budget\t{genome_budget}\n")
        fh.write(f"panel_fraction_pct\t{validation_metrics['panel_fraction_pct']:.6f}\n")
        fh.write(f"full_pairwise_count\t{validation_metrics['full_pairwise_count']}\n")
        fh.write(f"compression_ratio_vs_panel\t{validation_metrics['compression_ratio_vs_panel']:.6f}\n")
        fh.write(f"pearson_r_full_vs_panel_cluster_proportions\t{validation_metrics['pearson_r_full_vs_panel_cluster_proportions']:.6f}\n")
        fh.write(f"pearson_pvalue\t{validation_metrics['pearson_pvalue']:.6e}\n")
        fh.write("selection_scheme\tall_references_kept; genomes allocated proportionally across Leiden clusters; within-cluster centroid-first then farthest-point coverage\n")
        fh.write(f"representative_panel_fasta\t{panel_fasta}\n")
        fh.write(f"representative_panel_manifest\t{manifest_tsv}\n")
        fh.write(f"source_cluster_quotas\t{quota_tsv}\n")
