from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from .clustering import fit_agglomerative, normalize_cluster_labels
from .heatmap import maybe_heatmap, write_html_heatmap
from .io_utils import write_tsv
from .ordination import pcoa, run_precomputed_umap
from .palettes import PALETTE
from .validation import choose_closest_cluster_threshold, nmi, pcoa_variance_explained


def cluster_sort_key(cid: str) -> tuple:
    head = "".join(ch for ch in cid if ch.isalpha())
    tail = "".join(ch for ch in cid if ch.isdigit())
    return (head, int(tail) if tail else 0, cid)


def read_square_matrix(path: Path) -> Tuple[List[str], np.ndarray]:
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        ids = header[1:]
        rows = []
        row_ids = []
        for row in reader:
            if not row:
                continue
            row_ids.append(row[0])
            rows.append([float(v) for v in row[1:]])
    if ids != row_ids:
        raise RuntimeError(f"Matrix ids mismatch in {path}")
    return ids, np.array(rows, dtype=float)


def read_manifest(path: Path) -> List[dict]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    return rows


def read_optional_json(path: Path) -> dict:
    if not path.is_file():
        return {}
    return json.loads(path.read_text())


def read_key_value_summary(path: Path) -> dict:
    if not path.is_file():
        return {}
    values = {}
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            key, value = line.split("\t", 1)
            values[key] = value
    return values


def parse_state_id(state_id: str) -> dict:
    match = re.fullmatch(r"k(?P<k>\d+)_n(?P<n>\d+)_d(?P<d>[0-9p]+)(?:_r(?P<r>.+))?", state_id)
    if not match:
        return {"state_id": state_id, "kmer": None, "neighbors": None, "min_dist": None, "leiden_resolution": None}
    d_text = match.group("d").replace("p", ".")
    r_text = match.group("r") or ""
    return {
        "state_id": state_id,
        "kmer": int(match.group("k")),
        "neighbors": int(match.group("n")),
        "min_dist": float(d_text),
        "leiden_resolution": float(r_text) if r_text else None,
    }


def state_sort_key(state: dict) -> tuple:
    return (
        state["kmer"] if state["kmer"] is not None else 10**9,
        state["neighbors"] if state["neighbors"] is not None else 10**9,
        state["min_dist"] if state["min_dist"] is not None else 10**9,
        state["leiden_resolution"] if state["leiden_resolution"] is not None else 10**9,
        state["state_id"],
    )


def read_af_states(af_widget_dir: Path, panel_ids: set[str], genome_ids: set[str], default_state_id: str) -> tuple[list[dict], int]:
    states_dir = af_widget_dir / "states"
    if not states_dir.is_dir():
        return [], 0
    states = []
    for state_dir in sorted([p for p in states_dir.iterdir() if p.is_dir()]):
        coords_path = state_dir / "coords.csv"
        assignments_path = state_dir / "assignments.tsv"
        if not coords_path.is_file() or not assignments_path.is_file():
            continue
        coords = {}
        with coords_path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                seq_id = row["ID"]
                if seq_id in panel_ids:
                    coords[seq_id] = [float(row["AF_UMAP1"]), float(row["AF_UMAP2"])]
        clusters = {}
        with assignments_path.open(newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                seq_id = row["id"]
                if seq_id in panel_ids:
                    clusters[seq_id] = row["cluster"]
        if not coords or not clusters:
            continue
        genome_counts = Counter(cid for seq_id, cid in clusters.items() if seq_id in genome_ids)
        order = sorted(genome_counts, key=cluster_sort_key)
        colors = {cid: PALETTE[i % len(PALETTE)] for i, cid in enumerate(order)}
        state = {
            **parse_state_id(state_dir.name),
            "coords": coords,
            "clusters": clusters,
            "cluster_order": order,
            "cluster_colors": colors,
            "cluster_count": len(order),
            "panel_cluster_counts": dict(genome_counts),
        }
        states.append(state)
    states.sort(key=state_sort_key)
    default_index = 0
    for idx, state in enumerate(states):
        if state["state_id"] == default_state_id:
            default_index = idx
            break
    return states, default_index


def convex_hull(points: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts

    def cross(o: Tuple[float, float], a: Tuple[float, float], b: Tuple[float, float]) -> float:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: List[Tuple[float, float]] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper: List[Tuple[float, float]] = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def build_cluster_view(points: List[dict], coord_key: str, cluster_key: str, order: List[str], colors: Dict[str, str]) -> dict:
    grouped: Dict[str, List[Tuple[float, float]]] = defaultdict(list)
    genome_counts = Counter()
    for row in points:
        if row["item_class"] == "genome":
            grouped[row[cluster_key]].append(tuple(row[coord_key]))
            genome_counts[row[cluster_key]] += 1
    total_genomes = sum(genome_counts.values())
    labels = []
    hulls = []
    for idx, cid in enumerate(order):
        pts = grouped.get(cid, [])
        if pts:
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            labels.append(
                {
                    "cluster": cid,
                    "x": sum(xs) / len(xs),
                    "y": sum(ys) / len(ys),
                    "genomes": genome_counts[cid],
                    "pct": (100.0 * genome_counts[cid] / total_genomes) if total_genomes else 0.0,
                    "offset_x": ((idx % 3) - 1) * 18,
                    "offset_y": ((idx % 4) - 1.5) * 12,
                }
            )
            if len(pts) >= 3:
                hull = convex_hull(pts)
                if len(hull) >= 3:
                    hulls.append({"cluster": cid, "points": hull})
    return {
        "cluster_order": order,
        "cluster_colors": colors,
        "cluster_genome_counts": {cid: genome_counts.get(cid, 0) for cid in order},
        "cluster_pct": {cid: (100.0 * genome_counts.get(cid, 0) / total_genomes) if total_genomes else 0.0 for cid in order},
        "labels": labels,
        "hulls": hulls,
    }


def build_heatmaps(ids: List[str], ident: np.ndarray, manifest_rows: List[dict], outdir: Path) -> None:
    manifest = {row["id"]: row for row in manifest_rows}
    ref_ids = [row["id"] for row in manifest_rows if row["item_class"] == "reference"]
    genome_ids = [row["id"] for row in manifest_rows if row["item_class"] == "genome"]
    id_to_idx = {seq_id: idx for idx, seq_id in enumerate(ids)}

    genome_ids_sorted = sorted(genome_ids, key=lambda seq_id: (manifest[seq_id]["source_cluster"], manifest[seq_id]["label"]))
    ref_ids_sorted = ref_ids[:]
    square_ids = ref_ids_sorted + genome_ids_sorted

    genome_ref_matrix = [
        [float(ident[id_to_idx[g], id_to_idx[r]]) for r in ref_ids_sorted]
        for g in genome_ids_sorted
    ]
    square_matrix = [
        [float(ident[id_to_idx[rid], id_to_idx[cid]]) for cid in square_ids]
        for rid in square_ids
    ]

    genome_row_labels = [manifest[g]["label"] for g in genome_ids_sorted]
    ref_col_labels = [manifest[r]["label"] for r in ref_ids_sorted]
    square_labels = [manifest[s]["label"] for s in square_ids]

    genome_ref_tsv = outdir / "genome_x_reference_identity.tsv"
    square_tsv = outdir / "panel_square_identity.tsv"
    write_tsv(
        genome_ref_tsv,
        [{"id": gid, **{ref_col_labels[i]: f"{row[i]:.6f}" for i in range(len(ref_col_labels))}} for gid, row in zip(genome_ids_sorted, genome_ref_matrix)],
        ["id", *ref_col_labels],
    )
    write_tsv(
        square_tsv,
        [{"id": sid, **{square_labels[i]: f"{row[i]:.6f}" for i in range(len(square_labels))}} for sid, row in zip(square_ids, square_matrix)],
        ["id", *square_labels],
    )

    write_html_heatmap(
        outdir / "heatmap_genome_x_reference.html",
        genome_row_labels,
        ref_col_labels,
        genome_ref_matrix,
        title="Genomes vs references MUSCLE identity",
        row_export_ids=genome_ids_sorted,
        enable_row_select=True,
        row_axis_label="genomes",
        col_axis_label="references",
    )
    maybe_heatmap(
        genome_ref_matrix,
        genome_row_labels,
        ref_col_labels,
        outdir / "heatmap_genome_x_reference.png",
        title="Genomes vs references MUSCLE identity",
    )

    write_html_heatmap(
        outdir / "heatmap_panel_square.html",
        square_labels,
        square_labels,
        square_matrix,
        title="Representative panel MUSCLE identity",
        row_axis_label="panel sequences",
        col_axis_label="panel sequences",
    )
    maybe_heatmap(
        square_matrix,
        square_labels,
        square_labels,
        outdir / "heatmap_panel_square.png",
        title="Representative panel MUSCLE identity",
    )


def build_neighbor_joining_tree(ids: List[str], dist: np.ndarray, manifest_rows: List[dict], outdir: Path) -> tuple[dict, np.ndarray]:
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    safe_names = [f"T{i + 1:05d}" for i in range(len(ids))]
    safe_to_id = dict(zip(safe_names, ids))
    id_to_manifest = {row["id"]: row for row in manifest_rows}
    lower_triangle = [
        [max(0.0, float(dist[i, j])) for j in range(i + 1)]
        for i in range(len(ids))
    ]
    matrix = DistanceMatrix(safe_names, lower_triangle)
    tree = DistanceTreeConstructor().nj(matrix)

    Phylo.write(tree, outdir / "neighbor_joining_tree.nwk", "newick")
    write_tsv(
        outdir / "neighbor_joining_tree_tip_map.tsv",
        [
            {
                "tree_tip": safe,
                "id": seq_id,
                "label": id_to_manifest[seq_id]["label"],
                "item_class": id_to_manifest[seq_id]["item_class"],
                "source_cluster": id_to_manifest[seq_id]["source_cluster"],
            }
            for safe, seq_id in safe_to_id.items()
        ],
        ["tree_tip", "id", "label", "item_class", "source_cluster"],
    )

    terminals = tree.get_terminals()
    terminal_by_name = {clade.name: clade for clade in terminals}
    leaf_y = {clade.name: idx for idx, clade in enumerate(terminals)}
    coords: dict[int, tuple[float, float]] = {}
    edges: list[dict] = []

    def branch_len(clade) -> float:
        value = clade.branch_length
        if value is None or value < 0:
            return 0.0
        return float(value)

    def layout(clade, x_pos: float) -> tuple[float, float]:
        node_id = id(clade)
        if clade.is_terminal():
            y_pos = float(leaf_y[clade.name])
            coords[node_id] = (x_pos, y_pos)
            return x_pos, y_pos

        child_positions = []
        for child in clade.clades:
            child_positions.append(layout(child, x_pos + branch_len(child)))
        y_values = [pos[1] for pos in child_positions]
        y_pos = sum(y_values) / len(y_values) if y_values else 0.0
        coords[node_id] = (x_pos, y_pos)

        if child_positions:
            edges.append(
                {
                    "kind": "connector",
                    "x1": x_pos,
                    "y1": min(y_values),
                    "x2": x_pos,
                    "y2": max(y_values),
                    "branch_length": 0.0,
                    "distance_from_root": x_pos,
                }
            )
        for child in clade.clades:
            child_x, child_y = coords[id(child)]
            child_id = ""
            if child.is_terminal() and child.name in safe_to_id:
                child_id = safe_to_id[child.name]
            edges.append(
                {
                    "kind": "branch",
                    "x1": x_pos,
                    "y1": child_y,
                    "x2": child_x,
                    "y2": child_y,
                    "branch_length": max(0.0, child_x - x_pos),
                    "distance_from_root": child_x,
                    "child_id": child_id,
                }
            )
        return x_pos, y_pos

    layout(tree.root, 0.0)
    leaves = []
    for clade in terminals:
        seq_id = safe_to_id[clade.name]
        row = id_to_manifest[seq_id]
        x_pos, y_pos = coords[id(clade)]
        leaves.append(
            {
                "id": seq_id,
                "label": row["label"],
                "item_class": row["item_class"],
                "length_bp": int(row.get("length_bp") or 0),
                "leiden_cluster": row["source_cluster"],
                "x": x_pos,
                "y": y_pos,
            }
        )
    max_x = max([leaf["x"] for leaf in leaves] + [edge["x1"] for edge in edges] + [edge["x2"] for edge in edges] + [0.0])
    max_y = max([leaf["y"] for leaf in leaves] + [0.0])
    if max_x <= 0:
        max_x = 1.0
    if max_y <= 0:
        max_y = 1.0
    patristic = np.zeros((len(ids), len(ids)), dtype=float)
    for i, safe_i in enumerate(safe_names):
        clade_i = terminal_by_name[safe_i]
        for j in range(i + 1, len(ids)):
            clade_j = terminal_by_name[safe_names[j]]
            value = float(tree.distance(clade_i, clade_j))
            patristic[i, j] = value
            patristic[j, i] = value
    return {"edges": edges, "leaves": leaves, "max_x": max_x, "max_y": max_y}, patristic


def html_template(payload_json: str, title: str) -> str:
    template = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>__TITLE__</title>
  <style>
    body {
      margin: 0;
      font-family: Helvetica, Arial, sans-serif;
      background: #f3f1eb;
      color: #1d1d1d;
    }
    .wrap {
      max-width: 1500px;
      margin: 0 auto;
      padding: 18px 20px 28px;
    }
    h1 {
      margin: 0 0 8px;
      font-size: 28px;
      line-height: 1.1;
    }
    .sub {
      margin: 0 0 14px;
      color: #555;
      font-size: 14px;
    }
    .toolbar {
      display: flex;
      gap: 10px;
      margin: 0 0 14px;
      flex-wrap: wrap;
    }
    .toggle {
      border: 1px solid #d0c8b9;
      background: white;
      color: #333;
      border-radius: 999px;
      padding: 8px 14px;
      cursor: pointer;
      font-size: 13px;
    }
    .toggle.active {
      background: #222;
      color: white;
      border-color: #222;
    }
    .slider-wrap {
      background: white;
      border: 1px solid #ddd6ca;
      border-radius: 12px;
      padding: 12px 14px;
      margin: 0 0 14px;
      font-size: 13px;
    }
    .slider-wrap label {
      display: block;
      font-weight: 700;
      margin-bottom: 6px;
    }
    .slider-wrap input {
      width: 100%;
    }
    .legend {
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin: 10px 0 14px;
      font-size: 13px;
    }
    .legend-item {
      display: flex;
      align-items: center;
      gap: 6px;
      background: white;
      border: 1px solid #ddd6ca;
      border-radius: 999px;
      padding: 5px 10px;
    }
    .swatch {
      width: 12px;
      height: 12px;
      border-radius: 999px;
      display: inline-block;
      border: 1px solid rgba(0,0,0,0.15);
    }
    .meta {
      margin: 0 0 12px;
      color: #555;
      font-size: 12px;
    }
    .metrics {
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 10px 12px;
      margin: 0 0 14px;
    }
    .metric {
      background: white;
      border: 1px solid #ddd6ca;
      border-radius: 10px;
      padding: 10px 12px;
    }
    .metric .label {
      display: block;
      color: #666;
      font-size: 11px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.03em;
      margin-bottom: 4px;
    }
    .metric .value {
      color: #1d1d1d;
      font-size: 18px;
      font-weight: 700;
    }
    .agreement-panel {
      display: none;
      background: white;
      border: 1px solid #ddd6ca;
      border-radius: 12px;
      padding: 12px 14px;
      margin: 0 0 14px;
      font-size: 12px;
      overflow-x: auto;
    }
    .agreement-panel.active {
      display: block;
    }
    .agreement-panel table {
      border-collapse: collapse;
      width: 100%;
      min-width: 980px;
    }
    .agreement-panel th,
    .agreement-panel td {
      border: 1px solid #e3ded4;
      padding: 5px 7px;
      text-align: right;
      white-space: nowrap;
    }
    .agreement-panel th:first-child,
    .agreement-panel td:first-child {
      text-align: left;
      position: sticky;
      left: 0;
      background: white;
      z-index: 1;
    }
    .agreement-note {
      color: #555;
      margin: 0 0 8px;
      line-height: 1.35;
    }
    .agreement-cell {
      cursor: pointer;
      font-weight: 700;
      transition: outline 0.12s ease, transform 0.12s ease;
    }
    .agreement-cell:hover {
      outline: 2px solid #222;
      outline-offset: -2px;
      transform: scale(1.02);
    }
    .agreement-best {
      color: #333;
      background: #f8f2e6;
      border: 1px solid #e2d2b8;
      border-radius: 10px;
      display: flex;
      gap: 10px;
      align-items: center;
      justify-content: space-between;
      margin: 8px 0 10px;
      padding: 9px 11px;
    }
    .agreement-best button {
      flex: 0 0 auto;
      padding: 6px 11px;
      font-size: 12px;
    }
    .exportbar {
      display: flex;
      gap: 10px;
      align-items: flex-end;
      flex-wrap: wrap;
      margin: 0 0 14px;
      padding: 12px 14px;
      background: white;
      border: 1px solid #ddd6ca;
      border-radius: 12px;
    }
    .exportbar label {
      display: block;
      font-size: 13px;
      font-weight: 700;
      margin-bottom: 6px;
    }
    .multiselect-wrap {
      position: relative;
      min-width: 200px;
    }
    .multiselect-toggle {
      width: 100%;
      text-align: left;
      padding: 6px 8px;
      border: 1px solid #ccc;
      border-radius: 4px;
      background: white;
      cursor: pointer;
      font-size: 13px;
    }
    .multiselect-menu {
      display: none;
      position: absolute;
      z-index: 100;
      top: 100%;
      left: 0;
      min-width: 240px;
      max-height: 260px;
      overflow-y: auto;
      background: white;
      border: 1px solid #ccc;
      border-radius: 6px;
      box-shadow: 0 4px 16px rgba(0,0,0,0.12);
      padding: 4px 0;
    }
    .multiselect-menu.open { display: block; }
    .multiselect-item, .multiselect-all {
      display: flex;
      align-items: center;
      gap: 8px;
      padding: 5px 12px;
      font-size: 13px;
      cursor: pointer;
      user-select: none;
    }
    .multiselect-item:hover, .multiselect-all:hover { background: #f5f2ec; }
    .multiselect-all {
      border-bottom: 1px solid #eee;
      font-weight: 700;
    }
    .exportbar button {
      padding: 8px 12px;
      border: 1px solid #d0c8b9;
      border-radius: 999px;
      background: white;
      cursor: pointer;
      font-size: 13px;
    }
    .exportbar .hint {
      color: #666;
      font-size: 12px;
      padding-bottom: 2px;
    }
    canvas {
      width: 100%;
      height: auto;
      background: #fffdf8;
      border: 1px solid #d8d2c6;
      border-radius: 16px;
      box-shadow: 0 8px 24px rgba(0,0,0,0.06);
    }
  </style>
</head>
<body>
  <div class="wrap">
    <h1>__TITLE__</h1>
    <p class="sub">Use the threshold sliders to inspect complete-linkage agglomerative clustering on the aligned panel and tree-distance clustering from the neighbor-joining tree. Toggle between alignment-based PCoA, alignment-based UMAP, AF UMAP, best-vs-second-best reference identity, and tree views.</p>
    <div class="toolbar">
      <button class="toggle active" id="btn-pcoa">PCoA</button>
      <button class="toggle" id="btn-umap">Distance UMAP</button>
      <button class="toggle" id="btn-afumap">AF UMAP</button>
      <button class="toggle" id="btn-best">Best vs 2nd</button>
      <button class="toggle" id="btn-tree">NJ Tree</button>
      <button class="toggle" id="btn-agreement">Cluster Agreement</button>
    </div>
    <div class="toolbar">
      <button class="toggle active" id="btn-ref-include">Refs include self</button>
      <button class="toggle" id="btn-ref-exclude">Refs exclude self</button>
    </div>
    <div class="slider-wrap">
      <label for="novel-threshold-input">Best-vs-2nd best cutoff (%)</label>
      <input id="novel-threshold-input" type="number" min="0" max="100" step="0.1" style="width:120px;">
    </div>
    <div class="toolbar">
      <button class="toggle active" id="btn-leiden">Color by Leiden</button>
      <button class="toggle" id="btn-agg">Color by Agglomerative</button>
      <button class="toggle" id="btn-tree-cluster">Color by Tree</button>
    </div>
    <div class="toolbar">
      <button class="toggle" id="btn-scale-full">Full scale</button>
      <button class="toggle active" id="btn-scale-robust">Robust scale</button>
    </div>
    <div class="slider-wrap">
      <label for="threshold">Complete-linkage threshold: <span id="threshold-value"></span></label>
      <input id="threshold" type="range" min="0" max="0" step="1">
    </div>
    <div class="slider-wrap">
      <label for="tree-threshold">NJ tree-distance threshold: <span id="tree-threshold-value"></span></label>
      <input id="tree-threshold" type="range" min="0" max="0" step="1">
    </div>
    <div class="slider-wrap" id="af-state-wrap">
      <label>AF / Leiden state: <span id="af-state-value"></span></label>
      <div style="display:grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap:12px;">
        <label for="af-kmer" style="font-weight:700; margin:0;">k-mer: <span id="af-kmer-value"></span>
          <input id="af-kmer" type="range" min="0" max="0" step="1">
        </label>
        <label for="af-neighbors" style="font-weight:700; margin:0;">n_neighbors: <span id="af-neighbors-value"></span>
          <input id="af-neighbors" type="range" min="0" max="0" step="1">
        </label>
        <label for="af-min-dist" style="font-weight:700; margin:0;">min_dist: <span id="af-min-dist-value"></span>
          <input id="af-min-dist" type="range" min="0" max="0" step="1">
        </label>
      </div>
    </div>
    <div class="slider-wrap">
      <label>Length filter (nt)</label>
      <div style="display:flex; gap:14px; align-items:center; flex-wrap:wrap;">
        <label for="min-length" style="font-weight:400; margin:0;">min:
          <input id="min-length" type="number" step="1" style="margin-left:6px; width:110px;">
        </label>
        <label for="max-length" style="font-weight:400; margin:0;">max:
          <input id="max-length" type="number" step="1" style="margin-left:6px; width:110px;">
        </label>
      </div>
    </div>
    <div class="exportbar">
      <div>
        <label>Export Clusters</label>
        <div class="multiselect-wrap">
          <button type="button" class="multiselect-toggle" id="review-export-toggle">
            <span id="review-export-label">Select clusters ▾</span>
          </button>
          <div class="multiselect-menu" id="review-export-menu">
            <label class="multiselect-all"><input type="checkbox" id="review-export-all"> All clusters</label>
            <div id="review-export-items"></div>
          </div>
        </div>
      </div>
      <button id="review-export-ids">Download IDs</button>
      <button id="review-export-cmd">Download FASTA Command</button>
      <div class="hint">Exports the active clustering scheme at the current threshold.</div>
    </div>
    <div class="meta" id="meta"></div>
    <div class="metrics" id="metrics"></div>
    <div class="agreement-panel" id="agreement-panel">
      <p class="agreement-note">Agreement compares every precomputed AF/Leiden state against every distance-derived clustering state on sampled genomes only. Higher NMI indicates greater information overlap between partitions.</p>
      <div class="toolbar">
        <button class="toggle active" id="btn-agree-agg">Against Agglomerative</button>
        <button class="toggle" id="btn-agree-tree">Against Tree</button>
      </div>
      <div id="agreement-content"></div>
    </div>
    <div class="legend" id="legend"></div>
    <canvas id="plot" width="1450" height="980"></canvas>
  </div>
  <script>
    const payload = __PAYLOAD__;
    const canvas = document.getElementById('plot');
    const ctx = canvas.getContext('2d');
    const legend = document.getElementById('legend');
    const meta = document.getElementById('meta');
    const metrics = document.getElementById('metrics');
    const btnPcoa = document.getElementById('btn-pcoa');
    const btnUmap = document.getElementById('btn-umap');
    const btnAfUmap = document.getElementById('btn-afumap');
    const btnBest = document.getElementById('btn-best');
    const btnTree = document.getElementById('btn-tree');
    const btnAgreement = document.getElementById('btn-agreement');
    const agreementPanel = document.getElementById('agreement-panel');
    const agreementContent = document.getElementById('agreement-content');
    const btnAgreeAgg = document.getElementById('btn-agree-agg');
    const btnAgreeTree = document.getElementById('btn-agree-tree');
    const btnRefInclude = document.getElementById('btn-ref-include');
    const btnRefExclude = document.getElementById('btn-ref-exclude');
    const btnLeiden = document.getElementById('btn-leiden');
    const btnAgg = document.getElementById('btn-agg');
    const btnTreeCluster = document.getElementById('btn-tree-cluster');
    const btnScaleFull = document.getElementById('btn-scale-full');
    const btnScaleRobust = document.getElementById('btn-scale-robust');
    const novelThresholdInput = document.getElementById('novel-threshold-input');
    const thresholdSlider = document.getElementById('threshold');
    const thresholdValue = document.getElementById('threshold-value');
    const treeThresholdSlider = document.getElementById('tree-threshold');
    const treeThresholdValue = document.getElementById('tree-threshold-value');
    const afStateWrap = document.getElementById('af-state-wrap');
    const afStateValue = document.getElementById('af-state-value');
    const afKmerSlider = document.getElementById('af-kmer');
    const afNeighborsSlider = document.getElementById('af-neighbors');
    const afMinDistSlider = document.getElementById('af-min-dist');
    const afKmerValue = document.getElementById('af-kmer-value');
    const afNeighborsValue = document.getElementById('af-neighbors-value');
    const afMinDistValue = document.getElementById('af-min-dist-value');
    const minLengthInput = document.getElementById('min-length');
    const maxLengthInput = document.getElementById('max-length');
    const reviewExportToggle = document.getElementById('review-export-toggle');
    const reviewExportMenu = document.getElementById('review-export-menu');
    const reviewExportLabel = document.getElementById('review-export-label');
    const reviewExportAll = document.getElementById('review-export-all');
    const reviewExportItems = document.getElementById('review-export-items');
    const reviewExportIdsBtn = document.getElementById('review-export-ids');
    const reviewExportCmdBtn = document.getElementById('review-export-cmd');
    const thresholds = payload.thresholds;
    const afStates = payload.af_states || [];
    const defaultAfState = afStates[payload.default_af_state_index || 0] || null;
    const afKmerValues = Array.from(new Set(afStates.map(s => Number(s.kmer)))).sort((a, b) => a - b);
    const afNeighborsValues = Array.from(new Set(afStates.map(s => Number(s.neighbors)))).sort((a, b) => a - b);
    const afMinDistValues = Array.from(new Set(afStates.map(s => Number(s.min_dist)))).sort((a, b) => a - b);
    novelThresholdInput.value = String(payload.novel_threshold);
    thresholdSlider.max = String(thresholds.length - 1);
    thresholdSlider.value = String(payload.default_threshold_index);
    treeThresholdSlider.max = String(thresholds.length - 1);
    treeThresholdSlider.value = String(payload.default_tree_threshold_index || payload.default_threshold_index);
    afKmerSlider.max = String(Math.max(0, afKmerValues.length - 1));
    afNeighborsSlider.max = String(Math.max(0, afNeighborsValues.length - 1));
    afMinDistSlider.max = String(Math.max(0, afMinDistValues.length - 1));
    afKmerSlider.value = String(Math.max(0, afKmerValues.findIndex(v => defaultAfState && v === Number(defaultAfState.kmer))));
    afNeighborsSlider.value = String(Math.max(0, afNeighborsValues.findIndex(v => defaultAfState && v === Number(defaultAfState.neighbors))));
    afMinDistSlider.value = String(Math.max(0, afMinDistValues.findIndex(v => defaultAfState && v === Number(defaultAfState.min_dist))));
    if (!afStates.length) {
      afStateWrap.style.display = 'none';
      btnAfUmap.style.display = 'none';
    }
    const allLengths = payload.points.map(p => Number(p.length_bp || 0)).filter(v => Number.isFinite(v) && v > 0);
    const globalMinLength = allLengths.length ? Math.min(...allLengths) : 0;
    const globalMaxLength = allLengths.length ? Math.max(...allLengths) : 0;
    minLengthInput.value = String(globalMinLength);
    maxLengthInput.value = String(globalMaxLength);

    let mode = 'pcoa';
    let referenceBestMode = 'include';
    let scheme = 'leiden_cluster';
    let scaleMode = 'robust';
    let agreementVisible = false;
    let agreementTarget = 'agglomerative';
    let selectedAfStateId = defaultAfState ? defaultAfState.state_id : null;
    let hover = null;
    let hitboxes = [];
    const box = {x: 90, y: 90, w: 1220, h: 800};
    const pointById = new Map(payload.points.map(p => [p.id, p]));

    function currentThreshold() {
      return thresholds[Number(thresholdSlider.value)];
    }
    function currentTreeThreshold() {
      return thresholds[Number(treeThresholdSlider.value)];
    }
    function currentNovelThreshold() {
      const value = Number(novelThresholdInput.value);
      if (!Number.isFinite(value)) return Number(payload.novel_threshold);
      return Math.max(0, Math.min(100, value));
    }
    function currentAfValues() {
      return {
        kmer: afKmerValues[Number(afKmerSlider.value)],
        neighbors: afNeighborsValues[Number(afNeighborsSlider.value)],
        min_dist: afMinDistValues[Number(afMinDistSlider.value)],
      };
    }
    function currentAfState() {
      const vals = currentAfValues();
      const selectedState = afStates.find(s => s.state_id === selectedAfStateId);
      if (selectedState &&
          Number(selectedState.kmer) === vals.kmer &&
          Number(selectedState.neighbors) === vals.neighbors &&
          Number(selectedState.min_dist) === vals.min_dist) {
        return selectedState;
      }
      return afStates.find(s =>
        Number(s.kmer) === vals.kmer &&
        Number(s.neighbors) === vals.neighbors &&
        Number(s.min_dist) === vals.min_dist
      ) || null;
    }
    function pointCluster(point) {
      if (scheme === 'leiden_cluster') {
        const state = currentAfState();
        if (state && state.clusters[point.id]) return state.clusters[point.id];
        return point.leiden_cluster;
      }
      if (scheme === 'tree_cluster') return point.tree_clusters[currentTreeThreshold()];
      return point.agglomerative[currentThreshold()];
    }
    function currentView() {
      if (scheme === 'leiden_cluster') return payload.cluster_views.leiden_cluster[mode];
      if (scheme === 'tree_cluster') return payload.cluster_views.tree[currentTreeThreshold()][mode];
      return payload.cluster_views.agglomerative[currentThreshold()][mode];
    }
    function currentColorView() {
      if (scheme === 'leiden_cluster') {
        const state = currentAfState();
        if (state) return {cluster_order: state.cluster_order, cluster_colors: state.cluster_colors};
        return payload.cluster_views.leiden_cluster.pcoa;
      }
      if (scheme === 'tree_cluster') return payload.cluster_views.tree[currentTreeThreshold()].pcoa;
      return payload.cluster_views.agglomerative[currentThreshold()].pcoa;
    }
    function currentLengthBounds() {
      const minVal = Number(minLengthInput.value);
      const maxVal = Number(maxLengthInput.value);
      const minLen = Number.isFinite(minVal) ? minVal : globalMinLength;
      const maxLen = Number.isFinite(maxVal) ? maxVal : globalMaxLength;
      return {minLen, maxLen};
    }
    function filteredPoints() {
      const {minLen, maxLen} = currentLengthBounds();
      return payload.points.filter(p => {
        const len = Number(p.length_bp || 0);
        return len >= minLen && len <= maxLen;
      });
    }
    function pointCoords(point) {
      if (mode === 'tree') return point.pcoa;
      if (mode === 'afumap') {
        const state = currentAfState();
        if (state && state.coords[point.id]) return state.coords[point.id];
        return [NaN, NaN];
      }
      if (mode !== 'best2') return point[mode];
      if (point.item_class !== 'reference') return point.best2;
      return referenceBestMode === 'include' ? point.best2_ref_include : point.best2_ref_exclude;
    }
    function drawStar(x, y, fill) {
      ctx.save();
      ctx.translate(x, y);
      ctx.fillStyle = fill;
      ctx.strokeStyle = '#111';
      ctx.lineWidth = 0.8;
      ctx.beginPath();
      for (let k = 0; k < 5; k++) {
        const outer = -Math.PI / 2 + k * (2 * Math.PI / 5);
        const inner = outer + Math.PI / 5;
        const ox = Math.cos(outer) * 8;
        const oy = Math.sin(outer) * 8;
        const ix = Math.cos(inner) * 3.8;
        const iy = Math.sin(inner) * 3.8;
        if (k === 0) ctx.moveTo(ox, oy);
        else ctx.lineTo(ox, oy);
        ctx.lineTo(ix, iy);
      }
      ctx.closePath();
      ctx.fill();
      ctx.stroke();
      ctx.restore();
    }
    function quantile(sortedVals, q) {
      if (!sortedVals.length) return 0;
      const pos = (sortedVals.length - 1) * q;
      const lo = Math.floor(pos);
      const hi = Math.ceil(pos);
      if (lo === hi) return sortedVals[lo];
      const frac = pos - lo;
      return sortedVals[lo] * (1 - frac) + sortedVals[hi] * frac;
    }
    function bounds(points) {
      const coords = points.map(p => pointCoords(p)).filter(c => Number.isFinite(c[0]) && Number.isFinite(c[1]));
      if (!coords.length) return {minX: 0, maxX: 1, minY: 0, maxY: 1};
      const xs = coords.map(c => c[0]).sort((a, b) => a - b);
      const ys = coords.map(c => c[1]).sort((a, b) => a - b);
      const full = {minX: xs[0], maxX: xs[xs.length - 1], minY: ys[0], maxY: ys[ys.length - 1]};
      if (scaleMode === 'full' || xs.length < 10) return full;
      const robust = {
        minX: quantile(xs, 0.01),
        maxX: quantile(xs, 0.99),
        minY: quantile(ys, 0.01),
        maxY: quantile(ys, 0.99),
      };
      if (robust.minX === robust.maxX) {
        robust.minX = full.minX;
        robust.maxX = full.maxX;
      }
      if (robust.minY === robust.maxY) {
        robust.minY = full.minY;
        robust.maxY = full.maxY;
      }
      return robust;
    }
    function sx(x, b) {
      const clipped = Math.min(b.maxX, Math.max(b.minX, x));
      return box.x + ((clipped - b.minX) / (b.maxX - b.minX || 1)) * box.w;
    }
    function sy(y, b) {
      const clipped = Math.min(b.maxY, Math.max(b.minY, y));
      return box.y + box.h - ((clipped - b.minY) / (b.maxY - b.minY || 1)) * box.h;
    }
    function convexHull(pts) {
      const points = Array.from(new Set(pts.map(p => `${p[0]}|${p[1]}`))).map(s => s.split('|').map(Number));
      points.sort((a, b) => (a[0] - b[0]) || (a[1] - b[1]));
      if (points.length <= 2) return points;
      function cross(o, a, b) {
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
      }
      const lower = [];
      for (const p of points) {
        while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) lower.pop();
        lower.push(p);
      }
      const upper = [];
      for (let i = points.length - 1; i >= 0; i--) {
        const p = points[i];
        while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) upper.pop();
        upper.push(p);
      }
      lower.pop();
      upper.pop();
      return lower.concat(upper);
    }
    function computeView(points) {
      const colorView = currentColorView();
      const baseOrder = colorView.cluster_order || [];
      const baseColors = colorView.cluster_colors || {};
      const grouped = new Map();
      let visibleReferenceCount = 0;
      let totalGenomes = 0;
      for (const p of points) {
        if (p.item_class === 'reference') {
          visibleReferenceCount += 1;
          continue;
        }
        const coord = pointCoords(p);
        if (!Number.isFinite(coord[0]) || !Number.isFinite(coord[1])) continue;
        const cid = pointCluster(p);
        if (!grouped.has(cid)) grouped.set(cid, []);
        grouped.get(cid).push(coord);
        totalGenomes += 1;
      }
      const cluster_order = baseOrder.length
        ? baseOrder.filter(cid => grouped.has(cid))
        : Array.from(grouped.keys()).sort();
      const cluster_genome_counts = {};
      const cluster_pct = {};
      const labels = [];
      const hulls = [];
      cluster_order.forEach((cid, idx) => {
        const pts = grouped.get(cid) || [];
        const count = pts.length;
        cluster_genome_counts[cid] = count;
        cluster_pct[cid] = totalGenomes ? (100.0 * count / totalGenomes) : 0.0;
        if (pts.length) {
          const xs = pts.map(p => p[0]);
          const ys = pts.map(p => p[1]);
          labels.push({
            cluster: cid,
            x: xs.reduce((a, b) => a + b, 0) / xs.length,
            y: ys.reduce((a, b) => a + b, 0) / ys.length,
            genomes: count,
            pct: cluster_pct[cid],
            offset_x: ((idx % 3) - 1) * 18,
            offset_y: ((idx % 4) - 1.5) * 12,
          });
          if (pts.length >= 3) {
            const hull = convexHull(pts);
            if (hull.length >= 3) hulls.push({cluster: cid, points: hull});
          }
        }
      });
      return {cluster_order, cluster_colors: baseColors, cluster_genome_counts, cluster_pct, labels, hulls, totalGenomes, visibleReferenceCount};
    }
    function setLegend() {
      legend.innerHTML = '';
      const points = filteredPoints();
      const view = computeView(points);
      for (const cid of view.cluster_order) {
        const item = document.createElement('div');
        item.className = 'legend-item';
        const count = view.cluster_genome_counts[cid] || 0;
        const pct = view.cluster_pct[cid] || 0;
        item.innerHTML = `<span class="swatch" style="background:${view.cluster_colors[cid]}"></span>${cid} (${count}, ${pct.toFixed(1)}%)`;
        legend.appendChild(item);
      }
      const ref = document.createElement('div');
      ref.className = 'legend-item';
      ref.innerHTML = `<span style="font-size:14px;color:#000;">★</span> references (${view.visibleReferenceCount})`;
      legend.appendChild(ref);
      const thr = currentThreshold();
      const treeThr = currentTreeThreshold();
      let detail = `complete-linkage agglomerative clusters @ ${thr}`;
      if (scheme === 'leiden_cluster') detail = 'AF-widget Leiden clusters';
      if (scheme === 'tree_cluster') detail = `NJ tree-distance clusters @ ${treeThr}`;
      const afDetail = currentAfState() ? ` | AF state ${currentAfState().state_id} (${currentAfState().cluster_count} Leiden clusters)` : '';
      const refDetail = mode === 'best2' ? (referenceBestMode === 'include' ? 'refs use self as best' : 'refs use nearest non-self refs') : 'reference best-vs-second convention hidden';
      const {minLen, maxLen} = currentLengthBounds();
      const scaleDetail = `${scaleMode} scale`;
      meta.textContent = `${mode.toUpperCase()} | ${detail}${afDetail} | ${scaleDetail} | ${refDetail} | best-ref cutoff ${currentNovelThreshold()}% | visible points ${points.length}/${payload.points.length} | length ${minLen}-${maxLen} nt`;
    }
    function metricCard(label, value, tooltip) {
      const tip = tooltip ? ` title="${tooltip}"` : '';
      return `<div class="metric"${tip}><span class="label">${label}</span><span class="value">${value}</span></div>`;
    }
    function setMetrics() {
      const thr = currentThreshold();
      const treeThr = currentTreeThreshold();
      const thrMetrics = payload.validation.thresholds[thr] || {};
      const treeMetrics = payload.validation.tree_thresholds[treeThr] || {};
      const sample = payload.validation.sample || {};
      const align = payload.validation.align || {};
      const review = payload.validation.review || {};
      const afState = currentAfState();
      const visibleNovelGenomes = filteredPoints().filter(p =>
        p.item_class === 'genome' && Number(p.best_identity) < currentNovelThreshold()
      ).length;
      const minLineageSize = (sample.N && sample.sample_size_total)
        ? Math.ceil(sample.N / sample.sample_size_total) : null;
      metrics.innerHTML = [
        metricCard('Panel Fraction', sample.panel_fraction_pct !== undefined ? `${sample.panel_fraction_pct.toFixed(2)}%` : 'NA'),
        metricCard('Min Lineage Size', minLineageSize !== null ? `≥${minLineageSize} seqs` : 'NA',
          'Minimum number of sequences a lineage needs in the full dataset to be reliably captured in the panel (= ⌈N / panel_size⌉). Lineages below this threshold may be missed by proportional sampling.'),
        metricCard('Pearson r', sample.pearson_r_full_vs_panel_cluster_proportions !== undefined ? sample.pearson_r_full_vs_panel_cluster_proportions.toFixed(4) : 'NA'),
        metricCard('Runtime', align.runtime_seconds !== undefined ? `${(align.runtime_seconds / 3600).toFixed(2)} h` : 'NA'),
        metricCard('PCoA Var', (align.pcoa_axis1_variance_explained_pct !== undefined && align.pcoa_axis2_variance_explained_pct !== undefined) ? `${align.pcoa_axis1_variance_explained_pct.toFixed(1)} / ${align.pcoa_axis2_variance_explained_pct.toFixed(1)}%` : 'NA'),
        metricCard('Agg Clusters', thrMetrics.cluster_count !== undefined ? thrMetrics.cluster_count : 'NA'),
        metricCard('Tree Clusters', treeMetrics.cluster_count !== undefined ? treeMetrics.cluster_count : 'NA'),
        metricCard('Leiden Clusters', afState ? afState.cluster_count : 'NA'),
        metricCard('NMI', thrMetrics.nmi !== undefined ? thrMetrics.nmi.toFixed(3) : 'NA'),
        metricCard('Novel < Cutoff', review.novel_genome_count !== undefined ? visibleNovelGenomes : 'NA'),
      ].join('');
    }
    function agreementColor(value) {
      if (!Number.isFinite(value)) return '#fff';
      const v = Math.max(0, Math.min(1, value));
      return `rgba(76, 120, 168, ${0.10 + 0.70 * v})`;
    }
    function thresholdIndexFor(value) {
      const numeric = Number(value);
      const idx = thresholds.findIndex(t => Math.abs(Number(t) - numeric) < 1e-9);
      return idx >= 0 ? idx : 0;
    }
    function setAfSlidersFromState(stateId) {
      const state = afStates.find(s => s.state_id === stateId);
      if (!state) return false;
      const kIdx = afKmerValues.findIndex(v => v === Number(state.kmer));
      const nIdx = afNeighborsValues.findIndex(v => v === Number(state.neighbors));
      const dIdx = afMinDistValues.findIndex(v => Math.abs(v - Number(state.min_dist)) < 1e-9);
      if (kIdx >= 0) afKmerSlider.value = String(kIdx);
      if (nIdx >= 0) afNeighborsSlider.value = String(nIdx);
      if (dIdx >= 0) afMinDistSlider.value = String(dIdx);
      selectedAfStateId = state.state_id;
      return true;
    }
    function agreementRecord(target, stateId, threshold) {
      const numericThreshold = Number(threshold);
      return (payload.cluster_agreement || []).find(r =>
        r.target === target &&
        r.af_state === stateId &&
        Math.abs(Number(r.threshold) - numericThreshold) < 1e-9
      ) || null;
    }
    function bestAgreementRecord(target) {
      const rows = (payload.cluster_agreement || []).filter(r => r.target === target);
      if (!rows.length) return null;
      return rows.reduce((best, row) => Number(row.nmi) > Number(best.nmi) ? row : best, rows[0]);
    }
    function applyAgreementRecord(record) {
      if (!record) return;
      setAfSlidersFromState(record.af_state);
      if (record.target === 'tree') {
        treeThresholdSlider.value = String(thresholdIndexFor(record.threshold));
        setScheme('tree_cluster');
      } else {
        thresholdSlider.value = String(thresholdIndexFor(record.threshold));
        setScheme('agglomerative');
      }
      setMode('afumap');
    }
    function setAgreementPanel() {
      agreementPanel.classList.toggle('active', agreementVisible);
      btnAgreement.classList.toggle('active', agreementVisible);
      btnAgreeAgg.classList.toggle('active', agreementTarget === 'agglomerative');
      btnAgreeTree.classList.toggle('active', agreementTarget === 'tree');
      if (!agreementVisible) return;
      const rows = (payload.cluster_agreement || []).filter(r => r.target === agreementTarget);
      if (!rows.length) {
        agreementContent.innerHTML = '<p>No cluster-agreement metrics are available for this widget.</p>';
        return;
      }
      const thresholdsForRows = Array.from(new Set(rows.map(r => r.threshold)));
      const stateOrder = afStates.map(s => s.state_id);
      const hasResolution = rows.some(r => r.leiden_resolution !== undefined && r.leiden_resolution !== null);
      const byState = new Map();
      for (const row of rows) {
        if (!byState.has(row.af_state)) byState.set(row.af_state, []);
        byState.get(row.af_state).push(row);
      }
      const html = [];
      const best = bestAgreementRecord(agreementTarget);
      if (best) {
        html.push(`<div class="agreement-best"><span><strong>Best ${agreementTarget} NMI:</strong> ${Number(best.nmi).toFixed(3)} at AF state ${best.af_state}, threshold ${best.threshold}, ${best.target_cluster_count} ${agreementTarget} clusters. Click any NMI cell to apply that state to the graph.</span><button class="toggle agreement-apply" data-target="${best.target}" data-af-state="${best.af_state}" data-threshold="${best.threshold}">Apply best</button></div>`);
      }
      html.push('<table><thead><tr>');
      html.push('<th>AF state</th><th>k</th><th>n</th><th>dist</th>');
      if (hasResolution) html.push('<th>res</th>');
      html.push('<th>Leiden K</th>');
      for (const thr of thresholdsForRows) html.push(`<th>${thr}</th>`);
      html.push('</tr></thead><tbody>');
      for (const stateId of stateOrder) {
        const values = byState.get(stateId) || [];
        if (!values.length) continue;
        const first = values[0];
        const byThreshold = new Map(values.map(v => [v.threshold, v]));
        html.push('<tr>');
        html.push(`<td>${stateId}</td><td>${first.kmer}</td><td>${first.neighbors}</td><td>${first.min_dist}</td>`);
        if (hasResolution) html.push(`<td>${first.leiden_resolution !== undefined && first.leiden_resolution !== null ? first.leiden_resolution : 'NA'}</td>`);
        html.push(`<td>${first.af_cluster_count}</td>`);
        for (const thr of thresholdsForRows) {
          const rec = byThreshold.get(thr);
          const value = rec ? Number(rec.nmi) : NaN;
          const isBest = best && rec && rec.af_state === best.af_state && Math.abs(Number(rec.threshold) - Number(best.threshold)) < 1e-9;
          const title = rec ? `Apply AF ${rec.af_state} against ${rec.target} threshold ${rec.threshold}; target clusters: ${rec.target_cluster_count}` : '';
          const attrs = rec ? ` class="agreement-cell${isBest ? ' best-nmi' : ''}" data-target="${rec.target}" data-af-state="${rec.af_state}" data-threshold="${rec.threshold}"` : '';
          const border = isBest ? '; box-shadow: inset 0 0 0 2px #222' : '';
          html.push(`<td${attrs} title="${title}" style="background:${agreementColor(value)}${border}">${Number.isFinite(value) ? value.toFixed(3) : 'NA'}</td>`);
        }
        html.push('</tr>');
      }
      html.push('</tbody></table>');
      agreementContent.innerHTML = html.join('');
    }
    function drawTooltip(item, x, y) {
      const thr = currentThreshold();
      const lines = [
        item.label || item.id,
        `id: ${item.id}`,
        `type: ${item.item_class}`,
        `length: ${item.length_bp} bp`,
        `Leiden: ${item.leiden_cluster}`,
        `Agg (${thr}): ${item.agglomerative[thr]}`,
        `Tree (${currentTreeThreshold()}): ${item.tree_clusters[currentTreeThreshold()]}`,
        (mode === 'tree' && item.tree_root_distance !== undefined) ? `root distance: ${item.tree_root_distance.toFixed(5)} (${(item.tree_root_distance * 100).toFixed(2)}%)` : '',
        item.best_ref ? `best: ${item.best_ref} (${item.best_identity.toFixed(2)}%)` : '',
        (mode === 'best2' && item.item_class === 'reference') ? (referenceBestMode === 'include' ? 'refs plotted with self-match on x' : 'refs plotted without self-match') : '',
      ].filter(Boolean);
      drawTooltipBox(lines, x, y);
    }
    function drawEdgeTooltip(edge, x, y) {
      const lines = [
        edge.kind === 'branch' ? 'NJ branch' : 'NJ connector',
        edge.child_id ? `to: ${edge.child_id}` : '',
        `branch length: ${edge.branch_length.toFixed(5)} (${(edge.branch_length * 100).toFixed(2)}%)`,
        `root distance: ${edge.distance_from_root.toFixed(5)} (${(edge.distance_from_root * 100).toFixed(2)}%)`,
      ].filter(Boolean);
      drawTooltipBox(lines, x, y);
    }
    function drawTooltipBox(lines, x, y) {
      const pad = 8;
      ctx.font = '12px Helvetica';
      const width = Math.max(...lines.map(line => ctx.measureText(line).width)) + pad * 2;
      const height = lines.length * 16 + pad * 2;
      let tx = x + 14;
      let ty = y - 10;
      if (tx + width > canvas.width - 20) tx = x - width - 14;
      if (ty + height > canvas.height - 20) ty = canvas.height - height - 20;
      if (ty < 20) ty = 20;
      ctx.fillStyle = 'rgba(255,255,255,0.96)';
      ctx.strokeStyle = 'rgba(0,0,0,0.18)';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.roundRect(tx, ty, width, height, 10);
      ctx.fill();
      ctx.stroke();
      ctx.fillStyle = '#111';
      lines.forEach((line, i) => ctx.fillText(line, tx + pad, ty + pad + 12 + i * 16));
    }
    function pointSegmentDistance(px, py, x1, y1, x2, y2) {
      const dx = x2 - x1;
      const dy = y2 - y1;
      const denom = dx * dx + dy * dy;
      if (denom === 0) {
        const ex = px - x1;
        const ey = py - y1;
        return Math.sqrt(ex * ex + ey * ey);
      }
      let t = ((px - x1) * dx + (py - y1) * dy) / denom;
      t = Math.max(0, Math.min(1, t));
      const cx = x1 + t * dx;
      const cy = y1 + t * dy;
      const ex = px - cx;
      const ey = py - cy;
      return Math.sqrt(ex * ex + ey * ey);
    }
    function renderTree(points, view) {
      thresholdValue.textContent = currentThreshold();
      const visible = new Set(points.map(p => p.id));
      hitboxes = [];
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.fillStyle = '#fffdf8';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.strokeStyle = '#ddd6ca';
      ctx.lineWidth = 1;
      ctx.strokeRect(box.x, box.y, box.w, box.h);
      const tree = payload.tree;
      const left = box.x + 20;
      const top = box.y + 16;
      const width = box.w - 170;
      const height = box.h - 32;
      function tx(x) {
        return left + (x / (tree.max_x || 1)) * width;
      }
      function ty(y) {
        return top + (y / (tree.max_y || 1)) * height;
      }
      ctx.strokeStyle = '#b9b3aa';
      ctx.globalAlpha = 0.55;
      ctx.lineWidth = 0.75;
      for (const edge of tree.edges) {
        const x1 = tx(edge.x1);
        const y1 = ty(edge.y1);
        const x2 = tx(edge.x2);
        const y2 = ty(edge.y2);
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
        if (edge.kind === 'branch') {
          hitboxes.push({type: 'edge', x1, y1, x2, y2, r: 5, edge});
        }
      }
      ctx.globalAlpha = 1;
      ctx.font = '10px Helvetica';
      for (const leaf of tree.leaves) {
        if (!visible.has(leaf.id)) continue;
        const point = pointById.get(leaf.id);
        if (!point) continue;
        const cid = pointCluster(point);
        const color = view.cluster_colors[cid] || '#A0A0A0';
        const x = tx(leaf.x);
        const y = ty(leaf.y);
        if (leaf.item_class === 'reference') {
          drawStar(x, y, '#000');
          point.tree_root_distance = leaf.x;
          hitboxes.push({type: 'point', x, y, r: 10, item: point});
        } else {
          ctx.beginPath();
          ctx.fillStyle = color;
          ctx.globalAlpha = 0.9;
          ctx.arc(x, y, 2.7, 0, Math.PI * 2);
          ctx.fill();
          ctx.globalAlpha = 1;
          point.tree_root_distance = leaf.x;
          hitboxes.push({type: 'point', x, y, r: 5, item: point});
        }
      }
      ctx.fillStyle = '#666';
      ctx.font = '13px Helvetica';
      ctx.fillText('Neighbor-joining tree from MUSCLE pairwise distance matrix. Hover tips or branches for tree distances.', box.x + 20, box.y + box.h + 46);
      if (hover && hover.type === 'edge') drawEdgeTooltip(hover.edge, hover.x, hover.y);
      if (hover && hover.type !== 'edge') drawTooltip(hover.item, hover.x, hover.y);
    }
    function render() {
      const points = filteredPoints();
      const view = computeView(points);
      setMetrics();
      setAgreementPanel();
      if (afStates.length) {
        const state = currentAfState();
        const vals = currentAfValues();
        afKmerValue.textContent = String(vals.kmer);
        afNeighborsValue.textContent = String(vals.neighbors);
        afMinDistValue.textContent = String(vals.min_dist);
        afStateValue.textContent = state ? `${state.state_id} (${state.cluster_count} clusters)` : 'state not available';
      }
      if (mode === 'tree') {
        treeThresholdValue.textContent = currentTreeThreshold();
        renderTree(points, view);
        return;
      }
      if (!points.length) {
        thresholdValue.textContent = currentThreshold();
        treeThresholdValue.textContent = currentTreeThreshold();
        hitboxes = [];
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = '#fffdf8';
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        ctx.strokeStyle = '#ddd6ca';
        ctx.lineWidth = 1;
        ctx.strokeRect(box.x, box.y, box.w, box.h);
        ctx.fillStyle = '#444';
        ctx.font = '18px Helvetica';
        ctx.fillText('No points pass the current length filter.', box.x + 30, box.y + 40);
        return;
      }
      const b = bounds(points);
      thresholdValue.textContent = currentThreshold();
      treeThresholdValue.textContent = currentTreeThreshold();
      hitboxes = [];
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.fillStyle = '#fffdf8';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.strokeStyle = '#ddd6ca';
      ctx.lineWidth = 1;
      ctx.strokeRect(box.x, box.y, box.w, box.h);
      if (mode === 'best2') {
        const cutoff = currentNovelThreshold();
        const thrX = sx(cutoff, b);
        ctx.strokeStyle = '#666';
        ctx.setLineDash([6, 6]);
        ctx.beginPath();
        ctx.moveTo(thrX, box.y);
        ctx.lineTo(thrX, box.y + box.h);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.strokeStyle = '#333';
        ctx.beginPath();
        ctx.moveTo(thrX, box.y + box.h);
        ctx.lineTo(thrX, box.y + box.h + 9);
        ctx.stroke();
        const cutoffLabel = `cutoff ${cutoff.toFixed(1)}%`;
        ctx.font = '12px Helvetica';
        const labelWidth = ctx.measureText(cutoffLabel).width;
        const labelX = Math.max(box.x + 4, Math.min(box.x + box.w - labelWidth - 4, thrX - labelWidth / 2));
        ctx.fillStyle = '#222';
        ctx.fillText(cutoffLabel, labelX, box.y + box.h + 24);
      }
      for (const hull of view.hulls) {
        const color = view.cluster_colors[hull.cluster];
        ctx.fillStyle = color + '18';
        ctx.strokeStyle = color + 'AA';
        ctx.lineWidth = 1.1;
        ctx.beginPath();
        hull.points.forEach((p, idx) => {
          const x = sx(p[0], b);
          const y = sy(p[1], b);
          if (idx === 0) ctx.moveTo(x, y);
          else ctx.lineTo(x, y);
        });
        ctx.closePath();
        ctx.fill();
        ctx.stroke();
      }
      for (const p of points) {
        const coord = pointCoords(p);
        if (!Number.isFinite(coord[0]) || !Number.isFinite(coord[1])) continue;
        const x = sx(coord[0], b);
        const y = sy(coord[1], b);
        const cid = pointCluster(p);
        const color = view.cluster_colors[cid] || '#A0A0A0';
        if (p.item_class === 'reference') {
          drawStar(x, y, '#000');
          hitboxes.push({x, y, r: 10, item: p});
        } else {
          ctx.beginPath();
          ctx.fillStyle = color;
          ctx.globalAlpha = 0.75;
          ctx.arc(x, y, 3.1, 0, Math.PI * 2);
          ctx.fill();
          ctx.globalAlpha = 1;
          hitboxes.push({type: 'point', x, y, r: 6, item: p});
        }
      }
      ctx.font = '12px Helvetica';
      ctx.fillStyle = '#222';
      for (const label of view.labels) {
        const x = sx(label.x, b) + label.offset_x;
        const y = sy(label.y, b) + label.offset_y;
        ctx.fillText(`${label.cluster} (${label.genomes}, ${label.pct.toFixed(1)}%)`, x, y);
      }
      ctx.fillStyle = '#666';
      ctx.font = '13px Helvetica';
      const xlab = mode === 'best2' ? 'Best reference identity (%)' : (mode === 'pcoa' ? payload.pcoa_axis_labels.x : (mode === 'afumap' ? 'AF UMAP 1' : 'UMAP 1'));
      const ylab = mode === 'best2' ? 'Second-best reference identity (%)' : (mode === 'pcoa' ? payload.pcoa_axis_labels.y : (mode === 'afumap' ? 'AF UMAP 2' : 'UMAP 2'));
      ctx.fillText(xlab, box.x + box.w / 2 - 40, box.y + box.h + 46);
      ctx.save();
      ctx.translate(28, box.y + box.h / 2 + 30);
      ctx.rotate(-Math.PI / 2);
      ctx.fillText(ylab, 0, 0);
      ctx.restore();
      if (hover && hover.type === 'edge') drawEdgeTooltip(hover.edge, hover.x, hover.y);
      if (hover && hover.type !== 'edge') drawTooltip(hover.item, hover.x, hover.y);
    }
    canvas.addEventListener('mousemove', ev => {
      const rect = canvas.getBoundingClientRect();
      const x = (ev.clientX - rect.left) * (canvas.width / rect.width);
      const y = (ev.clientY - rect.top) * (canvas.height / rect.height);
      hover = null;
      for (const hb of hitboxes) {
        if (hb.type !== 'edge') {
          const dx = x - hb.x;
          const dy = y - hb.y;
          if ((dx * dx) + (dy * dy) <= hb.r * hb.r) {
            hover = {type: 'point', item: hb.item, x, y};
            break;
          }
        }
      }
      if (!hover) {
        for (const hb of hitboxes) {
          if (hb.type === 'edge' && pointSegmentDistance(x, y, hb.x1, hb.y1, hb.x2, hb.y2) <= hb.r) {
            hover = {type: 'edge', edge: hb.edge, x, y};
            break;
          }
        }
      }
      render();
    });
    canvas.addEventListener('mouseleave', () => { hover = null; render(); });
    function setMode(next) {
      mode = next;
      btnPcoa.classList.toggle('active', next === 'pcoa');
      btnUmap.classList.toggle('active', next === 'umap');
      btnAfUmap.classList.toggle('active', next === 'afumap');
      btnBest.classList.toggle('active', next === 'best2');
      btnTree.classList.toggle('active', next === 'tree');
      btnAgreement.classList.toggle('active', next === 'agreement');
      agreementVisible = next === 'agreement';
      setLegend();
      render();
    }
    function setScheme(next) {
      scheme = next;
      btnLeiden.classList.toggle('active', next === 'leiden_cluster');
      btnAgg.classList.toggle('active', next === 'agglomerative');
      btnTreeCluster.classList.toggle('active', next === 'tree_cluster');
      setLegend();
      render();
    }
    function setReferenceBestMode(next) {
      referenceBestMode = next;
      btnRefInclude.classList.toggle('active', next === 'include');
      btnRefExclude.classList.toggle('active', next === 'exclude');
      setLegend();
      render();
    }
    function setScaleMode(next) {
      scaleMode = next;
      btnScaleFull.classList.toggle('active', next === 'full');
      btnScaleRobust.classList.toggle('active', next === 'robust');
      setLegend();
      render();
    }
    btnPcoa.addEventListener('click', () => setMode('pcoa'));
    btnUmap.addEventListener('click', () => setMode('umap'));
    btnAfUmap.addEventListener('click', () => setMode('afumap'));
    btnBest.addEventListener('click', () => setMode('best2'));
    btnTree.addEventListener('click', () => setMode('tree'));
    btnAgreement.addEventListener('click', () => {
      agreementVisible = !agreementVisible;
      btnAgreement.classList.toggle('active', agreementVisible);
      setAgreementPanel();
    });
    btnAgreeAgg.addEventListener('click', () => { agreementTarget = 'agglomerative'; setAgreementPanel(); });
    btnAgreeTree.addEventListener('click', () => { agreementTarget = 'tree'; setAgreementPanel(); });
    agreementContent.addEventListener('click', (event) => {
      const el = event.target.closest('.agreement-cell, .agreement-apply');
      if (!el) return;
      const record = agreementRecord(el.dataset.target, el.dataset.afState, el.dataset.threshold);
      applyAgreementRecord(record);
    });
    btnRefInclude.addEventListener('click', () => setReferenceBestMode('include'));
    btnRefExclude.addEventListener('click', () => setReferenceBestMode('exclude'));
    btnLeiden.addEventListener('click', () => setScheme('leiden_cluster'));
    btnAgg.addEventListener('click', () => setScheme('agglomerative'));
    btnTreeCluster.addEventListener('click', () => setScheme('tree_cluster'));
    btnScaleFull.addEventListener('click', () => setScaleMode('full'));
    btnScaleRobust.addEventListener('click', () => setScaleMode('robust'));
    thresholdSlider.addEventListener('input', () => { setLegend(); render(); });
    treeThresholdSlider.addEventListener('input', () => { setLegend(); render(); });
    novelThresholdInput.addEventListener('input', () => { setLegend(); render(); });
    afKmerSlider.addEventListener('input', () => { selectedAfStateId = null; setLegend(); render(); });
    afNeighborsSlider.addEventListener('input', () => { selectedAfStateId = null; setLegend(); render(); });
    afMinDistSlider.addEventListener('input', () => { selectedAfStateId = null; setLegend(); render(); });
    minLengthInput.addEventListener('input', () => { setLegend(); render(); });
    maxLengthInput.addEventListener('input', () => { setLegend(); render(); });

    function downloadText(filename, text) {
      const blob = new Blob([text], {type: 'text/plain;charset=utf-8'});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = filename;
      document.body.appendChild(a);
      a.click();
      a.remove();
      URL.revokeObjectURL(url);
    }
    function updateReviewExportLabel() {
      const sel = selectedReviewClusters();
      if (sel.length === 0) reviewExportLabel.textContent = 'Select clusters ▾';
      else if (sel.length === 1) reviewExportLabel.textContent = sel[0] + ' ▾';
      else reviewExportLabel.textContent = sel.length + ' clusters selected ▾';
    }
    function selectedReviewClusters() {
      return [...reviewExportItems.querySelectorAll('input[type=checkbox]:checked')].map(cb => cb.value);
    }
    function syncReviewExportOptions() {
      const prior = selectedReviewClusters();
      reviewExportItems.innerHTML = '';
      const clusterSet = new Set();
      payload.points.forEach(p => {
        const c = pointCluster(p);
        if (c != null) clusterSet.add(String(c));
      });
      const clusterNames = Array.from(clusterSet).sort((a, b) => {
        const na = a.match(/\d+/), nb = b.match(/\d+/);
        if (na && nb) return Number(na[0]) - Number(nb[0]);
        return a.localeCompare(b);
      });
      clusterNames.forEach(name => {
        const lbl = document.createElement('label');
        lbl.className = 'multiselect-item';
        const cb = document.createElement('input');
        cb.type = 'checkbox';
        cb.value = name;
        if (prior.includes(name)) cb.checked = true;
        lbl.appendChild(cb);
        const count = payload.points.filter(p => String(pointCluster(p)) === name).length;
        lbl.appendChild(document.createTextNode(name + ' (' + count + ')'));
        reviewExportItems.appendChild(lbl);
      });
      const all = reviewExportItems.querySelectorAll('input[type=checkbox]');
      reviewExportAll.checked = all.length > 0 && [...all].every(cb => cb.checked);
      updateReviewExportLabel();
    }
    function currentReviewClusterMembers(clusterNames) {
      const nameSet = new Set(clusterNames);
      return payload.points
        .filter(p => nameSet.has(String(pointCluster(p))))
        .map(p => p.id);
    }
    function currentReviewExportCommand(clusterNames) {
      const em = payload.export_meta;
      const thresh = currentThreshold();
      const treeThresh = currentTreeThreshold();
      const joinedNames = clusterNames.join('_');
      const outDir = em.default_export_dir;
      if (scheme === 'agglomerative') {
        return [
          'cd seqscape',
          [
            'PYTHONPATH=src python -m seqscape.cli export-clusters',
            "--sequence-fasta '" + em.panel_fasta + "'",
            "--cluster-tsv '" + em.agglomerative_tsv + "'",
            '--cluster-column agglomerative_cluster',
            '--filter-column threshold',
            "--filter-value '" + thresh + "'",
            "--cluster-list '" + clusterNames.join(',') + "'",
            "--filename-prefix 'agg_" + thresh + '_' + joinedNames + "'",
            "--outdir '" + outDir + "'",
          ].join(' '),
        ].join('\n');
      }
      if (scheme === 'tree_cluster') {
        return [
          'cd seqscape',
          [
            'PYTHONPATH=src python -m seqscape.cli export-clusters',
            "--sequence-fasta '" + em.panel_fasta + "'",
            "--cluster-tsv '" + em.tree_tsv + "'",
            '--cluster-column tree_cluster',
            '--filter-column threshold',
            "--filter-value '" + treeThresh + "'",
            "--cluster-list '" + clusterNames.join(',') + "'",
            "--filename-prefix 'tree_" + treeThresh + '_' + joinedNames + "'",
            "--outdir '" + outDir + "'",
          ].join(' '),
        ].join('\n');
      }
      return null;
    }
    reviewExportToggle.addEventListener('click', () => {
      syncReviewExportOptions();
      reviewExportMenu.classList.toggle('open');
    });
    document.addEventListener('click', evt => {
      if (!reviewExportToggle.contains(evt.target) && !reviewExportMenu.contains(evt.target)) {
        reviewExportMenu.classList.remove('open');
      }
    });
    reviewExportAll.addEventListener('change', () => {
      const checked = reviewExportAll.checked;
      reviewExportItems.querySelectorAll('input[type=checkbox]').forEach(cb => { cb.checked = checked; });
      updateReviewExportLabel();
    });
    reviewExportItems.addEventListener('change', () => {
      const all = reviewExportItems.querySelectorAll('input[type=checkbox]');
      reviewExportAll.checked = [...all].every(cb => cb.checked);
      updateReviewExportLabel();
    });
    reviewExportIdsBtn.addEventListener('click', () => {
      const names = selectedReviewClusters();
      if (names.length === 0) { alert('Select at least one cluster.'); return; }
      const ids = currentReviewClusterMembers(names);
      const tag = scheme === 'tree_cluster' ? ('tree_' + currentTreeThreshold()) : (scheme === 'agglomerative' ? ('agg_' + currentThreshold()) : 'leiden');
      downloadText(tag + '_' + names.join('_') + '_ids.txt', ids.join('\n') + '\n');
    });
    reviewExportCmdBtn.addEventListener('click', () => {
      const names = selectedReviewClusters();
      if (names.length === 0) { alert('Select at least one cluster.'); return; }
      const cmd = currentReviewExportCommand(names);
      if (!cmd) { alert('Command export is only available for Agglomerative and Tree clustering schemes.'); return; }
      const tag = scheme === 'tree_cluster' ? ('tree_' + currentTreeThreshold()) : ('agg_' + currentThreshold());
      downloadText(tag + '_' + names.join('_') + '_export.sh', cmd + '\n');
    });

    setLegend();
    render();
  </script>
</body>
</html>
"""
    return template.replace("__PAYLOAD__", payload_json).replace("__TITLE__", title)


def run(args) -> None:
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    ids, ident = read_square_matrix(Path(args.identity_matrix).resolve())
    _, dist = read_square_matrix(Path(args.distance_matrix).resolve())
    manifest_rows = read_manifest(Path(args.manifest_tsv).resolve())
    manifest = {row["id"]: row for row in manifest_rows}
    if set(ids) != set(manifest.keys()):
        raise RuntimeError("Manifest ids do not match matrix ids")
    ordered_manifest = [manifest[seq_id] for seq_id in ids]

    pcoa_coords_raw, pcoa_eigvals = pcoa(dist)
    pcoa_coords = np.zeros((len(ids), 2), dtype=float)
    if pcoa_coords_raw.shape[1] >= 1:
        pcoa_coords[:, 0] = pcoa_coords_raw[:, 0]
    if pcoa_coords_raw.shape[1] >= 2:
        pcoa_coords[:, 1] = pcoa_coords_raw[:, 1]
    pcoa_pct1, pcoa_pct2 = pcoa_variance_explained(pcoa_eigvals)
    umap_coords = run_precomputed_umap(
        dist,
        neighbors=args.umap_neighbors,
        min_dist=args.umap_min_dist,
        spread=args.umap_spread,
        seed=args.umap_seed,
    )

    with (outdir / "pcoa_coords.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["ID", "PCoA1", "PCoA2"])
        for seq_id, (x, y) in zip(ids, pcoa_coords):
            writer.writerow([seq_id, f"{float(x):.6f}", f"{float(y):.6f}"])
    with (outdir / "distance_umap_coords.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["ID", "UMAP1", "UMAP2"])
        for seq_id, (x, y) in zip(ids, umap_coords):
            writer.writerow([seq_id, f"{float(x):.6f}", f"{float(y):.6f}"])

    ref_ids = [row["id"] for row in ordered_manifest if row["item_class"] == "reference"]
    ref_idx = [ids.index(rid) for rid in ref_ids]
    best_rows: List[dict] = []
    points: List[dict] = []
    for i, row in enumerate(ordered_manifest):
        values = ident[i, ref_idx]
        if row["item_class"] == "reference":
            local_order = [j for j in np.argsort(values)[::-1].tolist() if ref_ids[j] != row["id"]]
            if not local_order:
                best_pos = second_pos = 0
            elif len(local_order) == 1:
                best_pos = second_pos = local_order[0]
            else:
                best_pos, second_pos = local_order[0], local_order[1]
        else:
            local_order = np.argsort(values)[::-1].tolist()
            best_pos = local_order[0]
            second_pos = local_order[1] if len(local_order) > 1 else best_pos
        best_ref_id = ref_ids[best_pos]
        second_ref_id = ref_ids[second_pos]
        best_identity = float(values[best_pos])
        second_identity = float(values[second_pos])
        best_rows.append(
            {
                "id": row["id"],
                "label": row["label"],
                "item_class": row["item_class"],
                "source_cluster": row["source_cluster"],
                "best_ref_id": best_ref_id,
                "best_ref_label": manifest[best_ref_id]["label"],
                "best_identity": f"{best_identity:.6f}",
                "second_ref_id": second_ref_id,
                "second_ref_label": manifest[second_ref_id]["label"],
                "second_identity": f"{second_identity:.6f}",
                "candidate_novel": "1" if row["item_class"] == "genome" and best_identity < args.novel_threshold else "0",
            }
        )
        points.append(
            {
                "id": row["id"],
                "label": row["label"],
                "item_class": row["item_class"],
                "length_bp": int(row.get("length_bp") or 0),
                "leiden_cluster": row["source_cluster"],
                "best_ref": manifest[best_ref_id]["label"],
                "best_identity": best_identity,
                "best2": [best_identity, second_identity],
                "best2_ref_include": [100.0, best_identity] if row["item_class"] == "reference" else [best_identity, second_identity],
                "best2_ref_exclude": [best_identity, second_identity],
                "pcoa": [float(pcoa_coords[i, 0]), float(pcoa_coords[i, 1])],
                "umap": [float(umap_coords[i, 0]), float(umap_coords[i, 1])],
                "agglomerative": {},
                "tree_clusters": {},
            }
        )
    write_tsv(
        outdir / "best_reference_summary.tsv",
        best_rows,
        [
            "id",
            "label",
            "item_class",
            "source_cluster",
            "best_ref_id",
            "best_ref_label",
            "best_identity",
            "second_ref_id",
            "second_ref_label",
            "second_identity",
            "candidate_novel",
        ],
    )

    thresholds = [float(x.strip()) for x in args.thresholds.split(",") if x.strip()]
    if not thresholds:
        raise RuntimeError("At least one threshold is required")

    leiden_order = sorted({row["source_cluster"] for row in ordered_manifest if row["item_class"] == "genome"}, key=cluster_sort_key)
    leiden_colors = {cid: PALETTE[i % len(PALETTE)] for i, cid in enumerate(leiden_order)}
    cluster_views = {
        "leiden_cluster": {
            "pcoa": build_cluster_view(points, "pcoa", "leiden_cluster", leiden_order, leiden_colors),
            "umap": build_cluster_view(points, "umap", "leiden_cluster", leiden_order, leiden_colors),
            "best2": build_cluster_view(points, "best2", "leiden_cluster", leiden_order, leiden_colors),
        },
        "agglomerative": {},
        "tree": {},
    }
    threshold_rows = []
    long_agglom_rows = []
    threshold_validation_rows = []
    genome_leiden = {row["id"]: row["source_cluster"] for row in ordered_manifest if row["item_class"] == "genome"}
    agglomerative_genome_clusters: dict[str, dict[str, str]] = {}
    for threshold in thresholds:
        raw = fit_agglomerative(dist, "complete", threshold)
        cluster_names = normalize_cluster_labels(raw, prefix="A")
        agglom_order = sorted(set(cluster_names), key=lambda cid: cid)
        genome_agg = {}
        for point, cid in zip(points, cluster_names):
            point["agglomerative"][str(threshold)] = cid
            if point["item_class"] == "genome":
                genome_agg[point["id"]] = cid
        agglomerative_genome_clusters[f"{threshold:.6f}"] = genome_agg
        genome_counts = Counter(cid for row, cid in zip(ordered_manifest, cluster_names) if row["item_class"] == "genome")
        common_ids = sorted(set(genome_leiden) & set(genome_agg))
        nmi_value = nmi([genome_leiden[i] for i in common_ids], [genome_agg[i] for i in common_ids])
        threshold_rows.append(
            {
                "threshold": f"{threshold:.6f}",
                "cluster_count": len(agglom_order),
                "largest_genome_cluster": max(genome_counts.values()) if genome_counts else 0,
            }
        )
        threshold_validation_rows.append(
            {
                "threshold": f"{threshold:.6f}",
                "cluster_count": len(agglom_order),
                "nmi": f"{nmi_value:.6f}",
            }
        )
        for seq_id, cid in zip(ids, cluster_names):
            long_agglom_rows.append({"id": seq_id, "threshold": f"{threshold:.6f}", "agglomerative_cluster": cid})

    for threshold in thresholds:
        threshold_key = str(threshold)
        for point in points:
            point["current_agglomerative"] = point["agglomerative"][threshold_key]
        order = sorted({point["current_agglomerative"] for point in points})
        colors = {cid: PALETTE[i % len(PALETTE)] for i, cid in enumerate(order)}
        cluster_views["agglomerative"][threshold_key] = {
            "pcoa": build_cluster_view(points, "pcoa", "current_agglomerative", order, colors),
            "umap": build_cluster_view(points, "umap", "current_agglomerative", order, colors),
            "best2": build_cluster_view(points, "best2", "current_agglomerative", order, colors),
        }

    write_tsv(outdir / "agglomerative_threshold_summary.tsv", threshold_rows, ["threshold", "cluster_count", "largest_genome_cluster"])
    write_tsv(outdir / "agglomerative_threshold_assignments.tsv", long_agglom_rows, ["id", "threshold", "agglomerative_cluster"])
    write_tsv(outdir / "agglomerative_validation.tsv", threshold_validation_rows, ["threshold", "cluster_count", "nmi"])
    build_heatmaps(ids, ident, ordered_manifest, outdir)
    tree_payload, tree_dist = build_neighbor_joining_tree(ids, dist, ordered_manifest, outdir)
    tree_threshold_rows = []
    tree_long_rows = []
    tree_validation_rows = []
    tree_genome_clusters: dict[str, dict[str, str]] = {}
    for threshold in thresholds:
        raw = fit_agglomerative(tree_dist, "complete", threshold)
        cluster_names = normalize_cluster_labels(raw, prefix="T")
        tree_order = sorted(set(cluster_names), key=cluster_sort_key)
        genome_tree = {}
        for point, cid in zip(points, cluster_names):
            point["tree_clusters"][str(threshold)] = cid
            if point["item_class"] == "genome":
                genome_tree[point["id"]] = cid
        tree_genome_clusters[f"{threshold:.6f}"] = genome_tree
        genome_counts = Counter(cid for row, cid in zip(ordered_manifest, cluster_names) if row["item_class"] == "genome")
        common_ids = sorted(set(genome_leiden) & set(genome_tree))
        nmi_value = nmi([genome_leiden[i] for i in common_ids], [genome_tree[i] for i in common_ids])
        tree_threshold_rows.append(
            {
                "threshold": f"{threshold:.6f}",
                "cluster_count": len(tree_order),
                "largest_genome_cluster": max(genome_counts.values()) if genome_counts else 0,
            }
        )
        tree_validation_rows.append(
            {
                "threshold": f"{threshold:.6f}",
                "cluster_count": len(tree_order),
                "nmi": f"{nmi_value:.6f}",
            }
        )
        for seq_id, cid in zip(ids, cluster_names):
            tree_long_rows.append({"id": seq_id, "threshold": f"{threshold:.6f}", "tree_cluster": cid})
        for point in points:
            point["current_tree_cluster"] = point["tree_clusters"][str(threshold)]
        colors = {cid: PALETTE[i % len(PALETTE)] for i, cid in enumerate(tree_order)}
        cluster_views["tree"][str(threshold)] = {
            "pcoa": build_cluster_view(points, "pcoa", "current_tree_cluster", tree_order, colors),
            "umap": build_cluster_view(points, "umap", "current_tree_cluster", tree_order, colors),
            "best2": build_cluster_view(points, "best2", "current_tree_cluster", tree_order, colors),
        }
    write_tsv(outdir / "tree_threshold_summary.tsv", tree_threshold_rows, ["threshold", "cluster_count", "largest_genome_cluster"])
    write_tsv(outdir / "tree_threshold_assignments.tsv", tree_long_rows, ["id", "threshold", "tree_cluster"])
    write_tsv(outdir / "tree_validation.tsv", tree_validation_rows, ["threshold", "cluster_count", "nmi"])

    chosen_threshold_row = choose_closest_cluster_threshold(threshold_rows, len(leiden_order))
    threshold_metrics = {
        row["threshold"]: {
            "cluster_count": int(row["cluster_count"]),
            "nmi": float(vrow["nmi"]),
        }
        for row, vrow in zip(threshold_rows, threshold_validation_rows)
    }
    tree_threshold_metrics = {
        row["threshold"]: {
            "cluster_count": int(row["cluster_count"]),
            "nmi": float(vrow["nmi"]),
        }
        for row, vrow in zip(tree_threshold_rows, tree_validation_rows)
    }
    chosen_tree_threshold_row = choose_closest_cluster_threshold(tree_threshold_rows, len(leiden_order))
    novel_genome_count = sum(1 for row in best_rows if row["item_class"] == "genome" and row["candidate_novel"] == "1")
    sample_validation = read_optional_json(Path(args.manifest_tsv).resolve().parent / "validation_metrics.json")
    align_validation = read_optional_json(Path(args.distance_matrix).resolve().parent / "validation_metrics.json")
    sample_summary = read_key_value_summary(Path(args.manifest_tsv).resolve().parent / "summary.txt")
    af_widget_dir_text = getattr(args, "af_widget_dir", "") or sample_summary.get("af_widget_dir", "")
    chosen_state_id = sample_summary.get("chosen_state", "")
    af_states, default_af_state_index = ([], 0)
    genome_ids = {row["id"] for row in ordered_manifest if row["item_class"] == "genome"}
    if af_widget_dir_text:
        af_states, default_af_state_index = read_af_states(Path(af_widget_dir_text).resolve(), set(ids), genome_ids, chosen_state_id)
    chosen_tree_threshold = float(chosen_tree_threshold_row["threshold"])
    default_tree_threshold_index = min(range(len(thresholds)), key=lambda idx: abs(thresholds[idx] - chosen_tree_threshold))
    cluster_agreement_rows = []
    threshold_count_by_key = {row["threshold"]: int(row["cluster_count"]) for row in threshold_rows}
    tree_threshold_count_by_key = {row["threshold"]: int(row["cluster_count"]) for row in tree_threshold_rows}
    for state in af_states:
        af_genome_clusters = {
            seq_id: cid
            for seq_id, cid in state["clusters"].items()
            if seq_id in genome_ids
        }
        for threshold_key, target_clusters in agglomerative_genome_clusters.items():
            common_ids = sorted(set(af_genome_clusters) & set(target_clusters))
            nmi_value = nmi([af_genome_clusters[i] for i in common_ids], [target_clusters[i] for i in common_ids])
            cluster_agreement_rows.append(
                {
                    "target": "agglomerative",
                    "af_state": state["state_id"],
                    "kmer": state["kmer"],
                    "neighbors": state["neighbors"],
                    "min_dist": state["min_dist"],
                    "leiden_resolution": state["leiden_resolution"],
                    "af_cluster_count": state["cluster_count"],
                    "threshold": threshold_key,
                    "target_cluster_count": threshold_count_by_key[threshold_key],
                    "nmi": float(nmi_value),
                }
            )
        for threshold_key, target_clusters in tree_genome_clusters.items():
            common_ids = sorted(set(af_genome_clusters) & set(target_clusters))
            nmi_value = nmi([af_genome_clusters[i] for i in common_ids], [target_clusters[i] for i in common_ids])
            cluster_agreement_rows.append(
                {
                    "target": "tree",
                    "af_state": state["state_id"],
                    "kmer": state["kmer"],
                    "neighbors": state["neighbors"],
                    "min_dist": state["min_dist"],
                    "leiden_resolution": state["leiden_resolution"],
                    "af_cluster_count": state["cluster_count"],
                    "threshold": threshold_key,
                    "target_cluster_count": tree_threshold_count_by_key[threshold_key],
                    "nmi": float(nmi_value),
                }
            )
    if cluster_agreement_rows:
        write_tsv(
            outdir / "cluster_agreement_matrix.tsv",
            [
                {
                    **row,
                    "nmi": f"{row['nmi']:.6f}",
                }
                for row in cluster_agreement_rows
            ],
            [
                "target",
                "af_state",
                "kmer",
                "neighbors",
                "min_dist",
                "leiden_resolution",
                "af_cluster_count",
                "threshold",
                "target_cluster_count",
                "nmi",
            ],
        )
    review_validation = {
        "closest_threshold_to_leiden": chosen_threshold_row["threshold"],
        "closest_threshold_cluster_count": int(chosen_threshold_row["cluster_count"]),
        "closest_tree_threshold_to_leiden": chosen_tree_threshold_row["threshold"],
        "closest_tree_threshold_cluster_count": int(chosen_tree_threshold_row["cluster_count"]),
        "novel_genome_count": novel_genome_count,
        "novel_threshold_pct": args.novel_threshold,
        "pcoa_axis1_variance_explained_pct": pcoa_pct1,
        "pcoa_axis2_variance_explained_pct": pcoa_pct2,
    }
    validation_payload = {
        "sample": sample_validation,
        "align": align_validation,
        "review": review_validation,
        "thresholds": threshold_metrics,
        "tree_thresholds": tree_threshold_metrics,
    }
    (outdir / "validation_metrics.json").write_text(json.dumps(validation_payload, indent=2))

    payload = {
        "thresholds": [str(x) for x in thresholds],
        "export_meta": {
            "panel_fasta": str(Path(args.panel_fasta).resolve()),
            "agglomerative_tsv": str(outdir / "agglomerative_threshold_assignments.tsv"),
            "tree_tsv": str(outdir / "tree_threshold_assignments.tsv"),
            "default_export_dir": str(outdir / "exports"),
        },
        "default_threshold_index": min(len(thresholds) - 1, thresholds.index(0.06) if 0.06 in thresholds else 0),
        "novel_threshold": args.novel_threshold,
        "pcoa_axis_labels": {
            "x": f"PCoA 1 ({pcoa_pct1:.1f}%)",
            "y": f"PCoA 2 ({pcoa_pct2:.1f}%)",
        },
        "points": points,
        "cluster_views": cluster_views,
        "validation": validation_payload,
        "tree": tree_payload,
        "af_states": af_states,
        "default_af_state_index": default_af_state_index,
        "default_tree_threshold_index": default_tree_threshold_index,
        "cluster_agreement": cluster_agreement_rows,
    }
    (outdir / "panel_review_widget.html").write_text(
        html_template(json.dumps(payload, separators=(",", ":")), "Aligned Panel Review"),
        encoding="utf-8",
    )
    with (outdir / "summary.txt").open("w") as fh:
        fh.write(f"panel_fasta\t{Path(args.panel_fasta).resolve()}\n")
        fh.write(f"identity_matrix\t{Path(args.identity_matrix).resolve()}\n")
        fh.write(f"distance_matrix\t{Path(args.distance_matrix).resolve()}\n")
        fh.write(f"manifest_tsv\t{Path(args.manifest_tsv).resolve()}\n")
        if af_widget_dir_text:
            fh.write(f"af_widget_dir\t{Path(af_widget_dir_text).resolve()}\n")
            fh.write(f"af_states_loaded\t{len(af_states)}\n")
        fh.write(f"thresholds\t{','.join(str(x) for x in thresholds)}\n")
        fh.write(f"novel_threshold\t{args.novel_threshold}\n")
        fh.write(f"novel_genome_count\t{novel_genome_count}\n")
        fh.write(f"closest_threshold_to_leiden\t{chosen_threshold_row['threshold']}\n")
        fh.write(f"closest_threshold_cluster_count\t{chosen_threshold_row['cluster_count']}\n")
        fh.write(f"closest_tree_threshold_to_leiden\t{chosen_tree_threshold_row['threshold']}\n")
        fh.write(f"closest_tree_threshold_cluster_count\t{chosen_tree_threshold_row['cluster_count']}\n")
        fh.write(f"pcoa_axis1_variance_explained_pct\t{pcoa_pct1:.6f}\n")
        fh.write(f"pcoa_axis2_variance_explained_pct\t{pcoa_pct2:.6f}\n")
        fh.write(f"umap_neighbors\t{args.umap_neighbors}\n")
        fh.write(f"umap_min_dist\t{args.umap_min_dist}\n")
        fh.write(f"umap_spread\t{args.umap_spread}\n")
