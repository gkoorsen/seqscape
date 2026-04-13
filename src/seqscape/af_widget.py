from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import igraph as ig
import leidenalg
import numpy as np
from sklearn.neighbors import NearestNeighbors

from .io_utils import load_references, write_tsv
from .kmer import build_kmer_matrix, load_fasta
from .labels import load_cg_map
from .ordination import run_euclidean_umap
from .palettes import PALETTE


def status(msg: str) -> None:
    print(msg, flush=True)


def parse_ints(text: str) -> List[int]:
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def parse_floats(text: str) -> List[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def state_slug(kmer: int, neighbors: int, min_dist: float) -> str:
    return f"k{kmer}_n{neighbors}_d{str(min_dist).replace('.', 'p')}"


def build_leiden_graph(genome_coords: np.ndarray, neighbor_k: int) -> Tuple[ig.Graph, List[float]]:
    nn = NearestNeighbors(n_neighbors=min(neighbor_k + 1, len(genome_coords)), metric="euclidean")
    nn.fit(genome_coords)
    dists, nbrs = nn.kneighbors(genome_coords)
    positive = dists[:, 1:]
    sigma = float(np.median(positive[positive > 0])) if np.any(positive > 0) else 1.0
    edges: Dict[Tuple[int, int], float] = {}
    for i in range(len(genome_coords)):
        for dist, j in zip(dists[i, 1:], nbrs[i, 1:]):
            a, b = (i, int(j)) if i < int(j) else (int(j), i)
            if a == b:
                continue
            w = math.exp(-((float(dist) ** 2) / (2.0 * sigma * sigma))) if sigma > 0 else 1.0
            prev = edges.get((a, b))
            if prev is None or w > prev:
                edges[(a, b)] = w
    g = ig.Graph()
    g.add_vertices(len(genome_coords))
    g.add_edges(list(edges.keys()))
    return g, list(edges.values())


def normalize_cluster_labels(labels: Sequence[int], prefix: str = "C") -> List[str]:
    mapping: Dict[int, str] = {}
    out: List[str] = []
    next_idx = 1
    for raw in labels:
        if raw not in mapping:
            mapping[raw] = f"{prefix}{next_idx:02d}"
            next_idx += 1
        out.append(mapping[raw])
    return out


def assign_reference_clusters(
    genome_coords: np.ndarray,
    genome_clusters: List[str],
    ref_coords: np.ndarray,
) -> List[str]:
    nn = NearestNeighbors(n_neighbors=min(25, len(genome_coords)), metric="euclidean")
    nn.fit(genome_coords)
    _, idxs = nn.kneighbors(ref_coords)
    out: List[str] = []
    for neigh in idxs:
        votes = [genome_clusters[i] for i in neigh]
        out.append(Counter(votes).most_common(1)[0][0])
    return out


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


def html_template(payload_json: str, title: str) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{title}</title>
  <style>
    body {{
      margin: 0;
      font-family: Helvetica, Arial, sans-serif;
      background: #f3f1eb;
      color: #1d1d1d;
    }}
    .wrap {{
      max-width: 1500px;
      margin: 0 auto;
      padding: 18px 20px 28px;
    }}
    h1 {{
      margin: 0 0 8px;
      font-size: 28px;
      line-height: 1.1;
    }}
    .sub {{
      margin: 0 0 16px;
      color: #555;
      font-size: 14px;
    }}
    .controls {{
      display: grid;
      grid-template-columns: 1fr 1fr 1fr;
      gap: 14px 20px;
      margin-bottom: 16px;
      padding: 14px 16px;
      background: white;
      border: 1px solid #d8d2c6;
      border-radius: 12px;
    }}
    .control label {{
      display: block;
      font-size: 13px;
      font-weight: 700;
      margin-bottom: 6px;
    }}
    .control input {{
      width: 100%;
    }}
    .value {{
      font-weight: 700;
      color: #7d3f00;
    }}
    .legend {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin: 10px 0 14px;
      font-size: 13px;
    }}
    .legend-item {{
      display: flex;
      align-items: center;
      gap: 6px;
      background: white;
      border: 1px solid #ddd6ca;
      border-radius: 999px;
      padding: 5px 10px;
    }}
    .swatch {{
      width: 12px;
      height: 12px;
      border-radius: 999px;
      display: inline-block;
      border: 1px solid rgba(0,0,0,0.12);
    }}
    canvas {{
      width: 100%;
      height: auto;
      background: #fffdf8;
      border: 1px solid #d8d2c6;
      border-radius: 16px;
      box-shadow: 0 8px 24px rgba(0,0,0,0.06);
    }}
    .footer {{
      margin-top: 10px;
      color: #666;
      font-size: 12px;
    }}
    .tooltip {{
      position: fixed;
      display: none;
      pointer-events: none;
      z-index: 10;
      background: rgba(255,255,255,0.97);
      color: #111;
      border: 1px solid #444;
      border-radius: 10px;
      box-shadow: 0 10px 24px rgba(0,0,0,0.14);
      padding: 8px 10px;
      font-size: 12px;
      line-height: 1.35;
      white-space: nowrap;
    }}
    .tooltip .head {{
      font-weight: 700;
      margin-bottom: 4px;
    }}
  </style>
</head>
<body>
  <div class="wrap">
    <h1>{title}</h1>
    <p class="sub">Alignment-free UMAP over genomes plus the reference panel. Each slider state is independently clustered with Leiden, and the per-state assignments are written to disk for downstream representative sampling.</p>
    <div class="controls">
      <div class="control">
        <label for="kmer">`kmer`: <span id="kmer-value" class="value"></span></label>
        <input id="kmer" type="range" min="0" max="0" step="1">
      </div>
      <div class="control">
        <label for="neighbors">`n_neighbors`: <span id="neighbors-value" class="value"></span></label>
        <input id="neighbors" type="range" min="0" max="0" step="1">
      </div>
      <div class="control">
        <label for="min-dist">`min_dist`: <span id="min-dist-value" class="value"></span></label>
        <input id="min-dist" type="range" min="0" max="0" step="1">
      </div>
    </div>
    <div class="legend" id="legend"></div>
    <canvas id="plot" width="1400" height="980"></canvas>
    <div class="tooltip" id="tooltip"></div>
    <div class="footer" id="footer"></div>
  </div>
  <script>
    const payload = {payload_json};
    const plot = document.getElementById('plot');
    const ctx = plot.getContext('2d');
    const tooltip = document.getElementById('tooltip');
    const legend = document.getElementById('legend');
    const footer = document.getElementById('footer');
    const kmerSlider = document.getElementById('kmer');
    const neighborsSlider = document.getElementById('neighbors');
    const minDistSlider = document.getElementById('min-dist');
    const kmerValue = document.getElementById('kmer-value');
    const neighborsValue = document.getElementById('neighbors-value');
    const minDistValue = document.getElementById('min-dist-value');
    const W = plot.width, H = plot.height;
    const M = {{ left: 72, right: 34, top: 34, bottom: 62 }};

    const kmers = payload.kmers;
    const neighbors = payload.neighbors;
    const minDists = payload.min_dists;
    const points = payload.points;
    const states = payload.states;
    let kmerIndex = payload.default_kmer_index;
    let neighborIndex = payload.default_neighbor_index;
    let minDistIndex = payload.default_min_dist_index;
    let pointHitboxes = [];

    kmerSlider.max = String(kmers.length - 1);
    kmerSlider.value = String(kmerIndex);
    neighborsSlider.max = String(neighbors.length - 1);
    neighborsSlider.value = String(neighborIndex);
    minDistSlider.max = String(minDists.length - 1);
    minDistSlider.value = String(minDistIndex);

    function key() {{
      return `${{kmers[kmerIndex]}}|${{neighbors[neighborIndex]}}|${{minDists[minDistIndex]}}`;
    }}

    function currentState() {{
      return states[key()];
    }}

    function drawFilledStar(ctx2, x, y, fill, stroke) {{
      ctx2.save();
      ctx2.translate(x, y);
      ctx2.fillStyle = fill;
      ctx2.strokeStyle = stroke;
      ctx2.lineWidth = 1.1;
      ctx2.beginPath();
      const spikes = 5, outerR = 7, innerR = 3.2;
      for (let i = 0; i < spikes * 2; i++) {{
        const r = (i % 2 === 0) ? outerR : innerR;
        const ang = -Math.PI / 2 + (i * Math.PI / spikes);
        const px = Math.cos(ang) * r;
        const py = Math.sin(ang) * r;
        if (i === 0) ctx2.moveTo(px, py); else ctx2.lineTo(px, py);
      }}
      ctx2.closePath();
      ctx2.fill();
      ctx2.stroke();
      ctx2.restore();
    }}

    function drawLegend(state) {{
      legend.innerHTML = '';
      for (const cl of state.clusters) {{
        const item = document.createElement('div');
        item.className = 'legend-item';
        item.innerHTML = `<span class="swatch" style="background:${{cl.color}}"></span>${{cl.name}} (${{cl.genomes}}, ${{cl.pct.toFixed(1)}}%)`;
        legend.appendChild(item);
      }}
      const ref = document.createElement('div');
      ref.className = 'legend-item';
      ref.innerHTML = `<span style="font-size:14px;color:#1F77B4;">★</span> references (${{payload.reference_count}})`;
      legend.appendChild(ref);
    }}

    function showTooltip(evt, lines) {{
      tooltip.innerHTML = `<div class="head">${{lines[0]}}</div>${{lines.slice(1).map(t => `<div>${{t}}</div>`).join('')}}`;
      tooltip.style.display = 'block';
      tooltip.style.left = `${{Math.min(window.innerWidth - 260, evt.clientX + 14)}}px`;
      tooltip.style.top = `${{Math.max(10, evt.clientY - 12)}}px`;
    }}
    function hideTooltip() {{ tooltip.style.display = 'none'; }}

    function draw() {{
      const state = currentState();
      const coords = state.coords;
      kmerValue.textContent = kmers[kmerIndex];
      neighborsValue.textContent = neighbors[neighborIndex];
      minDistValue.textContent = minDists[minDistIndex];
      footer.textContent = `${{payload.genome_count}} genomes + ${{payload.reference_count}} references | Leiden clusters: ${{state.cluster_count}} | state: ${{state.state_slug}}`;

      const xs = coords.map(v => v[0]);
      const ys = coords.map(v => v[1]);
      const xmin = Math.min(...xs), xmax = Math.max(...xs);
      const ymin = Math.min(...ys), ymax = Math.max(...ys);
      const xspan = (xmax - xmin) || 1;
      const yspan = (ymax - ymin) || 1;

      function sx(x) {{ return M.left + ((x - xmin) / xspan) * (W - M.left - M.right); }}
      function sy(y) {{ return H - M.bottom - ((y - ymin) / yspan) * (H - M.top - M.bottom); }}

      ctx.clearRect(0, 0, W, H);
      ctx.fillStyle = "#fffdf8";
      ctx.fillRect(0, 0, W, H);
      ctx.strokeStyle = "#d9d1c4";
      ctx.lineWidth = 1;
      ctx.strokeRect(M.left, M.top, W - M.left - M.right, H - M.top - M.bottom);

      for (let i = 0; i <= 5; i++) {{
        const gx = M.left + (i / 5) * (W - M.left - M.right);
        const gy = M.top + (i / 5) * (H - M.top - M.bottom);
        ctx.strokeStyle = "rgba(0,0,0,0.05)";
        ctx.beginPath(); ctx.moveTo(gx, M.top); ctx.lineTo(gx, H - M.bottom); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(M.left, gy); ctx.lineTo(W - M.right, gy); ctx.stroke();
      }}

      for (const hull of state.hulls) {{
        const color = hull.color;
        ctx.beginPath();
        hull.points.forEach((p, idx) => {{
          const x = sx(p[0]), y = sy(p[1]);
          if (idx === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        }});
        ctx.closePath();
        ctx.fillStyle = color + "14";
        ctx.strokeStyle = color + "99";
        ctx.lineWidth = 1;
        ctx.fill();
        ctx.stroke();
      }}

      pointHitboxes = [];
      for (let i = 0; i < points.length; i++) {{
        const p = points[i];
        const c = coords[i];
        const x = sx(c[0]), y = sy(c[1]);
        const color = state.point_colors[i];
        if (p.item_class === "reference") {{
          drawFilledStar(ctx, x, y, color, "#111");
          pointHitboxes.push({{ x, y, r: 10, point: p, cluster: state.point_cluster_names[i] }});
        }} else {{
          ctx.beginPath();
          ctx.fillStyle = color;
          ctx.globalAlpha = 0.76;
          ctx.arc(x, y, 2.8, 0, Math.PI * 2);
          ctx.fill();
          ctx.globalAlpha = 1.0;
          pointHitboxes.push({{ x, y, r: 5.5, point: p, cluster: state.point_cluster_names[i] }});
        }}
      }}

      ctx.fillStyle = "#222";
      ctx.font = "12px Helvetica";
      for (const label of state.labels) {{
        const x = sx(label.x) + label.offset_x;
        const y = sy(label.y) + label.offset_y;
        ctx.fillText(`${{label.cluster}} (${{label.genomes}}, ${{label.pct.toFixed(1)}}%)`, x, y);
      }}

      drawLegend(state);
    }}

    function syncFromSliders() {{
      kmerIndex = Number(kmerSlider.value);
      neighborIndex = Number(neighborsSlider.value);
      minDistIndex = Number(minDistSlider.value);
      draw();
    }}

    plot.addEventListener('mousemove', evt => {{
      const rect = plot.getBoundingClientRect();
      const mx = (evt.clientX - rect.left) * (plot.width / rect.width);
      const my = (evt.clientY - rect.top) * (plot.height / rect.height);
      let hovered = null;
      for (const hb of pointHitboxes) {{
        const dx = mx - hb.x;
        const dy = my - hb.y;
        if ((dx * dx + dy * dy) <= hb.r * hb.r) {{
          hovered = hb;
          break;
        }}
      }}
      if (!hovered) {{
        hideTooltip();
        return;
      }}
      showTooltip(evt, [
        hovered.point.label,
        `id: ${{hovered.point.id}}`,
        `type: ${{hovered.point.item_class}}`,
        `length: ${{hovered.point.length_bp}} bp`,
        `cluster: ${{hovered.cluster}}`,
      ]);
    }});

    plot.addEventListener('mouseleave', () => hideTooltip());
    kmerSlider.addEventListener('input', syncFromSliders);
    neighborsSlider.addEventListener('input', syncFromSliders);
    minDistSlider.addEventListener('input', syncFromSliders);
    draw();
  </script>
</body>
</html>
"""


def run(args) -> None:
    outdir = Path(args.outdir).resolve()
    states_dir = outdir / "states"
    states_dir.mkdir(parents=True, exist_ok=True)

    kmers = parse_ints(args.kmer_values)
    neighbors_values = parse_ints(args.neighbors_values)
    min_dist_values = parse_floats(args.min_dist_values)
    if not kmers or not neighbors_values or not min_dist_values:
        raise RuntimeError("At least one kmer, neighbors, and min_dist value is required")

    genome_ids, genome_seqs = load_fasta(Path(args.input_fasta).resolve())
    cg_map = load_cg_map(Path(args.cg_map_tsv).resolve()) if args.cg_map_tsv else {}
    refs = load_references(Path(args.reference_fasta).resolve())
    ref_ids = [ref["ref_id"] for ref in refs]
    ref_seqs = [ref["sequence"] for ref in refs]
    ref_labels = {ref["ref_id"]: ref["ref_label"] for ref in refs}

    all_ids = genome_ids + ref_ids
    all_seqs = genome_seqs + ref_seqs
    all_points = []
    for seq_id in genome_ids:
        all_points.append(
            {
                "id": seq_id,
                "label": cg_map.get(seq_id, seq_id),
                "item_class": "genome",
                "length_bp": len(next(seq for gid, seq in zip(genome_ids, genome_seqs) if gid == seq_id)),
            }
        )
    for seq_id, seq in zip(ref_ids, ref_seqs):
        all_points.append(
            {
                "id": seq_id,
                "label": ref_labels.get(seq_id, seq_id),
                "item_class": "reference",
                "length_bp": len(seq),
            }
        )

    states_payload: Dict[str, dict] = {}
    kmer_cache: Dict[int, np.ndarray] = {}

    for kmer in kmers:
        status(f"build_kmer_matrix\t{kmer}")
        matrix = kmer_cache.get(kmer)
        if matrix is None:
            matrix = build_kmer_matrix(all_seqs, kmer)
            kmer_cache[kmer] = matrix
        for neighbors in neighbors_values:
            for min_dist in min_dist_values:
                slug = state_slug(kmer, neighbors, min_dist)
                status(f"state\t{slug}")
                coords = run_euclidean_umap(matrix, neighbors, min_dist, args.seed)
                genome_coords = coords[: len(genome_ids), :]
                ref_coords = coords[len(genome_ids) :, :]

                graph, weights = build_leiden_graph(
                    genome_coords, min(args.leiden_neighbor_k, max(15, len(genome_coords) - 1))
                )
                partition = leidenalg.find_partition(
                    graph,
                    leidenalg.RBConfigurationVertexPartition,
                    weights=weights,
                    resolution_parameter=args.leiden_resolution,
                    seed=args.seed,
                )
                genome_clusters = normalize_cluster_labels(partition.membership, prefix="C")
                reference_clusters = assign_reference_clusters(genome_coords, genome_clusters, ref_coords)
                point_cluster_names = genome_clusters + reference_clusters

                cluster_order = sorted(set(genome_clusters))
                cluster_colors = {cid: PALETTE[i % len(PALETTE)] for i, cid in enumerate(cluster_order)}
                point_colors = [cluster_colors[cid] for cid in point_cluster_names]

                genome_counts = Counter(genome_clusters)
                total_genomes = len(genome_clusters)
                labels = []
                hulls = []
                grouped_pts: Dict[str, List[Tuple[float, float]]] = defaultdict(list)
                for idx, cid in enumerate(genome_clusters):
                    grouped_pts[cid].append((float(genome_coords[idx, 0]), float(genome_coords[idx, 1])))

                for idx, cid in enumerate(cluster_order):
                    pts = grouped_pts[cid]
                    xs = [p[0] for p in pts]
                    ys = [p[1] for p in pts]
                    labels.append(
                        {
                            "cluster": cid,
                            "x": float(sum(xs) / len(xs)),
                            "y": float(sum(ys) / len(ys)),
                            "genomes": genome_counts[cid],
                            "pct": (100.0 * genome_counts[cid] / total_genomes) if total_genomes else 0.0,
                            "offset_x": ((idx % 3) - 1) * 18,
                            "offset_y": ((idx % 4) - 1.5) * 12,
                        }
                    )
                    if len(pts) >= 3:
                        hull = convex_hull(pts)
                        if len(hull) >= 3:
                            hulls.append({"cluster": cid, "color": cluster_colors[cid], "points": hull})

                state_path = states_dir / slug
                state_path.mkdir(parents=True, exist_ok=True)
                with (state_path / "coords.csv").open("w", newline="") as fh:
                    writer = csv.writer(fh)
                    writer.writerow(["ID", "AF_UMAP1", "AF_UMAP2"])
                    for seq_id, (x, y) in zip(all_ids, coords):
                        writer.writerow([seq_id, f"{float(x):.6f}", f"{float(y):.6f}"])
                assignment_rows = []
                for point, cluster_name in zip(all_points, point_cluster_names):
                    assignment_rows.append(
                        {
                            "id": point["id"],
                            "label": point["label"],
                            "item_class": point["item_class"],
                            "cluster": cluster_name,
                        }
                    )
                write_tsv(state_path / "assignments.tsv", assignment_rows, ["id", "label", "item_class", "cluster"])
                cluster_rows = [
                    {
                        "cluster": cid,
                        "genomes": genome_counts[cid],
                        "pct_genomes": f"{(100.0 * genome_counts[cid] / total_genomes) if total_genomes else 0.0:.6f}",
                    }
                    for cid in cluster_order
                ]
                write_tsv(state_path / "cluster_counts.tsv", cluster_rows, ["cluster", "genomes", "pct_genomes"])

                states_payload[f"{kmer}|{neighbors}|{min_dist}"] = {
                    "state_slug": slug,
                    "coords": [[round(float(x), 6), round(float(y), 6)] for x, y in coords],
                    "point_cluster_names": point_cluster_names,
                    "point_colors": point_colors,
                    "clusters": [
                        {
                            "name": cid,
                            "color": cluster_colors[cid],
                            "genomes": genome_counts[cid],
                            "pct": (100.0 * genome_counts[cid] / total_genomes) if total_genomes else 0.0,
                        }
                        for cid in cluster_order
                    ],
                    "cluster_count": len(cluster_order),
                    "labels": labels,
                    "hulls": hulls,
                }

    payload = {
        "kmers": kmers,
        "neighbors": neighbors_values,
        "min_dists": min_dist_values,
        "points": all_points,
        "states": states_payload,
        "genome_count": len(genome_ids),
        "reference_count": len(ref_ids),
        "default_kmer_index": kmers.index(args.default_kmer) if args.default_kmer in kmers else 0,
        "default_neighbor_index": neighbors_values.index(args.default_neighbors) if args.default_neighbors in neighbors_values else 0,
        "default_min_dist_index": min_dist_values.index(args.default_min_dist) if args.default_min_dist in min_dist_values else 0,
    }

    (outdir / "af_leiden_parameter_widget.html").write_text(
        html_template(json.dumps(payload, separators=(",", ":")), "Alignment-Free UMAP Parameter Explorer"),
        encoding="utf-8",
    )
    with (outdir / "summary.txt").open("w") as fh:
        fh.write(f"input_fasta\t{Path(args.input_fasta).resolve()}\n")
        fh.write(f"reference_fasta\t{Path(args.reference_fasta).resolve()}\n")
        fh.write(f"genomes\t{len(genome_ids)}\n")
        fh.write(f"references\t{len(ref_ids)}\n")
        fh.write(f"kmer_values\t{','.join(str(x) for x in kmers)}\n")
        fh.write(f"neighbors_values\t{','.join(str(x) for x in neighbors_values)}\n")
        fh.write(f"min_dist_values\t{','.join(str(x) for x in min_dist_values)}\n")
        fh.write(f"leiden_neighbor_k\t{args.leiden_neighbor_k}\n")
        fh.write(f"leiden_resolution\t{args.leiden_resolution}\n")
        fh.write(f"states_dir\t{states_dir}\n")
