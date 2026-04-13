from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .clustering import fit_agglomerative, normalize_cluster_labels
from .heatmap import maybe_heatmap, write_html_heatmap
from .io_utils import write_tsv
from .ordination import pcoa, run_precomputed_umap
from .palettes import PALETTE


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


def html_template(payload_json: str, title: str) -> str:
    template = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>__TITLE__</title>
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
      margin: 0 0 14px;
      color: #555;
      font-size: 14px;
    }}
    .toolbar {{
      display: flex;
      gap: 10px;
      margin: 0 0 14px;
      flex-wrap: wrap;
    }}
    .toggle {{
      border: 1px solid #d0c8b9;
      background: white;
      color: #333;
      border-radius: 999px;
      padding: 8px 14px;
      cursor: pointer;
      font-size: 13px;
    }}
    .toggle.active {{
      background: #222;
      color: white;
      border-color: #222;
    }}
    .slider-wrap {{
      background: white;
      border: 1px solid #ddd6ca;
      border-radius: 12px;
      padding: 12px 14px;
      margin: 0 0 14px;
      font-size: 13px;
    }}
    .slider-wrap label {{
      display: block;
      font-weight: 700;
      margin-bottom: 6px;
    }}
    .slider-wrap input {{
      width: 100%;
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
      border: 1px solid rgba(0,0,0,0.15);
    }}
    .meta {{
      margin: 0 0 12px;
      color: #555;
      font-size: 12px;
    }}
    canvas {{
      width: 100%;
      height: auto;
      background: #fffdf8;
      border: 1px solid #d8d2c6;
      border-radius: 16px;
      box-shadow: 0 8px 24px rgba(0,0,0,0.06);
    }}
  </style>
</head>
<body>
  <div class="wrap">
    <h1>__TITLE__</h1>
    <p class="sub">Use the threshold slider to inspect complete-linkage agglomerative clustering on the aligned panel. Toggle between alignment-based PCoA, alignment-based UMAP, and best-vs-second-best reference identity, and color by either the original Leiden clusters or the selected agglomerative clusters.</p>
    <div class="toolbar">
      <button class="toggle active" id="btn-pcoa">PCoA</button>
      <button class="toggle" id="btn-umap">Distance UMAP</button>
      <button class="toggle" id="btn-best">Best vs 2nd</button>
    </div>
    <div class="toolbar">
      <button class="toggle active" id="btn-ref-include">Refs include self</button>
      <button class="toggle" id="btn-ref-exclude">Refs exclude self</button>
    </div>
    <div class="toolbar">
      <button class="toggle active" id="btn-leiden">Color by Leiden</button>
      <button class="toggle" id="btn-agg">Color by Agglomerative</button>
    </div>
    <div class="slider-wrap">
      <label for="threshold">Complete-linkage threshold: <span id="threshold-value"></span></label>
      <input id="threshold" type="range" min="0" max="0" step="1">
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
    <div class="meta" id="meta"></div>
    <div class="legend" id="legend"></div>
    <canvas id="plot" width="1450" height="980"></canvas>
  </div>
  <script>
    const payload = __PAYLOAD__;
    const canvas = document.getElementById('plot');
    const ctx = canvas.getContext('2d');
    const legend = document.getElementById('legend');
    const meta = document.getElementById('meta');
    const btnPcoa = document.getElementById('btn-pcoa');
    const btnUmap = document.getElementById('btn-umap');
    const btnBest = document.getElementById('btn-best');
    const btnRefInclude = document.getElementById('btn-ref-include');
    const btnRefExclude = document.getElementById('btn-ref-exclude');
    const btnLeiden = document.getElementById('btn-leiden');
    const btnAgg = document.getElementById('btn-agg');
    const thresholdSlider = document.getElementById('threshold');
    const thresholdValue = document.getElementById('threshold-value');
    const minLengthInput = document.getElementById('min-length');
    const maxLengthInput = document.getElementById('max-length');
    const thresholds = payload.thresholds;
    thresholdSlider.max = String(thresholds.length - 1);
    thresholdSlider.value = String(payload.default_threshold_index);
    const allLengths = payload.points.map(p => Number(p.length_bp || 0)).filter(v => Number.isFinite(v) && v > 0);
    const globalMinLength = allLengths.length ? Math.min(...allLengths) : 0;
    const globalMaxLength = allLengths.length ? Math.max(...allLengths) : 0;
    minLengthInput.value = String(globalMinLength);
    maxLengthInput.value = String(globalMaxLength);

    let mode = 'pcoa';
    let referenceBestMode = 'include';
    let scheme = 'leiden_cluster';
    let hover = null;
    let hitboxes = [];
    const box = {x: 90, y: 90, w: 1220, h: 800};

    function currentThreshold() {
      return thresholds[Number(thresholdSlider.value)];
    }
    function pointCluster(point) {
      if (scheme === 'leiden_cluster') return point.leiden_cluster;
      return point.agglomerative[currentThreshold()];
    }
    function currentView() {
      if (scheme === 'leiden_cluster') return payload.cluster_views.leiden_cluster[mode];
      return payload.cluster_views.agglomerative[currentThreshold()][mode];
    }
    function currentColorView() {
      if (scheme === 'leiden_cluster') return payload.cluster_views.leiden_cluster.pcoa;
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
    function bounds(points) {
      const xs = points.map(p => pointCoords(p)[0]);
      const ys = points.map(p => pointCoords(p)[1]);
      return {minX: Math.min(...xs), maxX: Math.max(...xs), minY: Math.min(...ys), maxY: Math.max(...ys)};
    }
    function sx(x, b) { return box.x + ((x - b.minX) / (b.maxX - b.minX || 1)) * box.w; }
    function sy(y, b) { return box.y + box.h - ((y - b.minY) / (b.maxY - b.minY || 1)) * box.h; }
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
        const cid = pointCluster(p);
        if (!grouped.has(cid)) grouped.set(cid, []);
        grouped.get(cid).push(pointCoords(p));
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
      const detail = scheme === 'leiden_cluster' ? 'original Leiden clusters' : `complete-linkage agglomerative clusters @ ${thr}`;
      const refDetail = mode === 'best2' ? (referenceBestMode === 'include' ? 'refs use self as best' : 'refs use nearest non-self refs') : 'reference best-vs-second convention hidden';
      const {minLen, maxLen} = currentLengthBounds();
      meta.textContent = `${mode.toUpperCase()} | ${detail} | ${refDetail} | best-ref cutoff ${payload.novel_threshold}% | visible points ${points.length}/${payload.points.length} | length ${minLen}-${maxLen} nt`;
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
        item.best_ref ? `best: ${item.best_ref} (${item.best_identity.toFixed(2)}%)` : '',
        (mode === 'best2' && item.item_class === 'reference') ? (referenceBestMode === 'include' ? 'refs plotted with self-match on x' : 'refs plotted without self-match') : '',
      ].filter(Boolean);
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
    function render() {
      const points = filteredPoints();
      const view = computeView(points);
      if (!points.length) {
        thresholdValue.textContent = currentThreshold();
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
      hitboxes = [];
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.fillStyle = '#fffdf8';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.strokeStyle = '#ddd6ca';
      ctx.lineWidth = 1;
      ctx.strokeRect(box.x, box.y, box.w, box.h);
      if (mode === 'best2') {
        const thrX = sx(payload.novel_threshold, b);
        ctx.strokeStyle = '#666';
        ctx.setLineDash([6, 6]);
        ctx.beginPath();
        ctx.moveTo(thrX, box.y);
        ctx.lineTo(thrX, box.y + box.h);
        ctx.stroke();
        ctx.setLineDash([]);
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
          hitboxes.push({x, y, r: 6, item: p});
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
      const xlab = mode === 'best2' ? 'Best reference identity (%)' : (mode === 'pcoa' ? payload.pcoa_axis_labels.x : 'UMAP 1');
      const ylab = mode === 'best2' ? 'Second-best reference identity (%)' : (mode === 'pcoa' ? payload.pcoa_axis_labels.y : 'UMAP 2');
      ctx.fillText(xlab, box.x + box.w / 2 - 40, box.y + box.h + 46);
      ctx.save();
      ctx.translate(28, box.y + box.h / 2 + 30);
      ctx.rotate(-Math.PI / 2);
      ctx.fillText(ylab, 0, 0);
      ctx.restore();
      if (hover) drawTooltip(hover.item, hover.x, hover.y);
    }
    canvas.addEventListener('mousemove', ev => {
      const rect = canvas.getBoundingClientRect();
      const x = (ev.clientX - rect.left) * (canvas.width / rect.width);
      const y = (ev.clientY - rect.top) * (canvas.height / rect.height);
      hover = null;
      for (const hb of hitboxes) {
        const dx = x - hb.x;
        const dy = y - hb.y;
        if ((dx * dx) + (dy * dy) <= hb.r * hb.r) {
          hover = {item: hb.item, x, y};
          break;
        }
      }
      render();
    });
    canvas.addEventListener('mouseleave', () => { hover = null; render(); });
    function setMode(next) {
      mode = next;
      btnPcoa.classList.toggle('active', next === 'pcoa');
      btnUmap.classList.toggle('active', next === 'umap');
      btnBest.classList.toggle('active', next === 'best2');
      setLegend();
      render();
    }
    function setScheme(next) {
      scheme = next;
      btnLeiden.classList.toggle('active', next === 'leiden_cluster');
      btnAgg.classList.toggle('active', next === 'agglomerative');
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
    btnPcoa.addEventListener('click', () => setMode('pcoa'));
    btnUmap.addEventListener('click', () => setMode('umap'));
    btnBest.addEventListener('click', () => setMode('best2'));
    btnRefInclude.addEventListener('click', () => setReferenceBestMode('include'));
    btnRefExclude.addEventListener('click', () => setReferenceBestMode('exclude'));
    btnLeiden.addEventListener('click', () => setScheme('leiden_cluster'));
    btnAgg.addEventListener('click', () => setScheme('agglomerative'));
    thresholdSlider.addEventListener('input', () => { setLegend(); render(); });
    minLengthInput.addEventListener('input', () => { setLegend(); render(); });
    maxLengthInput.addEventListener('input', () => { setLegend(); render(); });
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
    total_pcoa = float(np.sum(pcoa_eigvals)) if pcoa_eigvals.size else 0.0
    pcoa_pct1 = (100.0 * float(pcoa_eigvals[0]) / total_pcoa) if total_pcoa > 0.0 and pcoa_eigvals.size >= 1 else 0.0
    pcoa_pct2 = (100.0 * float(pcoa_eigvals[1]) / total_pcoa) if total_pcoa > 0.0 and pcoa_eigvals.size >= 2 else 0.0
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
    }
    threshold_rows = []
    long_agglom_rows = []
    for threshold in thresholds:
        raw = fit_agglomerative(dist, "complete", threshold)
        cluster_names = normalize_cluster_labels(raw, prefix="A")
        agglom_order = sorted(set(cluster_names), key=lambda cid: cid)
        for point, cid in zip(points, cluster_names):
            point["agglomerative"][str(threshold)] = cid
        genome_counts = Counter(cid for row, cid in zip(ordered_manifest, cluster_names) if row["item_class"] == "genome")
        threshold_rows.append(
            {
                "threshold": f"{threshold:.6f}",
                "cluster_count": len(agglom_order),
                "largest_genome_cluster": max(genome_counts.values()) if genome_counts else 0,
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
    build_heatmaps(ids, ident, ordered_manifest, outdir)

    payload = {
        "thresholds": [str(x) for x in thresholds],
        "default_threshold_index": min(len(thresholds) - 1, thresholds.index(0.06) if 0.06 in thresholds else 0),
        "novel_threshold": args.novel_threshold,
        "pcoa_axis_labels": {
            "x": f"PCoA 1 ({pcoa_pct1:.1f}%)",
            "y": f"PCoA 2 ({pcoa_pct2:.1f}%)",
        },
        "points": points,
        "cluster_views": cluster_views,
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
        fh.write(f"thresholds\t{','.join(str(x) for x in thresholds)}\n")
        fh.write(f"novel_threshold\t{args.novel_threshold}\n")
        fh.write(f"umap_neighbors\t{args.umap_neighbors}\n")
        fh.write(f"umap_min_dist\t{args.umap_min_dist}\n")
        fh.write(f"umap_spread\t{args.umap_spread}\n")
