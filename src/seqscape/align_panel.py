from __future__ import annotations

import csv
import shutil
from pathlib import Path

import numpy as np

from .alignment import build_identity_matrix, identity_to_distance, resolve_muscle_exe
from .io_utils import load_fasta_records
from .ordination import pcoa


def write_matrix_tsv(path: Path, ids: list[str], matrix: np.ndarray, fmt: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["ID", *ids])
        for seq_id, row in zip(ids, matrix):
            writer.writerow([seq_id, *[fmt.format(float(v)) for v in row]])


def run(args) -> None:
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    panel_records = load_fasta_records(Path(args.panel_fasta).resolve())
    with Path(args.manifest_tsv).resolve().open(newline="") as fh:
        manifest_rows = list(csv.DictReader(fh, delimiter="\t"))
    items = []
    for row in manifest_rows:
        rec = panel_records.get(row["id"])
        if rec is None:
            raise RuntimeError(f"Manifest id missing from panel FASTA: {row['id']}")
        item = dict(row)
        item["sequence"] = str(rec.seq).upper()
        items.append(item)

    muscle_exe = resolve_muscle_exe(args.muscle)
    ids, ident = build_identity_matrix(items, args.aligner, args.mafft, muscle_exe, args.jobs, args.progress_every)
    dist = identity_to_distance(ident)
    write_matrix_tsv(outdir / "matrix_identity.tsv", ids, ident, "{:.6f}")
    write_matrix_tsv(outdir / "matrix_distance.tsv", ids, dist, "{:.6f}")
    coords_raw, _ = pcoa(dist)
    coords = np.zeros((len(ids), 2), dtype=float)
    if coords_raw.shape[1] >= 1:
        coords[:, 0] = coords_raw[:, 0]
    if coords_raw.shape[1] >= 2:
        coords[:, 1] = coords_raw[:, 1]
    with (outdir / "pcoa_coords.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["ID", "PCoA1", "PCoA2"])
        for seq_id, (x, y) in zip(ids, coords):
            writer.writerow([seq_id, f"{float(x):.6f}", f"{float(y):.6f}"])
    shutil.copy2(Path(args.manifest_tsv).resolve(), outdir / "representative_panel_manifest.tsv")
    shutil.copy2(Path(args.panel_fasta).resolve(), outdir / "representative_panel.fasta")
    with (outdir / "summary.txt").open("w") as fh:
        fh.write(f"panel_fasta\t{Path(args.panel_fasta).resolve()}\n")
        fh.write(f"manifest_tsv\t{Path(args.manifest_tsv).resolve()}\n")
        fh.write(f"aligner\t{args.aligner}\n")
        fh.write(f"jobs\t{args.jobs}\n")
        fh.write(f"matrix_identity\t{outdir / 'matrix_identity.tsv'}\n")
        fh.write(f"matrix_distance\t{outdir / 'matrix_distance.tsv'}\n")
        fh.write(f"pcoa_coords\t{outdir / 'pcoa_coords.csv'}\n")
