from __future__ import annotations

import importlib.util
import os
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
TEST_TMP = Path(tempfile.gettempdir()) / "seqscape_test_runtime"
TEST_TMP.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(TEST_TMP / "mplconfig"))
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _has_af_stack() -> bool:
    needed = ["igraph", "leidenalg", "umap"]
    return all(importlib.util.find_spec(name) is not None for name in needed)


def _write_fasta(path: Path, records: list[tuple[str, str, str]]) -> None:
    with path.open("w") as fh:
        for rec_id, desc, seq in records:
            header = f">{rec_id} {desc}".rstrip()
            fh.write(header + "\n")
            fh.write(seq + "\n")


@unittest.skipUnless(_has_af_stack(), "AF stack not installed")
class PipelineSmokeTests(unittest.TestCase):
    def test_synthetic_pipeline(self) -> None:
        from seqscape import umap_explorer, align_panel, distance_explorer, sampling

        with tempfile.TemporaryDirectory(prefix="seqscape_test_") as td:
            root = Path(td)
            genome_fa = root / "genomes.fasta"
            ref_fa = root / "refs.fasta"
            af_out = root / "af"
            panel_out = root / "panel"
            align_out = root / "aligned"
            review_out = root / "review"

            r1 = "AAAAACCCCCAAAAACCCCCAAAA"
            r2 = "GGGGGTTTTTGGGGGTTTTTGGGG"
            _write_fasta(
                genome_fa,
                [
                    ("g1", "genome", r1),
                    ("g2", "genome", "AAAAACCCCCAAAATCCCCCAAAA"),
                    ("g3", "genome", r2),
                    ("g4", "genome", "GGGGGTTTTTGGGGATTTTTGGGG"),
                ],
            )
            _write_fasta(
                ref_fa,
                [
                    ("r1", "reference", r1),
                    ("r2", "reference", r2),
                ],
            )

            umap_explorer.run(
                Namespace(
                    input_fasta=str(genome_fa),
                    reference_fasta=str(ref_fa),
                    cg_map_tsv="",
                    kmer_values="3",
                    neighbors_values="3",
                    min_dist_values="0.3",
                    seed=42,
                    leiden_neighbor_k=3,
                    leiden_resolution=0.01,
                    leiden_resolution_values="0.01",
                    default_kmer=3,
                    default_neighbors=3,
                    default_min_dist=0.3,
                    default_leiden_resolution=0.01,
                    validation_subsample_size=4,
                    validation_replicates=1,
                    validation_neighbors="2",
                    outdir=str(af_out),
                )
            )
            self.assertTrue((af_out / "umap_explorer.html").is_file())
            self.assertTrue((af_out / "states" / "k3_n3_d0p3_r1e-02" / "assignments.tsv").is_file())
            self.assertTrue((af_out / "states" / "k3_n3_d0p3_r1e-02" / "validation_metrics.json").is_file())

            sampling.run(
                Namespace(
                    input_fasta=str(genome_fa),
                    reference_fasta=str(ref_fa),
                    umap_explorer_dir=str(af_out),
                    chosen_kmer=3,
                    chosen_neighbors=3,
                    chosen_min_dist=0.3,
                    chosen_leiden_resolution=0.01,
                    sample_size=5,
                    outdir=str(panel_out),
                )
            )
            self.assertTrue((panel_out / "representative_panel.fasta").is_file())
            self.assertTrue((panel_out / "representative_panel_manifest.tsv").is_file())
            self.assertTrue((panel_out / "summary.txt").is_file())
            self.assertTrue((panel_out / "validation_metrics.json").is_file())

            align_panel.run(
                Namespace(
                    panel_fasta=str(panel_out / "representative_panel.fasta"),
                    manifest_tsv=str(panel_out / "representative_panel_manifest.tsv"),
                    aligner="pairwisealigner",
                    mafft="",
                    muscle="",
                    jobs=1,
                    progress_every=0,
                    outdir=str(align_out),
                )
            )
            self.assertTrue((align_out / "matrix_identity.tsv").is_file())
            self.assertTrue((align_out / "matrix_distance.tsv").is_file())
            self.assertTrue((align_out / "pcoa_coords.csv").is_file())
            self.assertTrue((align_out / "summary.txt").is_file())
            self.assertTrue((align_out / "validation_metrics.json").is_file())

            distance_explorer.run(
                Namespace(
                    panel_fasta=str(panel_out / "representative_panel.fasta"),
                    identity_matrix=str(align_out / "matrix_identity.tsv"),
                    distance_matrix=str(align_out / "matrix_distance.tsv"),
                    manifest_tsv=str(panel_out / "representative_panel_manifest.tsv"),
                    umap_explorer_dir="",
                    thresholds="0.05,0.10",
                    novel_threshold=92.0,
                    umap_neighbors=3,
                    umap_min_dist=0.3,
                    umap_spread=0.8,
                    umap_seed=42,
                    outdir=str(review_out),
                )
            )
            self.assertTrue((review_out / "distance_explorer.html").is_file())
            self.assertTrue((review_out / "heatmap_genome_x_reference.html").is_file())
            self.assertTrue((review_out / "heatmap_panel_square.html").is_file())
            self.assertTrue((review_out / "summary.txt").is_file())
            self.assertTrue((review_out / "validation_metrics.json").is_file())


if __name__ == "__main__":
    unittest.main()
