"""Tests for the recombination-screen fork, `collapse-panel`.

The load-bearing guarantees are (a) focal ids are never dropped regardless of
threshold, (b) the emitted panel satisfies the same contract `align-panel`
expects of `sample-panel` output, and (c) there is no per-cluster quota.
"""
from __future__ import annotations

import csv
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from seqscape.collapse_panel import collapse, jaccard, jaccard_to_identity, kmer_set, run


class Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


def _mutate(seq: str, n: int, rng: random.Random) -> str:
    s = list(seq)
    for pos in rng.sample(range(len(s)), n):
        s[pos] = rng.choice([b for b in "ACGT" if b != s[pos]])
    return "".join(s)


def _build_state(root: Path) -> tuple[Path, list[str], list[str]]:
    """Two tight groups of near-identical genomes plus a divergent focal pair."""
    rng = random.Random(0)
    base = "".join(rng.choice("ACGT") for _ in range(2000))
    other = "".join(rng.choice("ACGT") for _ in range(2000))

    seqs: dict[str, str] = {}
    clusters: dict[str, str] = {}
    for i in range(8):                       # cluster 0: near-identical to base
        sid = f"G{i:02d}"
        seqs[sid] = _mutate(base, 4, rng)
        clusters[sid] = "0"
    for i in range(8, 14):                   # cluster 1: near-identical to other
        sid = f"G{i:02d}"
        seqs[sid] = _mutate(other, 4, rng)
        clusters[sid] = "1"
    focal = ["F00", "F01"]                   # focal: in cluster 0, near-identical
    for sid in focal:                        # to base AND to each other
        seqs[sid] = _mutate(base, 6, rng)
        clusters[sid] = "0"
    refs = ["R00"]
    seqs["R00"] = _mutate(other, 10, rng)
    clusters["R00"] = "1"

    state = root / "explorer" / "states" / "k6_n15_d0p1_r1e00"
    state.mkdir(parents=True)
    with (state / "assignments.tsv").open("w") as fh:
        fh.write("id\tlabel\titem_class\tcluster\n")
        for sid in seqs:
            cls = "reference" if sid in refs else "genome"
            fh.write(f"{sid}\t{sid}\t{cls}\t{clusters[sid]}\n")
    with (state / "coords.csv").open("w") as fh:
        fh.write("ID,AF_UMAP1,AF_UMAP2\n")
        for i, sid in enumerate(seqs):
            fh.write(f"{sid},{float(i)},{float(i)}\n")

    with (root / "genomes.fasta").open("w") as fh:
        for sid, s in seqs.items():
            if sid not in refs:
                fh.write(f">{sid}\n{s}\n")
    with (root / "refs.fasta").open("w") as fh:
        for sid in refs:
            fh.write(f">{sid}\n{seqs[sid]}\n")
    (root / "focal.txt").write_text("\n".join(focal) + "\n")
    return root / "explorer", focal, refs


def _run(root: Path, explorer: Path, threshold: float, focal_file: bool = True) -> Path:
    outdir = root / f"out_{threshold}_{focal_file}"
    run(Args(
        input_fasta=str(root / "genomes.fasta"),
        reference_fasta=str(root / "refs.fasta"),
        umap_explorer_dir=str(explorer),
        chosen_kmer=6, chosen_neighbors=15, chosen_min_dist=0.1,
        chosen_leiden_resolution=1.0,
        identity_threshold=threshold, jaccard_k=15,
        focal_ids=str(root / "focal.txt") if focal_file else "",
        outdir=str(outdir),
    ))
    return outdir


class CollapsePanelTests(unittest.TestCase):
    def test_jaccard_to_identity_monotonic_and_bounded(self) -> None:
        self.assertEqual(jaccard_to_identity(0.0, 15), 0.0)
        self.assertEqual(jaccard_to_identity(1.0, 15), 100.0)
        vals = [jaccard_to_identity(j, 15) for j in (0.1, 0.3, 0.6, 0.9)]
        self.assertEqual(vals, sorted(vals))

    def test_identical_sequences_score_100(self) -> None:
        s = "ACGT" * 100
        self.assertAlmostEqual(jaccard_to_identity(jaccard(kmer_set(s, 15), kmer_set(s, 15)), 15), 100.0)

    def test_protected_ids_never_absorbed(self) -> None:
        """Even mutually near-identical focal ids must all survive."""
        rng = random.Random(1)
        base = "".join(rng.choice("ACGT") for _ in range(1000))
        seqs = {f"x{i}": _mutate(base, 2, rng) for i in range(5)}
        kms = {k: kmer_set(v, 15) for k, v in seqs.items()}
        prot = {"x0", "x1"}
        kept, members = collapse(list(seqs), kms, 90.0, 15, protected=prot)
        self.assertTrue(prot <= set(kept))
        # a protected id must never stand in for another genome
        for p in prot:
            self.assertEqual(members.get(p, []), [])

    def test_focal_retained_at_every_threshold(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            explorer, focal, _ = _build_state(root)
            for thr in (90.0, 95.0, 99.0, 99.9):
                outdir = _run(root, explorer, thr)
                m = json.loads((outdir / "validation_metrics.json").read_text())
                self.assertEqual(m["focal_retained"], len(focal),
                                 f"focal genome dropped at threshold {thr}")

    def test_panel_contract_matches_sample_panel(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            explorer, focal, refs = _build_state(root)
            outdir = _run(root, explorer, 95.0)

            with (outdir / "representative_panel_manifest.tsv").open(newline="") as fh:
                man = list(csv.DictReader(fh, delimiter="\t"))
            fasta_ids = {ln[1:].split()[0] for ln in
                         (outdir / "representative_panel.fasta").read_text().splitlines()
                         if ln.startswith(">")}
            # align-panel resolves every manifest id against the panel FASTA
            self.assertEqual({r["id"] for r in man}, fasta_ids)
            for col in ("id", "label", "item_class", "source_cluster"):
                self.assertIn(col, man[0])
            # references always survive
            self.assertTrue(set(refs) <= fasta_ids)

    def test_collapse_actually_dereplicates(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            explorer, focal, _ = _build_state(root)
            outdir = _run(root, explorer, 95.0)
            m = json.loads((outdir / "validation_metrics.json").read_text())
            # 14 near-identical genomes in two groups collapse hard; focal are exempt
            self.assertGreater(m["genomes_absorbed"], 8)
            self.assertLess(m["panel_size_genomes"], m["N"])

    def test_no_per_cluster_quota(self) -> None:
        """Retention must track cluster CONTENT, not cluster SIZE.

        Cluster 0 has more genomes than cluster 1 but is not more diverse, so a
        proportional scheme would give it more slots. Collapse must not.
        """
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            explorer, focal, _ = _build_state(root)
            outdir = _run(root, explorer, 95.0, focal_file=False)
            with (outdir / "source_cluster_retention.tsv").open(newline="") as fh:
                rows = {r["cluster"]: r for r in csv.DictReader(fh, delimiter="\t")}
            c0, c1 = rows["0"], rows["1"]
            self.assertGreater(int(c0["genomes_in_cluster"]), int(c1["genomes_in_cluster"]))
            # both clusters are internally near-identical -> both collapse to ~1
            self.assertLessEqual(int(c0["retained"]), 3)
            self.assertLessEqual(int(c1["retained"]), 3)

    def test_unknown_focal_id_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            explorer, _, _ = _build_state(root)
            (root / "focal.txt").write_text("NOT_A_REAL_ID\n")
            with self.assertRaises(RuntimeError):
                _run(root, explorer, 95.0)


if __name__ == "__main__":
    unittest.main()
