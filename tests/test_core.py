from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from seqscape.clustering import normalize_cluster_labels
from seqscape.kmer import all_kmers, build_kmer_matrix
from seqscape.sampling import largest_remainder_alloc


class CoreTests(unittest.TestCase):
    def test_largest_remainder_alloc_preserves_total(self) -> None:
        alloc = largest_remainder_alloc({"A": 7, "B": 2, "C": 1}, total=10)
        self.assertEqual(sum(alloc.values()), 10)
        self.assertGreaterEqual(alloc["A"], alloc["B"])
        self.assertGreaterEqual(alloc["B"], alloc["C"])

    def test_normalize_cluster_labels_is_stable(self) -> None:
        labels = normalize_cluster_labels([3, 3, 8, 8, 2], prefix="A")
        self.assertEqual(labels, ["A001", "A001", "A002", "A002", "A003"])

    def test_build_kmer_matrix_shape(self) -> None:
        kmers = all_kmers(3)
        mat = build_kmer_matrix(["AAACCCGGG", "AAAGGGCCC"], 3)
        self.assertEqual(len(kmers), 64)
        self.assertEqual(mat.shape, (2, 64))
        self.assertTrue(np.all(mat >= 0))


if __name__ == "__main__":
    unittest.main()
