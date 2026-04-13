from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
from Bio import SeqIO


def load_fasta(path: Path) -> tuple[list[str], list[str]]:
    ids, seqs = [], []
    with path.open() as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            ids.append(str(rec.id))
            seqs.append(str(rec.seq).upper().replace("-", ""))
    if not ids:
        raise RuntimeError(f"No sequences found in {path}")
    return ids, seqs


def all_kmers(k: int) -> list[str]:
    kmers = [""]
    for _ in range(k):
        kmers = [prefix + base for prefix in kmers for base in "ACGT"]
    return kmers


def kmer_vector(seq: str, k: int, kmers: list[str]) -> np.ndarray:
    clean = "".join(base for base in seq.upper() if base in "ACGT")
    counts = Counter(clean[i : i + k] for i in range(len(clean) - k + 1))
    vec = np.array([counts.get(km, 0) for km in kmers], dtype=np.float32)
    total = float(vec.sum())
    return vec / total if total > 0 else vec


def build_kmer_matrix(seqs: list[str], k: int) -> np.ndarray:
    kmers = all_kmers(k)
    return np.vstack([kmer_vector(seq, k, kmers) for seq in seqs])
