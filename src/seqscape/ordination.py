from __future__ import annotations

import os
import tempfile
from pathlib import Path
import sys

import numpy as np


def _add_vendor_umap_to_path() -> None:
    vendor_root = Path(__file__).resolve().parents[2] / "third_party" / "seqspace_vendor"
    if vendor_root.is_dir():
        sys.path.insert(0, str(vendor_root))


def _prepare_runtime_env() -> None:
    tmp = Path(tempfile.gettempdir())
    numba_cache = tmp / "seqscape_numba_cache"
    mpl_dir = tmp / "seqscape_mplconfig"
    numba_cache.mkdir(parents=True, exist_ok=True)
    mpl_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("NUMBA_CACHE_DIR", str(numba_cache))
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))


def pcoa(dist: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = dist.shape[0]
    d2 = dist**2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    pos = eigvals > 1e-12
    eigvals = eigvals[pos]
    eigvecs = eigvecs[:, pos]
    coords = eigvecs * np.sqrt(eigvals)
    return coords, eigvals


def run_precomputed_umap(
    dist: np.ndarray,
    neighbors: int,
    min_dist: float,
    spread: float,
    seed: int,
) -> np.ndarray:
    _add_vendor_umap_to_path()
    _prepare_runtime_env()
    import umap  # type: ignore

    reducer = umap.UMAP(
        n_neighbors=neighbors,
        min_dist=min_dist,
        spread=spread,
        n_components=2,
        metric="precomputed",
        random_state=seed,
    )
    return reducer.fit_transform(dist)


def run_euclidean_umap(matrix: np.ndarray, neighbors: int, min_dist: float, seed: int) -> np.ndarray:
    _add_vendor_umap_to_path()
    _prepare_runtime_env()
    import umap  # type: ignore

    reducer = umap.UMAP(
        n_neighbors=neighbors,
        min_dist=min_dist,
        n_components=2,
        metric="euclidean",
        random_state=seed,
    )
    return reducer.fit_transform(matrix)
