from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np
from sklearn.cluster import AgglomerativeClustering


def normalize_cluster_labels(labels: Sequence[int] | np.ndarray, prefix: str = "A") -> List[str]:
    mapping: Dict[int, str] = {}
    out: List[str] = []
    next_idx = 1
    iterable = labels.tolist() if hasattr(labels, "tolist") else list(labels)
    for label in iterable:
        if label not in mapping:
            mapping[label] = f"{prefix}{next_idx:03d}"
            next_idx += 1
        out.append(mapping[label])
    return out


def fit_agglomerative(dist: np.ndarray, linkage: str, threshold: float) -> np.ndarray:
    kwargs = {
        "n_clusters": None,
        "linkage": linkage,
        "distance_threshold": threshold,
    }
    try:
        model = AgglomerativeClustering(metric="precomputed", **kwargs)
    except TypeError:
        model = AgglomerativeClustering(affinity="precomputed", **kwargs)
    return model.fit_predict(dist)
