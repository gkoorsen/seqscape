from __future__ import annotations

import math
import random
from collections import Counter
from typing import Dict, Iterable, Sequence

import numpy as np
from scipy.stats import pearsonr
from sklearn.manifold import trustworthiness
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def full_pairwise_count(n: int) -> int:
    return n * (n - 1) // 2


def panel_fraction_pct(panel_size: int, total_size: int) -> float:
    return (100.0 * panel_size / total_size) if total_size else 0.0


def compression_ratio(total_size: int, panel_size: int) -> float:
    panel_pairs = full_pairwise_count(panel_size)
    return (full_pairwise_count(total_size) / panel_pairs) if panel_pairs else float("nan")


def pcoa_variance_explained(eigvals: np.ndarray) -> tuple[float, float]:
    pos = eigvals[eigvals > 1e-12]
    total = float(np.sum(pos)) if pos.size else 0.0
    if total <= 0.0:
        return 0.0, 0.0
    axis1 = (100.0 * float(pos[0]) / total) if pos.size >= 1 else 0.0
    axis2 = (100.0 * float(pos[1]) / total) if pos.size >= 2 else 0.0
    return axis1, axis2


def pearson_cluster_proportion_preservation(
    full_cluster_counts: Dict[str, int], panel_cluster_counts: Dict[str, int]
) -> tuple[float, float]:
    total_full = sum(full_cluster_counts.values())
    total_panel = sum(panel_cluster_counts.values())
    clusters = sorted(set(full_cluster_counts) | set(panel_cluster_counts))
    full_vals = [
        (100.0 * full_cluster_counts.get(cluster, 0) / total_full) if total_full else 0.0
        for cluster in clusters
    ]
    panel_vals = [
        (100.0 * panel_cluster_counts.get(cluster, 0) / total_panel) if total_panel else 0.0
        for cluster in clusters
    ]
    if len(clusters) < 2:
        return 1.0, 0.0
    return pearsonr(full_vals, panel_vals)


def choose_closest_cluster_threshold(summary_rows: Sequence[dict], target_cluster_count: int) -> dict:
    return sorted(
        summary_rows,
        key=lambda row: (abs(int(row["cluster_count"]) - target_cluster_count), float(row["threshold"])),
    )[0]


def nmi(leiden_labels: Sequence[str], distance_labels: Sequence[str]) -> float:
    return float(normalized_mutual_info_score(leiden_labels, distance_labels))


def ari(source_labels: Sequence[str], distance_labels: Sequence[str]) -> float:
    return float(adjusted_rand_score(source_labels, distance_labels))


def _stratified_ids(
    ids_by_cluster: Dict[str, Sequence[int]],
    total_pick: int,
    seed: int,
    replicate: int,
) -> list[int]:
    rng = random.Random(seed + replicate)
    total = sum(len(v) for v in ids_by_cluster.values())
    quotas: Dict[str, int] = {}
    remainders: list[tuple[float, str]] = []
    floor_total = 0
    for cluster, ids in ids_by_cluster.items():
        exact = total_pick * len(ids) / total
        base = min(len(ids), int(math.floor(exact)))
        quotas[cluster] = base
        floor_total += base
        remainders.append((exact - math.floor(exact), cluster))
    remaining = total_pick - floor_total
    for _, cluster in sorted(remainders, reverse=True):
        if remaining <= 0:
            break
        if quotas[cluster] < len(ids_by_cluster[cluster]):
            quotas[cluster] += 1
            remaining -= 1
    selected: list[int] = []
    for cluster, ids in ids_by_cluster.items():
        q = quotas[cluster]
        if q <= 0:
            continue
        if q >= len(ids):
            picked = list(ids)
        else:
            picked = rng.sample(list(ids), q)
        selected.extend(picked)
    return selected[:total_pick]


def trustworthiness_summary(
    feature_matrix: np.ndarray,
    coords: np.ndarray,
    cluster_labels: Sequence[str],
    neighbors_list: Sequence[int],
    subsample_size: int = 1000,
    replicates: int = 1,
    seed: int = 42,
) -> dict[str, dict]:
    n = feature_matrix.shape[0]
    if n == 0:
        return {}
    labels = list(cluster_labels)
    cluster_to_indices: Dict[str, list[int]] = {}
    for idx, cid in enumerate(labels):
        cluster_to_indices.setdefault(cid, []).append(idx)
    pick_n = min(subsample_size, n)
    scores: Dict[int, list[float]] = {int(k): [] for k in neighbors_list}
    for rep in range(max(1, replicates)):
        picked = _stratified_ids(cluster_to_indices, pick_n, seed, rep)
        sub_matrix = feature_matrix[picked]
        sub_coords = coords[picked]
        for neigh in neighbors_list:
            max_valid = max(1, (len(picked) - 1) // 2)
            k = min(int(neigh), max_valid)
            if k < 1:
                continue
            score = trustworthiness(sub_matrix, sub_coords, n_neighbors=k, metric="euclidean")
            scores[int(neigh)].append(float(score))
    out: dict[str, dict] = {}
    for neigh, vals in scores.items():
        if not vals:
            continue
        out[str(neigh)] = {
            "mean": float(np.mean(vals)),
            "std": float(np.std(vals, ddof=0)),
            "replicates": vals,
        }
    return out


def cluster_size_metrics(cluster_labels: Sequence[str]) -> dict[str, float | int]:
    counts = Counter(cluster_labels)
    total = sum(counts.values())
    largest = max(counts.values()) if counts else 0
    singletons = sum(1 for v in counts.values() if v == 1)
    small_clusters = sum(1 for v in counts.values() if v <= 10)
    return {
        "cluster_count": len(counts),
        "largest_cluster_genomes": largest,
        "largest_cluster_pct": (100.0 * largest / total) if total else 0.0,
        "singleton_cluster_count": singletons,
        "small_cluster_count_le_10": small_clusters,
        "genomes_counted": total,
    }
