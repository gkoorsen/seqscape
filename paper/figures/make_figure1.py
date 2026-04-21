"""Generate Figure 1: SeqScape workflow overview and interactive outputs.

Panels:
  1A: Four-stage pipeline schematic
  1B: umap-explorer UMAP preview (4,700-seq norovirus run, chosen state)
  1C: distance-explorer preview (Chhabra 302 PCoA + NJ tree)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

FIG_DIR = Path(__file__).parent
OUT_FIG = FIG_DIR / "figure1.png"
OUT_SVG = FIG_DIR / "figure1.svg"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
})


def panel_1a(ax):
    """Four-stage pipeline schematic."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.2)
    ax.axis("off")

    stages = [
        {
            "x": 0.2, "title": "1  umap-explorer",
            "subtitle": "k-mer + UMAP + Leiden",
            "in_label": "FASTA\n(N seqs)",
            "out_label": "UMAP states\n(interactive HTML)",
            "color": "#4e79a7",
        },
        {
            "x": 2.7, "title": "2  sample-panel",
            "subtitle": "proportional,\ncentroid + farthest-point",
            "in_label": "chosen state\n+ panel size",
            "out_label": "representative\npanel FASTA",
            "color": "#f28e2b",
        },
        {
            "x": 5.2, "title": "3  align-panel",
            "subtitle": "MUSCLE / MAFFT\nall-vs-all",
            "in_label": "panel FASTA\n(n seqs)",
            "out_label": "identity +\ndistance matrices",
            "color": "#59a14f",
        },
        {
            "x": 7.7, "title": "4  distance-explorer",
            "subtitle": "PCoA + NJ tree +\nthreshold explorer",
            "in_label": "matrices\n+ references",
            "out_label": "interactive HTML\n(clusters, novelty)",
            "color": "#b07aa1",
        },
    ]

    box_w, box_h = 2.1, 1.5
    y_center = 1.55

    for i, s in enumerate(stages):
        x = s["x"]
        box = FancyBboxPatch(
            (x, y_center - box_h / 2), box_w, box_h,
            boxstyle="round,pad=0.05,rounding_size=0.08",
            facecolor=s["color"], edgecolor="black", alpha=0.85, lw=1.0,
        )
        ax.add_patch(box)
        ax.text(x + box_w / 2, y_center + 0.45, s["title"],
                ha="center", va="center", fontsize=10, weight="bold",
                color="white")
        ax.text(x + box_w / 2, y_center - 0.05, s["subtitle"],
                ha="center", va="center", fontsize=8,
                color="white", style="italic")

        # Input label (above)
        ax.text(x + box_w / 2, y_center + box_h / 2 + 0.35, s["in_label"],
                ha="center", va="center", fontsize=7.5, color="#333333")
        ax.annotate("", xy=(x + box_w / 2, y_center + box_h / 2 + 0.02),
                    xytext=(x + box_w / 2, y_center + box_h / 2 + 0.22),
                    arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8))

        # Output label (below)
        ax.text(x + box_w / 2, y_center - box_h / 2 - 0.35, s["out_label"],
                ha="center", va="center", fontsize=7.5, color="#333333")
        ax.annotate("", xy=(x + box_w / 2, y_center - box_h / 2 - 0.22),
                    xytext=(x + box_w / 2, y_center - box_h / 2 - 0.02),
                    arrowprops=dict(arrowstyle="<-", color="#333333", lw=0.8))

        # Arrow to next stage
        if i < len(stages) - 1:
            nx = stages[i + 1]["x"]
            arr = FancyArrowPatch(
                (x + box_w, y_center), (nx, y_center),
                arrowstyle="-|>", mutation_scale=16, color="#444444", lw=1.5,
            )
            ax.add_patch(arr)

    ax.set_title("A  SeqScape staged workflow", loc="left", weight="bold",
                 fontsize=11, pad=10)


def panel_1b(ax):
    """umap-explorer: UMAP of 4,700-seq norovirus at chosen state, coloured by Leiden."""
    state_dir = Path("examples/norovirus_vp1/runs/af/states/k7_n100_d0p1_r1e00")
    coords = pd.read_csv(state_dir / "coords.csv")
    assign = pd.read_csv(state_dir / "assignments.tsv", sep="\t")
    df = coords.merge(assign[["id", "cluster", "item_class"]],
                      left_on="ID", right_on="id", how="left")

    # Point colour by cluster
    clusters = sorted(df["cluster"].dropna().unique())
    n_clusters = len(clusters)
    cmap = plt.get_cmap("tab20", n_clusters)
    cluster_to_color = {c: cmap(i) for i, c in enumerate(clusters)}

    genomes = df[df["item_class"] == "genome"]
    refs = df[df["item_class"] == "reference"]

    colors = genomes["cluster"].map(cluster_to_color).tolist()
    ax.scatter(genomes["AF_UMAP1"], genomes["AF_UMAP2"],
               s=8, alpha=0.55, c=colors, edgecolors="none")

    # References overlaid
    if len(refs) > 0:
        ax.scatter(refs["AF_UMAP1"], refs["AF_UMAP2"],
                   s=60, alpha=1.0, c="black", marker="*",
                   edgecolors="white", lw=0.5,
                   label=f"References (n={len(refs)})")

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(
        f"B  umap-explorer — norovirus VP1 (N=4,700), k=7, {n_clusters} Leiden clusters",
        loc="left", weight="bold",
    )
    ax.legend(loc="lower right", frameon=False)


def panel_1c(ax_left, ax_right):
    """distance-explorer preview: PCoA (left) + NJ tree (right) for Chhabra 302."""
    import re
    review = Path("examples/norovirus_vp1/runs/chhabra302/review")

    # PCoA coords + genotype colouring
    pcoa = pd.read_csv(review / "pcoa_coords.csv")
    pcoa = pcoa.rename(columns={"ID": "id", "PCoA1": "PC1", "PCoA2": "PC2"})

    def genogroup(seq_id):
        m = re.match(r"^(GI{1,5}V?|GIX|GX)[.]", seq_id)
        return m.group(1) if m else "other"

    pcoa["genogroup"] = pcoa["id"].map(genogroup)
    groups = sorted(pcoa["genogroup"].unique())
    cmap = plt.get_cmap("tab10", len(groups))
    for i, g in enumerate(groups):
        sub = pcoa[pcoa["genogroup"] == g]
        ax_left.scatter(sub["PC1"], sub["PC2"], s=18, alpha=0.85,
                        color=cmap(i), edgecolors="white", lw=0.3,
                        label=f"{g} (n={len(sub)})")

    ax_left.set_xlabel("PCoA 1")
    ax_left.set_ylabel("PCoA 2")
    ax_left.set_title(
        "C  distance-explorer PCoA — Chhabra 302, coloured by genogroup",
        loc="left", weight="bold",
    )
    ax_left.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0),
                   frameon=False, fontsize=7, ncol=1)

    # Threshold/cluster count curve — the core distance-explorer output
    import json
    with open(review / "validation_metrics.json") as f:
        vm = json.load(f)

    thresholds = sorted(float(t) for t in vm["tree_thresholds"].keys())
    cluster_counts = [vm["tree_thresholds"][f"{t:.6f}"]["cluster_count"]
                      for t in thresholds]
    nmis = [vm["tree_thresholds"][f"{t:.6f}"]["nmi"] for t in thresholds]

    ax_right.plot(thresholds, cluster_counts, "o-", color="#333333", lw=2,
                  label="cluster count")
    ax_right.axhline(53, color="#d62728", ls=":", lw=1.2,
                     label="Chhabra = 53 genotypes")
    ax_right.axvline(0.15, color="#cccccc", ls="--", lw=1)
    ax_right.scatter([0.15], [53], s=140, facecolors="none",
                     edgecolors="red", lw=2, zorder=5)
    ax_right.annotate("slider → 0.15\n→ 53 clusters",
                      xy=(0.15, 53), xytext=(0.17, 90),
                      fontsize=8, color="red",
                      arrowprops=dict(arrowstyle="->", color="red", lw=0.8))

    ax_right.set_xlabel("Agglomerative threshold (slider)")
    ax_right.set_ylabel("Cluster count")
    ax_right.set_title("D  distance-explorer threshold explorer",
                       loc="left", weight="bold")
    ax_right.legend(loc="upper right", frameon=False)


def main():
    fig = plt.figure(figsize=(14, 11))
    gs = fig.add_gridspec(3, 2, width_ratios=[1, 1], height_ratios=[0.9, 1.1, 1],
                          hspace=0.45, wspace=0.35)

    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :])
    ax_c = fig.add_subplot(gs[2, 0])
    ax_d = fig.add_subplot(gs[2, 1])

    panel_1a(ax_a)
    panel_1b(ax_b)
    panel_1c(ax_c, ax_d)

    fig.suptitle(
        "Figure 1. SeqScape staged workflow with interactive HTML widgets "
        "at stages 1 and 4.",
        y=0.995, x=0.01, ha="left", fontsize=11, weight="bold",
    )

    fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    print(f"Saved {OUT_FIG}")
    print(f"Saved {OUT_SVG}")


if __name__ == "__main__":
    main()
