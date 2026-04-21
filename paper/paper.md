---
title: "SeqScape: sequence-space exploration, representative sampling, and alignment-based review"
tags:
  - Python
  - bioinformatics
  - sequence space
  - clustering
  - UMAP
  - representative sampling
authors:
  - name: Gerrit Koorsen
    orcid: "https://orcid.org/0000-0001-5159-1806"
    affiliation: 1
affiliations:
  - name: "Department of Biochemistry, University of Johannesburg, Johannesburg, South Africa"
    index: 1
date: 2026-04-13
bibliography: paper.bib
---

# Summary

SeqScape is a staged software workflow for exploring large sequence collections in alignment-free sequence space, selecting representative subsets for tractable downstream analysis, and reviewing the resulting alignment-based structure in interactive visualizations. The software provides four connected stages: `umap-explorer`, `sample-panel`, `align-panel`, and `distance-explorer`. Together, these stages support a pragmatic analysis path in which large sequence collections are first organized in alignment-free space, then sampled proportionally, then validated and interpreted using pairwise identity matrices, ordination, agglomerative clustering, and reference-space summaries. The workflow combines UMAP-based embedding [@mcinnes2018umap], Leiden clustering [@traag2019leiden], pairwise sequence identity estimation using Biopython [@cock2009biopython] or external aligners such as MUSCLE [@edgar2004muscle] and MAFFT [@katoh2013mafft], and interactive HTML review outputs. The workflow was developed in the context of viral genomics, but is intentionally general enough to support other compact sequence collections such as genes, ORFs, or similar homologous sequence sets.

# Statement of need

Large sequence datasets are increasingly common, but direct all-versus-all alignment and distance analysis across full collections is often computationally expensive and operationally awkward. This is especially visible in viral genomics, where compact genomes and ORFs are often available in large numbers and exploratory clustering must be balanced against tractable downstream validation. At the same time, purely alignment-free embeddings are useful for exploration but do not on their own provide the same level of interpretability as explicit alignment-based distances, reference comparisons, or thresholded cluster assignments.

SeqScape addresses this gap by connecting alignment-free exploration to alignment-based follow-up in a single staged workflow. It is designed for users who need to:

- inspect large sequence collections in an interactive embedding
- compare clustering behavior across parameter settings
- select representative panels from the full collection instead of aligning every sequence
- validate the sampled panel with pairwise identity or distance calculations
- review the sampled panel in ordination, thresholded agglomerative clustering, and reference-space views

The software is especially well suited to compact genomes, genes, or ORFs where large sequence collections can be meaningfully explored in k-mer space and then reviewed in a smaller alignment-based panel. Although the motivating use case was viral genome and ORF analysis, the workflow is not limited to viral sequences.

# Core functionality

SeqScape provides four connected stages.

## 1. Alignment-free sequence-space widget

The `umap-explorer` stage builds an interactive alignment-free parameter explorer over genomes plus a reference panel. Sequences are represented in k-mer space, embedded with UMAP [@mcinnes2018umap], and clustered with Leiden [@traag2019leiden]. Multiple parameter states can be precomputed and compared through a single HTML widget.

## 2. Representative panel sampling

The `sample-panel` stage selects a reduced representative panel from a chosen alignment-free state. References are retained explicitly, while genome slots are allocated proportionally across Leiden clusters. Within each cluster, sampling prioritizes a centroid-like sequence and then expands coverage toward more peripheral points.

## 3. Pairwise alignment panel

The `align-panel` stage computes pairwise identity and distance matrices for the sampled panel using Biopython pairwise alignment [@cock2009biopython] or external aligners such as MUSCLE [@edgar2004muscle] or MAFFT [@katoh2013mafft]. This stage writes reusable identity, distance, and ordination outputs for downstream review.

## 4. Post-alignment review widget

The `distance-explorer` stage builds an interactive review bundle from the sampled panel and its pairwise identity or distance matrices. Views include PCoA, alignment-distance UMAP, agglomerative clustering with threshold controls, best-versus-second-best reference identity, and heatmaps for genome-versus-reference and full panel similarity.

# Availability

SeqScape is implemented in Python and distributed as an open-source package. Source code is available at:

`https://github.com/gkoorsen/seqscape`

# Current manuscript notes

Before submission, the following still need to be finalized in this paper draft:

- release DOI after tagging and Zenodo archiving
- any additional domain-specific citations that should be included in `paper.bib`
