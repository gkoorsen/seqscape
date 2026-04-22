# SeqScape: interactive software for alignment-free sequence-space exploration, representative sampling, and alignment-based review of large compact sequence collections

Gerrit Koorsen  
Department of Biochemistry, University of Johannesburg, Johannesburg, South Africa  
ORCID: 0000-0001-5159-1806  
Correspondence: gkoorsen@uj.ac.za

## Abstract

**Summary:** SeqScape is a Python workflow for exploring large compact-sequence collections in alignment-free space, selecting proportionally representative subsets, and reviewing the resulting structure through interactive HTML visualizations. Four stages chain k-mer UMAP embedding and Leiden clustering to pairwise alignment, PCoA, agglomerative clustering, and reference-based novelty screening. This makes all-versus-all pairwise alignment tractable at scales where direct full-collection analysis is infeasible — routine at current viral sequencing depth. We benchmark the workflow on the norovirus VP1 taxonomy of Chhabra et al. (2019). SeqScape recovers all 53 Chhabra genotype groups from 275 curated VP1 sequences at amino-acid distance 0.15 (NMI = 0.984, ARI = 0.982), and compresses a 4,700-sequence current GenBank collection 88.5-fold while preserving cluster proportions (Pearson r = 0.9998). The resulting panel is directly usable for distance-based phylogenetics, species demarcation, recombination detection, and diversity estimation.

**Availability and implementation:** SeqScape is implemented in Python (>=3.10) and distributed under the MIT license. Source code, documentation, and installation instructions are available at https://github.com/gkoorsen/seqscape. The software is archived at [ZENODO DOI — to be added before submission].

**Keywords:** alignment-free sequence analysis; representative subsampling; interactive visualization; viral taxonomy; UMAP.

## 1 Introduction

High-throughput sequencing has made it routine to assemble collections of tens of thousands of compact homologous sequences, including viral genomes, gene segments, open reading frames, and conserved gene families, from environmental, clinical, and surveillance contexts. Exploratory analysis of such collections faces a fundamental scalability constraint: all-versus-all pairwise alignment scales quadratically with collection size, rendering direct alignment-based distance analysis computationally infeasible at the scales now commonly encountered. A collection of 10,000 sequences requires approximately 50 million pairwise alignments; at 18,000 sequences, this rises to roughly 162 million. Even with aggressive parallelisation, this precludes the iterative, parameter-sensitive exploration that characterises early-stage diversity surveys. Yet pairwise identity matrices remain the currency of many essential downstream analyses, including distance-based phylogenetic reconstruction, ICTV-criteria species demarcation, recombination detection, haplotype deduplication, vaccine and diagnostic representative selection, and population genetic diversity estimation. The challenge is therefore not to avoid pairwise alignment, but to make it tractable by identifying the smallest set of sequences that faithfully represents the diversity of the full collection. At the same time, alignment-free embeddings are exploratory tools that do not on their own provide the interpretive precision of explicit pairwise identity matrices, reference-based classification, or thresholded cluster assignments.

Several existing tools address parts of this problem but none span the full interactive exploration-to-validation arc that SeqScape provides (Table 1). iDeLUCS (Arias et al., 2023) applies alignment-free k-mer UMAP clustering to DNA sequences with an interactive viewer, but is a single-stage classifier and does not produce a proportional subsample for downstream pairwise analysis. TARDiS (Morales-Arce et al., 2023) uses a genetic algorithm to select diverse, temporally distributed subsets for phylogenetic inference, but operates on pre-aligned input and does not provide interactive exploration of alignment-free embeddings. Vclust (Zielezinski et al., 2025) delivers Lempel-Ziv–based ANI clustering for viral genomes at scale but is a batch classifier without subsampling or interactive review. PopPUNK (Lees et al., 2019) uses variable-length k-mer sketches for scalable bacterial lineage assignment within an established reference frame, not for de novo exploration of new collections. vConTACT2 (Bin Jang et al., 2019) builds protein-sharing networks for phage taxonomy and is restricted to that domain. SeqScape addresses this gap by combining interactive alignment-free parameter exploration, proportional representative sampling, and interactive alignment-based review within a single generalizable workflow applicable to arbitrary compact sequence collections. It treats alignment-free embedding as a navigational tool rather than a final result: users explore the sequence landscape interactively, select a stable parameter state for representative sampling, and then validate and interpret a tractably sized panel through a second interactive review bundle built from explicit pairwise alignment distances. We demonstrate the approach on the norovirus VP1 genotyping benchmark of Chhabra et al. (2019), which provides well-defined ground-truth classifications against which SeqScape's output can be quantitatively validated.

**Table 1.** Positioning of SeqScape relative to existing alignment-free and representative-subsampling tools for sequence collection analysis.

| Tool          | Alignment-free exploration | Proportional subsampling | Interactive review widget | Reference-based novelty screening | Domain |
|---------------|:---:|:---:|:---:|:---:|---|
| iDeLUCS       | ✓ | – | – | – | DNA clustering |
| TARDiS        | – | ✓ (diversity + time) | – | – | phylogenetic subsampling |
| Vclust        | ✓ (ANI) | – | – | – | viral genome clustering |
| PopPUNK       | ✓ (sketch) | – | – | ✓ | bacterial lineage assignment |
| vConTACT2     | – | – | – | ✓ | phage taxonomy |
| **SeqScape**  | **✓** | **✓** | **✓** | **✓** | **general compact sequences** |

## 2 Software description

SeqScape is implemented in Python and exposes a four-stage command-line interface. The workflow is designed to be run sequentially, with each stage writing output files consumed by the next.

**Stage 1: Alignment-free sequence-space widget (`umap-explorer`).** Input sequences and an optional reference panel are represented as normalized k-mer frequency vectors. UMAP is applied to embed sequences in two dimensions using Euclidean distance on the k-mer matrix. Leiden clustering is performed on a k-nearest-neighbour graph of the UMAP coordinates. Multiple k-mer sizes, UMAP `n_neighbors` values, and `min_dist` values can be specified in a single run, and the resulting states are embedded into a single interactive HTML widget with slider controls.

**Stage 2: Representative panel sampling (`sample-panel`).** A user-selected parameter state is used to sample a proportionally representative panel from the full genome collection. References are retained. Genome slots are allocated across Leiden clusters using the largest-remainder method to preserve cluster proportions exactly under integer quotas. Within each cluster, selection follows a centroid-first, farthest-point coverage strategy: the sequence closest to the cluster centroid in UMAP space is selected first, and subsequent selections are chosen to maximize coverage of the remaining cluster extent.

**Stage 3: Pairwise alignment panel (`align-panel`).** All-versus-all pairwise global alignment is computed for the sampled panel using Biopython's `PairwiseAligner`, MUSCLE, or MAFFT. Identity is computed as the fraction of matched positions in ungapped aligned columns. A distance matrix is derived as `d = 1 − (identity/100)`. Principal coordinates analysis (PCoA) is computed from the distance matrix.

**Stage 4: Post-alignment review widget (`distance-explorer`).** The review stage builds an interactive HTML bundle from the sampled panel and its alignment-derived matrices. The bundle provides four coordinated views: PCoA of the alignment distances; UMAP embedding on the precomputed distance matrix; a best-versus-second-best reference identity scatter plot with a configurable novelty threshold line; and interactive agglomerative clustering with a complete-linkage threshold slider. Coloring can be toggled between the original Leiden cluster assignments and the agglomerative cluster assignments at any threshold.

## 3 Applications and results

### 3.1 Norovirus VP1 genotype recovery

Norovirus taxonomy provides a well-defined benchmark: Chhabra et al. (2019), corrected by a 2020 corrigendum, established 48 confirmed genotypes across genogroups GI–GX using maximum-likelihood phylogenies (PhyML) of ClustalW-aligned VP1 amino acid sequences, with cluster boundaries defined by a 2×standard-deviation criterion on patristic distance distributions (computed in Patristic and R); GII.4 sub-variants were designated below the genotype rank on the basis of phylogenetic clustering combined with documented epidemic circulation in at least two geographically diverse locations. We applied SeqScape to test (a) whether the distance-explorer can recover this taxonomy when working from the curated reference set, and (b) whether the staged workflow scales to the current depth of GenBank norovirus VP1 deposits.

**Validation on the Chhabra reference set.** We extracted VP1 coding sequences from the 302 accessions listed in Chhabra et al. Supplementary Table S1, translating to amino acids, yielding 275 sequences representing 53 distinct genotype groups (48 confirmed plus five tentative assignments present in the supplementary data). Twenty-seven accessions corresponded to complete genomes with non-standard CDS annotations and were recovered by searching for VP1 product keywords in the GenBank feature table; 26 accessions encoding other ORFs or lacking sufficient annotation were excluded. Pairwise MUSCLE alignment of the 275-sequence set required 37,675 comparisons and completed in ~5 min on 8 CPU threads. The distance-explorer NJ tree at an amino acid distance threshold of 0.15 partitioned the 275 sequences into exactly 53 clusters, matching the 53 Chhabra genotype groups with NMI = 0.984 and ARI = 0.982 (Figure 2A). A biologist navigating the distance-explorer would locate this threshold by observing the plateau in cluster count between 0.12 and 0.20 (53–61 clusters), which coincides with the known genotype-level divergence in norovirus VP1: the distance-explorer provides the threshold explorer; the biologist selects the threshold consistent with biological knowledge.

**Scalability on the full GenBank collection.** We retrieved all GenBank norovirus records annotated as VP1 or extractable from complete genome sequences (2,500–10,000 nt; 400–1,100 nt for direct VP1 deposits), deduplicated by sequence hash, and retained 4,700 sequences. Alignment-free UMAP embedding at k=7, `n_neighbors=100`, `min_dist=0.1`, Leiden resolution 1.0 resolved 21 clusters. A representative panel of 500 sequences (10.6% of the collection; 13 CDC typing-tool references retained in full) was sampled proportionally across Leiden clusters, preserving cluster proportions almost exactly (Pearson r = 0.9998, p < 10⁻³³; Figure 2B). The panel reduces the alignment burden from 11,042,650 to 124,750 comparisons — an 88.5-fold compression.

MUSCLE alignment of the 500-sequence panel completed in ~14 min on 8 CPU threads. PCoA resolved 32.7% and 18.5% of variance on axes 1 and 2. At the Leiden-derived agglomerative threshold (0.10), 59 clusters were recovered; against Chhabra genotype labels this yields NMI = 0.876 and ARI = 0.795. The lower concordance relative to the validation benchmark reflects genuine biological signal: GII.4, comprising 46.6% of the panel (228/489 typed sequences), is split into 10 sub-clusters corresponding to recognised GII.4 variants (Sydney 2012, DenHaag 2006, New Orleans 2009, etc.) that Chhabra's single "GII.4" label collapses. This sub-structure is visible in the distance-explorer PCoA and is interpretable by a biologist familiar with GII.4 epidemiology (Figure 2C). Among the 37 genotypes represented in the panel (of 48 confirmed), all had ≥10 sequences in the full dataset — above the minimum lineage size ⌈4700/500⌉ = 10 required for reliable proportional capture. The 11 confirmed genotypes absent from the panel each had fewer than 10 sequences in the current GenBank VP1 collection, consistent with this threshold.

A total of 204 panel sequences (41.9% of genomes) fell below 92% amino acid identity to the nearest CDC reference, flagged by the distance-explorer novelty indicator as representing lineages not well-anchored by the current reference set.

![**Figure 1.** SeqScape staged workflow with interactive HTML widgets at stages 1 and 4.](figures/figure1.png)

![**Figure 2.** SeqScape quantitative validation on the norovirus VP1 benchmark. (A) Chhabra 302 cluster count, NMI, and ARI vs amino acid distance threshold. (B) Panel cluster proportion vs full-collection cluster proportion (Pearson r = 0.9998). (C) GII.4 dendrogram at threshold 0.10.](figures/figure2.png)

### 3.2 Workflow reproducibility

Complete workflows including dataset retrieval, filtering, reference selection, and SeqScape CLI commands are provided in the repository examples at https://github.com/gkoorsen/seqscape/tree/main/examples.

## 4 Conclusion

SeqScape provides a practical, general-purpose workflow for navigating large collections of compact homologous sequences through interactive alignment-free exploration, proportional representative sampling, and alignment-based validation. The staged design makes explicit pairwise alignment tractable at scales where full all-versus-all comparison is computationally infeasible. On the norovirus VP1 benchmark, 11,042,650 theoretical pairwise alignments in the 4,700-sequence collection were reduced to 124,750 within a 500-sequence representative panel — an 88.5-fold compression — while preserving the cluster diversity structure of the full collection (Pearson r = 0.9998). Applied directly to the curated Chhabra VP1 reference set (275 sequences representing 53 genotype groups), the distance-explorer recovers all 53 groups with NMI = 0.984 and ARI = 0.982 at an amino acid distance threshold of 0.15, demonstrating that the tool replicates expert taxonomy when the reference set is well-curated; on the 500-sequence representative panel drawn from 4,700 GenBank VP1 sequences, agreement against Chhabra labels is necessarily lower (NMI = 0.876, ARI = 0.795) because the panel resolves biologically real GII.4 sub-variant structure that the reference taxonomy collapses. The representative panel produced by SeqScape is directly suitable as input to downstream analyses including distance-based phylogenetic reconstruction, ICTV-criteria species demarcation, recombination detection, population genetic diversity estimation, and representative selection for vaccines or diagnostics. The interactive HTML outputs are self-contained and require no server infrastructure, making them suitable for sharing as supplementary files or for collaborative exploratory analysis.

## Funding

[Add funding statement here if applicable, or state: This work received no specific grant funding.]

## References

1. Arias, P.M., Butorac, A., Suzek, B.E., et al. (2023) iDeLUCS: a deep learning interactive tool for alignment-free clustering of DNA sequences. *Bioinformatics*, **39**, btad508.
2. Bin Jang, H., Bolduc, B., Zablocki, O., et al. (2019) Taxonomic assignment of uncultivated prokaryotic virus genomes is enabled by gene-sharing networks. *Nat. Biotechnol.*, **37**, 632–639.
3. Chhabra, P., de Graaf, M., Parra, G.I., et al. (2019) Updated classification of norovirus genogroups and genotypes. *J. Gen. Virol.*, **100**, 1393–1406. [Corrigendum: (2020) *J. Gen. Virol.*, **101**, 893.]
4. Cock, P.J.A., Antao, T., Chang, J.T., et al. (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, **25**, 1422–1423.
5. Edgar, R.C. (2022) MUSCLE v5 enables improved estimates of phylogenetic tree confidence by ensemble bootstrapping. *Nat. Commun.*, **13**, 6968.
6. Katoh, K. and Standley, D.M. (2013) MAFFT multiple sequence alignment software version 7. *Mol. Biol. Evol.*, **30**, 772–780.
7. Lees, J.A., Harris, S.R., Tonkin-Hill, G., et al. (2019) Fast and flexible bacterial genomic epidemiology with PopPUNK. *Genome Res.*, **29**, 304–316.
8. McInnes, L., Healy, J., and Melville, J. (2018) UMAP: Uniform Manifold Approximation and Projection for dimension reduction. arXiv:1802.03426.
9. Morales-Arce, A.Y., Johri, P., and Jensen, J.D. (2023) TARDiS: a tool for diversity- and time-aware subsampling of genomic datasets. *Bioinformatics*, **39**, btad352.
10. Saitou, N. and Nei, M. (1987) The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Mol. Biol. Evol.*, **4**, 406–425.
11. Sayers, E.W., Bolton, E.E., Brister, J.R., et al. (2022) Database resources of the National Center for Biotechnology Information. *Nucleic Acids Res.*, **50**, D20–D26.
12. Traag, V.A., Waltman, L., and van Eck, N.J. (2019) From Louvain to Leiden: guaranteeing well-connected communities. *Sci. Rep.*, **9**, 5233.
13. Vinh, N.X., Epps, J., and Bailey, J. (2010) Information theoretic measures for clusterings comparison. *J. Mach. Learn. Res.*, **11**, 2837–2854.
14. Zielezinski, A., Gudyś, A., Barylski, J., et al. (2025) Ultrafast and accurate sequence alignment and clustering of viral genomes. *Nat. Methods*, **22**, 654–663.
