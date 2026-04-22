# Widgets

SeqScape produces two interactive, self-contained HTML widgets — one per
interactive stage. Both are single `.html` files that require no server,
render in any modern browser, and are suitable for sharing as supplementary
files.

---

## 1. UMAP Explorer (`umap_explorer.html`)

Written to `{outdir}/umap_explorer.html` by `seqscape umap-explorer`.

The widget scans a grid of alignment-free parameter states (k-mer length,
UMAP `n_neighbors`, UMAP `min_dist`, Leiden `resolution`) and lets the user
switch between states interactively. Each state is rendered from
precomputed coordinates and cluster assignments stored in
`{outdir}/states/{state_slug}/`.

### Panels

- **Main UMAP scatter.** Genomes shown as filled circles; user-supplied
  reference sequences shown as stars. Points coloured by Leiden cluster.
  Convex hulls drawn around each cluster.
- **Legend.** Cluster labels (C01, C02, …) with genome counts and
  percentage of the total dataset.
- **Metrics strip** (eight cards above the scatter):
  - `Genomes` — number of genome sequences embedded
  - `References` — number of references embedded
  - `Leiden Clusters` — cluster count for the currently selected state
  - `Largest Cluster` — size of the largest cluster
  - `Trustworthiness n=10` and `Trustworthiness n=30` — local embedding
    quality at two neighbourhood sizes
  - `Small Clusters` — number of clusters with fewer than 10 genomes
  - `Validation Sample` — subsample size used for trustworthiness

### Interactive controls

- **Parameter sliders** (four range sliders): `k-mer`, `n_neighbors`,
  `min_dist`, `Leiden resolution`. Each slider snaps to the grid values
  supplied to the CLI. Moving any slider swaps to the corresponding
  precomputed state.
- **Export Clusters dropdown.** A multi-select checkbox dropdown labelled
  `Select clusters ▾`. Contains an `All clusters` master checkbox plus
  one checkbox per Leiden cluster.
- **Sample Size** (number input, default `500`). Sets the genome budget
  used by the generated `sample-panel` command.
- **Export buttons:**
  - `Download IDs` — writes `{state_slug}_{cluster_names}_ids.txt`, a
    newline-separated list of the IDs in the selected clusters.
  - `Download FASTA Command` — writes a bash script that calls
    `seqscape export-clusters` with the correct arguments to extract
    those clusters as per-cluster FASTA files.
  - `Download Next-Step Command` — writes a bash script that calls
    `seqscape sample-panel` with the currently selected parameter state
    and sample size.

### Files written per state

`{outdir}/states/{state_slug}/` contains:

- `assignments.tsv` — columns `id`, `label`, `item_class` (`genome` or
  `reference`), `cluster`
- `cluster_counts.tsv` — cluster name, size, percentage
- `coords.csv` — UMAP x/y coordinates
- `validation_metrics.json` — trustworthiness and cluster diagnostics

State slugs encode the parameter quadruple:
`k{k}_n{neighbors}_d{min_dist_p}_r{resolution}`, e.g. `k5_n100_d0p3_r1e00`.

### How to read it

1. Start with the default state and scan the metric strip. If
   trustworthiness is low (<0.9) or the cluster count looks implausibly
   large, try coarser `n_neighbors` and larger `min_dist`.
2. Sweep the `Leiden resolution` slider until the cluster count
   stabilises — the plateau is a good operating point.
3. Visually check the main scatter: well-separated clusters with intact
   reference stars near same-clade genomes indicate a sensible
   embedding.
4. Select the state that gives the cleanest combination of
   trustworthiness and interpretable cluster structure, then use
   `Download Next-Step Command` to produce the panel sampling command.

---

## 2. Distance Explorer (`distance_explorer.html`)

Written to `{outdir}/distance_explorer.html` by `seqscape distance-explorer`.

Built after pairwise alignment of a representative panel, this widget
shows the alignment-based geometry of the panel and compares distance-based
cluster partitions against the alignment-free Leiden partitions carried
over from the UMAP Explorer.

### Panels (switchable via view-mode buttons)

- **PCoA** — Principal coordinates projection of the panel distance
  matrix. Axes show % variance explained.
- **Distance UMAP** — UMAP embedding computed directly on the alignment
  distance matrix.
- **AF UMAP** — the alignment-free UMAP carried over from the UMAP
  Explorer (requires `--umap-explorer-dir`).
- **Best vs 2nd** — scatter of each panel genome's best reference
  identity (x) against its second-best reference identity (y). A flagged
  novelty region highlights candidate lineages without a close reference
  match.
- **NJ Tree** — canvas-rendered neighbor-joining tree of the panel
  distance matrix.
- **Cluster Agreement** — table of Normalised Mutual Information between
  every UMAP-Explorer state and every distance-derived partition
  (agglomerative and tree-based). Cells are clickable: clicking a cell
  switches the colour scheme to that AF state and sets the distance
  threshold accordingly.

### Interactive controls

- **Colour scheme** (three buttons):
  - `Color by Leiden` — colour by AF/Leiden cluster from the UMAP
    Explorer
  - `Color by Agglomerative` — colour by complete-linkage agglomerative
    cluster at the current threshold
  - `Color by Tree` — colour by NJ-tree distance cluster at the current
    threshold
- **Complete-linkage threshold** (range slider) — controls the
  agglomerative cut-off. Snaps to the discrete values supplied via
  `--thresholds`.
- **NJ tree-distance threshold** (range slider) — controls the tree-cut
  threshold.
- **AF state selectors** (three sliders: `k-mer`, `n_neighbors`,
  `min_dist`) — switch among alignment-free parameter states carried
  over from the UMAP Explorer.
- **Best-vs-2nd best cutoff (%)** (number input, default `92.0`) — sets
  the novelty threshold. Any panel genome whose best reference identity
  is below this value is flagged as a candidate novel lineage in the
  `Best vs 2nd` view and in exports.
- **Reference self-handling** (two buttons):
  - `Refs include self` — plot references at (100%, best_identity)
  - `Refs exclude self` — plot references against their best
    non-self match
- **Scale** (`Full scale` / `Robust scale`) — axis scaling for the
  current plot.
- **Length filter (nt)** (min/max inputs) — hide panel members outside a
  nucleotide-length window.
- **Export Clusters dropdown.** Same pattern as the UMAP Explorer:
  multi-select checkboxes over the currently coloured clusters with a
  master `All clusters` toggle.
- **Export buttons:**
  - `Download IDs` — writes `{tag}_{cluster_names}_ids.txt`
  - `Download FASTA Command` — writes a bash script calling
    `seqscape export-clusters` with the correct `--cluster-tsv`,
    `--cluster-column`, and `--filter-column`/`--filter-value` for the
    currently selected colour scheme and threshold.
- **Agreement panel mode** (shown when the agreement table is open):
  `Against Agglomerative` or `Against Tree` — switches the reference
  partition the NMI is computed against.

### Tabular outputs written alongside the widget

- `agglomerative_threshold_assignments.tsv` — columns `id`, `threshold`,
  `agglomerative_cluster`. One row per genome per threshold.
- `tree_threshold_assignments.tsv` — same layout with `tree_cluster`.
- `agglomerative_threshold_summary.tsv` and
  `tree_threshold_summary.tsv` — cluster counts and size distributions
  at each threshold.
- `cluster_agreement_matrix.tsv` — NMI values with columns `target`,
  `af_state`, `kmer`, `neighbors`, `min_dist`, `leiden_resolution`,
  `af_cluster_count`, `threshold`, `target_cluster_count`, `nmi`.
- `best_reference_summary.tsv` — per-genome best and second-best
  reference identities used by the `Best vs 2nd` view.
- `genome_x_reference_identity.tsv` — full genome-by-reference identity
  matrix.
- `panel_square_identity.tsv` — panel-by-panel identity matrix.
- `pcoa_coords.csv`, `distance_umap_coords.csv` — coordinates for the
  PCoA and distance-UMAP panels.
- `neighbor_joining_tree.nwk` and `neighbor_joining_tree_tip_map.tsv`
  — NJ tree in Newick format and the mapping from panel IDs to tree
  tip labels.
- `heatmap_genome_x_reference.{html,png}` and
  `heatmap_panel_square.{html,png}` — static heatmap renderings.

### How to read it

1. Start in **PCoA** with `Color by Leiden`. If AF/Leiden clusters
   separate cleanly in PCoA, the alignment-free embedding is
   consistent with the alignment geometry.
2. Switch to **Cluster Agreement** and scan the NMI table. The
   highest-NMI cell identifies the distance threshold that best
   matches an AF state. Clicking it applies that configuration.
3. Use **NJ Tree** with `Color by Tree` to inspect tree-cut partitions
   directly; the tree-distance threshold slider exposes the plateau.
4. In **Best vs 2nd**, watch for genomes far below the novelty
   cutoff — these are candidate new lineages worth flagging for
   independent review.
5. For any final cluster set of interest, use the export dropdown to
   pull per-cluster FASTA via `export-clusters`.
