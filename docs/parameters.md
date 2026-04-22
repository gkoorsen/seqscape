# Parameters

Reference for every CLI flag exposed by `seqscape`. Defaults and
semantics are taken directly from `src/seqscape/cli.py`; use
`seqscape {stage} --help` for the same information at the terminal.

Invocation forms:

```bash
# installed console script
seqscape {stage} ...

# from a checkout
PYTHONPATH=src python -m seqscape.cli {stage} ...
```

Shared conventions:

- Comma-separated lists (e.g. `--kmer-values 5,7,9`) declare grids of
  parameter values. Do **not** include spaces inside the list.
- `min_dist` float values are encoded in the on-disk state slug by
  replacing the decimal point with `p` (e.g. `0.3 → d0p3`).
- Leiden resolution is encoded in scientific notation (e.g.
  `1e-06 → r1e-06`).

---

## `seqscape umap-explorer`

Build the alignment-free UMAP + Leiden parameter explorer.

### Inputs

| Flag | Default | Description |
|---|---|---|
| `--input-fasta` | *required* | Nucleotide FASTA of genomes to embed |
| `--reference-fasta` | *required* | Nucleotide FASTA of reference anchors (rendered as stars) |
| `--cg-map-tsv` | `""` | Optional 2-column TSV mapping sequence IDs to human-readable labels |

### Alignment-free embedding grid

| Flag | Default | Description |
|---|---|---|
| `--kmer-values` | `5` | Comma-separated k-mer lengths to sweep |
| `--neighbors-values` | `50,100,150,200` | Comma-separated UMAP `n_neighbors` values |
| `--min-dist-values` | `0.05,0.1,0.3,0.5,0.8` | Comma-separated UMAP `min_dist` values |
| `--seed` | `42` | Random seed for UMAP |

### Leiden clustering

| Flag | Default | Description |
|---|---|---|
| `--leiden-neighbor-k` | `300` | k for the shared-nearest-neighbour graph Leiden runs on |
| `--leiden-resolution` | `1e-6` | Single Leiden resolution (used when `--leiden-resolution-values` is empty) |
| `--leiden-resolution-values` | `""` | Comma-separated resolutions to sweep (overrides the single value) |

### Default state shown on widget open

| Flag | Default | Description |
|---|---|---|
| `--default-kmer` | `5` | k-mer of the initially displayed state |
| `--default-neighbors` | `100` | `n_neighbors` of the initially displayed state |
| `--default-min-dist` | `0.3` | `min_dist` of the initially displayed state |
| `--default-leiden-resolution` | unset | Resolution of the initially displayed state |

### Trustworthiness validation

| Flag | Default | Description |
|---|---|---|
| `--validation-subsample-size` | `1000` | Sample size for trustworthiness computation |
| `--validation-replicates` | `1` | Number of replicates (averaged) |
| `--validation-neighbors` | `10,30,100` | Neighbourhood sizes at which to report trustworthiness |

### Output

| Flag | Default | Description |
|---|---|---|
| `--outdir` | *required* | Output directory; receives `umap_explorer.html`, `summary.txt`, and `states/{slug}/` subdirectories |

---

## `seqscape sample-panel`

Sample a proportional representative panel from a chosen UMAP-Explorer
state.

| Flag | Default | Description |
|---|---|---|
| `--input-fasta` | *required* | Same genome FASTA used upstream |
| `--reference-fasta` | *required* | Same reference FASTA used upstream; all references are retained in the panel |
| `--umap-explorer-dir` | *required* | Directory produced by `seqscape umap-explorer` |
| `--chosen-kmer` | *required* | k-mer of the selected state |
| `--chosen-neighbors` | *required* | `n_neighbors` of the selected state |
| `--chosen-min-dist` | *required* | `min_dist` of the selected state |
| `--chosen-leiden-resolution` | `1e-6` | Leiden resolution of the selected state |
| `--sample-size` | *required* | Target total panel size (references + genomes) |
| `--outdir` | *required* | Output directory; receives `representative_panel.fasta`, `representative_panel_manifest.tsv`, `source_cluster_quotas.tsv`, `summary.txt`, `validation_metrics.json` |

The genome budget is `sample-size − references_retained`. Sampling is
proportional to Leiden cluster size in the selected state.

---

## `seqscape align-panel`

Run all-versus-all pairwise alignment on the panel.

| Flag | Default | Description |
|---|---|---|
| `--panel-fasta` | *required* | Panel FASTA from `sample-panel` |
| `--manifest-tsv` | *required* | Panel manifest from `sample-panel` |
| `--aligner` | `muscle` | `muscle`, `mafft`, or `pairwisealigner` (Biopython) |
| `--mafft` | `/opt/homebrew/bin/mafft` | Path to the MAFFT binary |
| `--muscle` | `muscle` | Path or command name for the MUSCLE binary |
| `--jobs` | `8` | Parallel alignment workers |
| `--progress-every` | `500` | Log progress every N completed pairs |
| `--outdir` | *required* | Output directory; receives `matrix_identity.tsv`, `matrix_distance.tsv`, `pcoa_coords.csv`, and a copy of the panel FASTA and manifest |

Identity matrix values are percent identity (0–100). Distance matrix
values are `1 − identity/100` (range 0–1) and feed directly into the
distance-explorer clustering.

---

## `seqscape distance-explorer`

Build the post-alignment distance explorer.

### Inputs

| Flag | Default | Description |
|---|---|---|
| `--panel-fasta` | *required* | Panel FASTA from `align-panel` |
| `--identity-matrix` | *required* | `matrix_identity.tsv` from `align-panel` |
| `--distance-matrix` | *required* | `matrix_distance.tsv` from `align-panel` |
| `--manifest-tsv` | *required* | Panel manifest |
| `--umap-explorer-dir` | `""` | Optional upstream UMAP-Explorer directory; enables the `AF UMAP` view and the cluster-agreement table |

### Clustering and novelty

| Flag | Default | Description |
|---|---|---|
| `--thresholds` | `0.03,0.04,0.05,0.06,0.07,0.08,0.10` | Distance thresholds at which to precompute agglomerative and tree partitions |
| `--novel-threshold` | `92.0` | Best-vs-second-best reference identity (%) below which a panel genome is flagged as candidate novel |

### Distance-UMAP settings

| Flag | Default | Description |
|---|---|---|
| `--umap-neighbors` | `30` | `n_neighbors` for the distance-UMAP panel |
| `--umap-min-dist` | `0.5` | `min_dist` for the distance-UMAP panel |
| `--umap-spread` | `0.8` | `spread` for the distance-UMAP panel |
| `--umap-seed` | `42` | Random seed |

### Output

| Flag | Default | Description |
|---|---|---|
| `--outdir` | *required* | Output directory; receives `distance_explorer.html` and the tabular outputs listed in [`widgets.md`](widgets.md#tabular-outputs-written-alongside-the-widget) |

---

## `seqscape export-clusters`

Export cluster members to per-cluster FASTA files. Not a required
workflow stage — a reusable utility driven by any cluster-assignment
TSV produced by the pipeline (UMAP-Explorer `assignments.tsv`, panel
`representative_panel_manifest.tsv`, or distance-explorer
`agglomerative_threshold_assignments.tsv` / `tree_threshold_assignments.tsv`).

### Inputs

| Flag | Default | Description |
|---|---|---|
| `--sequence-fasta` | *required* | FASTA containing the sequences to export |
| `--reference-fasta` | `""` | Optional additional FASTA of reference sequences to merge in |
| `--cluster-tsv` | *required* | TSV with one row per sequence and a cluster-label column |
| `--manifest-tsv` | `""` | Optional manifest used when cluster IDs need to be resolved to panel IDs |
| `--id-column` | `id` | Name of the sequence-ID column in `--cluster-tsv` |
| `--cluster-column` | `cluster` | Name of the cluster-label column in `--cluster-tsv` |

### Filtering

| Flag | Default | Description |
|---|---|---|
| `--filter-column` | `""` | Optional column to filter rows on (e.g. `threshold`) |
| `--filter-value` | `""` | Required value in `--filter-column` (as printed in the TSV — typically 6-decimal float, e.g. `0.050000`) |
| `--cluster-list` | `""` | Optional comma-separated whitelist of cluster names to export (if empty, exports all) |

### Output

| Flag | Default | Description |
|---|---|---|
| `--filename-prefix` | `""` | Prefix prepended to every exported FASTA filename |
| `--outdir` | *required* | Output directory; receives one FASTA per cluster |

See [`workflow.md`](workflow.md) for three common usage patterns
(full-set Leiden, panel Leiden, agglomerative at a chosen threshold).
