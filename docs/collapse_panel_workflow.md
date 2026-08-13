# Forking after `umap-explorer`: two panels, two objectives

`umap-explorer` is the shared entry point. It answers "what is in this sample?"
-- composition, cluster structure, where the focal genomes sit. That first glance
is useful regardless of what you do next, and nothing below changes it.

What follows depends on the question:

```
                    umap-explorer            <- inspect composition, choose a state
                          |
              +-----------+-----------+
              |                       |
        sample-panel            collapse-panel
    (proportional)              (dereplicated)
      "describe the             "screen for
       diversity"                recombination"
              |                       |
              +-----------+-----------+
                          |
                     align-panel
                          |
                   distance-explorer
                          |
                    recombination         <- OpenRDP
```

Both selectors read the **same** chosen explorer state and emit the **same** two
files -- `representative_panel.fasta` and `representative_panel_manifest.tsv` --
so everything downstream accepts either without modification. Verified: the
manifest carries `id`, `label`, `item_class`, `source_cluster`, every manifest id
resolves in the panel FASTA, and `align-panel` consumes it unchanged.

## Why two selectors rather than one flag

The two have opposite objectives, and an objective is not a parameter.

`sample-panel` allocates a fixed budget across Leiden clusters **in proportion to
cluster size** and is scored on `pearson_cluster_proportion_preservation`. It is
built to reproduce the composition of the input. For a descriptive panel that is
exactly right.

For a recombination screen it is exactly wrong, because **the composition is an
artefact**. Public sequence databases record outbreak resequencing effort, not
diversity. On the 1,792-genome ToCSV comparator pull the raw imbalance is 57:1
(1,761 TYLCV-cluster records against 31 ToCSV). A faithful sample faithfully
reproduces that bias.

Measured, running seqscape's own `largest_remainder_alloc` / `farthest_point_order`
/ `select_genome_panel` at budget 173:

| | ToCSV retained | near-duplicate pairs among picks |
|---|---|---|
| `sample-panel` (proportional) | **2 / 31** | 73% |
| `collapse-panel` (threshold)  | **31 / 31** | 0% |

All 31 ToCSV genomes fall in one Leiden cluster (46 genomes: 31 ToCSV + 15
others). Its proportional quota is 46/1792 x 173 = 4.44 -> 4 slots, shared with
the 15 non-ToCSV members. Two survive.

Note this is **not a defect in `sample-panel`**. Farthest-point ordering does
spread picks -- lower within-panel identity than random draws in 21 of 26
sizeable clusters, mean advantage 0.83 pp. The quota is what truncates it, and
the quota is the correct behaviour for the objective `sample-panel` was built for.

## Why dereplication specifically, for recombination

Triplet-based detectors (RDP, GENECONV, MaxChi, Chimaera, BootScan, SiScan) scan
parent-parent-child triplets. Near-identical duplicates add triplets without
adding donor diversity, and the multiple-testing correction applied across those
triplets then costs power at the events that matter. One representative per
distinct lineage is the standard recommendation for RDP-class input
(Martin et al. 2015, *Virus Evolution* 1:vev003).

Concretely, on this pull: 955,889,535 triplets across all 1,791 genomes versus
833,340 across the 171 representatives -- **1,147x fewer tests** for a panel that
retains more of the diversity that matters and all of the focal set.

## Usage

```bash
seqscape collapse-panel \
  --input-fasta genomes.fasta \
  --reference-fasta refs.fasta \
  --umap-explorer-dir out/umap_explorer \
  --chosen-kmer 6 --chosen-neighbors 15 --chosen-min-dist 0.1 \
  --chosen-leiden-resolution 1.0 \
  --identity-threshold 98.1 \
  --focal-ids tocsv_candidates.txt \
  --outdir out/collapse_panel
```

Then continue down the shared path exactly as with `sample-panel`:

```bash
seqscape align-panel \
  --panel-fasta out/collapse_panel/representative_panel.fasta \
  --manifest-tsv out/collapse_panel/representative_panel_manifest.tsv \
  --aligner mafft --outdir out/aligned
```

### Selection rule

Strict priority, no per-cluster budget:

1. `item_class == "reference"` -- always retained (same as `sample-panel`)
2. `--focal-ids` -- always retained, never collapsed, and never allowed to act as
   a representative for another genome
3. everything else -- greedy collapse at `--identity-threshold`, longest genome
   in each group becoming its representative

**Panel size is an output, not an input.** It is set by how many distinct lineages
the data contains. That is the point: a budget is what reintroduces the bias.

### Choosing `--identity-threshold`

Default 98.1% is the ICTV begomovirus **strain** demarcation threshold, so groups
collapse to roughly one representative per strain. Measured sweep, same pull:

| threshold | panel genomes | focal retained | triplets |
|---|---|---|---|
| 95.0% | 60 | 31/31 | 35,990 |
| 97.0% | 110 | 31/31 | 221,815 |
| **98.1%** | **171** | **31/31** | **833,340** |
| 99.0% | 316 | 31/31 | 5,259,030 |
| 99.5% | 632 | 31/31 | 42,072,556 |

The focal set is retained at every threshold -- `--focal-ids` is protected by
construction, so the threshold trades panel size against donor resolution and
nothing else. Lower it if OpenRDP runtime becomes limiting; raise it if you
suspect the true donor is a close relative of a retained representative.

### Identity estimation

Selection uses 15-mer Jaccard converted to identity by the Mash relation
(Ondov et al. 2016, *Genome Biol* 17:132), so it stays alignment-free and cheap.
Exact identities are computed afterwards by `align-panel` on the much smaller
panel. Do not substitute offset/rotation arithmetic for alignment here -- on
circular begomovirus genomes it reports ~28% identity for 98%-identical pairs.

## Outputs beyond the shared two

- `collapse_groups.tsv` -- every representative and the ids it absorbed, so any
  dropped genome can be traced to the one standing in for it
- `source_cluster_retention.tsv` -- per Leiden cluster: size, retained, focal
  count; this is where you check the focal cluster was not thinned
- `validation_metrics.json` -- panel sizes, compression ratio, and both triplet
  counts

## Reproducibility

Leiden's spectral initialisation falls back to random on this data (small
eigengap), so cluster count and per-cluster retention shift somewhat between
seeds. The structural result does not: collapse has no per-cluster quota, so
focal retention is independent of the partition entirely.
