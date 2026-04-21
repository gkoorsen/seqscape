# Quickstart

This is the minimal staged workflow for SeqScape.

## 1. Build the UMAP Explorer

```bash
cd seqscape
PYTHONPATH=src python -m seqscape.cli umap-explorer \
  --input-fasta path/to/phased_filtered.fasta \
  --reference-fasta path/to/reference_panel.fasta \
  --kmer-values 5 \
  --neighbors-values 100 \
  --min-dist-values 0.3 \
  --outdir runs/example_af
```

## 2. Sample a representative panel

```bash
cd seqscape
PYTHONPATH=src python -m seqscape.cli sample-panel \
  --input-fasta path/to/phased_filtered.fasta \
  --reference-fasta path/to/reference_panel.fasta \
  --umap-explorer-dir runs/example_af \
  --chosen-kmer 5 \
  --chosen-neighbors 100 \
  --chosen-min-dist 0.3 \
  --sample-size 500 \
  --outdir runs/example_panel
```

## 3. Align the panel

```bash
cd seqscape
PYTHONPATH=src python -m seqscape.cli align-panel \
  --panel-fasta runs/example_panel/representative_panel.fasta \
  --manifest-tsv runs/example_panel/representative_panel_manifest.tsv \
  --aligner muscle \
  --jobs 8 \
  --outdir runs/example_aligned
```

## 4. Build the Distance Explorer

```bash
cd seqscape
PYTHONPATH=src python -m seqscape.cli distance-explorer \
  --panel-fasta runs/example_aligned/representative_panel.fasta \
  --identity-matrix runs/example_aligned/matrix_identity.tsv \
  --distance-matrix runs/example_aligned/matrix_distance.tsv \
  --manifest-tsv runs/example_aligned/representative_panel_manifest.tsv \
  --thresholds 0.05,0.06 \
  --outdir runs/example_review
```

## Notes

- `umap-explorer` and `distance-explorer` require the UMAP/Leiden stack
- `align-panel` can run with `muscle`, `mafft`, or `pairwisealigner`
- the repository test suite includes a synthetic end-to-end smoke run
