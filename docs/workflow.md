# Workflow

SeqScape exposes a four-stage workflow:

1. `af-widget`
2. `sample-panel`
3. `align-panel`
4. `review-widget`

Auxiliary command:

5. `export-clusters`

`export-clusters` is not a required stage in the main pipeline. It is a reusable export utility that can write one FASTA per cluster from:

- a full-set AF Leiden state `assignments.tsv`
- a sampled-panel `representative_panel_manifest.tsv` using `source_cluster`
- a review-stage `agglomerative_threshold_assignments.tsv` filtered to a chosen threshold

Typical usage patterns:

- Full AF Leiden clusters:
```bash
PYTHONPATH=src python -m seqscape.cli export-clusters \
  --sequence-fasta genomes.fasta \
  --reference-fasta refs.fasta \
  --cluster-tsv af_widget/states/k5_n100_d0p3/assignments.tsv \
  --cluster-column cluster \
  --filename-prefix leiden_full \
  --outdir exports/full_leiden
```

- Sampled-panel Leiden clusters:
```bash
PYTHONPATH=src python -m seqscape.cli export-clusters \
  --sequence-fasta panel/representative_panel.fasta \
  --cluster-tsv panel/representative_panel_manifest.tsv \
  --manifest-tsv panel/representative_panel_manifest.tsv \
  --cluster-column source_cluster \
  --filename-prefix leiden_panel \
  --outdir exports/panel_leiden
```

- Agglomerative clusters at a chosen threshold:
```bash
PYTHONPATH=src python -m seqscape.cli export-clusters \
  --sequence-fasta panel/representative_panel.fasta \
  --cluster-tsv review/agglomerative_threshold_assignments.tsv \
  --manifest-tsv panel/representative_panel_manifest.tsv \
  --cluster-column agglomerative_cluster \
  --filter-column threshold \
  --filter-value 0.050000 \
  --filename-prefix agg005 \
  --outdir exports/agglomerative_0p05
```

Canonical entry points:

- package CLI:
  - `PYTHONPATH=src python -m seqscape.cli ...`
- installed console script:
  - `seqscape ...`

Repository compatibility entry point:

- from the legacy repo root:
  - `python sequence_space_pipeline.py ...`

The root-level `sequence_space_pipeline.py` is now only a thin wrapper around the package CLI. The package implementation in `src/seqscape/` is the authoritative code path.
