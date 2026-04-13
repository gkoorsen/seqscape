# Workflow

SeqScape exposes a four-stage workflow:

1. `af-widget`
2. `sample-panel`
3. `align-panel`
4. `review-widget`

Canonical entry points:

- package CLI:
  - `PYTHONPATH=src python -m seqscape.cli ...`
- installed console script:
  - `seqscape ...`

Repository compatibility entry point:

- from the legacy repo root:
  - `python sequence_space_pipeline.py ...`

The root-level `sequence_space_pipeline.py` is now only a thin wrapper around the package CLI. The package implementation in `src/seqscape/` is the authoritative code path.
