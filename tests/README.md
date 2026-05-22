# Tests

This directory contains lightweight unit tests and a synthetic end-to-end smoke test for the extracted SeqScape package.

The smoke test is the smallest self-contained dataset for trying the software.
It creates temporary FASTA files with four genome-like sequences and two
reference sequences, then runs all four SeqScape stages:

1. `umap-explorer`
2. `sample-panel`
3. `align-panel`
4. `distance-explorer`

Run from the `seqscape/` directory:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

To run only the synthetic end-to-end example:

```bash
PYTHONPATH=src python -m unittest tests.test_pipeline_smoke -v
```

The test data are generated inside `tests/test_pipeline_smoke.py` and written
to a temporary directory during the run. No external FASTA files, NCBI download,
or manuscript example data are required.

Environment notes:
- the end-to-end smoke test needs the AF stack available:
  - `igraph`
  - `leidenalg`
  - `umap-learn`
- `align-panel` smoke uses `pairwisealigner`, not `muscle`, to keep the test fast and self-contained
