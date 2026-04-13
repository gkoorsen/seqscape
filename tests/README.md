# Tests

This directory contains lightweight unit tests and a synthetic end-to-end smoke test for the extracted SeqScape package.

Run from the `seqscape/` directory:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

Environment notes:
- the end-to-end smoke test needs the AF stack available:
  - `igraph`
  - `leidenalg`
  - `umap-learn`
- `align-panel` smoke uses `pairwisealigner`, not `muscle`, to keep the test fast and self-contained
