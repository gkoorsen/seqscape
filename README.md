# SeqScape

`SeqScape` is a staged workflow for sequence-space exploration, representative sampling, and alignment-based review.

It is being extracted from the larger exploratory repository into a smaller publication-facing software package for the planned JOSS and BMC Bioinformatics papers.

Repository:
- `https://github.com/gkoorsen/seqscape`

## Current Status

- package CLI implemented in [src/seqscape](/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/seqscape/src/seqscape)
- four pipeline stages available:
  - `umap-explorer`
  - `sample-panel`
  - `align-panel`
  - `distance-explorer`
- auxiliary export command available:
  - `export-clusters`
- package test suite present
- synthetic end-to-end smoke test passing
- under active development

## Installation

Recommended environment:

```bash
cd seqscape
conda env create -f environment.yml
conda activate seqscape
```

The `align-panel` stage can use `muscle`, `mafft`, or Biopython `pairwisealigner`.

External tools:
- `MUSCLE` is recommended for production pairwise alignment panels
- `MAFFT` is optional if you want the `mafft` alignment mode

## CLI

Package entry point:

```bash
cd seqscape
PYTHONPATH=src python -m seqscape.cli --help
```

Compatibility entry point from the parent repo:

```bash
python sequence_space_pipeline.py --help
```

## Workflow

1. `umap-explorer`
2. `sample-panel`
3. `align-panel`
4. `distance-explorer`

Auxiliary command:
- `export-clusters`

See:
- [docs/workflow.md](/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/seqscape/docs/workflow.md)
- [docs/widgets.md](/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/seqscape/docs/widgets.md)
- [docs/parameters.md](/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/seqscape/docs/parameters.md)

## Quickstart

Run the package tests, including the synthetic end-to-end smoke test:

```bash
cd seqscape
PYTHONPATH=src python -m unittest discover -s tests -v
```

For the staged command-line workflow, see:
- [examples/quickstart.md](/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/seqscape/examples/quickstart.md)

## Testing

```bash
cd seqscape
PYTHONPATH=src python -m unittest discover -s tests -v
```

The current suite includes:
- core utility tests
- a synthetic end-to-end smoke test covering all four pipeline stages

## Contributing

Contribution guidance is in:
- [CONTRIBUTING.md](/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/seqscape/CONTRIBUTING.md)

## License

SeqScape is distributed under the MIT license:
- [LICENSE](/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/seqscape/LICENSE)
