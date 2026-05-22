# SeqScape

`SeqScape` is a staged workflow for sequence-space exploration, representative sampling, and alignment-based review.

It is being extracted from the larger exploratory repository into a smaller publication-facing software package for the planned JOSS and BMC Bioinformatics papers.

Repository:
- `https://github.com/gkoorsen/seqscape`

## Current Status

- package CLI implemented in [src/seqscape](src/seqscape)
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
python -m pip install -e .
```

The `align-panel` stage can use `muscle`, `mafft`, or Biopython `pairwisealigner`.

External tools:
- `MUSCLE` is recommended for production pairwise alignment panels
- `MAFFT` is optional if you want the `mafft` alignment mode

## CLI

Package entry point:

```bash
cd seqscape
seqscape --help
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
- [docs/workflow.md](docs/workflow.md)
- [docs/widgets.md](docs/widgets.md)
- [docs/parameters.md](docs/parameters.md)

## Interactive Example Widgets

The norovirus VP1 example includes rendered browser widgets that can be opened
without rerunning the workflow. They are static HTML files and do not need a web
server.

Important: GitHub displays HTML files as source code. For interactive use,
download the repository or release archive, unzip it, and open the files below
directly in a browser such as Chrome, Edge, Firefox, or Safari.

Open these files from the downloaded `seqscape/` folder:

- `examples/norovirus_vp1/runs/af/af_leiden_parameter_widget.html` — alignment-free Leiden parameter explorer
- `examples/norovirus_vp1/runs/review/panel_review_widget.html` — full norovirus VP1 panel review widget
- `examples/norovirus_vp1/runs/review/heatmap_panel_square.html` — panel-by-panel identity heatmap
- `examples/norovirus_vp1/runs/review/heatmap_genome_x_reference.html` — panel-by-reference identity heatmap
- `examples/norovirus_vp1/runs/chhabra302/review/panel_review_widget.html` — Chhabra reference-set validation widget
- `examples/norovirus_vp1/runs/chhabra302/review/heatmap_panel_square.html` — Chhabra reference-set heatmap

The same list is maintained in [examples/README.md](examples/README.md), with
links that work when browsing the repository.

## Quickstart

Run the package tests, including the synthetic end-to-end smoke test:

```bash
cd seqscape
PYTHONPATH=src python -m unittest discover -s tests -v
```

The end-to-end smoke test uses small synthetic FASTA records generated inside
the test itself, so it does not require any external data downloads. It runs all
four stages and writes temporary outputs for the alignment-free explorer,
representative panel, alignment matrices, and distance explorer.

For the staged command-line workflow, see:
- [examples/quickstart.md](examples/quickstart.md)

For full worked biological case studies, see:
- [examples/README.md](examples/README.md)
- [examples/norovirus_vp1/README.md](examples/norovirus_vp1/README.md)
- [examples/begomovirus_dna_a/README.md](examples/begomovirus_dna_a/README.md)

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
- [CONTRIBUTING.md](CONTRIBUTING.md)

## License

SeqScape is distributed under the MIT license:
- [LICENSE](LICENSE)
