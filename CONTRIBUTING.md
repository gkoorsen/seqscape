# Contributing

SeqScape is under active development as a publication-facing extraction from the larger exploratory repository.

## Scope

Contributions should stay aligned with the four-stage SeqScape workflow:

1. `af-widget`
2. `sample-panel`
3. `align-panel`
4. `review-widget`

Avoid adding exploratory one-off analyses or unrelated legacy workflows to this package.

## Development

Recommended setup:

```bash
cd seqscape
conda env create -f environment.yml
conda activate seqscape
```

Run tests before proposing changes:

```bash
cd seqscape
PYTHONPATH=src python -m unittest discover -s tests -v
```

## Coding Expectations

- keep the package surface small and publication-facing
- prefer changes inside `src/seqscape`
- keep the CLI behavior stable unless a deliberate interface change is being made
- add or update tests when changing stage behavior
- update docs when changing parameters, outputs, or workflow assumptions

## Reporting Issues

When reporting a problem, include:
- command used
- environment details
- input sizes
- exact error output
- whether the issue occurs in `af-widget`, `sample-panel`, `align-panel`, or `review-widget`

## Maintainer

Current maintainer:
- Gerrit Koorsen
