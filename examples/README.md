# SeqScape examples

Two fully worked case studies are included. Each contains fetch scripts, a
prepared reference set, and a reproducible SeqScape workflow with its output
under `runs/`.

| Case study | Anchor taxonomy | Full-collection size | Panel size |
|---|---|---:|---:|
| [`norovirus_vp1/`](norovirus_vp1/) | Chhabra et al. 2019 — 10 genogroups, 48 confirmed VP1 genotypes | ~4,700 VP1 sequences | 500 |
| [`begomovirus_dna_a/`](begomovirus_dna_a/) | Brown et al. 2015 — 39 begomovirus groups, 91%/94% species/strain cutoffs | ~15,700 DNA-A sequences | 500 |

Each case-study directory contains:

- `README.md` — dataset provenance, anchor paper, and workflow notes
- `fetch_dataset.py` — NCBI Entrez retrieval of the full sequence pool
- `fetch_*_references.py` / `parse_*.py` — reference-set assembly from the
  anchor paper's supplementary material
- `data/` — fetched FASTA and manifest files
- `runs/` — SeqScape pipeline output (UMAP explorer, panel, alignment,
  distance explorer)
- auxiliary scripts for genotype assignment and comparison against the
  anchor taxonomy

## Tracked output widgets

The norovirus VP1 example includes selected rendered HTML outputs for quick
inspection without rerunning the full workflow:

- [AF Leiden parameter widget](norovirus_vp1/runs/af/af_leiden_parameter_widget.html)
- [Panel review widget](norovirus_vp1/runs/review/panel_review_widget.html)
- [Panel square heatmap](norovirus_vp1/runs/review/heatmap_panel_square.html)
- [Genome x reference heatmap](norovirus_vp1/runs/review/heatmap_genome_x_reference.html)
- [Chhabra 302 panel review widget](norovirus_vp1/runs/chhabra302/review/panel_review_widget.html)
- [Chhabra 302 panel square heatmap](norovirus_vp1/runs/chhabra302/review/heatmap_panel_square.html)

These files are self-contained browser widgets. When using a downloaded copy of
the repository or release archive, open the `.html` files directly in a browser;
no command-line step or local web server is required. GitHub's web interface
shows HTML files as source code, so the interactive view is best accessed from a
downloaded copy.

Suggested first views:

- `norovirus_vp1/runs/af/af_leiden_parameter_widget.html` for the
  alignment-free parameter search and Leiden cluster states
- `norovirus_vp1/runs/review/panel_review_widget.html` for the full
  500-sequence norovirus VP1 panel review
- `norovirus_vp1/runs/chhabra302/review/panel_review_widget.html` for the
  Chhabra reference-set validation benchmark

## Reproduction

The `data/` directories and most `runs/` outputs are **not** tracked in the
repository — the input FASTA pools (~4,700 norovirus and ~15,700 begomovirus
sequences) are regenerated directly from NCBI Entrez so that the examples
always reflect the current GenBank state rather than a frozen snapshot. To
reproduce a case study end-to-end:

```bash
cd examples/norovirus_vp1          # or examples/begomovirus_dna_a
python fetch_dataset.py --email you@example.org --outdir data/
# (norovirus only) download chhabra_s1.xlsx manually and place it in data/
python parse_chhabra_s1.py         # writes data/chhabra_s1_accessions.txt
python fetch_cdc_references.py     # norovirus references
# then follow examples/quickstart.md with --input-fasta data/*.fasta
```

A full end-to-end run takes ~20 min on 8 CPU threads for a 500-sequence
panel. Estimated disk footprint is ~500 MB per case study for the full
pipeline outputs.
