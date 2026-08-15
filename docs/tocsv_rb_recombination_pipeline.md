# ToCSV RB Recombination Pipeline

This pipeline screens the clean phased ToCSV genomes for candidate
resistance-breaking recombination signatures in the origin/IR region.

The key design choice is to collapse whole-genome redundancy only after
protecting origin-window diversity. That prevents near-identical genomes from
inflating the recombination test set, while avoiding the obvious failure mode:
discarding a genome that is almost identical overall but carries a distinct
origin/IR tract.

## Run

From the SeqScape repository:

```bash
scripts/run_tocsv_rb_recombination_pipeline.sh --skip-recombination
```

On this workstation the script automatically uses:

```text
/Users/gerritkoorsen/opt/anaconda3/envs/seqscape/bin/python
```

when `PYTHON` is not set. It also sets writable Matplotlib and numba cache
directories inside the run folder, which avoids UMAP/numba cache errors from
the base shell environment.

An existing `umap-explorer` directory is reused only when its `summary.txt`
matches the same clean FASTA and reference FASTA. If the comparator FASTA differs
from the old UMAP run, the script rebuilds UMAP in the new output directory.

That produces the protected-ID set, dereplicated panel, origin-window FASTA, and
MAFFT alignment. Remove `--skip-recombination` when OpenRDP is installed and you
are ready to run the recombination screen:

```bash
scripts/run_tocsv_rb_recombination_pipeline.sh
```

The defaults point to the current local full ToCSV run:

```text
clean genomes:
/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/runs/tocsv_260203_full_aws_20260612/full_run/analysis/seqscape_clean_phased_20260625/clean_phased_2400_3300_unique.fasta

support table:
/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/runs/tocsv_260203_full_aws_20260612/full_run/analysis/seqscape_clean_phased_20260625/clean_phased_2400_3300_unique_support.tsv

comparators:
/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/runs/tocsv_260203_full_aws_20260612/full_run/comparators/comparator_panel.fasta

anchor references:
/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/runs/tocsv_260203_full_aws_20260612/full_run/analysis/phylo_compare_20260612_223925/references.fasta
```

## Stages

1. Extract a 400 nt circular-genome window centred on the conserved
   begomovirus origin motif `TAATATTAC`.
2. Collapse those local origin windows at 98.1 percent identity and select one
   well-supported local representative per regional haplotype group.
3. Combine the clean local genomes with the comparator FASTA as collapsible
   whole-genome input. Anchor references are kept separate.
4. Run or reuse `umap-explorer` for whole-genome context.
5. Run `collapse-panel` at 95 percent whole-genome identity, using the selected
   origin-window representatives as protected focal IDs. Anchor references are
   always retained; comparator genomes are dereplicated unless they are distinct.
6. Extract the same origin-centred window from the dereplicated genome panel.
7. Align the regional panel with MAFFT.
8. Run `seqscape recombination` with a fast OpenRDP method set.

## Important Outputs

```text
protected_origin_ids.txt
protected_origin_groups.tsv
collapse_panel/representative_panel.fasta
collapse_panel/representative_panel_manifest.tsv
regions/panel_origin_core_400nt.fasta
regions/panel_origin_core_400nt_aligned.fasta
recombination_origin_core_400nt/recombination_events.tsv
recombination_origin_core_400nt/recombination_consensus.tsv
```

Use `recombination_events.tsv` for screening and `recombination_consensus.tsv`
for higher-confidence events. With large panels, `--min-methods 1` is a triage
setting, not a final biological call. Any candidate should be rechecked as a
specific triplet or small focused panel before being described as a likely
resistance-breaking recombinant.

## Tuning

For a faster first pass:

```bash
COLLAPSE_IDENTITY_THRESHOLD=95.0 \
REGION_IDENTITY_THRESHOLD=98.1 \
scripts/run_tocsv_rb_recombination_pipeline.sh
```

For a more conservative regional protection set:

```bash
REGION_IDENTITY_THRESHOLD=99.0 \
scripts/run_tocsv_rb_recombination_pipeline.sh
```

For a broader origin window:

```bash
WINDOW_SIZE=1000 \
scripts/run_tocsv_rb_recombination_pipeline.sh
```

The anomalous ToCSV FASTA should not be mixed into this primary run. Run it as a
separate rescue/QC branch after filtering for plausible complete circular
genomes or clearly interpretable origin-region candidates.
