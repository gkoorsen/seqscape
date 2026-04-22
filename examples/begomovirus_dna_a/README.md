# Case study 2 — Begomovirus DNA-A

Anchor papers:

- Muhire et al. 2014, *PLOS ONE* — introduced SDT; demonstrated all-vs-all
  pairwise analysis on 1,000 begomovirus genomes (499,500 alignments;
  ~62 h serial / 1.6 h on 40 cores).
- Brown et al. 2015, *Arch Virol* — revised begomovirus taxonomy from 3,123
  DNA-A genomes; established the 91 %/94 % species/strain identity thresholds.

## Dataset plan

- **Full collection** (this directory): Begomovirus nuccore records, 2500–3200 nt
  (the DNA-A range), with satellites, DNA-B segments, and clearly partial
  records excluded, ambiguity ≤1 %, deduplicated by sequence hash. Current
  GenBank pool ≈ 15,730 records before filtering.
- **Reference panel (anchors):** Brown 2015 supplementary material (MOESM1)
  contains 39 per-group SDT pairwise-identity matrices — the accessions used
  as group exemplars across all 39 groups are listed in
  `brown_2015_group_exemplar_accessions.txt` (252 unique accessions extracted
  directly from the MOESM1 xls sheet headers). Pass this list to SeqScape via
  `--reference-fasta`. The full 3,123-genome Brown dataset itself is not
  redistributed by the authors — use the 252 exemplars as anchors and the
  current GenBank pull as the exploration collection.
- **Ground-truth / comparator:** run SDT on the 500-panel and compare the
  91 % species cutoff against af-Leiden clusters and against complete-linkage
  agglomerative clustering at distance 0.09.

## Fetch

```bash
python fetch_dataset.py \
  --email you@example.org \
  --outdir data/
```

Outputs `data/begomovirus_dna_a.fasta` and `data/begomovirus_dna_a_manifest.tsv`.

## SeqScape workflow

Follow [`../quickstart.md`](../quickstart.md) with
`--input-fasta data/begomovirus_dna_a.fasta` and the Brown 2015 accessions as
`--reference-fasta`. Target panel size 500.
