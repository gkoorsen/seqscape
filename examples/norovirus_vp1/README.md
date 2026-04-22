# Case study 1 — Norovirus VP1

Anchor paper: Chhabra et al. 2019, *J Gen Virol* — established the current
10-genogroup / 49-genotype norovirus taxonomy from pairwise-distance analysis
of 305 curated complete VP1 reference sequences.

## Dataset plan

- **Full collection** (this directory): Norovirus nuccore records, 1500–2500 nt,
  filtered to complete VP1 CDSs 1500–1800 nt, ambiguity ≤1 %, deduplicated by
  sequence hash. Current GenBank pool ≈ 5,906 records before filtering.
- **Reference panel (anchors):** Chhabra 2019 Table S1 accessions (305 seqs).
  PMC holds this paper outside the Open Access bulk set, so Table S1 has to
  be downloaded manually from the journal page
  (https://www.microbiologyresearch.org/content/journal/jgv/10.1099/jgv.0.001318
  → Supplementary Material → Excel). Save as `chhabra_s1.xlsx`, then run
  `parse_chhabra_s1.py` to extract the accession list.
- **Ground-truth labels:** CDC [Human Calicivirus Typing tool]
  (https://calicivirustypingtool.cdc.gov/) genotype calls on the 500-panel.

## Fetch

```bash
python fetch_dataset.py \
  --email you@example.org \
  --outdir data/
```

Outputs `data/norovirus_vp1.fasta` and `data/norovirus_vp1_manifest.tsv`.

## SeqScape workflow

Follow [`../quickstart.md`](../quickstart.md) with
`--input-fasta data/norovirus_vp1.fasta` and the Chhabra S1 accessions as
`--reference-fasta`. Target panel size 500.
