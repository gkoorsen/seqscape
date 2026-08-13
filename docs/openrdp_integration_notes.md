# OpenRDP integration: install route and validation notes

Written while validating the `collapse-panel` -> `recombination` fork against the
IS76 positive control (Panno et al. 2018; TYLCV backbone with a TYLCSV-derived
tract spanning the intergenic region).

## 1. OpenRDP is not on PyPI

`pip install openrdp` fails: `No matching distribution found for openrdp`.
Checked `openrdp`, `OpenRDP`, `open-rdp`, `rdp5` -- all HTTP 404 on PyPI. It
installs only from source:

```sh
pip install git+https://github.com/PoonLab/OpenRDP.git
```

### Do not install it into the seqscape environment

OpenRDP pins `numpy<2.0.0`. The seqscape env here runs numpy 2.4.6, and umap/numba
are built against it, so installing OpenRDP alongside downgrades numpy and breaks
`umap-explorer`. Because `resolve_openrdp()` looks for the `openrdp` console
script on `PATH` before trying `python -m openrdp`, a separate environment works
with no shared interpreter:

```sh
conda create -n openrdp python=3.11 'numpy=1.26' 'scipy=1.11' h5py pip
conda activate openrdp && pip install git+https://github.com/PoonLab/OpenRDP.git
PATH="$PATH:/path/to/envs/openrdp/bin" seqscape recombination ...
```

APPEND to `PATH`, do not prepend: prepending shadows the seqscape interpreter and
every `Bio` import fails. The SeqScape wrapper executes the `openrdp` script via
the adjacent Python interpreter when one exists, so the OpenRDP package can live
in its own environment even though the SeqScape command is run from the SeqScape
environment.

`3seq` and `geneconv` ship as bundled x86_64 Mach-O binaries. On Apple Silicon
they execute under Rosetta -- verified working, no action needed.

## 2. Five defects, all of which made the positive control silently negative

The control returned `raw events: 0  consensus events: 0` with exit code 0. Five
separate causes, found in order:

### (a) Repeated `-m` overwrote the method list

The wrapper built `-m rdp -m geneconv -m maxchi ...`. OpenRDP's `-m` is
`nargs="+"`, so repeating the flag REPLACES rather than accumulates: only the
last method (`siscan`) ran, found nothing, and exited 0. Methods must be passed
space-delimited after a single flag: `-m rdp geneconv maxchi ...`.

This is the highest-impact defect: it silently reduced a six-method consensus
screen to one method, so `--min-methods 2` could never be satisfied by
construction, on any dataset.

### (b) Parent columns were dropped

The parser matched header substrings `"major"` / `"minor"`, but OpenRDP emits
`Parent1` / `Parent2`. Every parent field came back empty. Now matches both.

Related: OpenRDP writes `-` for an unassigned parent; that is now mapped to empty
rather than carried through as a literal dash.

### (c) `Parent1`/`Parent2` were assumed to mean major/minor

They are POSITIONAL columns with no major/minor semantics. On the IS76 control,
`Parent1` is the *minor* (tract donor) parent -- the reverse of the assumption.
Reporting them as `parent_major`/`parent_minor` therefore asserted a backbone
that was wrong.

Fixed in two parts: the columns are now named `parent_a`/`parent_b`, and the
backbone is determined from the data by `polarise_event()`, which counts
phylogenetically informative sites (positions where the two candidate parents
differ) inside the called tract versus in the flanks. A genuine recombinant
follows one parent in the flanks and the other in the tract; the direction of
that flip identifies backbone and donor. When the pattern does not invert, no
call is made -- which is itself informative.

On IS76 the polarisation is unambiguous:

| region | informative sites | favours TYLCV | favours TYLCSV |
|---|---|---|---|
| tract 1343-1492 | 36 | 16 | 20 |
| flanks | 592 | 582 | 10 |

-> backbone TYLCV_IL_AM409201, donor TYLCSV_NC_003828, matching the published
event. The conserved nonanucleotide origin (TAATATTAC) sits at alignment column
1398, inside the called interval.

### (d) Consensus grouped on an exact parent tuple

Events were keyed on `(recombinant, major, minor)`. Methods legitimately disagree
about parents for the same event: GENECONV reported one parent where RDP reported
two. The exact-tuple key therefore split one real event into two single-method
groups, and `--min-methods 2` rejected both.

An event is defined by its recombinant and its breakpoint region, not by the
parent pair. Grouping is now per-recombinant with single-linkage merging on
breakpoint proximity; parents are reported as the union across supporting methods.

### (e) Reciprocal overlap alone is the wrong merge criterion

The first fix used a 0.5 reciprocal-overlap rule. The control still failed:
GENECONV (1382-1492) and RDP (1343-1420) describe the same event but reciprocally
overlap 0.4935 -- just under the cutoff.

Rather than tune the threshold until the control passed, the criterion was
corrected. What methods estimate is a pair of BREAKPOINT POSITIONS; the interval
between them is a by-product, and intervals of differing length can bracket the
same breakpoint while failing any fixed overlap fraction. Events now merge when
EITHER intervals overlap substantially OR corresponding breakpoints agree within
`breakpoint_tolerance` (default 100 nt).

### (f) Out-of-range p-values -- and the wrong fix for them

OpenRDP's RDP method emitted `36.363636...` as a "Pvalue" for the IS76 event.
The first fix here flagged such values (`pvalue_out_of_range`) and EXEMPTED them
from the `--pvalue` threshold, reasoning that comparing a non-probability to a
probability cutoff is meaningless.

**That was wrong, and it was wrong in the dangerous direction.** Reading
`openrdp/rdp.py:178`, the reported value is

```
pvalue = G * (L/N) * binom.sf(M-1, n=N, p=p)
```

where `G` is the number of triplets and `L/N` the alignment-length over
recombinant-region-length term. This is a Bonferroni-type correction, applied
uncapped. A corrected value ABOVE 1 therefore means NOT SIGNIFICANT -- the
correction has consumed all the evidence. Exempting those values inverted the
statistic: it admitted the least significant calls as though they had passed.

The scale of the error, on the 61-genome panel: of 54,300 RDP rows, **43,069
(79.3%) exceed 1** and would all have been admitted as significant support.

Out-of-range values are now clamped to 1 and then compared against the threshold
normally, so they fail any cutoff below 1. `pvalue_out_of_range` is retained as
an audit flag and such rows are still excluded from `min_pvalue`.

#### What this does to the positive control

On the control the arithmetic is exact and worth stating in full:

| term | value |
|---|---|
| G (triplets, 3 sequences) | C(3,3) = 1 |
| L (alignment length) | 2800 |
| N (recombinant region, 1420-1343) | 77 |
| L/N | 36.363636... |
| reported "Pvalue" | 36.363636... |
| **implied `binom.sf`** | **exactly 1.0** |

RDP's binomial test contributes NO evidence for this event. The reported number
is the correction factor alone. The earlier claim in these notes that the control
was recovered by "geneconv + rdp" was an artifact of the exemption bug.

**The IS76 control is recovered by GENECONV alone.** It is a genuine detection
(GENECONV p = 0.0 and 0.03171 on the two fragments, correctly polarised) but it
is a ONE-METHOD detection, and the wrapper's `--min-methods 2` default rejects
it. This is a property of the data and the methods, not something to tune away:
of the six methods, only GENECONV finds a published, experimentally validated
recombination event. Run a control-bearing screen at `--min-methods 1` and treat
method count as a ranking signal rather than a gate.

## 3. IUPAC ambiguity codes must be recoded

On the real 172-genome panel OpenRDP exits 1 with:

```
Alignment contains invalid characters KS.
Sequences can only contain A,T,G,C,-,N.
```

The panel contains 5 ambiguous bases (R, Y, S, W, K -- one each) across 4 of 172
records: 0.00067% of characters. These arise from consensus calling at
low-coverage or genuinely polymorphic sites.

`sanitise_alignment_for_openrdp()` recodes them to N and writes
`openrdp_input.fasta` alongside the results, leaving the user's alignment
untouched. Substitution is 1:1 so column coordinates -- and therefore breakpoint
positions -- are preserved. The count is printed and recorded in `summary.txt`
(`iupac_sites_recoded_to_n`) rather than applied silently. N is already OpenRDP's
missing-data character and an ambiguous base carries no signal for a triplet
scan, so this loses nothing.

## 4. Validation status

- IS76 positive control recovered: 1 consensus event at `--min-methods 1`, by
  GENECONV alone (see 2(f) -- RDP's binomial term is exactly 1.0 and contributes
  nothing), backbone and donor correctly polarised by `polarise_event()`,
  breakpoints bracketing the intergenic origin.
- 27 regression tests in `tests/test_recombination_consensus.py`, including the
  parser defects, p-value handling, role polarisation, bounded consensus merge,
  parentless-row filtering, paired GENECONV fragment handling, IUPAC recoding,
  gap-preserving alignment reads, and isolated OpenRDP script resolution.
- Mutation-tested: each fix reverted individually in an isolated tree, all five
  mutants killed by the suite. A test that passes with the bug reintroduced would
  be worthless, so this was checked rather than assumed.
- Full suite in the `seqscape` conda environment: 48 tests, all passing.

### Reproduce the IS76 positive-control validation

From the SeqScape repository root:

```sh
python -m venv /private/tmp/openrdp-validation-venv
/private/tmp/openrdp-validation-venv/bin/python -m pip install \
  'git+https://github.com/PoonLab/OpenRDP.git'

PATH="$PATH:/private/tmp/openrdp-validation-venv/bin" \
PYTHONPATH=src \
/Users/gerritkoorsen/opt/anaconda3/envs/seqscape/bin/python -m seqscape.cli \
  recombination \
  --alignment-fasta /Users/gerritkoorsen/ciderseq-mono/cideseq-mono/runs/tocsv_260203_full_aws_20260612/full_run/comparators/is76_control_alignment.fasta \
  --methods rdp,geneconv,maxchi,chimaera,threeseq,bootscan,siscan \
  --min-methods 1 \
  --pvalue 0.05 \
  --outdir runs/tylcv_is76_recombination_validation
```

Expected strict consensus:

| field | value |
|---|---|
| `recombinant_resolved` | `IS76_LN812978` |
| `backbone_parent` | `TYLCV_IL_AM409201` |
| `donor_parent` | `TYLCSV_NC_003828` |
| `methods` | `geneconv` |
| `breakpoint_start` / `breakpoint_end` | `1382` / `1492` |
| `min_pvalue` | `0.0` |

The significant GENECONV rows are complementary fragments for one event: a short
donor-like tract plus the longer backbone-like block in the same recombinant.
They are not two independent recombination events.

## 5. Caveat on the p-value semantics

RDP's reported value is Bonferroni-corrected over triplets, so it is strongly
panel-size dependent: the same triplet that is significant on its own becomes
non-significant inside a large panel purely through the `G` term (on 61 genomes
G = 35,990, so any raw survival value above 2.78e-05 corrects to above 1). This
is defensible as multiple-testing control, but it means RDP support is not
comparable across panels of different size, and a large panel will suppress RDP
almost entirely.

Practical consequence: `--min-methods 2` is NOT a safe default on a large panel,
because the positive control itself does not meet it. Screen at
`--min-methods 1`, rank by method count and by GENECONV significance, and confirm
anything of interest by re-testing the specific triplet in RDP5 directly. The
wrapper is a triage tool over a dereplicated panel, not a
substitute for per-event scrutiny.

Two of the seven methods (`geneconv`, `threeseq`) are distributed under
non-commercial licences; `--methods` can exclude them.

## 6. Runtime: which methods you can actually afford

Per-method cost measured on a 12-genome subset of the real panel (4,326 columns),
then projected by triplet count. OpenRDP enumerates C(n,3) triplets and runs every
method on each, so cost is cubic in panel size.

| method | ms/triplet | projected, n=172 (833,340 triplets) |
|---|---|---|
| siscan | 6209.04 | **1437 h (60 days)** |
| bootscan | 168.13 | 38.9 h |
| chimaera | 32.99 | 7.6 h |
| geneconv | 25.18 | 5.8 h |
| rdp | 22.20 | 5.1 h |
| maxchi | 6.72 | 1.6 h (crashes, see below) |

This is why the first full-panel run appeared to hang: `siscan` is 37x slower than
the next-slowest method and would never have completed.

On the IS76 control, `siscan` and `bootscan` each returned ZERO events -- the
detection came entirely from `geneconv` (see 2(f): RDP's apparent support was an
artifact of the p-value exemption bug). They cost 99% of the runtime and
contributed nothing to the one event we can independently verify. `siscan` should
be off by default; `bootscan` is affordable only on a small panel.

### Panel size is the real lever

`collapse-panel --identity-threshold` controls panel size, and all 31 focal
genomes are retained at EVERY threshold in the sweep -- so a lower threshold costs
focal coverage nothing, only background resolution.

| threshold | panel | triplets | 4 fast methods | + bootscan |
|---|---|---|---|---|
| 95.0 | 61 | 35,990 | **0.9 h** | 2.6 h |
| 97.0 | 111 | 221,815 | 5.4 h | 15.7 h |
| 98.1 | 172 | 833,340 | 20.2 h | 59.1 h |
| 99.0 | 317 | 5,259,030 | 127 h | 373 h |

Recommended: screen at 95% first (under an hour, all focal genomes present), then
re-run any triplet of interest at 98.1% or in RDP5 directly for the definitive
call. Treat the wide-panel run as confirmation, not discovery.

### maxchi crashes on real data

`ZeroDivisionError: float division by zero` at `openrdp/maxchi.py:192`
(`0.05/l` where `l` is the count of chi-squared windows, zero for triplets that
yield none). This is an upstream bug, not a wrapper bug. Since one crashing method
was aborting the entire run and discarding the results of methods that had already
finished, the wrapper now invokes methods ONE PER SUBPROCESS, records per-method
status in `summary.txt` (`method_status.<method>`), and continues past failures --
exiting only if every method fails.

### Reference mode does not help

OpenRDP's `-r` reduces enumeration from C(n,3) to n_query x C(n_ref,2)
(`openrdp/__init__.py:262`), which would map neatly onto "focal genomes as
candidate recombinants, panel as candidate parents". It is not usable: `rdp`
crashes on it (`ValueError: math domain error`, `rdp.py:178`), and `geneconv`,
`maxchi` and `chimaera` all recover NOTHING from the IS76 control in this mode.
A mode that fails the positive control cannot be used regardless of its speed
advantage. All-vs-all on a smaller panel is the supported path.

## 7. A seventh defect: OpenRDP's recombinant label is not stable

Found while validating the per-method runner. Running `-m rdp` alone versus all
six methods gives the same triplet and the same breakpoints, with the ROLES
PERMUTED:

| methods run | Recombinant | Parent1 | Parent2 | breakpoints |
|---|---|---|---|---|
| all six | IS76 | TYLCSV | TYLCV | 1343-1420 |
| rdp alone, or rdp+geneconv | **TYLCSV** | IS76 | TYLCV | 1343-1420 |

Deterministic across seeds 3, 7 and 42 -- so this is not a permutation-seed
effect but a dependence of the reported role on the method set. It means the
`Recombinant` column cannot be used as a grouping key or reported as a finding.

Two consequences for the wrapper:

1. Consensus now groups on the SET OF SEQUENCES INVOLVED plus breakpoint
   agreement (two events merge when they share >= 2 sequences and their
   breakpoints agree). OpenRDP's own labels are preserved in
   `openrdp_recombinant_labels` for audit, not used for grouping.
2. `assign_roles()` recovers the recombinant from the alignment: each triplet
   member is tried as the putative recombinant, and the assignment whose
   informative-site pattern inverts between tract and flanks is kept. Reported in
   `recombinant_resolved` with `role_support`. Triplets where no assignment
   inverts get an empty call rather than a guess.

Verified: the control now yields the same consensus event under `rdp,geneconv`
and under `rdp,geneconv,maxchi,chimaera` -- `recombinant_resolved=IS76_LN812978`,
backbone TYLCV, donor TYLCSV, `role_support=628`, while
`openrdp_recombinant_labels` records both conflicting labels.

## 8. An eighth defect, visible only at panel scale: transitive merge chaining

The consensus grouper merges two events when they share sequences AND their
breakpoints agree. On the 3-sequence positive control this is correct and cannot
misbehave. On the real 61-genome panel it collapses:

| min_methods | consensus events | largest event |
|---|---|---|
| 1 | 62 | -- |
| 2 | 8 | **59 of 61 sequences, spanning columns 0-3760 of 3689** |
| 3 | 2 | same |

A recombination event involves three sequences. An "event" containing 59 is not
an event. The cause is transitive chaining: A merges with B, B merges with C, so
A, B and C land in one group even when A and C share no interval. Across 142,985
raw events on a dense panel, the chain closes over almost the whole panel and
every reported event spans the entire genome.

**This is why the positive control is necessary but NOT sufficient as a
validation.** A 3-sequence control exercises the parsing, the polarisation and
the p-value handling, but it cannot exercise the grouping at scale -- the failure
mode needs a dense panel to appear at all. Any future change to the merge rule
must be checked against BOTH the control and a panel-scale run.

FIXED. The merge is now bounded at `max_members=3` -- the biological definition
of an event (a recombinant and two parents), not a tuning knob. Candidate
clusters are found through a pair index rather than a linear scan, since a merge
requires >= 2 shared sequences; grouping the 142,985 panel events went from
~10 min to 6 s.

| min_methods | before | after | largest group |
|---|---|---|---|
| 1 | 62 | 62,862 | 3 |
| 2 | 8 | 4,856 | 3 |
| 3 | 2 | 262 | 3 |

Covered by `TestMergeBound` in `tests/test_recombination_consensus.py`, which
builds a synthetic chain (A~B, B~C, C~D at identical breakpoints) because the
3-sequence control cannot exercise this path. Four tests: the bound holds, the
chain does not collapse to one event, RAISING the bound restores chaining
(proving the bound is what prevents it), and the control still passes. Removing
the bound kills two of them.

### A second consequence: merging widens short tracts

Consensus breakpoints are medians across the supporting methods, so a precise
short tract merged with wide bootscan/rdp intervals is reported wide. The 193 nt
GENECONV call on MT878433 (cols 2489-2682) is absorbed into a
bootscan+geneconv+rdp group reported as 100-2664, i.e. 2,564 nt. Short-tract
screening must therefore run on the PER-TRIPLET table
(`recombination_events.tsv`), with the consensus used to ask whether a candidate
has independent multi-method support -- not to measure its extent.

## 9. Screen result on the 61-genome panel

Completed run, 2.4 h wall clock: rdp 54,300 events / 632 s; geneconv 481 / 600 s;
chimaera 0 / 2,841 s; bootscan 88,204 / 3,999 s; maxchi crashed after 40 s.

Using GENECONV (the only method that recovers the control) and its per-triplet
calls rather than the chained consensus:

- 478 of 481 GENECONV events significant at p <= 0.05
- 213 significant events involve >= 1 of the 31 focal genomes
- **all 31 focal genomes are implicated in at least one event**

The last point is the substantive finding: on begomoviruses, "carries detectable
recombination" discriminates nothing. Applying the IS76 signature from
Panno et al. (2018) -- a SHORT tract at the IR/C1 junction, 77 nt in the control --
narrows 213 events to 11 with a tract <= 400 nt, and to ONE with a tract <= 400 nt
and a breakpoint near the IR: MT878433, 193 nt at reference nt 2233-2405
(TYLCSV NC_003828 coordinates), p = 0.00714, partner GU951759.

That candidate is IR-PROXIMAL, not clearly IR-spanning, and sits close enough to
the window edge that widening or narrowing the IR definition changes the call. It
is a lead to confirm in RDP5 directly, not a result. Genome rotation was checked
and is consistent (55/61 carry TAATATTAC spanning the linearisation junction at a
fixed offset), so the breakpoints are not rotation artifacts.

## 10. Re-analysis with the bounded consensus

Re-derived from the same raw per-method CSVs (no rescan needed):

- 4,856 consensus events at >= 2 methods, every one a proper triplet
- 4,099 involve >= 1 focal genome, spanning 29 of 31 focal genomes
- 1,007 have a tract <= 400 nt; 473 of those are also IR-proximal

But all 473 IS76-like consensus events are `bootscan,rdp`. **Not one includes
GENECONV** -- the only method that recovers the positive control. Methods with no
demonstrated sensitivity on a known event agreeing with each other is not
corroboration; bootscan and rdp are also the two highest-volume methods (88,204
and 54,300 raw events), so their agreement is expected under noise.

Conversely the GENECONV per-triplet candidate (MT878433 / GU951759, 193 nt at
reference nt 2233-2405) has no short-tract multi-method backing: the seven
consensus events containing both sequences all report tracts of 240 nt or wider,
and the group that absorbs the actual call reports 2,564 nt.

**Conclusion: the screen does not identify a resistance-breaking sequence.** It
produces one lead (MT878433) supported by a single method at a single locus,
which is below the bar for reporting. What the screen does establish, and this is
worth having: all 31 focal genomes are recombinant, so recombination per se does
not discriminate, and the IS76 signature does not appear cleanly in this panel.
