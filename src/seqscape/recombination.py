"""Formal recombination detection for SeqScape panels.

WHAT THIS IS
------------
A thin, reproducible wrapper around OpenRDP -- an open-source re-implementation
of the RDP4/RDP5 recombination detection program -- plus optional mapping of
detected breakpoints onto ORF/intergenic coordinates.

WHY A WRAPPER AND NOT AN IMPLEMENTATION
---------------------------------------
Everything else in SeqScape is pairwise identity -> distance -> ordination /
clustering. Breakpoint detection is a different class of statistics
(phylogenetic incongruence, substitution-pattern scanning, permutation tests).
Re-implementing those tests here would be neither defensible nor citable.
OpenRDP runs the established method suite (RDP, GENECONV, MaxChi, Chimaera,
3Seq, BootScan, SiScan) against a multiple sequence alignment and reports
breakpoints with putative parents and p-values, which is what a reviewer will
expect to see.

The consensus convention in the field -- and the one applied here -- is that a
recombination event is only called when it is supported by several independent
methods; `--min-methods` defaults to 3 for that reason. Single-method hits are
retained in the raw output but are not promoted to the consensus table.

LICENCE NOTE
------------
Two of the methods OpenRDP invokes ship as external binaries with non-commercial
/ academic-use terms (3Seq is CC BY-NC-SA; GENECONV is academic-use). That is
fine for research use, but the restriction is real -- do not redistribute
outputs commercially. Use `--methods` to exclude them if required.

INPUT
-----
An ALIGNED FASTA. `seqscape align-panel` output, or any MAFFT/MUSCLE alignment.
For resistance-breaking work the intended input is a region alignment produced
by `extract_genome_regions.py` (e.g. the CP_REP window), not a whole-genome
alignment -- an IS76-type event is ~76 nt and is diluted at genome scale.

OUTPUTS
-------
  recombination_raw.csv        OpenRDP output, unmodified
  recombination_events.tsv     parsed per-method events
  recombination_consensus.tsv  events passing --min-methods, with parents
  breakpoint_regions.tsv       breakpoints mapped to ORF/IR coordinates (if
                               --region-manifest given)
  summary.txt                  run parameters and counts
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

# The method suite OpenRDP exposes. 3seq/geneconv are external binaries.
ALL_METHODS = ["rdp", "geneconv", "maxchi", "chimaera", "threeseq", "bootscan", "siscan"]


def resolve_openrdp() -> list[str] | None:
    """Return a runnable OpenRDP command, or None if it is not installed.

    Tries the console script first, then `python -m openrdp`, so the wrapper
    works with either a pip install or a source checkout on PYTHONPATH.
    """
    found = shutil.which("openrdp")
    if found:
        # OpenRDP installs its console script with `#!/usr/bin/env python3` in
        # some environments. Execute it with the adjacent interpreter when that
        # exists, otherwise the script can resolve to SeqScape's Python and fail
        # to import OpenRDP from its own isolated environment.
        adjacent_python = Path(found).with_name("python")
        if adjacent_python.exists() and adjacent_python.is_file():
            return [str(adjacent_python), found]
        return [found]
    probe = subprocess.run([sys.executable, "-c", "import openrdp"], capture_output=True)
    if probe.returncode == 0:
        return [sys.executable, "-m", "openrdp"]
    return None


def parse_openrdp_csv(path: Path) -> list[dict]:
    """Parse OpenRDP's CSV output into event dicts.

    OpenRDP's column naming has varied between releases, so columns are matched
    case-insensitively by substring rather than by exact header, and unmatched
    rows are skipped rather than crashing the run.
    """
    if not path.exists():
        return []
    rows: list[dict] = []
    with path.open(newline="") as fh:
        for raw in csv.DictReader(fh):
            norm = { (k or "").strip().lower(): (v or "").strip() for k, v in raw.items() }

            def pick(*fragments, default=""):
                for frag in fragments:
                    for key, val in norm.items():
                        if frag in key:
                            # OpenRDP writes "-" for an unassigned parent
                            return "" if val == "-" else val
                return default

            start, end = pick("start"), pick("end")
            if not start or not end:
                continue
            try:
                start_i, end_i = int(float(start)), int(float(end))
            except ValueError:
                continue
            pval_raw = pick("p-value", "pvalue", "p_value")
            try:
                pval = float(pval_raw)
            except ValueError:
                pval = float("nan")
            # OpenRDP's RDP method can emit a score rather than a probability
            # (values >1 have been observed). Such a value cannot be compared
            # against a p-value threshold; mark it so the event is reported as
            # unfiltered rather than silently failing the cutoff.
            pval_out_of_range = pval == pval and not (0.0 <= pval <= 1.0)
            rows.append({
                # OpenRDP emits Parent1/Parent2, which are POSITIONAL columns and
                # carry no major/minor semantics: on the IS76 control, Parent1 is
                # the minor (tract donor) parent. Only true RDP4-style Major/Minor
                # headers are semantic. Do not relabel positional parents as
                # major/minor -- determine which is the backbone from informative
                # sites in the flanking regions instead.
                "recombinant": pick("recombinant"),
                "parent_a": pick("major", "parent1"),
                "parent_b": pick("minor", "parent2"),
                "start": start_i,
                "end": end_i,
                "method": pick("method", "test").lower(),
                "pvalue": pval,
                "pvalue_out_of_range": pval_out_of_range,
            })
    return rows


def read_alignment(path: Path) -> dict[str, str]:
    """Read an aligned FASTA preserving gap columns.

    `io_utils.load_fasta_records` strips '-', which is correct for unaligned
    input but destroys the column coordinates that OpenRDP breakpoints index
    into. Alignment-aware consumers must use this instead.
    """
    seqs: dict[str, list[str]] = {}
    current = None
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current is not None:
                seqs[current].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


IUPAC_AMBIGUITY = set("RYSWKMBDHV")


def sanitise_alignment_for_openrdp(alignment: dict[str, str], path: Path) -> dict:
    """Recode IUPAC ambiguity codes to N and write an OpenRDP-safe alignment.

    OpenRDP rejects any alignment containing characters outside A,T,G,C,-,N with
    "Alignment contains invalid characters", which real begomovirus panels trip:
    consensus calling leaves occasional R/Y/S/W/K at low-coverage or genuinely
    polymorphic sites. Recoding them to N is the conventional handling -- an
    ambiguous base carries no phylogenetic signal for a triplet scan, and N is
    already OpenRDP's missing-data character.

    Column coordinates are preserved (substitution is 1:1, no re-alignment), so
    breakpoints remain comparable to the input alignment. Returns a report so the
    substitution is recorded rather than silent.
    """
    recoded: dict[str, str] = {}
    per_record: dict[str, int] = {}
    codes: dict[str, int] = {}
    for name, seq in alignment.items():
        out = []
        n = 0
        for ch in seq:
            if ch in IUPAC_AMBIGUITY:
                codes[ch] = codes.get(ch, 0) + 1
                n += 1
                out.append("N")
            else:
                out.append(ch)
        recoded[name] = "".join(out)
        if n:
            per_record[name] = n

    with path.open("w") as fh:
        for name, seq in recoded.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    return {
        "n_sites_recoded": sum(per_record.values()),
        "n_records_affected": len(per_record),
        "codes": dict(sorted(codes.items())),
        "records": dict(sorted(per_record.items())),
    }


def assign_roles(alignment: dict[str, str], candidates: list[str],
                 start: int, end: int) -> dict:
    """Decide which member of a triplet is the recombinant, from the alignment.

    OpenRDP's `Recombinant` column cannot be trusted for this: which member of a
    triplet it labels the recombinant depends on which methods were run (see
    build_consensus). The role is instead recovered from the data by trying each
    member as the putative recombinant and keeping the assignment whose
    informative-site pattern inverts between the tract and the flanks -- the
    defining signature of a recombinant. Ties and non-inverting triplets yield no
    call, which is reported rather than guessed.
    """
    members = [c for c in dict.fromkeys(candidates) if c and c in alignment]
    if len(members) < 3:
        # Two-sequence calls (e.g. GENECONV) cannot be polarised on their own.
        base = polarise_event(alignment, members[0] if members else "",
                              members[1:] if len(members) > 1 else [], start, end)
        base["recombinant_resolved"] = ""
        base["role_support"] = 0
        return base

    best, best_support = None, -1
    for i, rec in enumerate(members):
        parents = [m for j, m in enumerate(members) if j != i]
        out = polarise_event(alignment, rec, parents, start, end)
        if not out["backbone_parent"]:
            continue
        # Support = strength of the inversion: informative sites backing the
        # tract call plus those backing the flank call.
        support = out["informative_sites_tract"] + out["informative_sites_flank"]
        if support > best_support:
            best, best_support = dict(out, recombinant_resolved=rec,
                                      role_support=support), support
    if best is None:
        empty = polarise_event(alignment, members[0], members[1:], start, end)
        empty["recombinant_resolved"] = ""
        empty["role_support"] = 0
        return empty
    return best


def polarise_event(alignment: dict[str, str], recombinant: str, parents: list[str],
                   start: int, end: int) -> dict:
    """Decide which parent is the backbone and which donated the tract.

    OpenRDP's Parent1/Parent2 columns are positional, not semantic, so the
    backbone cannot be read off them. It is determined here the same way the
    IS76 event is established in the literature: count phylogenetically
    informative sites where the two candidate parents differ, and see which
    parent the recombinant follows INSIDE the called tract versus OUTSIDE it.
    A genuine recombinant follows one parent in the flanks and the other within
    the tract; the sign of that flip is the evidence.

    Returns counts plus `backbone_parent` / `donor_parent`, or empty strings when
    the pattern does not invert (which is itself diagnostic -- it means the call
    is not supported by informative sites).
    """
    usable = [p for p in parents if p and p in alignment]
    if recombinant not in alignment or len(usable) != 2:
        return {"backbone_parent": "", "donor_parent": "", "informative_sites_tract": 0,
                "informative_sites_flank": 0, "tract_supports": "", "flank_supports": ""}
    rec = alignment[recombinant]
    p1, p2 = alignment[usable[0]], alignment[usable[1]]
    n = min(len(rec), len(p1), len(p2))

    def tally(lo: int, hi: int) -> tuple[int, int]:
        c1 = c2 = 0
        for i in range(max(0, lo), min(n, hi)):
            x, y, z = rec[i], p1[i], p2[i]
            if x == "-" or y == "-" or z == "-" or y == z:
                continue
            if x == y:
                c1 += 1
            elif x == z:
                c2 += 1
        return c1, c2

    in1, in2 = tally(start, end)
    f1a, f2a = tally(0, start)
    f1b, f2b = tally(end, n)
    fl1, fl2 = f1a + f1b, f2a + f2b

    tract_winner = usable[0] if in1 > in2 else (usable[1] if in2 > in1 else "")
    flank_winner = usable[0] if fl1 > fl2 else (usable[1] if fl2 > fl1 else "")
    inverted = bool(tract_winner) and bool(flank_winner) and tract_winner != flank_winner
    return {
        "backbone_parent": flank_winner if inverted else "",
        "donor_parent": tract_winner if inverted else "",
        "informative_sites_tract": in1 + in2,
        "informative_sites_flank": fl1 + fl2,
        "tract_supports": tract_winner,
        "flank_supports": flank_winner,
    }


def build_consensus(events: list[dict], min_methods: int, pvalue_max: float,
                    min_overlap: float = 0.5,
                    breakpoint_tolerance: int = 100,
                    max_members: int = 3) -> list[dict]:
    """Collapse per-method events into consensus events.

    An event is defined by its RECOMBINANT and its BREAKPOINT REGION, not by the
    parent pair. Methods legitimately disagree on parent assignment for the same
    event -- GENECONV reports a single parent where RDP reports major and minor,
    and a method may leave a parent unassigned entirely. Grouping on an exact
    (recombinant, major, minor) tuple therefore splits one real event into
    several single-method groups, and `--min-methods` then rejects all of them.
    That is how the IS76 positive control was missed.

    Events are instead grouped by recombinant and merged by breakpoint overlap
    (single linkage, `min_overlap` reciprocal). Parents are reported as the union
    of what the supporting methods assigned.
    """
    def members(ev: dict) -> set:
        return {ev["recombinant"], ev["parent_a"], ev["parent_b"]} - {""}

    def event_parents(ev: dict) -> list[str]:
        return [p for p in (ev["parent_a"], ev["parent_b"]) if p]

    kept: list[dict] = []
    for ev in events:
        # RDP reports a Bonferroni-type corrected value, not a raw probability:
        #   pvalue = G * (L/N) * binom.sf(...)        [openrdp/rdp.py:178]
        # where G is the triplet count and L/N the window-position term. The
        # product is NOT capped at 1, so values above 1 are routine (79% of RDP
        # rows on a 61-genome panel). A Bonferroni-corrected value above 1 means
        # NOT SIGNIFICANT -- the correction has consumed all the evidence -- so
        # such rows are clamped to 1 and then fail any threshold below 1.
        #
        # Do not "exempt" out-of-range values from the threshold on the grounds
        # that they are not probabilities. That inverts the meaning of the
        # statistic and admits the least significant calls as though they were
        # the most significant.
        pv = ev["pvalue"]
        if pv == pv and ev.get("pvalue_out_of_range") and pv > 1.0:
            pv = 1.0
        if pv == pv and pv > pvalue_max:
            continue
        # Rows with no parent assignment are OpenRDP bookkeeping/fragments, not
        # interpretable recombination events. Keeping them promotes singleton
        # "events" in the consensus table.
        if len(members(ev)) < 2:
            continue
        kept.append(ev)

    # GENECONV can report the two halves of one recombinant as separate
    # single-parent fragments: a short donor-like tract plus a longer
    # backbone-like block in the same recombinant. Collapse that pattern into a
    # triplet event so role assignment has both parents available. This is the
    # IS76 control pattern after rotating the circular genome so the IR is
    # internal.
    consumed: set[int] = set()
    paired: list[dict] = []

    def interval_len(ev: dict) -> int:
        return max(0, ev["end"] - ev["start"])

    def reciprocal_overlap(a: dict, b: dict) -> float:
        lo, hi = max(a["start"], b["start"]), min(a["end"], b["end"])
        inter = max(0, hi - lo)
        shorter = min(interval_len(a), interval_len(b))
        return inter / shorter if shorter else 0.0

    for i, a in enumerate(kept):
        if i in consumed or len(event_parents(a)) != 1:
            continue
        for j in range(i + 1, len(kept)):
            b = kept[j]
            if j in consumed or len(event_parents(b)) != 1:
                continue
            if a["method"] != b["method"] or a["method"] != "geneconv":
                continue
            if a["recombinant"] != b["recombinant"] or not a["recombinant"]:
                continue
            if event_parents(a)[0] == event_parents(b)[0]:
                continue
            short, long = (a, b) if interval_len(a) <= interval_len(b) else (b, a)
            if interval_len(long) < 2 * interval_len(short):
                continue
            if reciprocal_overlap(short, long) > 0.1:
                continue
            merged = dict(short)
            merged["parent_a"] = event_parents(short)[0]
            merged["parent_b"] = event_parents(long)[0]
            merged["paired_single_parent_fragments"] = True
            paired.append(merged)
            consumed.update({i, j})
            break
    kept = paired + [ev for i, ev in enumerate(kept) if i not in consumed]

    def overlaps(a: dict, b: dict) -> bool:
        """Do two per-method events describe the same recombination event?

        Reciprocal-overlap alone is the wrong test here. What methods estimate is
        a pair of BREAKPOINT POSITIONS, and they disagree on those by tens to
        hundreds of nucleotides; the interval between them is a by-product. Two
        methods can bracket the same breakpoint and still fall below any fixed
        overlap fraction when their intervals differ in length -- on the IS76
        control, GENECONV (1382-1492) and RDP (1343-1420) describe one event and
        reciprocally overlap by 0.4935.

        An event pair is therefore merged when EITHER the intervals overlap
        substantially, OR the corresponding breakpoints agree to within
        `breakpoint_tolerance` nucleotides.
        """
        lo, hi = max(a["start"], b["start"]), min(a["end"], b["end"])
        inter = max(0, hi - lo)
        shorter = min(a["end"] - a["start"], b["end"] - b["start"])
        if shorter > 0 and inter / shorter >= min_overlap:
            return True
        return (abs(a["start"] - b["start"]) <= breakpoint_tolerance
                and abs(a["end"] - b["end"]) <= breakpoint_tolerance)

    # Events are grouped by the SET OF SEQUENCES INVOLVED, not by which one
    # OpenRDP called the recombinant. That column is not stable: on the IS76
    # control, RDP names IS76 the recombinant when all six methods run and names
    # TYLCSV the recombinant when run alone -- same triplet, same breakpoints
    # 1343-1420, roles permuted (deterministic across seeds 3/7/42). Keying on it
    # therefore splits one event by method set. Two events join when they share at
    # least two sequences and their breakpoints agree, which also lets GENECONV's
    # two-sequence call merge with RDP's three-sequence call.
    # The merge is BOUNDED at `max_members` sequences. Single-linkage grouping
    # with no bound chains transitively -- A joins B, B joins C, so A, B and C
    # land in one group even when A and C share no interval. On the 61-genome
    # panel (142,985 raw events) that collapsed to 8 groups at min_methods=2,
    # the largest holding 59 of 61 sequences and spanning the whole alignment.
    # A recombination event involves a recombinant and two parents; a group of
    # 59 is not an event. The bound is the biological definition, not a tuning
    # knob. Raising it restores chaining and is only useful for testing.
    #
    # NOTE: this failure mode is invisible on the 3-sequence positive control --
    # chaining needs a dense panel to appear at all. Changes here must be checked
    # against BOTH the control and a panel-scale run.
    #
    # Candidate clusters are looked up through a pair index rather than scanned:
    # a merge requires >= 2 shared sequences, so only clusters already containing
    # one of this event's member-pairs can qualify. Without it the scan is
    # quadratic in the number of clusters (~10 min on the 61-genome panel).
    groups: list[list[dict]] = []
    group_members: list[set] = []
    pair_index: dict[frozenset, set] = {}

    def pairs_of(ms: set):
        ms = sorted(ms)
        return [frozenset((ms[i], ms[j]))
                for i in range(len(ms)) for j in range(i + 1, len(ms))]

    for ev in sorted(kept, key=lambda e: (e["start"], e["end"])):
        ev_ms = members(ev)
        candidates: set = set()
        for pr in pairs_of(ev_ms):
            candidates |= pair_index.get(pr, set())
        for gi in sorted(candidates):
            if len(group_members[gi] | ev_ms) > max_members:
                continue
            if any(overlaps(ev, m) and len(ev_ms & members(m)) >= 2
                   for m in groups[gi]):
                groups[gi].append(ev)
                group_members[gi] |= ev_ms
                for pr in pairs_of(group_members[gi]):
                    pair_index.setdefault(pr, set()).add(gi)
                break
        else:
            gi = len(groups)
            groups.append([ev])
            group_members.append(set(ev_ms))
            for pr in pairs_of(ev_ms):
                pair_index.setdefault(pr, set()).add(gi)

    consensus = []
    for evs in groups:
        methods = sorted({e["method"] for e in evs if e["method"]})
        if len(methods) < min_methods:
            continue
        # Retain OpenRDP's own labels for audit, but do not treat them as roles.
        all_members = set()
        for e in evs:
            all_members |= members(e)
        rec_labels = ";".join(sorted({e["recombinant"] for e in evs if e["recombinant"]}))
        starts = sorted(e["start"] for e in evs)
        ends = sorted(e["end"] for e in evs)
        pvals = [e["pvalue"] for e in evs
                 if e["pvalue"] == e["pvalue"] and not e.get("pvalue_out_of_range")]
        consensus.append({
            "members": ";".join(sorted(all_members)),
            "openrdp_recombinant_labels": rec_labels,
            "n_methods": len(methods),
            "methods": ",".join(methods),
            "breakpoint_start": starts[len(starts) // 2],
            "breakpoint_end": ends[len(ends) // 2],
            "start_min": starts[0],
            "end_max": ends[-1],
            "min_pvalue": min(pvals) if pvals else float("nan"),
        })
    consensus.sort(key=lambda r: (-r["n_methods"], r["min_pvalue"]))
    return consensus


def load_region_manifest(path: Path) -> dict[str, tuple[int, int]]:
    """Load {id: (start, end)} from an extract_genome_regions.py manifest."""
    spans: dict[str, tuple[int, int]] = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                spans[row["id"]] = (int(row["start"]), int(row["end"]))
            except (KeyError, ValueError):
                continue
    return spans


def annotate_breakpoints(consensus: list[dict], feature_tsv: Path | None) -> list[dict]:
    """Flag whether each consensus breakpoint falls in a feature of interest.

    `feature_tsv` is a simple TSV of: name<TAB>start<TAB>end (alignment
    coordinates). For RB work the feature of interest is the IR / CP-Rep
    junction, which is where the IS76-type breakpoint lies.
    """
    if feature_tsv is None or not feature_tsv.exists():
        return []
    feats = []
    with feature_tsv.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0].lower() in ("name", "feature"):
                continue
            try:
                feats.append((parts[0], int(parts[1]), int(parts[2])))
            except ValueError:
                continue
    out = []
    for ev in consensus:
        hit = [name for name, lo, hi in feats
               if not (ev["breakpoint_end"] < lo or ev["breakpoint_start"] > hi)]
        row = dict(ev)
        row["overlapping_features"] = ",".join(hit) if hit else "none"
        out.append(row)
    return out


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run(args) -> None:
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    aln = Path(args.alignment_fasta).resolve()
    if not aln.exists():
        sys.exit(f"Alignment not found: {aln}")

    raw_csv = outdir / "recombination_raw.csv"
    cmd = resolve_openrdp()
    methods = [m.strip().lower() for m in args.methods.split(",") if m.strip()]

    # OpenRDP accepts only A,T,G,C,-,N. Recode IUPAC ambiguity codes to N and run
    # on the sanitised copy, leaving the user's alignment untouched.
    alignment = read_alignment(aln)
    openrdp_input = outdir / "openrdp_input.fasta"
    recode = sanitise_alignment_for_openrdp(alignment, openrdp_input)
    if recode["n_sites_recoded"]:
        print(f"recoded {recode['n_sites_recoded']} IUPAC ambiguity site(s) to N "
              f"across {recode['n_records_affected']} record(s): {recode['codes']}",
              flush=True)

    if cmd is None:
        sys.exit(
            "OpenRDP is not installed.\n"
            "  pip install openrdp\n"
            "or install from https://github.com/PoonLab/OpenRDP\n"
            "Note: the 3Seq and GENECONV components are external binaries with "
            "non-commercial/academic-use licences."
        )

    # Methods are run ONE PER INVOCATION rather than all in one call, because a
    # crash in any single method aborts the whole run and loses the results of
    # every method that already completed. OpenRDP's maxchi raises
    # ZeroDivisionError on triplets that yield no valid chi-squared window
    # (openrdp/maxchi.py, `0.05/l` with l == 0), which is reachable on real
    # panels. Per-method invocation degrades to a partial result instead.
    #
    # NB: OpenRDP's -m takes nargs="+", so a REPEATED -m overwrites rather than
    # accumulates and only the last method actually runs (silently: exit 0, empty
    # output table). Methods must be passed space-delimited after a single -m.
    events: list[dict] = []
    method_status: dict[str, str] = {}
    stdout_parts, stderr_parts = [], []
    for m in methods:
        per_csv = outdir / f"raw_{m}.csv"
        full_cmd = cmd + [str(openrdp_input), "-o", str(per_csv), "-m", m]
        if args.cfg:
            full_cmd += ["-c", str(Path(args.cfg).resolve())]
        print(f"running {m}: {' '.join(full_cmd)}", flush=True)
        t0 = time.time()
        proc = subprocess.run(full_cmd, capture_output=True, text=True)
        dt = time.time() - t0
        stdout_parts.append(f"### {m}\n{proc.stdout or ''}")
        stderr_parts.append(f"### {m}\n{proc.stderr or ''}")
        if proc.returncode != 0:
            method_status[m] = f"failed(exit {proc.returncode})"
            print(f"  {m} FAILED after {dt:.0f}s -- continuing with other methods",
                  flush=True)
            continue
        found = parse_openrdp_csv(per_csv)
        events.extend(found)
        method_status[m] = "ok"
        print(f"  {m}: {len(found)} event(s) in {dt:.0f}s", flush=True)

    (outdir / "openrdp_stdout.log").write_text("\n".join(stdout_parts))
    (outdir / "openrdp_stderr.log").write_text("\n".join(stderr_parts))
    if not any(v == "ok" for v in method_status.values()):
        sys.exit(f"Every method failed; see {outdir/'openrdp_stderr.log'}")

    write_tsv(raw_csv, events,
              ["method", "start", "end", "recombinant", "parent_a", "parent_b", "pvalue"])
    write_tsv(outdir / "recombination_events.tsv", events,
              ["recombinant", "parent_a", "parent_b", "start", "end", "method", "pvalue"])

    consensus = build_consensus(events, args.min_methods, args.pvalue)

    # Positional Parent1/Parent2 do not say which parent is the backbone. Resolve
    # that from informative sites in the alignment itself.
    # `alignment` was read verbatim above (gaps preserved, pre-recoding), which is
    # what the informative-site comparison needs.
    for ev in consensus:
        ev.update(assign_roles(alignment, ev["members"].split(";"),
                               ev["start_min"], ev["end_max"]))

    cons_cols = ["recombinant_resolved", "backbone_parent", "donor_parent", "members",
                 "openrdp_recombinant_labels",
                 "n_methods", "methods", "breakpoint_start", "breakpoint_end",
                 "start_min", "end_max", "min_pvalue",
                 "informative_sites_tract", "informative_sites_flank",
                 "tract_supports", "flank_supports", "role_support"]
    write_tsv(outdir / "recombination_consensus.tsv", consensus, cons_cols)

    feature_tsv = Path(args.feature_tsv).resolve() if args.feature_tsv else None
    annotated = annotate_breakpoints(consensus, feature_tsv)
    if annotated:
        write_tsv(outdir / "breakpoint_regions.tsv", annotated, cons_cols + ["overlapping_features"])

    with (outdir / "summary.txt").open("w") as fh:
        fh.write(f"alignment_fasta\t{aln}\n")
        fh.write(f"methods\t{','.join(methods)}\n")
        fh.write(f"min_methods\t{args.min_methods}\n")
        fh.write(f"pvalue_max\t{args.pvalue}\n")
        fh.write(f"iupac_sites_recoded_to_n\t{recode['n_sites_recoded']}\n")
        fh.write(f"iupac_records_affected\t{recode['n_records_affected']}\n")
        for m in methods:
            fh.write(f"method_status.{m}\t{method_status.get(m, 'not_run')}\n")
        fh.write(f"n_events_raw\t{len(events)}\n")
        fh.write(f"n_consensus_events\t{len(consensus)}\n")
        fh.write(f"raw_csv\t{raw_csv}\n")

    print(f"raw events: {len(events)}   consensus events (>= {args.min_methods} methods): {len(consensus)}")
    print(f"outputs: {outdir}")
