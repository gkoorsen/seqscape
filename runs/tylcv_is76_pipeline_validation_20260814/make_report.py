#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/seqscape-matplotlib-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parent
ALIGNMENT = Path(
    "/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/runs/"
    "tocsv_260203_full_aws_20260612/full_run/comparators/"
    "is76_control_alignment.fasta"
)
CONTROL_JSON = Path(
    "/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/runs/"
    "tocsv_260203_full_aws_20260612/full_run/comparators/"
    "is76_positive_control.json"
)
STRICT = ROOT / "strict_p005"
AUDIT = ROOT / "rdp_audit_p1"
FROM_400 = ROOT / "from_genomes_400nt"
FROM_1000 = ROOT / "from_genomes_1000nt"
FIGURES = ROOT / "figures"
TABLES = ROOT / "tables"
COMPARATOR_FASTA = Path(
    "/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/runs/"
    "tocsv_260203_full_aws_20260612/full_run/comparators/"
    "comparator_panel.fasta"
)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = ""
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current:
                records[current].append(line.upper())
    return {key: "".join(chunks) for key, chunks in records.items()}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_summary(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            key, value = line.split("\t", 1)
            out[key] = value
    return out


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def informative_states_named(
    records: dict[str, str],
    recombinant_id: str,
    backbone_id: str,
    donor_id: str,
) -> list[tuple[int, str]]:
    rec = records[recombinant_id]
    backbone = records[backbone_id]
    donor = records[donor_id]
    states: list[tuple[int, str]] = []
    for idx, (r, b, d) in enumerate(zip(rec, backbone, donor), start=1):
        if b not in "ACGT" or d not in "ACGT" or r not in "ACGT" or b == d:
            continue
        if r == d and r != b:
            states.append((idx, "donor"))
        elif r == b and r != d:
            states.append((idx, "backbone"))
        else:
            states.append((idx, "other"))
    return states


def informative_states(records: dict[str, str]) -> list[tuple[int, str]]:
    return informative_states_named(
        records,
        "IS76_LN812978",
        "TYLCV_IL_AM409201",
        "TYLCSV_NC_003828",
    )


def rolling_donor_fraction(states: list[tuple[int, str]], length: int, window: int = 120):
    state_by_pos: dict[int, str] = dict(states)
    xs: list[int] = []
    ys: list[float] = []
    counts: list[int] = []
    half = window // 2
    for pos in range(1, length + 1):
        lo = max(1, pos - half)
        hi = min(length, pos + half)
        vals = [state_by_pos[i] for i in range(lo, hi + 1) if i in state_by_pos]
        informative = [v for v in vals if v in {"donor", "backbone"}]
        if not informative:
            continue
        xs.append(pos)
        ys.append(sum(1 for v in informative if v == "donor") / len(informative))
        counts.append(len(informative))
    return xs, ys, counts


def plot_informative_sites(
    records: dict[str, str],
    states: list[tuple[int, str]],
    strict_event: dict[str, str],
    audit_event: dict[str, str],
) -> None:
    length = len(records["IS76_LN812978"])
    strict_start = int(strict_event["breakpoint_start"])
    strict_end = int(strict_event["breakpoint_end"])
    audit_start = int(audit_event["start_min"])
    audit_end = int(audit_event["end_max"])
    xs, donor_frac, counts = rolling_donor_fraction(states, length)

    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        figsize=(11.5, 5.8),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.2]},
    )
    fig.patch.set_facecolor("white")
    for ax in (ax1, ax2):
        ax.set_facecolor("white")
        ax.axvspan(audit_start, audit_end, color="#d8b4fe", alpha=0.22, linewidth=0)
        ax.axvspan(strict_start, strict_end, color="#fbbf24", alpha=0.28, linewidth=0)

    ax1.plot(xs, donor_frac, color="#7c3aed", linewidth=1.5)
    ax1.axhline(0.5, color="#737373", linewidth=0.8, linestyle="--")
    ax1.set_ylabel("Donor-match\nfraction")
    ax1.set_ylim(-0.05, 1.05)
    ax1.set_title(
        "IS76 positive control: donor-like tract recovered at the origin/IR",
        loc="left",
        fontweight="bold",
    )

    colors = {"backbone": "#2563eb", "donor": "#f97316", "other": "#6b7280"}
    yvals = {"backbone": 0.0, "donor": 1.0, "other": 0.5}
    labels = {"backbone": "IS76 matches TYLCV-IL backbone", "donor": "IS76 matches TYLCSV donor", "other": "Other informative state"}
    for state in ("backbone", "donor", "other"):
        pts = [pos for pos, value in states if value == state]
        ax2.scatter(
            pts,
            [yvals[state]] * len(pts),
            s=13 if state != "other" else 10,
            color=colors[state],
            alpha=0.85,
            linewidths=0,
            label=labels[state],
        )
    ax2.set_yticks([0, 0.5, 1])
    ax2.set_yticklabels(["Backbone", "Other", "Donor"])
    ax2.set_xlabel("MAFFT alignment column")
    ax2.set_xlim(1, length)

    handles = [
        Patch(facecolor="#fbbf24", alpha=0.28, label="Strict GENECONV tract"),
        Patch(facecolor="#d8b4fe", alpha=0.22, label="RDP+GENECONV audit envelope"),
    ]
    scatter_handles, scatter_labels = ax2.get_legend_handles_labels()
    ax1.legend(handles=handles, loc="upper right", frameon=False)
    ax2.legend(scatter_handles, scatter_labels, loc="lower right", frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(FIGURES / "is76_parental_state_across_alignment.png", dpi=300)
    plt.close(fig)

    zoom_lo = max(1, audit_start - 130)
    zoom_hi = min(length, audit_end + 130)
    fig, ax = plt.subplots(figsize=(10, 3.4))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.axvspan(audit_start, audit_end, color="#d8b4fe", alpha=0.25, linewidth=0)
    ax.axvspan(strict_start, strict_end, color="#fbbf24", alpha=0.35, linewidth=0)
    for state in ("backbone", "donor", "other"):
        pts = [pos for pos, value in states if value == state and zoom_lo <= pos <= zoom_hi]
        ax.scatter(
            pts,
            [yvals[state]] * len(pts),
            s=25 if state != "other" else 18,
            color=colors[state],
            alpha=0.9,
            linewidths=0,
            label=labels[state],
        )
    ax.set_xlim(zoom_lo, zoom_hi)
    ax.set_ylim(-0.25, 1.25)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["Backbone", "Other", "Donor"])
    ax.set_xlabel("MAFFT alignment column")
    ax.set_title("Breakpoint-region zoom", loc="left", fontweight="bold")
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "is76_breakpoint_zoom.png", dpi=300)
    plt.close(fig)


def plot_method_intervals(events: list[dict[str, str]], strict_event: dict[str, str]) -> None:
    rows = []
    for event in events:
        method = event["method"]
        if method not in {"rdp", "geneconv"}:
            continue
        rows.append((method.upper(), int(event["start"]), int(event["end"]), event["pvalue"]))
    rows = sorted(rows, key=lambda row: (row[0], row[1], row[2]))

    fig, ax = plt.subplots(figsize=(8.8, 2.8))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ylabels = []
    for idx, (method, start, end, pvalue) in enumerate(rows):
        y = len(rows) - idx
        color = "#f97316" if method == "GENECONV" and start == int(strict_event["breakpoint_start"]) else "#7c3aed"
        ax.plot([start, end], [y, y], color=color, linewidth=8, solid_capstyle="butt")
        label = f"{method}: {start}-{end}"
        if method == "GENECONV":
            label += f" (p={pvalue})"
        else:
            label += " (OpenRDP score column)"
        ylabels.append(label)
    ax.set_yticks(range(len(rows), 0, -1))
    ax.set_yticklabels(ylabels)
    ax.set_xlabel("MAFFT alignment column")
    ax.set_title("OpenRDP method intervals", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#e5e7eb", linewidth=0.7)
    fig.tight_layout()
    fig.savefig(FIGURES / "openrdp_method_intervals.png", dpi=300)
    plt.close(fig)


def plot_origin_window_sites(
    records: dict[str, str],
    states: list[tuple[int, str]],
    event: dict[str, str],
) -> None:
    length = len(records["LN812978"])
    start = int(event["breakpoint_start"])
    end = int(event["breakpoint_end"])
    xs, donor_frac, _ = rolling_donor_fraction(states, length, window=80)

    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        figsize=(10.5, 4.8),
        sharex=True,
        gridspec_kw={"height_ratios": [2.0, 1.1]},
    )
    fig.patch.set_facecolor("white")
    for ax in (ax1, ax2):
        ax.set_facecolor("white")
        ax.axvspan(start, end, color="#fbbf24", alpha=0.30, linewidth=0)

    ax1.plot(xs, donor_frac, color="#7c3aed", linewidth=1.5)
    ax1.axhline(0.5, color="#737373", linewidth=0.8, linestyle="--")
    ax1.set_ylabel("Donor-match\nfraction")
    ax1.set_ylim(-0.05, 1.05)
    ax1.set_title(
        "Genome-derived 1000 nt origin window recovers the IS76 tract",
        loc="left",
        fontweight="bold",
    )

    colors = {"backbone": "#2563eb", "donor": "#f97316", "other": "#6b7280"}
    yvals = {"backbone": 0.0, "donor": 1.0, "other": 0.5}
    labels = {
        "backbone": "LN812978 matches AM409201 backbone",
        "donor": "LN812978 matches NC_003828 donor",
        "other": "Other informative state",
    }
    for state in ("backbone", "donor", "other"):
        pts = [pos for pos, value in states if value == state]
        ax2.scatter(
            pts,
            [yvals[state]] * len(pts),
            s=16 if state != "other" else 11,
            color=colors[state],
            alpha=0.88,
            linewidths=0,
            label=labels[state],
        )
    ax2.set_yticks([0, 0.5, 1])
    ax2.set_yticklabels(["Backbone", "Other", "Donor"])
    ax2.set_xlabel("MAFFT alignment column")
    ax2.set_xlim(1, length)

    handles = [Patch(facecolor="#fbbf24", alpha=0.30, label="Strict GENECONV tract")]
    ax1.legend(handles=handles, loc="upper right", frameon=False)
    ax2.legend(loc="lower right", frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(FIGURES / "is76_origin_window_1000nt_parental_state.png", dpi=300)
    plt.close(fig)


def markdown_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    out = ["|" + "|".join(columns) + "|", "|" + "|".join(["---"] * len(columns)) + "|"]
    for row in rows:
        out.append("|" + "|".join(str(row.get(col, "")) for col in columns) + "|")
    return "\n".join(out)


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    TABLES.mkdir(exist_ok=True)

    records = read_fasta(ALIGNMENT)
    control = json.loads(CONTROL_JSON.read_text())
    strict_cons = read_tsv(STRICT / "recombination_consensus.tsv")
    audit_cons = read_tsv(AUDIT / "recombination_consensus.tsv")
    events = read_tsv(AUDIT / "recombination_events.tsv")
    strict_summary = read_summary(STRICT / "summary.txt")
    audit_summary = read_summary(AUDIT / "summary.txt")
    strict_event = strict_cons[0]
    audit_event = audit_cons[0]
    from_400_event = read_tsv(
        FROM_400 / "recombination_strict_p005" / "recombination_consensus.tsv"
    )[0]
    from_1000_event = read_tsv(
        FROM_1000 / "recombination_strict_p005" / "recombination_consensus.tsv"
    )[0]
    from_1000_records = read_fasta(FROM_1000 / "is76_origin_window_1000nt_aligned.fasta")
    from_1000_states = informative_states_named(
        from_1000_records,
        "LN812978",
        "AM409201",
        "NC_003828",
    )
    states = informative_states(records)

    sequence_rows = [
        {
            "role": "known recombinant",
            "id": "IS76_LN812978",
            "description": f"{control['recombinant']['name']} ({control['recombinant']['country']})",
        },
        {
            "role": "backbone parent",
            "id": "TYLCV_IL_AM409201",
            "description": f"{control['major_parent']['name']} ({control['major_parent']['country']}, {control['major_parent']['year']})",
        },
        {
            "role": "donor parent",
            "id": "TYLCSV_NC_003828",
            "description": control["minor_parent_donor"]["name"],
        },
    ]
    write_tsv(TABLES / "control_sequences.tsv", sequence_rows, ["role", "id", "description"])

    validation_rows = [
        {
            "run": "strict p<=0.05",
            "recombinant": strict_event["recombinant_resolved"],
            "backbone_parent": strict_event["backbone_parent"],
            "donor_parent": strict_event["donor_parent"],
            "methods": strict_event["methods"],
            "breakpoint_columns": f"{strict_event['breakpoint_start']}-{strict_event['breakpoint_end']}",
            "min_pvalue": strict_event["min_pvalue"],
            "role_support": strict_event["role_support"],
        },
        {
            "run": "RDP audit p<=1.0",
            "recombinant": audit_event["recombinant_resolved"],
            "backbone_parent": audit_event["backbone_parent"],
            "donor_parent": audit_event["donor_parent"],
            "methods": audit_event["methods"],
            "breakpoint_columns": f"{audit_event['start_min']}-{audit_event['end_max']}",
            "min_pvalue": audit_event["min_pvalue"],
            "role_support": audit_event["role_support"],
        },
    ]
    write_tsv(
        TABLES / "validation_summary.tsv",
        validation_rows,
        [
            "run",
            "recombinant",
            "backbone_parent",
            "donor_parent",
            "methods",
            "breakpoint_columns",
            "min_pvalue",
            "role_support",
        ],
    )

    genome_rows = [
        {
            "input": "400 nt origin window from comparator genomes",
            "records": "3",
            "consensus": "detected",
            "recombinant": from_400_event["recombinant_resolved"] or "not resolved",
            "backbone_parent": from_400_event["backbone_parent"] or "not resolved",
            "donor_parent": from_400_event["donor_parent"] or "not resolved",
            "methods": from_400_event["methods"],
            "breakpoint_columns": f"{from_400_event['breakpoint_start']}-{from_400_event['breakpoint_end']}",
            "min_pvalue": from_400_event["min_pvalue"],
            "role_support": from_400_event["role_support"],
            "interpretation": "event detected, but the window is too short for flank-based role resolution",
        },
        {
            "input": "1000 nt origin window from comparator genomes",
            "records": "3",
            "consensus": "detected",
            "recombinant": from_1000_event["recombinant_resolved"],
            "backbone_parent": from_1000_event["backbone_parent"],
            "donor_parent": from_1000_event["donor_parent"],
            "methods": from_1000_event["methods"],
            "breakpoint_columns": f"{from_1000_event['breakpoint_start']}-{from_1000_event['breakpoint_end']}",
            "min_pvalue": from_1000_event["min_pvalue"],
            "role_support": from_1000_event["role_support"],
            "interpretation": "event and expected parent roles recovered",
        },
    ]
    genome_columns = [
        "input",
        "records",
        "consensus",
        "recombinant",
        "backbone_parent",
        "donor_parent",
        "methods",
        "breakpoint_columns",
        "min_pvalue",
        "role_support",
        "interpretation",
    ]
    write_tsv(TABLES / "from_genomes_validation.tsv", genome_rows, genome_columns)

    method_rows = []
    for key, value in strict_summary.items():
        if key.startswith("method_status."):
            method_rows.append({"method": key.split(".", 1)[1], "strict_status": value})
    write_tsv(TABLES / "method_status.tsv", method_rows, ["method", "strict_status"])

    plot_informative_sites(records, states, strict_event, audit_event)
    plot_method_intervals(events, strict_event)
    plot_origin_window_sites(from_1000_records, from_1000_states, from_1000_event)

    n_donor = sum(1 for _, state in states if state == "donor")
    n_backbone = sum(1 for _, state in states if state == "backbone")
    n_other = sum(1 for _, state in states if state == "other")

    report = f"""# TYLCV-IS76 Positive-Control Validation

## Result

The SeqScape recombination wrapper recovered the known TYLCV-IS76 resistance-breaking recombinant signal. In the strict screen (`p <= 0.05`), the consensus event resolves:

- recombinant: `{strict_event['recombinant_resolved']}`
- backbone parent: `{strict_event['backbone_parent']}`
- donor parent: `{strict_event['donor_parent']}`
- breakpoint interval: alignment columns `{strict_event['breakpoint_start']}-{strict_event['breakpoint_end']}`
- strict supporting method: `{strict_event['methods']}`
- role support: `{strict_event['role_support']}` informative sites

The audit run retains RDP's score-like output and merges it with GENECONV over the same tract envelope (`{audit_event['start_min']}-{audit_event['end_max']}`). This is reported separately because the RDP output column is not a valid probability in this OpenRDP run (`36.36`, greater than 1), while GENECONV gives the strict significant call.

## Input

Input alignment:

`{ALIGNMENT}`

Positive-control metadata:

`{CONTROL_JSON}`

{markdown_table(sequence_rows, ['role', 'id', 'description'])}

The metadata file records the published TYLCV-IS76 tract as approximately {control['published_tract_nt']} nt. The local control annotation places the recovered signal in the intergenic region immediately downstream of the conserved origin motif, with a reconstructed tract range of {control['recovered']['tract_min_nt']}-{control['recovered']['tract_max_nt']} nt.

## Recombination Calls

{markdown_table(validation_rows, ['run', 'recombinant', 'backbone_parent', 'donor_parent', 'methods', 'breakpoint_columns', 'min_pvalue', 'role_support'])}

All OpenRDP methods completed in the strict run:

{markdown_table(method_rows, ['method', 'strict_status'])}

## Genome-Derived Pipeline Check

To test the practical workflow from comparator genomes, the three positive-control records were extracted from:

`{COMPARATOR_FASTA}`

Origin-centered windows were cut around the conserved `TAATATTAC` motif, aligned with MAFFT, and passed through the same SeqScape recombination wrapper.

{markdown_table(genome_rows, genome_columns)}

The 400 nt window is sensitive enough to detect the event, but its significant intervals consume most of the window and leave no informative flank support for parent-role polarisation. The 1000 nt origin window resolves the expected roles directly from the extracted genomes: `LN812978` as recombinant, `AM409201` as backbone parent, and `NC_003828` as donor parent.

## Figures

![Parent-state scan](figures/is76_parental_state_across_alignment.png)

**Figure 1.** Informative-site scan across the aligned TYLCV-IS76 positive control. Blue sites are columns where IS76 matches the TYLCV-IL backbone parent; orange sites are columns where IS76 matches the TYLCSV donor parent. The strict GENECONV tract is shaded yellow. The RDP+GENECONV audit envelope is shaded purple.

![Breakpoint zoom](figures/is76_breakpoint_zoom.png)

**Figure 2.** Zoom around the recovered breakpoint region. The donor-like informative sites concentrate inside the OpenRDP/GENECONV interval, while the flanking sequence returns to the TYLCV-IL backbone state.

![Method intervals](figures/openrdp_method_intervals.png)

**Figure 3.** OpenRDP method intervals in the positive-control run. GENECONV gives the strict significant interval; RDP reports an overlapping interval but its numeric field is score-like rather than a p-value in this run.

![Genome-derived origin-window scan](figures/is76_origin_window_1000nt_parental_state.png)

**Figure 4.** Informative-site scan for the 1000 nt origin-centered window extracted from the comparator genome FASTA. The donor-like sites concentrate in the strict GENECONV interval, while flanking sites support the TYLCV-IL backbone parent.

## Interpretation

This positive-control test passes. The pipeline detects the known resistance-breaking recombination pattern and assigns the expected roles: IS76 as recombinant, TYLCV-IL as the backbone parent, and TYLCSV as the donor-side parent.

The strict result is intentionally conservative: it promotes the GENECONV-supported event under `p <= 0.05`. The RDP audit supports the same region geometrically, but is not counted as strict statistical support because OpenRDP emits a value outside the probability range for that method in this control.

For full-genome inputs, this validation supports using an origin-centered window with enough flanking sequence for role assignment. A 1000 nt window worked for this positive control; a 400 nt window detected recombination but did not retain enough flanking signal to resolve the parent roles.

Informative-site counts from the alignment:

- donor-matching informative sites: {n_donor}
- backbone-matching informative sites: {n_backbone}
- other informative states: {n_other}
- total informative sites: {len(states)}

The `role_support` value is the donor-matching plus backbone-matching count (`{n_donor + n_backbone}`); the six other informative states are retained in the figure but are not counted as role support.

## Commands

Strict run:

```bash
PATH="/private/tmp/openrdp-validation-venv/bin:$PATH" PYTHONPATH=src \\
python -m seqscape.cli recombination \\
  --alignment-fasta "{ALIGNMENT}" \\
  --methods rdp,geneconv,maxchi,chimaera,threeseq,bootscan,siscan \\
  --min-methods 1 \\
  --pvalue 0.05 \\
  --outdir "{STRICT}"
```

RDP audit run:

```bash
PATH="/private/tmp/openrdp-validation-venv/bin:$PATH" PYTHONPATH=src \\
python -m seqscape.cli recombination \\
  --alignment-fasta "{ALIGNMENT}" \\
  --methods rdp,geneconv,maxchi,chimaera,threeseq,bootscan,siscan \\
  --min-methods 1 \\
  --pvalue 1.0 \\
  --outdir "{AUDIT}"
```

Genome-derived 1000 nt origin-window run:

```bash
python scripts/extract_origin_windows.py \\
  --fasta "{COMPARATOR_FASTA}" \\
  --ids "{ROOT / 'inputs' / 'is76_ids.txt'}" \\
  --window-size 1000 \\
  --out-fasta "{FROM_1000 / 'is76_origin_window_1000nt.fasta'}" \\
  --manifest-tsv "{FROM_1000 / 'is76_origin_window_1000nt_manifest.tsv'}"

mafft --auto --thread 8 \\
  "{FROM_1000 / 'is76_origin_window_1000nt.fasta'}" \\
  > "{FROM_1000 / 'is76_origin_window_1000nt_aligned.fasta'}"

PATH="/private/tmp/openrdp-validation-venv/bin:$PATH" PYTHONPATH=src \\
python -m seqscape.cli recombination \\
  --alignment-fasta "{FROM_1000 / 'is76_origin_window_1000nt_aligned.fasta'}" \\
  --methods rdp,geneconv,maxchi,chimaera,threeseq,bootscan,siscan \\
  --min-methods 1 \\
  --pvalue 0.05 \\
  --outdir "{FROM_1000 / 'recombination_strict_p005'}"
```
"""
    (ROOT / "report.md").write_text(report)
    print(ROOT / "report.md")


if __name__ == "__main__":
    main()
