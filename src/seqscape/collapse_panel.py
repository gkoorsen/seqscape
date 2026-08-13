"""Redundancy-collapse panel selection -- the recombination-screen fork.

Sibling of `sampling.py`. Both consume the SAME UMAP Explorer state and emit the
SAME two files (`representative_panel.fasta` + `representative_panel_manifest.tsv`),
so `align-panel`, `distance-explorer` and `export-clusters` accept either without
modification. The scientist inspects composition in `umap-explorer`, then forks:

    umap-explorer ──┬── sample-panel    (proportional; describe diversity)
                    └── collapse-panel  (dereplicated; screen for recombination)

Why a separate selector rather than a flag on `sample-panel`: the two have
opposite objectives, and the objective is not a parameter.

`sample-panel` allocates a fixed budget across Leiden clusters in proportion to
cluster size, and is scored on `pearson_cluster_proportion_preservation` -- it is
built to REPRODUCE the composition of the input. That is correct for a
first-glance overview panel.

A recombination screen needs the opposite. Public sequence databases reflect
outbreak resequencing effort, not diversity, so composition is an artefact. A
proportional panel reproduces the artefact: measured on the 1,792-genome ToCSV
comparator pull, a budget of 173 retains 2 of 31 focal ToCSV genomes, because all
31 sit in one 46-genome Leiden cluster whose proportional quota is 4 slots shared
with 15 non-focal members. It also leaves 73% of retained pairs above the collapse
threshold, since a quota larger than the number of distinct lineages in a cluster
cannot help but keep near-duplicates.

Recombination detection is additionally *sensitive* to this in a way that
diversity description is not. Triplet-based methods (RDP, GENECONV, MaxChi,
BootScan, SiScan) test parent-parent-child triplets. Near-identical duplicates
inflate the triplet count without adding donor diversity, and multiple-testing
correction across those triplets then costs power at the events that matter.
One representative per distinct lineage is the standard recommendation for RDP-class
input (Martin et al. 2015, Virus Evolution 1:vev003).

Selection rule here, in strict priority order:

  1. references (`item_class == "reference"`)  -- always retained, as in sample-panel
  2. focal ids (`--focal-ids`)                 -- always retained, never collapsed
  3. everything else                           -- greedy collapse at a fixed
                                                  identity threshold, no quota

There is no per-cluster budget. Panel size is an OUTPUT, set by how many distinct
lineages the data actually contains -- which is the point.

Identity is estimated by k-mer Jaccard rather than alignment, so selection stays
alignment-free and O(n^2) in cheap set operations; `align-panel` then computes
exact identities on the (much smaller) selected panel. The Jaccard->identity
conversion is the standard Mash-style relation (Ondov et al. 2016, Genome Biol
17:132).
"""
from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .io_utils import load_fasta_records, read_assignments, read_coords, write_fasta, write_tsv
from .umap_explorer import format_resolution
from .validation import compression_ratio, full_pairwise_count, panel_fraction_pct

MANIFEST_COLUMNS = [
    "id", "label", "item_class", "source_cluster", "selection_rank",
    "selection_mode", "length_bp", "af_umap1", "af_umap2", "centroid_distance",
    "collapsed_members", "represents",
]


def kmer_set(seq: str, k: int) -> set[str]:
    s = seq.upper()
    return {s[i:i + k] for i in range(len(s) - k + 1) if "N" not in s[i:i + k]}


def jaccard(a: set[str], b: set[str]) -> float:
    if not a or not b:
        return 0.0
    inter = len(a & b)
    return inter / (len(a) + len(b) - inter)


def jaccard_to_identity(j: float, k: int) -> float:
    """Mash-style conversion: j -> approximate nucleotide identity (percent)."""
    if j <= 0.0:
        return 0.0
    if j >= 1.0:
        return 100.0
    return 100.0 * (1.0 + (1.0 / k) * math.log(2.0 * j / (1.0 + j)))


def collapse(
    ordered_ids: list[str],
    kmers: dict[str, set[str]],
    threshold: float,
    k: int,
    protected: set[str],
) -> tuple[list[str], dict[str, list[str]]]:
    """Greedy single-linkage-style collapse in priority order.

    `ordered_ids` is consumed in order, so callers control which genome becomes
    the representative of its group. Protected ids are always retained and never
    absorb-or-are-absorbed, so a focal genome can never be dropped, and can never
    silently stand in for an unrelated genome.
    """
    kept: list[str] = []
    members: dict[str, list[str]] = defaultdict(list)
    for sid in ordered_ids:
        if sid in protected:
            kept.append(sid)
            continue
        absorbed_by = None
        for rep in kept:
            if rep in protected:
                continue
            if jaccard_to_identity(jaccard(kmers[sid], kmers[rep]), k) >= threshold:
                absorbed_by = rep
                break
        if absorbed_by is None:
            kept.append(sid)
        else:
            members[absorbed_by].append(sid)
    return kept, members


def run(args) -> None:
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    chosen_slug = (
        f"k{args.chosen_kmer}_n{args.chosen_neighbors}_d{str(args.chosen_min_dist).replace('.', 'p')}"
        f"_r{format_resolution(args.chosen_leiden_resolution)}"
    )
    state_dir = Path(args.umap_explorer_dir).resolve() / "states" / chosen_slug
    assignments_tsv = state_dir / "assignments.tsv"
    coords_csv = state_dir / "coords.csv"
    if not assignments_tsv.is_file() or not coords_csv.is_file():
        legacy_slug = f"k{args.chosen_kmer}_n{args.chosen_neighbors}_d{str(args.chosen_min_dist).replace('.', 'p')}"
        legacy_state_dir = Path(args.umap_explorer_dir).resolve() / "states" / legacy_slug
        if (legacy_state_dir / "assignments.tsv").is_file() and (legacy_state_dir / "coords.csv").is_file():
            chosen_slug = legacy_slug
            state_dir = legacy_state_dir
            assignments_tsv = state_dir / "assignments.tsv"
            coords_csv = state_dir / "coords.csv"
        else:
            raise RuntimeError(f"Chosen AF state not found: {state_dir}")

    genome_records = load_fasta_records(Path(args.input_fasta).resolve())
    ref_records = load_fasta_records(Path(args.reference_fasta).resolve()) if args.reference_fasta else {}
    assignments = read_assignments(assignments_tsv)
    coords_map = read_coords(coords_csv)

    focal: set[str] = set()
    if args.focal_ids:
        focal_path = Path(args.focal_ids).resolve()
        if not focal_path.is_file():
            raise RuntimeError(f"--focal-ids file not found: {focal_path}")
        focal = {ln.strip() for ln in focal_path.read_text().splitlines()
                 if ln.strip() and not ln.startswith("#")}

    genome_rows = [r for r in assignments if r["item_class"] == "genome"]
    ref_rows = [r for r in assignments if r["item_class"] == "reference"]
    known = {r["id"] for r in genome_rows}
    missing_focal = sorted(focal - known)
    if missing_focal:
        raise RuntimeError(
            f"{len(missing_focal)} --focal-ids not present as genomes in the chosen state, "
            f"e.g. {missing_focal[:5]}"
        )

    k = args.jaccard_k
    kmers: dict[str, set[str]] = {}
    for r in genome_rows:
        rec = genome_records.get(r["id"])
        if rec is None:
            raise RuntimeError(f"Genome {r['id']} in assignments but absent from --input-fasta")
        kmers[r["id"]] = kmer_set(str(rec.seq), k)

    # Priority order: focal first (protected), then remaining genomes. Within the
    # non-focal set, order by descending ungapped length so the most complete
    # genome becomes its group's representative.
    non_focal = [r["id"] for r in genome_rows if r["id"] not in focal]
    non_focal.sort(key=lambda sid: (-len(str(genome_records[sid].seq).replace("-", "")), sid))
    ordered = sorted(focal) + non_focal

    kept, members = collapse(ordered, kmers, args.identity_threshold, k, protected=focal)

    id_to_row = {r["id"]: r for r in genome_rows}
    panel_items: list[dict] = []
    panel_records: list[SeqRecord] = []

    for rank, sid in enumerate(kept, start=1):
        row = id_to_row[sid]
        rec = genome_records[sid]
        seq = str(rec.seq).upper()
        x, y = coords_map.get(sid, (float("nan"), float("nan")))
        grp = members.get(sid, [])
        panel_items.append({
            "id": sid,
            "label": row["label"],
            "item_class": "genome",
            "source_cluster": row["cluster"],
            "selection_rank": rank,
            "selection_mode": "focal" if sid in focal else "collapse_representative",
            "length_bp": len(seq),
            "af_umap1": f"{x:.6f}", "af_umap2": f"{y:.6f}",
            "centroid_distance": "",
            "collapsed_members": len(grp),
            "represents": ";".join(sorted(grp)),
        })
        panel_records.append(SeqRecord(Seq(seq), id=sid, name="", description=""))

    for row in ref_rows:
        rec = ref_records.get(row["id"]) or genome_records.get(row["id"])
        if rec is None:
            raise RuntimeError(f"Reference {row['id']} not found in --reference-fasta")
        seq = str(rec.seq).upper()
        x, y = coords_map.get(row["id"], (float("nan"), float("nan")))
        panel_items.append({
            "id": row["id"], "label": row["label"], "item_class": "reference",
            "source_cluster": row["cluster"], "selection_rank": 0,
            "selection_mode": "reference", "length_bp": len(seq),
            "af_umap1": f"{x:.6f}", "af_umap2": f"{y:.6f}",
            "centroid_distance": "", "collapsed_members": 0, "represents": "",
        })
        panel_records.append(SeqRecord(Seq(seq), id=row["id"], name="", description=""))

    manifest_tsv = outdir / "representative_panel_manifest.tsv"
    panel_fasta = outdir / "representative_panel.fasta"
    groups_tsv = outdir / "collapse_groups.tsv"
    summary_txt = outdir / "summary.txt"

    write_fasta(panel_records, panel_fasta)
    write_tsv(
        manifest_tsv,
        sorted(panel_items, key=lambda r: (r["item_class"] != "reference",
                                           r["selection_mode"] != "focal",
                                           int(r["selection_rank"]))),
        MANIFEST_COLUMNS,
    )
    write_tsv(
        groups_tsv,
        [{"representative": rep, "n_absorbed": len(m), "absorbed_ids": ";".join(sorted(m))}
         for rep, m in sorted(members.items(), key=lambda kv: -len(kv[1]))],
        ["representative", "n_absorbed", "absorbed_ids"],
    )

    # Per-cluster retention, so the fork's effect on the focal cluster is auditable
    per_cluster: dict[str, dict] = defaultdict(lambda: {"genomes_in_cluster": 0, "retained": 0, "focal": 0})
    for r in genome_rows:
        per_cluster[r["cluster"]]["genomes_in_cluster"] += 1
        if r["id"] in focal:
            per_cluster[r["cluster"]]["focal"] += 1
    for sid in kept:
        per_cluster[id_to_row[sid]["cluster"]]["retained"] += 1
    write_tsv(
        outdir / "source_cluster_retention.tsv",
        [{"cluster": c, **v} for c, v in sorted(per_cluster.items())],
        ["cluster", "genomes_in_cluster", "retained", "focal"],
    )

    n_genomes = len(genome_rows)
    panel_total = len(panel_items)
    metrics = {
        "N": n_genomes,
        "panel_size_total": panel_total,
        "panel_size_genomes": len(kept),
        "panel_size_references": len(ref_rows),
        "focal_supplied": len(focal),
        "focal_retained": sum(1 for sid in kept if sid in focal),
        "identity_threshold": args.identity_threshold,
        "jaccard_k": k,
        "genomes_absorbed": sum(len(m) for m in members.values()),
        "panel_fraction_pct": panel_fraction_pct(panel_total, n_genomes),
        "full_pairwise_count": full_pairwise_count(n_genomes),
        "compression_ratio_vs_panel": compression_ratio(n_genomes, panel_total),
        "full_triplet_count": math.comb(n_genomes, 3) if n_genomes >= 3 else 0,
        "panel_triplet_count": math.comb(panel_total, 3) if panel_total >= 3 else 0,
    }
    (outdir / "validation_metrics.json").write_text(json.dumps(metrics, indent=2))

    with summary_txt.open("w") as fh:
        fh.write(f"input_fasta\t{Path(args.input_fasta).resolve()}\n")
        fh.write(f"umap_explorer_dir\t{Path(args.umap_explorer_dir).resolve()}\n")
        fh.write(f"chosen_state\t{chosen_slug}\n")
        fh.write(f"identity_threshold\t{args.identity_threshold}\n")
        fh.write(f"jaccard_k\t{k}\n")
        fh.write(f"focal_supplied\t{metrics['focal_supplied']}\n")
        fh.write(f"focal_retained\t{metrics['focal_retained']}\n")
        fh.write(f"panel_size_genomes\t{len(kept)}\n")
        fh.write(f"panel_size_references\t{len(ref_rows)}\n")
        fh.write(f"genomes_absorbed\t{metrics['genomes_absorbed']}\n")
        fh.write(f"compression_ratio_vs_panel\t{metrics['compression_ratio_vs_panel']:.6f}\n")
        fh.write(f"full_triplet_count\t{metrics['full_triplet_count']}\n")
        fh.write(f"panel_triplet_count\t{metrics['panel_triplet_count']}\n")
        fh.write("selection_scheme\tall_references_kept; focal ids protected; "
                 "remaining genomes collapsed greedily at a fixed k-mer identity "
                 "threshold; NO per-cluster quota\n")
        fh.write(f"representative_panel_fasta\t{panel_fasta}\n")
        fh.write(f"representative_panel_manifest\t{manifest_tsv}\n")
        fh.write(f"collapse_groups\t{groups_tsv}\n")

    print(f"genomes {n_genomes} -> {len(kept)} representatives "
          f"(+{len(ref_rows)} references) | focal {metrics['focal_retained']}/{metrics['focal_supplied']} "
          f"| absorbed {metrics['genomes_absorbed']} at >={args.identity_threshold}% identity")
