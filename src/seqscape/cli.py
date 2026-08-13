from __future__ import annotations

import argparse


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Staged SeqScape pipeline: build the UMAP Explorer, sample a representative "
            "panel, compute pairwise alignments, then build the Distance Explorer."
        )
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_af = sub.add_parser("umap-explorer", help="Build the alignment-free UMAP+Leiden parameter explorer")
    ap_af.add_argument("--input-fasta", required=True)
    ap_af.add_argument("--reference-fasta", required=True)
    ap_af.add_argument("--cg-map-tsv", default="")
    ap_af.add_argument("--kmer-values", default="5")
    ap_af.add_argument("--neighbors-values", default="50,100,150,200")
    ap_af.add_argument("--min-dist-values", default="0.05,0.1,0.3,0.5,0.8")
    ap_af.add_argument("--seed", type=int, default=42)
    ap_af.add_argument("--leiden-neighbor-k", type=int, default=300)
    ap_af.add_argument("--leiden-resolution", type=float, default=1e-6)
    ap_af.add_argument("--leiden-resolution-values", default="")
    ap_af.add_argument("--default-kmer", type=int, default=5)
    ap_af.add_argument("--default-neighbors", type=int, default=100)
    ap_af.add_argument("--default-min-dist", type=float, default=0.3)
    ap_af.add_argument("--default-leiden-resolution", type=float, default=None)
    ap_af.add_argument("--validation-subsample-size", type=int, default=1000)
    ap_af.add_argument("--validation-replicates", type=int, default=1)
    ap_af.add_argument("--validation-neighbors", default="10,30,100")
    ap_af.add_argument("--outdir", required=True)

    ap_sample = sub.add_parser("sample-panel", help="Sample a proportional representative panel from a chosen UMAP Explorer state")
    ap_sample.add_argument("--input-fasta", required=True)
    ap_sample.add_argument("--reference-fasta", required=True)
    ap_sample.add_argument("--umap-explorer-dir", required=True)
    ap_sample.add_argument("--chosen-kmer", type=int, required=True)
    ap_sample.add_argument("--chosen-neighbors", type=int, required=True)
    ap_sample.add_argument("--chosen-min-dist", type=float, required=True)
    ap_sample.add_argument("--chosen-leiden-resolution", type=float, default=1e-6)
    ap_sample.add_argument("--sample-size", type=int, required=True)
    ap_sample.add_argument("--outdir", required=True)

    ap_collapse = sub.add_parser(
        "collapse-panel",
        help="Dereplicated panel from a chosen UMAP Explorer state (recombination-screen fork of sample-panel)",
    )
    ap_collapse.add_argument("--input-fasta", required=True)
    ap_collapse.add_argument("--reference-fasta", default="")
    ap_collapse.add_argument("--umap-explorer-dir", required=True)
    ap_collapse.add_argument("--chosen-kmer", type=int, required=True)
    ap_collapse.add_argument("--chosen-neighbors", type=int, required=True)
    ap_collapse.add_argument("--chosen-min-dist", type=float, required=True)
    ap_collapse.add_argument("--chosen-leiden-resolution", type=float, default=1e-6)
    ap_collapse.add_argument("--identity-threshold", type=float, default=98.1,
                             help="Collapse genomes at or above this %% nucleotide identity "
                                  "(default: 98.1, the ICTV begomovirus strain demarcation level)")
    ap_collapse.add_argument("--jaccard-k", type=int, default=15,
                             help="k-mer size for the alignment-free identity estimate (default: 15)")
    ap_collapse.add_argument("--focal-ids", default="",
                             help="Optional file of one id per line, never collapsed and never dropped "
                                  "(e.g. your ToCSV candidates)")
    ap_collapse.add_argument("--outdir", required=True)

    ap_align = sub.add_parser("align-panel", help="Run all-vs-all identities on the sampled panel")
    ap_align.add_argument("--panel-fasta", required=True)
    ap_align.add_argument("--manifest-tsv", required=True)
    ap_align.add_argument("--aligner", choices=["pairwisealigner", "mafft", "muscle"], default="muscle")
    ap_align.add_argument("--mafft", default="mafft")
    ap_align.add_argument("--muscle", default="muscle")
    ap_align.add_argument("--jobs", type=int, default=8)
    ap_align.add_argument("--progress-every", type=int, default=500)
    ap_align.add_argument("--outdir", required=True)

    ap_review = sub.add_parser("distance-explorer", help="Build the post-alignment distance explorer with threshold clustering, reference scatter, and heatmaps")
    ap_review.add_argument("--panel-fasta", required=True)
    ap_review.add_argument("--identity-matrix", required=True)
    ap_review.add_argument("--distance-matrix", required=True)
    ap_review.add_argument("--manifest-tsv", required=True)
    ap_review.add_argument("--umap-explorer-dir", default="")
    ap_review.add_argument("--thresholds", default="0.03,0.04,0.05,0.06,0.07,0.08,0.10")
    ap_review.add_argument("--novel-threshold", type=float, default=92.0)
    ap_review.add_argument("--umap-neighbors", type=int, default=30)
    ap_review.add_argument("--umap-min-dist", type=float, default=0.5)
    ap_review.add_argument("--umap-spread", type=float, default=0.8)
    ap_review.add_argument("--umap-seed", type=int, default=42)
    ap_review.add_argument("--outdir", required=True)

    ap_recomb = sub.add_parser(
        "recombination",
        help="Detect recombination breakpoints in an aligned panel via OpenRDP (RDP/GENECONV/MaxChi/Chimaera/3Seq/BootScan/SiScan)",
    )
    ap_recomb.add_argument("--alignment-fasta", required=True,
                           help="ALIGNED FASTA (e.g. a region alignment from extract_genome_regions.py + MAFFT/MUSCLE)")
    ap_recomb.add_argument("--methods", default="rdp,geneconv,maxchi,chimaera,threeseq,bootscan,siscan",
                           help="Comma-separated OpenRDP methods to run")
    ap_recomb.add_argument("--min-methods", type=int, default=3,
                           help="Minimum number of independent methods supporting an event for the consensus table (default: 3)")
    ap_recomb.add_argument("--pvalue", type=float, default=0.05, help="Maximum p-value for an event to count (default: 0.05)")
    ap_recomb.add_argument("--feature-tsv", default="",
                           help="Optional TSV (name<TAB>start<TAB>end, alignment coords) to map breakpoints onto, e.g. the IR / CP-Rep window")
    ap_recomb.add_argument("--cfg", default="", help="Optional OpenRDP config file")
    ap_recomb.add_argument("--outdir", required=True)

    ap_export = sub.add_parser("export-clusters", help="Export cluster members to per-cluster FASTA files")
    ap_export.add_argument("--sequence-fasta", required=True)
    ap_export.add_argument("--reference-fasta", default="")
    ap_export.add_argument("--cluster-tsv", required=True)
    ap_export.add_argument("--manifest-tsv", default="")
    ap_export.add_argument("--id-column", default="id")
    ap_export.add_argument("--cluster-column", default="cluster")
    ap_export.add_argument("--filter-column", default="")
    ap_export.add_argument("--filter-value", default="")
    ap_export.add_argument("--cluster-list", default="")
    ap_export.add_argument("--filename-prefix", default="")
    ap_export.add_argument("--outdir", required=True)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    if args.cmd == "umap-explorer":
        from . import umap_explorer

        umap_explorer.run(args)
    elif args.cmd == "sample-panel":
        from . import sampling

        sampling.run(args)
    elif args.cmd == "collapse-panel":
        from . import collapse_panel

        collapse_panel.run(args)
    elif args.cmd == "align-panel":
        from . import align_panel

        align_panel.run(args)
    elif args.cmd == "distance-explorer":
        from . import distance_explorer

        distance_explorer.run(args)
    elif args.cmd == "recombination":
        from . import recombination

        recombination.run(args)
    elif args.cmd == "export-clusters":
        from . import export_clusters

        export_clusters.run(args)


if __name__ == "__main__":
    main()
