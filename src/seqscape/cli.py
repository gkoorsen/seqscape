from __future__ import annotations

import argparse


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Staged SeqScape pipeline: build AF Leiden widget, sample a representative "
            "panel, compute pairwise alignments, then build the post-alignment review bundle."
        )
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_af = sub.add_parser("af-widget", help="Build the AF UMAP+Leiden parameter widget")
    ap_af.add_argument("--input-fasta", required=True)
    ap_af.add_argument("--reference-fasta", required=True)
    ap_af.add_argument("--cg-map-tsv", default="")
    ap_af.add_argument("--kmer-values", default="5")
    ap_af.add_argument("--neighbors-values", default="50,100,150,200")
    ap_af.add_argument("--min-dist-values", default="0.05,0.1,0.3,0.5,0.8")
    ap_af.add_argument("--seed", type=int, default=42)
    ap_af.add_argument("--leiden-neighbor-k", type=int, default=300)
    ap_af.add_argument("--leiden-resolution", type=float, default=1e-6)
    ap_af.add_argument("--default-kmer", type=int, default=5)
    ap_af.add_argument("--default-neighbors", type=int, default=100)
    ap_af.add_argument("--default-min-dist", type=float, default=0.3)
    ap_af.add_argument("--outdir", required=True)

    ap_sample = sub.add_parser("sample-panel", help="Sample a proportional representative panel from a chosen AF state")
    ap_sample.add_argument("--input-fasta", required=True)
    ap_sample.add_argument("--reference-fasta", required=True)
    ap_sample.add_argument("--af-widget-dir", required=True)
    ap_sample.add_argument("--chosen-kmer", type=int, required=True)
    ap_sample.add_argument("--chosen-neighbors", type=int, required=True)
    ap_sample.add_argument("--chosen-min-dist", type=float, required=True)
    ap_sample.add_argument("--sample-size", type=int, required=True)
    ap_sample.add_argument("--outdir", required=True)

    ap_align = sub.add_parser("align-panel", help="Run all-vs-all identities on the sampled panel")
    ap_align.add_argument("--panel-fasta", required=True)
    ap_align.add_argument("--manifest-tsv", required=True)
    ap_align.add_argument("--aligner", choices=["pairwisealigner", "mafft", "muscle"], default="muscle")
    ap_align.add_argument("--mafft", default="/opt/homebrew/bin/mafft")
    ap_align.add_argument("--muscle", default="muscle")
    ap_align.add_argument("--jobs", type=int, default=8)
    ap_align.add_argument("--progress-every", type=int, default=500)
    ap_align.add_argument("--outdir", required=True)

    ap_review = sub.add_parser("review-widget", help="Build the post-alignment threshold explorer, reference scatter, and heatmaps")
    ap_review.add_argument("--panel-fasta", required=True)
    ap_review.add_argument("--identity-matrix", required=True)
    ap_review.add_argument("--distance-matrix", required=True)
    ap_review.add_argument("--manifest-tsv", required=True)
    ap_review.add_argument("--thresholds", default="0.03,0.04,0.05,0.06,0.07,0.08,0.10")
    ap_review.add_argument("--novel-threshold", type=float, default=92.0)
    ap_review.add_argument("--umap-neighbors", type=int, default=30)
    ap_review.add_argument("--umap-min-dist", type=float, default=0.5)
    ap_review.add_argument("--umap-spread", type=float, default=0.8)
    ap_review.add_argument("--umap-seed", type=int, default=42)
    ap_review.add_argument("--outdir", required=True)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    if args.cmd == "af-widget":
        from . import af_widget

        af_widget.run(args)
    elif args.cmd == "sample-panel":
        from . import sampling

        sampling.run(args)
    elif args.cmd == "align-panel":
        from . import align_panel

        align_panel.run(args)
    elif args.cmd == "review-widget":
        from . import review_widget

        review_widget.run(args)


if __name__ == "__main__":
    main()
