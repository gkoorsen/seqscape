#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Run the ToCSV resistance-breaking recombination screen.

Default inputs point to the current local ToCSV full run. Override any path with
flags or environment variables.

Required tools for the full run:
  - python
  - seqscape, or run from this repo so PYTHONPATH=src can be used
  - mafft for the regional multiple-sequence alignment
  - OpenRDP for the final recombination stage

Examples:
  scripts/run_tocsv_rb_recombination_pipeline.sh --skip-recombination

  WINDOW_SIZE=400 REGION_IDENTITY_THRESHOLD=98.1 \
    scripts/run_tocsv_rb_recombination_pipeline.sh --outdir runs/rb_screen_400nt

Options:
  --clean-fasta PATH          Clean, phased, exact-deduplicated local genomes
  --comparator-fasta PATH     Collapsible public/local comparator genomes
  --reference-fasta PATH      Always-kept anchor references for UMAP context
  --support-tsv PATH          Support table for clean local genomes
  --umap-dir PATH             Reuse an existing umap-explorer directory
  --outdir PATH               Output directory
  --window-size N             Origin window size, default 400
  --region-identity PCT       Regional protection threshold, default 98.1
  --collapse-identity PCT     Whole-genome collapse threshold, default 95.0
  --methods LIST              OpenRDP methods, default geneconv,rdp,chimaera
  --min-methods N             Consensus method threshold, default 1
  --jobs N                    Threads for UMAP/MAFFT where supported, default 8
  --force-umap                Recompute umap-explorer instead of reusing it
  --skip-recombination        Stop after producing the aligned regional panel
  --dry-run                   Print commands without running them
  -h, --help                  Show this help
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

PYTHON="${PYTHON:-}"
RUN_ROOT="${RUN_ROOT:-/Users/gerritkoorsen/ciderseq-mono/cideseq-mono/runs/tocsv_260203_full_aws_20260612/full_run}"
CLEAN_FASTA="${CLEAN_FASTA:-}"
COMPARATOR_FASTA="${COMPARATOR_FASTA:-}"
REFERENCE_FASTA="${REFERENCE_FASTA:-}"
SUPPORT_TSV="${SUPPORT_TSV:-}"
UMAP_DIR="${UMAP_DIR:-}"
OUTDIR="${OUTDIR:-}"

WINDOW_SIZE="${WINDOW_SIZE:-400}"
REGION_IDENTITY_THRESHOLD="${REGION_IDENTITY_THRESHOLD:-98.1}"
COLLAPSE_IDENTITY_THRESHOLD="${COLLAPSE_IDENTITY_THRESHOLD:-95.0}"
METHODS="${METHODS:-geneconv,rdp,chimaera}"
MIN_METHODS="${MIN_METHODS:-1}"
PVALUE="${PVALUE:-0.05}"
JOBS="${JOBS:-8}"

CHOSEN_KMER="${CHOSEN_KMER:-6}"
CHOSEN_NEIGHBORS="${CHOSEN_NEIGHBORS:-50}"
CHOSEN_MIN_DIST="${CHOSEN_MIN_DIST:-0.3}"
CHOSEN_LEIDEN_RESOLUTION="${CHOSEN_LEIDEN_RESOLUTION:-0.001}"

FORCE_UMAP=0
SKIP_RECOMBINATION=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --clean-fasta) CLEAN_FASTA="$2"; shift 2 ;;
    --comparator-fasta) COMPARATOR_FASTA="$2"; shift 2 ;;
    --reference-fasta) REFERENCE_FASTA="$2"; shift 2 ;;
    --support-tsv) SUPPORT_TSV="$2"; shift 2 ;;
    --umap-dir) UMAP_DIR="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --window-size) WINDOW_SIZE="$2"; shift 2 ;;
    --region-identity) REGION_IDENTITY_THRESHOLD="$2"; shift 2 ;;
    --collapse-identity) COLLAPSE_IDENTITY_THRESHOLD="$2"; shift 2 ;;
    --methods) METHODS="$2"; shift 2 ;;
    --min-methods) MIN_METHODS="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --force-umap) FORCE_UMAP=1; shift ;;
    --skip-recombination) SKIP_RECOMBINATION=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 2 ;;
  esac
done

: "${CLEAN_FASTA:=${RUN_ROOT}/analysis/seqscape_clean_phased_20260625/clean_phased_2400_3300_unique.fasta}"
: "${COMPARATOR_FASTA:=${RUN_ROOT}/comparators/comparator_panel.fasta}"
: "${REFERENCE_FASTA:=${RUN_ROOT}/analysis/phylo_compare_20260612_223925/references.fasta}"
: "${SUPPORT_TSV:=${RUN_ROOT}/analysis/seqscape_clean_phased_20260625/clean_phased_2400_3300_unique_support.tsv}"
: "${UMAP_DIR:=${RUN_ROOT}/analysis/seqscape_clean_phased_20260625/umap_explorer}"
: "${OUTDIR:=runs/tocsv_rb_recombination_$(date +%Y%m%d_%H%M%S)}"

if [[ -z "${PYTHON}" ]]; then
  if [[ -x "/Users/gerritkoorsen/opt/anaconda3/envs/seqscape/bin/python" ]]; then
    PYTHON="/Users/gerritkoorsen/opt/anaconda3/envs/seqscape/bin/python"
  else
    PYTHON="python"
  fi
fi

export MPLCONFIGDIR="${MPLCONFIGDIR:-${OUTDIR}/.matplotlib}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${OUTDIR}/.numba_cache}"

mkdir -p "${OUTDIR}/logs" "${OUTDIR}/regions" "${MPLCONFIGDIR}" "${NUMBA_CACHE_DIR}"
COMMAND_LOG="${OUTDIR}/commands.log"
: > "${COMMAND_LOG}"

log() {
  printf '%s\n' "$*" | tee -a "${COMMAND_LOG}"
}

quote_cmd() {
  local quoted=()
  local arg
  for arg in "$@"; do
    quoted+=("$(printf '%q' "${arg}")")
  done
  printf '%s ' "${quoted[@]}"
}

run() {
  log "+ $(quote_cmd "$@")"
  if [[ "${DRY_RUN}" == "1" ]]; then
    return 0
  fi
  "$@"
}

require_file() {
  if [[ ! -f "$1" ]]; then
    echo "Missing required file: $1" >&2
    exit 1
  fi
}

require_file "${CLEAN_FASTA}"
require_file "${COMPARATOR_FASTA}"
require_file "${REFERENCE_FASTA}"
require_file "${SUPPORT_TSV}"

if [[ -n "${SEQSCAPE_CMD:-}" ]]; then
  read -r -a SEQSCAPE <<< "${SEQSCAPE_CMD}"
elif command -v seqscape >/dev/null 2>&1; then
  SEQSCAPE=(seqscape)
else
  export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"
  SEQSCAPE=("${PYTHON}" -m seqscape.cli)
fi

LOCAL_ORIGIN_FASTA="${OUTDIR}/regions/local_origin_core_${WINDOW_SIZE}nt.fasta"
LOCAL_ORIGIN_MANIFEST="${OUTDIR}/regions/local_origin_core_${WINDOW_SIZE}nt_manifest.tsv"
ANALYSIS_FASTA="${OUTDIR}/analysis_genomes_plus_comparators.fasta"
ANALYSIS_DUPLICATES="${OUTDIR}/analysis_genomes_plus_comparators_duplicates.tsv"
PROTECTED_IDS="${OUTDIR}/protected_origin_ids.txt"
PROTECTED_GROUPS="${OUTDIR}/protected_origin_groups.tsv"
COLLAPSE_DIR="${OUTDIR}/collapse_panel"
PANEL_ORIGIN_FASTA="${OUTDIR}/regions/panel_origin_core_${WINDOW_SIZE}nt.fasta"
PANEL_ORIGIN_MANIFEST="${OUTDIR}/regions/panel_origin_core_${WINDOW_SIZE}nt_manifest.tsv"
ALIGNED_FASTA="${OUTDIR}/regions/panel_origin_core_${WINDOW_SIZE}nt_aligned.fasta"
RECOMB_DIR="${OUTDIR}/recombination_origin_core_${WINDOW_SIZE}nt"

log "Output directory: ${OUTDIR}"
log "Clean genomes: ${CLEAN_FASTA}"
log "Collapsible comparators: ${COMPARATOR_FASTA}"
log "Anchor references: ${REFERENCE_FASTA}"

run "${PYTHON}" scripts/combine_fasta_unique.py \
  --fasta "${CLEAN_FASTA}" \
  --fasta "${COMPARATOR_FASTA}" \
  --exclude-fasta "${REFERENCE_FASTA}" \
  --out-fasta "${ANALYSIS_FASTA}" \
  --duplicates-tsv "${ANALYSIS_DUPLICATES}"

if [[ "${FORCE_UMAP}" != "1" && -f "${UMAP_DIR}/summary.txt" ]]; then
  UMAP_INPUT="$(awk '$1 == "input_fasta" {print $2}' "${UMAP_DIR}/summary.txt")"
  UMAP_REFERENCE="$(awk '$1 == "reference_fasta" {print $2}' "${UMAP_DIR}/summary.txt")"
  if [[ "${UMAP_INPUT}" != "${ANALYSIS_FASTA}" || "${UMAP_REFERENCE}" != "${REFERENCE_FASTA}" ]]; then
    log "Existing UMAP Explorer uses different inputs; rebuilding in ${OUTDIR}/af"
    log "  existing input: ${UMAP_INPUT}"
    log "  existing refs:  ${UMAP_REFERENCE}"
    UMAP_DIR="${OUTDIR}/af"
    FORCE_UMAP=1
  fi
fi

run "${PYTHON}" scripts/extract_origin_windows.py \
  --fasta "${CLEAN_FASTA}" \
  --window-size "${WINDOW_SIZE}" \
  --out-fasta "${LOCAL_ORIGIN_FASTA}" \
  --manifest-tsv "${LOCAL_ORIGIN_MANIFEST}"

run "${PYTHON}" scripts/select_region_protected_ids.py \
  --region-fasta "${LOCAL_ORIGIN_FASTA}" \
  --support-tsv "${SUPPORT_TSV}" \
  --identity-threshold "${REGION_IDENTITY_THRESHOLD}" \
  --out-ids "${PROTECTED_IDS}" \
  --out-groups-tsv "${PROTECTED_GROUPS}"

if [[ "${FORCE_UMAP}" == "1" || ! -d "${UMAP_DIR}/states" ]]; then
  UMAP_DIR="${OUTDIR}/af"
  run "${SEQSCAPE[@]}" umap-explorer \
    --input-fasta "${ANALYSIS_FASTA}" \
    --reference-fasta "${REFERENCE_FASTA}" \
    --kmer-values "${CHOSEN_KMER}" \
    --neighbors-values "${CHOSEN_NEIGHBORS}" \
    --min-dist-values "${CHOSEN_MIN_DIST}" \
    --leiden-resolution-values "${CHOSEN_LEIDEN_RESOLUTION}" \
    --default-kmer "${CHOSEN_KMER}" \
    --default-neighbors "${CHOSEN_NEIGHBORS}" \
    --default-min-dist "${CHOSEN_MIN_DIST}" \
    --default-leiden-resolution "${CHOSEN_LEIDEN_RESOLUTION}" \
    --outdir "${UMAP_DIR}"
else
  log "Reusing existing UMAP Explorer: ${UMAP_DIR}"
fi

run "${SEQSCAPE[@]}" collapse-panel \
  --input-fasta "${ANALYSIS_FASTA}" \
  --reference-fasta "${REFERENCE_FASTA}" \
  --umap-explorer-dir "${UMAP_DIR}" \
  --chosen-kmer "${CHOSEN_KMER}" \
  --chosen-neighbors "${CHOSEN_NEIGHBORS}" \
  --chosen-min-dist "${CHOSEN_MIN_DIST}" \
  --chosen-leiden-resolution "${CHOSEN_LEIDEN_RESOLUTION}" \
  --identity-threshold "${COLLAPSE_IDENTITY_THRESHOLD}" \
  --focal-ids "${PROTECTED_IDS}" \
  --outdir "${COLLAPSE_DIR}"

run "${PYTHON}" scripts/extract_origin_windows.py \
  --fasta "${COLLAPSE_DIR}/representative_panel.fasta" \
  --window-size "${WINDOW_SIZE}" \
  --out-fasta "${PANEL_ORIGIN_FASTA}" \
  --manifest-tsv "${PANEL_ORIGIN_MANIFEST}"

if [[ "${DRY_RUN}" == "1" ]]; then
  log "+ mafft --auto --thread ${JOBS} ${PANEL_ORIGIN_FASTA} > ${ALIGNED_FASTA}"
else
  if ! command -v mafft >/dev/null 2>&1; then
    echo "mafft is required to produce the OpenRDP multiple sequence alignment" >&2
    exit 1
  fi
  log "+ mafft --auto --thread ${JOBS} ${PANEL_ORIGIN_FASTA} > ${ALIGNED_FASTA}"
  mafft --auto --thread "${JOBS}" "${PANEL_ORIGIN_FASTA}" \
    > "${ALIGNED_FASTA}" \
    2> "${OUTDIR}/logs/mafft_origin_core_${WINDOW_SIZE}nt.log"
fi

if [[ "${SKIP_RECOMBINATION}" == "1" ]]; then
  log "Skipping recombination stage. Aligned FASTA: ${ALIGNED_FASTA}"
else
  run "${SEQSCAPE[@]}" recombination \
    --alignment-fasta "${ALIGNED_FASTA}" \
    --methods "${METHODS}" \
    --min-methods "${MIN_METHODS}" \
    --pvalue "${PVALUE}" \
    --outdir "${RECOMB_DIR}"
fi

cat > "${OUTDIR}/README.txt" <<EOF
ToCSV RB recombination screen

clean_fasta	${CLEAN_FASTA}
analysis_fasta	${ANALYSIS_FASTA}
comparator_fasta	${COMPARATOR_FASTA}
reference_fasta	${REFERENCE_FASTA}
support_tsv	${SUPPORT_TSV}
umap_dir	${UMAP_DIR}
window_size	${WINDOW_SIZE}
region_identity_threshold	${REGION_IDENTITY_THRESHOLD}
collapse_identity_threshold	${COLLAPSE_IDENTITY_THRESHOLD}
methods	${METHODS}
min_methods	${MIN_METHODS}

Key outputs:
protected_ids	${PROTECTED_IDS}
collapse_manifest	${COLLAPSE_DIR}/representative_panel_manifest.tsv
regional_panel	${PANEL_ORIGIN_FASTA}
regional_alignment	${ALIGNED_FASTA}
recombination_dir	${RECOMB_DIR}
commands	${COMMAND_LOG}
EOF

log "Pipeline complete: ${OUTDIR}"
