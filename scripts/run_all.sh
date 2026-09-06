#!/usr/bin/env bash
# =============================================================================
#  run_all.sh - end-to-end Ribo-seq analysis for any dataset
# =============================================================================
#  Runs every stage in order, skipping work already done, and finishes with an
#  HTML report. Each stage is resumable, so a failure costs only that stage.
#
#  Usage:
#      bash run_all.sh                       # everything, interactive config
#      bash run_all.sh --threads 12 --ram 32 # non-interactive
#      bash run_all.sh --from te             # skip to downstream analysis
#      bash run_all.sh --only report
#      bash run_all.sh --skip-xtail          # xtail is the slow stage
#      bash run_all.sh --bins 2000           # faster, coarser xtail p-values
#
#  Stages: configure install genome fastq pipeline te expression de report
# =============================================================================
set -euo pipefail
GUIDE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
SCRIPTS="${GUIDE_DIR}/scripts"

ALL_STAGES=(configure install genome fastq pipeline te expression de report)
FROM=""; ONLY=""; SKIP_XTAIL=0; BINS=""; PASSTHRU=()
ALPHA="0.05"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)        FROM="$2"; shift 2 ;;
        --only)        ONLY="$2"; shift 2 ;;
        --skip-xtail)  SKIP_XTAIL=1; shift ;;
        --bins)        BINS="$2"; shift 2 ;;
        --alpha)       ALPHA="$2"; shift 2 ;;
        --threads|--ram|--project|--accession)
                       PASSTHRU+=("$1" "$2"); shift 2 ;;
        --yes|-y)      PASSTHRU+=("--yes"); shift ;;
        -h|--help)     sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

want() {  # should this stage run?
    local stage="$1"
    [[ -n "$ONLY" ]] && { [[ "$stage" == "$ONLY" ]]; return; }
    if [[ -n "$FROM" ]]; then
        local seen=0
        for s in "${ALL_STAGES[@]}"; do
            [[ "$s" == "$FROM" ]] && seen=1
            [[ "$s" == "$stage" ]] && { (( seen )); return; }
        done
        return 1
    fi
    return 0
}

banner() { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }
START=$(date +%s)

# ---- 1. configure -----------------------------------------------------------
if want configure; then
    banner "configure"
    bash "${SCRIPTS}/00_configure.sh" "${PASSTHRU[@]}"
fi
[[ -f "${GUIDE_DIR}/project.conf" ]] || { echo "[error] no project.conf - run the configure stage" >&2; exit 1; }
# shellcheck disable=SC1090
source "${GUIDE_DIR}/project.conf"
mkdir -p "$LOG_DIR"

run_r() {  # run_r <script> [args...] - inside te_analysis, from PROJECT_DIR
    local script="$1"; shift
    ( cd "$PROJECT_DIR" && Rscript "${SCRIPTS}/${script}" "$@" )
}

# ---- 2-5. install, genome, fastq, pipeline ----------------------------------
want install  && { banner "install";        bash "${SCRIPTS}/01_install.sh"; }
want genome   && { banner "genome";         bash "${SCRIPTS}/02_prepare_genome.sh"; }
want fastq    && { banner "fastq download"; bash "${SCRIPTS}/03_download_fastq.sh"; }
want pipeline && { banner "nextflow pipeline"; bash "${SCRIPTS}/04_run_pipeline.sh"; }

# ---- 6. downstream analysis (needs the te_analysis env) ---------------------
if want te || want expression || want de || want report; then
    CONDA_ROOT="${CONDA_ROOT:-$HOME/miniconda3}"
    # shellcheck disable=SC1091
    source "${CONDA_ROOT}/etc/profile.d/conda.sh"
    conda activate "$TE_ENV_NAME"
    export TE_THREADS
    # The R scripts source te_common.R and resolve inputs relative to the
    # working directory, so they are run from PROJECT_DIR with SCRIPTS on hand.
    cp -f "${SCRIPTS}/te_common.R" "${PROJECT_DIR}/te_common.R"
fi

COUNTS="${RESULTS_DIR}/quantification/salmon/salmon.merged.gene_counts.tsv"
DE_OUT="${RESULTS_DIR}/differential_expression"

if want te; then
    banner "translational efficiency"
    if [[ ! -f "$COUNTS" ]]; then
        echo "[warn] ${COUNTS} not found - the pipeline stage has not produced counts yet."
        echo "[warn] skipping TE, expression and DE stages."
    else
        run_r 05_te_deseq2.R  --counts "$COUNTS" --outdir "$DE_OUT" \
            2>&1 | tee "${LOG_DIR}/te_deseq2.log"  | grep -E '^\[' || true

        if (( SKIP_XTAIL )); then
            echo "[skip] xtail (--skip-xtail)"
        else
            XA=(--counts "$COUNTS" --outdir "$DE_OUT" --threads "$TE_THREADS")
            [[ -n "$BINS" ]] && XA+=(--bins "$BINS")
            echo "[te] xtail is the slow stage (~20 min at bins=10000);"
            echo "     --bins 2000 is ~5x faster with coarser p-values."
            run_r 06_te_xtail.R "${XA[@]}" \
                2>&1 | tee "${LOG_DIR}/te_xtail.log" | grep -E '^\[' || true
        fi

        run_r 07_te_riborex.R --counts "$COUNTS" --outdir "$DE_OUT" \
            2>&1 | tee "${LOG_DIR}/te_riborex.log" | grep -E '^\[' || true
    fi
fi

if want expression && [[ -f "$COUNTS" ]]; then
    banner "expression matrices"
    run_r 08_expression_matrices.R --counts "$COUNTS" --outdir "$DE_OUT" \
        2>&1 | tee "${LOG_DIR}/expression.log" | grep -E '^\[' || true
fi

if want de && [[ -f "$COUNTS" ]]; then
    banner "differential expression + master table"
    run_r 09_de_rna_ribo.R --counts "$COUNTS" --outdir "$DE_OUT" --alpha "$ALPHA" \
        2>&1 | tee "${LOG_DIR}/de.log" | grep -E '^\[|^    ' || true
fi

if want report; then
    banner "report"
    PROJECT_CONF="${GUIDE_DIR}/project.conf" \
        run_r 10_report.R --outdir "$DE_OUT" --alpha "$ALPHA" \
        2>&1 | tee "${LOG_DIR}/report.log" | grep -E '^\[' || true
fi

ELAPSED=$(( $(date +%s) - START ))
printf '\n\033[1m=== finished in %dh %dm %ds ===\033[0m\n' \
       $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60))
echo "  results: ${RESULTS_DIR}"
[[ -f "${DE_OUT}/riboseq_report.html" ]] && echo "  report : ${DE_OUT}/riboseq_report.html"
