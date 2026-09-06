#!/usr/bin/env bash
# =============================================================================
#  lib_common.sh - shared helpers sourced by every 0N_*.sh stage script
# =============================================================================
#  Every 0N_*.sh stage script sources this file the same way:
#    STEP=<stage>; source ".../lib_common.sh"
#    mkdirs                      # create project dirs from project.conf
#    export_caches               # (04_run_pipeline.sh only) NXF_* env vars
#    activate_env "$ENV_NAME"    # (04_run_pipeline.sh only) conda activate
#    log "message"   / die "message"
#  and reads/writes the same project.conf variables documented alongside
#  00_configure.sh. If real behaviour ever needs to diverge from what's
#  implemented here, that divergence is confined to this one file.
# =============================================================================
set -euo pipefail

_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
GUIDE_DIR="$(cd "${_LIB_DIR}/.." && pwd)"
CONF="${GUIDE_DIR}/project.conf"

[[ -f "$CONF" ]] || {
    echo "[error] ${CONF} not found - run scripts/00_configure.sh first" >&2
    exit 1
}
# shellcheck disable=SC1090
source "$CONF"

STEP="${STEP:-guide}"

log() { printf '[%s] %s\n' "$STEP" "$*"; }
die() { printf '[%s] [error] %s\n' "$STEP" "$*" >&2; exit 1; }

# Create every directory project.conf points at. Safe to call repeatedly.
mkdirs() {
    mkdir -p \
        "$PROJECT_DIR" "$WORK_DIR" "$CACHE_DIR" "$CACHE_DIR/singularity" \
        "$RAW_DIR" "$RESULTS_DIR" "$LOG_DIR"
}

# Environment variables Nextflow itself should see. Keeps every cache under
# the project directory instead of $HOME, and selects the legacy v1 config
# parser: nf-core/riboseq 1.2.0's modules.config mixes `if` statements with
# plain config statements, which Nextflow >=26's strict v2 parser rejects
# ("If statements cannot be mixed with config statements"). v1 is already the
# default on Nextflow 25.x, so exporting it is a no-op there and required here.
export_caches() {
    export NXF_SYNTAX_PARSER=v1
    export NXF_HOME="${CACHE_DIR}/nextflow"
    export NXF_SINGULARITY_CACHEDIR="${CACHE_DIR}/singularity"
    export NXF_ASSETS="${CACHE_DIR}/nextflow-assets"
    mkdir -p "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR" "$NXF_ASSETS"
}

# activate_env <name> - activate a conda env by name, from CONDA_ROOT.
activate_env() {
    local env_name="$1"
    local root="${CONDA_ROOT:-$HOME/miniconda3}"
    [[ -f "${root}/etc/profile.d/conda.sh" ]] || die "miniconda not found at ${root} - set CONDA_ROOT"
    # shellcheck disable=SC1091
    source "${root}/etc/profile.d/conda.sh"
    conda env list | awk '{print $1}' | grep -qx "$env_name" \
        || die "conda env '${env_name}' does not exist - run scripts/01_install.sh"
    conda activate "$env_name"
}
