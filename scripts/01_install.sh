#!/usr/bin/env bash
# =============================================================================
#  01_install.sh - Nextflow, Singularity and the two conda environments
# =============================================================================
#  Two separate environments on purpose:
#    nextflow      drives the pipeline
#    te_analysis   R + DESeq2 + Xtail + riborex for downstream analysis
#  Xtail and riborex are unmaintained and compiled from GitHub; if a build
#  breaks against a newer R, only te_analysis is affected, not the pipeline.
# =============================================================================
STEP=install
source "$(dirname "${BASH_SOURCE[0]:-$0}")/lib_common.sh"

CONDA_ROOT="${CONDA_ROOT:-$HOME/miniconda3}"
[[ -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ]] || die "install miniconda first, or set CONDA_ROOT"
# shellcheck disable=SC1091
source "${CONDA_ROOT}/etc/profile.d/conda.sh"

have_env() { conda env list | awk '{print $1}' | grep -qx "$1"; }

# -----------------------------------------------------------------------------
# 1. Singularity - must come from the system, not conda
# -----------------------------------------------------------------------------
if command -v singularity >/dev/null || command -v apptainer >/dev/null; then
    log "container runtime: $(command -v singularity || command -v apptainer)"
else
    log "Singularity/Apptainer not found. Install it with ONE of:"
    log "    sudo apt-get update && sudo apt-get install -y singularity-container"
    log "    sudo apt-get install -y apptainer"
    die "container runtime required - the pipeline runs every tool in a container"
fi

# -----------------------------------------------------------------------------
# 2. Pipeline environment
# -----------------------------------------------------------------------------
if have_env "$CONDA_ENV_NAME"; then
    log "conda env '${CONDA_ENV_NAME}' exists - skipping"
else
    log "creating conda env '${CONDA_ENV_NAME}' ..."
    mamba create -n "$CONDA_ENV_NAME" -y -c conda-forge -c bioconda nextflow openjdk
fi

conda activate "$CONDA_ENV_NAME"
log "nextflow: $(nextflow -v 2>&1)"

log "pulling ${PIPELINE} revision ${PIPELINE_REV} ..."
nextflow pull "$PIPELINE" -r "$PIPELINE_REV" || log "pull failed - will retry at run time"

# -----------------------------------------------------------------------------
# 3. Downstream analysis environment
# -----------------------------------------------------------------------------
if have_env "$TE_ENV_NAME"; then
    log "conda env '${TE_ENV_NAME}' exists - skipping"
else
    log "creating conda env '${TE_ENV_NAME}' (several minutes) ..."
    # R pinned to 4.3: Xtail's Rcpp code does not compile against R >= 4.4.
    mamba create -n "$TE_ENV_NAME" -y -c conda-forge -c bioconda \
        r-base=4.3 bioconductor-deseq2 bioconductor-edger bioconductor-biocparallel \
        r-remotes r-rcpp r-rcpparmadillo r-ggplot2 r-data.table r-fdrtool \
        r-matrixstats r-locfit r-knitr r-rmarkdown compilers make
fi

conda activate "$TE_ENV_NAME"

# --- Xtail: DESCRIPTION must be repaired before R will install it ------------
# The authors obfuscated their address against spam harvesters:
#     person("Zhengtao","Xiao", email="zhengtao.xiao[at]xjtu.edu.cn", ...)
# R derives Maintainer from Authors@R and demands a real email, so "[at]" makes
# it reject the package with "Malformed maintainer field" - at BOTH `R CMD
# build` and `R CMD INSTALL`, which is why remotes' build=FALSE is not enough.
if Rscript -e 'quit(status = !requireNamespace("xtail", quietly = TRUE))'; then
    log "xtail already installed"
else
    log "installing xtail (patching DESCRIPTION) ..."
    XT="$(mktemp -d)"; trap 'rm -rf "$XT"' EXIT
    curl -sSL https://github.com/xryanglab/xtail/archive/refs/heads/master.tar.gz -o "${XT}/x.tgz"
    tar xzf "${XT}/x.tgz" -C "$XT"
    sed -i 's/\[at\]/@/g' "${XT}/xtail-master/DESCRIPTION"
    R CMD INSTALL "${XT}/xtail-master"
fi

Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (requireNamespace("riborex", quietly = TRUE)) {
    cat("[install] riborex already installed\n")
} else {
    cat("[install] installing riborex ...\n")
    remotes::install_github("smithlabcode/riborex", upgrade = "never", quiet = FALSE)
}'

# -----------------------------------------------------------------------------
# 4. Verify - a package that installs but will not load is still broken
# -----------------------------------------------------------------------------
log "verifying R packages ..."
Rscript -e '
pkgs <- c("DESeq2","xtail","riborex","ggplot2","data.table","rmarkdown")
bad <- character()
for (p in pkgs) {
  ok <- suppressWarnings(suppressPackageStartupMessages(
        require(p, character.only = TRUE, quietly = TRUE)))
  cat(sprintf("  %-12s %s\n", p, if (ok) paste("OK", as.character(packageVersion(p))) else "FAILED"))
  if (!ok) bad <- c(bad, p)
}
if (length(bad)) { cat("\n[error] failed to load:", paste(bad, collapse=", "), "\n"); quit(status=1) }
cat("\n[install] all R packages available\n")'

log "done. Next:  bash scripts/02_prepare_genome.sh"
