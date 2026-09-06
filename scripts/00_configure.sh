#!/usr/bin/env bash
# =============================================================================
#  00_configure.sh - detect the machine, ask the user, write project.conf
# =============================================================================
#  Everything downstream reads project.conf. Nothing else hardcodes a thread
#  count, a memory ceiling or a path, so the same pipeline runs on a 12-thread
#  workstation and a 128-thread server without editing a single script.
#
#  Usage:
#      bash 00_configure.sh                      # interactive
#      bash 00_configure.sh --threads 12 --ram 32 --yes
#      bash 00_configure.sh --yes                # accept all detected defaults
# =============================================================================
set -euo pipefail

GUIDE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
CONF="${GUIDE_DIR}/project.conf"

THREADS=""; RAM_GB=""; ASSUME_YES=0
PROJECT_DIR=""; GENOME_ACC=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)   THREADS="$2";     shift 2 ;;
        --ram)       RAM_GB="$2";      shift 2 ;;
        --project)   PROJECT_DIR="$2"; shift 2 ;;
        --accession) GENOME_ACC="$2";  shift 2 ;;
        --yes|-y)    ASSUME_YES=1;     shift ;;
        -h|--help)   sed -n '2,14p' "$0"; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

# -----------------------------------------------------------------------------
# 1. Detect what the machine actually has
# -----------------------------------------------------------------------------
DET_THREADS="$(nproc 2>/dev/null || echo 4)"
# MemTotal is what the kernel sees; `free -g` rounds down and under-reports.
DET_RAM_GB="$(awk '/MemTotal/ {printf "%d", $2/1024/1024}' /proc/meminfo 2>/dev/null || echo 8)"
DET_DISK_GB="$(df -BG --output=avail "${GUIDE_DIR}" 2>/dev/null | tail -1 | tr -dc '0-9' || echo 0)"

echo "=============================================================="
echo " riboseq guide - machine configuration"
echo "=============================================================="
echo "  detected CPU threads : ${DET_THREADS}"
echo "  detected RAM         : ${DET_RAM_GB} GB"
echo "  free disk here       : ${DET_DISK_GB} GB"
echo

ask() {  # ask <prompt> <default> <varname>
    local prompt="$1" default="$2" __var="$3" reply
    if [[ $ASSUME_YES -eq 1 || ! -t 0 ]]; then
        printf -v "$__var" '%s' "$default"; return
    fi
    read -r -p "  ${prompt} [${default}]: " reply || reply=""
    printf -v "$__var" '%s' "${reply:-$default}"
}

[[ -n "$THREADS"     ]] || ask "CPU threads to use"         "$DET_THREADS"    THREADS
[[ -n "$RAM_GB"      ]] || ask "RAM to use (GB)"            "$DET_RAM_GB"     RAM_GB
[[ -n "$PROJECT_DIR" ]] || ask "Project/working directory"  "$(dirname "$GUIDE_DIR")/riboseq_run" PROJECT_DIR
[[ -n "$GENOME_ACC"  ]] || ask "NCBI genome accession"      "GCF_000007065.1" GENOME_ACC

# -----------------------------------------------------------------------------
# 2. Validate - a wrong number here wastes hours later
# -----------------------------------------------------------------------------
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "[error] threads must be an integer" >&2; exit 1; }
[[ "$RAM_GB"  =~ ^[0-9]+$ ]] || { echo "[error] RAM must be an integer (GB)" >&2; exit 1; }
(( THREADS >= 1 )) || { echo "[error] need at least 1 thread" >&2; exit 1; }
(( RAM_GB  >= 4 )) || { echo "[error] need at least 4 GB RAM" >&2; exit 1; }

if (( THREADS > DET_THREADS )); then
    echo "[warn] asking for ${THREADS} threads but only ${DET_THREADS} exist."
    echo "       Oversubscribing makes jobs slower, not faster."
fi
if (( RAM_GB > DET_RAM_GB )); then
    echo "[warn] asking for ${RAM_GB} GB but the machine has ${DET_RAM_GB} GB."
    echo "       Tasks will be OOM-killed. Lower it unless you have swap to burn."
fi
if (( DET_DISK_GB < 50 )); then
    echo "[warn] only ${DET_DISK_GB} GB free. A small riboseq run needs ~50 GB"
    echo "       for FASTQs, the work directory and results."
fi

# -----------------------------------------------------------------------------
# 3. Derive the resource ladder from the two numbers
# -----------------------------------------------------------------------------
#  Reserve headroom: the OS, the JVM driver and Singularity all need memory that
#  is not available to tasks. A task requesting every last GB gets OOM-killed.
#  ~12%, floor 2 GB, cap 8 GB.
RESERVE=$(( RAM_GB * 12 / 100 )); (( RESERVE < 2 )) && RESERVE=2; (( RESERVE > 8 )) && RESERVE=8
MAX_MEM=$(( RAM_GB - RESERVE ))

#  Per-label requests, deliberately BELOW the ceiling: if process_high took
#  every core, tasks would run one at a time and the machine would idle
#  whenever that task was single-threaded.
HIGH_CPUS=$(( THREADS * 2 / 3 ));   (( HIGH_CPUS < 1 )) && HIGH_CPUS=1
MED_CPUS=$((  THREADS / 3 ));       (( MED_CPUS  < 1 )) && MED_CPUS=1
LOW_CPUS=$((  THREADS / 6 ));       (( LOW_CPUS  < 1 )) && LOW_CPUS=1

HIGH_MEM=$(( MAX_MEM * 2 / 3 ));    (( HIGH_MEM < 4 )) && HIGH_MEM=4
MED_MEM=$((  MAX_MEM / 3 ));        (( MED_MEM  < 2 )) && MED_MEM=2
LOW_MEM=$((  MAX_MEM / 6 ));        (( LOW_MEM  < 2 )) && LOW_MEM=2
HIGHMEM=$((  MAX_MEM * 85 / 100 )); (( HIGHMEM < 4 )) && HIGHMEM=4

#  The floors above can exceed the ceiling on a very small machine (a 4 GB box
#  leaves 2 GB usable, but the process_high floor is 4 GB). resourceLimits would
#  clamp it at run time, but the generated config would advertise a request the
#  machine cannot honour - so clamp here and keep config and reality identical.
for v in HIGH_MEM MED_MEM LOW_MEM HIGHMEM; do
    if (( ${!v} > MAX_MEM )); then
        printf -v "$v" '%s' "$MAX_MEM"
    fi
done
(( HIGH_CPUS > THREADS )) && HIGH_CPUS=$THREADS
(( MED_CPUS  > THREADS )) && MED_CPUS=$THREADS
(( LOW_CPUS  > THREADS )) && LOW_CPUS=$THREADS

#  JVM driver heap: the Nextflow process itself, not the tasks. It only tracks
#  state, so it never needs much - but starving it stalls very large runs.
JVM_HEAP=$(( RAM_GB / 8 )); (( JVM_HEAP < 2 )) && JVM_HEAP=2; (( JVM_HEAP > 8 )) && JVM_HEAP=8

#  Parallel FASTQ downloads: network-bound, not CPU-bound. Beyond ~4 streams
#  this rarely helps and ENA starts throttling.
DL_JOBS=$(( THREADS / 3 )); (( DL_JOBS < 1 )) && DL_JOBS=1; (( DL_JOBS > 4 )) && DL_JOBS=4

mkdir -p "$PROJECT_DIR"
PROJECT_DIR="$(cd "$PROJECT_DIR" && pwd)"

cat > "$CONF" <<EOF
# =============================================================================
#  project.conf - generated by 00_configure.sh on $(date -Iseconds)
#  Re-run 00_configure.sh to change these.
# =============================================================================

# ---- hardware ----
THREADS=${THREADS}
RAM_GB=${RAM_GB}
MAX_MEM_GB=${MAX_MEM}          # ceiling handed to Nextflow (RAM minus headroom)
RESERVE_GB=${RESERVE}          # left for OS + JVM + Singularity
JVM_HEAP_GB=${JVM_HEAP}        # Nextflow driver heap

# ---- derived per-label resources ----
HIGH_CPUS=${HIGH_CPUS}
MED_CPUS=${MED_CPUS}
LOW_CPUS=${LOW_CPUS}
HIGH_MEM_GB=${HIGH_MEM}
MED_MEM_GB=${MED_MEM}
LOW_MEM_GB=${LOW_MEM}
HIGHMEM_GB=${HIGHMEM}

# ---- job concurrency ----
DL_JOBS=${DL_JOBS}             # parallel FASTQ downloads (network-bound)
TE_THREADS=${HIGH_CPUS}        # threads for Xtail

# ---- paths ----
PROJECT_DIR=${PROJECT_DIR}
GUIDE_DIR=${GUIDE_DIR}
WORK_DIR=${PROJECT_DIR}/work
CACHE_DIR=${PROJECT_DIR}/cache
RAW_DIR=${PROJECT_DIR}/raw_data
RESULTS_DIR=${PROJECT_DIR}/results
LOG_DIR=${PROJECT_DIR}/logs

# ---- reference genome ----
GENOME_ACC=${GENOME_ACC}
GTF_PAD=50                     # synthetic UTR padding, nt

# ---- pipeline ----
PIPELINE=nf-core/riboseq
PIPELINE_REV=1.2.0
CONDA_ENV_NAME=nextflow
TE_ENV_NAME=te_analysis
EOF

echo
echo "  Wrote ${CONF}"
echo "  ------------------------------------------------------------"
printf "  %-22s %s\n" "usable for tasks:"  "${MAX_MEM} GB of ${RAM_GB} GB (${RESERVE} GB reserved)"
printf "  %-22s %s\n" "process_high:"      "${HIGH_CPUS} cpus / ${HIGH_MEM} GB"
printf "  %-22s %s\n" "process_medium:"    "${MED_CPUS} cpus / ${MED_MEM} GB"
printf "  %-22s %s\n" "process_low:"       "${LOW_CPUS} cpus / ${LOW_MEM} GB"
printf "  %-22s %s\n" "high_memory:"       "${HIGHMEM} GB"
printf "  %-22s %s\n" "JVM driver heap:"   "${JVM_HEAP} GB"
printf "  %-22s %s\n" "download jobs:"     "${DL_JOBS}"
printf "  %-22s %s\n" "project dir:"       "${PROJECT_DIR}"
echo "  ------------------------------------------------------------"
echo "  Next:  bash scripts/01_install.sh"
