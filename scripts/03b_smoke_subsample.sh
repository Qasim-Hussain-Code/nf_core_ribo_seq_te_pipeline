#!/usr/bin/env bash
# =============================================================================
#  03b_smoke_subsample.sh - fast, low-RAM smoke test before the full run
# =============================================================================
#  Added for this analysis because the target machine has only ~7.6 GB RAM.
#  Instead of downloading the
#  full ~4.6 GB of FASTQs (as 03_download_fastq.sh does) and only then finding
#  out whether the genome/config/pipeline plumbing works, this script streams
#  just the first N reads of each ENA run directly (curl | zcat | head | gzip),
#  so nothing beyond those reads is ever transferred or written to disk.
#
#  It writes raw_data_smoke/<sample>.fastq.gz and samplesheet.smoke.csv
#  (a copy of samplesheet.csv with fastq_1 pointed at the subsampled files).
#  Run the pipeline against samplesheet.smoke.csv first; once it completes
#  end-to-end, run 03_download_fastq.sh + 04_run_pipeline.sh for the real
#  full-data analysis against the original samplesheet.csv.
#
#  Usage:
#      bash 03b_smoke_subsample.sh                # 1,000,000 reads/sample
#      bash 03b_smoke_subsample.sh --reads 200000  # smaller/faster
# =============================================================================
STEP=smoke
source "$(dirname "${BASH_SOURCE[0]:-$0}")/lib_common.sh"
mkdirs
cd "$PROJECT_DIR"

READS=1000000
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads) READS="$2"; shift 2 ;;
        *) die "unknown option: $1" ;;
    esac
done
LINES=$(( READS * 4 ))

SAMPLESHEET="${PROJECT_DIR}/samplesheet.csv"
SMOKE_SHEET="${PROJECT_DIR}/samplesheet.smoke.csv"
SMOKE_DIR="${PROJECT_DIR}/raw_data_smoke"
mkdir -p "$SMOKE_DIR"

[[ -f "$SAMPLESHEET" ]] || die "samplesheet not found: ${SAMPLESHEET}"

COL_SAMPLE=$(head -1 "$SAMPLESHEET" | tr ',' '\n' | grep -n '^sample$' | cut -d: -f1)
COL_SRR=$(head -1 "$SAMPLESHEET" | tr ',' '\n' | grep -n '^srr$' | cut -d: -f1)
[[ -n "$COL_SAMPLE" && -n "$COL_SRR" ]] || die "samplesheet needs 'sample' and 'srr' columns"

ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"

log "subsampling ${READS} reads/sample from ENA (streamed, nothing full is downloaded) ..."

awk -F, -v s="$COL_SAMPLE" -v r="$COL_SRR" 'NR>1 && $r != "" {print $s","$r}' "$SAMPLESHEET" | \
while IFS=, read -r sample srr; do
    dest="${SMOKE_DIR}/${sample}.fastq.gz"
    if [[ -s "$dest" ]]; then
        log "${sample}: already subsampled - skipping"
        continue
    fi
    ftp=$(curl -sf --retry 3 --retry-delay 5 \
        "${ENA_API}?accession=${srr}&result=read_run&fields=fastq_ftp&format=tsv" \
        | awk 'NR==2')
    [[ -n "$ftp" ]] || die "ENA returned no FTP path for ${srr} (${sample})"
    url="ftp://${ftp}"
    log "${sample} <- ${srr} (first ${READS} reads) ..."
    # SIGPIPE from `head` closing the pipe early stops curl/zcat before the
    # remote file is fully transferred - only the requested reads are fetched.
    set +o pipefail
    curl -sSL "$url" | zcat | head -n "$LINES" | gzip -1 > "${dest}.tmp"
    set -o pipefail
    mv "${dest}.tmp" "$dest"
    n=$(( $(zcat "$dest" | wc -l) / 4 ))
    log "${sample}: wrote ${n} reads -> ${dest}"
done

python3 - "$SAMPLESHEET" "$SMOKE_SHEET" "$SMOKE_DIR" <<'PY'
import csv, sys, os
sheet, out, raw = sys.argv[1], sys.argv[2], sys.argv[3]
with open(sheet, newline='') as fh:
    rows = list(csv.DictReader(fh)); fields = list(rows[0].keys())
for r in rows:
    r['fastq_1'] = os.path.join(raw, f"{r['sample']}.fastq.gz")
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
print(f"[smoke] wrote {out}")
PY

log "done. Next: bash scripts/04_run_pipeline.sh --input ${SMOKE_SHEET}  (or see run_all.sh notes)"
