#!/usr/bin/env bash
# =============================================================================
#  03_download_fastq.sh - fetch raw reads from the ENA FTP mirror
# =============================================================================
#  ENA hosts the same runs as SRA but serves ready-made .fastq.gz over FTP, so
#  there is no fasterq-dump conversion - far faster than the SRA toolkit.
#
#  Accessions come from the samplesheet's `srr` column, so the two files cannot
#  drift apart. Safe to re-run: verified files are skipped, partial files
#  resume, corrupt files are deleted so the next run starts clean.
# =============================================================================
STEP=fastq
source "$(dirname "${BASH_SOURCE[0]:-$0}")/lib_common.sh"
mkdirs
cd "$PROJECT_DIR"

SAMPLESHEET="${1:-${PROJECT_DIR}/samplesheet.csv}"
MANIFEST="${RAW_DIR}/ena_manifest.tsv"
ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"

[[ -f "$SAMPLESHEET" ]] || die "samplesheet not found: ${SAMPLESHEET}
    It needs columns: sample,fastq_1,fastq_2,strandedness,type,condition,replicate,pair,srr"

# Locate the sample and srr columns by NAME, not position - samplesheets differ.
COL_SAMPLE=$(head -1 "$SAMPLESHEET" | tr ',' '\n' | grep -n '^sample$' | cut -d: -f1)
COL_SRR=$(head -1 "$SAMPLESHEET"    | tr ',' '\n' | grep -n '^srr$'    | cut -d: -f1)
[[ -n "$COL_SAMPLE" ]] || die "samplesheet has no 'sample' column"
[[ -n "$COL_SRR" ]] || die "samplesheet has no 'srr' column - add run accessions, or place FASTQs manually"

# -----------------------------------------------------------------------------
# 1. Build sample -> FTP url + md5 manifest from ENA
# -----------------------------------------------------------------------------
if [[ ! -s "$MANIFEST" ]]; then
    log "querying ENA for FTP paths ..."
    : > "${MANIFEST}.tmp"
    awk -F, -v s="$COL_SAMPLE" -v r="$COL_SRR" 'NR>1 && $r != "" {print $s"\t"$r}' "$SAMPLESHEET" |
    while IFS=$'\t' read -r sample srr; do
        row=$(curl -sf --retry 3 --retry-delay 5 \
              "${ENA_API}?accession=${srr}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" \
              | awk 'NR==2')
        ftp=$(echo "$row" | cut -f2); md5=$(echo "$row" | cut -f3)
        [[ -n "$ftp" ]] || { echo "[error] ENA returned no FTP path for ${srr} (${sample})" >&2; exit 1; }
        [[ "$ftp" == *";"* ]] && { echo "[error] ${srr} is paired-end; this script handles single-end" >&2; exit 1; }
        printf '%s\t%s\tftp://%s\t%s\n' "$sample" "$srr" "$ftp" "$md5" >> "${MANIFEST}.tmp"
        echo "[fastq]   ${sample}  <-  ${srr}"
    done
    mv "${MANIFEST}.tmp" "$MANIFEST"
fi
log "$(wc -l < "$MANIFEST") runs -> ${RAW_DIR}  (${DL_JOBS} parallel streams)"

# -----------------------------------------------------------------------------
# 2. Download, resume, verify
# -----------------------------------------------------------------------------
fetch_one() {
    local sample srr url md5 dest
    IFS=$'\t' read -r sample srr url md5 <<< "$1"
    dest="${RAW_DIR}/${sample}.fastq.gz"

    if [[ -s "$dest" ]] && [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$md5" ]]; then
        echo "[fastq] ${sample}: verified - skipping"; return 0
    fi
    echo "[fastq] ${sample}: downloading ${srr} ..."
    # -nv not --show-progress: several parallel wgets interleave their progress
    # bars into unreadable noise and bury real errors.
    wget -c -nv --tries=10 --waitretry=15 -O "$dest" "$url"

    if [[ "$(md5sum "$dest" | cut -d' ' -f1)" != "$md5" ]]; then
        # Discard it. A half-written file makes `wget -c` resume onto corrupt
        # bytes forever, and the pipeline only finds out much later with a
        # cryptic "corrupt deflate stream" from fq lint.
        rm -f "$dest"
        echo "[error] ${sample}: md5 mismatch - removed, re-run to fetch again" >&2
        return 1
    fi
    echo "[fastq] ${sample}: OK"
}
export -f fetch_one
export RAW_DIR

status=0
xargs -a "$MANIFEST" -d '\n' -I{} -P "$DL_JOBS" bash -c 'fetch_one "$@"' _ {} || status=$?
(( status == 0 )) || die "one or more downloads failed - re-run to resume"

# -----------------------------------------------------------------------------
# 3. Point the samplesheet at the downloaded files
# -----------------------------------------------------------------------------
python3 - "$SAMPLESHEET" "$RAW_DIR" <<'PY'
import csv, sys, os
sheet, raw = sys.argv[1], sys.argv[2]
with open(sheet, newline='') as fh:
    rows = list(csv.DictReader(fh)); fields = list(rows[0].keys())
changed = 0
for r in rows:
    want = os.path.join(raw, f"{r['sample']}.fastq.gz")
    if r.get('fastq_1') != want:
        r['fastq_1'] = want; changed += 1
with open(sheet, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
print(f"[fastq] samplesheet: {changed} path(s) updated")
PY

log "done ($(du -sh "$RAW_DIR" | cut -f1)). Next:  bash scripts/04_run_pipeline.sh"
