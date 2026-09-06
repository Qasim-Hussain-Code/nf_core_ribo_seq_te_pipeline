#!/usr/bin/env bash
# =============================================================================
#  02_prepare_genome.sh - fetch the reference and make the GTF usable
# =============================================================================
#  Resolves any NCBI assembly accession to its FTP path, downloads the FASTA
#  and GTF, then repairs the GTF. That repair is not cosmetic - see fix_gtf.py.
# =============================================================================
STEP=genome
source "$(dirname "${BASH_SOURCE[0]:-$0}")/lib_common.sh"
mkdirs
cd "$PROJECT_DIR"

ACC="${1:-$GENOME_ACC}"
log "assembly accession: ${ACC}"

# -----------------------------------------------------------------------------
# 1. Resolve the accession to an FTP directory via the NCBI datasets API
# -----------------------------------------------------------------------------
# NCBI's layout is  .../GCF/000/007/065/GCF_000007065.1_ASM706v1/  - the digits
# are split into groups of three, and the suffix is the assembly NAME, which
# cannot be derived from the accession alone. Ask the API for it.
if [[ ! -s genome.fasta || ! -s genome.gtf ]]; then
    log "resolving FTP path from NCBI ..."
    PREFIX="${ACC%%_*}"                      # GCF or GCA
    DIGITS="${ACC#*_}"; DIGITS="${DIGITS%%.*}"
    P1="${DIGITS:0:3}"; P2="${DIGITS:3:3}"; P3="${DIGITS:6:3}"
    BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/${PREFIX}/${P1}/${P2}/${P3}"

    ASM_DIR="$(curl -sSL "${BASE}/" | grep -o "${ACC}_[^\"/]*" | head -1)"
    [[ -n "$ASM_DIR" ]] || die "could not resolve ${ACC} under ${BASE}/ - check the accession"
    URL="${BASE}/${ASM_DIR}/${ASM_DIR}"
    log "resolved: ${ASM_DIR}"
fi

fetch() {  # fetch <suffix> <outfile>
    local suffix="$1" out="$2"
    [[ -s "$out" ]] && { log "${out} present - skipping"; return; }
    log "downloading ${out} ..."
    curl -sSL --retry 5 --retry-delay 5 "${URL}_${suffix}" -o "${out}.gz" \
        || die "download failed: ${URL}_${suffix}"
    gunzip -f "${out}.gz"
}

fetch "genomic.fna.gz" genome.fasta
fetch "genomic.gtf.gz" genome.gtf

[[ -s genome.fasta ]] || die "genome.fasta empty"
[[ -s genome.gtf   ]] || die "genome.gtf empty"

GENOME_BP=$(grep -v '^>' genome.fasta | tr -d '\n' | wc -c)
log "genome: $(grep -c '^>' genome.fasta) sequence(s), ${GENOME_BP} bp"

# -----------------------------------------------------------------------------
# 2. Repair the GTF
# -----------------------------------------------------------------------------
if [[ ! -s genome.exons.gtf || genome.gtf -nt genome.exons.gtf ]]; then
    python3 "${GUIDE_DIR}/scripts/fix_gtf.py" genome.gtf genome.exons.gtf \
        --pad "${GTF_PAD}" --fasta genome.fasta
else
    log "genome.exons.gtf up to date - skipping"
fi

# -----------------------------------------------------------------------------
# 3. STAR's suffix-array parameter depends on genome size
# -----------------------------------------------------------------------------
# STAR's default --genomeSAindexNbases 14 suits a mammalian genome. On a 4 Mb
# bacterial genome it over-allocates and the index is unusable. The formula is
# min(14, log2(genomeLength)/2 - 1).
SA_N=$(python3 -c "import math;print(min(14,max(4,int(math.log2($GENOME_BP)/2-1))))")
log "STAR --genomeSAindexNbases ${SA_N} (from ${GENOME_BP} bp)"
grep -q '^SA_INDEX_NBASES=' "$CONF" && sed -i "/^SA_INDEX_NBASES=/d" "$CONF"
echo "SA_INDEX_NBASES=${SA_N}" >> "$CONF"

log "done. Next:  bash scripts/03_download_fastq.sh"
