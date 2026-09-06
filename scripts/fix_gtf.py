#!/usr/bin/env python3
"""
Add exon (and transcript) features to an NCBI prokaryotic/archaeal GTF.

NCBI does not emit `exon` records for protein-coding genes in its bacterial and
archaeal GTFs - only CDS, start_codon and stop_codon. Most tools cope, but
riboWaltz builds a TxDb via GenomicFeatures::makeTxDbFromGRanges, which requires
every CDS to sit inside an exon of the same transcript. Without them it dies with:

    Error in .make_splicings(exons, cds, stop_codons) :
      some CDS parts cannot be mapped to an exon

Prokaryotes/archaea have no introns, so a transcript's single exon is simply the
span of its own features. The stop codon is included: NCBI's CDS ends just before
it, and riboWaltz needs the stop inside the transcript to place the CDS end.

The exon is also padded by --pad nt on each side to create synthetic UTRs. This
is not cosmetic - riboWaltz's psite() selects the reads it calibrates P-site
offsets from with:

    site_dist_end5 = end5 - cds_start
    site_sub <- dt[site_dist_end5 <= -flanking & site_dist_end3 >= flanking - 1]

so a read only counts if its 5' end lies at least `flanking` (default 6) nt
UPSTREAM of the CDS start. With an exon equal to the CDS, cds_start is 1 for
every transcript, no read can qualify, and psite() dies on an empty table with
"argument is of length zero". The pad must comfortably exceed a footprint's
5' overhang (~12 nt for the P-site), hence the 50 nt default.

Padding is applied only to synthesised exons; real tRNA/rRNA exons already in
the annotation are left untouched. Coordinates are clamped to the sequence ends.

Finally, spaces in the source column are collapsed to underscores. STAR's GTF
parser tokenises on whitespace, so NCBI's "Protein Homology" source shifts the
feature-type field and every such line is SILENTLY skipped - STAR builds a
transcriptome of 330 instead of 3517 transcripts, and salmon then quantifies
only ~9% of genes with no error anywhere.

Usage: fix_gtf.py <in.gtf> <out.gtf> [--pad N] [--fasta genome.fasta]
"""
import re
import sys
from collections import OrderedDict

TX_RE = re.compile(r'transcript_id "([^"]+)"')
GENE_RE = re.compile(r'gene_id "([^"]+)"')

# Order features within a transcript the way tools expect to read them.
RANK = {"transcript": 0, "exon": 1, "CDS": 2, "start_codon": 3, "stop_codon": 4}


def seq_lengths(fasta):
    """Sequence name -> length, so padded exons can be clamped to the contig."""
    lengths, name, n = {}, None, 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def main(src, dst, pad=50, fasta=None):
    limits = seq_lengths(fasta) if fasta else {}
    header = []
    records = []                 # (seqname, start, end, feature, transcript, line)
    tx = OrderedDict()           # transcript_id -> collected info
    saw_exon, saw_transcript = set(), set()
    despaced = 0

    with open(src) as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue

            # STAR's GTF parser tokenises on whitespace, so a source field
            # containing a space (NCBI uses "Protein Homology") shifts the
            # feature-type column and the line is silently skipped - STAR reads
            # 330 of 3517 transcripts from this annotation. Collapse it.
            if " " in f[1]:
                f[1] = f[1].replace(" ", "_")
                line = "\t".join(f) + "\n"
                despaced += 1

            seqname, source, feature, start, end, _, strand, _, attrs = f
            start, end = int(start), int(end)

            m = TX_RE.search(attrs)
            if not m:
                records.append((seqname, start, end, feature, "", line))
                continue
            tid = m.group(1)

            if feature == "exon":
                saw_exon.add(tid)
            elif feature == "transcript":
                saw_transcript.add(tid)

            info = tx.get(tid)
            if info is None:
                gene = GENE_RE.search(attrs)
                tx[tid] = {
                    "seqname": seqname,
                    "source": source,
                    "strand": strand,
                    "gene": gene.group(1) if gene else tid,
                    "start": start,
                    "end": end,
                }
            else:
                info["start"] = min(info["start"], start)
                info["end"] = max(info["end"], end)

            records.append((seqname, start, end, feature, tid, line))

    added_exon = added_tx = 0
    for tid, i in tx.items():
        if tid in saw_exon and tid in saw_transcript:
            continue
        attrs = f'gene_id "{i["gene"]}"; transcript_id "{tid}";'
        # Synthetic UTRs, clamped to the contig ends.
        pstart = max(1, i["start"] - pad)
        pend = i["end"] + pad
        limit = limits.get(i["seqname"])
        if limit:
            pend = min(pend, limit)
        base = (i["seqname"], i["source"], pstart, pend, i["strand"])

        if tid not in saw_transcript:
            records.append((
                i["seqname"], pstart, pend, "transcript", tid,
                "{}\t{}\ttranscript\t{}\t{}\t.\t{}\t.\t{}\n".format(
                    base[0], base[1], base[2], base[3], base[4], attrs)
            ))
            added_tx += 1

        if tid not in saw_exon:
            records.append((
                i["seqname"], pstart, pend, "exon", tid,
                '{}\t{}\texon\t{}\t{}\t.\t{}\t.\t{} exon_number "1";\n'.format(
                    base[0], base[1], base[2], base[3], base[4], attrs)
            ))
            added_exon += 1
        i["start"], i["end"] = pstart, pend

    # Keep each transcript's records contiguous and coordinate-ordered.
    tx_start = {t: i["start"] for t, i in tx.items()}
    records.sort(key=lambda r: (
        r[0],
        tx_start.get(r[4], r[1]),
        r[4],
        RANK.get(r[3], 9),
        r[1],
        r[2],
    ))

    with open(dst, "w") as out:
        out.writelines(header)
        out.writelines(r[5] for r in records)

    print(f"[gtf] transcripts seen      : {len(tx)}")
    print(f"[gtf] exon records added    : {added_exon}")
    print(f"[gtf] transcript records add: {added_tx}")
    print(f"[gtf] UTR pad each side     : {pad} nt")
    print(f"[gtf] source fields despaced: {despaced}")
    print(f"[gtf] wrote {dst}")


if __name__ == "__main__":
    argv = sys.argv[1:]
    pad, fasta = 50, None
    for flag, cast in (("--pad", int), ("--fasta", str)):
        if flag in argv:
            k = argv.index(flag)
            val = cast(argv[k + 1])
            if flag == "--pad":
                pad = val
            else:
                fasta = val
            del argv[k:k + 2]
    if len(argv) != 2:
        sys.exit(__doc__)
    main(argv[0], argv[1], pad=pad, fasta=fasta)
