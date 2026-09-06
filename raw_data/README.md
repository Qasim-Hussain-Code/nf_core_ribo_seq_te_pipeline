## Raw sequencing reads

This directory holds the eight single-end FASTQ files listed in
`config/samplesheet.csv` (four RNA-seq libraries, four Ribo-seq libraries,
two conditions, two replicates each). They are not tracked in this
repository because of their size (~600 MB per file).

The reads are public and were retrieved from the ENA FTP mirror of
BioProject PRJNA1004312. Running `scripts/03_download_fastq.sh` will fetch
them here directly from the run accessions already present in the `srr`
column of the samplesheet, verify each file against its published md5, and
populate `config/samplesheet.csv` with the resulting local paths.
