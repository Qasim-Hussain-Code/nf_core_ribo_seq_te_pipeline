## Results

This directory holds the analysis output that is small enough to version
directly and useful to inspect without re-running anything.

| Directory | Contents |
|---|---|
| `differential_expression/` | Translational-efficiency tables from all three methods, RNA-only and Ribo-only differential expression, normalised expression matrices, and the self-contained HTML report |
| `quantification/` | Merged Salmon gene- and transcript-level count and TPM matrices across all eight libraries |
| `orf_predictions/` | Genome-wide open reading frame calls from RiboTish |
| `ribotricer/` | Actively translated ORF calls from RiboTricer |
| `riboseq_qc/` | Per-library ribosome-profiling quality metrics (read-length distribution, metagene coverage, RiboTish quality scores) |

A handful of larger, purely derivative directories (raw alignments,
reference indices, per-library QC reports, the full riboWaltz output) are
excluded for size; each has its own short placeholder explaining what
belongs there and how to regenerate it.
