#!/usr/bin/env Rscript
#
#  Normalised expression matrices for every RNA-seq and Ribo-seq library.
#
#  Which normalisation to use for what:
#
#    CPM     - depth only. Fine for a quick look, WRONG for comparing different
#              genes to each other because it ignores gene length.
#    RPKM    - depth + length. Its per-sample totals are not equal, so RPKM
#    /FPKM     values are not strictly comparable BETWEEN samples. For our
#              single-end libraries RPKM and FPKM are numerically identical:
#              "fragments" and "reads" are the same thing when reads are not
#              paired. Both are written so downstream tools find the name they
#              expect, but they are the same numbers.
#    TPM     - length first, then depth. Columns sum to 1e6, so TPM IS
#              comparable between samples. Prefer this for expression levels.
#    DESeq2  - median-of-ratios size factors. What the DE models actually use.
#    VST     - variance-stabilised, homoscedastic. For PCA, clustering and
#              heatmaps - NOT for reporting expression levels.
#
#  None of these are for differential testing: DE runs on raw counts (see
#  09_de_rna_ribo.R). Normalised matrices are for visualising and reporting.
#
#  Usage:
#      Rscript 08_expression_matrices.R
#      Rscript 08_expression_matrices.R --outdir <dir> --min_count 20
#

source("te_common.R")
suppressPackageStartupMessages({
    library(DESeq2)
})

opts <- te_args()   # defaults to results/differential_expression - one flat folder
d <- load_expression_data(opts)

counts <- d$counts
if (is.null(d$lengths)) {
    stop("gene lengths not found - RPKM/FPKM/TPM need ", opts$lengths, call. = FALSE)
}
# Effective length from salmon; guard against zeros before dividing.
len <- d$lengths
len[len <= 0] <- NA

EXPR_DIR <- tool_dir(opts, "expression")
p <- function(f) file.path(EXPR_DIR, f)
nm <- function(f) sprintf("%s.%s", d$name, f)

# ---------------------------------------------------------------------------
# CPM - depth normalised only
# ---------------------------------------------------------------------------
lib <- colSums(counts)
cpm <- t(t(counts) / lib) * 1e6
write_matrix(round(cpm, 4), p(nm("gene_CPM.tsv")))

# ---------------------------------------------------------------------------
# RPKM / FPKM - depth then length
# ---------------------------------------------------------------------------
rpkm <- cpm / (len / 1000)
write_matrix(round(rpkm, 4), p(nm("gene_RPKM.tsv")))
# Identical for single-end data; written under both names for convenience.
write_matrix(round(rpkm, 4), p(nm("gene_FPKM.tsv")))

# ---------------------------------------------------------------------------
# TPM - length then depth, so columns sum to 1e6
# ---------------------------------------------------------------------------
rate <- counts / (len / 1000)          # reads per kilobase
tpm  <- t(t(rate) / colSums(rate, na.rm = TRUE)) * 1e6
write_matrix(round(tpm, 4), p(nm("gene_TPM.tsv")))

stopifnot(all(abs(colSums(tpm, na.rm = TRUE) - 1e6) < 1))
message("[expr] TPM columns sum to 1e6 - OK")

# ---------------------------------------------------------------------------
# DESeq2 size-factor normalised counts, and VST for plotting
# ---------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(round(counts), d$meta, ~ assay + condition)
dds <- estimateSizeFactors(dds)
write_matrix(round(counts(dds, normalized = TRUE), 4), p(nm("gene_DESeq2_normalised.tsv")))

sf <- data.frame(sample = names(sizeFactors(dds)),
                 size_factor = round(as.numeric(sizeFactors(dds)), 4),
                 library_size = lib[names(sizeFactors(dds))])
write.table(sf, p(nm("size_factors.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
message("[expr] wrote ", basename(p(nm("size_factors.tsv"))))

# blind = TRUE: the transform must not learn from the design, or downstream
# clustering looks better than the data warrants.
vsd <- vst(dds, blind = TRUE)
write_matrix(round(assay(vsd), 4), p(nm("gene_VST.tsv")))

# ---------------------------------------------------------------------------
# QC plots on the VST matrix
# ---------------------------------------------------------------------------
pdf(p(nm("expression_qc.pdf")), width = 8, height = 7)

pca <- plotPCA(vsd, intgroup = c("assay", "condition"), returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"))
print(
    ggplot2::ggplot(pca, ggplot2::aes(PC1, PC2, colour = condition, shape = assay)) +
        ggplot2::geom_point(size = 4) +
        ggplot2::labs(x = sprintf("PC1: %d%% variance", pv[1]),
                      y = sprintf("PC2: %d%% variance", pv[2]),
                      title = "PCA of VST-transformed counts") +
        ggplot2::theme_bw(base_size = 12)
)

d_mat <- as.matrix(dist(t(assay(vsd))))
heatmap(d_mat, symm = TRUE, margins = c(12, 12),
        main = "Sample-to-sample distance (VST)")

# Ribo:RNA ratio per replicate - a crude per-sample TE sanity check.
# Match through `pair` explicitly: positional matching here would silently
# divide one culture's footprints by another culture's mRNA.
ribo_meta <- d$meta[d$meta$assay == "Ribo", ]
rna_meta  <- d$meta[d$meta$assay == "RNA",  ]
rna_for_ribo <- rna_meta$sample[match(ribo_meta$pair, rna_meta$pair)]
stopifnot(!anyNA(rna_for_ribo))

te_proxy <- log2((tpm[, ribo_meta$sample, drop = FALSE] + 1) /
                 (tpm[, rna_for_ribo,     drop = FALSE] + 1))
colnames(te_proxy) <- sprintf("%s/%s", ribo_meta$sample, rna_for_ribo)
boxplot(te_proxy, las = 2, ylab = "log2 (Ribo TPM + 1) / (RNA TPM + 1)",
        main = "Per-replicate TE proxy", cex.axis = 0.55)
abline(h = 0, lty = 2, col = "grey40")

# Library size and per-assay expression distribution
barplot(lib / 1e6, las = 2, cex.names = 0.6, ylab = "million assigned reads",
        main = "Library size", col = ifelse(d$meta$assay == "Ribo", "firebrick", "steelblue"))
legend("topright", c("Ribo", "RNA"), fill = c("firebrick", "steelblue"), bty = "n")

boxplot(log2(tpm + 1), las = 2, cex.axis = 0.55, ylab = "log2(TPM + 1)",
        main = "TPM distribution per library",
        col = ifelse(d$meta$assay == "Ribo", "firebrick", "steelblue"))

# TPM vs RPKM on one library, to make the difference concrete
s1 <- colnames(tpm)[1]
plot(pmax(rpkm[, s1], 1e-2), pmax(tpm[, s1], 1e-2), log = "xy", pch = 19, cex = 0.3,
     col = "grey40", xlab = sprintf("RPKM (%s)", s1), ylab = sprintf("TPM (%s)", s1),
     main = "TPM vs RPKM - same ranking, different scaling")
abline(0, 1, lty = 2, col = "firebrick")

invisible(dev.off())
message("[expr] wrote ", basename(p(nm("expression_qc.pdf"))))

message(sprintf("[expr] done - matrices in %s", EXPR_DIR))
