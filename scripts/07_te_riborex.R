#!/usr/bin/env Rscript
#
#  Translational efficiency by riborex.
#
#  riborex (Li et al. 2017, Bioinformatics) is the lightest of the three. It
#  pools the RNA and RPF libraries into one generalised linear model and tests
#  the same interaction term DESeq2 does - but it hands the fit to an existing
#  engine (DESeq2, edgeR, or a modified edgeR) rather than reimplementing it.
#
#  So riborex and 05_te_deseq2.R --design simple are close to the same model.
#  Where they differ is the pairing: riborex has no way to express the per-
#  replicate `pair` term, so it always fits the unpaired form. Treat agreement
#  between them as a sanity check on the fit, not as independent evidence.
#
#  Usage:
#      Rscript 07_te_riborex.R
#      Rscript 07_te_riborex.R --engine edgeR
#      Rscript 07_te_riborex.R --counts <file> --outdir <dir>
#

source("te_common.R")
.te_require("riborex")
suppressPackageStartupMessages({
    library(riborex)
})

opts <- te_args()

# --engine is riborex-specific, so it is parsed here rather than in te_common.R
a <- commandArgs(trailingOnly = TRUE)
engine <- if ("--engine" %in% a) a[which(a == "--engine") + 1] else "DESeq2"
if (!engine %in% c("DESeq2", "edgeR", "edgeRD", "Voom")) {
    stop("--engine must be one of DESeq2, edgeR, edgeRD, Voom")
}

d <- load_te_data(opts)

# Pass the FACTOR, not as.character(). riborex feeds the condition vector
# straight into a data.frame; a character vector gets factorised with R's
# default alphabetical levels, which here makes NitrogenMinus the reference and
# silently reports the contrast BACKWARDS ("NitrogenPlus vs NitrogenMinus").
# The factor from load_te_data() already carries levels c(reference, target).
cond <- d$condition
message(sprintf("[riborex] engine: %s", engine))
message(sprintf("[riborex] levels: %s (reference first)",
                paste(levels(cond), collapse = " < ")))

stopifnot(identical(rownames(d$rna), rownames(d$ribo)))

# riborex wants data.frames, and it matches the two count tables positionally -
# already guaranteed by the pair-matching in load_te_data().
res <- riborex(rnaCntTable  = as.data.frame(d$rna),
               riboCntTable = as.data.frame(d$ribo),
               rnaCond      = cond,
               riboCond     = cond,
               engine       = engine)

# Confirm the direction actually tested matches contrasts.csv. Without this a
# future change to riborex or to the level order flips every sign unnoticed.
if (engine == "DESeq2") {
    desc <- mcols(res)$description
    desc <- desc[grep("log2 fold", desc)]
    want <- sprintf("%s vs %s", d$contrast$target, d$contrast$reference)
    if (length(desc) && !grepl(want, desc[1], fixed = TRUE)) {
        stop(sprintf("riborex tested '%s' but contrasts.csv asks for '%s'",
                     desc[1], want), call. = FALSE)
    }
    message(sprintf("[riborex] verified contrast: %s", want))
}

# Each engine returns its own object shape; normalise to one table.
tab <- as.data.frame(res)
lfc_col  <- grep("log2FoldChange|logFC",     colnames(tab), value = TRUE)[1]
p_col    <- grep("^pvalue$|^PValue$",        colnames(tab), value = TRUE)[1]
padj_col <- grep("^padj$|^FDR$|adj",         colnames(tab), value = TRUE)[1]
if (is.na(lfc_col) || is.na(p_col)) {
    stop("unexpected riborex output columns: ", paste(colnames(tab), collapse = ", "))
}
if (is.na(padj_col)) {
    tab$padj <- p.adjust(tab[[p_col]], method = "BH")
    padj_col <- "padj"
}

out <- data.frame(
    gene_id            = rownames(tab),
    log2FC_TE_riborex  = tab[[lfc_col]],
    pvalue_TE_riborex  = tab[[p_col]],
    padj_TE_riborex    = tab[[padj_col]],
    row.names          = NULL
)
if ("baseMean" %in% colnames(tab)) out$baseMean_counts <- tab$baseMean

e <- load_expression_data(opts)
ann <- expr_annotation(e$counts, e$lengths, e$meta)
out <- annotate_results(out, ann, "riborex")
out <- out[order(out$padj_TE_riborex, -abs(out$log2FC_TE_riborex)), ]

write_flat(out, opts, sprintf("%s.riborex_TE.tsv", d$name), tool = "riborex")
summarise_te(out$padj_TE_riborex, out$log2FC_TE_riborex, "riborex")
message(sprintf("[riborex] log2FC basis: %s", LOG2FC_BASIS[["riborex"]]))

saveRDS(res, file.path(tool_dir(opts, "riborex"),
                       sprintf("%s.riborex_TE_%s.rds", d$name, engine)))

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
f <- open_pdf(opts, sprintf("%s.riborex_TE_plots.pdf", d$name), tool = "riborex")

plot_volcano(out$log2FC_TE_riborex, out$padj_TE_riborex,
             sprintf("riborex (%s) TE volcano - %s", engine, d$name),
             xlab = "log2 fold change in TE")

if (!is.null(out$baseMean_counts)) {
    plot_ma(out$baseMean_counts, out$log2FC_TE_riborex, out$padj_TE_riborex,
            sprintf("riborex (%s) TE MA", engine), xlab = "mean normalised count")
}

ribo_tpm <- rowMeans(out[, grep("^TPM_Ribo_", names(out)), drop = FALSE], na.rm = TRUE)
plot_lfc_vs_expr(ribo_tpm, out$log2FC_TE_riborex, out$padj_TE_riborex,
                 "TE change vs Ribo expression", xlab = "mean Ribo TPM")

invisible(dev.off())
message(sprintf("[out] wrote %s", basename(f)))
message("[riborex] done")
