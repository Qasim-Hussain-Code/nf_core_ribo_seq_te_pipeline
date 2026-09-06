#!/usr/bin/env Rscript
#
#  Translational efficiency by DESeq2 interaction model.
#
#  TE is the ratio of ribosome footprints to mRNA. A change in TE between two
#  conditions is therefore a change in that ratio - which is exactly an
#  interaction term. Stacking the Ribo and RNA libraries into one matrix and
#  fitting
#
#      ~ condition + condition:ind + assay + condition:assay
#
#  makes `conditionTarget.assayRibo` the log2 fold change in TE. This is
#  DESeq2's "individuals nested within groups" design, with `ind` the replicate
#  number WITHIN each condition - it pairs each culture's Ribo library to its
#  own RNA library instead of the group mean.
#
#  The obvious `~ pair + assay + assay:condition` is NOT full rank: without the
#  `condition` main effect R emits one interaction column per assay level, and
#  pair3 + pair4 equals their sum exactly. --design simple drops the pairing
#  entirely, which is what Xtail and riborex fit internally.
#
#  log2FC comes from RAW COUNTS via DESeq2's median-of-ratios normalisation -
#  never from RPKM/FPKM/TPM. Those are attached to the output as descriptive
#  expression levels only; see LOG2FC_BASIS in te_common.R.
#
#  Usage:
#      Rscript 05_te_deseq2.R
#      Rscript 05_te_deseq2.R --design simple
#      Rscript 05_te_deseq2.R --counts <file> --outdir <dir>
#

source("te_common.R")
suppressPackageStartupMessages({
    library(DESeq2)
})

opts <- te_args()
if (is.null(opts$design)) opts$design <- "paired"
d <- load_te_data(opts)

# ---------------------------------------------------------------------------
# Stack Ribo and RNA into a single matrix, one column per library
# ---------------------------------------------------------------------------
counts <- cbind(d$ribo, d$rna)
colnames(counts) <- c(paste0("Ribo_", d$ribo_meta$pair),
                      paste0("RNA_",  d$rna_meta$pair))

coldata <- data.frame(
    assay     = factor(rep(c("Ribo", "RNA"), each = ncol(d$ribo)), levels = c("RNA", "Ribo")),
    condition = factor(rep(as.character(d$condition), 2), levels = levels(d$condition)),
    pair      = factor(rep(as.character(d$pair), 2)),
    # Replicate number WITHIN condition (1,2 in each), not the global pair id.
    # This is what makes the nested paired design below estimable.
    ind       = factor(rep(as.character(d$ribo_meta$replicate), 2)),
    row.names = colnames(counts)
)

if (opts$design == "paired") {
    # DESeq2's "individuals nested within groups" form. The obvious
    # `~ pair + assay + assay:condition` is NOT full rank here: dropping the
    # `condition` main effect makes R emit an interaction column per assay
    # level, and pair3 + pair4 equals their sum exactly (rank 6 of 7). Nesting
    # the within-condition replicate `ind` instead gives 6 columns of rank 6.
    design <- ~ condition + condition:ind + assay + condition:assay
} else if (opts$design == "simple") {
    # Unpaired - the model riborex and Xtail fit internally.
    design <- ~ condition + assay + condition:assay
} else {
    stop("--design must be 'paired' or 'simple'")
}
message(sprintf("[deseq2] design: %s", paste(deparse(design), collapse = "")))
print(coldata)

# Check rank BEFORE handing the design to DESeq2, so a bad design produces an
# actionable message instead of DESeq2's generic checkFullRank error.
mm <- model.matrix(design, coldata)
if (qr(mm)$rank < ncol(mm)) {
    stop(sprintf(paste0("design '%s' is not full rank (%d columns, rank %d).\n",
                        "  Re-run with --design simple for an unpaired fit."),
                 opts$design, ncol(mm), qr(mm)$rank), call. = FALSE)
}
message(sprintf("[deseq2] %d samples, %d coefficients, %d residual df",
                nrow(coldata), ncol(mm), nrow(coldata) - ncol(mm)))

dds <- DESeqDataSetFromMatrix(counts, coldata, design)

dds <- DESeq(dds)

# The interaction coefficient IS the TE change. Pick it by name so a change in
# factor levels can never silently select the wrong contrast.
target <- levels(d$condition)[2]
wanted <- grep(sprintf("^assayRibo\\.condition%s$|^condition%s\\.assayRibo$",
                       target, target),
               resultsNames(dds), value = TRUE)
if (length(wanted) != 1) {
    stop("could not identify the interaction coefficient; available: ",
         paste(resultsNames(dds), collapse = ", "))
}
message(sprintf("[deseq2] TE coefficient: %s", wanted))

res <- results(dds, name = wanted, alpha = 0.05)
summary(res)

out <- data.frame(
    gene_id            = rownames(res),
    log2FC_TE_deseq2   = res$log2FoldChange,
    lfcSE_TE_deseq2    = res$lfcSE,
    stat_TE_deseq2     = res$stat,
    pvalue_TE_deseq2   = res$pvalue,
    padj_TE_deseq2     = res$padj,
    baseMean_counts    = res$baseMean,
    row.names          = NULL
)

# Attach the expression levels the fold change should be read against, plus an
# explicit record of what the log2FC was actually computed from.
e <- load_expression_data(opts)
ann <- expr_annotation(e$counts, e$lengths, e$meta)
out <- annotate_results(out, ann, "deseq2")
out <- out[order(out$padj_TE_deseq2, -abs(out$log2FC_TE_deseq2)), ]

write_flat(out, opts, sprintf("%s.deseq2_TE.tsv", d$name), tool = "deseq2")
summarise_te(out$padj_TE_deseq2, out$log2FC_TE_deseq2, "deseq2")
message(sprintf("[deseq2] log2FC basis: %s", LOG2FC_BASIS[["deseq2"]]))

saveRDS(dds, file.path(tool_dir(opts, "deseq2"), sprintf("%s.deseq2_TE_dds.rds", d$name)))

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
f <- open_pdf(opts, sprintf("%s.deseq2_TE_plots.pdf", d$name), tool = "deseq2")

plot_volcano(out$log2FC_TE_deseq2, out$padj_TE_deseq2,
             sprintf("DESeq2 TE volcano - %s", d$name),
             xlab = "log2 fold change in TE")
plot_ma(out$baseMean_counts, out$log2FC_TE_deseq2, out$padj_TE_deseq2,
        sprintf("DESeq2 TE MA - %s", d$name), xlab = "mean normalised count")

# TE change against the Ribo abundance it rests on - a large fold change on a
# gene with almost no footprints is not the same finding as one on an abundant
# gene.
ribo_tpm <- rowMeans(out[, grep("^TPM_Ribo_", names(out)), drop = FALSE], na.rm = TRUE)
plot_lfc_vs_expr(ribo_tpm, out$log2FC_TE_deseq2, out$padj_TE_deseq2,
                 "TE change vs Ribo expression", xlab = "mean Ribo TPM")

plotDispEsts(dds, main = "DESeq2 dispersion estimates")
invisible(dev.off())
message(sprintf("[out] wrote %s", basename(f)))
message("[deseq2] done")
