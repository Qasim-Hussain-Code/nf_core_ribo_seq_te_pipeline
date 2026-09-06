#!/usr/bin/env Rscript
#
#  Translational efficiency by Xtail.
#
#  Xtail (Xiao et al. 2016, Nat Commun) does not fit a single interaction term.
#  It evaluates TE change two ways and keeps the more conservative:
#
#    1. the difference between the mRNA log2FC and the RPF log2FC
#       (how differently the two assays respond to the condition), and
#    2. the difference between the within-condition RPF:mRNA log2 ratios
#       (how the TE itself differs between conditions)
#
#  It takes whichever gives the LARGER p-value, which is why Xtail is usually
#  the strictest of the three methods here. That conservatism is the reason it
#  is recommended at low replication - it is built to avoid the false positives
#  a 2-vs-2 comparison invites.
#
#  Xtail is the slow one. Two knobs control that, both exposed below:
#
#    --threads N  parallelises the per-gene loops. Setup is serial, so scaling
#                 is sub-linear - do not expect 2x cores to halve the runtime.
#    --bins N     resolution of the numerical integration xtail runs per gene,
#                 and the dominant cost. Runtime scales roughly linearly with
#                 it. 10000 is the package default and the value used in the
#                 paper; 2000 is ~5x faster and moves log2FC_TE barely at all,
#                 though borderline padj values shift. Use a low value while
#                 iterating and 10000 for numbers you intend to publish.
#
#  Usage:
#      Rscript 06_te_xtail.R
#      Rscript 06_te_xtail.R --threads 10
#      Rscript 06_te_xtail.R --threads 10 --bins 2000     # fast pass
#      Rscript 06_te_xtail.R --counts <file> --outdir <dir>
#

source("te_common.R")
.te_require("xtail")
suppressPackageStartupMessages({
    library(xtail)
})

opts <- te_args()

# --threads / --bins are xtail-specific, so they are parsed here rather than in
# te_common.R. TE_THREADS stays supported as a fallback for existing scripts.
a <- commandArgs(trailingOnly = TRUE)
arg_int <- function(flag, fallback) {
    if (flag %in% a) {
        v <- suppressWarnings(as.integer(a[which(a == flag) + 1]))
        if (is.na(v) || v < 1) stop(sprintf("%s needs a positive integer", flag), call. = FALSE)
        v
    } else fallback
}

threads <- arg_int("--threads", as.integer(Sys.getenv("TE_THREADS", "4")))
bins    <- arg_int("--bins", 10000)

ncores <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
if (!is.na(ncores) && threads > ncores) {
    warning(sprintf("--threads %d exceeds the %d cores available; using %d",
                    threads, ncores, ncores), call. = FALSE)
    threads <- ncores
}
message(sprintf("[xtail] threads: %d%s", threads,
                if (!is.na(ncores)) sprintf(" of %d cores", ncores) else ""))
message(sprintf("[xtail] bins: %d%s", bins,
                if (bins < 10000) "  (reduced - faster, coarser p-values)" else ""))

d <- load_te_data(opts)

cond <- as.character(d$condition)

# xtail's `baseLevel` defaults to condition[1] - whichever condition happens to
# sit in the first column. Our columns are ordered with the TARGET first, so the
# default would make the target the base and report every fold change backwards.
# Pin it to the reference from contrasts.csv instead.
base_level <- d$contrast$reference
stopifnot(base_level %in% cond)
message(sprintf("[xtail] condition vector: %s", paste(cond, collapse = ", ")))
message(sprintf("[xtail] baseLevel (reference): %s -> log2FC_TE is %s vs %s",
                base_level, d$contrast$target, d$contrast$reference))

# xtail() wants mRNA and RPF matrices with identical gene order and columns in
# the same sample order as `condition`. load_te_data() has already matched them
# through the `pair` column, but assert it rather than trust it.
stopifnot(identical(rownames(d$rna), rownames(d$ribo)))
stopifnot(ncol(d$rna) == length(cond), ncol(d$ribo) == length(cond))

res <- xtail(mrna      = d$rna,
             rpf       = d$ribo,
             condition = cond,
             baseLevel = base_level,
             bins      = bins,
             minMeanCount = 1,
             threads   = threads)

tab <- resultsTable(res, log2FCs = TRUE, log2Rs = TRUE)

out <- data.frame(
    gene_id           = rownames(tab),
    log2FC_TE_xtail   = tab[["log2FC_TE_final"]],
    pvalue_TE_xtail   = tab[["pvalue_final"]],
    padj_TE_xtail     = tab[["pvalue.adjust"]],
    log2FC_mRNA_xtail = tab[[grep("^mRNA_log2FC$", colnames(tab))[1]]],
    log2FC_RPF_xtail  = tab[[grep("^RPF_log2FC$",  colnames(tab))[1]]],
    row.names         = NULL
)

e <- load_expression_data(opts)
ann <- expr_annotation(e$counts, e$lengths, e$meta)
out <- annotate_results(out, ann, "xtail")
out <- out[order(out$padj_TE_xtail, -abs(out$log2FC_TE_xtail)), ]

write_flat(out, opts, sprintf("%s.xtail_TE.tsv", d$name), tool = "xtail")
summarise_te(out$padj_TE_xtail, out$log2FC_TE_xtail, "xtail")
message(sprintf("[xtail] log2FC basis: %s", LOG2FC_BASIS[["xtail"]]))

saveRDS(res, file.path(tool_dir(opts, "xtail"), sprintf("%s.xtail_TE.rds", d$name)))

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
f <- open_pdf(opts, sprintf("%s.xtail_TE_plots.pdf", d$name), tool = "xtail")

plot_volcano(out$log2FC_TE_xtail, out$padj_TE_xtail,
             sprintf("Xtail TE volcano - %s", d$name),
             xlab = "log2 fold change in TE")

# Xtail's own view: mRNA fold change against RPF fold change. Points off the
# diagonal are where transcription and translation disagree.
tryCatch(plotFCs(res),
         error = function(e) message("[xtail] plotFCs skipped: ", conditionMessage(e)))

# plotRs() raises "subscript contains invalid names" on this data: xtail names
# BOTH per-condition log2TE columns after the base level, so resultsTable gives
# "<ref>_log2TE" and "<ref>_log2TE.1" and plotRs cannot find the second. The
# numbers are correct, only the labels are wrong - draw it directly instead.
te_cols <- grep("_log2TE", colnames(tab), value = TRUE)
if (length(te_cols) >= 2) {
    ok <- is.finite(tab[[te_cols[1]]]) & is.finite(tab[[te_cols[2]]])
    plot(tab[[te_cols[1]]][ok], tab[[te_cols[2]]][ok], pch = 19, cex = 0.4,
         col = "grey40",
         xlab = sprintf("log2 TE in %s", d$contrast$reference),
         ylab = sprintf("log2 TE in %s", d$contrast$target),
         main = "Within-condition log2 TE (xtail)")
    abline(0, 1, lty = 2, col = "firebrick")
} else {
    message("[xtail] no per-condition log2TE columns - R plot skipped")
}

# mRNA vs RPF fold change, drawn from our own columns so it renders even if
# xtail's plotFCs fails.
ok <- is.finite(out$log2FC_mRNA_xtail) & is.finite(out$log2FC_RPF_xtail)
sig <- !is.na(out$padj_TE_xtail) & out$padj_TE_xtail < 0.05
plot(out$log2FC_mRNA_xtail[ok], out$log2FC_RPF_xtail[ok], pch = 19, cex = 0.4,
     col = ifelse(sig[ok], "firebrick", "grey75"),
     xlab = "mRNA log2FC", ylab = "RPF log2FC",
     main = sprintf("mRNA vs RPF fold change - %s", d$name))
abline(0, 1, lty = 2, col = "grey40"); abline(h = 0, v = 0, lty = 3, col = "grey70")
legend("topleft", bty = "n", cex = 0.75, pch = 19, col = c("firebrick", "grey75"),
       legend = c("TE padj < 0.05", "ns"))

ribo_tpm <- rowMeans(out[, grep("^TPM_Ribo_", names(out)), drop = FALSE], na.rm = TRUE)
plot_lfc_vs_expr(ribo_tpm, out$log2FC_TE_xtail, out$padj_TE_xtail,
                 "TE change vs Ribo expression", xlab = "mean Ribo TPM")

invisible(dev.off())
message(sprintf("[out] wrote %s", basename(f)))
message("[xtail] done")
