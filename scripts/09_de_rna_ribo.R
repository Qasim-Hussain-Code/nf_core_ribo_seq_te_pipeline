#!/usr/bin/env Rscript
#
#  Differential expression WITHIN each assay, then the combined picture.
#
#  Three questions, three answers:
#
#    RNA  DE  - did transcript abundance change?          (transcriptional response)
#    Ribo DE  - did ribosome occupancy change?            (translational output)
#    TE       - did the Ribo:RNA ratio change?            (translational control)
#
#  TE alone is not interpretable on its own. A gene can show no TE change while
#  its expression doubles at both levels, and a gene can show a large TE change
#  purely because its mRNA collapsed. Running RNA and Ribo DE alongside the TE
#  results is what separates those cases - and the classification below is
#  exactly the regulatory-mode call that anota2seq would have produced if it
#  could run at n=2:
#
#    Transcriptional  RNA and Ribo move together, TE unchanged
#    Translational    Ribo moves, RNA does not (TE change drives it)
#    Buffered         RNA moves, Ribo does not (translation absorbs it)
#    Intensified      both move the same way AND TE amplifies it
#
#  Usage:
#      Rscript 09_de_rna_ribo.R
#      Rscript 09_de_rna_ribo.R --te results/translational_efficiency/deseq2/<file>.tsv
#      Rscript 09_de_rna_ribo.R --outdir <dir> --alpha 0.1
#

source("te_common.R")
suppressPackageStartupMessages({
    library(DESeq2)
})

opts <- te_args()   # defaults to results/differential_expression - one flat folder

a <- commandArgs(trailingOnly = TRUE)
alpha  <- if ("--alpha" %in% a) as.numeric(a[which(a == "--alpha") + 1]) else 0.05
te_file <- if ("--te" %in% a) a[which(a == "--te") + 1] else NULL

d <- load_expression_data(opts)
dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# DE within one assay
# ---------------------------------------------------------------------------
run_de <- function(assay_name) {
    keep <- d$meta$assay == assay_name
    meta <- droplevels(d$meta[keep, , drop = FALSE])
    cts  <- round(d$counts[, meta$sample, drop = FALSE])

    # Drop genes that are silent within THIS assay. A gene with no footprints
    # at all should not carry a dispersion estimate in the Ribo model just
    # because its mRNA is abundant.
    cts <- cts[rowSums(cts) > opts$min_count, , drop = FALSE]

    message(sprintf("\n[de:%s] %d genes, %d samples (%s)", assay_name, nrow(cts),
                    nrow(meta), paste(levels(meta$condition), collapse = " -> ")))

    dds <- DESeq(DESeqDataSetFromMatrix(cts, meta, ~ condition))

    coef <- sprintf("condition_%s_vs_%s", d$contrast$target, d$contrast$reference)
    if (!coef %in% resultsNames(dds)) {
        stop(sprintf("expected coefficient '%s'; got: %s", coef,
                     paste(resultsNames(dds), collapse = ", ")), call. = FALSE)
    }
    res <- results(dds, name = coef, alpha = alpha)

    tag <- tolower(assay_name)
    out <- data.frame(gene_id = rownames(res),
                      log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
                      pvalue = res$pvalue, padj = res$padj,
                      baseMean_counts = res$baseMean, row.names = NULL)
    # Tag every statistic with the assay it came from, so the columns stay
    # unambiguous once everything is merged into the master table.
    names(out)[-1] <- sprintf("%s_%s", names(out)[-1], assay_name)

    # RPKM/FPKM/TPM for THIS assay only - the levels the fold change describes.
    ann <- expr_annotation(d$counts, d$lengths, d$meta, assays = assay_name)
    out <- annotate_results(out, ann, tag)
    names(out)[names(out) == "log2FC_basis"] <- sprintf("log2FC_basis_%s", assay_name)

    lfc <- out[[sprintf("log2FC_%s", assay_name)]]
    padj <- out[[sprintf("padj_%s", assay_name)]]
    out <- out[order(padj, -abs(lfc)), ]

    write_flat(out, opts, sprintf("%s.%s_DE.tsv", d$name, tag), tool = "rna_ribo_de")

    sig <- !is.na(padj) & padj < alpha
    message(sprintf("[de:%s] %d/%d significant at padj < %g (%d up, %d down)",
                    assay_name, sum(sig), nrow(out), alpha,
                    sum(sig & lfc > 0, na.rm = TRUE), sum(sig & lfc < 0, na.rm = TRUE)))
    message(sprintf("[de:%s] log2FC basis: %s", assay_name, LOG2FC_BASIS[[tag]]))

    saveRDS(dds, file.path(tool_dir(opts, "rna_ribo_de"),
                           sprintf("%s.%s_DE_dds.rds", d$name, tag)))
    list(res = out, dds = dds)
}

rna  <- run_de("RNA")
ribo <- run_de("Ribo")

# ---------------------------------------------------------------------------
# MASTER TABLE - every method, every statistic, every expression measure
# ---------------------------------------------------------------------------
# One row per gene, carrying: RNA DE, Ribo DE, TE from all three methods,
# RPKM/FPKM/TPM per assay x condition, and the regulatory-mode call. This is
# the table to open when comparing methods against each other.

# Each TE method writes into its own sub-directory, named after the tool.
te_methods <- list(deseq2  = file.path("deseq2",  sprintf("%s.deseq2_TE.tsv",  d$name)),
                   xtail   = file.path("xtail",   sprintf("%s.xtail_TE.tsv",   d$name)),
                   riborex = file.path("riborex", sprintf("%s.riborex_TE.tsv", d$name)))
if (!is.null(te_file)) te_methods <- list(custom = te_file)

master <- merge(rna$res, ribo$res, by = "gene_id", all = TRUE)

found <- character()
for (meth in names(te_methods)) {
    rel <- te_methods[[meth]]
    f <- if (file.exists(rel)) rel else file.path(opts$outdir, rel)
    if (!file.exists(f)) {
        message(sprintf("[master] %s TE table not found - skipping (run 0%s first)",
                        meth, switch(meth, deseq2 = "5_te_deseq2.R", xtail = "6_te_xtail.R",
                                     riborex = "7_te_riborex.R", "?")))
        next
    }
    te <- read.delim(f)
    # Keep only this method's own statistics; the expression columns are
    # identical across methods and are added once, below.
    cols <- c("gene_id", grep("_TE_", names(te), value = TRUE))
    master <- merge(master, te[, cols, drop = FALSE], by = "gene_id", all = TRUE)
    found <- c(found, meth)
}
message(sprintf("\n[master] TE methods merged: %s",
                if (length(found)) paste(found, collapse = ", ") else "none"))

# ---------------------------------------------------------------------------
# Regulatory-mode classification, driven by the primary TE method
# ---------------------------------------------------------------------------
primary <- if ("deseq2" %in% found) "deseq2" else if (length(found)) found[1] else NA

if (!is.na(primary)) {
    lfc_te  <- master[[sprintf("log2FC_TE_%s", primary)]]
    padj_te <- master[[sprintf("padj_TE_%s",  primary)]]

    sig <- function(p) !is.na(p) & p < alpha
    s_rna  <- sig(master$padj_RNA)
    s_ribo <- sig(master$padj_Ribo)
    s_te   <- sig(padj_te)
    same_dir <- !is.na(master$log2FC_RNA) & !is.na(master$log2FC_Ribo) &
                sign(master$log2FC_RNA) == sign(master$log2FC_Ribo)

    master$mode <- "Unchanged"
    master$mode[s_rna &  s_ribo & !s_te & same_dir] <- "Transcriptional"
    master$mode[!s_rna & s_ribo &  s_te]            <- "Translational"
    master$mode[s_rna & !s_ribo &  s_te]            <- "Buffered"
    master$mode[s_rna &  s_ribo &  s_te & same_dir] <- "Intensified"
    # TE moved but neither assay reached significance alone - real at low n,
    # but weaker evidence than the four named modes.
    master$mode[!s_rna & !s_ribo & s_te]            <- "TE-only"
    master$mode_basis <- sprintf("TE from %s, alpha = %g", primary, alpha)

    # How many methods call each gene significant - the practical confidence
    # score when three methods disagree as much as they do at n=2.
    padj_cols <- grep("^padj_TE_", names(master), value = TRUE)
    master$n_TE_methods_sig <- rowSums(
        sapply(padj_cols, function(cc) !is.na(master[[cc]]) & master[[cc]] < alpha))

    master <- master[order(-master$n_TE_methods_sig, padj_te, master$padj_Ribo), ]
} else {
    message("[master] no TE table available - mode classification skipped")
    master <- master[order(master$padj_Ribo), ]
}

write_flat(master, opts, sprintf("%s.MASTER_all_results.tsv", d$name), tool = "master")
message(sprintf("[master] %d genes x %d columns", nrow(master), ncol(master)))

if (!is.na(primary)) {
    message(sprintf("\n[master] regulatory mode at padj < %g (TE from %s):", alpha, primary))
    tb <- sort(table(master$mode), decreasing = TRUE)
    for (n in names(tb)) message(sprintf("    %-16s %d", n, tb[[n]]))
}

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
f <- open_pdf(opts, sprintf("%s.DE_plots.pdf", d$name), tool = "master")

# --- per-assay DE ---
for (a_name in c("RNA", "Ribo")) {
    lfc  <- master[[sprintf("log2FC_%s", a_name)]]
    padj <- master[[sprintf("padj_%s",   a_name)]]
    base <- master[[sprintf("baseMean_counts_%s", a_name)]]
    plot_volcano(lfc, padj, sprintf("%s-seq DE volcano - %s", a_name, d$name), alpha)
    plot_ma(base, lfc, padj, sprintf("%s-seq DE MA - %s", a_name, d$name), alpha,
            xlab = "mean normalised count")
    tpm_cols <- grep(sprintf("^TPM_%s_", a_name), names(master), value = TRUE)
    if (length(tpm_cols)) {
        plot_lfc_vs_expr(rowMeans(master[, tpm_cols, drop = FALSE], na.rm = TRUE),
                         lfc, padj, sprintf("%s log2FC vs expression", a_name), alpha,
                         xlab = sprintf("mean %s TPM", a_name))
    }
}

# --- the key riboseq figure: transcription vs translation ---
if (!is.na(primary)) {
    ok <- is.finite(master$log2FC_RNA) & is.finite(master$log2FC_Ribo)
    cols <- c(Unchanged = "grey80", Transcriptional = "steelblue",
              Translational = "firebrick", Buffered = "darkorange",
              Intensified = "purple", `TE-only` = "grey45")
    plot(master$log2FC_RNA[ok], master$log2FC_Ribo[ok],
         col = cols[master$mode[ok]], pch = 19, cex = 0.5,
         xlab = "RNA log2FC", ylab = "Ribo log2FC",
         main = sprintf("Transcription vs translation - %s", d$name))
    abline(0, 1, lty = 2, col = "grey40"); abline(h = 0, v = 0, lty = 3, col = "grey70")
    legend("topleft", legend = names(cols), col = cols, pch = 19, cex = 0.7, bty = "n")

    bp <- barplot(sort(table(master$mode), decreasing = TRUE), las = 2, cex.names = 0.75,
                  col = "steelblue", ylab = "genes", main = "Regulatory mode")
    tb2 <- sort(table(master$mode), decreasing = TRUE)
    text(bp, tb2, labels = tb2, pos = 3, cex = 0.75, xpd = NA)

    # --- method comparison ---
    lfc_cols <- grep("^log2FC_TE_", names(master), value = TRUE)
    if (length(lfc_cols) >= 2) {
        pairs(master[, lfc_cols], pch = 19, cex = 0.25, col = "#00000030",
              main = "log2FC in TE - method vs method")
        for (i in seq_len(length(lfc_cols) - 1)) {
            for (j in (i + 1):length(lfc_cols)) {
                ok2 <- is.finite(master[[lfc_cols[i]]]) & is.finite(master[[lfc_cols[j]]])
                if (sum(ok2) < 10) next
                r <- cor(master[[lfc_cols[i]]][ok2], master[[lfc_cols[j]]][ok2])
                plot(master[[lfc_cols[i]]][ok2], master[[lfc_cols[j]]][ok2],
                     pch = 19, cex = 0.35, col = "#00000040",
                     xlab = lfc_cols[i], ylab = lfc_cols[j],
                     main = sprintf("%s vs %s  (r = %.3f)",
                                    sub("log2FC_TE_", "", lfc_cols[i]),
                                    sub("log2FC_TE_", "", lfc_cols[j]), r))
                abline(0, 1, lty = 2, col = "firebrick")
            }
        }
    }
    if (length(padj_cols) >= 2) {
        barplot(sapply(padj_cols, function(cc) sum(!is.na(master[[cc]]) & master[[cc]] < alpha)),
                names.arg = sub("padj_TE_", "", padj_cols), col = "firebrick",
                ylab = "significant genes", main = sprintf("TE hits per method (padj < %g)", alpha))
        barplot(table(factor(master$n_TE_methods_sig, levels = 0:length(padj_cols))),
                col = "steelblue", xlab = "number of methods calling the gene significant",
                ylab = "genes", main = "Method agreement")
    }

    # --- heatmap of top TE genes ---
    top <- head(master$gene_id[!is.na(padj_te)][order(padj_te[!is.na(padj_te)])], 40)
    if (length(top) > 2) {
        norm <- compute_norm(d$counts, d$lengths)
        hm <- log2(norm$TPM[top, , drop = FALSE] + 1)
        hm <- hm[rowSums(is.finite(hm)) == ncol(hm), , drop = FALSE]
        if (nrow(hm) > 2) {
            heatmap(hm, scale = "row", margins = c(12, 7), cexRow = 0.5, cexCol = 0.8,
                    main = sprintf("Top %d TE genes (log2 TPM, row-scaled)", nrow(hm)))
        }
    }
}

invisible(dev.off())
message(sprintf("[out] wrote %s", basename(f)))
message(sprintf("\n[de] done - all results in %s", opts$outdir))
