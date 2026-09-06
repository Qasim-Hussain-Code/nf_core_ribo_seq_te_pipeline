#  Shared input handling for the three TE scripts.
#
#  Every method needs the same thing: a Ribo count matrix and an RNA count
#  matrix with identical gene order and matched sample order, plus the
#  condition each column belongs to. Doing that once here keeps the three
#  methods genuinely comparable - any difference in their output is the
#  method, not the input handling.
#

#  These scripts only work inside the `te_analysis` conda env. Without it
#  Rscript resolves to the system R, which has none of these packages, and the
#  failure ("there is no package called 'DESeq2'") does not hint at the cause.
#  Say it plainly instead.
.te_require <- function(pkgs) {
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
        stop(sprintf(paste0(
            "missing R package(s): %s\n",
            "  Using R at: %s (%s)\n\n",
            "  The TE packages live in the `te_analysis` conda env. Activate it first:\n",
            "      source ~/miniconda3/etc/profile.d/conda.sh\n",
            "      conda activate te_analysis\n\n",
            "  If the env does not exist yet, build it with:\n",
            "      bash 04_install_te_tools.sh"),
            paste(missing, collapse = ", "),
            R.home("bin"), R.version.string), call. = FALSE)
    }
}
.te_require("DESeq2")

suppressPackageStartupMessages({
    library(DESeq2)
})

PROJECT_DIR <- getwd()

TE_DEFAULTS <- list(
    counts      = file.path(PROJECT_DIR, "results/quantification/salmon/salmon.merged.gene_counts.tsv"),
    lengths     = file.path(PROJECT_DIR, "results/quantification/salmon/salmon.merged.gene_lengths.tsv"),
    samplesheet = file.path(PROJECT_DIR, "samplesheet.csv"),
    contrasts   = file.path(PROJECT_DIR, "contrasts.csv"),
    # Single flat output folder for every script - TE, DE and matrices alike.
    outdir      = file.path(PROJECT_DIR, "results/differential_expression"),
    min_count   = 10   # genes must exceed this summed over all samples
)

#' Parse `--key value` pairs off the command line, falling back to TE_DEFAULTS.
te_args <- function() {
    a <- commandArgs(trailingOnly = TRUE)
    opts <- TE_DEFAULTS
    i <- 1
    while (i <= length(a)) {
        key <- sub("^--", "", a[i])
        if (key %in% names(opts) && i < length(a)) {
            opts[[key]] <- a[i + 1]
            i <- i + 2
        } else if (key == "design" && i < length(a)) {
            opts$design <- a[i + 1]
            i <- i + 2
        } else {
            i <- i + 1
        }
    }
    opts$min_count <- as.numeric(opts$min_count)
    opts
}

#' Load counts + metadata and return matched Ribo/RNA matrices.
#'
#' salmon emits fractional estimated counts; every method here models raw
#' integer counts, so they are rounded. Genes with almost no signal are
#' dropped up front - with n=2 they only destabilise dispersion estimates.
load_te_data <- function(opts) {
    stopifnot(file.exists(opts$counts), file.exists(opts$samplesheet))

    raw <- read.delim(opts$counts, check.names = FALSE)
    rownames(raw) <- raw$gene_id
    mat <- as.matrix(raw[, setdiff(colnames(raw), c("gene_id", "gene_name")), drop = FALSE])
    mode(mat) <- "numeric"
    mat <- round(mat)

    ss <- read.csv(opts$samplesheet, stringsAsFactors = FALSE)
    ss <- ss[ss$sample %in% colnames(mat), , drop = FALSE]
    if (!nrow(ss)) stop("no samplesheet rows match the count matrix columns")

    ribo <- ss[ss$type == "riboseq", , drop = FALSE]
    rna  <- ss[ss$type == "rnaseq",  , drop = FALSE]
    if (!nrow(ribo) || !nrow(rna)) stop("need both riboseq and rnaseq rows in the samplesheet")

    # Match RNA to Ribo through the `pair` column so column i of both matrices
    # is the same biological replicate. Ordering by name would silently
    # mis-pair if the two assays were named differently.
    ribo <- ribo[order(ribo$condition, ribo$pair), , drop = FALSE]
    rna  <- rna[match(ribo$pair, rna$pair), , drop = FALSE]
    if (anyNA(rna$sample)) {
        stop("every riboseq sample needs an rnaseq sample with the same `pair` value")
    }

    ribo_mat <- mat[, ribo$sample, drop = FALSE]
    rna_mat  <- mat[, rna$sample,  drop = FALSE]

    keep <- (rowSums(ribo_mat) + rowSums(rna_mat)) > opts$min_count
    ribo_mat <- ribo_mat[keep, , drop = FALSE]
    rna_mat  <- rna_mat[keep, , drop = FALSE]

    ctr <- read.csv(opts$contrasts, stringsAsFactors = FALSE)[1, ]
    cond <- factor(ribo$condition, levels = c(ctr$reference, ctr$target))

    message(sprintf("[te] %d genes kept of %d (sum > %g)",
                    nrow(ribo_mat), nrow(mat), opts$min_count))
    message(sprintf("[te] contrast: %s (target) vs %s (reference)",
                    ctr$target, ctr$reference))
    message(sprintf("[te] Ribo: %s", paste(ribo$sample, collapse = ", ")))
    message(sprintf("[te] RNA : %s", paste(rna$sample,  collapse = ", ")))

    list(ribo = ribo_mat, rna = rna_mat, condition = cond,
         pair = factor(ribo$pair), ribo_meta = ribo, rna_meta = rna,
         contrast = ctr,
         name = paste0(ctr$target, "_vs_", ctr$reference))
}

#' Write a results table sorted by significance.
write_te <- function(df, opts, method, name) {
    dir.create(file.path(opts$outdir, method), recursive = TRUE, showWarnings = FALSE)
    out <- file.path(opts$outdir, method, sprintf("%s.%s_TE.tsv", name, method))
    write.table(df, out, sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("[te] wrote %s", out))
    out
}

#' Consistent one-line summary so the three methods can be compared at a glance.
summarise_te <- function(padj, lfc, method, alpha = 0.05) {
    sig <- !is.na(padj) & padj < alpha
    message(sprintf(
        "[%s] %d/%d genes with TE padj < %.2f  (%d up, %d down)",
        method, sum(sig), length(padj), alpha,
        sum(sig & lfc > 0, na.rm = TRUE), sum(sig & lfc < 0, na.rm = TRUE)))
}

#' Load the full 8-sample count matrix plus gene lengths and metadata.
#'
#' Unlike load_te_data(), this keeps every library as its own column - what
#' expression matrices and within-assay differential expression need. Rows are
#' filtered on the same threshold so the numbers stay comparable to the TE
#' tables.
load_expression_data <- function(opts) {
    stopifnot(file.exists(opts$counts), file.exists(opts$samplesheet))

    read_matrix <- function(path) {
        x <- read.delim(path, check.names = FALSE)
        rownames(x) <- x$gene_id
        m <- as.matrix(x[, setdiff(colnames(x), c("gene_id", "gene_name")), drop = FALSE])
        mode(m) <- "numeric"
        m
    }

    counts <- read_matrix(opts$counts)

    ss <- read.csv(opts$samplesheet, stringsAsFactors = FALSE)
    ss <- ss[ss$sample %in% colnames(counts), , drop = FALSE]
    # assay/condition first so groups sit together in plots and heatmaps
    ss <- ss[order(ss$type, ss$condition, ss$replicate), , drop = FALSE]
    counts <- counts[, ss$sample, drop = FALSE]

    lengths <- NULL
    if (!is.null(opts$lengths) && file.exists(opts$lengths)) {
        lengths <- read_matrix(opts$lengths)[rownames(counts), ss$sample, drop = FALSE]
    }

    keep <- rowSums(counts) > opts$min_count
    counts <- counts[keep, , drop = FALSE]
    if (!is.null(lengths)) lengths <- lengths[keep, , drop = FALSE]

    ctr <- read.csv(opts$contrasts, stringsAsFactors = FALSE)[1, ]

    meta <- data.frame(
        sample    = ss$sample,
        assay     = factor(ifelse(ss$type == "riboseq", "Ribo", "RNA"),
                           levels = c("RNA", "Ribo")),
        condition = factor(ss$condition, levels = c(ctr$reference, ctr$target)),
        replicate = factor(ss$replicate),
        pair      = factor(ss$pair),
        row.names = ss$sample,
        stringsAsFactors = FALSE
    )

    message(sprintf("[expr] %d genes x %d samples (sum > %g)",
                    nrow(counts), ncol(counts), opts$min_count))

    list(counts = counts, lengths = lengths, meta = meta, contrast = ctr,
         name = paste0(ctr$target, "_vs_", ctr$reference))
}

#' Write a matrix with a leading gene_id column.
write_matrix <- function(m, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    df <- data.frame(gene_id = rownames(m), m, check.names = FALSE)
    write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("[expr] wrote %s", basename(path)))
}

# ===========================================================================
#  Normalised expression, provenance, and plotting - shared by 05-09
# ===========================================================================

#  IMPORTANT - what log2FC is actually computed from.
#
#  Every log2FC and p-value in these outputs is derived from RAW COUNTS, using
#  each tool's own internal normalisation (DESeq2 median-of-ratios for DESeq2
#  and riborex; xtail's own size factors). RPKM/FPKM/TPM are NEVER used for
#  testing, and must not be: they are already normalised, which destroys the
#  count-level mean-variance relationship the negative-binomial model needs.
#  Handing TPM to DESeq2 yields confident-looking and wrong p-values.
#
#  RPKM/FPKM/TPM are attached to every table as descriptive expression levels
#  so a fold change can be read next to the abundance it came from. The
#  `log2FC_basis` column on every table records this explicitly.
LOG2FC_BASIS <- c(
    deseq2  = "raw_counts / DESeq2 median-of-ratios",
    riborex = "raw_counts / DESeq2 median-of-ratios (via riborex)",
    xtail   = "raw_counts / xtail internal size factors",
    rna     = "raw_counts / DESeq2 median-of-ratios",
    ribo    = "raw_counts / DESeq2 median-of-ratios"
)

#' CPM / RPKM / FPKM / TPM from a count matrix and matching length matrix.
#'
#' RPKM and FPKM are identical for single-end libraries - a "fragment" is a
#' "read" when reads are not paired. Both are returned so tables can carry the
#' name a reader expects, but they are the same numbers.
compute_norm <- function(counts, lengths) {
    if (is.null(lengths)) stop("gene lengths required for RPKM/FPKM/TPM", call. = FALSE)
    len <- lengths
    len[len <= 0] <- NA

    cpm  <- t(t(counts) / colSums(counts)) * 1e6
    rpkm <- cpm / (len / 1000)
    rate <- counts / (len / 1000)
    tpm  <- t(t(rate) / colSums(rate, na.rm = TRUE)) * 1e6

    list(CPM = cpm, RPKM = rpkm, FPKM = rpkm, TPM = tpm)
}

#' Mean expression per assay x condition group, as a gene_id-keyed data.frame.
#'
#' `assays` restricts which assays appear (e.g. only "RNA" for the RNA DE
#' table). Column names read like TPM_Ribo_NitrogenPlus.
expr_annotation <- function(counts, lengths, meta, assays = c("RNA", "Ribo"),
                            measures = c("TPM", "RPKM", "FPKM")) {
    norm <- compute_norm(counts, lengths)
    out <- data.frame(gene_id = rownames(counts), stringsAsFactors = FALSE)

    for (msr in measures) {
        m <- norm[[msr]]
        for (a in assays) {
            for (cond in levels(meta$condition)) {
                sel <- meta$assay == a & meta$condition == cond
                if (!any(sel)) next
                out[[sprintf("%s_%s_%s", msr, a, cond)]] <-
                    round(rowMeans(m[, meta$sample[sel], drop = FALSE], na.rm = TRUE), 3)
            }
        }
    }
    out
}

#' Attach expression columns and the log2FC provenance to a results table.
annotate_results <- function(df, ann, method) {
    df$log2FC_basis <- unname(LOG2FC_BASIS[method])
    merge(df, ann, by = "gene_id", all.x = TRUE, sort = FALSE)
}

# ---------------------------------------------------------------------------
#  Plots
# ---------------------------------------------------------------------------

#' Volcano: effect size against significance.
plot_volcano <- function(lfc, padj, title, alpha = 0.05, xlab = "log2 fold change") {
    ok <- is.finite(lfc) & !is.na(padj)
    if (!any(ok)) { plot.new(); title(paste(title, "- no data")); return(invisible()) }
    sig <- padj[ok] < alpha
    plot(lfc[ok], -log10(padj[ok]), pch = 19, cex = 0.4,
         col = ifelse(sig, ifelse(lfc[ok] > 0, "firebrick", "steelblue"), "grey75"),
         xlab = xlab, ylab = "-log10 adjusted p", main = title)
    abline(h = -log10(alpha), lty = 2, col = "grey40")
    abline(v = 0, lty = 3, col = "grey70")
    legend("topleft", bty = "n", cex = 0.75, pch = 19,
           col = c("firebrick", "steelblue", "grey75"),
           legend = c(sprintf("up (%d)",   sum(sig & lfc[ok] > 0)),
                      sprintf("down (%d)", sum(sig & lfc[ok] < 0)),
                      sprintf("ns (%d)",   sum(!sig))))
}

#' MA: effect size against expression level.
plot_ma <- function(base, lfc, padj, title, alpha = 0.05,
                    xlab = "mean expression") {
    ok <- is.finite(lfc) & is.finite(base) & base > 0
    if (!any(ok)) { plot.new(); title(paste(title, "- no data")); return(invisible()) }
    sig <- !is.na(padj[ok]) & padj[ok] < alpha
    plot(base[ok], lfc[ok], log = "x", pch = 19, cex = 0.4,
         col = ifelse(sig, "firebrick", "grey75"),
         xlab = xlab, ylab = "log2 fold change", main = title)
    abline(h = 0, lty = 2, col = "grey40")
}

#' log2FC against a normalised expression level, to expose fold changes that
#' rest on almost no signal.
plot_lfc_vs_expr <- function(expr, lfc, padj, title, alpha = 0.05,
                             xlab = "mean TPM") {
    ok <- is.finite(lfc) & is.finite(expr)
    if (!any(ok)) { plot.new(); title(paste(title, "- no data")); return(invisible()) }
    sig <- !is.na(padj[ok]) & padj[ok] < alpha
    plot(pmax(expr[ok], 1e-2), lfc[ok], log = "x", pch = 19, cex = 0.4,
         col = ifelse(sig, "firebrick", "grey75"),
         xlab = xlab, ylab = "log2 fold change", main = title)
    abline(h = 0, lty = 2, col = "grey40")
}

#  Output layout: one top-level folder, one sub-directory per tool, so a file's
#  provenance is visible from its path alone.
#
#      results/differential_expression/
#        deseq2/       05 - TE by DESeq2 interaction model
#        xtail/        06 - TE by Xtail
#        riborex/      07 - TE by riborex
#        expression/   08 - CPM / RPKM / FPKM / TPM / VST matrices
#        rna_ribo_de/  09 - DE within RNA-seq and within Ribo-seq
#        master/       09 - merged master table + method comparison

#' Path to a tool's sub-directory, created on demand.
tool_dir <- function(opts, tool) {
    d <- if (is.null(tool) || !nzchar(tool)) opts$outdir else file.path(opts$outdir, tool)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    d
}

#' Open a PDF inside a tool's sub-directory.
open_pdf <- function(opts, name, tool = NULL, width = 8, height = 7) {
    f <- file.path(tool_dir(opts, tool), name)
    pdf(f, width = width, height = height)
    f
}

#' Write a results table into a tool's sub-directory.
write_flat <- function(df, opts, filename, tool = NULL) {
    d <- tool_dir(opts, tool)
    f <- file.path(d, filename)
    write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("[out] wrote %s/%s",
                    if (is.null(tool)) basename(opts$outdir) else tool, filename))
    f
}
