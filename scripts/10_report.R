#!/usr/bin/env Rscript
#
#  10_report.R - self-contained HTML report from whatever results exist.
#
#  Deliberately defensive: every section checks for its inputs and reports what
#  is missing rather than failing. A run that stopped after quantification
#  still produces a readable report saying so.
#
#  Usage:
#      Rscript 10_report.R
#      Rscript 10_report.R --outdir <results dir> --alpha 0.1
#

source("te_common.R")
opts <- te_args()
a <- commandArgs(trailingOnly = TRUE)
alpha <- if ("--alpha" %in% a) as.numeric(a[which(a == "--alpha") + 1]) else 0.05

de_dir <- opts$outdir
rpt <- file.path(de_dir, "riboseq_report.html")
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

esc <- function(x) { x <- as.character(x)
    x <- gsub("&","&amp;",x,fixed=TRUE); x <- gsub("<","&lt;",x,fixed=TRUE)
    gsub(">","&gt;",x,fixed=TRUE) }

H <- c()
add <- function(...) H <<- c(H, ...)

add('<!doctype html><meta charset="utf-8"><title>Ribo-seq report</title>',
    '<style>
     body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
          max-width:1100px;margin:2rem auto;padding:0 1.5rem;line-height:1.6;color:#1a1a1a}
     h1{border-bottom:3px solid #2c5aa0;padding-bottom:.3rem}
     h2{margin-top:2.5rem;border-bottom:1px solid #ddd;padding-bottom:.2rem}
     table{border-collapse:collapse;width:100%;margin:1rem 0;font-size:.9rem}
     th{background:#2c5aa0;color:#fff;text-align:left;padding:.5rem}
     td{border-bottom:1px solid #eee;padding:.4rem .5rem}
     tr:nth-child(even){background:#fafafa}
     .miss{background:#fff4e5;border-left:4px solid #e08600;padding:.8rem 1rem;margin:1rem 0}
     .ok{background:#eaf6ea;border-left:4px solid #2e7d32;padding:.8rem 1rem;margin:1rem 0}
     code{background:#f4f4f4;padding:.1rem .3rem;border-radius:3px}
     details{margin:.5rem 0;background:#fafafa;border:1px solid #e0e0e0;border-radius:4px;padding:.5rem .8rem}
     summary{cursor:pointer;font-weight:600}
     .num{text-align:right;font-variant-numeric:tabular-nums}
     .foot{margin-top:3rem;padding-top:1rem;border-top:1px solid #ddd;
           font-size:.87rem;color:#555;text-align:center;line-height:1.9}
     .foot a{color:#0f5257;text-decoration:none} .foot a:hover{text-decoration:underline}
     .foot .nm{font-weight:600;color:#0f5257;font-size:.95rem}
     </style>')
add(sprintf(paste0("<h1>Ribo-seq analysis report</h1>",
                   "<p><em>generated %s</em></p>"),
            format(Sys.time(), "%Y-%m-%d %H:%M")))

note_missing <- function(what, how)
    add(sprintf('<div class="miss"><strong>%s not available.</strong> %s</div>', esc(what), esc(how)))

tbl <- function(df, n = 25) {
    df <- head(df, n)
    num <- vapply(df, is.numeric, logical(1))
    df[num] <- lapply(df[num], function(x) signif(x, 4))
    add("<table><tr>", paste0("<th>", esc(names(df)), "</th>", collapse = ""), "</tr>")
    for (i in seq_len(nrow(df))) {
        cells <- vapply(seq_along(df), function(j)
            sprintf('<td%s>%s</td>', if (num[j]) ' class="num"' else "", esc(df[i, j])), character(1))
        add("<tr>", paste0(cells, collapse = ""), "</tr>")
    }
    add("</table>")
}

# ---------------------------------------------------------------- run summary
add("<h2>1. Run configuration</h2>")
conf <- file.path(dirname(dirname(de_dir)), "riboseq_nextflow_guide", "project.conf")
if (!file.exists(conf)) conf <- Sys.getenv("PROJECT_CONF", "")
if (nzchar(conf) && file.exists(conf)) {
    kv <- read.dcf(textConnection(gsub("^([A-Z_]+)=", "\\1: ",
          grep("^[A-Z_]+=", readLines(conf), value = TRUE))))
    d <- data.frame(setting = colnames(kv), value = as.character(kv[1, ]))
    d <- d[d$setting %in% c("THREADS","RAM_GB","MAX_MEM_GB","HIGH_CPUS","HIGH_MEM_GB",
                            "GENOME_ACC","PIPELINE","PIPELINE_REV","GTF_PAD"), ]
    tbl(d, 20)
} else {
    note_missing("project.conf", "Run 00_configure.sh to record the hardware profile.")
}

# ------------------------------------------------------------- sample overview
add("<h2>2. Samples</h2>")
if (file.exists(opts$samplesheet)) {
    ss <- read.csv(opts$samplesheet, stringsAsFactors = FALSE)
    keep <- intersect(c("sample","type","condition","replicate","pair","srr"), names(ss))
    tbl(ss[, keep, drop = FALSE], 50)
    add(sprintf("<p>%d libraries: %s</p>", nrow(ss),
                paste(sprintf("%d %s", table(ss$type), names(table(ss$type))), collapse = ", ")))
} else note_missing("samplesheet.csv", "Expected at the project root.")

# ------------------------------------------------------------------ expression
add("<h2>3. Expression matrices</h2>")
expr_dir <- file.path(de_dir, "expression")
mats <- if (dir.exists(expr_dir)) list.files(expr_dir, pattern = "\\.tsv$") else character()
if (length(mats)) {
    add("<div class='ok'>Normalised matrices written. <strong>Use TPM to compare samples</strong> - ",
        "TPM columns sum to 1e6, RPKM columns do not, so RPKM is not comparable between libraries. ",
        "RPKM and FPKM are identical for single-end data.</div>")
    add("<ul>", paste0("<li><code>", esc(mats), "</code></li>", collapse = ""), "</ul>")
} else note_missing("Expression matrices", "Run 08_expression_matrices.R.")

# ----------------------------------------------------------------- DE + master
add("<h2>4. Differential expression and translational efficiency</h2>")
master_f <- list.files(file.path(de_dir, "master"), pattern = "MASTER_all_results\\.tsv$", full.names = TRUE)
if (length(master_f)) {
    m <- read.delim(master_f[1])
    add(sprintf("<p>Master table: <code>%s</code> - %d genes x %d columns.</p>",
                esc(basename(master_f[1])), nrow(m), ncol(m)))

    add("<div class='ok'><strong>Where log2FC comes from.</strong> Every fold change and p-value ",
        "is computed from <em>raw counts</em> using each tool's own normalisation ",
        "(DESeq2 median-of-ratios; xtail size factors). RPKM/FPKM/TPM are never used for testing - ",
        "they are already normalised, which breaks the count-level mean-variance relationship the ",
        "negative-binomial model needs. They appear as descriptive expression levels only.</div>")

    padj_cols <- grep("^padj_TE_", names(m), value = TRUE)
    if (length(padj_cols)) {
        counts <- vapply(padj_cols, function(cc) sum(!is.na(m[[cc]]) & m[[cc]] < alpha), integer(1))
        tbl(data.frame(method = sub("padj_TE_", "", padj_cols),
                       significant = as.integer(counts),
                       tested = vapply(padj_cols, function(cc) sum(!is.na(m[[cc]])), integer(1))), 10)
    }
    for (an in c("RNA","Ribo")) {
        pc <- sprintf("padj_%s", an)
        if (pc %in% names(m))
            add(sprintf("<p><strong>%s-seq DE:</strong> %d of %d genes at padj &lt; %g</p>",
                        an, sum(!is.na(m[[pc]]) & m[[pc]] < alpha), sum(!is.na(m[[pc]])), alpha))
    }
    if ("mode" %in% names(m)) {
        add("<h3>Regulatory modes</h3>")
        add("<p>What each gene's behaviour implies about where control is exerted:</p>")
        tb <- as.data.frame(sort(table(m$mode), decreasing = TRUE))
        names(tb) <- c("mode", "genes"); tbl(tb, 10)
        add("<details><summary>How these are defined</summary><ul>",
            "<li><strong>Transcriptional</strong> - RNA and Ribo both move, TE unchanged: output follows transcription.</li>",
            "<li><strong>Translational</strong> - Ribo moves, RNA does not: translation acts alone.</li>",
            "<li><strong>Buffered</strong> - RNA moves, Ribo does not: translation absorbs the mRNA change.</li>",
            "<li><strong>Intensified</strong> - both move the same way and TE amplifies it.</li>",
            "<li><strong>TE-only</strong> - TE moved but neither assay reached significance alone; weaker evidence.</li>",
            "</ul></details>")
    }
    if ("n_TE_methods_sig" %in% names(m)) {
        top <- m[order(-m$n_TE_methods_sig, m[[padj_cols[1]]]), ]
        show <- intersect(c("gene_id","n_TE_methods_sig","mode","log2FC_RNA","log2FC_Ribo",
                            grep("^log2FC_TE_", names(m), value = TRUE)), names(top))
        add("<h3>Top genes by method agreement</h3>",
            "<p>Sorted by how many TE methods call the gene significant - the practical ",
            "confidence score when methods disagree at low replication.</p>")
        tbl(top[, show, drop = FALSE], 25)
    }
} else note_missing("Master results table", "Run 05-07 then 09_de_rna_ribo.R.")

# ---------------------------------------------------------------------- plots
add("<h2>5. Figures</h2>")
pdfs <- list.files(de_dir, pattern = "\\.pdf$", recursive = TRUE)
if (length(pdfs)) {
    add("<ul>", paste0("<li><code>", esc(pdfs), "</code></li>", collapse = ""), "</ul>")
} else note_missing("Figures", "Each analysis script writes a PDF beside its table.")

# ------------------------------------------------------------------- provenance
add("<h2>6. Provenance</h2><details><summary>Session info</summary><pre>",
    esc(paste(capture.output(sessionInfo()), collapse = "\n")), "</pre></details>")

add('<div class="foot">',
    '<span class="nm">Ribo-seq translational efficiency analysis</span>',
    '</div>')

writeLines(H, rpt)
message(sprintf("[report] wrote %s", rpt))
