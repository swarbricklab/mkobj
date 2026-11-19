#! /usr/bin/env Rscript
# workflow/scripts/detect_doublets_scdblfinder.R
# Run scDblFinder on every called cell of one capture and write long-format
# calls (engine = scdblfinder) in the doublets/<capture>.csv schema (#41).
#
# Reads the Cell Ranger filtered matrix directly (Gene Expression features only),
# i.e. the same cell set as per_capture_raw/<capture>.h5ad after #39.

# Open the log first so a failed library load is recorded there too.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con); sink(log_con, type = "message")

suppressPackageStartupMessages({
    library(Matrix)
    library(SingleCellExperiment)
    library(scDblFinder)
})

capture    <- snakemake@wildcards[["capture"]]
matrix_dir <- snakemake@input[["matrix_dir"]]
out_csv    <- snakemake@output[["calls"]]
seed       <- as.integer(snakemake@params[["seed"]])
known_csv  <- snakemake@input[["known"]]           # NULL when demux is not configured
rate       <- read.csv(snakemake@input[["params_csv"]], stringsAsFactors = FALSE)  # #43
clus_df    <- read.csv(snakemake@input[["clusters"]], colClasses = "character")    # #47
known_use  <- snakemake@params[["known_use"]]
exclude_genes    <- unlist(snakemake@params[["exclude_genes"]])
exclude_patterns <- unlist(snakemake@params[["exclude_patterns"]])

message("scDblFinder ", as.character(packageVersion("scDblFinder")),
        " on capture ", capture, " (", matrix_dir, ")")

# ── Read the Cell Ranger matrix (features x cells) ──────────────────────────────
features <- read.delim(file.path(matrix_dir, "features.tsv.gz"), header = FALSE,
                       col.names = c("gene_id", "gene_name", "feature_type"),
                       colClasses = "character")
barcodes <- read.delim(file.path(matrix_dir, "barcodes.tsv.gz"), header = FALSE,
                       colClasses = "character")[[1]]
counts <- as(readMM(file.path(matrix_dir, "matrix.mtx.gz")), "CsparseMatrix")
gex <- features$feature_type == "Gene Expression"
counts <- counts[gex, , drop = FALSE]
gene_names <- features$gene_name[gex]
rownames(counts) <- make.unique(features$gene_id[gex])
colnames(counts) <- barcodes
message("Read ", ncol(counts), " cells x ", nrow(counts), " GEX features")

# ── Known (genotype-detected) doublets (#42) ────────────────────────────────────
known <- NULL
if (!is.null(known_csv) && length(known_csv) && nzchar(known_csv)) {
    kd <- read.csv(known_csv, colClasses = c(barcode = "character"))
    known <- barcodes %in% kd$barcode[as.logical(kd$known_doublet)]
    message("Known doublets: ", sum(known), " (", kd$known_source[1], "), knownUse = ",
            known_use)
}

# With knownUse = "positive" the known doublets are training positives. They are
# mostly cross-donor, so drop sex-chromosome and HLA genes from the expression
# space; otherwise the classifier can learn "two donors" instead of "two cells".
# (scDblFinder's trainingFeatures argument selects classifier features such as
# cxds_score, not genes, so the genes are removed from the matrix instead.)
if (!is.null(known) && known_use == "positive") {
    drop <- gene_names %in% exclude_genes
    for (p in exclude_patterns) drop <- drop | grepl(p, gene_names)
    message("knownUse=positive: excluding ", sum(drop), " genes from features: ",
            paste(sort(unique(gene_names[drop])), collapse = ", "))
    counts <- counts[!drop, , drop = FALSE]
    gene_names <- gene_names[!drop]
}

sce <- SingleCellExperiment(assays = list(counts = counts))
rowData(sce)$gene_name <- gene_names

# ── Run scDblFinder ─────────────────────────────────────────────────────────────
args <- list(sce)
if (!is.null(known)) {
    args$knownDoublets <- known
    args$knownUse <- known_use
}
# Clusters for artificial doublets (#47): lineage labels, or none (random).
if (nrow(clus_df) > 0) {
    cl <- clus_df$cluster[match(colnames(sce), clus_df$barcode)]
    n_missing <- sum(is.na(cl))
    # Lineage labels cover every called cell (#45); many misses means the barcodes
    # don't line up (e.g. a capture prefix), which would silently give one cluster.
    if (n_missing > 0.05 * length(cl))
        stop(n_missing, "/", length(cl), " cells have no lineage label in ",
             snakemake@input[["clusters"]], "; check that barcodes match")
    cl[is.na(cl)] <- "unlabelled"
    args$clusters <- cl
    message("clusters=", rate$clusters, ": passing ", length(unique(cl)), " clusters (",
            n_missing, " cells without a label -> unlabelled)")
} else {
    message("clusters=none: random artificial doublets")
}
# Expected doublet rate (#43): see doublet_rate.py / scdblfinder_params.py.
args$dbr.per1k <- rate$dbr_per1k
if (isTRUE(as.logical(rate$pass_dbr)))    args$dbr    <- rate$dbr_used
if (isTRUE(as.logical(rate$pass_dbr_sd))) args$dbr.sd <- rate$dbr_sd_used
fallback <- as.character(rate$fallback_reason)
message("dbr_source=", rate$dbr_source, " -> ", rate$dbr_source_used,
        if (!is.na(fallback) && nzchar(fallback)) paste0(" (", fallback, ")") else "",
        "; dbr.per1k=", rate$dbr_per1k,
        if (!is.null(args[["dbr"]])) paste0(", dbr=", signif(args[["dbr"]], 3)) else "",
        if (!is.null(args[["dbr.sd"]])) paste0(", dbr.sd=", signif(args[["dbr.sd"]], 3)) else "")
set.seed(seed)
sce <- do.call(scDblFinder, args)
# Audit the features actually used: count the configured exclusions (names and
# patterns) still present. Expect 0 with knownUse = "positive".
still <- rowData(sce)$gene_name %in% exclude_genes
for (p in exclude_patterns) still <- still | grepl(p, rowData(sce)$gene_name)
message("Features used (", nrow(sce), "): excluded genes/patterns present = ", sum(still))

calls <- data.frame(
    capture_id    = capture,
    barcode       = colnames(sce),
    engine        = "scdblfinder",
    doublet_score = sce$scDblFinder.score,
    doublet_call  = sce$scDblFinder.class == "doublet",
    stringsAsFactors = FALSE
)
message("Called ", sum(calls$doublet_call), " / ", nrow(calls), " doublets (",
        sprintf("%.2f", 100 * mean(calls$doublet_call)), "%)")

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
# write.csv keeps only 15 significant digits, so scores would not round-trip to
# the float64 scDblFinder computed. 17 significant digits always do;
# combine_doublet_calls re-writes them in shortest round-trip form.
calls$doublet_score <- sprintf("%.17g", calls$doublet_score)
write.csv(calls, out_csv, row.names = FALSE, quote = FALSE)
message("Wrote ", out_csv)
