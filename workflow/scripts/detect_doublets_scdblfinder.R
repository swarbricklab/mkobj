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
rownames(counts) <- make.unique(features$gene_id[gex])
colnames(counts) <- barcodes
message("Read ", ncol(counts), " cells x ", nrow(counts), " GEX features")

sce <- SingleCellExperiment(assays = list(counts = counts))

# ── Run scDblFinder (package defaults; #42/#43/#47 add knownDoublets, dbr, clusters)
set.seed(seed)
sce <- scDblFinder(sce)

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
write.csv(calls, out_csv, row.names = FALSE, quote = FALSE)
message("Wrote ", out_csv)
