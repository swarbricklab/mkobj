#! /usr/bin/env Rscript
# workflow/scripts/merge_captures.R
# Merge individual capture Seurat objects into a single object
# Optionally filters cells to a subset defined by samples.csv

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(readr)
    library(dplyr)
    library(qs)
})

# Set up logging
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con); sink(log_con, type = "message")

# Get input and output paths
rds_files <- snakemake@input[["rds_files"]]
output_file <- snakemake@output[["merged"]]
samples_file <- snakemake@params[["samples_file"]]

# QC filtering parameters
min_counts <- snakemake@params[["min_counts"]]
max_counts <- snakemake@params[["max_counts"]]
min_features <- snakemake@params[["min_features"]]
max_features <- snakemake@params[["max_features"]]
max_pct_mt <- snakemake@params[["max_pct_mt"]]
max_pct_ribo <- snakemake@params[["max_pct_ribo"]]
filter_doublets <- snakemake@params[["filter_doublets"]]
min_cells_per_gene <- snakemake@params[["min_cells_per_gene"]]

message("Merging ", length(rds_files), " capture objects...")
message("Output file: ", output_file)

# ── Load samples filter (if any) ────────────────────────────────────────────────
samples_filter <- NULL
if (!is.null(samples_file) && samples_file != "" && file.exists(samples_file)) {
    message("Loading samples filter from: ", samples_file)
    samples_filter <- read_csv(samples_file, show_col_types = FALSE)
    message("  Filtering to ", nrow(samples_filter), " sample-capture pairs")
}

# ── Helper: filter a Seurat object to a subset ──────────────────────────────────
filter_object <- function(obj, samples_filter) {
    if (is.null(samples_filter)) return(obj)

    # Determine capture from the first barcode prefix
    capture <- obj$capture[1]
    valid_samples <- samples_filter %>%
        filter(capture_id == capture) %>%
        pull(sample_id) %>%
        unique()

    if (length(valid_samples) == 0) {
        message("  No samples listed for capture ", capture, " — keeping all cells")
        return(obj)
    }

    n_before <- ncol(obj)
    cell_meta <- obj[[]]

    keep <- cell_meta$sample_id %in% valid_samples |
            (!is.na(cell_meta$status) & cell_meta$status == "doublet") |
            is.na(cell_meta$sample_id)

    obj <- obj[, keep]
    message("  Filtered ", capture, ": ", n_before, " -> ", ncol(obj), " cells",
            " (kept ", sum(cell_meta$sample_id %in% valid_samples), " cohort singlets, ",
            sum(!is.na(cell_meta$status) & cell_meta$status == "doublet"), " doublets, ",
            sum(is.na(cell_meta$sample_id)), " unassigned)")
    return(obj)
}

# ── Load and (optionally) filter objects ─────────────────────────────────────────
objects <- list()
for (i in seq_along(rds_files)) {
    message("Loading: ", rds_files[i])
    obj <- readRDS(rds_files[i])
    obj <- filter_object(obj, samples_filter)
    if (ncol(obj) > 0) {
        objects[[length(objects) + 1]] <- obj
    } else {
        message("  Skipping — no cells remaining after filter")
    }
}

if (length(objects) == 0) {
    stop("No cells remaining after filtering — cannot create merged object")
}

# Merge objects
if (length(objects) == 1) {
    message("Only one capture - no merging needed")
    merged_object <- objects[[1]]
} else {
    message("Merging captures...")
    merged_object <- merge(
        x = objects[[1]],
        y = objects[2:length(objects)]
    )
    
    # Join layers for RNA assay
    message("Joining RNA layers...")
    merged_object[["RNA"]] <- JoinLayers(merged_object[["RNA"]])
    
    # Join layers for AB assay if present
    if ("AB" %in% names(merged_object@assays)) {
        message("Joining AB layers...")
        merged_object[["AB"]] <- JoinLayers(merged_object[["AB"]])
    }
}

message("Merged object has ", ncol(merged_object), " cells total")

# ── QC cell filtering ────────────────────────────────────────────────────────────
n_before <- ncol(merged_object)
keep <- rep(TRUE, n_before)
meta <- merged_object[[]]

if (!is.null(filter_doublets) && isTRUE(filter_doublets) && "predicted_doublet" %in% colnames(meta)) {
    is_doublet <- !is.na(meta$predicted_doublet) & meta$predicted_doublet
    message("  Doublet filter: flagged ", sum(is_doublet), " cells")
    keep <- keep & !is_doublet
}

if (!is.null(min_counts) && "nCount_RNA" %in% colnames(meta)) {
    fail <- meta$nCount_RNA < min_counts
    message("  min_counts (", min_counts, "): flagged ", sum(fail, na.rm = TRUE), " cells")
    keep <- keep & !fail
}

if (!is.null(max_counts) && "nCount_RNA" %in% colnames(meta)) {
    fail <- meta$nCount_RNA > max_counts
    message("  max_counts (", max_counts, "): flagged ", sum(fail, na.rm = TRUE), " cells")
    keep <- keep & !fail
}

if (!is.null(min_features) && "nFeature_RNA" %in% colnames(meta)) {
    fail <- meta$nFeature_RNA < min_features
    message("  min_features (", min_features, "): flagged ", sum(fail, na.rm = TRUE), " cells")
    keep <- keep & !fail
}

if (!is.null(max_features) && "nFeature_RNA" %in% colnames(meta)) {
    fail <- meta$nFeature_RNA > max_features
    message("  max_features (", max_features, "): flagged ", sum(fail, na.rm = TRUE), " cells")
    keep <- keep & !fail
}

if (!is.null(max_pct_mt) && "percent.mt" %in% colnames(meta)) {
    fail <- meta$percent.mt > max_pct_mt
    message("  max_pct_mt (", max_pct_mt, "): flagged ", sum(fail, na.rm = TRUE), " cells")
    keep <- keep & !fail
}

if (!is.null(max_pct_ribo) && "percent.ribo" %in% colnames(meta)) {
    fail <- meta$percent.ribo > max_pct_ribo
    message("  max_pct_ribo (", max_pct_ribo, "): flagged ", sum(fail, na.rm = TRUE), " cells")
    keep <- keep & !fail
}

merged_object <- merged_object[, keep]
message("QC cell filtering: ", n_before, " -> ", ncol(merged_object), " cells (removed ", n_before - ncol(merged_object), ")")

# ── Gene QC: flag low-expression genes (do NOT drop) ─────────────────────────────
if (!is.null(min_cells_per_gene)) {
    counts_mat <- GetAssayData(merged_object, assay = "RNA", layer = "counts")
    cells_per_gene <- Matrix::rowSums(counts_mat > 0)
    merged_object[["RNA"]]@meta.data$n_cells <- cells_per_gene
    merged_object[["RNA"]]@meta.data$is_filtered <- cells_per_gene < min_cells_per_gene
    n_filtered_genes <- sum(cells_per_gene < min_cells_per_gene)
    message("Gene QC: flagged ", n_filtered_genes, "/", nrow(merged_object),
            " genes as is_filtered (expressed in < ", min_cells_per_gene, " cells)")
}

# Save merged object
message("Saving merged object to: ", output_file)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
qsave(merged_object, file = output_file)

message("Done merging captures")

sink(type = "message"); sink()
close(log_con)
