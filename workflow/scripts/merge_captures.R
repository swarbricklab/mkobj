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

# Save merged object
message("Saving merged object to: ", output_file)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
qsave(merged_object, file = output_file)

message("Done merging captures")

sink(type = "message"); sink()
close(log_con)
