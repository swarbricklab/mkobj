# Plan: Integrating Preprocessing Steps from `brca_mega_atlas` into `mkobj`

This document outlines the plan to align the `mkobj` AnnData/Seurat creation workflow with the more comprehensive preprocessing steps currently used in the `brca_mega_atlas` repository.

## 1. Comparison of Current Approaches

### Similarities
- Both workflows aim to take raw counts and generate individual and merged AnnData/Seurat objects for downstream analysis.
- Both aggregate multiple captures/samples.

### Incompatibilities and Divergences
- **Scope**: `mkobj` is strictly focused on object creation and merging (a lightweight structural pipeline). In contrast, `brca_mega_atlas/analysis/preprocessing` acts as a full preprocessing pipeline (handling filtering, doublet detection, normalization, dimensionality reduction, and preliminary annotation).
- **Filtering**: `mkobj` defers most filtering to downstream steps. `brca_mega_atlas` performs comprehensive filtering (e.g., visualising QC metrics to set thresholds for counts, features, mitochondrial, and ribosomal content).
- **Doublet Detection**: `mkobj` currently does not perform explicit doublet detection (it may rely on demultiplexing flags). `brca_mega_atlas` runs a robust dual-method doublet detection:
  1. Computational prediction using **Scrublet**.
  2. Biology-driven **Cross-lineage** detection (identifying cells expressing strong marker genes from multiple incompatible lineages).
- **Data Transformations**: `brca_mega_atlas` applies normalization, log-transformation, batch-aware highly variable gene (HVG) selection, and produces PCA/UMAP embeddings early in the process. `mkobj` keeps the data mostly raw.
- **Output States**: `brca_mega_atlas` maintains intermediate `.h5ad` files at various filtering stages and clearly separates the raw counts into a `.layers["counts"]` before transforming `adata.X`.

## 2. Plan for `mkobj` Integration

To incorporate the advanced preprocessing steps into `mkobj`, we will add new rules and scripts to the Snakemake workflow. This ensures `mkobj` produces high-quality, clean objects ready for downstream biological analysis.

### Step 2.1: Enhance Initial Object Creation (QC Metrics)
**Goal:** Calculate comprehensive QC metrics during the initial AnnData/Seurat object creation.
**Actions:**
- Update `create_anndata.py` and `create_seurat.R` to calculate both mitochondrial (`percent.mt`) and ribosomal (`percent.ribo`) content.
- Compute standard metrics (e.g., `nFeature`, `nCount`).
- Save these metrics in the objects' metadata/obs to facilitate downstream visualization and filtering.

### Step 2.2: Implement Doublet Detection
**Goal:** Add a dedicated doublet detection rule to identify and flag (or remove) doublets before merging.
**Actions:**
- Create a new script (e.g., `detect_doublets.py`).
- **Method A (Scrublet):** Implement Scrublet to score cells and flag them based on a defined threshold (e.g., > 0.2).
- **Method B (Cross-lineage):** Implement the cross-lineage marker detection logic. This will require an input configuration defining major lineage marker genes.
- Combine the results of both methods into a single `predicted_doublet` boolean column in `adata.obs`.
- Add a corresponding Snakemake rule (`doublet_detection.smk`) to execute this script on each capture prior to, or immediately following, the merge step.

### Step 2.3: Implement QC Filtering and Merging Logic
**Goal:** Filter out low-quality cells and true doublets, while preserving genes for CELLxGENE schema.
**Actions:**
- Update the merge scripts (`merge_anndata.py` / `merge_captures.R`) or create a dedicated filtering script to apply minimum/maximum thresholds for `nCount`, `nFeature`, `percent.mt`, and `percent.ribo` to remove poor quality cells.
- Remove cells flagged as `predicted_doublet`.
- **Gene Retention rule:** Do NOT drop genes. Instead of removing genes that do not meet minimum expression criteria, flag them by adding a boolean `is_filtered` column in `adata.var` (or equivalent in Seurat).
- Optionally, output QC summary plots at this stage to document the effect of filtering.

### Step 2.4: Standardize Data Transformations and Projection
**Goal:** Provide normalized data, embeddings, and clustering as part of the core output while preserving the full gene set.
**Actions:**
- Add a processing script (`process_merged_data.py`) that takes the merged, filtered object.
- Back up raw counts to `adata.layers["counts"]`.
- While processing, subset to highly variable genes for internal clustering/PCA/UMAP computations, but **re-attach** the resulting embeddings and cluster identifications back onto the full object containing all genes.
- Perform total count normalization and log1p transformation.
- Compute Highly Variable Genes (accounting for batch effects if standard).
- Run PCA and UMAP.
- This ensures the output of `mkobj` is immediately usable for visualization and clustering in downstream atlasing, without losing features required by CELLxGENE.

### Step 2.5: Format Metadata for CELLxGENE
**Goal:** Ensure the final dataset structure contains standard metadata columns required by the CELLxGENE schema.
**Actions:**
- Add a step to map existing metadata to the standardized columns expected by CELLxGENE (e.g. standardizing sample origins, assay types, organism annotations).
- Validate the object against basic schema assumptions before completing the pipeline.

## 3. Implementation Phasing

1. **Phase 1: QC calculation and Doublet Flagging.** Retrofit the existing `create_*` scripts to calculate all QC metrics. Add the Scrublet approach to a new rule. This provides immediate value without disrupting the pipeline structure.
2. **Phase 2: Filtering and Reporting.** Add explicit filtering rules, mark non-target genes in `.var` as `is_filtered` instead of dropping them, and generate pre/post-filter summary statistics.
3. **Phase 3: Cross-lineage & Transformations.** Integrate the cross-lineage doublet detection (requiring marker gene definitions) and finalize the standard transformations (Normalization, HVG, PCA) on the merged dataset, keeping all genes on the final object.
4. **Phase 4: CELLxGENE Standardization.** Format the resulting metadata to match the standardized keys for the CELLxGENE upload schema.