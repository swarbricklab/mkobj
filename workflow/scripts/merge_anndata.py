#!/usr/bin/env python3
"""
workflow/scripts/merge_anndata.py
Merge individual capture AnnData objects into a single object.
Optionally filters cells to a subset defined by samples.csv.

Parallels the logic in merge_captures.R.
"""

import sys
import logging
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path

# Set up logging — write to snakemake log file AND stderr
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(snakemake.log[0]),
        logging.StreamHandler(sys.stderr),
    ],
)
logger = logging.getLogger(__name__)

# Ensure uncaught exceptions are logged to the log file
def _excepthook(exc_type, exc_value, exc_tb):
    logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_tb))
    sys.__excepthook__(exc_type, exc_value, exc_tb)

sys.excepthook = _excepthook

# Get input and output paths
h5ad_files = snakemake.input.h5ad_files
output_file = snakemake.output.merged
samples_file = snakemake.params.get("samples_file", "")

# QC filtering parameters
min_counts = snakemake.params.get("min_counts", None)
max_counts = snakemake.params.get("max_counts", None)
min_features = snakemake.params.get("min_features", None)
max_features = snakemake.params.get("max_features", None)
max_pct_mt = snakemake.params.get("max_pct_mt", None)
max_pct_ribo = snakemake.params.get("max_pct_ribo", None)
filter_doublets = snakemake.params.get("filter_doublets", True)
min_cells_per_gene = snakemake.params.get("min_cells_per_gene", 3)

logger.info(f"Merging {len(h5ad_files)} capture objects...")
logger.info(f"Output file: {output_file}")


# ── Load samples filter (if any) ────────────────────────────────────────────────
samples_filter = None
if samples_file and Path(samples_file).exists():
    logger.info(f"Loading samples filter from: {samples_file}")
    samples_filter = pd.read_csv(samples_file)
    logger.info(f"  Filtering to {len(samples_filter)} sample-capture pairs")


def filter_adata(adata: ad.AnnData, samples_filter: pd.DataFrame | None) -> ad.AnnData:
    """Filter an AnnData object to cells belonging to the subset."""
    if samples_filter is None:
        return adata

    capture = adata.obs['capture'].iloc[0]
    valid_samples = samples_filter.loc[
        samples_filter['capture_id'] == capture, 'sample_id'
    ].unique().tolist()

    if not valid_samples:
        logger.info(f"  No samples listed for capture {capture} — keeping all cells")
        return adata

    n_before = adata.n_obs

    is_cohort = adata.obs['sample_id'].isin(valid_samples)
    is_doublet = adata.obs['status'].eq('doublet') & adata.obs['status'].notna()
    is_unassigned = adata.obs['sample_id'].isna()

    keep = is_cohort | is_doublet | is_unassigned
    adata = adata[keep].copy()

    logger.info(
        f"  Filtered {capture}: {n_before} -> {adata.n_obs} cells "
        f"(kept {is_cohort.sum()} cohort singlets, "
        f"{is_doublet.sum()} doublets, "
        f"{is_unassigned.sum()} unassigned)"
    )
    return adata


# ── Load and (optionally) filter objects ─────────────────────────────────────────
objects = []
for i, h5ad_file in enumerate(h5ad_files):
    logger.info(f"Loading: {h5ad_file}")
    adata = ad.read_h5ad(h5ad_file)
    adata = filter_adata(adata, samples_filter)
    if adata.n_obs > 0:
        objects.append(adata)
    else:
        logger.info("  Skipping — no cells remaining after filter")

if not objects:
    raise ValueError("No cells remaining after filtering — cannot create merged object")

# Merge objects
if len(objects) == 1:
    logger.info("Only one capture - no merging needed")
    merged = objects[0]
else:
    logger.info("Merging captures...")
    # Concatenate along obs axis (cells)
    merged = ad.concat(
        objects,
        join='outer',  # Keep all genes from all objects
        merge='same',  # Only merge identical elements
        label='batch',
        keys=[Path(f).stem for f in h5ad_files],
        index_unique=None  # Barcodes already prefixed with capture
    )
    
    # Handle AB modality if present in obsm
    # ad.concat handles obsm automatically, but we need to ensure alignment
    ab_present = all('AB' in obj.obsm for obj in objects)
    if ab_present:
        logger.info("AB modality present in all captures")
        # Check if all AB features are the same
        ab_features = [obj.uns.get('AB_features', []) for obj in objects]
        if all(f == ab_features[0] for f in ab_features):
            merged.uns['AB_features'] = ab_features[0]
        else:
            logger.warning("AB features differ across captures - storing only feature names")
            # Store unique union of AB features
            all_ab_features = []
            for f in ab_features:
                all_ab_features.extend(f)
            merged.uns['AB_features'] = list(dict.fromkeys(all_ab_features))

logger.info(f"Merged object has {merged.n_obs} cells total")

# ── QC cell filtering ────────────────────────────────────────────────────────────
n_before = merged.n_obs
keep = np.ones(merged.n_obs, dtype=bool)

if filter_doublets and 'predicted_doublet' in merged.obs.columns:
    is_doublet = merged.obs['predicted_doublet'].astype(bool)
    n_doublet = is_doublet.sum()
    keep &= ~is_doublet
    logger.info(f"  Doublet filter: flagged {n_doublet} cells")

if min_counts is not None and 'nCount_RNA' in merged.obs.columns:
    fail = merged.obs['nCount_RNA'] < min_counts
    logger.info(f"  min_counts ({min_counts}): flagged {fail.sum()} cells")
    keep &= ~fail

if max_counts is not None and 'nCount_RNA' in merged.obs.columns:
    fail = merged.obs['nCount_RNA'] > max_counts
    logger.info(f"  max_counts ({max_counts}): flagged {fail.sum()} cells")
    keep &= ~fail

if min_features is not None and 'nFeature_RNA' in merged.obs.columns:
    fail = merged.obs['nFeature_RNA'] < min_features
    logger.info(f"  min_features ({min_features}): flagged {fail.sum()} cells")
    keep &= ~fail

if max_features is not None and 'nFeature_RNA' in merged.obs.columns:
    fail = merged.obs['nFeature_RNA'] > max_features
    logger.info(f"  max_features ({max_features}): flagged {fail.sum()} cells")
    keep &= ~fail

if max_pct_mt is not None and 'percent.mt' in merged.obs.columns:
    fail = merged.obs['percent.mt'] > max_pct_mt
    logger.info(f"  max_pct_mt ({max_pct_mt}): flagged {fail.sum()} cells")
    keep &= ~fail

if max_pct_ribo is not None and 'percent.ribo' in merged.obs.columns:
    fail = merged.obs['percent.ribo'] > max_pct_ribo
    logger.info(f"  max_pct_ribo ({max_pct_ribo}): flagged {fail.sum()} cells")
    keep &= ~fail

merged = merged[keep].copy()
logger.info(f"QC cell filtering: {n_before} -> {merged.n_obs} cells (removed {n_before - merged.n_obs})")

# ── Gene QC: flag low-expression genes (do NOT drop) ────────────────────────────
if min_cells_per_gene is not None:
    from scipy.sparse import issparse
    if issparse(merged.X):
        cells_per_gene = np.asarray((merged.X > 0).sum(axis=0)).flatten()
    else:
        cells_per_gene = np.asarray((merged.X > 0).sum(axis=0)).flatten()
    merged.var['n_cells'] = cells_per_gene
    merged.var['is_filtered'] = cells_per_gene < min_cells_per_gene
    n_filtered_genes = merged.var['is_filtered'].sum()
    logger.info(f"Gene QC: flagged {n_filtered_genes}/{merged.n_vars} genes as is_filtered (expressed in < {min_cells_per_gene} cells)")

# Save merged object
logger.info(f"Saving merged object to: {output_file}")
Path(output_file).parent.mkdir(parents=True, exist_ok=True)
merged.write_h5ad(output_file)

logger.info("Done merging captures")
