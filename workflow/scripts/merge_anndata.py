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

# Save merged object
logger.info(f"Saving merged object to: {output_file}")
Path(output_file).parent.mkdir(parents=True, exist_ok=True)
merged.write_h5ad(output_file)

logger.info("Done merging captures")
