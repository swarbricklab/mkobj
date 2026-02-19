#!/usr/bin/env python3
"""
workflow/scripts/attach_metadata_anndata.py
Attach experimental metadata to merged AnnData object.

Parallels the logic in attach_metadata.R.
"""

import logging
import pandas as pd
import anndata as ad
from pathlib import Path

# Set up logging
logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Get input and output paths
merged_file = snakemake.input.merged
metadata_file = snakemake.input.metadata
output_file = snakemake.output.annotated
join_column = snakemake.params.join_column

logger.info("Attaching metadata to merged AnnData object")
logger.info(f"Input: {merged_file}")
logger.info(f"Metadata: {metadata_file}")
logger.info(f"Metadata join column: {join_column}")

# Load merged object
logger.info("Loading merged object...")
adata = ad.read_h5ad(merged_file)
logger.info(f"Loaded object with {adata.n_obs} cells")

# Load metadata
logger.info("Loading metadata...")
meta = pd.read_csv(metadata_file)
logger.info(f"Loaded metadata with {len(meta)} rows and {len(meta.columns)} columns")

# Rename the metadata join column to 'sample_id' to match the AnnData object
# The AnnData object always uses 'sample_id', but the metadata CSV may use a different column name
if join_column != "sample_id":
    if join_column not in meta.columns:
        raise ValueError(
            f"Join column '{join_column}' not found in metadata. "
            f"Available columns: {', '.join(meta.columns)}"
        )
    logger.info(f"Renaming metadata column '{join_column}' to 'sample_id' for join")
    meta = meta.rename(columns={join_column: 'sample_id'})

# Remove capture_id column if present (already stored in AnnData as 'capture')
# and deduplicate to get one row per sample_id
if 'capture_id' in meta.columns:
    logger.info("Removing 'capture_id' column (already in AnnData as 'capture')")
    meta = meta.drop(columns=['capture_id'])

n_before = len(meta)
meta = meta.drop_duplicates()
logger.info(f"Deduplicated metadata: {n_before} -> {len(meta)} rows")

# Set sample_id as index for joining
meta = meta.set_index('sample_id')

# Join metadata to cell metadata on 'sample_id'
logger.info("Joining metadata on 'sample_id' column...")
n_cells = adata.n_obs

# Get current obs and preserve index
current_obs = adata.obs.copy()
current_obs['_cell_barcode'] = current_obs.index

# Join metadata
for col in meta.columns:
    # Skip columns that already exist in obs
    if col in current_obs.columns:
        logger.info(f"  Skipping column '{col}' (already exists in obs)")
        continue
    # Map metadata values based on sample_id
    current_obs[col] = current_obs['sample_id'].map(meta[col])

# Restore index
current_obs = current_obs.set_index('_cell_barcode')
current_obs.index.name = None

# Verify join didn't create issues
if len(current_obs) != n_cells:
    raise ValueError(
        f"Join issue: {len(current_obs)} rows vs {n_cells} cells"
    )

# Check for successful join
meta_cols = [c for c in meta.columns if c in current_obs.columns]
if meta_cols:
    n_matched = current_obs[meta_cols[0]].notna().sum()
    logger.info(f"Matched metadata for {n_matched} / {len(current_obs)} cells")

# Update obs in AnnData
adata.obs = current_obs

# Save annotated object
logger.info(f"Saving annotated object to: {output_file}")
Path(output_file).parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(output_file)

logger.info("Done attaching metadata")
