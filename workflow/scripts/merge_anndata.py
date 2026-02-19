#!/usr/bin/env python3
"""
workflow/scripts/merge_anndata.py
Merge individual capture AnnData objects into a single object.

Parallels the logic in merge_captures.R.
"""

import logging
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
h5ad_files = snakemake.input.h5ad_files
output_file = snakemake.output.merged

logger.info(f"Merging {len(h5ad_files)} capture objects...")
logger.info(f"Output file: {output_file}")

# Load all objects
objects = []
for i, h5ad_file in enumerate(h5ad_files):
    logger.info(f"Loading: {h5ad_file}")
    adata = ad.read_h5ad(h5ad_file)
    logger.info(f"  - Cells: {adata.n_obs}")
    objects.append(adata)

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
