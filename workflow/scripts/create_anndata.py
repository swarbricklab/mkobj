#!/usr/bin/env python3
"""
workflow/scripts/create_anndata.py
Create an AnnData object from Cell Ranger output for a single capture.

Parallels the logic in create_seurat.R.
"""

import sys
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
from scipy.io import mmread
from scipy.sparse import csr_matrix

# Set up logging
logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Extract parameters from snakemake
matrix_dir = Path(snakemake.input.matrix_dir)
output_h5ad = snakemake.output.h5ad
assignment_root = snakemake.params.assignment_root or ""
annotation_root = snakemake.params.annotation_root or ""
ambient_root = snakemake.params.ambient_root or ""
samples_file = snakemake.params.samples_file or ""
modality = snakemake.params.modality or "auto"
capture = snakemake.wildcards.capture

logger.info(f"Creating AnnData object for capture: {capture}")
logger.info(f"Matrix directory: {matrix_dir}")


def get_valid_samples(capture: str, samples_file: str) -> list | None:
    """Get valid sample_ids for this capture from samples.csv."""
    if not samples_file:
        logger.info("No samples file configured. No sample filtering will be applied.")
        return None
    
    samples_path = Path(samples_file)
    if not samples_path.exists():
        logger.info(f"Samples file not found: {samples_file}. No sample filtering will be applied.")
        return None
    
    samples = pd.read_csv(samples_path)
    
    if 'sample_id' not in samples.columns or 'capture_id' not in samples.columns:
        logger.info("Samples file missing required columns (sample_id, capture_id). No filtering applied.")
        return None
    
    # Filter to this capture and get unique sample_ids
    valid_samples = samples[samples['capture_id'] == capture]['sample_id'].dropna().unique().tolist()
    
    if not valid_samples:
        logger.info(f"No samples found for capture {capture} in samples file. No filtering applied.")
        return None
    
    logger.info(f"Valid samples from samples.csv: {', '.join(valid_samples)}")
    return valid_samples


def load_assignments(capture: str, barcodes: list, assignment_root: str) -> pd.DataFrame:
    """Load sample assignments for cells."""
    if not assignment_root:
        logger.info("No assignment root configured. Using capture as sample ID.")
        return pd.DataFrame({
            'status': ['singlet'] * len(barcodes),
            'sample_id': [capture] * len(barcodes)
        }, index=barcodes)
    
    assignment_path = Path(assignment_root) / capture / "cell_assignment.tsv"
    logger.info(f"Looking for assignments at: {assignment_path}")
    
    if not assignment_path.exists():
        logger.info("No sample assignments found. Using capture as sample ID.")
        return pd.DataFrame({
            'status': ['singlet'] * len(barcodes),
            'sample_id': [capture] * len(barcodes)
        }, index=barcodes)
    
    assignments = pd.read_csv(assignment_path, sep='\t')
    assignments = assignments.rename(columns={'assignment': 'sample_id'})
    assignments = assignments.set_index('barcode')
    logger.info(f"Loaded assignments for {len(assignments)} cells")
    return assignments


def load_annotations(capture: str, annotation_root: str) -> pd.DataFrame:
    """Load cell type annotations."""
    if not annotation_root:
        logger.info("No annotation root configured.")
        return pd.DataFrame()
    
    annotation_path = Path(annotation_root) / capture / "cell_types.csv"
    logger.info(f"Looking for annotations at: {annotation_path}")
    
    if not annotation_path.exists():
        logger.info("No cell annotations found.")
        return pd.DataFrame()
    
    annotations = pd.read_csv(annotation_path)
    annotations = annotations.set_index('barcode')
    logger.info(f"Loaded annotations for {len(annotations)} cells")
    return annotations


def load_ambient(capture: str, ambient_root: str) -> pd.DataFrame:
    """Load ambient RNA profiles."""
    if not ambient_root:
        logger.info("No ambient root configured.")
        return pd.DataFrame()
    
    ambient_path = Path(ambient_root) / capture / "ambient_summary.csv"
    logger.info(f"Looking for ambient profiles at: {ambient_path}")
    
    if not ambient_path.exists():
        logger.info("No ambient profiles found.")
        return pd.DataFrame()
    
    ambient = pd.read_csv(ambient_path)
    ambient = ambient.set_index('barcode')
    logger.info(f"Loaded ambient profiles for {len(ambient)} cells")
    return ambient


def read_10x_mtx_multimodal(matrix_dir: Path) -> dict:
    """
    Read 10X Cell Ranger output, handling multimodal data.
    Returns a dict with modality names as keys.
    """
    features_file = matrix_dir / "features.tsv.gz"
    barcodes_file = matrix_dir / "barcodes.tsv.gz"
    matrix_file = matrix_dir / "matrix.mtx.gz"
    
    # Read features
    features = pd.read_csv(features_file, sep='\t', header=None,
                           names=['gene_id', 'gene_name', 'feature_type'])
    
    # Read barcodes
    barcodes = pd.read_csv(barcodes_file, sep='\t', header=None, names=['barcode'])
    
    # Read matrix
    matrix = mmread(matrix_file).T.tocsr()  # Transpose: cells x features
    
    # Check for multimodal data
    feature_types = features['feature_type'].unique()
    
    result = {}
    for ftype in feature_types:
        mask = features['feature_type'] == ftype
        ftype_features = features[mask]
        ftype_matrix = matrix[:, mask.values]
        
        # Create AnnData for this modality
        adata = ad.AnnData(X=ftype_matrix)
        adata.obs_names = barcodes['barcode'].values
        adata.var_names = ftype_features['gene_id'].values
        adata.var['gene_name'] = ftype_features['gene_name'].values
        adata.var['feature_type'] = ftype_features['feature_type'].values
        
        result[ftype] = adata
    
    return result


# Read Cell Ranger output
logger.info("Reading Cell Ranger output...")
modalities = read_10x_mtx_multimodal(matrix_dir)

available_modalities = list(modalities.keys())
logger.info(f"Available modalities: {', '.join(available_modalities)}")

# Handle multimodal data
if len(modalities) > 1:
    logger.info("Detected multimodal data")
    
    # Determine which modality to use for main object
    if modality == "auto" or modality == "Gene Expression":
        if "Gene Expression" not in modalities:
            raise ValueError("Gene Expression modality not found in multimodal data")
        adata = modalities["Gene Expression"]
    else:
        adata = modalities[modality]
    
    # Add antibody capture as additional modality if present
    if "Antibody Capture" in modalities:
        logger.info("Adding Antibody Capture as additional modality...")
        ab_data = modalities["Antibody Capture"]
        # Filter to common cells
        common_cells = adata.obs_names.intersection(ab_data.obs_names)
        ab_data = ab_data[common_cells, :]
        # Store in obsm as a dense matrix and feature names in uns
        adata.obsm['AB'] = ab_data[adata.obs_names, :].X.toarray()
        adata.uns['AB_features'] = ab_data.var_names.tolist()
        logger.info(f"Added AB assay with {ab_data.n_vars} features")
else:
    # Unimodal data
    logger.info("Detected unimodal data (gene expression only)")
    adata = list(modalities.values())[0]

logger.info(f"Created AnnData object with {adata.n_obs} cells and {adata.n_vars} features")

# Get barcodes for metadata attachment
barcodes = adata.obs_names.tolist()

# Load and attach metadata
assignments = load_assignments(capture, barcodes, assignment_root)
annotations = load_annotations(capture, annotation_root)
ambient = load_ambient(capture, ambient_root)

# Merge all metadata
for df in [assignments, annotations, ambient]:
    if not df.empty:
        # Reindex to match adata.obs_names
        df_aligned = df.reindex(adata.obs_names)
        for col in df_aligned.columns:
            adata.obs[col] = df_aligned[col].values

# Subset cells based on samples.csv if provided
valid_samples = get_valid_samples(capture, samples_file)

if valid_samples is not None:
    logger.info(f"\nFiltering cells for capture: {capture}")
    n_input = adata.n_obs
    
    # Determine which cells to keep:
    # 1. Cells with sample_id in valid_samples (cohort singlets)
    # 2. Cells with status == "doublet"
    # 3. Cells with NA sample_id (unassigned)
    is_cohort_singlet = adata.obs['sample_id'].isin(valid_samples)
    is_doublet = adata.obs['status'].fillna('') == 'doublet'
    is_unassigned = adata.obs['sample_id'].isna()
    
    keep_cells = is_cohort_singlet | is_doublet | is_unassigned
    
    # Log filtering stats
    n_cohort_singlets = is_cohort_singlet.sum()
    n_doublets = is_doublet.sum()
    n_unassigned = is_unassigned.sum()
    n_kept = keep_cells.sum()
    n_filtered = n_input - n_kept
    
    logger.info(f"  Valid samples from samples.csv: {', '.join(valid_samples)}")
    logger.info(f"  Input: {n_input} cells")
    logger.info(f"  Kept: {n_kept} cells")
    logger.info(f"    - Cohort singlets: {n_cohort_singlets}")
    logger.info(f"    - Doublets: {n_doublets}")
    logger.info(f"    - Unassigned: {n_unassigned}")
    logger.info(f"  Filtered out: {n_filtered} cells (samples not in cohort)")
    
    # Subset the AnnData object
    adata = adata[keep_cells, :].copy()

# Add capture ID to cell barcodes (for uniqueness when merging)
prefixed_barcodes = [f"{capture}_{bc}" for bc in adata.obs_names]
adata.obs_names = prefixed_barcodes

# Add capture metadata
adata.obs['capture'] = capture

# Save AnnData object
logger.info(f"Saving AnnData object to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(output_h5ad)
logger.info(f"Done creating AnnData object for capture: {capture}")
