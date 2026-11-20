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

from capture_io import get_valid_samples, load_assignments, read_10x_mtx_multimodal

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


def get_sample_donor_map(samples_file: str) -> dict:
    """sample_id -> donor_id from samples.csv, used to map a cell's assignment to
    its prop_ambient donor column (which is keyed by donor_id)."""
    if not samples_file or not Path(samples_file).exists():
        return {}
    samples = pd.read_csv(samples_file, dtype=str)
    if 'sample_id' not in samples.columns or 'donor_id' not in samples.columns:
        return {}
    return dict(zip(samples['sample_id'], samples['donor_id']))


def _resolve_donor_column(sample_id: str, donor_cols: list, sample_to_donor: dict):
    """Map a demux assignment label to its prop_ambient donor column.

    prop_ambient columns are genotype donor_ids; the assignment carries the
    sample_id. Try the label directly, then its donor_id (samples.csv), then the
    bare id before a '-N' suffix (e.g. 4270-1 -> 4270). PBMC-6 -> PBMC_Donor_6
    only resolves via the samples.csv donor_id. Returns None if unmappable.
    """
    if sample_id in donor_cols:
        return sample_id
    donor = sample_to_donor.get(sample_id)
    if donor and donor in donor_cols:
        return donor
    bare = sample_id.split('-')[0]
    if bare in donor_cols:
        return bare
    return None


def load_ambient(capture: str, ambient_root: str,
                 assignments: pd.DataFrame, sample_to_donor: dict) -> pd.DataFrame:
    """Load per-cell ambient RNA metadata: ambient_other and logLik_ratio.

    logLik_ratio (vireo's per-cell singlet-vs-doublet log-likelihood ratio) comes
    from prop_ambient.tsv, which carries it for EVERY cell -- including the doublet/
    unassigned cells that ambient_summary.csv leaves blank.

    ambient_other (a cell's read-share from donors OTHER than the one it is assigned
    to) is recomputed as 1 - prop_ambient[cell][assigned_donor], taking the assigned
    donor from the demuxafy CONSENSUS assignment (status == 'singlet') rather than
    vireo's self-assignment. The per-donor proportions sum to 1, so this reproduces
    ambient_summary.csv exactly wherever vireo had assigned the cell, and additionally
    fills the value for the demux singlets vireo left ambiguous. It stays NaN for
    doublets/unassigned, which have no single consensus donor (nor a kept-cell role).

    Falls back to reading ambient_summary.csv directly when prop_ambient.tsv is absent
    (no per-donor proportions available).
    """
    if not ambient_root:
        logger.info("No ambient root configured.")
        return pd.DataFrame()

    prop_path = Path(ambient_root) / capture / "prop_ambient.tsv"
    if not prop_path.exists():
        summary_path = Path(ambient_root) / capture / "ambient_summary.csv"
        if not summary_path.exists():
            logger.info("No ambient profiles found.")
            return pd.DataFrame()
        logger.info(f"No prop_ambient.tsv for {capture}; falling back to {summary_path}")
        return pd.read_csv(summary_path).set_index('barcode')

    prop = pd.read_csv(prop_path, sep='\t').set_index('cell')
    donor_cols = [c for c in prop.columns if c != 'logLik_ratio']

    result = pd.DataFrame(index=prop.index)
    result['logLik_ratio'] = prop['logLik_ratio']

    ambient_other = pd.Series(np.nan, index=prop.index, dtype=float)
    asg = assignments.reindex(prop.index)
    if {'sample_id', 'status'}.issubset(asg.columns):
        is_singlet = asg['status'].eq('singlet')
        for sid in asg.loc[is_singlet, 'sample_id'].dropna().unique():
            col = _resolve_donor_column(sid, donor_cols, sample_to_donor)
            if col is None:
                logger.warning(f"ambient_other: no prop_ambient donor column for "
                               f"assignment {sid!r} in {capture}; leaving NaN.")
                continue
            mask = is_singlet & asg['sample_id'].eq(sid)
            ambient_other.loc[mask] = 1.0 - prop.loc[mask, col].astype(float)
    result['ambient_other'] = ambient_other

    logger.info(f"Loaded ambient profiles for {len(result)} cells "
                f"({int(ambient_other.notna().sum())} with ambient_other)")
    return result[['ambient_other', 'logLik_ratio']]


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

# Determine valid samples up front so load_assignments can use it as a fallback
# for single-sample captures with no demux output.
valid_samples = get_valid_samples(capture, samples_file)

# Load and attach metadata
assignments = load_assignments(capture, barcodes, assignment_root, valid_samples)
annotations = load_annotations(capture, annotation_root)
# ambient_other is resolved against each cell's consensus assignment, so load the
# assignments first and pass them (plus the sample_id -> donor_id map) in.
sample_to_donor = get_sample_donor_map(samples_file)
ambient = load_ambient(capture, ambient_root, assignments, sample_to_donor)

# Merge all metadata
for df in [assignments, annotations, ambient]:
    if not df.empty:
        # Reindex to match adata.obs_names
        df_aligned = df.reindex(adata.obs_names)
        for col in df_aligned.columns:
            adata.obs[col] = df_aligned[col].values

# No cohort subsetting here (#39). The per-capture object keeps EVERY called cell,
# so detect_doublets sees the whole capture: passenger cells (other projects'
# singlets) form doublets with cohort cells, and doublet simulation and expected
# rates refer to the full capture. Cohort subsetting (keep cohort singlets,
# doublets and unassigned cells) is applied by merge_anndata.py, after doublet
# detection, using the same samples sheet.
if valid_samples is not None:
    n_cohort = int(adata.obs['sample_id'].isin(valid_samples).sum())
    logger.info(f"Keeping all {adata.n_obs} called cells for {capture} "
                f"({n_cohort} cohort singlets; cohort subsetting happens at merge)")

# Add capture ID to cell barcodes (for uniqueness when merging)
prefixed_barcodes = [f"{capture}_{bc}" for bc in adata.obs_names]
adata.obs_names = prefixed_barcodes

# Add capture metadata
adata.obs['capture'] = capture

# Calculate QC metrics
logger.info("Computing QC metrics (nFeature, nCount, percent.mt, percent.ribo)...")
adata.var['mt'] = adata.var['gene_name'].str.upper().str.startswith("MT-")
adata.var['ribo'] = adata.var['gene_name'].str.upper().str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

# Keep the scanpy QC names only (total_counts / n_genes_by_counts / pct_counts_mt /
# pct_counts_ribo). The duplicate Seurat-named copies (nCount_RNA / nFeature_RNA /
# percent.mt / percent.ribo) were dropped -- identical values, and the publication
# repo recomputes QC anyway. Downstream reads the scanpy names.

# Convert string columns with NA values to proper string type for h5ad compatibility
# h5py can't handle mixed types or np.nan in string arrays
for col in adata.obs.columns:
    if adata.obs[col].dtype == 'object' or pd.api.types.is_string_dtype(adata.obs[col]):
        # Fill NA with "NA" string and convert to string type
        adata.obs[col] = adata.obs[col].fillna("NA").astype(str)

# Save AnnData object
logger.info(f"Saving AnnData object to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(output_h5ad)
logger.info(f"Done creating AnnData object for capture: {capture}")
