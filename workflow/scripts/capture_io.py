"""
workflow/scripts/capture_io.py
Per-capture readers shared by create_anndata.py and the annotation-stage rules
(detect_doublets, annotate_lineages, extract_known_doublets).

The annotation stage reads its inputs from here (Cell Ranger matrix, demux
files) rather than from the temp per_capture_raw/<capture>.h5ad, so rerunning
create_anndata for the merge does not invalidate the persisted annotation CSVs
(#49). Logic moved verbatim from create_anndata.py.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread

logger = logging.getLogger(__name__)


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


def load_assignments(capture: str, barcodes: list, assignment_root: str,
                     valid_samples: list | None = None) -> pd.DataFrame:
    """Load sample assignments for cells.

    When no demux assignment file is available, fall back to the sole sample
    listed for this capture in samples.csv (single-sample capture); otherwise
    fall back to the capture name.
    """
    if valid_samples is not None and len(valid_samples) == 1:
        fallback_sample = valid_samples[0]
    else:
        fallback_sample = capture

    if not assignment_root:
        logger.info(f"No assignment root configured. Using {fallback_sample} as sample ID.")
        return pd.DataFrame({
            'status': ['singlet'] * len(barcodes),
            'sample_id': [fallback_sample] * len(barcodes)
        }, index=barcodes)
    
    assignment_path = Path(assignment_root) / capture / "cell_assignment.tsv"
    logger.info(f"Looking for assignments at: {assignment_path}")

    if not assignment_path.exists():
        # demux IS configured but this capture has no assignment file. That is only
        # legitimate for a declared single-sample capture; otherwise refuse to
        # silently label every cell a singlet named after the capture — that is how a
        # failed run with no demux (e.g. Neomet_Pool029a) slipped into the merge as
        # named singletons.
        if valid_samples is not None and len(valid_samples) == 1:
            logger.info(f"No demux file; single-sample capture -> using {fallback_sample}.")
            return pd.DataFrame({
                'status': ['singlet'] * len(barcodes),
                'sample_id': [fallback_sample] * len(barcodes)
            }, index=barcodes)
        raise FileNotFoundError(
            f"Capture {capture!r}: demux root is configured but {assignment_path} is "
            f"missing, and it is not a declared single-sample capture "
            f"(valid_samples={valid_samples}). Refusing to treat all cells as singlets "
            f"named after the capture; register its demux output or its single sample."
        )
    
    assignments = pd.read_csv(assignment_path, sep='\t', dtype={'assignment': str, 'status': str})
    assignments = assignments.rename(columns={'assignment': 'sample_id'})
    assignments = assignments.set_index('barcode')

    # Resolve the non-singlet status into 'doublet' vs 'unassigned'. cell_assignment.tsv
    # records only 'singlet' or 'NA', and its singlet call is the demuxafy
    # AtLeastHalfSinglet consensus -- so the same consensus resolves the rest: a
    # non-singlet the consensus calls a doublet is a 'doublet'; otherwise it is a
    # singlet droplet with no confident donor -> 'unassigned'. (Note a SAME-DONOR
    # doublet is indistinguishable from a singlet to SNP demux, so 'singlet' never
    # fully excludes a doublet -- see docs/run.md.)
    non_singlet = assignments['status'] != 'singlet'
    sb_path = Path(assignment_root) / capture / "sample_barcode_assignments.tsv"
    if sb_path.exists() and non_singlet.any():
        consensus = pd.read_csv(
            sb_path, sep='\t',
            usecols=['Barcode', 'AtLeastHalfSinglet_DropletType'], dtype=str
        ).set_index('Barcode')['AtLeastHalfSinglet_DropletType'].reindex(assignments.index)
        assignments.loc[non_singlet, 'status'] = 'unassigned'
        assignments.loc[non_singlet & consensus.eq('doublet'), 'status'] = 'doublet'
    else:
        assignments.loc[non_singlet, 'status'] = 'unassigned'

    logger.info(f"Loaded assignments for {len(assignments)} cells "
                f"(status: {assignments['status'].value_counts().to_dict()})")
    return assignments


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
    #
    # scipy.io.mmread returns int64 for an integer MatrixMarket file — 8 bytes
    # per non-zero to hold UMI counts that never exceed a few tens of thousands.
    # Downcast at the door, before .tocsr() allocates a second int64 array:
    #
    #   * it halves the data array in every object downstream. On the 1.04M-cell
    #     merged atlas (2.78e9 non-zeros) that is 20.7GB -> 10.4GB, in the
    #     per-capture objects, both merges, and process_merged_data's
    #     `layers["counts"] = adata.X.copy()`.
    #   * float32 is what scanpy normalises into anyway, and what the CELLxGENE
    #     schema requires of the raw matrix — so format_cellxgene's promotion to
    #     raw.X becomes a no-op instead of a full conversion.
    #
    # float32 holds every integer exactly up to 2**24 (16,777,216); the largest
    # count observed across this registry is ~47k. Warn rather than round
    # silently if that ever stops being true.
    matrix = mmread(matrix_file)  # COO, features x cells
    if matrix.dtype != np.float32:
        largest = int(matrix.data.max()) if matrix.nnz else 0
        if largest > 2 ** 24:
            logger.warning(
                f"Largest count {largest:,} exceeds the float32 exact-integer limit "
                f"({2 ** 24:,}) — counts would be rounded; keeping {matrix.dtype}"
            )
        else:
            logger.info(
                f"Downcasting counts {matrix.dtype} -> float32 "
                f"({matrix.nnz:,} non-zeros, max {largest:,})"
            )
            matrix.data = matrix.data.astype(np.float32)
    matrix = matrix.T.tocsr()  # Transpose: cells x features

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


def read_gex(matrix_dir, capture: str) -> ad.AnnData:
    """Gene Expression counts of every called cell, as create_anndata builds them
    (float32, var_names = gene ids, var['gene_name']), with obs_names prefixed
    '<capture>_' like the per-capture objects."""
    modalities = read_10x_mtx_multimodal(Path(matrix_dir))
    adata = modalities["Gene Expression"] if "Gene Expression" in modalities \
        else list(modalities.values())[0]
    adata.obs_names = [f"{capture}_{bc}" for bc in adata.obs_names]
    return adata
