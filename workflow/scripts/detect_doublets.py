#!/usr/bin/env python3
"""
workflow/scripts/detect_doublets.py
Detect doublets in a single capture using Scrublet and optionally cross-lineage markers.
Writes long-format calls only (doublets/engines/scrublet/<capture>.csv).
"""

import sys
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path

from capture_io import read_gex
from doublet_io import scrublet_calls, write_calls

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(snakemake.log[0]),
        logging.StreamHandler(sys.stderr),
    ],
)
logger = logging.getLogger(__name__)

def _excepthook(exc_type, exc_value, exc_tb):
    logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_tb))
    sys.__excepthook__(exc_type, exc_value, exc_tb)

sys.excepthook = _excepthook

# Inputs
matrix_dir = snakemake.input.matrix_dir
output_calls = snakemake.output.calls
markers_file = snakemake.params.get("lineage_markers", "")
scrublet_threshold = snakemake.params.get("scrublet_threshold", 0.2)
cross_lineage_threshold = snakemake.params.get("cross_lineage_threshold", 0.7)

capture = snakemake.wildcards.capture

logger.info(f"Detecting doublets for capture: {capture}")
# Every called cell, read straight from the Cell Ranger matrix exactly as
# create_anndata builds it (not from the temp per-capture object; see #49).
adata = read_gex(matrix_dir, capture)

# 1. Scrublet
logger.info("Running Scrublet...")
sc.pp.scrublet(adata)

# Ensure score vector is float and threshold
adata.obs['doublet_score'] = adata.obs['doublet_score'].astype(float)
adata.obs['predicted_scrublet'] = adata.obs['doublet_score'] > scrublet_threshold
logger.info(f"Found {adata.obs['predicted_scrublet'].sum()} Scrublet doublets out of {adata.n_obs} cells.")

# 2. Cross-lineage
adata.obs['predicted_cross_lineage_doublet'] = False
lineage_cols = []
if markers_file and Path(markers_file).exists():
    logger.info(f"Running cross-lineage detection using {markers_file}...")
    cepo_genes = pd.read_csv(markers_file, index_col=0).to_dict(orient='list')
    
    top_n_genes = 10
    cepo_genes2 = {}
    for cty, genes_of_interest in cepo_genes.items():
        # genes_of_interest = genes_of_interest[:20] is done in brca atlas
        genes_of_interest_set = set(genes_of_interest[:20])
        var_names_set = set(adata.var_names)
        intersection_genes = list(genes_of_interest_set & var_names_set)
        cepo_genes2[cty] = [g for g in genes_of_interest if g in intersection_genes][:top_n_genes]
    
    # Calculate proportion expressed
    for cty, genes in cepo_genes2.items():
        if not genes:
            adata.obs[cty] = 0.0
            continue
        gene_indices = [adata.var_names.get_loc(g) for g in genes]
        X_subset = adata.X[:, gene_indices] > 0
        prop_exprs = np.sum(X_subset, axis=1) / top_n_genes
        # handle np matrix conversions carefully
        adata.obs[cty] = prop_exprs.A1 if hasattr(prop_exprs, "A1") else np.squeeze(np.asarray(prop_exprs))

    # A cell is a cross-lineage doublet if it expresses >threshold for >1 lineage
    # Note: cepo_genes2.keys() may contain lineage names as columns in adata.obs
    lineage_cols = list(cepo_genes2.keys())
    doublet_counts = np.sum(adata.obs[lineage_cols] > cross_lineage_threshold, axis=1)
    adata.obs['predicted_cross_lineage_doublet'] = doublet_counts > 1
    adata.obs['doublet_cross_lineage_score'] = doublet_counts / 10.0
    logger.info(f"Found {adata.obs['predicted_cross_lineage_doublet'].sum()} cross-lineage doublets.")
else:
    logger.info("No lineage markers file provided or found. Skipping cross-lineage detection.")

# 3. Combine
adata.obs['predicted_doublet'] = adata.obs['predicted_scrublet'] | adata.obs['predicted_cross_lineage_doublet']
total_predicted = adata.obs['predicted_doublet'].sum()
logger.info(f"Total predicted doublets combined: {total_predicted} ({(total_predicted/adata.n_obs)*100:.2f}%)")

# Persist the Scrublet calls in the long-format schema (#40); combined across
# engines into doublets/<capture>.csv by combine_doublet_calls. This CSV is the
# only output: the merge restores the obs columns from it (#49), so the
# annotation stage and the assembly stage can run separately.
calls = scrublet_calls(adata.obs, capture, lineage_cols)
Path(output_calls).parent.mkdir(parents=True, exist_ok=True)
write_calls(calls, output_calls)
logger.info(f"Wrote {len(calls)} Scrublet call rows to {output_calls}")
logger.info("Done.")
