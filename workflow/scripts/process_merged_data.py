#!/usr/bin/env python3
"""
workflow/scripts/process_merged_data.py
Normalize, select HVGs, compute PCA/UMAP, and cluster the merged AnnData object.

All computations use a temporary HVG-only subset. The resulting embeddings
(PCA, UMAP) and cluster labels are transferred back to the full-gene object
so that no genes are lost — critical for CELLxGENE compatibility.
"""

import sys
import logging
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path

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

# ── Inputs / outputs / params ────────────────────────────────────────────────────
input_h5ad = snakemake.input.merged
output_h5ad = snakemake.output.processed

n_top_genes = int(snakemake.params.get("n_top_genes", 2000))
hvg_batch_key = snakemake.params.get("hvg_batch_key", None) or None
n_pcs = int(snakemake.params.get("n_pcs", 50))
leiden_resolutions = snakemake.params.get("leiden_resolutions", [0.1, 0.3, 0.5, 1.0])

logger.info(f"Processing merged object: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad)
logger.info(f"Input shape: {adata.n_obs} cells x {adata.n_vars} genes")

# ── Back up raw counts ───────────────────────────────────────────────────────────
logger.info("Backing up raw counts to adata.layers['counts']")
adata.layers["counts"] = adata.X.copy()

# ── Normalize on the full object ─────────────────────────────────────────────────
logger.info("Normalizing (total-count + log1p)...")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# ── HVG selection ────────────────────────────────────────────────────────────────
logger.info(f"Selecting {n_top_genes} highly variable genes (batch_key={hvg_batch_key})...")
# Use only non-filtered genes for HVG selection if is_filtered exists
if 'is_filtered' in adata.var.columns:
    # Temporarily mark filtered genes as not highly variable
    mask_good = ~adata.var['is_filtered']
    adata_good = adata[:, mask_good].copy()
    sc.pp.highly_variable_genes(adata_good, n_top_genes=n_top_genes, batch_key=hvg_batch_key)
    # Transfer HVG annotations back to full object
    adata.var['highly_variable'] = False
    adata.var.loc[mask_good, 'highly_variable'] = adata_good.var['highly_variable'].values
    for col in ['means', 'dispersions', 'dispersions_norm']:
        if col in adata_good.var.columns:
            adata.var[col] = np.nan
            adata.var.loc[mask_good, col] = adata_good.var[col].values
    del adata_good
else:
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=hvg_batch_key)

n_hvg = adata.var['highly_variable'].sum()
logger.info(f"Selected {n_hvg} highly variable genes")

# ── PCA on HVG subset ───────────────────────────────────────────────────────────
logger.info(f"Running PCA (n_pcs={n_pcs}) on HVG subset...")
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.tl.pca(adata_hvg, n_comps=n_pcs)

# Transfer PCA back to full object
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.uns['pca'] = adata_hvg.uns['pca']
# Store which genes were used for PCA
adata.varm['PCs'] = np.zeros((adata.n_vars, n_pcs))
hvg_idx = np.where(adata.var['highly_variable'])[0]
adata.varm['PCs'][hvg_idx, :] = adata_hvg.varm['PCs']

# ── Neighbors + UMAP ────────────────────────────────────────────────────────────
logger.info("Computing neighbors and UMAP...")
# Use PCA from the subset but operate on the full object's obsm
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.umap(adata)

# ── Leiden clustering at multiple resolutions ────────────────────────────────────
for res in leiden_resolutions:
    key = f"leiden_res_{res:4.2f}"
    logger.info(f"Leiden clustering at resolution {res} -> {key}")
    sc.tl.leiden(adata, resolution=res, key_added=key, flavor="igraph", n_iterations=2)

# Clean up temporary subset
del adata_hvg

logger.info(f"Output shape: {adata.n_obs} cells x {adata.n_vars} genes")
logger.info(f"Layers: {list(adata.layers.keys())}")
logger.info(f"obsm keys: {list(adata.obsm.keys())}")

# ── Save ─────────────────────────────────────────────────────────────────────────
logger.info(f"Saving processed object to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(output_h5ad)
logger.info("Done.")
