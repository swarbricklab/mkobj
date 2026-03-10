#!/usr/bin/env python3
"""
workflow/scripts/format_cellxgene.py
Format a processed AnnData object for CELLxGENE upload.

Ensures the object has the required obs/var columns per the CELLxGENE schema
(https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/).

This script maps existing metadata to schema-required fields and fills in
configurable defaults so the output is ready (or nearly ready) for submission.
"""

import sys
import logging
import pandas as pd
import numpy as np
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
input_h5ad = snakemake.input.processed
output_h5ad = snakemake.output.cellxgene

# CELLxGENE defaults — these can be overridden via config
cellxgene_cfg = snakemake.params.get("cellxgene", {}) or {}

# Column mappings: schema field -> existing obs column name
# If the obs column already exists under the schema name, no mapping is needed.
column_map = cellxgene_cfg.get("column_map", {}) or {}

# Default values for required schema fields when no mapping exists
defaults = cellxgene_cfg.get("defaults", {}) or {}

logger.info(f"Formatting for CELLxGENE: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad)
logger.info(f"Input shape: {adata.n_obs} cells x {adata.n_vars} genes")

# ── Required obs fields per CELLxGENE schema ─────────────────────────────────────
# https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md
REQUIRED_OBS_FIELDS = [
    "assay_ontology_term_id",
    "cell_type_ontology_term_id",
    "development_stage_ontology_term_id",
    "disease_ontology_term_id",
    "donor_id",
    "is_primary_data",
    "organism_ontology_term_id",
    "self_reported_ethnicity_ontology_term_id",
    "sex_ontology_term_id",
    "suspension_type",
    "tissue_ontology_term_id",
]

# ── Apply column mappings ────────────────────────────────────────────────────────
for schema_field, obs_column in column_map.items():
    if obs_column in adata.obs.columns:
        logger.info(f"  Mapping obs['{obs_column}'] -> obs['{schema_field}']")
        adata.obs[schema_field] = adata.obs[obs_column].values
    else:
        logger.warning(f"  column_map: obs column '{obs_column}' not found for schema field '{schema_field}'")

# ── Fill defaults for missing required fields ────────────────────────────────────
for field in REQUIRED_OBS_FIELDS:
    if field not in adata.obs.columns:
        if field in defaults:
            val = defaults[field]
            logger.info(f"  Setting obs['{field}'] = '{val}' (from config default)")
            adata.obs[field] = val
        else:
            logger.warning(f"  Required field '{field}' is missing and no default provided — setting to 'unknown'")
            adata.obs[field] = "unknown"

# ── Ensure is_primary_data is boolean ────────────────────────────────────────────
if "is_primary_data" in adata.obs.columns:
    adata.obs["is_primary_data"] = adata.obs["is_primary_data"].astype(bool)

# ── var: rename is_filtered -> feature_is_filtered (CELLxGENE convention) ────────
if "is_filtered" in adata.var.columns and "feature_is_filtered" not in adata.var.columns:
    logger.info("  Renaming var['is_filtered'] -> var['feature_is_filtered']")
    adata.var["feature_is_filtered"] = adata.var["is_filtered"]

# Ensure feature_is_filtered exists
if "feature_is_filtered" not in adata.var.columns:
    logger.info("  Setting var['feature_is_filtered'] = False (no gene filtering info found)")
    adata.var["feature_is_filtered"] = False

# ── var: ensure feature_biotype column exists ────────────────────────────────────
if "feature_biotype" not in adata.var.columns:
    logger.info("  Setting var['feature_biotype'] = 'gene'")
    adata.var["feature_biotype"] = "gene"

# ── var: ensure feature_reference column exists ──────────────────────────────────
if "feature_reference" not in adata.var.columns:
    ref = cellxgene_cfg.get("feature_reference", "NCBITaxon:9606")
    logger.info(f"  Setting var['feature_reference'] = '{ref}'")
    adata.var["feature_reference"] = ref

# ── Validation summary ───────────────────────────────────────────────────────────
logger.info("Validation summary:")
missing = []
unknown_fields = []
for field in REQUIRED_OBS_FIELDS:
    if field not in adata.obs.columns:
        missing.append(field)
    elif (adata.obs[field] == "unknown").all():
        unknown_fields.append(field)

if missing:
    logger.warning(f"  MISSING required obs fields: {missing}")
else:
    logger.info("  All required obs fields present.")

if unknown_fields:
    logger.warning(f"  Fields set to 'unknown' (need real values before upload): {unknown_fields}")

logger.info(f"  var['feature_is_filtered']: {adata.var['feature_is_filtered'].sum()} / {adata.n_vars} genes flagged")
logger.info(f"  layers: {list(adata.layers.keys())}")
logger.info(f"  obsm: {list(adata.obsm.keys())}")

# ── Save ─────────────────────────────────────────────────────────────────────────
logger.info(f"Saving CELLxGENE-formatted object to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(output_h5ad)
logger.info("Done.")
