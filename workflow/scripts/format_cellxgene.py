#!/usr/bin/env python3
"""
workflow/scripts/format_cellxgene.py
Format a processed AnnData object for CELLxGENE upload.

Targets schema 7.1.0
(https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.1.0/schema.md).

Three sources fill the schema fields, applied in this order:

  1. `cellxgene.metadata` — a list of CSV join tables. Each is keyed on an
     existing obs column and contributes its remaining columns to obs. This
     is how per-sample clinical metadata, per-capture assay chemistry and
     per-cell cell-type ontology maps get in.
  2. `cellxgene.column_map` — rename an obs column onto a schema field.
  3. `cellxgene.defaults` — a constant, used both for an absent column and to
     fill rows a join did not match.

Anything still unset lands on "unknown", which the schema permits for every
ontology field that can be unavailable.
"""

import sys
import logging
import numpy as np
import pandas as pd
import scipy.sparse as sp
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

SCHEMA_VERSION = "7.1.0"

# ── Inputs / outputs / params ────────────────────────────────────────────────────
input_h5ad = snakemake.input.processed
output_h5ad = snakemake.output.cellxgene

cellxgene_cfg = snakemake.params.get("cellxgene", {}) or {}
subset = snakemake.params.get("subset", "") or ""

metadata_tables = cellxgene_cfg.get("metadata", []) or []
column_map = cellxgene_cfg.get("column_map", {}) or {}
defaults = cellxgene_cfg.get("defaults", {}) or {}
uns_cfg = cellxgene_cfg.get("uns", {}) or {}

# The cellxgene block is global, but uns['title'] must distinguish the datasets
# within a collection — and with `subsets:` one run emits several. Allow the
# literal {subset} in any uns string so a single config can title them all.
# Plain str.replace, not str.format, so other braces pass through untouched.
if subset:
    uns_cfg = {
        k: v.replace("{subset}", subset) if isinstance(v, str) else v
        for k, v in uns_cfg.items()
    }

logger.info(f"Formatting for CELLxGENE schema {SCHEMA_VERSION}: {input_h5ad}")
adata = ad.read_h5ad(input_h5ad)
logger.info(f"Input shape: {adata.n_obs} cells x {adata.n_vars} genes")

# ── Schema field inventory ───────────────────────────────────────────────────────
# obs fields the curator MUST annotate. organism_ontology_term_id is NOT here:
# it moved from obs to uns in schema 6.0.0.
#
# Deliberately omitted:
#   experimental_condition_ontology_term_id — MUST NOT be present when every
#       observation would be "na", which is the case for an unperturbed dataset.
#   genetic_perturbation_id / _strategy     — likewise, gated on uns['genetic_perturbations'].
#   array_col / array_row / in_tissue       — Visium only.
REQUIRED_OBS_FIELDS = [
    "assay_ontology_term_id",
    "cell_type_ontology_term_id",
    "development_stage_ontology_term_id",
    "disease_ontology_term_id",
    "donor_id",
    "is_primary_data",
    "self_reported_ethnicity_ontology_term_id",
    "sex_ontology_term_id",
    "suspension_type",
    "tissue_ontology_term_id",
    "tissue_type",
]

# uns fields the curator MUST annotate.
REQUIRED_UNS_FIELDS = [
    "organism_ontology_term_id",
    "title",
]

# var columns CELLxGENE Discover annotates itself — a curator MUST NOT supply
# them, so we strip any that earlier steps (or an older mkobj) wrote.
DISCOVER_VAR_COLUMNS = [
    "feature_biotype",
    "feature_length",
    "feature_name",
    "feature_reference",
    "feature_type",
]

BOOL_OBS_FIELDS = {"is_primary_data"}

# ── 1. Join metadata tables ──────────────────────────────────────────────────────
# Each entry is {path, key}, where `key` names an existing obs column. Every
# other column in the CSV is written to obs under its own name. Rows that do
# not match are left NaN and picked up by the defaults pass below.
#
# Tables are applied in order, so a later table overrides an earlier one on a
# shared column name.
for spec in metadata_tables:
    path = spec.get("path")
    key = spec.get("key")
    if not path or not key:
        logger.warning(f"  metadata: skipping malformed entry {spec!r} (need 'path' and 'key')")
        continue

    if key not in adata.obs.columns:
        logger.warning(f"  metadata: obs has no column '{key}' — skipping {path}")
        continue

    table = pd.read_csv(path, dtype=str)
    if key not in table.columns:
        logger.warning(f"  metadata: {path} has no column '{key}' — skipping")
        continue

    n_dup = table.duplicated(subset=[key]).sum()
    if n_dup:
        logger.warning(f"  metadata: {path} has {n_dup} duplicate '{key}' rows — keeping the first of each")
        table = table.drop_duplicates(subset=[key])

    # Compare as strings on both sides; obs keys are often categorical.
    left = adata.obs[key].astype(str).where(adata.obs[key].notna())
    value_columns = [c for c in table.columns if c != key]

    matched = left.isin(set(table[key]))
    n_null_key = int(left.isna().sum())
    unmatched_keys = sorted(left[~matched].dropna().unique())
    logger.info(
        f"  metadata: {path} on obs['{key}'] — "
        f"{matched.sum():,}/{adata.n_obs:,} cells matched, "
        f"contributing {len(value_columns)} columns"
    )
    if n_null_key:
        logger.info(f"    {n_null_key:,} cells have no '{key}' value to join on")
    if unmatched_keys:
        shown = unmatched_keys[:10]
        logger.warning(
            f"    {len(unmatched_keys)} '{key}' values absent from the table: {shown}"
            + (" ..." if len(unmatched_keys) > len(shown) else "")
        )

    for col in value_columns:
        adata.obs[col] = left.map(dict(zip(table[key], table[col])))

# ── 2. Apply column mappings ─────────────────────────────────────────────────────
for schema_field, obs_column in column_map.items():
    if obs_column in adata.obs.columns:
        logger.info(f"  Mapping obs['{obs_column}'] -> obs['{schema_field}']")
        adata.obs[schema_field] = adata.obs[obs_column].values
    else:
        logger.warning(f"  column_map: obs column '{obs_column}' not found for schema field '{schema_field}'")

# ── 3. Fill required obs fields ──────────────────────────────────────────────────
# A field can be absent entirely, or present with gaps left by an unmatched
# join. Both get the configured default, or "unknown" if none is configured.
for field in REQUIRED_OBS_FIELDS:
    fallback = defaults.get(field, "unknown")

    if field not in adata.obs.columns:
        source = "config default" if field in defaults else "no value available"
        logger.info(f"  Setting obs['{field}'] = {fallback!r} for all cells ({source})")
        adata.obs[field] = fallback
        continue

    if field in BOOL_OBS_FIELDS:
        continue

    # Categorical columns reject an unseen fill value, so widen first.
    if isinstance(adata.obs[field].dtype, pd.CategoricalDtype):
        adata.obs[field] = adata.obs[field].astype(object)

    n_missing = adata.obs[field].isna().sum()
    if n_missing:
        logger.info(
            f"  obs['{field}']: filling {n_missing:,}/{adata.n_obs:,} unset cells with {fallback!r}"
        )
        adata.obs[field] = adata.obs[field].fillna(fallback)

# Schema wants str categoricals for these; is_primary_data is the one bool.
for field in REQUIRED_OBS_FIELDS:
    if field in BOOL_OBS_FIELDS:
        adata.obs[field] = adata.obs[field].astype(bool)
    else:
        adata.obs[field] = adata.obs[field].astype(str).astype("category")

# ── 4. uns fields ────────────────────────────────────────────────────────────────
for field in REQUIRED_UNS_FIELDS:
    if field in uns_cfg:
        adata.uns[field] = uns_cfg[field]
        logger.info(f"  Setting uns['{field}'] = {uns_cfg[field]!r}")
    elif field in adata.uns:
        logger.info(f"  uns['{field}'] already set to {adata.uns[field]!r}")
    else:
        logger.warning(f"  Required uns field '{field}' is unset — add it under cellxgene.uns")

# organism_ontology_term_id moved obs -> uns in schema 6.0.0. Migrate a value
# left in obs by an older config rather than silently dropping it.
if "organism_ontology_term_id" in adata.obs.columns:
    obs_values = set(adata.obs["organism_ontology_term_id"].astype(str).unique())
    if "organism_ontology_term_id" not in adata.uns and len(obs_values) == 1:
        adata.uns["organism_ontology_term_id"] = obs_values.pop()
        logger.info(
            f"  Moved organism_ontology_term_id from obs to uns "
            f"(= {adata.uns['organism_ontology_term_id']!r}); it is a uns field from schema 6.0.0"
        )
    else:
        logger.info("  Dropping obs['organism_ontology_term_id'] — it belongs in uns")
    del adata.obs["organism_ontology_term_id"]

# ── 5. raw counts ────────────────────────────────────────────────────────────────
# For UMI scRNA-seq the raw matrix is REQUIRED, and MUST live in raw.X when a
# normalized matrix is in X. mkobj carries counts in layers['counts'], so
# promote them. The layer stays put — extra layers are allowed.
if adata.raw is None and "counts" in adata.layers:
    counts = adata.layers["counts"]
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    counts = counts.astype(np.float32)
    adata.raw = ad.AnnData(X=counts, obs=adata.obs[[]].copy(), var=adata.var[[]].copy())
    logger.info("  Populated raw.X from layers['counts'] as float32 (schema requires raw in raw.X)")
elif adata.raw is None:
    logger.warning("  No raw counts available — schema REQUIRES a raw matrix for UMI scRNA-seq")

# ── 6. var ───────────────────────────────────────────────────────────────────────
# 10x's own feature_type ("Gene Expression" / "Antibody Capture") collides with
# a Discover-reserved column name, so rename rather than lose it.
if "feature_type" in adata.var.columns and "cellranger_feature_type" not in adata.var.columns:
    adata.var["cellranger_feature_type"] = adata.var["feature_type"]
    logger.info("  Renamed var['feature_type'] -> var['cellranger_feature_type'] (reserved name)")

dropped = [c for c in DISCOVER_VAR_COLUMNS if c in adata.var.columns]
if dropped:
    adata.var = adata.var.drop(columns=dropped)
    logger.info(f"  Dropped Discover-annotated var columns: {dropped}")

# feature_is_filtered means "zeroed out in X but retained in raw", which lets
# Explorer mask the gene. mkobj's own `is_filtered` flags low-expression genes
# WITHOUT zeroing them, so the two are not interchangeable — reusing it here
# would claim masking that has not happened. Report only genes that really are
# all-zero in X while non-zero in raw.
if adata.raw is not None:
    x_nonzero = np.asarray(abs(adata.X).sum(axis=0)).ravel() > 0
    raw_nonzero = np.asarray(abs(adata.raw.X).sum(axis=0)).ravel() > 0
    feature_is_filtered = (~x_nonzero) & raw_nonzero
else:
    # "When a raw matrix is not present, the value for all features MUST be False."
    feature_is_filtered = np.zeros(adata.n_vars, dtype=bool)

adata.var["feature_is_filtered"] = feature_is_filtered
logger.info(
    f"  var['feature_is_filtered']: {int(feature_is_filtered.sum())}/{adata.n_vars} genes "
    f"zeroed in X but present in raw"
)
if "is_filtered" in adata.var.columns:
    logger.info(
        f"  var['is_filtered'] retained as a non-schema flag "
        f"({int(adata.var['is_filtered'].sum())} low-expression genes, not masked)"
    )

# ── 7. Validation summary ────────────────────────────────────────────────────────
logger.info("Validation summary:")

unknown_fields = [
    f for f in REQUIRED_OBS_FIELDS
    if f not in BOOL_OBS_FIELDS and (adata.obs[f].astype(str) == "unknown").all()
]
partial_fields = [
    (f, int((adata.obs[f].astype(str) == "unknown").sum()))
    for f in REQUIRED_OBS_FIELDS
    if f not in BOOL_OBS_FIELDS and f not in unknown_fields
    and (adata.obs[f].astype(str) == "unknown").any()
]

logger.info(f"  All {len(REQUIRED_OBS_FIELDS)} required obs fields present.")
if unknown_fields:
    logger.warning(f"  Entirely 'unknown' (need real values before upload): {unknown_fields}")
for field, n in partial_fields:
    logger.info(f"  obs['{field}']: {n:,}/{adata.n_obs:,} cells 'unknown'")

missing_uns = [f for f in REQUIRED_UNS_FIELDS if f not in adata.uns]
if missing_uns:
    logger.warning(f"  MISSING required uns fields: {missing_uns}")

# Cheap structural checks the validator would otherwise fail on later.
if adata.raw is not None:
    raw_counts = adata.raw.X
    empty_cells = int((np.asarray(raw_counts.sum(axis=1)).ravel() == 0).sum())
    if empty_cells:
        logger.warning(f"  {empty_cells:,} cells have an all-zero raw matrix (schema requires >= 1 non-zero)")
    sample = raw_counts[:1000] if raw_counts.shape[0] > 1000 else raw_counts
    if not np.allclose(sample.data, np.rint(sample.data)):
        logger.warning("  raw.X contains non-integer values (schema requires de-duplicated molecule counts)")

logger.info(f"  layers: {list(adata.layers.keys())}")
logger.info(f"  obsm: {list(adata.obsm.keys())}")
logger.info(f"  uns: {sorted(adata.uns.keys())}")

# ── Save ─────────────────────────────────────────────────────────────────────────
logger.info(f"Saving CELLxGENE-formatted object to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(output_h5ad)
