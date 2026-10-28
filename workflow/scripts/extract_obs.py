#!/usr/bin/env python3
"""
workflow/scripts/extract_obs.py
Stage 1 of the config-gated keep-list filter (see rules/filter_deposit.smk).

Dump an AnnData's .obs to a flat CSV so the downstream DuckDB step can resolve a
keep-list without ever loading the counts matrix. The object is opened backed, so
this is cheap (obs only) even for a 50 GB raw-deposit merged.h5ad.

The row key is the obs index (obs_names) written as `cell_id`. For a merged
multi-capture object those names are capture-prefixed barcodes and so are unique
across captures — this is the join key the keep-list is resolved on, NOT the bare
16 bp cell barcode (which repeats across captures). Uniqueness is asserted here.
"""

import sys
import logging

import anndata as ad

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)

in_h5ad = snakemake.input.h5ad
out_csv = snakemake.output.obs_csv
# Optional column allow-list; default dumps every obs column.
obs_columns = snakemake.params.get("obs_columns") or None

logger.info(f"Reading obs (backed) from: {in_h5ad}")
adata = ad.read_h5ad(in_h5ad, backed="r")
obs = adata.obs.copy()
logger.info(f"  n_obs={adata.n_obs} n_vars={adata.n_vars}; obs columns: {list(obs.columns)}")

# obs_names is the unique per-cell key (capture-prefixed barcode).
if obs.index.has_duplicates:
    n_dup = int(obs.index.duplicated().sum())
    raise ValueError(
        f"obs index (cell_id) is not unique: {n_dup} duplicate names. "
        "The keep-list join requires a unique capture-prefixed cell id."
    )

if obs_columns:
    missing = [c for c in obs_columns if c not in obs.columns]
    if missing:
        raise KeyError(f"obs_columns requested but absent from .obs: {missing}")
    obs = obs[obs_columns]

obs.index.name = "cell_id"
# Cast everything to str so categorical / nullable-string columns round-trip as
# plain text for DuckDB (NaN -> empty, handled as NULL by read_csv_auto).
out = obs.astype(object).reset_index()
out.to_csv(out_csv, index=False)
logger.info(f"Wrote {len(out)} obs rows x {out.shape[1]} cols -> {out_csv}")
