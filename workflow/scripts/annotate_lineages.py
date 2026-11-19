#!/usr/bin/env python3
"""
workflow/scripts/annotate_lineages.py
Label every called cell of one capture with the reference lineage classifier
and write lineages/<capture_id>.csv (#45):

  capture_id, barcode, lineage_major, lineage_fine,
  confidence_major, confidence_fine, label_source

barcode is the raw Cell Ranger barcode (no capture prefix), matching
doublets/<capture_id>.csv.
"""

import sys
import logging

import anndata as ad
import pandas as pd

from doublet_io import strip_capture_prefix
from lineage_model import load_model, predict

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)

LINEAGE_COLUMNS = ["capture_id", "barcode", "lineage_major", "lineage_fine",
                   "confidence_major", "confidence_fine", "label_source"]

capture = snakemake.wildcards.capture
model = load_model(snakemake.params.model_dir)
info = model["info"]
logger.info(f"Model {info.get('version')} ({snakemake.params.model_dir}); levels "
            f"{ {k: len(v['classes']) for k, v in info['levels'].items()} }")

adata = ad.read_h5ad(snakemake.input.h5ad)
symbols = adata.var["gene_name"].values if "gene_name" in adata.var else adata.var_names.values
pred, n_matched = predict(model, adata.X, symbols, adata.var_names.values)
n_panel = len(info["genes"])
logger.info(f"{n_matched}/{n_panel} panel genes matched")
if n_matched < 0.8 * n_panel:
    logger.warning(f"Only {n_matched / n_panel:.0%} of panel genes matched; check the gene space")

out = pd.DataFrame({
    "capture_id": capture,
    "barcode": strip_capture_prefix(adata.obs_names, capture),
    **{c: pred[c].values for c in pred.columns},
})
for level in ("major", "fine"):
    if f"lineage_{level}" not in out:
        out[f"lineage_{level}"] = pd.NA
        out[f"confidence_{level}"] = pd.NA
out["label_source"] = "model"
out = out[LINEAGE_COLUMNS]

for level in ("major", "fine"):
    vc = out[f"lineage_{level}"].value_counts()
    logger.info(f"{level}: " + ", ".join(f"{k}={v}" for k, v in vc.items()))
out.to_csv(snakemake.output.csv, index=False)
logger.info(f"Wrote {len(out)} rows to {snakemake.output.csv}")
