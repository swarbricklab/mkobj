#!/usr/bin/env python3
"""
workflow/scripts/annotate_lineages.py
Label every called cell of one capture and write lineages/<capture_id>.csv:

  capture_id, barcode, lineage_major, lineage_fine,
  confidence_major, confidence_fine, label_source
  [+ lineage_major_crosswalk when source = both]

params.lineages.source:
  model      reference classifier (#45)
  crosswalk  cellranger annotate label -> major lineage via the reference's crosswalk
             (#46); fine labels and confidences are NA, label_source = crosswalk
  both       model labels, plus lineage_major_crosswalk for comparison and a
             model-vs-crosswalk agreement table (lineages/agreement/<capture>.csv)

barcode is the raw Cell Ranger barcode (no capture prefix), matching
doublets/<capture_id>.csv.
"""

import sys
import logging
from pathlib import Path

import anndata as ad
import pandas as pd

from doublet_io import strip_capture_prefix
from lineage_crosswalk import agreement_table, apply_crosswalk, load_crosswalk
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
source = snakemake.params.source
adata = ad.read_h5ad(snakemake.input.h5ad)
barcodes = strip_capture_prefix(adata.obs_names, capture)
out = pd.DataFrame({"capture_id": capture, "barcode": barcodes})
model_classes = None

# ── Model labels ─────────────────────────────────────────────────────────────────
if source in ("model", "both"):
    model = load_model(snakemake.params.model_dir)
    info = model["info"]
    model_classes = set(info["levels"].get("major", {}).get("classes", []))
    logger.info(f"Model {info.get('version')} ({snakemake.params.model_dir}); levels "
                f"{ {k: len(v['classes']) for k, v in info['levels'].items()} }")
    symbols = adata.var["gene_name"].values if "gene_name" in adata.var else adata.var_names.values
    pred, n_matched = predict(model, adata.X, symbols, adata.var_names.values)
    n_panel = len(info["genes"])
    logger.info(f"{n_matched}/{n_panel} panel genes matched")
    if n_matched < 0.8 * n_panel:
        logger.warning(f"Only {n_matched / n_panel:.0%} of panel genes matched; check the gene space")
    for c in pred.columns:
        out[c] = pred[c].values
    out["label_source"] = "model"

# ── Crosswalk labels ─────────────────────────────────────────────────────────────
if source in ("crosswalk", "both"):
    xw = load_crosswalk(snakemake.input.crosswalk)
    label_col = snakemake.params.label_column
    ann_path = Path(snakemake.params.annotation)
    if snakemake.params.annotation and ann_path.exists():
        ann = pd.read_csv(ann_path, dtype=str)
        if label_col not in ann.columns:
            raise ValueError(f"{ann_path} has no {label_col!r} column "
                             f"(params.lineages.crosswalk_label_column); has {list(ann.columns)}")
        labels = ann.set_index("barcode")[label_col].reindex(barcodes).reset_index(drop=True)
    else:
        logger.warning(f"No cellranger annotate labels at {ann_path}; crosswalk labels are NA")
        labels = pd.Series([pd.NA] * len(barcodes), dtype="string")
    mapped, unmapped = apply_crosswalk(labels, xw)
    logger.info(f"Crosswalk ({snakemake.input.crosswalk}, column {label_col}): "
                f"{int(mapped.notna().sum())}/{len(mapped)} cells mapped")
    if len(unmapped):
        logger.warning("Unmapped cellranger labels (-> NA): "
                       + "; ".join(f"{k} ({v})" for k, v in unmapped.items()))
    if model_classes:
        foreign = sorted(set(xw.values()) - model_classes
                         - set(snakemake.params.collapse.values()))
        if foreign:
            logger.warning(f"Crosswalk lineages not in the model's major vocabulary: {foreign}")

    if source == "crosswalk":
        out["lineage_major"] = mapped.values
        out["lineage_fine"] = pd.NA
        out["confidence_major"] = pd.NA
        out["confidence_fine"] = pd.NA
        out["label_source"] = "crosswalk"
    else:
        out["lineage_major_crosswalk"] = mapped.values
        tab = agreement_table(capture, out["lineage_major"], mapped,
                              collapse=snakemake.params.collapse)
        n_cmp = tab.loc[tab.lineage_major_crosswalk != "NA", "n_cells"].sum()
        n_ok = tab.loc[tab.agree, "n_cells"].sum()
        logger.info(f"Model vs crosswalk agreement: {n_ok}/{n_cmp} "
                    f"({n_ok / max(n_cmp, 1):.1%}) of crosswalk-labelled cells")
        Path(snakemake.output.agreement).parent.mkdir(parents=True, exist_ok=True)
        tab.to_csv(snakemake.output.agreement, index=False)

for level in ("major", "fine"):
    for prefix in ("lineage", "confidence"):
        if f"{prefix}_{level}" not in out:
            out[f"{prefix}_{level}"] = pd.NA
    vc = out[f"lineage_{level}"].value_counts(dropna=False)
    logger.info(f"{level}: " + ", ".join(f"{k}={v}" for k, v in vc.items()))
extra = [c for c in out.columns if c not in LINEAGE_COLUMNS]
out[LINEAGE_COLUMNS + extra].to_csv(snakemake.output.csv, index=False)
logger.info(f"Wrote {len(out)} rows to {snakemake.output.csv}")
