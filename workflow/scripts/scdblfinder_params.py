#!/usr/bin/env python3
"""
workflow/scripts/scdblfinder_params.py
Per-capture scDblFinder parameters (#43): the expected doublet rate (dbr /
dbr.sd) chosen by params.doublets.dbr_source, recorded as
doublets/scdblfinder_params/<capture>.csv and read by
detect_doublets_scdblfinder.R. See doublet_rate.py for the rules.
"""

import gzip
import sys
import logging
from pathlib import Path

import pandas as pd

from doublet_clusters import cluster_labels
from doublet_rate import PARAM_COLUMNS, choose_dbr

CLUSTER_COLUMNS = ["clusters", "n_clusters", "n_fine_merged", "min_cluster_size"]

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)

capture = snakemake.wildcards.capture
source = snakemake.params.dbr_source

with gzip.open(Path(snakemake.input.matrix_dir) / "barcodes.tsv.gz", "rt") as fh:
    n_cells = sum(1 for line in fh if line.strip())

chip = None
captures_file = snakemake.params.captures_file
if captures_file and Path(captures_file).exists():
    caps = pd.read_csv(captures_file, dtype=str)
    if "chip" in caps.columns:
        hit = caps.loc[caps["capture"] == capture, "chip"]
        chip = hit.iloc[0] if len(hit) else None

summary = None
demux_root = snakemake.params.demux_root
if source == "genotype" and demux_root:
    path = Path(demux_root) / capture / "doublet_summary.csv"
    if path.exists():
        df = pd.read_csv(path, dtype={"capture_id": str})
        if "capture_id" in df.columns:
            df = df[df["capture_id"] == capture]
        summary = df.iloc[0].to_dict() if len(df) else None
        logger.info(f"doublet_summary: {path}")

row = choose_dbr(capture, n_cells, source, chip=chip, summary=summary,
                 per1k_map=snakemake.params.per1k)
if row["fallback_reason"]:
    logger.warning(f"{capture}: dbr_source=genotype -> loading ({row['fallback_reason']})")

# ── Clusters for artificial-doublet simulation (#47) ────────────────────────────
mode = snakemake.params.clusters
min_size = int(snakemake.params.min_cluster_size)
lineages = None
if mode != "none":
    lineages = pd.read_csv(snakemake.input.lineages, dtype={"barcode": str})
labels, summary = cluster_labels(lineages, mode, min_size) if lineages is not None \
    else cluster_labels(pd.DataFrame(), "none")
row.update(summary, min_cluster_size=min_size)
if labels is not None:
    pd.DataFrame({"barcode": lineages["barcode"], "cluster": labels.values}).to_csv(
        snakemake.output.clusters, index=False)
    logger.info(f"clusters={mode}: {summary['n_clusters']} clusters "
                f"({summary['n_fine_merged']} small fine types folded into major): "
                + ", ".join(f"{k}={v}" for k, v in labels.value_counts().items()))
else:
    pd.DataFrame(columns=["barcode", "cluster"]).to_csv(snakemake.output.clusters, index=False)

cols = PARAM_COLUMNS + CLUSTER_COLUMNS
logger.info(", ".join(f"{k}={row[k]}" for k in cols))
pd.DataFrame([row], columns=cols).to_csv(snakemake.output.csv, index=False)
