#!/usr/bin/env python3
"""
workflow/scripts/doublet_concordance.py
Per-capture doublet concordance report (#48): engine calls vs genotype doublets,
by lineage pair. See doublet_report.py for definitions.
"""

import sys
import logging
from pathlib import Path

import pandas as pd

from doublet_io import read_calls
from doublet_report import capture_pairs, concordance

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)

capture = snakemake.wildcards.capture
calls = read_calls(snakemake.input.calls)
known = pd.read_csv(snakemake.input.known, dtype={"barcode": str, "demux_status": str})

lineages = None
if snakemake.input.get("lineages"):
    lineages = pd.read_csv(snakemake.input.lineages, dtype={"barcode": str})
else:
    logger.warning("No lineage labels: every known doublet is pair_category 'unknown'")
pairs = capture_pairs(calls, lineages, min_second_prob=snakemake.params.min_second_prob,
                      pair_source=snakemake.params.pair_source)

d_within = None
demux_root = snakemake.params.demux_root
summary = Path(demux_root) / capture / "doublet_summary.csv" if demux_root else None
if summary is not None and summary.exists():
    s = pd.read_csv(summary)
    if "d_within" in s.columns and pd.notna(s["d_within"].iloc[0]):
        d_within = float(s["d_within"].iloc[0])
logger.info(f"d_within = {d_within} ({summary if d_within is not None else 'not available'})")

trained = ("scdblfinder",) if snakemake.params.known_use == "positive" else ()
if trained:
    logger.warning("scDblFinder ran with knownUse=positive: its recall_known is ~1 by "
                   "construction (known doublets were training positives)")
rep = concordance(capture, calls, known, pairs, d_within=d_within,
                  ratio_flag=snakemake.params.ratio_flag, trained_on_known=trained)
logger.info("\n" + rep[["engine", "n_called", "n_known", "recall_known",
                        "recall_known_cross_lineage", "recall_known_same_lineage",
                        "n_called_singlets", "singlet_flag"]].to_string(index=False))
Path(snakemake.output.csv).parent.mkdir(parents=True, exist_ok=True)
rep.to_csv(snakemake.output.csv, index=False)
