#!/usr/bin/env python3
"""
workflow/scripts/doublet_concordance_cohort.py
Cohort roll-up of the per-capture doublet concordance reports (#48).
"""

import sys
import logging

import pandas as pd

from doublet_report import cohort_rollup

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)

per_capture = pd.concat([pd.read_csv(p) for p in snakemake.input.reports], ignore_index=True)
cohort = cohort_rollup(per_capture, ratio_flag=snakemake.params.ratio_flag)
logger.info("\n" + cohort.to_string(index=False))
cohort.to_csv(snakemake.output.csv, index=False)
