#!/usr/bin/env python3
"""
workflow/scripts/cell_qc_features.py
Per-cell QC features of every called cell (#52): log total counts, log genes,
percent mito, percent ribo. Read from the Cell Ranger matrix (annotation stage).
"""

import sys
import logging

from capture_io import read_gex
from doublet_io import strip_capture_prefix
from within_lineage import qc_features

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
)
capture = snakemake.wildcards.capture
adata = read_gex(snakemake.input.matrix_dir, capture)
f = qc_features(adata.X, adata.var["gene_name"].values)
f.insert(0, "barcode", strip_capture_prefix(adata.obs_names, capture))
f.insert(0, "capture_id", capture)
f.to_csv(snakemake.output.csv, index=False)
logging.info(f"{capture}: QC features for {len(f)} cells")
