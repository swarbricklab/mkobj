#!/usr/bin/env python3
"""
workflow/scripts/extract_known_doublets.py
Per-capture genotype-detected ("known") doublets for scDblFinder's knownDoublets
and the known_doublet column of doublets/<capture>.csv (#42).

Source, in order of preference:
  1. <demux>/<capture>/doublet_posteriors.csv (swarbricklab/snp_demux#53):
     consensus_status == 'doublet'
  2. the per-capture object's `status` column (create_anndata's resolution of
     cell_assignment.tsv + the AtLeastHalfSinglet consensus)
"""

import sys
import logging
from pathlib import Path

import anndata as ad
import pandas as pd

from doublet_io import known_from_posteriors, known_from_status, posteriors_resolve_doublets

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)

capture = snakemake.wildcards.capture
demux_root = snakemake.params.demux_root
post_path = Path(demux_root) / capture / "doublet_posteriors.csv" if demux_root else None

obs = ad.read_h5ad(snakemake.input.h5ad, backed="r").obs
post = None
fallback_reason = ""
if post_path is not None and post_path.exists():
    post = pd.read_csv(post_path, dtype=str)
    usable, why = posteriors_resolve_doublets(post)
    if not usable:
        logger.warning(f"Not using {post_path} for known doublets: {why}. "
                       "Falling back to the resolved demux `status` column.")
        post = None
        fallback_reason = f"{post_path} unusable"
else:
    fallback_reason = f"no {post_path}" if post_path else "no demux root"
if post is not None:
    known = known_from_posteriors(post, capture)
    logger.info(f"Known doublets from {post_path}")
    called = set(obs.index.str.replace(f"^{capture}_", "", regex=True))
    n_missing = int((~known["barcode"].isin(called)).sum())
    if n_missing:
        logger.warning(f"{n_missing} posterior barcodes are not called cells; ignored downstream")
else:
    known = known_from_status(obs, capture)
    logger.info(f"Known doublets from the demux `status` column ({fallback_reason})")

logger.info(f"{capture}: {int(known['known_doublet'].sum())} known doublets "
            f"of {len(known)} barcodes")
known.to_csv(snakemake.output.csv, index=False)
