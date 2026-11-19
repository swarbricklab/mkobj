#!/usr/bin/env python3
"""
workflow/scripts/combine_doublet_calls.py
Row-bind each engine's per-capture calls into the persisted long-format
doublets/<capture_id>.csv (#40). Schema: see doublet_io.py.
"""

import sys
import logging

import pandas as pd

from doublet_io import attach_known, combine_calls, read_calls, write_calls

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)

frames = []
for path in snakemake.input.calls:
    df = read_calls(path)
    logger.info(f"{path}: {len(df)} rows, engines {sorted(df['engine'].unique())}")
    frames.append(df)

combined = combine_calls(frames)
if snakemake.input.get("known"):
    known = pd.read_csv(snakemake.input.known, dtype={"capture_id": str, "barcode": str})
    combined = attach_known(combined, known)
    logger.info(f"Attached known_doublet ({int(known['known_doublet'].sum())} known)")
counts = combined.groupby("engine")["doublet_call"].agg(["size", "sum"])
for engine, row in counts.iterrows():
    logger.info(f"{engine}: {int(row['sum'])} / {int(row['size'])} called doublets")
write_calls(combined, snakemake.output.csv)
logger.info(f"Wrote {len(combined)} rows to {snakemake.output.csv}")
