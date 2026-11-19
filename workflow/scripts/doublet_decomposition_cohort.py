#!/usr/bin/env python3
"""
workflow/scripts/doublet_decomposition_cohort.py
Cohort-level doublet decomposition (#51; exploratory). Writes
doublets/decomposition/{per_capture,pair_enrichment,regression}.csv and a short
report.html. Method: docs/doublet_decomposition.md.
"""

import sys
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from doublet_decomposition import capture_summary, ols, pair_enrichment
from doublet_io import read_calls
from doublet_report import lineage_pairs

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)
P = snakemake.params
out_dir = Path(snakemake.output.per_capture).parent
out_dir.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(snakemake.input.capture_metadata, dtype={"capture_id": str}) \
    if snakemake.input.get("capture_metadata") else pd.DataFrame(columns=["capture_id"])
meta = meta.rename(columns={"capture": "capture_id"}).set_index("capture_id")

rows, groups = [], []
for calls_csv, known_csv, lin_csv in zip(snakemake.input.calls, snakemake.input.known,
                                         snakemake.input.lineages):
    capture = Path(calls_csv).stem
    calls, known = read_calls(calls_csv), pd.read_csv(known_csv, dtype={"barcode": str})
    if P.engine not in set(calls["engine"]):
        raise ValueError(f"{capture}: engine {P.engine!r} not in {calls_csv}")
    lin = pd.read_csv(lin_csv, dtype={"barcode": str})
    pairs = lineage_pairs(lin, pd.Index(calls["barcode"].unique()),
                          min_second_prob=P.min_second_prob)
    d_within = f_x = None
    summ = Path(P.demux_root) / capture / "doublet_summary.csv" if P.demux_root else None
    if summ is not None and summ.exists():
        s = pd.read_csv(summ).iloc[0]
        d_within = float(s["d_within"]) if pd.notna(s.get("d_within")) else None
        for col in (P.f_x_column, "f_x_posterior", "f_x_hard"):
            if col in s and pd.notna(s[col]):
                f_x = float(s[col])
                break
    row, grp = capture_summary(capture, calls, known, pairs, P.engine, d_within, f_x)
    row["f_x_source"] = "doublet_summary" if f_x is not None else "n_known / n_cells"
    row["dissociation_protocol"] = meta["dissociation_protocol"].get(capture, np.nan) \
        if "dissociation_protocol" in meta else np.nan
    loaded = meta["cells_loaded"].get(capture, np.nan) if "cells_loaded" in meta else np.nan
    row["cells_loaded"] = loaded if pd.notna(loaded) else row["n_cells"]
    row["cells_loaded_source"] = "capture_metadata" if pd.notna(loaded) else "called cells (proxy)"
    grp["dissociation_protocol"] = row["dissociation_protocol"]
    rows.append(row)
    groups.append(grp)

per_capture = pd.DataFrame(rows)
per_capture.to_csv(snakemake.output.per_capture, index=False)
groups = pd.concat(groups, ignore_index=True)

enr = [pair_enrichment(groups, P.adjacent_pairs).assign(protocol="all")]
for proto, g in groups.dropna(subset=["dissociation_protocol"]).groupby("dissociation_protocol"):
    enr.append(pair_enrichment(g, P.adjacent_pairs).assign(protocol=proto))
enrichment = pd.concat(enr, ignore_index=True)
enrichment.to_csv(snakemake.output.pair_enrichment, index=False)

# f_x ~ cells loaded (+ protocol): excess cross-donor doublets in co-dissociated
# captures beyond what loading explains point to cross-sample re-aggregation.
per_capture["cells_loaded_k"] = per_capture["cells_loaded"] / 1000
factors = ["dissociation_protocol"] if per_capture["dissociation_protocol"].nunique() > 1 else []
if not factors:
    logger.warning("dissociation_protocol missing or constant: regressing f_x on loading only")
reg = ols(per_capture, "f_x", ["cells_loaded_k"], factors)
reg.insert(0, "model", "f_x ~ cells_loaded_k" + (" + dissociation_protocol" if factors else ""))
reg.to_csv(snakemake.output.regression, index=False)
logger.info("\n" + reg.to_string(index=False))

# Short HTML summary.
fmt = lambda df: df.to_html(index=False, float_format=lambda x: f"{x:.4g}", border=0)
html = f"""<!doctype html><html><head><meta charset="utf-8"><title>Doublet decomposition</title>
<style>body{{font-family:system-ui,sans-serif;max-width:1100px;margin:2em auto;padding:0 1em}}
table{{border-collapse:collapse;font-size:13px;margin-bottom:2em}}td,th{{padding:3px 8px;border-bottom:1px solid #ddd;text-align:right}}
th{{background:#f4f4f4}}</style></head><body>
<h1>Doublet decomposition (exploratory)</h1>
<p>Engine: <b>{P.engine}</b>; {len(per_capture)} captures. Method and assumptions:
<code>docs/doublet_decomposition.md</code> in mkobj.</p>
<h2>Regression: cross-donor doublet fraction vs loading and protocol</h2>{fmt(reg)}
<h2>Lineage-pair enrichment of within-donor calls vs cross-donor template</h2>
{fmt(enrichment[enrichment.protocol == "all"].head(30))}
<h2>Per capture</h2>{fmt(per_capture.drop(columns=["cells_loaded_k"]))}
</body></html>"""
Path(snakemake.output.html).write_text(html)
logger.info(f"Wrote {out_dir}")
