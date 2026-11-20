#!/usr/bin/env python3
"""
workflow/scripts/within_lineage_score.py
Cohort fit of the exploratory within-lineage doublet score (#52). Writes
doublets/within_lineage/separability.csv (per lineage: n, CV AUC, PU label
frequency c, estimated prevalence, coefficients) and per_capture.csv (predicted
vs expected within-donor same-lineage doublets). See within_lineage.py.
"""

import sys
import logging
from pathlib import Path

import pandas as pd

from doublet_io import read_calls
from doublet_report import capture_pairs
from within_lineage import fit_lineage, normalise_within, per_capture_estimates

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)
P = snakemake.params

cells, d_within, sources = [], {}, []
for qc_csv, known_csv, lin_csv, calls_csv in zip(snakemake.input.qc, snakemake.input.known,
                                                 snakemake.input.lineages, snakemake.input.calls):
    capture = Path(qc_csv).stem
    qc = pd.read_csv(qc_csv, dtype={"barcode": str})
    kn = pd.read_csv(known_csv, dtype={"barcode": str, "demux_status": str})
    lin = pd.read_csv(lin_csv, dtype={"barcode": str})
    # Same pairing as the concordance report (#48): scDblFinder's origin where
    # available (needs doublets.clusters), else the classifier's top two classes.
    pairs = capture_pairs(read_calls(calls_csv), lin, min_second_prob=P.min_second_prob,
                          pair_source=P.pair_source).reindex(qc["barcode"])
    kb = kn.loc[kn["known_doublet"].astype(bool), "barcode"]
    sources.extend(pairs["pair_source"].reindex(kb).dropna().tolist())
    df = (qc.merge(kn[["barcode", "known_doublet", "demux_status"]], on="barcode", how="left")
            .merge(lin[["barcode", "lineage_major"]], on="barcode", how="left"))
    df["pair_category"] = pairs["pair_category"].reindex(df["barcode"]).values
    cells.append(df)
    summ = Path(P.demux_root) / capture / "doublet_summary.csv" if P.demux_root else None
    if summ is not None and summ.exists():
        s = pd.read_csv(summ).iloc[0]
        if pd.notna(s.get("d_within")):
            d_within[capture] = float(s["d_within"])
src = pd.Series(sources).value_counts()
logger.info(f"Known-doublet lineage pairs by source: {src.to_dict()}")
if src.get("lineage_top2", 0) > src.get("scdblfinder_origin", 0):
    logger.warning("Most known-doublet pairs come from the classifier's top two classes, "
                   "which label nearly every real doublet same-lineage. Run scDblFinder with "
                   "params.doublets.clusters: lineage_fine so the labelled set is right.")
cells = normalise_within(pd.concat(cells, ignore_index=True), "capture_id")
cells["known_doublet"] = cells["known_doublet"].fillna(False).astype(bool)

sep_rows, scored = [], []
for lineage in sorted(cells["lineage_major"].dropna().unique()):
    row, sc = fit_lineage(cells, lineage, min_pos=P.min_pos)
    sep_rows.append(row)
    if sc is not None:
        scored.append(sc)
sep = pd.DataFrame(sep_rows)
sep["pair_source_used"] = ",".join(f"{k}={v}" for k, v in src.items()) or "none"
sep.to_csv(snakemake.output.separability, index=False)
logger.info("\n" + sep.to_string(index=False))

# Expected within-donor same-lineage doublets per capture x lineage:
# n_singlets x d_within x share of L|L among cross-donor known doublets (cohort
# template of random pairing).
known_pairs = cells[cells["known_doublet"]]
same = known_pairs[known_pairs["pair_category"] == "same_lineage"]
share = (same.groupby("lineage_major").size() / max(len(known_pairs), 1)).to_dict()
n_single = cells[cells["demux_status"].eq("singlet")].groupby("capture_id").size()
expected = pd.DataFrame([
    {"capture_id": c, "lineage": l,
     "expected_within_same_lineage": n_single.get(c, 0) * d_within[c] * share.get(l, 0.0)}
    for c in d_within for l in cells["lineage_major"].dropna().unique()
]) if d_within else None
est = per_capture_estimates(pd.concat(scored) if scored else None, expected)
est.to_csv(snakemake.output.per_capture, index=False)
logger.info("\n" + est.to_string(index=False))
