"""
workflow/scripts/doublet_rate.py
Choose scDblFinder's expected doublet rate (dbr / dbr.sd) per capture (#43).

dbr_source:
  loading   dbr.per1k from the capture's chip (captures.csv `chip` column;
            standard -> 0.008, ht -> 0.004 by default). scDblFinder then sets
            dbr = dbr.per1k * n_cells / 1000 and dbr.sd = 40% of dbr.
  genotype  dbr = d_total, dbr.sd = d_total_sd from snp_demux's
            doublet_summary.csv (swarbricklab/snp_demux#55). A lower bound.
            Falls back to `loading` for single-sample captures, a missing
            summary or NA values; the reason is recorded.
  none      dbr.sd = 1: no expectation on the number of doublets.

Pure functions only, so they can be unit-tested.
"""

from __future__ import annotations

import math

DEFAULT_PER1K = {"standard": 0.008, "ht": 0.004}
SCDBLFINDER_SD_FRACTION = 0.4   # scDblFinder's default dbr.sd = 40% of dbr
SOURCES = ("genotype", "loading", "none")

PARAM_COLUMNS = [
    "capture_id", "n_cells", "chip", "dbr_source", "dbr_source_used", "dbr_per1k",
    "dbr_used", "dbr_sd_used", "pass_dbr", "pass_dbr_sd", "fallback_reason",
]


def _num(x):
    """float(x), or None for missing / NA / non-numeric."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def per1k_for_chip(chip, per1k_map: dict | None = None) -> float:
    per1k_map = {**DEFAULT_PER1K, **(per1k_map or {})}
    key = str(chip).strip().lower() if chip is not None and str(chip).strip() else "standard"
    if key in ("nan", "na", "none"):
        key = "standard"
    if key not in per1k_map:
        raise ValueError(f"unknown chip {chip!r}; known: {sorted(per1k_map)}")
    return float(per1k_map[key])


def choose_dbr(capture: str, n_cells: int, source: str, chip=None,
               summary: dict | None = None, per1k_map: dict | None = None) -> dict:
    """Return the per-capture scDblFinder rate parameters as a PARAM_COLUMNS dict.

    pass_dbr / pass_dbr_sd say whether dbr / dbr.sd are passed explicitly; when
    False scDblFinder derives them (from dbr.per1k, and 40% of dbr).
    """
    if source not in SOURCES:
        raise ValueError(f"dbr_source={source!r}; use one of {SOURCES}")
    per1k = per1k_for_chip(chip, per1k_map)
    loading_dbr = per1k * n_cells / 1000.0
    row = dict(capture_id=capture, n_cells=int(n_cells),
               chip=str(chip) if chip is not None else "standard",
               dbr_source=source, dbr_per1k=per1k, fallback_reason="")

    if source == "genotype":
        reason = ""
        if not summary:
            reason = "no doublet_summary.csv"
        else:
            n_samples = _num(summary.get("n_samples"))
            d_total = _num(summary.get("d_total"))
            if n_samples is not None and n_samples <= 1:
                reason = "single-sample capture"
            elif d_total is None:
                reason = "d_total is NA"
        if not reason:
            sd = _num(summary.get("d_total_sd"))
            return {**row, "dbr_source_used": "genotype", "dbr_used": d_total,
                    "dbr_sd_used": sd if sd is not None else SCDBLFINDER_SD_FRACTION * d_total,
                    "pass_dbr": True, "pass_dbr_sd": sd is not None}
        row["fallback_reason"] = reason
        source = "loading"

    if source == "loading":
        return {**row, "dbr_source_used": "loading", "dbr_used": loading_dbr,
                "dbr_sd_used": SCDBLFINDER_SD_FRACTION * loading_dbr,
                "pass_dbr": False, "pass_dbr_sd": False}

    # none
    return {**row, "dbr_source_used": "none", "dbr_used": loading_dbr, "dbr_sd_used": 1.0,
            "pass_dbr": False, "pass_dbr_sd": True}
