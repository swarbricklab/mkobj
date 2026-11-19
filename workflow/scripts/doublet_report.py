"""
workflow/scripts/doublet_report.py
Doublet concordance: expression-engine calls vs genotype doublets (#48).

Per capture x engine:
  recall on known (genotype) doublets: overall, cross-lineage, same-lineage
  calls among genotype singlets vs the expected n_singlets x d_within
  (swarbricklab/snp_demux#55), with a ratio and a flag
  agreement with every other engine (Jaccard, Cohen's kappa)

Lineage pair of a barcode (pair_source):
  lineage_top2 (default)  the classifier's top two major classes
       (lineage_major, lineage_major_2) when the runner-up has probability
       >= min_second_prob, otherwise (lineage_major, lineage_major). Falls back
       to the scDblFinder origin when the labels have no runner-up column
       (e.g. lineages.source: crosswalk).
  scdblfinder_origin  scDblFinder's most likely origin ("A+B", clusters mapped
       to their major lineage; needs doublets.clusters), else lineage_top2.
       On the synthetic fixture the origin got the category right for 68-81%
       of known doublets vs 100% for lineage_top2, hence the default.
  cross_lineage when the two majors differ, same_lineage when they match,
  unknown without labels.

Pure functions only, so they can be unit-tested.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

CATEGORIES = ("cross_lineage", "same_lineage", "unknown")


PAIR_SOURCES = ("lineage_top2", "scdblfinder_origin")


def lineage_pairs(lineages: pd.DataFrame | None, barcodes: pd.Index,
                  origin: pd.Series | None = None, min_second_prob: float = 0.2,
                  cluster_to_major: dict | None = None,
                  pair_source: str = "lineage_top2") -> pd.DataFrame:
    """Per barcode: lineage_pair ("A|B", sorted) and pair_category."""
    if pair_source not in PAIR_SOURCES:
        raise ValueError(f"pair_source={pair_source!r}; use one of {PAIR_SOURCES}")
    out = pd.DataFrame(index=barcodes)
    out["lineage_pair"] = pd.NA
    has_top2 = False
    if lineages is not None and len(lineages):
        lin = lineages.set_index("barcode").reindex(barcodes)
        a = lin["lineage_major"].astype("string")
        b = a.copy()
        has_top2 = "lineage_major_2" in lin and "confidence_major_2" in lin
        if has_top2:
            second = lin["confidence_major_2"].astype(float) >= min_second_prob
            b = b.where(~second, lin["lineage_major_2"].astype("string"))
        out["lineage_pair"] = _pair(a, b)
    use_origin = origin is not None and (pair_source == "scdblfinder_origin" or not has_top2)
    if use_origin:
        o = origin.reindex(barcodes).astype("string")
        has = o.notna() & o.str.contains(r"\+", regex=True)
        if has.any():
            parts = o[has].str.split("+", n=1, expand=True)
            m = cluster_to_major or {}
            pa = parts[0].map(lambda x: m.get(x, x)).astype("string")
            pb = parts[1].map(lambda x: m.get(x, x)).astype("string")
            out.loc[has[has].index, "lineage_pair"] = _pair(pa, pb).values
    pair = out["lineage_pair"].astype("string")
    split = pair.str.split("|", n=1, expand=True) if pair.notna().any() else None
    cat = pd.Series("unknown", index=barcodes, dtype=object)
    if split is not None:
        known = pair.notna().values
        cat.loc[known] = np.where(split[0].values[known] == split[1].values[known],
                                  "same_lineage", "cross_lineage")
    out["pair_category"] = cat
    return out


def _pair(a: pd.Series, b: pd.Series) -> pd.Series:
    both = a.notna() & b.notna()
    res = pd.Series(pd.NA, index=a.index, dtype="string")
    res.loc[both.values] = ["|".join(sorted((x, y))) for x, y in zip(a[both], b[both])]
    return res


def _kappa(x: np.ndarray, y: np.ndarray) -> float:
    n = len(x)
    if n == 0:
        return float("nan")
    po = np.mean(x == y)
    pe = np.mean(x) * np.mean(y) + np.mean(~x) * np.mean(~y)
    return float("nan") if pe == 1 else float((po - pe) / (1 - pe))


def concordance(capture: str, calls: pd.DataFrame, known: pd.DataFrame,
                pairs: pd.DataFrame, d_within: float | None = None,
                ratio_flag: float = 2.0) -> pd.DataFrame:
    """One row per engine for one capture. See the module docstring."""
    kn = known.set_index("barcode")
    rows = []
    wide = calls.pivot(index="barcode", columns="engine", values="doublet_call").astype(bool)
    for engine, g in calls.groupby("engine", sort=True):
        g = g.set_index("barcode")
        called = g["doublet_call"].astype(bool)
        is_known = kn["known_doublet"].reindex(g.index).fillna(False).astype(bool)
        status = kn["demux_status"].reindex(g.index).astype("string") \
            if "demux_status" in kn else pd.Series(pd.NA, index=g.index, dtype="string")
        cat = pairs["pair_category"].reindex(g.index).fillna("unknown")
        r = {"capture_id": capture, "engine": engine, "n_cells": len(g),
             "n_called": int(called.sum()), "n_known": int(is_known.sum()),
             "recall_known": _rate(called & is_known, is_known)}
        for c in CATEGORIES:
            sel = is_known & cat.eq(c)
            r[f"n_known_{c}"] = int(sel.sum())
            r[f"n_called_known_{c}"] = int((called & sel).sum())
            r[f"recall_known_{c}"] = _rate(called & sel, sel)
        singlet = status.eq("singlet").fillna(False)
        r["n_singlets"] = int(singlet.sum())
        r["n_called_singlets"] = int((called & singlet).sum())
        r["d_within"] = d_within if d_within is not None else np.nan
        exp = r["n_singlets"] * d_within if d_within is not None else np.nan
        r["expected_within_donor"] = exp
        r["ratio_called_vs_expected"] = r["n_called_singlets"] / exp if exp and exp > 0 else np.nan
        ratio = r["ratio_called_vs_expected"]
        r["singlet_flag"] = ("NA" if np.isnan(ratio) else
                             "excess" if ratio > ratio_flag else
                             "deficit" if ratio < 1 / ratio_flag else "ok")
        for other in wide.columns:
            if other == engine:
                continue
            x, y = wide[engine].values, wide[other].values
            union = (x | y).sum()
            r[f"jaccard_vs_{other}"] = float((x & y).sum() / union) if union else np.nan
            r[f"kappa_vs_{other}"] = _kappa(x, y)
        rows.append(r)
    return pd.DataFrame(rows)


def _rate(num: pd.Series, den: pd.Series) -> float:
    d = int(den.sum())
    return float(num.sum()) / d if d else np.nan


def cohort_rollup(per_capture: pd.DataFrame, ratio_flag: float = 2.0) -> pd.DataFrame:
    """Sum counts across captures per engine and recompute the rates."""
    count_cols = [c for c in per_capture.columns
                  if c.startswith("n_") or c == "expected_within_donor"]
    agg = per_capture.groupby("engine")[count_cols].sum(min_count=1).reset_index()
    agg.insert(0, "capture_id", "cohort")
    agg["n_captures"] = per_capture.groupby("engine").size().values
    agg["recall_known"] = agg["n_called_known_cross_lineage"].add(
        agg["n_called_known_same_lineage"]).add(agg["n_called_known_unknown"]) / agg["n_known"]
    for c in CATEGORIES:
        agg[f"recall_known_{c}"] = agg[f"n_called_known_{c}"] / agg[f"n_known_{c}"].replace(0, np.nan)
    # Only captures with a d_within expectation enter the singlet ratio.
    has_exp = per_capture["expected_within_donor"].notna()
    called_exp = per_capture[has_exp].groupby("engine")["n_called_singlets"].sum()
    agg["n_called_singlets_with_expected"] = agg["engine"].map(called_exp)
    exp = agg["expected_within_donor"]
    agg["ratio_called_vs_expected"] = agg["n_called_singlets_with_expected"] / exp.replace(0, np.nan)
    agg["singlet_flag"] = np.select(
        [agg["ratio_called_vs_expected"].isna(), agg["ratio_called_vs_expected"] > ratio_flag,
         agg["ratio_called_vs_expected"] < 1 / ratio_flag], ["NA", "excess", "deficit"], "ok")
    for col in [c for c in per_capture.columns if c.startswith(("jaccard_vs_", "kappa_vs_"))]:
        w = {}
        for engine, g in per_capture.groupby("engine"):
            ok = g[col].notna()
            w[engine] = np.average(g.loc[ok, col], weights=g.loc[ok, "n_cells"]) if ok.any() else np.nan
        agg[f"mean_{col}"] = agg["engine"].map(w)  # cell-weighted mean over captures
    return agg
