"""
workflow/scripts/doublet_decomposition.py
Decompose doublets into random co-encapsulation vs aggregates (#51; exploratory).
See docs/doublet_decomposition.md for the method and its assumptions.

Pure functions only, so they can be unit-tested.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def capture_summary(capture: str, calls: pd.DataFrame, known: pd.DataFrame,
                    pairs: pd.DataFrame, engine: str, d_within: float | None,
                    f_x: float | None) -> tuple[dict, pd.DataFrame]:
    """Per-capture excess within-donor calls and per-barcode pair labels.

    Returns (summary row, frame of barcodes in two groups with their lineage pair):
      group = within_donor_call   engine call among genotype singlets
      group = cross_donor_known   genotype (cross-donor) doublet: the random-pairing
                                  template for lineage-pair composition
    """
    g = calls[calls["engine"] == engine].set_index("barcode")
    kn = known.set_index("barcode").reindex(g.index)
    singlet = kn["demux_status"].astype("string").eq("singlet").fillna(False).values
    is_known = kn["known_doublet"].astype("boolean").fillna(False).values.astype(bool)
    called = g["doublet_call"].astype(bool).values
    n_singlets = int(singlet.sum())
    n_calls = int((called & singlet).sum())
    expected = n_singlets * d_within if d_within is not None else np.nan
    row = dict(capture_id=capture, engine=engine, n_cells=len(g), n_singlets=n_singlets,
               n_known=int(is_known.sum()),
               f_x=f_x if f_x is not None else is_known.mean() if len(g) else np.nan,
               d_within=d_within if d_within is not None else np.nan,
               within_donor_calls=n_calls, expected_within_donor=expected,
               excess_within_donor=n_calls - expected,
               excess_rate=(n_calls - expected) / n_singlets if n_singlets else np.nan)
    p = pairs.reindex(g.index)
    grp = pd.concat([
        pd.DataFrame({"group": "within_donor_call", "lineage_pair": p["lineage_pair"][called & singlet]}),
        pd.DataFrame({"group": "cross_donor_known", "lineage_pair": p["lineage_pair"][is_known]}),
    ])
    grp.insert(0, "capture_id", capture)
    return row, grp


def pair_enrichment(groups: pd.DataFrame, adjacent: list[str] | None = None,
                    pseudocount: float = 0.5) -> pd.DataFrame:
    """Lineage-pair composition of within-donor calls vs the cross-donor template.

    enrichment = share among within-donor calls / share among cross-donor known
    doublets (pseudocount-smoothed). Pairs are 'A|B' sorted. `adjacent` marks the
    spatially adjacent pairs of interest.
    """
    g = groups.dropna(subset=["lineage_pair"])
    tab = pd.crosstab(g["lineage_pair"], g["group"])
    for col in ("within_donor_call", "cross_donor_known"):
        if col not in tab:
            tab[col] = 0
    k = len(tab)
    share_w = (tab["within_donor_call"] + pseudocount) / (tab["within_donor_call"].sum() + pseudocount * k)
    share_x = (tab["cross_donor_known"] + pseudocount) / (tab["cross_donor_known"].sum() + pseudocount * k)
    out = pd.DataFrame({
        "n_within_donor_calls": tab["within_donor_call"],
        "n_cross_donor_known": tab["cross_donor_known"],
        "share_within_donor": share_w, "share_cross_donor": share_x,
        "enrichment": share_w / share_x,
    }).rename_axis("lineage_pair").reset_index()
    parts = out["lineage_pair"].str.split("|", n=1, expand=True)
    out["same_lineage"] = parts[0] == parts[1]
    adj = {"|".join(sorted(p.split("|"))) for p in (adjacent or [])}
    out["adjacent"] = out["lineage_pair"].isin(adj)
    return out.sort_values("enrichment", ascending=False, ignore_index=True)


def ols(df: pd.DataFrame, y: str, numeric: list[str], factors: list[str]) -> pd.DataFrame:
    """Ordinary least squares with treatment-coded factors; returns
    term, estimate, se, t, p, n. Rows with NA in any used column are dropped."""
    d = df.dropna(subset=[y] + numeric + factors)
    X = [np.ones(len(d))]
    names = ["intercept"]
    for c in numeric:
        X.append(d[c].astype(float).values)
        names.append(c)
    for c in factors:
        levels = sorted(d[c].astype(str).unique())
        for lev in levels[1:]:
            X.append((d[c].astype(str) == lev).astype(float).values)
            names.append(f"{c}[{lev}]")
    X = np.column_stack(X)
    yv = d[y].astype(float).values
    n, k = X.shape
    if n <= k:
        return pd.DataFrame({"term": names, "estimate": np.nan, "se": np.nan, "t": np.nan,
                             "p": np.nan, "n": n})
    beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
    resid = yv - X @ beta
    sigma2 = resid @ resid / (n - k)
    se = np.sqrt(np.diag(sigma2 * np.linalg.pinv(X.T @ X)))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), df=n - k)
    return pd.DataFrame({"term": names, "estimate": beta, "se": se, "t": t, "p": p, "n": n})
