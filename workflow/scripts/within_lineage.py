"""
workflow/scripts/within_lineage.py
Within-lineage doublet score from QC features (#52; exploratory).

Same-lineage within-donor doublets (e.g. T + T from one donor) are invisible to
genotypes and nearly invisible to expression engines. Cross-donor same-lineage
known doublets are confirmed examples of the same physical event, so per major
lineage we learn how they differ from singlets of that lineage in QC features
(UMIs, genes, mito/ribo fractions), normalised within capture.

Per lineage:
  labelled   cross-donor known doublets whose lineage pair is L|L
  unlabelled genotype singlets labelled L. These are NOT clean negatives: they
             hide the within-donor same-lineage doublets we want to count.
This is positive-unlabelled (PU) learning. Following Elkan & Noto (2008):
  g(x) = P(labelled | x), fitted by logistic regression (out-of-fold)
  c    = mean g over held-out labelled cells = P(labelled | doublet)
  each unlabelled cell's probability of being a doublet:
       w(x) = (1 - c) / c * g(x) / (1 - g(x)), clipped to [0, 1]
  predicted within-donor same-lineage doublets = sum of w over a capture's
  unlabelled cells of L.
The key assumption (SCAR) is that genotype detection selects doublets
independently of their QC features: cross-donor and within-donor same-lineage
doublets are the same physical event. A threshold-based count corrected with
the false-positive rate on the same singlets is NOT used: it collapses to ~0
by construction because the hidden doublets inflate that rate.
AUC (labelled vs unlabelled, cross-validated) measures separability and is
conservative, since the unlabelled set contains positives.

Pure functions only (scikit-learn for fitting).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

FEATURES = ["log_total_counts", "log_n_genes", "pct_counts_mt", "pct_counts_ribo"]


def qc_features(X, gene_names) -> pd.DataFrame:
    """Raw per-cell QC features from a cells x genes count matrix."""
    from scipy import sparse
    X = sparse.csr_matrix(X)
    names = pd.Index([str(g).upper() for g in gene_names])
    total = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    mt = np.asarray(X[:, np.where(names.str.startswith("MT-"))[0]].sum(axis=1)).ravel()
    ribo = np.asarray(X[:, np.where(names.str.startswith(("RPS", "RPL")))[0]].sum(axis=1)).ravel()
    denom = np.where(total > 0, total, 1)
    return pd.DataFrame({"log_total_counts": np.log1p(total), "log_n_genes": np.log1p(n_genes),
                         "pct_counts_mt": 100 * mt / denom, "pct_counts_ribo": 100 * ribo / denom})


def normalise_within(df: pd.DataFrame, by: str, cols=FEATURES) -> pd.DataFrame:
    """z-score each feature within groups of `by` (e.g. capture)."""
    out = df.copy()
    g = df.groupby(by)[cols]
    out[cols] = (df[cols] - g.transform("mean")) / g.transform("std").replace(0, 1)
    return out


def fit_lineage(cells: pd.DataFrame, lineage: str, n_splits: int = 5, seed: int = 0,
                min_pos: int = 10):
    """Fit and cross-validate one lineage (Elkan-Noto PU). `cells` has FEATURES
    (normalised), lineage_major, pair_category, known_doublet, demux_status,
    capture_id.

    Returns (separability row, unlabelled cells of the lineage with g and w).
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold

    lin = cells["lineage_major"].astype(str) == lineage
    pos = lin & cells["known_doublet"].astype(bool) & cells["pair_category"].eq("same_lineage")
    unl = lin & cells["demux_status"].astype(str).eq("singlet")
    row = {"lineage": lineage, "n_labelled": int(pos.sum()), "n_unlabelled": int(unl.sum()),
           "auc_cv": np.nan, "c_hat": np.nan, "est_prevalence_unlabelled": np.nan}
    if pos.sum() < min_pos or unl.sum() < min_pos:
        row["note"] = f"fewer than {min_pos} labelled or unlabelled cells"
        return row, None

    train = cells[pos | unl]
    X = train[FEATURES].values
    y = pos[pos | unl].values.astype(int)
    oof = np.zeros(len(y))
    k = int(min(n_splits, y.sum(), (1 - y).sum()))
    for tr, te in StratifiedKFold(n_splits=k, shuffle=True, random_state=seed).split(X, y):
        m = LogisticRegression(max_iter=1000).fit(X[tr], y[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
    c = float(oof[y == 1].mean())
    g = np.clip(oof[y == 0], 0, 1 - 1e-9)
    w = np.clip((1 - c) / c * g / (1 - g), 0, 1)
    model = LogisticRegression(max_iter=1000).fit(X, y)
    row.update(auc_cv=float(roc_auc_score(y, oof)), c_hat=c,
               est_prevalence_unlabelled=float(w.mean()), note="",
               **{f"coef_{f}": float(v) for f, v in zip(FEATURES, model.coef_[0])})
    unlabelled = train[y == 0].copy()
    unlabelled["g"] = g
    unlabelled["w"] = w
    return row, unlabelled


def per_capture_estimates(scored: pd.DataFrame, expected: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per capture x lineage: unlabelled singlets and the PU estimate of
    within-donor same-lineage doublets among them (sum of w), plus the expected
    count when `expected` is given (capture_id, lineage,
    expected_within_same_lineage)."""
    if scored is None or not len(scored):
        return pd.DataFrame(columns=["capture_id", "lineage", "n_singlets",
                                     "predicted_within_same_lineage"])
    out = (scored.groupby(["capture_id", "lineage_major"])
           .agg(n_singlets=("w", "size"), predicted_within_same_lineage=("w", "sum"))
           .reset_index().rename(columns={"lineage_major": "lineage"}))
    if expected is not None and len(out):
        out = out.merge(expected, on=["capture_id", "lineage"], how="left")
        out["ratio_predicted_vs_expected"] = (out["predicted_within_same_lineage"]
                                              / out["expected_within_same_lineage"])
    return out
