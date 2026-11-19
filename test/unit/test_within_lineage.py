import numpy as np
import pandas as pd
import pytest

from within_lineage import FEATURES, fit_lineage, normalise_within, per_capture_estimates, qc_features


def test_qc_features():
    X = np.array([[5, 3, 2, 0], [0, 0, 0, 0]])
    f = qc_features(X, ["MT-CO1", "RPS3", "GENE", "RPL13"])
    assert f.loc[0, "pct_counts_mt"] == pytest.approx(50) and f.loc[0, "pct_counts_ribo"] == pytest.approx(30)
    assert f.loc[1, "log_total_counts"] == 0


def _cells(rng, n_singlet=1500, n_known=150, n_hidden=0, shift=2.0):
    """Two captures; known same-lineage doublets, plus `n_hidden` doublets per
    capture hidden among the genotype singlets (the PU setting)."""
    rows = []
    for cap in ("C1", "C2"):
        offset = rng.normal(0, 2, len(FEATURES))   # capture batch effect
        kinds = ["singlet"] * n_singlet + ["known"] * n_known + ["hidden"] * n_hidden
        for kind in kinds:
            dbl = kind != "singlet"
            f = rng.normal(0, 1, len(FEATURES)) + offset + (shift if dbl else 0) * np.array([1, 1, 0, 0])
            rows.append({**dict(zip(FEATURES, f)), "capture_id": cap, "lineage_major": "T",
                         "hidden": kind == "hidden", "known_doublet": kind == "known",
                         "pair_category": "same_lineage" if kind == "known" else "unknown",
                         "demux_status": "doublet" if kind == "known" else "singlet"})
    return normalise_within(pd.DataFrame(rows), "capture_id")


def test_pu_estimate_tracks_hidden_doublets():
    """Elkan-Noto has an upward floor when classes overlap (documented), but the
    estimate moves ~1:1 with the number of hidden doublets."""
    def est(hidden, seed):
        _, scored = fit_lineage(_cells(np.random.default_rng(seed), n_hidden=hidden), "T")
        return per_capture_estimates(scored)["predicted_within_same_lineage"].mean()
    base = np.mean([est(0, s) for s in range(3)])
    more = np.mean([est(60, s) for s in range(3)])
    assert 40 < more - base < 80          # truth: +60 per capture
    assert base < 60                      # floor at AUC ~0.97


def test_too_few_positives():
    rng = np.random.default_rng(1)
    row, scored = fit_lineage(_cells(rng, n_singlet=200, n_known=3), "T")
    assert scored is None and np.isnan(row["auc_cv"])
