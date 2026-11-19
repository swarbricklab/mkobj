import numpy as np
import pandas as pd
import pytest

from doublet_decomposition import capture_summary, ols, pair_enrichment


def test_capture_summary_excess_and_groups():
    bc = [f"b{i}" for i in range(8)]
    calls = pd.DataFrame({"capture_id": "C", "barcode": bc, "engine": "e",
                          "doublet_score": 0.0,
                          "doublet_call": [1, 1, 0, 1, 1, 0, 0, 0]}).astype({"doublet_call": bool})
    known = pd.DataFrame({"barcode": bc, "known_doublet": [1, 1, 0, 0, 0, 0, 0, 0],
                          "demux_status": ["doublet"] * 2 + ["singlet"] * 6}).astype({"known_doublet": bool})
    pairs = pd.DataFrame({"lineage_pair": ["A|B", "A|A", "A|A", "A|C", "A|C", "A|A", "B|B", "B|B"]},
                         index=bc)
    row, grp = capture_summary("C", calls, known, pairs, "e", d_within=0.25, f_x=0.2)
    assert row["n_singlets"] == 6 and row["within_donor_calls"] == 2
    assert row["expected_within_donor"] == pytest.approx(1.5)
    assert row["excess_within_donor"] == pytest.approx(0.5)
    assert grp.groupby("group").size().to_dict() == {"cross_donor_known": 2, "within_donor_call": 2}


def test_pair_enrichment_flags_adjacent():
    groups = pd.DataFrame({
        "group": ["within_donor_call"] * 4 + ["cross_donor_known"] * 4,
        "lineage_pair": ["A|C"] * 3 + ["A|A"] + ["A|C", "A|A", "A|B", "A|B"],
    })
    out = pair_enrichment(groups, adjacent=["C|A"]).set_index("lineage_pair")
    assert out.loc["A|C", "enrichment"] > 1 > out.loc["A|B", "enrichment"]
    assert out.loc["A|C", "adjacent"] and not out.loc["A|A", "adjacent"]
    assert out.loc["A|A", "same_lineage"]


def test_ols_recovers_coefficients():
    rng = np.random.default_rng(0)
    n = 200
    df = pd.DataFrame({"cells": rng.uniform(1, 20, n),
                       "protocol": rng.choice(["co_dissociated", "separate_then_pooled"], n)})
    df["f_x"] = 0.01 + 0.004 * df.cells + 0.02 * (df.protocol == "separate_then_pooled") \
        + rng.normal(0, 0.002, n)
    r = ols(df, "f_x", ["cells"], ["protocol"]).set_index("term")
    assert r.loc["cells", "estimate"] == pytest.approx(0.004, abs=5e-4)
    assert r.loc["protocol[separate_then_pooled]", "estimate"] == pytest.approx(0.02, abs=2e-3)
    assert r.loc["cells", "p"] < 1e-6
