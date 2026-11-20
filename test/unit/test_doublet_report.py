import numpy as np
import pandas as pd
import pytest

from doublet_report import cohort_rollup, concordance, lineage_pairs


def _fixture():
    # 10 barcodes: b0-b3 cross-lineage known doublets, b4-b5 same-lineage known
    # doublets, b6-b9 genotype singlets.
    bc = [f"b{i}" for i in range(10)]
    lineages = pd.DataFrame({
        "barcode": bc,
        "lineage_major": ["Epi", "Epi", "T", "T", "Epi", "T", "Epi", "Epi", "T", "T"],
        "lineage_major_2": ["T", "T", "Epi", "Epi", "T", "Epi", "T", "T", "Epi", "Epi"],
        "confidence_major_2": [0.4, 0.45, 0.3, 0.35, 0.01, 0.02, 0.0, 0.0, 0.0, 0.0],
    })
    known = pd.DataFrame({"barcode": bc, "known_doublet": [True] * 6 + [False] * 4,
                          "demux_status": ["doublet"] * 6 + ["singlet"] * 4})
    # engine A: finds 3/4 cross, 0/2 same, calls 1 singlet; engine B: finds all known
    a = [1, 1, 1, 0, 0, 0, 1, 0, 0, 0]
    b = [1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
    calls = pd.concat([
        pd.DataFrame({"capture_id": "C", "barcode": bc, "engine": "A",
                      "doublet_score": np.array(a, float), "doublet_call": np.array(a, bool)}),
        pd.DataFrame({"capture_id": "C", "barcode": bc, "engine": "B",
                      "doublet_score": np.array(b, float), "doublet_call": np.array(b, bool)}),
    ])
    return bc, lineages, known, calls


def test_recall_by_category():
    bc, lineages, known, calls = _fixture()
    pairs = lineage_pairs(lineages, pd.Index(bc))
    assert pairs["pair_category"].tolist()[:6] == ["cross_lineage"] * 4 + ["same_lineage"] * 2
    rep = concordance("C", calls, known, pairs, d_within=0.25).set_index("engine")
    assert rep.loc["A", "recall_known"] == pytest.approx(3 / 6)
    assert rep.loc["A", "recall_known_cross_lineage"] == pytest.approx(3 / 4)
    assert rep.loc["A", "recall_known_same_lineage"] == 0
    assert rep.loc["B", "recall_known_same_lineage"] == 1
    # singlets: 4 x 0.25 = 1 expected; A calls 1 -> ok, B calls 0 -> deficit
    assert rep.loc["A", "expected_within_donor"] == pytest.approx(1.0)
    assert rep.loc["A", "singlet_flag"] == "ok"
    assert rep.loc["B", "singlet_flag"] == "deficit"
    assert rep.loc["A", "jaccard_vs_B"] == pytest.approx(3 / 7)


def test_origin_overrides_top2_and_maps_clusters():
    bc, lineages, known, calls = _fixture()
    origin = pd.Series(["CD4+CD8"] + [None] * 9, index=bc)
    m = {"CD4": "T", "CD8": "T"}
    # lineage_top2 ignores the origin when the runner-up column exists
    assert lineage_pairs(lineages, pd.Index(bc), origin=origin, cluster_to_major=m,
                         pair_source="lineage_top2").loc["b0", "pair_category"] == "cross_lineage"
    # ... and falls back to it when there is no runner-up (crosswalk labels)
    no_top2 = lineages.drop(columns=["lineage_major_2", "confidence_major_2"])
    assert lineage_pairs(no_top2, pd.Index(bc), origin=origin, cluster_to_major=m,
                         pair_source="lineage_top2").loc["b0", "pair_category"] == "same_lineage"
    # default scdblfinder_origin: origin where present, top-2 for the other barcodes
    pairs = lineage_pairs(lineages, pd.Index(bc), origin=origin, cluster_to_major=m)
    assert pairs.loc["b0", "lineage_pair"] == "T|T"
    assert pairs.loc["b0", "pair_category"] == "same_lineage"
    assert pairs.loc["b1", "pair_category"] == "cross_lineage"
    assert pairs.loc["b0", "pair_source"] == "scdblfinder_origin"
    assert pairs.loc["b1", "pair_source"] == "lineage_top2"
    rep = concordance("C", calls, known, pairs).set_index("engine")
    assert rep.loc["A", "pair_source_used"] == "mixed"


def test_no_lineages_is_unknown_and_cohort_rollup():
    bc, _, known, calls = _fixture()
    pairs = lineage_pairs(None, pd.Index(bc))
    assert set(pairs["pair_category"]) == {"unknown"}
    rep = concordance("C", calls, known, pairs, d_within=None)
    assert rep["singlet_flag"].eq("NA").all()
    two = pd.concat([rep, rep.assign(capture_id="D")])
    co = cohort_rollup(two).set_index("engine")
    assert co.loc["A", "n_known"] == 12 and co.loc["A", "recall_known"] == pytest.approx(0.5)
    assert co.loc["A", "n_captures"] == 2


def test_cohort_ratio_uses_only_captures_with_expectation():
    bc, lineages, known, calls = _fixture()
    pairs = lineage_pairs(lineages, pd.Index(bc))
    with_exp = concordance("C", calls, known, pairs, d_within=0.25)
    without = concordance("D", calls, known, pairs, d_within=None)
    co = cohort_rollup(pd.concat([with_exp, without])).set_index("engine")
    # engine A: 1 singlet call in C (expected 1.0); D's call is excluded
    assert co.loc["A", "n_called_singlets"] == 2
    assert co.loc["A", "ratio_called_vs_expected"] == pytest.approx(1.0)


def test_known_used_in_training_is_flagged():
    bc, lineages, known, calls = _fixture()
    pairs = lineage_pairs(lineages, pd.Index(bc))
    rep = concordance("C", calls, known, pairs, trained_on_known=("B",)).set_index("engine")
    assert bool(rep.loc["B", "known_used_in_training"])
    assert not rep.loc["A", "known_used_in_training"]
    co = cohort_rollup(rep.reset_index()).set_index("engine")
    assert co.loc["B", "known_used_in_training"] and not co.loc["A", "known_used_in_training"]
