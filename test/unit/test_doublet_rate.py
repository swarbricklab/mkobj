import pytest

from doublet_rate import PARAM_COLUMNS, choose_dbr, per1k_for_chip


def test_multi_donor_capture_with_summary_uses_genotype():
    r = choose_dbr("C1", 5000, "genotype", chip="standard",
                   summary={"n_samples": 4, "d_total": 0.08, "d_total_sd": 0.01})
    assert r["dbr_source_used"] == "genotype"
    assert r["dbr_used"] == pytest.approx(0.08)
    assert r["dbr_sd_used"] == pytest.approx(0.01)
    assert r["pass_dbr"] and r["pass_dbr_sd"]
    assert set(r) == set(PARAM_COLUMNS)


def test_single_donor_capture_falls_back_to_loading():
    r = choose_dbr("C2", 5000, "genotype", chip="standard",
                   summary={"n_samples": 1, "d_total": "NA", "d_total_sd": "NA"})
    assert r["dbr_source_used"] == "loading"
    assert r["fallback_reason"] == "single-sample capture"
    assert r["dbr_used"] == pytest.approx(0.04)          # 0.008 * 5
    assert not r["pass_dbr"]


@pytest.mark.parametrize("summary,reason", [
    (None, "no doublet_summary.csv"),
    ({"n_samples": 3, "d_total": float("nan")}, "d_total is NA"),
])
def test_genotype_fallback_reasons(summary, reason):
    r = choose_dbr("C", 1000, "genotype", summary=summary)
    assert (r["dbr_source_used"], r["fallback_reason"]) == ("loading", reason)


def test_ht_chip_halves_per1k():
    r = choose_dbr("C3", 10000, "loading", chip="HT")
    assert r["dbr_per1k"] == pytest.approx(0.004)
    assert r["dbr_used"] == pytest.approx(0.04)
    assert r["dbr_sd_used"] == pytest.approx(0.016)


def test_none_disables_expectation():
    r = choose_dbr("C", 1000, "none")
    assert r["dbr_sd_used"] == 1.0 and r["pass_dbr_sd"] and not r["pass_dbr"]


def test_genotype_missing_sd_uses_package_default_fraction():
    r = choose_dbr("C", 1000, "genotype", summary={"n_samples": 2, "d_total": 0.05})
    assert r["dbr_sd_used"] == pytest.approx(0.02) and not r["pass_dbr_sd"]


def test_unknown_chip_and_source_rejected():
    assert per1k_for_chip(None) == 0.008 and per1k_for_chip("nan") == 0.008
    with pytest.raises(ValueError):
        per1k_for_chip("v4")
    with pytest.raises(ValueError):
        choose_dbr("C", 1, "magic")
