import pandas as pd
import pytest

from doublet_io import COLUMNS, combine_calls, scrublet_calls, strip_capture_prefix, validate_calls


def _obs(cross=False):
    obs = pd.DataFrame(
        {"doublet_score": [0.1, 0.5], "predicted_scrublet": [False, True]},
        index=["CAP_A_AAAC-1", "CAP_A_TTTG-1"],
    )
    if cross:
        obs["predicted_cross_lineage_doublet"] = [True, False]
        obs["doublet_cross_lineage_score"] = [0.2, 0.0]
    return obs


def test_strip_prefix_only_when_present():
    assert strip_capture_prefix(["CAP_A_AAAC-1", "TTTG-1"], "CAP_A") == ["AAAC-1", "TTTG-1"]


def test_scrublet_calls_schema_and_values():
    df = scrublet_calls(_obs(), "CAP_A")
    assert list(df.columns) == COLUMNS
    assert df["barcode"].tolist() == ["AAAC-1", "TTTG-1"]
    assert df["doublet_call"].tolist() == [False, True]
    assert (df["engine"] == "scrublet").all()


def test_cross_lineage_flag_only_when_scored():
    assert "cross_lineage_flag" not in scrublet_calls(_obs(), "CAP_A")
    assert scrublet_calls(_obs(cross=True), "CAP_A")["cross_lineage_flag"].tolist() == [True, False]


def test_combine_two_engines_long_format():
    a = scrublet_calls(_obs(cross=True), "CAP_A")
    b = a[COLUMNS].assign(engine="scdblfinder", doublet_score=[0.9, 0.05])
    df = combine_calls([a, b])
    assert len(df) == 4
    assert df.groupby("engine").size().to_dict() == {"scdblfinder": 2, "scrublet": 2}
    # optional column kept, NA for the engine that lacks it
    assert df.loc[df.engine == "scdblfinder", "cross_lineage_flag"].isna().all()


def test_validate_rejects_duplicates_and_missing():
    a = scrublet_calls(_obs(), "CAP_A")
    with pytest.raises(ValueError, match="duplicate"):
        combine_calls([a, a])
    with pytest.raises(ValueError, match="missing"):
        validate_calls(a.drop(columns="doublet_score"))


def test_write_read_round_trips_scores_exactly(tmp_path):
    from doublet_io import read_calls, write_calls
    obs = _obs(cross=True)
    obs["doublet_score"] = [0.17030567685589518, 1 / 3]
    df = scrublet_calls(obs, "CAP_A")
    path = tmp_path / "calls.csv"
    write_calls(df, path)
    back = read_calls(path)
    assert back["doublet_score"].tolist() == df["doublet_score"].tolist()
    assert back["doublet_call"].tolist() == df["doublet_call"].tolist()
    assert back["cross_lineage_flag"].tolist() == df["cross_lineage_flag"].tolist()


def test_known_from_status_and_posteriors():
    from doublet_io import known_from_posteriors, known_from_status
    obs = pd.DataFrame({"status": ["singlet", "doublet", "unassigned"]},
                       index=["CAP_A_A-1", "CAP_A_B-1", "CAP_A_C-1"])
    k = known_from_status(obs, "CAP_A")
    assert k["barcode"].tolist() == ["A-1", "B-1", "C-1"]
    assert k["known_doublet"].tolist() == [False, True, False]
    post = pd.DataFrame({"barcode": ["A-1", "B-1"], "consensus_status": ["Doublet", "singlet"]})
    assert known_from_posteriors(post, "CAP_A")["known_doublet"].tolist() == [True, False]
    with pytest.raises(ValueError):
        known_from_posteriors(post.drop(columns="consensus_status"), "CAP_A")


def test_attach_known_fills_false_for_unlisted():
    from doublet_io import attach_known
    calls = scrublet_calls(_obs(), "CAP_A")
    known = pd.DataFrame({"capture_id": ["CAP_A"], "barcode": ["TTTG-1"], "known_doublet": [True]})
    out = attach_known(calls, known)
    assert out["known_doublet"].tolist() == [False, True]
    assert list(out.columns[:5]) == COLUMNS


def test_posteriors_with_unresolved_vocabulary_are_rejected():
    from doublet_io import posteriors_resolve_doublets
    # snp_demux cell_assignment.tsv passes through as singlet / NA: no doublets visible
    raw = pd.DataFrame({"barcode": ["A-1", "B-1"], "consensus_status": ["singlet", None]})
    ok, why = posteriors_resolve_doublets(raw)
    assert not ok and "NA" in why
    resolved = pd.DataFrame({"barcode": ["A-1", "B-1", "C-1"],
                             "consensus_status": ["singlet", "doublet", "unassigned"]})
    assert posteriors_resolve_doublets(resolved)[0]
    odd = pd.DataFrame({"barcode": ["A-1"], "consensus_status": ["maybe"]})
    assert not posteriors_resolve_doublets(odd)[0]


def test_legacy_scrublet_obs_restores_columns_in_order():
    from doublet_io import legacy_scrublet_obs
    obs = _obs(cross=True)
    obs["Epi"] = [0.9, 0.1]
    calls = scrublet_calls(obs, "CAP_A", lineage_cols=["Epi"])
    # merge order differs from the CSV order; extra engine rows are ignored
    other = calls[COLUMNS].assign(engine="scdblfinder")
    names = ["CAP_A_TTTG-1", "CAP_A_AAAC-1"]
    out = legacy_scrublet_obs(pd.concat([calls, other]), names, "CAP_A")
    assert list(out.columns) == ["doublet_score", "predicted_doublet", "predicted_scrublet",
                                 "predicted_cross_lineage_doublet", "Epi",
                                 "doublet_cross_lineage_score"]
    assert out["doublet_score"].tolist() == [0.5, 0.1]
    assert out["predicted_doublet"].tolist() == [True, True]   # scrublet | cross-lineage
    assert out["Epi"].tolist() == [0.1, 0.9]


def test_legacy_scrublet_obs_rejects_missing_cells():
    from doublet_io import legacy_scrublet_obs
    calls = scrublet_calls(_obs(), "CAP_A")
    with pytest.raises(ValueError, match="stale"):
        legacy_scrublet_obs(calls, ["CAP_A_GGGG-1"], "CAP_A")
