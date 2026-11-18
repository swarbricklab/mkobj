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
