import pandas as pd

from attach_obs import doublet_columns, lineage_columns, read_table, somatic_columns


def test_lineage_columns_renamed_and_joined():
    lin = pd.DataFrame({"capture_id": "C", "barcode": ["a", "b"], "lineage_major": ["Epi", "T"],
                        "lineage_fine": ["Lum", "CD4"], "confidence_major": [0.9, 0.8],
                        "confidence_fine": [0.7, 0.6], "label_source": "model",
                        "lineage_major_2": ["T", "Epi"]})
    out, n = lineage_columns(lin, ["b", "a", "z"])
    assert n == 2
    assert list(out.columns) == ["lineage_major", "lineage_fine", "lineage_confidence_major",
                                 "lineage_confidence_fine", "lineage_source"]
    assert out["lineage_major"].tolist()[:2] == ["T", "Epi"] and pd.isna(out.loc["z", "lineage_major"])


def test_doublet_columns_wide_per_engine():
    calls = pd.DataFrame({"capture_id": "C", "barcode": ["a", "b", "a", "b"],
                          "engine": ["scrublet", "scrublet", "scdblfinder", "scdblfinder"],
                          "doublet_score": [0.1, 0.9, 0.2, 0.8],
                          "doublet_call": [False, True, False, True],
                          "known_doublet": [False, True, False, True]})
    out, n = doublet_columns(calls, ["a", "b"])
    assert n == 2
    assert set(out.columns) == {"doublet_score_scrublet", "doublet_call_scrublet",
                                "doublet_score_scdblfinder", "doublet_call_scdblfinder",
                                "known_doublet"}
    assert out["doublet_call_scdblfinder"].tolist() == [False, True]


def test_somatic_csv_then_tsv_fallback(tmp_path):
    (tmp_path / "cell_calls.tsv").write_text(
        "barcode\tdonor\tevaluable\tburden_matched\tmalignant\n"
        "a\t01\t1\t2\t1\nb\t02\t0\t\t\n")
    calls = read_table(tmp_path / "cell_calls.csv", tmp_path / "cell_calls.tsv")
    clones = pd.DataFrame({"barcode": ["a"], "donor": ["01"], "p_tumour": [0.95]})
    out, n = somatic_columns(calls, clones, ["a", "b", "c"])
    assert n == 2
    assert out.loc["a", "somatic_malignant"] == 1 and pd.isna(out.loc["b", "somatic_burden_matched"])
    assert out.loc["a", "somatic_p_tumour"] == 0.95 and pd.isna(out.loc["b", "somatic_p_tumour"])
    assert read_table(tmp_path / "nope.csv") is None


def test_somatic_r_na_tokens_are_missing_and_numeric(tmp_path):
    # cardelino's R write.table writes NA as "NA"; somatifind's calls use "".
    from attach_obs import read_table, somatic_columns
    calls = tmp_path / "cell_calls.tsv"
    calls.write_text("barcode\tdonor\tevaluable\tburden_matched\tmalignant\n"
                     "A-1\t4825\t1\t3\t1\nB-1\t\t0\t0\t\n")
    clones = tmp_path / "cell_clone_calls.tsv"
    clones.write_text("barcode\tdonor\tp_tumour\tmalignant\nA-1\t4825\t0.9\t1\nB-1\tNA\tNA\t\n")
    c, k = read_table(tmp_path / "x.csv", calls), read_table(tmp_path / "y.csv", clones)
    out, n = somatic_columns(c, k, ["A-1", "B-1", "C-1"])
    assert n == 2
    assert out["somatic_p_tumour"].dtype.kind == "f"
    assert out["somatic_p_tumour"].tolist()[0] == 0.9
    assert out["somatic_p_tumour"].isna().tolist() == [False, True, True]
    assert out["somatic_malignant"].isna().tolist() == [False, True, True]
    assert all(out[col].dtype.kind in "fi" for col in out.columns)


def test_doublet_calls_plain_bool_without_gaps_nullable_with_gaps():
    import pandas as pd
    from attach_obs import doublet_columns
    calls = pd.DataFrame({"barcode": ["A-1", "B-1"], "engine": "scrublet",
                          "doublet_score": [0.1, 0.9], "doublet_call": [False, True],
                          "known_doublet": [False, True]})
    full, _ = doublet_columns(calls, ["A-1", "B-1"])
    assert full["doublet_call_scrublet"].dtype == bool and full["known_doublet"].dtype == bool
    gap, _ = doublet_columns(calls, ["A-1", "B-1", "C-1"])
    assert str(gap["doublet_call_scrublet"].dtype) == "boolean"
    assert gap["doublet_call_scrublet"].isna().tolist() == [False, False, True]
