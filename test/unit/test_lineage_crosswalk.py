import pandas as pd
import pytest

from lineage_crosswalk import agreement_table, apply_crosswalk, load_crosswalk


def test_known_label_maps_and_unmapped_reported(tmp_path):
    p = tmp_path / "xw.csv"
    p.write_text("# reviewed: no\ncellranger_label,lineage_major\nT cell,T-cells\n"
                 " epithelial cell , Epithelial\n")
    xw = load_crosswalk(p)
    assert xw == {"T cell": "T-cells", "epithelial cell": "Epithelial"}
    labels = pd.Series(["T cell", "epithelial cell", "mystery", "mystery", None, "odd"])
    mapped, unmapped = apply_crosswalk(labels, xw)
    assert mapped.tolist()[:2] == ["T-cells", "Epithelial"]
    assert mapped.iloc[2:].isna().all()
    assert unmapped.to_dict() == {"mystery": 2, "odd": 1}


def test_duplicate_labels_rejected(tmp_path):
    p = tmp_path / "xw.csv"
    p.write_text("cellranger_label,lineage_major\nT cell,T-cells\nT cell,B-cells\n")
    with pytest.raises(ValueError, match="more than once"):
        load_crosswalk(p)


def test_agreement_with_collapse():
    model = pd.Series(["Cancer Epithelial", "Normal Epithelial", "T-cells", "T-cells"])
    xw = pd.Series(["Epithelial", "Epithelial", "T-cells", None])
    tab = agreement_table("C", model, xw, collapse={"Cancer Epithelial": "Epithelial",
                                                    "Normal Epithelial": "Epithelial"})
    got = {(r.lineage_major_model, r.lineage_major_crosswalk): (r.n_cells, r.agree)
           for r in tab.itertuples()}
    assert got == {("Epithelial", "Epithelial"): (2, True), ("T-cells", "T-cells"): (1, True),
                   ("T-cells", "NA"): (1, False)}
