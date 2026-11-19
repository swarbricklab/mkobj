import pandas as pd
import pytest

from doublet_clusters import cluster_labels


def _lin():
    return pd.DataFrame({
        "lineage_major": ["Epi"] * 6 + ["T"] * 4 + [None],
        "lineage_fine": ["Lum"] * 4 + ["Bas"] * 2 + ["CD4"] * 3 + [None, None],
    })


def test_none_returns_no_clusters():
    lab, s = cluster_labels(_lin(), "none")
    assert lab is None and s["clusters"] == "none"


def test_fine_small_types_fold_into_major():
    lab, s = cluster_labels(_lin(), "lineage_fine", min_cluster_size=2)
    # Bas (2) is kept; the NA-fine T cell folds into major T, a 1-cell remnant that
    # joins the largest kept T fine type (CD4).
    assert lab.tolist() == ["Lum"] * 4 + ["Bas"] * 2 + ["CD4"] * 3 + ["CD4", "unlabelled"]
    lab, s = cluster_labels(_lin(), "lineage_fine", min_cluster_size=3)
    # Bas (2) folds into Epi, a remnant of 2 < 3, which joins Lum.
    assert lab.tolist() == ["Lum"] * 6 + ["CD4"] * 4 + ["unlabelled"]
    assert s == {"clusters": "lineage_fine", "n_clusters": 3, "n_fine_merged": 1}


def test_fine_remnant_large_enough_stays_major():
    lin = pd.DataFrame({"lineage_major": ["Epi"] * 9, "lineage_fine": ["Lum"] * 5 + ["Bas"] * 2 + ["Sec"] * 2})
    lab, s = cluster_labels(lin, "lineage_fine", min_cluster_size=3)
    # Bas and Sec fold into Epi: a remnant of 4 >= 3 stays its own cluster.
    assert lab.tolist() == ["Lum"] * 5 + ["Epi"] * 4
    assert s["n_fine_merged"] == 2


def test_small_whole_lineage_stays_its_own_cluster():
    lin = pd.DataFrame({"lineage_major": ["Epi"] * 6 + ["Plasma"] * 2,
                        "lineage_fine": ["Lum"] * 6 + ["Plasmablast"] * 2})
    lab, _ = cluster_labels(lin, "lineage_fine", min_cluster_size=3)
    # No kept fine type in Plasma to join: the small lineage stays distinct.
    assert lab.tolist() == ["Lum"] * 6 + ["Plasma"] * 2


def test_major():
    lab, s = cluster_labels(_lin(), "lineage_major")
    assert lab.tolist() == ["Epi"] * 6 + ["T"] * 4 + ["unlabelled"]
    assert s["n_clusters"] == 3


def test_bad_mode():
    with pytest.raises(ValueError):
        cluster_labels(_lin(), "leiden")


def test_remnant_ignores_fine_types_of_other_lineages():
    # One B cell carries a stray PVL fine label; the B-cell remnant must not join PVL.
    lin = pd.DataFrame({"lineage_major": ["PVL"] * 5 + ["B"] * 3,
                        "lineage_fine": ["PVL Diff"] * 5 + ["B Mem", "B Naive", "PVL Diff"]})
    lab, _ = cluster_labels(lin, "lineage_fine", min_cluster_size=3)
    assert lab.tolist() == ["PVL Diff"] * 5 + ["B", "B", "PVL Diff"]
