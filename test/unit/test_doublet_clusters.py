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
    lab, s = cluster_labels(_lin(), "lineage_fine", min_cluster_size=3)
    assert lab.tolist() == ["Lum"] * 4 + ["Epi"] * 2 + ["CD4"] * 3 + ["T", "unlabelled"]
    assert s == {"clusters": "lineage_fine", "n_clusters": 5, "n_fine_merged": 1}


def test_major():
    lab, s = cluster_labels(_lin(), "lineage_major")
    assert lab.tolist() == ["Epi"] * 6 + ["T"] * 4 + ["unlabelled"]
    assert s["n_clusters"] == 3


def test_bad_mode():
    with pytest.raises(ValueError):
        cluster_labels(_lin(), "leiden")
