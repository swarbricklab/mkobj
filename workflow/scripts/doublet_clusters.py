"""
workflow/scripts/doublet_clusters.py
Cluster labels for scDblFinder's artificial-doublet simulation (#47).

params.doublets.clusters:
  none           random artificial doublets (scDblFinder default)
  lineage_fine   fine lineage labels from lineages/<capture>.csv; fine types with
                 fewer than min_cluster_size cells (or NA) fall back to their
                 major lineage. If what is folded into a major lineage is still
                 smaller than min_cluster_size, it joins the largest kept fine type
                 of that lineage, so no small remnant clusters are left. A lineage
                 that is small as a whole (e.g. a few plasma cells) stays its own
                 cluster: it is a distinct cell type.
  lineage_major  major lineage labels

Cells with no usable label go to an "unlabelled" cluster. Pure functions only.
"""

from __future__ import annotations

import pandas as pd

CLUSTER_MODES = ("none", "lineage_fine", "lineage_major")
UNLABELLED = "unlabelled"


def cluster_labels(lineages: pd.DataFrame, mode: str, min_cluster_size: int = 50):
    """Return (Series of cluster labels indexed like `lineages`, summary dict).

    summary: clusters (mode), n_clusters, n_fine_merged (fine types folded into
    their major lineage).
    """
    if mode not in CLUSTER_MODES:
        raise ValueError(f"params.doublets.clusters={mode!r}; use one of {CLUSTER_MODES}")
    if mode == "none":
        return None, {"clusters": "none", "n_clusters": pd.NA, "n_fine_merged": pd.NA}

    major = lineages["lineage_major"].astype("string")
    if mode == "lineage_major":
        lab = major.fillna(UNLABELLED)
        merged = 0
    else:
        fine = lineages["lineage_fine"].astype("string")
        sizes = fine.value_counts()
        small = set(sizes[sizes < min_cluster_size].index)
        use_major = fine.isna() | fine.isin(small)
        lab = fine.where(~use_major, major)
        merged = len(small)
        # Second pass: a major-lineage remnant (folded small fine types) that is still
        # small joins the largest kept fine type of the same lineage. Major and fine
        # are predicted separately and occasionally disagree, so a fine type belongs
        # to the lineage most of its cells are called, not to every lineage that one
        # stray cell carries.
        kept = pd.DataFrame({"major": major, "fine": fine})[~use_major & major.notna()]
        largest = pd.Series(dtype=str)
        if len(kept):
            home = kept.groupby("fine", observed=True)["major"].agg(lambda m: m.value_counts().idxmax())
            n_fine = kept["fine"].value_counts()
            ranked = pd.DataFrame({"home": home, "n": n_fine.reindex(home.index)})
            largest = ranked.sort_values("n", ascending=False).groupby("home").head(1)
            largest = pd.Series(largest.index, index=largest["home"])
        remnant = use_major & major.notna()
        rsize = major[remnant].value_counts()
        for m, n in rsize.items():
            if n < min_cluster_size and m in largest.index:
                lab[remnant & (major == m)] = largest[m]
        lab = lab.fillna(UNLABELLED)
    lab = lab.astype(str)
    return lab, {"clusters": mode, "n_clusters": int(lab.nunique()), "n_fine_merged": merged}
