"""
workflow/scripts/doublet_clusters.py
Cluster labels for scDblFinder's artificial-doublet simulation (#47).

params.doublets.clusters:
  none           random artificial doublets (scDblFinder default)
  lineage_fine   fine lineage labels from lineages/<capture>.csv; fine types with
                 fewer than min_cluster_size cells (or NA) fall back to their
                 major lineage, so there are no tiny clusters
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
        lab = fine.where(~use_major, major).fillna(UNLABELLED)
        merged = len(small)
    lab = lab.astype(str)
    return lab, {"clusters": mode, "n_clusters": int(lab.nunique()), "n_fine_merged": merged}
