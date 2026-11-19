"""
workflow/scripts/lineage_crosswalk.py
Crosswalk cellranger annotate labels to major lineages (#46).

The crosswalk CSV (deps.lineage_crosswalk) has columns
  cellranger_label, lineage_major
Lines starting with '#' are comments. Labels are matched exactly after
stripping whitespace; unmapped labels become NA and are reported.

Pure functions only, so they can be unit-tested.
"""

from __future__ import annotations

import pandas as pd


def load_crosswalk(path) -> dict[str, str]:
    xw = pd.read_csv(path, comment="#", dtype=str, skipinitialspace=True)
    missing = {"cellranger_label", "lineage_major"} - set(xw.columns)
    if missing:
        raise ValueError(f"crosswalk {path} missing columns {sorted(missing)}")
    xw = xw.dropna(subset=["cellranger_label"])
    xw["cellranger_label"] = xw["cellranger_label"].str.strip()
    dup = xw["cellranger_label"].duplicated(keep=False)
    if dup.any():
        raise ValueError(f"crosswalk {path} maps labels more than once: "
                         f"{sorted(xw.loc[dup, 'cellranger_label'].unique())}")
    return dict(zip(xw["cellranger_label"], xw["lineage_major"].str.strip()))


def apply_crosswalk(labels: pd.Series, crosswalk: dict[str, str]):
    """Map labels -> lineage_major. Returns (mapped Series with NA for unmapped,
    counts of unmapped labels as a Series, most frequent first). Missing input
    labels (NaN) stay NA and are not reported as unmapped."""
    lab = labels.astype("string").str.strip()
    mapped = lab.map(crosswalk).astype("string")
    unmapped = lab[lab.notna() & mapped.isna()].value_counts()
    return mapped, unmapped


def agreement_table(capture: str, model_major: pd.Series, crosswalk_major: pd.Series,
                    collapse: dict[str, str] | None = None) -> pd.DataFrame:
    """Long-format model-vs-crosswalk agreement for one capture.

    collapse maps model classes onto crosswalk names before comparing (e.g.
    Cancer/Normal Epithelial -> Epithelial when the crosswalk cannot tell them
    apart). Rows: capture_id, lineage_major_model, lineage_major_crosswalk,
    n_cells, agree.
    """
    m = model_major.astype("string").replace(collapse or {})
    x = crosswalk_major.astype("string")
    df = pd.DataFrame({"lineage_major_model": m.fillna("NA"),
                       "lineage_major_crosswalk": x.fillna("NA")})
    tab = df.value_counts().rename("n_cells").reset_index()
    tab.insert(0, "capture_id", capture)
    tab["agree"] = tab["lineage_major_model"] == tab["lineage_major_crosswalk"]
    return tab.sort_values(["lineage_major_model", "n_cells"], ascending=[True, False],
                           ignore_index=True)
