"""
workflow/scripts/attach_obs.py
Per-cell annotation columns attached to the merged objects (#50).

Each builder takes the per-capture source table(s) and the capture's raw
barcodes and returns a DataFrame indexed like those barcodes, plus the number
of barcodes matched (join coverage). Pure functions only.

  lineages   lineages/<capture>.csv  -> lineage_major, lineage_fine,
             lineage_confidence_major, lineage_confidence_fine, lineage_source
  doublets   doublets/<capture>.csv  -> doublet_score_<engine>,
             doublet_call_<engine>, known_doublet
  somatifind calls/<capture>/cell_calls.csv (or .tsv) -> somatic_evaluable,
             somatic_burden_matched, somatic_malignant; cardelino
             cell_clone_calls.csv (or .tsv) -> somatic_p_tumour
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

LINEAGE_MAP = {
    "lineage_major": "lineage_major",
    "lineage_fine": "lineage_fine",
    "confidence_major": "lineage_confidence_major",
    "confidence_fine": "lineage_confidence_fine",
    "label_source": "lineage_source",
}
SOMATIC_MAP = {
    "evaluable": "somatic_evaluable",
    "burden_matched": "somatic_burden_matched",
    "malignant": "somatic_malignant",
}


def _join(table: pd.DataFrame, barcodes, mapping: dict) -> tuple[pd.DataFrame, int]:
    t = table.drop_duplicates("barcode").set_index("barcode")
    cols = {src: dst for src, dst in mapping.items() if src in t.columns}
    out = t[list(cols)].rename(columns=cols).reindex(pd.Index(barcodes))
    return out, int(pd.Index(barcodes).isin(t.index).sum())


def lineage_columns(lineages: pd.DataFrame, barcodes):
    return _join(lineages, barcodes, LINEAGE_MAP)


def doublet_columns(calls: pd.DataFrame, barcodes):
    """Wide per-engine score/call columns plus known_doublet."""
    idx = pd.Index(barcodes)
    out = pd.DataFrame(index=idx)
    for engine, g in calls.groupby("engine", sort=True):
        g = g.drop_duplicates("barcode").set_index("barcode")
        out[f"doublet_score_{engine}"] = g["doublet_score"].reindex(idx).astype(float).values
        out[f"doublet_call_{engine}"] = g["doublet_call"].reindex(idx).astype("boolean").values
    if "known_doublet" in calls.columns:
        k = calls.drop_duplicates("barcode").set_index("barcode")["known_doublet"]
        out["known_doublet"] = k.reindex(idx).astype("boolean").values
    # Plain bool when every cell has a value (the usual case): nullable "boolean"
    # is only needed for gaps and is less portable (Seurat conversion, validators).
    for c in out.columns:
        if out[c].dtype == "boolean" and not out[c].isna().any():
            out[c] = out[c].astype(bool)
    return out, int(idx.isin(calls["barcode"]).sum())


def read_table(path_csv: Path, path_tsv: Path | None = None) -> pd.DataFrame | None:
    """Read somatifind's CSV, falling back to the legacy TSV. None if neither exists."""
    for path, sep in ((path_csv, ","), (path_tsv, "\t")):
        if path is not None and Path(path).exists():
            # Missing values: "" (somatifind's Python writer) and "NA" (cardelino's
            # R write.table). Donor ids stay strings.
            return pd.read_csv(path, sep=sep, dtype={"barcode": str, "donor": str},
                               keep_default_na=False, na_values=["", "NA"])
    return None


def somatic_columns(calls: pd.DataFrame | None, clones: pd.DataFrame | None, barcodes):
    idx = pd.Index(barcodes)
    out = pd.DataFrame(index=idx)
    n = 0
    if calls is not None:
        cols, n = _join(calls, barcodes, SOMATIC_MAP)
        out = out.join(cols)
    if clones is not None and "p_tumour" in clones.columns:
        cols, n_c = _join(clones, barcodes, {"p_tumour": "somatic_p_tumour"})
        out = out.join(cols)
        n = max(n, n_c)
    # All somatic columns are numeric (flags, counts, a probability). Coerce strictly,
    # so an unexpected token fails here instead of becoming a categorical column.
    for c in out.columns:
        out[c] = pd.to_numeric(out[c], errors="raise")
    return out, n
