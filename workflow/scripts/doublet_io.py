"""
workflow/scripts/doublet_io.py
Shared schema and helpers for the persisted per-capture doublet calls (#40).

doublets/<capture_id>.csv is long format: one row per called barcode per engine,
so several engines can coexist (e.g. during the Scrublet -> scDblFinder
transition). Barcodes are the raw Cell Ranger barcodes, without the
'<capture>_' prefix the per-capture objects carry; join on (capture_id, barcode).

Pure functions only (no snakemake globals), so they can be unit-tested.
"""

from __future__ import annotations

import pandas as pd

# Required columns, in order. Engines may add optional columns (see OPTIONAL).
COLUMNS = ["capture_id", "barcode", "engine", "doublet_score", "doublet_call"]
# Optional columns, written after COLUMNS when any engine provides them.
OPTIONAL = ["cross_lineage_flag", "known_doublet"]


def strip_capture_prefix(obs_names, capture: str) -> list[str]:
    """'<capture>_<barcode>' -> '<barcode>'; names without the prefix pass through."""
    prefix = f"{capture}_"
    return [n[len(prefix):] if n.startswith(prefix) else n for n in obs_names]


def scrublet_calls(obs: pd.DataFrame, capture: str) -> pd.DataFrame:
    """Long-format calls for the Scrublet engine from a detect_doublets obs frame.

    Expects 'doublet_score' and 'predicted_scrublet'; adds 'cross_lineage_flag'
    when cross-lineage detection produced a score column.
    """
    out = pd.DataFrame({
        "capture_id": capture,
        "barcode": strip_capture_prefix(obs.index, capture),
        "engine": "scrublet",
        "doublet_score": obs["doublet_score"].astype(float).values,
        "doublet_call": obs["predicted_scrublet"].astype(bool).values,
    })
    if "doublet_cross_lineage_score" in obs.columns:
        out["cross_lineage_flag"] = obs["predicted_cross_lineage_doublet"].astype(bool).values
    return out


def validate_calls(df: pd.DataFrame, source: str = "") -> pd.DataFrame:
    """Check the long-format schema and return the frame with columns ordered."""
    missing = [c for c in COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"doublet calls {source!r} missing columns {missing}")
    dup = df.duplicated(["capture_id", "barcode", "engine"])
    if dup.any():
        raise ValueError(
            f"doublet calls {source!r} have {int(dup.sum())} duplicate "
            f"(capture_id, barcode, engine) rows"
        )
    extra = [c for c in df.columns if c not in COLUMNS]
    ordered = [c for c in OPTIONAL if c in extra] + [c for c in extra if c not in OPTIONAL]
    return df[COLUMNS + ordered]


def combine_calls(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Row-bind per-engine call frames into one validated long-format frame."""
    if not frames:
        raise ValueError("no engine call frames to combine")
    combined = pd.concat([validate_calls(f) for f in frames], ignore_index=True, sort=False)
    return validate_calls(combined, "combined")


def read_calls(path) -> pd.DataFrame:
    """Read a doublets/<capture>.csv with stable dtypes (barcodes/ids as str)."""
    df = pd.read_csv(path, dtype={"capture_id": str, "barcode": str, "engine": str})
    return validate_calls(df, str(path))
