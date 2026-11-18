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


def _float_repr(x) -> str:
    """Shortest string that round-trips to the same float64 ('' for NA)."""
    return "" if pd.isna(x) else repr(float(x))


def write_calls(df: pd.DataFrame, path) -> None:
    """Write calls so float columns round-trip exactly.

    pandas' default to_csv float formatting is not guaranteed to round-trip
    (pandas 3 writes 16 significant digits, losing the last bit of e.g.
    0.17030567685589518), so the persisted scores would not equal the obs values.
    Format floats with Python's shortest round-trip repr instead.
    """
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = out[col].map(_float_repr)
    out.to_csv(path, index=False)


def read_calls(path) -> pd.DataFrame:
    """Read a doublets/<capture>.csv with stable dtypes (barcodes/ids as str)
    and exact float round-trip."""
    df = pd.read_csv(path, dtype={"capture_id": str, "barcode": str, "engine": str},
                     float_precision="round_trip")
    return validate_calls(df, str(path))


# ── Known (genotype-detected) doublets (#42) ─────────────────────────────────────
KNOWN_COLUMNS = ["capture_id", "barcode", "known_doublet", "known_source"]


def known_from_posteriors(post: pd.DataFrame, capture: str) -> pd.DataFrame:
    """Known doublets from snp_demux's doublet_posteriors.csv (swarbricklab/snp_demux#53):
    consensus_status == 'doublet'."""
    missing = {"barcode", "consensus_status"} - set(post.columns)
    if missing:
        raise ValueError(f"doublet_posteriors.csv for {capture} missing {sorted(missing)}")
    return pd.DataFrame({
        "capture_id": capture,
        "barcode": post["barcode"].astype(str).values,
        "known_doublet": post["consensus_status"].astype(str).str.lower().eq("doublet").values,
        "known_source": "doublet_posteriors",
    })


def known_from_status(obs: pd.DataFrame, capture: str) -> pd.DataFrame:
    """Known doublets from the per-capture object's demux `status` column
    (create_anndata's resolution: singlet / doublet / unassigned)."""
    return pd.DataFrame({
        "capture_id": capture,
        "barcode": strip_capture_prefix(obs.index, capture),
        "known_doublet": obs["status"].astype(str).eq("doublet").values,
        "known_source": "assignment_status",
    })


def attach_known(calls: pd.DataFrame, known: pd.DataFrame) -> pd.DataFrame:
    """Left-join known_doublet onto long-format calls by (capture_id, barcode).
    Barcodes absent from `known` get False."""
    k = known[["capture_id", "barcode", "known_doublet"]]
    out = calls.drop(columns=["known_doublet"], errors="ignore").merge(
        k, on=["capture_id", "barcode"], how="left", validate="many_to_one")
    out["known_doublet"] = out["known_doublet"].astype("boolean").fillna(False).astype(bool)
    return validate_calls(out)
