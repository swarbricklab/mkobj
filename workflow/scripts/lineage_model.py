"""
workflow/scripts/lineage_model.py
A lightweight, versioned per-cell lineage classifier (#44).

Model: multinomial logistic regression on a fixed gene panel, one model per
label level ("major", "fine"). Cells are normalised to 10k counts, log1p'd and
z-scaled with the reference means/SDs (clipped at +-10), as in CellTypist.

The trained model is stored framework-free so prediction needs only numpy and
never depends on the scikit-learn version used for training:

  <model_dir>/model.npz   per level: coef_<level> (classes x genes),
                          intercept_<level>; shared: mean, scale
  <model_dir>/model.json  version, genes (symbol + Ensembl id), classes per
                          level, fine->major map, training metadata, and
                          the sha256 of model.npz (checked on load)

Training uses scikit-learn (shipped with scanpy) and is done once, offline
(tools/train_lineage_classifier.py). Prediction runs in
the existing scanpy container. Pure functions only (no snakemake globals).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

LEVELS = ("major", "fine")
TARGET_SUM = 1e4
CLIP = 10.0


# ── Normalisation ────────────────────────────────────────────────────────────────
def log_normalise(X):
    """Counts -> log1p(CP10k), returned dense float32 (cells x genes)."""
    X = sparse.csr_matrix(X, dtype=np.float32)
    totals = np.asarray(X.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    X = sparse.diags(TARGET_SUM / totals) @ X
    X.data = np.log1p(X.data)
    return np.asarray(X.todense(), dtype=np.float32)


def align_genes(X, var_symbols, var_ids, genes: pd.DataFrame):
    """Reorder/expand the columns of X to the model gene panel.

    Matches on Ensembl id first (version-stripped), then on gene symbol. Genes
    missing from X become zero columns. Returns (X_panel, n_matched).
    """
    strip = lambda s: pd.Index([str(v).split(".")[0] for v in s])
    pos_id = {g: i for i, g in enumerate(strip(var_ids))} if var_ids is not None else {}
    pos_sym = {}
    for i, g in enumerate(var_symbols):
        pos_sym.setdefault(str(g), i)
    panel_ids = strip(genes["gene_id"].fillna(""))
    cols = np.array([pos_id.get(gid, -1) for gid in panel_ids])
    # Symbol fallback only for panel genes without an id match, and never onto a
    # column already claimed by an id match (no column feeds two panel genes).
    claimed = set(cols[cols >= 0])
    for k, sym in enumerate(genes["gene_symbol"]):
        if cols[k] < 0:
            j = pos_sym.get(str(sym), -1)
            if j >= 0 and j not in claimed:
                cols[k] = j
                claimed.add(j)
    hit = cols >= 0
    X = sparse.csr_matrix(X, dtype=np.float32)
    # Take column 0 as a placeholder for missing genes, then zero those columns.
    out = X[:, np.where(hit, cols, 0)] @ sparse.diags(hit.astype(np.float32))
    return sparse.csr_matrix(out), int(hit.sum())


# ── Training ─────────────────────────────────────────────────────────────────────
def train(X_counts, genes: pd.DataFrame, labels: dict[str, np.ndarray], *,
          C: float = 0.1, max_iter: int = 1000, seed: int = 0, meta: dict | None = None):
    """Fit one multinomial logistic regression per level.

    X_counts: cells x panel-genes counts (already in `genes` order).
    labels:   {"major": array, "fine": array}, one label per cell.
    Returns (arrays, info) for save_model.
    """
    from sklearn.linear_model import LogisticRegression

    Xn = log_normalise(X_counts)
    mean = Xn.mean(axis=0)
    scale = Xn.std(axis=0)
    scale[scale == 0] = 1.0
    Xs = np.clip((Xn - mean) / scale, -CLIP, CLIP)

    arrays = {"mean": mean.astype(np.float32), "scale": scale.astype(np.float32)}
    info = {"levels": {}, "meta": meta or {}}
    for level, y in labels.items():
        clf = LogisticRegression(C=C, max_iter=max_iter, random_state=seed)
        clf.fit(Xs, np.asarray(y).astype(str))
        coef, intercept = clf.coef_, clf.intercept_
        if coef.shape[0] == 1:  # binary: sklearn stores one row; expand to two classes
            coef = np.vstack([-coef, coef]) / 2
            intercept = np.array([-intercept[0], intercept[0]]) / 2
        arrays[f"coef_{level}"] = coef.astype(np.float32)
        arrays[f"intercept_{level}"] = intercept.astype(np.float32)
        info["levels"][level] = {"classes": [str(c) for c in clf.classes_]}
    if {"major", "fine"} <= set(labels):
        pairs = pd.DataFrame({"fine": np.asarray(labels["fine"]).astype(str),
                              "major": np.asarray(labels["major"]).astype(str)})
        info["fine_to_major"] = (pairs.groupby("fine")["major"]
                                 .agg(lambda s: s.value_counts().index[0]).to_dict())
    info["genes"] = genes[["gene_symbol", "gene_id"]].fillna("").to_dict(orient="list")
    return arrays, info


def save_model(model_dir, arrays: dict, info: dict, version: str) -> Path:
    model_dir = Path(model_dir)
    model_dir.mkdir(parents=True, exist_ok=True)
    npz = model_dir / "model.npz"
    np.savez_compressed(npz, **arrays)
    info = {**info, "version": version, "format": "mkobj-lineage-logreg/1",
            "npz_sha256": hashlib.sha256(npz.read_bytes()).hexdigest()}
    (model_dir / "model.json").write_text(json.dumps(info, indent=2) + "\n")
    return model_dir


# ── Prediction ───────────────────────────────────────────────────────────────────
def load_model(model_dir) -> dict:
    model_dir = Path(model_dir)
    info = json.loads((model_dir / "model.json").read_text())
    npz = model_dir / "model.npz"
    sha = hashlib.sha256(npz.read_bytes()).hexdigest()
    if sha != info.get("npz_sha256"):
        raise ValueError(f"{npz}: sha256 {sha} does not match model.json "
                         f"({info.get('npz_sha256')}); model files are inconsistent")
    with np.load(npz) as z:
        arrays = {k: z[k] for k in z.files}
    info["genes"] = pd.DataFrame(info["genes"])
    return {"info": info, **arrays}


def _softmax(z):
    z = z - z.max(axis=1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)


def predict(model: dict, X_counts, var_symbols, var_ids=None, chunk: int = 20000):
    """Predict labels + confidences for every cell.

    Returns (DataFrame with lineage_<level> and confidence_<level> per level,
    n_genes_matched). Confidence is the softmax probability of the chosen class.
    """
    info = model["info"]
    Xp, n_matched = align_genes(X_counts, var_symbols, var_ids, info["genes"])
    out = {f"{p}_{lvl}": [] for lvl in info["levels"] for p in ("lineage", "confidence")}
    for start in range(0, Xp.shape[0], chunk):
        Xs = np.clip((log_normalise(Xp[start:start + chunk]) - model["mean"]) / model["scale"],
                     -CLIP, CLIP)
        for level, lv in info["levels"].items():
            prob = _softmax(Xs @ model[f"coef_{level}"].T + model[f"intercept_{level}"])
            k = prob.argmax(axis=1)
            out[f"lineage_{level}"].append(np.asarray(lv["classes"])[k])
            out[f"confidence_{level}"].append(prob[np.arange(len(k)), k])
    df = pd.DataFrame({k: np.concatenate(v) if v else np.array([]) for k, v in out.items()})
    return df, n_matched
