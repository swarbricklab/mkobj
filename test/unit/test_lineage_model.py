import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from lineage_model import align_genes, load_model, predict, save_model, train


def _toy(rng, n=200):
    genes = pd.DataFrame({"gene_symbol": ["EPCAM", "KRT5", "PTPRC", "CD8A", "G1", "G2"],
                          "gene_id": ["E1", "E2", "E3", "E4", "E5", "E6"]})
    fine = rng.choice(["Lum", "Bas", "CD4", "CD8"], size=n)
    major = np.where(np.isin(fine, ["Lum", "Bas"]), "Epi", "T")
    mu = np.full((n, 6), 2.0)
    mu[major == "Epi", 0] += 30
    mu[fine == "Bas", 1] += 30
    mu[major == "T", 2] += 30
    mu[fine == "CD8", 3] += 30
    return sparse.csr_matrix(rng.poisson(mu)), genes, {"major": major, "fine": fine}


def test_train_save_load_predict_roundtrip(tmp_path):
    rng = np.random.default_rng(0)
    X, genes, labels = _toy(rng)
    arrays, info = train(X, genes, labels)
    save_model(tmp_path, arrays, info, version="t1")
    m = load_model(tmp_path)
    assert m["info"]["version"] == "t1"
    assert m["info"]["fine_to_major"]["Bas"] == "Epi"
    # Query with shuffled columns, symbols only, one extra gene
    order = [3, 1, 0, 5, 2, 4]
    Xq = sparse.hstack([X[:, order], sparse.csr_matrix(np.ones((X.shape[0], 1)))]).tocsr()
    syms = list(genes.gene_symbol.values[order]) + ["EXTRA"]
    pred, n = predict(m, Xq, syms)
    assert n == 6
    assert (pred["lineage_major"].values == labels["major"]).mean() > 0.98
    assert (pred["lineage_fine"].values == labels["fine"]).mean() > 0.95
    assert pred.filter(like="confidence").apply(lambda c: c.between(0, 1).all()).all()


def test_align_genes_prefers_id_and_zero_fills_missing():
    genes = pd.DataFrame({"gene_symbol": ["A", "B", "C"], "gene_id": ["ENSG1", "ENSG2", "ENSG3"]})
    X = sparse.csr_matrix(np.array([[1, 2], [3, 4]]))
    # column 0 has symbol 'B' but id ENSG1.5 -> matched to A by (version-stripped) id
    out, n = align_genes(X, ["B", "ZZZ"], ["ENSG1.5", "ENSG9"], genes)
    assert n == 1
    assert out.toarray().tolist() == [[1, 0, 0], [3, 0, 0]]


def test_load_rejects_tampered_npz(tmp_path):
    rng = np.random.default_rng(1)
    X, genes, labels = _toy(rng, 80)
    arrays, info = train(X, genes, labels)
    save_model(tmp_path, arrays, info, version="t")
    np.savez_compressed(tmp_path / "model.npz", **{**arrays, "mean": arrays["mean"] + 1})
    with pytest.raises(ValueError, match="sha256"):
        load_model(tmp_path)
