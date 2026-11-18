#!/usr/bin/env python3
"""
tools/train_lineage_classifier.py
Train the versioned per-capture lineage classifier (#44) once, offline.

Input: a reference AnnData with raw counts and two label columns, e.g. the
Wu et al. 2021 breast cancer mini-atlas (celltype_major / celltype_minor).
Output (--out): model.npz + model.json (see workflow/scripts/lineage_model.py),
metrics.csv (held-out accuracy per level), confusion_<level>.csv,
gene_space.csv (panel vs. target matrix gene space) and a README.md stub.

Typical use (on a normal node, no GPU):

  python tools/train_lineage_classifier.py \\
      --reference wu2021_miniatlas.h5ad --major-key celltype_major \\
      --fine-key celltype_minor \\
      --target-features refdata-gex-GRCh38-2020-A/features.tsv.gz \\
      --version wu2021-v1 --out resources/lineage_classifier/wu2021-v1
  dvc add resources/lineage_classifier/wu2021-v1

Then point params.lineages.model at the directory.
"""

import argparse
import sys
from datetime import date
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
from lineage_model import align_genes, load_model, predict, save_model, train  # noqa: E402


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--reference", required=True, help="reference .h5ad with raw counts")
    ap.add_argument("--counts-layer", default=None, help="layer holding raw counts (default: X)")
    ap.add_argument("--major-key", default="celltype_major")
    ap.add_argument("--fine-key", default="celltype_minor")
    ap.add_argument("--symbol-key", default=None, help="var column with gene symbols (default: var_names)")
    ap.add_argument("--gene-id-key", default=None, help="var column with Ensembl ids (optional)")
    ap.add_argument("--target-features", default=None,
                    help="features.tsv.gz of the target Cell Ranger reference; restricts the "
                         "panel to genes present there and reports gene-space compatibility")
    ap.add_argument("--n-genes", type=int, default=3000, help="highly variable genes in the panel")
    ap.add_argument("--n-markers", type=int, default=50,
                    help="top marker genes per fine class added to the panel (HVG selection "
                         "alone can miss strong, high-mean markers)")
    ap.add_argument("--min-cells-per-class", type=int, default=20,
                    help="fine classes with fewer cells are folded into their major class")
    ap.add_argument("--test-frac", type=float, default=0.2)
    ap.add_argument("--C", type=float, default=0.1, help="inverse L2 strength")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-refit", action="store_true",
                    help="ship the train-split model instead of refitting on all cells")
    ap.add_argument("--version", required=True)
    ap.add_argument("--out", required=True)
    return ap.parse_args()


def main():
    a = parse_args()
    out = Path(a.out)
    ref = ad.read_h5ad(a.reference)
    X = ref.layers[a.counts_layer] if a.counts_layer else ref.X
    X = sparse.csr_matrix(X)
    if X.data.size and not np.allclose(X.data[:10000], np.round(X.data[:10000])):
        sys.exit("reference matrix is not raw counts; pass --counts-layer")
    symbols = ref.var[a.symbol_key].astype(str).values if a.symbol_key else ref.var_names.astype(str).values
    ids = ref.var[a.gene_id_key].astype(str).values if a.gene_id_key else np.array([""] * ref.n_vars)

    major = ref.obs[a.major_key].astype(str).values
    fine = ref.obs[a.fine_key].astype(str).values
    counts = pd.Series(fine).value_counts()
    small = counts[counts < a.min_cells_per_class].index
    fine = np.where(pd.Series(fine).isin(small), major, fine)
    print(f"{ref.n_obs} cells; {len(set(major))} major, {len(set(fine))} fine classes "
          f"({len(small)} small fine classes folded into major)")

    # ── Gene space: restrict to genes present in the target reference ───────────
    keep = np.ones(ref.n_vars, dtype=bool)
    gene_space = {"reference_genes": ref.n_vars}
    if a.target_features:
        tf = pd.read_csv(a.target_features, sep="\t", header=None,
                         names=["gene_id", "gene_symbol", "feature_type"])
        tf = tf[tf.feature_type == "Gene Expression"]
        t_ids = set(tf.gene_id.str.split(".").str[0])
        t_sym = set(tf.gene_symbol)
        by_id = pd.Series(ids).str.split(".").str[0].isin(t_ids).values
        by_sym = pd.Series(symbols).isin(t_sym).values
        keep = by_id | by_sym
        gene_space.update(target_genes=len(tf), matched_by_id=int(by_id.sum()),
                          matched_by_symbol_only=int((by_sym & ~by_id).sum()),
                          unmatched=int((~keep).sum()))
        print(f"gene space: {gene_space}")

    # ── HVG panel on the shared genes ───────────────────────────────────────────
    tmp = ad.AnnData(X=X[:, keep].astype(np.float32))
    tmp.var_names = pd.Index(symbols[keep]).astype(str)
    tmp.var_names_make_unique()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    sc.pp.highly_variable_genes(tmp, n_top_genes=min(a.n_genes, tmp.n_vars), flavor="cell_ranger")
    in_panel = tmp.var["highly_variable"].values.copy()
    if a.n_markers:
        tmp.obs["fine"] = pd.Categorical(fine)
        sc.tl.rank_genes_groups(tmp, "fine", method="t-test", n_genes=a.n_markers)
        markers = set(pd.DataFrame(tmp.uns["rank_genes_groups"]["names"]).values.ravel())
        in_panel |= tmp.var_names.isin(markers)
        print(f"panel: {int(tmp.var['highly_variable'].sum())} HVGs + markers -> {int(in_panel.sum())} genes")
    panel_idx = np.where(keep)[0][in_panel]
    genes = pd.DataFrame({"gene_symbol": symbols[panel_idx], "gene_id": ids[panel_idx]})
    Xp = X[:, panel_idx]
    gene_space["panel_genes"] = len(panel_idx)

    # ── Held-out evaluation (stratified by fine label) ──────────────────────────
    rng = np.random.default_rng(a.seed)
    test = np.zeros(ref.n_obs, dtype=bool)
    for lab in np.unique(fine):
        idx = np.where(fine == lab)[0]
        n_test = int(round(len(idx) * a.test_frac))
        test[rng.choice(idx, size=n_test, replace=False)] = True
    labels = {"major": major, "fine": fine}
    meta = dict(reference=Path(a.reference).name, major_key=a.major_key, fine_key=a.fine_key,
                n_cells=int(ref.n_obs), C=a.C, seed=a.seed, trained=str(date.today()),
                gene_space=gene_space)
    arrays, info = train(Xp[~test], genes, {k: v[~test] for k, v in labels.items()},
                         C=a.C, seed=a.seed, meta={**meta, "split": "train"})
    save_model(out / "_heldout", arrays, info, version=f"{a.version}-heldout")
    pred, _ = predict(load_model(out / "_heldout"), Xp[test], genes.gene_symbol, genes.gene_id)

    rows = []
    for level in ("major", "fine"):
        y, p = labels[level][test], pred[f"lineage_{level}"].values
        per_class = pd.crosstab(pd.Series(y, name="true"), pd.Series(p, name="predicted"))
        per_class.to_csv(out / f"confusion_{level}.csv")
        recall = np.diag(per_class.reindex(index=per_class.index, columns=per_class.index,
                                           fill_value=0)) / per_class.sum(axis=1)
        rows.append(dict(level=level, n_test=int(test.sum()), accuracy=float((y == p).mean()),
                         balanced_accuracy=float(recall.mean()),
                         median_confidence=float(np.median(pred[f"confidence_{level}"]))))
    metrics = pd.DataFrame(rows)
    metrics.to_csv(out / "metrics.csv", index=False)
    pd.Series(gene_space).to_csv(out / "gene_space.csv", header=["value"])
    print(metrics.to_string(index=False))

    # ── Shipped model ───────────────────────────────────────────────────────────
    if not a.no_refit:
        arrays, info = train(Xp, genes, labels, C=a.C, seed=a.seed, meta={**meta, "split": "all"})
    info["meta"]["heldout_metrics"] = rows
    save_model(out, arrays, info, version=a.version)
    (out / "README.md").write_text(f"""# Lineage classifier {a.version}

Multinomial logistic regression (see `workflow/scripts/lineage_model.py`), trained
{date.today()} with `tools/train_lineage_classifier.py`.

- Reference: `{Path(a.reference).name}`, {ref.n_obs} cells; labels `{a.major_key}` (major)
  and `{a.fine_key}` (fine; classes with < {a.min_cells_per_class} cells folded into major).
- Gene panel: {len(panel_idx)} highly variable genes (`genes` in `model.json`).
- Gene space: {gene_space}
- Held-out accuracy ({a.test_frac:.0%} stratified split): see `metrics.csv`.
- {"Shipped model refit on all cells." if not a.no_refit else "Shipped model is the train-split model."}
- `model.npz` sha256 is in `model.json` and checked on load.
""")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
