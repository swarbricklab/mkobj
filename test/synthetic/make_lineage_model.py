#!/usr/bin/env python3
"""
test/synthetic/make_lineage_model.py
Train a toy lineage classifier on the synthetic fixture (CAP_A singlets) with
tools/train_lineage_classifier.py, for running annotate_lineages locally.
Run make_fixture.py first. Writes test/synthetic/data/lineage_model/.
"""
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io

here = Path(__file__).resolve().parent
data = here / "data"
cr = data / "cellranger" / "CAP_A"
X = io.mmread(cr / "matrix.mtx.gz").T.tocsr()
feats = pd.read_csv(cr / "features.tsv.gz", sep="\t", header=None, names=["id", "sym", "type"])
bcs = pd.read_csv(cr / "barcodes.tsv.gz", header=None)[0]
truth = pd.read_csv(data / "truth" / "CAP_A.csv").set_index("barcode").loc[bcs]
single = ~truth["is_doublet"].values
ref = ad.AnnData(X=X[single].astype(np.float32),
                 obs=truth.loc[single, ["major", "fine"]].reset_index(drop=True).astype(str))
ref.obs_names = bcs[single].values
ref.var_names = feats["sym"].values
ref.var["gene_ids"] = feats["id"].values
ref_path = data / "lineage_reference.h5ad"
ref.write_h5ad(ref_path)
subprocess.run([sys.executable, str(here.parents[1] / "tools" / "train_lineage_classifier.py"),
                "--reference", str(ref_path), "--major-key", "major", "--fine-key", "fine",
                "--gene-id-key", "gene_ids", "--target-features", str(cr / "features.tsv.gz"),
                "--n-genes", "200", "--n-markers", "20", "--version", "synthetic-v1",
                "--out", str(data / "lineage_model")], check=True)
