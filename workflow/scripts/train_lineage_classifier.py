#!/usr/bin/env python3
"""
workflow/scripts/train_lineage_classifier.py
Snakemake wrapper for the train_lineage_classifier rule: turns
params.lineages.train into lineage_train.py arguments and trains the model
into lineage_classifier/<version>/ under the results directory.
"""
import sys

from lineage_train import main

cfg = snakemake.params.train
out_dir = snakemake.params.out_dir

argv = ["--reference", snakemake.input.reference,
        "--major-key", cfg["major_key"], "--fine-key", cfg["fine_key"],
        "--version", cfg["version"], "--out", out_dir]
if snakemake.input.get("target_features"):
    argv += ["--target-features", snakemake.input.target_features]
for key in ("counts_layer", "symbol_key", "gene_id_key", "group_key",
            "n_genes", "n_markers", "min_cells_per_class", "test_frac", "C", "seed"):
    if cfg.get(key) is not None:
        argv += ["--" + key.replace("_", "-"), str(cfg[key])]
if cfg.get("use_raw"):
    argv.append("--use-raw")

with open(snakemake.log[0], "w") as log:
    sys.stdout = sys.stderr = log
    print("lineage_train.py " + " ".join(argv), flush=True)
    main(argv)
