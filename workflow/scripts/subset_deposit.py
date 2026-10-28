#!/usr/bin/env python3
"""
workflow/scripts/subset_deposit.py
Stage 3 of the config-gated keep-list filter (see rules/filter_deposit.smk).

Subset the merged object to the resolved keep-list and write the deposit object
with the same CLASSIC var/obs index encoding + gene-axis assertion the merge
stage uses, so older readers never silently lose the gene axis.
"""

import sys
import logging

import numpy as np
import pandas as pd
import anndata as ad
import h5py

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)

in_h5ad = snakemake.input.h5ad
keep_csv = snakemake.input.keep
out_h5ad = snakemake.output.deposit
expected_n_vars = snakemake.params.get("expected_n_vars")

keep = pd.read_csv(keep_csv)["cell_id"].astype(str).to_numpy()
keep_set = set(keep)
logger.info(f"Keep-list: {len(keep_set)} unique cell_ids")

logger.info(f"Reading {in_h5ad} (full, in-memory)...")
adata = ad.read_h5ad(in_h5ad)
logger.info(f"  loaded n_obs={adata.n_obs} n_vars={adata.n_vars}")

mask = adata.obs_names.isin(keep_set)
n_match = int(mask.sum())
if n_match != len(keep_set):
    missing = len(keep_set) - n_match
    raise ValueError(
        f"{missing} keep-list cell_ids not found in the object "
        f"(matched {n_match}/{len(keep_set)}). Keep-list and object are out of sync."
    )

sub = adata[mask].copy()
del adata
logger.info(f"Subset: n_obs={sub.n_obs} n_vars={sub.n_vars}")

# ── CLASSIC string encoding for the var/obs index (see merge_anndata.py) ─────────
sub.var.index = pd.Index(np.asarray(sub.var.index, dtype=object), name=sub.var.index.name)
sub.obs.index = pd.Index(np.asarray(sub.obs.index, dtype=object), name=sub.obs.index.name)

logger.info(f"Writing deposit object -> {out_h5ad}")
sub.write_h5ad(out_h5ad)


def _force_classic_index(path, axis):
    with h5py.File(path, "r+") as hf:
        node = hf[axis]["_index"]
        if not isinstance(node, h5py.Group):
            return
        vals = node["values"][:]
        del hf[axis]["_index"]
        ds = hf[axis].create_dataset(
            "_index", data=vals, dtype=h5py.string_dtype(encoding="utf-8")
        )
        ds.attrs["encoding-type"] = "string-array"
        ds.attrs["encoding-version"] = "0.2.0"
        logger.info(f"Rewrote {axis}/_index nullable-string-array -> classic string-array")


for _axis in ("var", "obs"):
    with h5py.File(out_h5ad, "r") as _hf:
        _rewrite = dict(_hf[_axis]["_index"].attrs).get("encoding-type") != "string-array"
    if _rewrite:
        _force_classic_index(out_h5ad, _axis)

# ── Post-write verification ──────────────────────────────────────────────────────
check = ad.read_h5ad(out_h5ad, backed="r")
with h5py.File(out_h5ad, "r") as _hf:
    _var_enc = dict(_hf["var"]["_index"].attrs).get("encoding-type")
logger.info(
    f"Re-read check: n_obs={check.n_obs} n_vars={check.n_vars} "
    f"var0={check.var_names[0]!r} var_index_enc={_var_enc}"
)
if expected_n_vars is not None:
    assert _var_enc == "string-array", f"var index encoding is {_var_enc}, not classic string-array."
    assert check.n_vars == int(expected_n_vars), (
        f"GENE AXIS LOST: re-read n_vars={check.n_vars}, expected {expected_n_vars}."
    )
    assert str(check.var_names[0]).startswith("ENSG"), (
        f"First var_name {check.var_names[0]!r} is not an Ensembl id."
    )
    logger.info(f"Gene-axis assertion passed: n_vars={check.n_vars}")
check.file.close()
logger.info("Done subsetting deposit")
