#!/usr/bin/env python3
"""
workflow/scripts/subset_deposit.py
Stage 3 of the config-gated keep-list filter (see rules/filter_deposit.smk).

Subset the merged object to the resolved keep-list and write the deposit object
with the same CLASSIC var/obs index encoding + gene-axis assertion the merge
stage uses, so older readers never silently lose the gene axis.

Optional, config-gated size reduction for the deposit (filter.x_dtype /
filter.compression_level): X is written by hand with h5py rather than by
anndata, so that
  * X/data can be cast losslessly to a compact integer dtype (raw counts arrive
    as float32; every value is checked to round-trip exactly before casting),
  * X/indices is written as int32 even when nnz > 2**31 -- scipy forces indices
    and indptr to share one dtype, so a large in-memory matrix always carries
    int64 indices, and anndata writes them as-is,
  * gzip is paired with the HDF5 byte-shuffle filter, which anndata's
    write_h5ad(compression=...) cannot enable and which is most of the gain on
    integer counts.
Only standard HDF5 filters (gzip + shuffle) are used, so any HDF5 reader
(h5py, anndata, rhdf5/zellkonverter, anndataR) can open the file without plugins.
"""

import sys
import logging

import numpy as np
import pandas as pd
import anndata as ad
import h5py
import scipy.sparse as sp

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
flag_low_expression_genes = bool(snakemake.params.get("flag_low_expression_genes", False))
min_cells_per_gene = snakemake.params.get("min_cells_per_gene", 3)
# Deposit size reduction (both default off -> plain anndata write, as before).
x_dtype = snakemake.params.get("x_dtype")                      # e.g. "uint32"
compression_level = snakemake.params.get("compression_level")  # gzip 0-9

if x_dtype is not None:
    x_dtype = np.dtype(x_dtype)
if compression_level is not None:
    compression_level = int(compression_level)
    if not 0 <= compression_level <= 9:
        raise ValueError(f"compression_level must be 0-9 (gzip), got {compression_level}")
dset_kw = (
    dict(compression="gzip", compression_opts=compression_level, shuffle=True)
    if compression_level is not None
    else {}
)
custom_x = x_dtype is not None or compression_level is not None
CHUNK = 50_000_000  # elements per streamed block when checking / writing X

keep_df = pd.read_csv(keep_csv, dtype=str)
keep = keep_df["cell_id"].astype(str).to_numpy()
keep_set = set(keep)
# resolve_keep carries the per-cell `usage` disposition (analysed/control) alongside
# the cell_id so the deposit object can expose it without re-joining.
usage_by_cell = keep_df.set_index("cell_id")["usage"] if "usage" in keep_df.columns else None
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

# Attach the per-cell usage disposition from the keep-list.
if usage_by_cell is not None:
    sub.obs["usage"] = usage_by_cell.reindex(sub.obs_names.astype(str)).astype(str).to_numpy()
    logger.info(f"usage breakdown: {sub.obs['usage'].value_counts().to_dict()}")

# Optional, config-gated gene-level low-expression flagging. Computed HERE, on the
# deposited subset, so var['n_cells'] counts only the cells actually in the deposit
# (not the pre-filter merge). Off by default: the raw deposit keeps the full,
# unannotated Cell Ranger gene axis and leaves the filter decision to consumers.
if flag_low_expression_genes:
    n_cells = np.asarray((sub.X > 0).sum(axis=0)).ravel()
    sub.var["n_cells"] = n_cells
    sub.var["is_filtered"] = n_cells < int(min_cells_per_gene)
    logger.info(
        f"Gene flagging: {int(sub.var['is_filtered'].sum())}/{sub.n_vars} genes "
        f"below min_cells_per_gene={min_cells_per_gene}"
    )

# ── CLASSIC string encoding for the var/obs index (see merge_anndata.py) ─────────
sub.var.index = pd.Index(np.asarray(sub.var.index, dtype=object), name=sub.var.index.name)
sub.obs.index = pd.Index(np.asarray(sub.obs.index, dtype=object), name=sub.obs.index.name)



def _check_lossless_cast(data, dtype):
    """Raise unless every value in `data` survives a cast to `dtype` exactly."""
    if data.dtype == dtype:
        return
    lo, hi = (np.iinfo(dtype).min, np.iinfo(dtype).max) if dtype.kind in "iu" else (None, None)
    for s in range(0, data.size, CHUNK):
        c = data[s : s + CHUNK]
        if lo is not None and (c.min() < lo or c.max() > hi):
            raise ValueError(
                f"X values [{c.min()}, {c.max()}] fall outside the {dtype} range [{lo}, {hi}]; "
                f"choose a wider filter.x_dtype."
            )
        if not np.array_equal(c.astype(dtype).astype(c.dtype), c):
            raise ValueError(
                f"X holds values that do not round-trip through {dtype} (non-integer counts?); "
                f"refusing a lossy cast. Unset filter.x_dtype to keep {data.dtype}."
            )


def _write_sparse_x(path, X, data_dtype, kw):
    """Stream a CSR/CSC matrix into `path`/X with compact dtypes + compression.

    Writes the standard anndata csr_matrix/csc_matrix (0.1.0) encoding, so the
    result reads back with ad.read_h5ad like any other object.
    """
    fmt = X.format
    if fmt not in ("csr", "csc"):
        X = X.tocsr()
        fmt = "csr"
    nnz = int(X.nnz)
    minor = X.shape[1] if fmt == "csr" else X.shape[0]
    idx_dtype = np.int32 if minor <= np.iinfo(np.int32).max else np.int64
    ptr_dtype = np.int32 if nnz <= np.iinfo(np.int32).max else np.int64

    def _dset(group, name, n, dtype):
        # Chunked + compressed when asked; empty arrays cannot be chunked.
        extra = dict(chunks=True, **kw) if (kw and n > 0) else {}
        return group.create_dataset(name, shape=(n,), dtype=dtype, **extra)

    with h5py.File(path, "r+") as hf:
        if "X" in hf:
            del hf["X"]
        g = hf.create_group("X")
        g.attrs["encoding-type"] = f"{fmt}_matrix"
        g.attrs["encoding-version"] = "0.1.0"
        g.attrs["shape"] = np.asarray(X.shape, dtype=np.int64)

        d_ptr = _dset(g, "indptr", X.indptr.size, ptr_dtype)
        d_ptr[...] = X.indptr.astype(ptr_dtype, copy=False)
        d_idx = _dset(g, "indices", nnz, idx_dtype)
        d_dat = _dset(g, "data", nnz, data_dtype)
        for s in range(0, nnz, CHUNK):
            e = min(s + CHUNK, nnz)
            d_idx[s:e] = X.indices[s:e].astype(idx_dtype, copy=False)
            d_dat[s:e] = X.data[s:e].astype(data_dtype, copy=False)
            logger.info(f"  X written {e / nnz * 100:5.1f}% ({e:,}/{nnz:,} nnz)")
    logger.info(
        f"Wrote X ({fmt}): data={np.dtype(data_dtype)} indices={np.dtype(idx_dtype)} "
        f"indptr={np.dtype(ptr_dtype)} "
        f"{'gzip-%d+shuffle' % kw['compression_opts'] if kw else 'uncompressed'}"
    )


logger.info(f"Writing deposit object -> {out_h5ad}")
if not custom_x:
    sub.write_h5ad(out_h5ad)
else:
    X = sub.X
    if not sp.issparse(X):
        raise ValueError("filter.x_dtype / compression_level need a sparse X; got a dense matrix.")
    out_x_dtype = x_dtype if x_dtype is not None else X.dtype
    _check_lossless_cast(X.data, out_x_dtype)
    # Write everything except X through anndata (obs/var/uns, compressed), then
    # add X by hand. X is detached first rather than written-then-replaced:
    # HDF5 never reclaims space from deleted datasets, so overwriting an
    # anndata-written X would leave its full uncompressed size as dead weight.
    sub.X = None
    if compression_level is not None:
        sub.write_h5ad(out_h5ad, compression="gzip", compression_opts=compression_level)
    else:
        sub.write_h5ad(out_h5ad)
    _write_sparse_x(out_h5ad, X, out_x_dtype, dset_kw)
    sub.X = X


def _force_classic_index(path, axis):
    with h5py.File(path, "r+") as hf:
        node = hf[axis]["_index"]
        if not isinstance(node, h5py.Group):
            return
        vals = node["values"][:]
        del hf[axis]["_index"]
        ds = hf[axis].create_dataset(
            "_index", data=vals, dtype=h5py.string_dtype(encoding="utf-8"),
            **(dict(chunks=True, **dset_kw) if (dset_kw and len(vals)) else {}),
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

if custom_x:
    with h5py.File(out_h5ad, "r") as _hf:
        _x = _hf["X"]
        for _name in ("data", "indices", "indptr"):
            _d = _x[_name]
            logger.info(
                f"X/{_name}: dtype={_d.dtype} compression={_d.compression} "
                f"level={_d.compression_opts} shuffle={_d.shuffle}"
            )
        assert _x["data"].dtype == out_x_dtype, f"X/data is {_x['data'].dtype}, expected {out_x_dtype}"
        if compression_level is not None and _x["data"].shape[0] > 0:
            assert _x["data"].compression == "gzip" and _x["data"].shuffle, "X/data not gzip+shuffle"
    # Lossless spot-check: compare a few row windows of the written X with the
    # in-memory matrix (exact, after the same dtype cast).
    assert check.X.shape == X.shape and check.n_obs == sub.n_obs
    rng = np.random.default_rng(0)
    win = min(1000, sub.n_obs)
    starts = {0, max(sub.n_obs - win, 0)} | set(rng.integers(0, max(sub.n_obs - win, 1), size=3).tolist())
    for a in sorted(starts):
        got = sp.csr_matrix(check.X[a : a + win])
        want = sp.csr_matrix(X[a : a + win]).astype(out_x_dtype)
        if got.shape != want.shape or (got != want).nnz != 0:
            raise AssertionError(f"Written X differs from the in-memory matrix in rows {a}:{a + win}")
    logger.info(f"Lossless spot-check passed on {len(starts)} row windows of {win}")
check.file.close()
logger.info("Done subsetting deposit")
