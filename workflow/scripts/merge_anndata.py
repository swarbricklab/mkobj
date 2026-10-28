#!/usr/bin/env python3
"""
workflow/scripts/merge_anndata.py
Merge individual capture AnnData objects into a single object.
Optionally filters cells to a subset defined by samples.csv.

Parallels the logic in merge_captures.R.
"""

import sys
import logging
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path

# Set up logging — write to snakemake log file AND stderr
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(snakemake.log[0]),
        logging.StreamHandler(sys.stderr),
    ],
)
logger = logging.getLogger(__name__)

# Ensure uncaught exceptions are logged to the log file
def _excepthook(exc_type, exc_value, exc_tb):
    logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_tb))
    sys.__excepthook__(exc_type, exc_value, exc_tb)

sys.excepthook = _excepthook

# Get input and output paths
h5ad_files = snakemake.input.h5ad_files
output_file = snakemake.output.merged
samples_file = snakemake.params.get("samples_file", "")

# QC filtering parameters
min_counts = snakemake.params.get("min_counts", None)
max_counts = snakemake.params.get("max_counts", None)
min_features = snakemake.params.get("min_features", None)
max_features = snakemake.params.get("max_features", None)
max_pct_mt = snakemake.params.get("max_pct_mt", None)
max_pct_ribo = snakemake.params.get("max_pct_ribo", None)
filter_doublets = snakemake.params.get("filter_doublets", True)
min_cells_per_gene = snakemake.params.get("min_cells_per_gene", 3)
# Optional post-write guard: when set (e.g. the Cell Ranger reference build's
# 36601), the object is re-read with a standard anndata and the gene axis is
# asserted — the backstop against the var encoding silently dropping all genes.
expected_n_vars = snakemake.params.get("expected_n_vars", None)

logger.info(f"Merging {len(h5ad_files)} capture objects...")
logger.info(f"Output file: {output_file}")


# ── Load samples filter (if any) ────────────────────────────────────────────────
samples_filter = None
if samples_file and Path(samples_file).exists():
    logger.info(f"Loading samples filter from: {samples_file}")
    samples_filter = pd.read_csv(samples_file)
    logger.info(f"  Filtering to {len(samples_filter)} sample-capture pairs")


def filter_adata(adata: ad.AnnData, samples_filter: pd.DataFrame | None) -> ad.AnnData:
    """Filter an AnnData object to cells belonging to the subset."""
    if samples_filter is None:
        return adata

    capture = adata.obs['capture'].iloc[0]
    valid_samples = samples_filter.loc[
        samples_filter['capture_id'] == capture, 'sample_id'
    ].unique().tolist()

    if not valid_samples:
        logger.info(f"  No samples listed for capture {capture} — keeping all cells")
        return adata

    n_before = adata.n_obs

    is_cohort = adata.obs['sample_id'].isin(valid_samples)
    is_doublet = adata.obs['status'].eq('doublet') & adata.obs['status'].notna()
    is_unassigned = adata.obs['sample_id'].isna()

    keep = is_cohort | is_doublet | is_unassigned
    adata = adata[keep].copy()

    logger.info(
        f"  Filtered {capture}: {n_before} -> {adata.n_obs} cells "
        f"(kept {is_cohort.sum()} cohort singlets, "
        f"{is_doublet.sum()} doublets, "
        f"{is_unassigned.sum()} unassigned)"
    )

    # Relabel onto the project's declared sample_id, and record the gussid.
    #
    # The samples sheet is a crosswalk: one row per (capture, obs_label), where
    # obs_label is the value create_anndata wrote into obs['sample_id'] — the
    # gussid for a pooled capture (from cell_assignment.tsv) or the capture name
    # for a non-pooled one (mkobj's fallback). Selection above matched on it.
    # Here we set obs['gussid'] and overwrite obs['sample_id'] with the declared
    # name the project uses as its merge key. Cells with no crosswalk match
    # (doublets, unassigned) keep their sentinel/NA sample_id and get a null
    # gussid. Guarded: an old-format sheet (sample_id + capture_id only) still
    # selects, it just does not relabel.
    xwalk = {"gussid", "declared_sample_id", "donor_id", "tissue_id"}
    if xwalk.issubset(samples_filter.columns):
        sub = samples_filter.loc[samples_filter["capture_id"] == capture]
        by_label = sub.drop_duplicates("sample_id").set_index("sample_id")
        obs_label = adata.obs["sample_id"].astype(object)
        adata.obs["gussid"] = obs_label.map(by_label["gussid"]).values
        adata.obs["donor_id"] = obs_label.map(by_label["donor_id"]).values
        adata.obs["tissue_id"] = obs_label.map(by_label["tissue_id"]).values
        declared = obs_label.map(by_label["declared_sample_id"])
        adata.obs["sample_id"] = declared.where(declared.notna(), obs_label).values
        logger.info(
            f"  Relabelled {capture} onto declared sample_id; set obs['gussid']"
        )

    return adata


# ── Load and (optionally) filter objects ─────────────────────────────────────────
objects = []
for i, h5ad_file in enumerate(h5ad_files):
    logger.info(f"Loading: {h5ad_file}")
    adata = ad.read_h5ad(h5ad_file)
    adata = filter_adata(adata, samples_filter)
    if adata.n_obs > 0:
        objects.append(adata)
    else:
        logger.info("  Skipping — no cells remaining after filter")

if not objects:
    raise ValueError("No cells remaining after filtering — cannot create merged object")

# Merge objects
if len(objects) == 1:
    logger.info("Only one capture - no merging needed")
    merged = objects[0]
else:
    logger.info("Merging captures...")
    # Concatenate along obs axis (cells)
    merged = ad.concat(
        objects,
        join='outer',  # Keep all genes from all objects
        merge='same',  # Only merge identical elements
        label='batch',
        keys=[Path(f).stem for f in h5ad_files],
        index_unique=None  # Barcodes already prefixed with capture
    )
    
    # Handle AB modality if present in obsm
    # ad.concat handles obsm automatically, but we need to ensure alignment
    ab_present = all('AB' in obj.obsm for obj in objects)
    if ab_present:
        logger.info("AB modality present in all captures")
        # Check if all AB features are the same
        ab_features = [list(obj.uns.get('AB_features', [])) for obj in objects]
        if all(f == ab_features[0] for f in ab_features):
            merged.uns['AB_features'] = ab_features[0]
        else:
            logger.warning("AB features differ across captures - storing only feature names")
            # Store unique union of AB features
            all_ab_features = []
            for f in ab_features:
                all_ab_features.extend(f)
            merged.uns['AB_features'] = list(dict.fromkeys(all_ab_features))

logger.info(f"Merged object has {merged.n_obs} cells total")

# ── QC cell filtering ────────────────────────────────────────────────────────────
n_before = merged.n_obs
keep = np.ones(merged.n_obs, dtype=bool)

if filter_doublets and 'predicted_doublet' in merged.obs.columns:
    is_doublet = merged.obs['predicted_doublet'].astype(bool)
    n_doublet = is_doublet.sum()
    keep &= ~is_doublet
    logger.info(f"  Doublet filter: flagged {n_doublet} cells")

if min_counts is not None and 'nCount_RNA' in merged.obs.columns:
    fail = merged.obs['nCount_RNA'] < min_counts
    logger.info(f"  min_counts ({min_counts}): flagged {fail.sum()} cells")
    keep &= ~fail

if max_counts is not None and 'nCount_RNA' in merged.obs.columns:
    fail = merged.obs['nCount_RNA'] > max_counts
    logger.info(f"  max_counts ({max_counts}): flagged {fail.sum()} cells")
    keep &= ~fail

if min_features is not None and 'nFeature_RNA' in merged.obs.columns:
    fail = merged.obs['nFeature_RNA'] < min_features
    logger.info(f"  min_features ({min_features}): flagged {fail.sum()} cells")
    keep &= ~fail

if max_features is not None and 'nFeature_RNA' in merged.obs.columns:
    fail = merged.obs['nFeature_RNA'] > max_features
    logger.info(f"  max_features ({max_features}): flagged {fail.sum()} cells")
    keep &= ~fail

if max_pct_mt is not None and 'percent.mt' in merged.obs.columns:
    fail = merged.obs['percent.mt'] > max_pct_mt
    logger.info(f"  max_pct_mt ({max_pct_mt}): flagged {fail.sum()} cells")
    keep &= ~fail

if max_pct_ribo is not None and 'percent.ribo' in merged.obs.columns:
    fail = merged.obs['percent.ribo'] > max_pct_ribo
    logger.info(f"  max_pct_ribo ({max_pct_ribo}): flagged {fail.sum()} cells")
    keep &= ~fail

merged = merged[keep].copy()
logger.info(f"QC cell filtering: {n_before} -> {merged.n_obs} cells (removed {n_before - merged.n_obs})")

# ── Gene QC: flag low-expression genes (do NOT drop) ────────────────────────────
if min_cells_per_gene is not None:
    from scipy.sparse import issparse
    if issparse(merged.X):
        cells_per_gene = np.asarray((merged.X > 0).sum(axis=0)).flatten()
    else:
        cells_per_gene = np.asarray((merged.X > 0).sum(axis=0)).flatten()
    merged.var['n_cells'] = cells_per_gene
    merged.var['is_filtered'] = cells_per_gene < min_cells_per_gene
    n_filtered_genes = merged.var['is_filtered'].sum()
    logger.info(f"Gene QC: flagged {n_filtered_genes}/{merged.n_vars} genes as is_filtered (expressed in < {min_cells_per_gene} cells)")

# ── Provenance join: donor_id + tissue_id from the samples sheet ─────────────────
# obs['sample_id'] holds the SNP-demux call (or the capture-name fallback). The
# samples sheet is a clean 1:1 sample_id -> donor_id / tissue_id crosswalk, so
# join on sample_id. Cells with no sample match (doublets, unassigned, the "NA"
# sentinel) get "NA" — note every capture here is pooled (multi-donor), so there
# is no unambiguous capture_id fallback for donor/tissue.
if samples_filter is not None and {"donor_id", "tissue_id"}.issubset(samples_filter.columns):
    by_sample = samples_filter.drop_duplicates("sample_id").set_index("sample_id")
    sid = merged.obs["sample_id"].astype(object)
    for col in ("donor_id", "tissue_id"):
        mapped = sid.map(by_sample[col])
        merged.obs[col] = mapped.where(mapped.notna(), "NA").astype(object).values
    n_matched = sid.isin(by_sample.index).sum()
    logger.info(
        f"Joined donor_id/tissue_id from samples sheet: {n_matched}/{merged.n_obs} "
        f"cells matched on sample_id"
    )
else:
    logger.info("No donor_id/tissue_id columns in samples sheet — skipping provenance join")

# ── Grouping keys: capture_id + pool_id ──────────────────────────────────────────
# create_anndata writes obs['capture']; ad.concat adds obs['batch'] from the
# per-capture file stems, which are identical to 'capture' here. Expose them
# under the deposit's grouping-key names. (pool_id == capture_id in this dataset:
# the 10x capture IS the demux pool, there is no coarser pool grouping.)
if "capture" in merged.obs.columns:
    merged.obs["capture_id"] = merged.obs["capture"].astype(object)
    merged.obs = merged.obs.drop(columns=["capture"])
if "batch" in merged.obs.columns:
    merged.obs["pool_id"] = merged.obs["batch"].astype(object)
    merged.obs = merged.obs.drop(columns=["batch"])

# ── CLASSIC string encoding for the var/obs index ────────────────────────────────
# ad.concat produces a pandas StringDtype (nullable) var index, which anndata
# writes as 'nullable-string-array' — older readers silently interpret that as an
# EMPTY gene axis (the bug that dropped all 36601 genes last time). Force a plain
# object-dtype index via DIRECT `.index =` assignment. (Reassigning the whole
# frame — `merged.var = df` — does NOT work: anndata's DataFrame setter re-infers
# a pandas StringDtype under pandas >= 2.x string inference, reintroducing the
# nullable encoding. String columns here are categorical, which every anndata
# version reads, so only the index needs coercing.)
merged.var.index = pd.Index(np.asarray(merged.var.index, dtype=object),
                            name=merged.var.index.name)
merged.obs.index = pd.Index(np.asarray(merged.obs.index, dtype=object),
                            name=merged.obs.index.name)

# Save merged object
logger.info(f"Saving merged object to: {output_file}")
Path(output_file).parent.mkdir(parents=True, exist_ok=True)
merged.write_h5ad(output_file)

# Belt-and-suspenders: if any anndata/pandas combination still wrote the index as
# nullable-string-array (a group with values+mask), rewrite it to the classic
# string-array dataset encoding in place.
import h5py
def _force_classic_index(path, axis):
    with h5py.File(path, "r+") as hf:
        node = hf[axis]["_index"]
        if not isinstance(node, h5py.Group):
            return  # already classic string-array
        vals = node["values"][:]
        del hf[axis]["_index"]
        ds = hf[axis].create_dataset("_index", data=vals,
                                     dtype=h5py.string_dtype(encoding="utf-8"))
        ds.attrs["encoding-type"] = "string-array"
        ds.attrs["encoding-version"] = "0.2.0"
        logger.info(f"Rewrote {axis}/_index nullable-string-array -> classic string-array")

for _axis in ("var", "obs"):
    with h5py.File(output_file, "r") as _hf:
        if dict(_hf[_axis]["_index"].attrs).get("encoding-type") != "string-array":
            _rewrite = True
        else:
            _rewrite = False
    if _rewrite:
        _force_classic_index(output_file, _axis)

# ── Post-write verification: re-read and assert the gene axis survived ───────────
check = ad.read_h5ad(output_file, backed="r")
with h5py.File(output_file, "r") as _hf:
    _var_enc = dict(_hf["var"]["_index"].attrs).get("encoding-type")
logger.info(f"Re-read check: n_obs={check.n_obs} n_vars={check.n_vars} "
            f"var0={check.var_names[0]!r} var_index_enc={_var_enc}")
if expected_n_vars is not None:
    assert _var_enc == "string-array", (
        f"var index encoding is {_var_enc}, not classic string-array."
    )
if expected_n_vars is not None:
    assert check.n_vars == int(expected_n_vars), (
        f"GENE AXIS LOST: re-read n_vars={check.n_vars}, expected {expected_n_vars}. "
        f"The var encoding dropped genes — fix before using this object."
    )
    v0 = str(check.var_names[0])
    assert v0.startswith("ENSG"), (
        f"First var_name {v0!r} is not an Ensembl id — var index encoding is wrong."
    )
    logger.info(f"Gene-axis assertion passed: n_vars={check.n_vars}, var0={v0}")
check.file.close()

logger.info("Done merging captures")
