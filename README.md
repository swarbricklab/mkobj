# mkobj

Create merged Seurat and AnnData objects from Cell Ranger filtered feature barcode matrices.

## Overview

This workflow processes 10X Chromium single-cell data from Cell Ranger output and creates both:
- **Seurat objects** (`.qs`) via the R/Bioconductor ecosystem
- **AnnData objects** (`.h5ad`) via the Python/Scanpy ecosystem

It supports:
- **Unimodal data**: Gene expression only
- **Multimodal data**: Gene expression + Antibody Capture (CITE-seq)

The workflow can optionally attach:
- Sample assignments from SNP-based demultiplexing
- Cell type annotations
- Ambient RNA profiles
- Sample-level subsetting (keep only cells from a specified cohort)

## Modes

This workflow can be run in either [standalone mode](#standalone-mode) or [module mode](#module-mode).

### Standalone mode

In "standalone" mode, the data is included in the same repo as the workflow.
This mode is used mainly for testing.

```bash
./run_test.sh
```

### Module mode

This workflow can be embedded into a dataset as a [git submodule](https://www.atlassian.com/git/tutorials/git-submodule).

To use in module mode:
1. Add this workflow as a submodule to your dataset
2. Copy and configure the config file
3. Run the workflow using `run_mod.sh`

```bash
# From the dataset root
git submodule add <repo-url> modules/mkobj
mkdir -p config/mkobj
cp modules/mkobj/config/template.yaml config/mkobj/config.yaml
# Edit config.yaml for your dataset
./modules/mkobj/run_mod.sh
```

## Workflow Structure

The workflow is organized into parallel Seurat and AnnData pipelines that run concurrently:

```
Cell Ranger matrices
        │
        ├──► create_seurat_object (per capture) ──► merge_captures ──► merged.qs
        │
        └──► create_anndata_object (per capture) ──► merge_anndata_captures ──► merged.h5ad
```

### Seurat Pipeline (R)
1. **create_seurat_object**: Create individual Seurat objects per capture
   - Reads Cell Ranger filtered matrices (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`)
   - Handles multimodal data (GEX + AB stored as a separate assay)
   - Attaches sample assignments, annotations, and ambient profiles to cell metadata
   - Optionally subsets cells to a specified cohort via `samples.csv`
   - Prefixes barcodes with capture ID for uniqueness across captures

2. **merge_captures**: Merge all per-capture objects into one
   - Combines all captures using Seurat's `merge()` function
   - Joins layers for proper integration
   - Output: `merged.qs`

### AnnData Pipeline (Python)
1. **create_anndata_object**: Create individual AnnData objects per capture
   - Reads Cell Ranger filtered matrices (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`)
   - Handles multimodal data — gene expression in `X`, antibody capture in `obsm['AB']` with feature names in `uns['AB_features']`
   - Attaches sample assignments, annotations, and ambient profiles to `obs`
   - Optionally subsets cells to a specified cohort via `samples.csv` (keeps cohort singlets, doublets, and unassigned cells)
   - Prefixes barcodes with capture ID for uniqueness across captures
   - Converts string columns with NA values to proper string type for h5ad compatibility

2. **merge_anndata_captures**: Merge all per-capture objects into one
   - Combines all captures using `anndata.concat()` with `join='outer'`
   - Preserves multimodal data in `obsm` across captures
   - Output: `merged.h5ad`

## Re-running with additional captures

Adding a capture is **not** an incremental operation. The per-capture objects are
`temp()` — `create_seurat_object` (`per_capture/{capture}.rds`),
`create_anndata_object` (`per_capture_raw/{capture}.h5ad`) and `detect_doublets`
(`per_capture/{capture}.h5ad`) — so they are deleted once a merge succeeds. On the
next run they are missing, and every merge that consumes them has to rebuild them.

Per-capture objects are **shared across subsets**, so the rebuild is not confined to
the subset you edited. Adding one capture to subset A also regenerates the captures
A shares with B, which in turn makes B's merged object out of date:

```
subset A: X, Y     add W to A only     rebuilds X, Y, W  and  Y, Z
subset B: Y, Z          ────────►      re-merges BOTH A and B
```

Snakemake reports the two cases differently — `Set of input files has changed since
last execution` for the subset you actually edited, and `Input files updated by
another job` for the ones dragged in with it. Expect work proportional to the union
of every subset that shares a capture with the one you changed, which in a
well-connected dataset is most of it.

Always check the size of the job before committing to it:

```bash
./modules/mkobj/run_mod.sh -n     # per-rule job counts
```

### Keeping the per-capture objects

If a dataset grows by a few captures at a time, the rebuild is usually worse than the
disk it saves. Snakemake's `--notemp` ignores the `temp()` declarations and leaves the
per-capture objects in place:

```bash
./modules/mkobj/run_mod.sh --notemp
```

Every later addition is then genuinely incremental: only the new capture is built, and
only the subsets whose capture list actually changed are re-merged — the cascade above
does not happen, because the shared captures are still on disk and unchanged.

The cost is one full set of per-capture objects. Measured on a 148-capture dataset:
roughly 12 GB of `per_capture/` plus 35 GB of `per_capture_raw/`. Weigh that against
recomputing three rules per capture across every connected subset.

`temp()` is the better default for a one-shot build, where those files are pure waste.
It is the wrong default for a dataset you keep adding to.

### Under DVC

Datasets embedding mkobj as a submodule usually wrap it in a DVC stage. DVC removes a
stage's outputs *before* re-running it, so with no further configuration the merged
objects for every subset are deleted first — nothing can be skipped even in principle,
and a failed run leaves neither the new results nor the old ones.

`persist: true` on the stage outputs prevents that deletion:

```yaml
outs:
  - data/objects:
      persist: true
```

**But check your `cache.type` first.** `persist` does not simply skip the removal — in
`Stage.remove_outs()` it substitutes `unprotect` for `remove`, and unprotecting a
symlinked or hardlinked file means byte-copying it back into the workspace. With
`cache.type = hardlink,symlink` the pre-run step becomes a full copy of every file in
the stage's outputs, and it is charged before any compute starts. On a stage with
multi-terabyte outputs this can cost far more than the rebuild it was meant to avoid,
and it duplicates data the links existed to share.

With the default `cache.type=copy` there is no link to break, so `persist` is cheap and
the above does not apply.

Where the unprotect pass is too expensive, run the workflow directly and register the
result afterwards:

```bash
./modules/mkobj/run_mod.sh --notemp     # snakemake decides what to skip
dvc commit <stage>                      # record the outputs DVC did not produce
```

This keeps the incremental behaviour at the cost of the stage no longer being
reproducible purely through `dvc repro` — a deliberate trade, worth a note in the
dataset repo when you make it.

## Configuration

See the [configuration guide](config/README.md) for detailed instructions.

Quick start:
```yaml
deps:
  cellranger: "data/counts"
  captures: "config/mkobj/captures.csv"
  demux: "data/demux"

outs:
  results: "data/objects"
  logs: "logs/mkobj"
```

## Outputs

| File | Description |
|------|-------------|
| `merged.qs` | Merged Seurat object (R, serialized with [qs](https://github.com/traversc/qs)) |
| `merged.h5ad` | Merged AnnData object (Python, [HDF5-backed](https://anndata.readthedocs.io/en/latest/fileformat-prose.html)) |

Both outputs contain identical data — the same cells, metadata, and (where applicable) multimodal assays — in their respective ecosystem formats.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) >= 8.0
- [snakemake-executor-plugin-cluster-generic](https://github.com/snakemake/snakemake-executor-plugin-cluster-generic)
- [qxub](https://github.com/swarbricklab/qxub) (for PBS cluster submission)
- Conda
- Apptainer (for containerized environments)

Rule-level conda environments are defined in `workflow/envs/` and installed automatically:

| Environment | Key packages |
|-------------|-------------|
| `seurat.yaml` | R, Seurat 5.1, SeuratObject, qs, tidyverse |
| `scanpy.yaml` | Python >= 3.10, scanpy >= 1.10, anndata >= 0.10, pandas, numpy, scipy |

## Authors

Originally developed as part of the Swarbrick Lab data processing pipelines.
