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
