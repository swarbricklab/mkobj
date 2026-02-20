# mkseurat

Create merged Seurat and AnnData objects from Cell Ranger filtered feature barcode matrices.

## Overview

This workflow processes 10X Chromium single-cell data from Cell Ranger output and creates both:
- **Seurat objects** (R/Bioconductor ecosystem)
- **AnnData objects** (Python/ScanPy ecosystem)

It supports:
- **Unimodal data**: Gene expression only
- **Multimodal data**: Gene expression + Antibody Capture (CITE-seq)

The workflow can optionally attach:
- Sample assignments from SNP-based demultiplexing
- Cell type annotations
- Ambient RNA profiles

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
git submodule add <repo-url> modules/mkseurat
mkdir -p config/mkseurat
cp modules/mkseurat/config/template.yaml config/mkseurat/config.yaml
# Edit config.yaml for your dataset
./modules/mkseurat/run_mod.sh
```

## Workflow Structure

The workflow is organized into parallel Seurat and AnnData pipelines:

### Seurat Pipeline (R)
1. **create_seurat_object**: Create individual Seurat objects per capture
   - Reads Cell Ranger filtered matrices
   - Handles multimodal data (GEX + AB)
   - Attaches sample assignments, annotations, and ambient profiles
   - Prefixes barcodes with capture ID

2. **merge_captures**: Merge all capture objects into one
   - Combines all captures using Seurat's merge function
   - Joins layers for proper integration

### AnnData Pipeline (Python)
1. **create_anndata_object**: Create individual AnnData objects per capture
   - Reads Cell Ranger filtered matrices
   - Handles multimodal data (GEX + AB stored in obsm)
   - Attaches sample assignments, annotations, and ambient profiles
   - Prefixes barcodes with capture ID

2. **merge_anndata_captures**: Merge all capture objects into one
   - Combines all captures using anndata.concat

## Configuration

See the [configuration guide](config/README.md) for detailed instructions.

Quick start:
```yaml
deps:
  cellranger: "data/sc/cellranger/count/filtered"
  captures: "config/chromium-preprocessing/captures.csv"
  demux: "data/snp/demux/assignment"

outs:
  results: "data/sc/seurat"
  logs: "logs/mkseurat"
```

## Outputs

| File | Description |
|------|-------------|
| `merged.qs` | Merged Seurat object in qs format |
| `merged.h5ad` | Merged AnnData object in h5ad format |

## Requirements

- Snakemake 7.32.4
- Conda/Mamba
- R packages: Seurat 5.1, qs, tidyverse
- Python packages: scanpy>=1.10, anndata>=0.10, pandas, numpy

## Authors

Originally developed as part of the Swarbrick Lab data processing pipelines.
