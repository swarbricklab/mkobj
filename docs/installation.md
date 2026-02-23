# Installation

## Prerequisites

- Conda
- Git (for submodule installation)
- DVC (for data tracking, optional)
- [qxub](https://github.com/swarbricklab/qxub) (for PBS cluster submission)

## Snakemake Environment

This workflow requires Snakemake 8+ with the cluster-generic executor plugin.

Create or activate the Snakemake environment:

```bash
conda activate snakemake_8.30.0
```

The environment should include:
- `snakemake >= 8.0`
- `snakemake-executor-plugin-cluster-generic`
- `qxub` / `qxtat` (PBS submission and status checking)

## Installing as a Submodule

To use this workflow in a dataset project:

```bash
cd /path/to/dataset
git submodule add <repo-url> modules/mkobj
git submodule update --init --recursive
```

## Rule Environments

Rule-level conda environments are defined in `workflow/envs/` and installed automatically
by Snakemake on the first run (via Apptainer/Conda software deployment).

### Seurat (R)

Defined in `workflow/envs/seurat.yaml`:
- R Seurat 5.1
- SeuratObject
- qs (fast serialization)
- tidyverse (readr, dplyr, tibble)

### Scanpy (Python)

Defined in `workflow/envs/scanpy.yaml`:
- Python >= 3.10
- scanpy >= 1.10
- anndata >= 0.10
- pandas, numpy, scipy, h5py
