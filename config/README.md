# Configuration

Configuring the workflow involves editing config files and preparing input data.

## Config file

To use this workflow as a module in a dataset project, copy the template config file:

```bash
cd top/of/project
mkdir -p config/mkobj
cp modules/mkobj/config/template.yaml config/mkobj/config.yaml
```

Edit the config file to match your dataset structure.

### Required dependencies

| Key | Description |
|-----|-------------|
| `deps.cellranger` | Path to Cell Ranger filtered feature barcode matrices. Each capture should be a subdirectory containing `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` |
| `deps.captures` | CSV file listing captures to process. Must have a `capture` column |

### Optional dependencies

| Key | Description |
|-----|-------------|
| `deps.samples` | CSV file for subsetting cells by sample. Must have `sample_id` and `capture_id` columns. For each capture, only cells assigned to listed samples are kept. Doublets and unassigned cells are always retained. |
| `deps.demux` | Path to sample assignment files from SNP demux. Each capture directory should contain `cell_assignment.tsv` |
| `deps.annotation` | Path to cell type annotations. Each capture directory should contain `cell_types.csv` |
| `deps.ambient` | Path to ambient RNA profiles. Each capture directory should contain `ambient_summary.csv` |

### Outputs

| Key | Description |
|-----|-------------|
| `outs.results` | Directory for output files (`merged.qs` and `merged.h5ad`) |
| `outs.logs` | Directory for log files |

### Parameters

| Key | Description | Default |
|-----|-------------|---------|
| `params.modality` | How to handle multimodal data (`auto`, `Gene Expression`) | `auto` |

When `modality` is `auto`, multimodal captures will use Gene Expression as the primary assay
and store Antibody Capture data as an additional modality (Seurat: separate assay; AnnData: `obsm['AB']`).

## Captures CSV

The captures CSV file must have at minimum a `capture` column listing the capture IDs to process.
Capture IDs should match the subdirectory names in `deps.cellranger`.

Example:
```csv
capture
NC001
NC002
NC003
```

## Samples file (for subsetting)

The optional samples CSV file (`deps.samples`) is used to subset cells to a specific cohort.
It must have `sample_id` and `capture_id` columns.

For each capture, cells are kept if:
- Their assigned `sample_id` is in the samples file for that capture, OR
- They are marked as `status == "doublet"`, OR
- Their `sample_id` is NA (unassigned cells)

This is useful when multiple studies are pooled in the same capture but you only want cells from your cohort.

Example:
```csv
sample_id,capture_id,tissue_id,donor_id
4063,Atlas_Pool_2a,4063,4063
4063,Atlas_Pool_2b,4063,4063
4218,Atlas_Pool_2a,4218,4218
```

This filtering logic is applied identically in both the Seurat and AnnData pipelines.

## Sample assignment files

If using SNP-based demultiplexing, each capture should have a `cell_assignment.tsv` file with columns:
- `barcode`: Cell barcode
- `status`: Assignment status (singlet, doublet, unassigned)
- `assignment`: Sample assignment (stored as `sample_id` in cell metadata)

If no demux data is configured, all cells default to `status = "singlet"` and `sample_id = {capture}`.
