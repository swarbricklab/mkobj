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
| `deps.lineage_markers` | Path to a lineage markers CSV used for cross-lineage doublet detection |

### Outputs

| Key | Description |
|-----|-------------|
| `outs.results` | Directory for output files (`merged.qs` and `merged.h5ad`) |
| `outs.logs` | Directory for log files |

### Parameters

| Key | Description | Default |
|-----|-------------|---------|
| `params.modality` | How to handle multimodal data (`auto`, `Gene Expression`) | `auto` |
| `params.scrublet_threshold` | Threshold for flagging standard computational scrublet doublets | `0.2` |
| `params.cross_lineage_threshold` | Ratio of markers expressed to flag crossing over an explicit lineage | `0.7` |
| `params.min_counts` | Minimum total UMI counts per cell | `500` |
| `params.max_counts` | Maximum total UMI counts per cell | `50000` |
| `params.min_features` | Minimum number of genes detected per cell | `200` |
| `params.max_features` | Maximum number of genes detected per cell | `null` (disabled) |
| `params.max_pct_mt` | Maximum percent mitochondrial content | `60` |
| `params.max_pct_ribo` | Maximum percent ribosomal content | `null` (disabled) |
| `params.filter_doublets` | Remove cells flagged as `predicted_doublet` | `true` |
| `params.min_cells_per_gene` | Minimum cells expressing a gene; genes below this are flagged `is_filtered` in `.var` rather than dropped | `3` |
| `params.n_top_genes` | Number of highly variable genes for PCA/UMAP | `2000` |
| `params.hvg_batch_key` | Batch key for HVG selection (e.g. `pool_id`, `capture`) | `null` |
| `params.n_pcs` | Number of principal components | `50` |
| `params.leiden_resolutions` | Leiden clustering resolutions to compute | `[0.1, 0.3, 0.5, 1.0]` |

### CELLxGENE formatting (optional)

The `cellxgene` config block controls how the final `cellxgene.h5ad` is formatted for upload.
It targets [schema 7.1.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.1.0/schema.md).

| Key | Description |
|-----|-------------|
| `cellxgene.metadata` | List of CSV join tables, each `{path, key}`. `key` names an existing `obs` column; every other column in the CSV is written to `obs` under its own name |
| `cellxgene.column_map` | Dict mapping CELLxGENE schema field names to existing `obs` column names |
| `cellxgene.defaults` | Dict of constants, used both when a column is absent and to fill rows a join did not match |
| `cellxgene.uns` | Dict of `uns` values — `organism_ontology_term_id` and `title` |

The three sources are applied in order (`metadata`, then `column_map`, then `defaults`),
so a later one wins. Anything still unset becomes `"unknown"`, which the schema permits
for every ontology field that can be unavailable.

#### Metadata joins

`cellxgene.metadata` is the general way to bring outside metadata in. Because each table
declares its own key, one mechanism covers joins at every grain:

```yaml
cellxgene:
  metadata:
    - path: "config/mkobj/cellxgene_sample.csv"    # per sample  — clinical/ontology terms
      key: "sample_id"
    - path: "config/mkobj/cellxgene_capture.csv"   # per capture — assay chemistry
      key: "capture"
    - path: "config/mkobj/cell_type_cl.csv"        # per cell    — label -> CL term
      key: "fine_cell_type"
```

Join tables are rule inputs, so editing one re-runs the format step.

#### Required fields

Required `obs` fields: `assay_ontology_term_id`, `cell_type_ontology_term_id`,
`development_stage_ontology_term_id`, `disease_ontology_term_id`, `donor_id`,
`is_primary_data`, `self_reported_ethnicity_ontology_term_id`, `sex_ontology_term_id`,
`suspension_type`, `tissue_ontology_term_id`, `tissue_type`.

Required `uns` fields: `organism_ontology_term_id`, `title`.

The `cellxgene` block is global, but a `subsets:` run emits several objects and `title`
must distinguish the datasets in a collection. Any `uns` string may contain the literal
`{subset}`, which is replaced with the subset name:

```yaml
cellxgene:
  uns:
    title: "BRCA consortium chromium — {subset}"
```

`organism_ontology_term_id` is a `uns` field, not `obs` — it moved in schema 6.0.0. A value
left in `obs` by an older config is migrated automatically.

`experimental_condition_ontology_term_id` and the `genetic_perturbation_*` fields are
deliberately not written: the schema says they MUST NOT be present unless the dataset
actually has perturbations. `array_col` / `array_row` / `in_tissue` are Visium-only.

#### var and raw

`feature_biotype`, `feature_length`, `feature_name`, `feature_reference` and `feature_type`
are annotated by CELLxGENE Discover; a curator MUST NOT supply them, so they are dropped
if present. Cell Ranger's own `feature_type` (`Gene Expression` / `Antibody Capture`) is
preserved as `cellranger_feature_type` to avoid the reserved name.

`var['feature_is_filtered']` means "zeroed in `X` but retained in `raw`", so Explorer can
mask the gene. It is computed from the matrices, not taken from mkobj's `is_filtered` —
that flags low-expression genes without zeroing them, and the two are not interchangeable.
`is_filtered` is kept as a non-schema column.

Raw counts are **moved** from `layers['counts']` to `raw.X` as `float32`, which is where
the schema requires them when `X` holds normalized values. Moved rather than copied:
keeping both would duplicate the count matrix in memory and on disk, which on a
million-cell object is tens of GB for no benefit. Read them from `raw.X` downstream.

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
