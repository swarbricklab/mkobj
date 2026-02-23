# Running the Workflow

## Prerequisites

1. [Install](installation.md) the workflow environment
2. [Configure](../config/README.md) the workflow for your dataset
3. Prepare input data (Cell Ranger outputs)

## Module Mode

When embedded in a dataset project:

```bash
cd /path/to/dataset
./modules/mkobj/run_mod.sh
```

### Common options

```bash
# Dry run (show what would be executed)
./modules/mkobj/run_mod.sh -n

# Force re-run all rules
./modules/mkobj/run_mod.sh --forceall

# Run specific rule
./modules/mkobj/run_mod.sh merge_captures

# Run only the AnnData pipeline
./modules/mkobj/run_mod.sh merge_anndata_captures
```

## Standalone Mode

For testing:

```bash
cd /path/to/mkobj
./run_test.sh
```

## Outputs

After successful completion, outputs will be in the configured `outs.results` directory:

| File | Description |
|------|-------------|
| `merged.qs` | Merged Seurat object (R) |
| `merged.h5ad` | Merged AnnData object (Python) |

Both files contain the same cells and metadata in their respective formats.

Per-capture intermediate files (`per_capture/{capture}.rds`, `per_capture/{capture}.h5ad`)
are marked as `temp()` and automatically removed after merging.

## Cluster Execution

In module mode, `run_mod.sh` uses the cluster profile at `profiles/cluster/config.yaml`,
which submits jobs to PBS via [qxub](https://github.com/swarbricklab/qxub).

Key cluster settings:
- Jobs are submitted with per-rule resource allocation (memory, CPUs, walltime)
- Job status is tracked via `qxtat check --snakemake`
- Platform-specific settings (project, storage volumes) are handled by qxub

## Logs

Logs for each rule are saved to the configured `outs.logs` directory:
- `create_seurat/{capture}.log` — Per-capture Seurat object creation
- `create_anndata/{capture}.log` — Per-capture AnnData object creation
- `merge_captures.log` — Seurat merge step
- `merge_anndata.log` — AnnData merge step

## Troubleshooting

### Common Issues

1. **Missing captures.csv**: Ensure `deps.captures` points to a valid CSV with a `capture` column
2. **Cell Ranger paths not found**: Check that `deps.cellranger/{capture}` directories exist with `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`
3. **Memory errors on merge**: Increase `mem_mb` for `merge_captures` / `merge_anndata_captures` in `profiles/cluster/config.yaml` under `set-resources`
4. **h5ad write errors with mixed types**: The AnnData pipeline converts NA values in string columns to `"NA"` strings for h5py compatibility. If you encounter type errors, check your metadata files for unexpected mixed types.
