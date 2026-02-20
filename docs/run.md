# Running the Workflow

## Prerequisites

1. [Install](installation.md) the workflow environment
2. [Configure](../config/README.md) the workflow for your dataset
3. Prepare input data (Cell Ranger outputs)

## Module Mode

When embedded in a dataset project:

```bash
cd /path/to/dataset
./modules/mkseurat/run_mod.sh
```

### Common options

```bash
# Dry run (show what would be executed)
./modules/mkseurat/run_mod.sh -n

# Force re-run all rules
./modules/mkseurat/run_mod.sh --forceall

# Run specific rule
./modules/mkseurat/run_mod.sh merge_captures

# Limit parallelism
./modules/mkseurat/run_mod.sh -j 4
```

## Standalone Mode

For testing:

```bash
cd /path/to/mkseurat
./run_test.sh
```

## Outputs

After successful completion, outputs will be in the configured `outs.results` directory:

- `merged.qs` - Merged Seurat object
- `merged.h5ad` - Merged AnnData object

## Logs

Logs for each rule are saved to the configured `outs.logs` directory.

## Troubleshooting

### Common Issues

1. **Missing captures.csv**: Ensure `deps.captures` points to a valid CSV with a `capture` column
2. **Cell Ranger paths not found**: Check that `deps.cellranger/{capture}` directories exist
3. **Memory errors on merge**: Increase memory allocation in `profiles/workflow/config.yaml`
