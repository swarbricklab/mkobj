#!/bin/bash
# Run mkobj workflow in module mode (as a submodule of a dataset)

set -e
eval "$(conda shell.bash hook)"
conda activate snakemake_8  # TODO: update environment name

module=mkobj
profile="--profile modules/$module/profiles/cluster"

mkdir -p logs/joblogs

snakemake $profile \
    --snakefile modules/$module/workflow/Snakefile \
    --configfile config/$module/config.yaml \
    "$@"
