#!/bin/bash
# Run mkobj workflow in module mode (as a submodule of a dataset)
# Snakemake orchestrates on a compute node; individual rules are dispatched
# to further compute nodes via the cluster profile (which also uses qxub).

set -e

module=mkobj
profile="--profile modules/$module/profiles/cluster"

mkdir -p logs/joblogs

# Pre-create conda envs on login node (compute nodes lack internet)
qx --env snakemake_8.30.0 --mem 16GB  --queue copyq -- \
    snakemake \
        --snakefile workflow/Snakefile \
        --configfile config/test.yaml \
        --software-deployment-method conda \
        --conda-create-envs-only

# Run workflow on a compute node
conda activate qxub
qx --env snakemake_8.30.0 --mem 4GB --cpus 1 --runtime 12h -- \
    snakemake $profile \
    --snakefile modules/$module/workflow/Snakefile \
    --configfile config/$module/config.yaml \
    "$@"
