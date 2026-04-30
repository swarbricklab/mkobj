#!/bin/bash
# Run mkobj workflow in standalone mode (for testing)
# Executes on a compute node via qxub

set -e

# Pre-create conda envs on login node (compute nodes lack internet)
qx --env snakemake_8.30.0 --mem 16GB  --queue copyq -- \
    snakemake \
        --snakefile workflow/Snakefile \
        --configfile config/test.yaml \
        --software-deployment-method conda \
        --conda-create-envs-only

# Run workflow on a compute node

qx --env snakemake_8.30.0 --mem 16GB --cpus 4 --runtime 1h -- \
    snakemake \
        --snakefile workflow/Snakefile \
        --configfile config/test.yaml \
        --software-deployment-method conda \
        --cores 4 \
        "$@"
