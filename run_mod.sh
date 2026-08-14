#!/bin/bash
# Run mkobj workflow in module mode (as a submodule of a dataset)
# Snakemake orchestrates on a compute node; individual rules are dispatched
# to further compute nodes via the cluster profile (which also uses qxub).
#
# Run from the DATASET root, not from inside the submodule:
#   ./modules/mkobj/run_mod.sh [snakemake args...]

set -e
eval "$(conda shell.bash hook)"
conda activate snakemake_8.30.0

module=mkobj
profile="--profile modules/$module/profiles/cluster"
snakefile="modules/$module/workflow/Snakefile"
configfile="config/$module/config.yaml"

mkdir -p logs/joblogs

# Pre-create conda envs on an internet-capable node (compute nodes lack
# outbound access). Passes the profile so envs land in its conda-prefix (the
# shared a56 cache) rather than ./.snakemake/conda; the profile also supplies
# software-deployment-method, so that flag is not repeated here.
qx --env snakemake_8.30.0 --mem 16GB --queue copyq -- \
    snakemake $profile \
        --snakefile "$snakefile" \
        --configfile "$configfile" \
        --conda-create-envs-only

# Run workflow on a compute node
qx --env snakemake_8.30.0 --mem 4GB --cpus 1 --runtime 12h -- \
    snakemake $profile \
        --snakefile "$snakefile" \
        --configfile "$configfile" \
        "$@"
