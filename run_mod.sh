#!/bin/bash
# Run mkseurat workflow in module mode (as a submodule of a dataset)

set -e
eval "$(conda shell.bash hook)"
conda activate snakemake_7.32.4

module=mkobj

if [[ "$(hostname)" == *"nci"* ]]; then
    echo "Running on NCI"
    global_profile="--profile modules/$module/profiles/global/nci"
    workflow_profile="--workflow-profile modules/$module/profiles/workflow"
    module load singularity
    mkdir -p logs/joblogs
else
    echo "WARNING: Unknown host: $(hostname)"
    echo "WARNING: No known global profile for this host"
    echo "WARNING: Running without global profile"
    echo "See https://github.com/swarbricklab/snakemake_config/blob/main/README.md"
    global_profile=""
    workflow_profile="--workflow-profile modules/$module/profiles/workflow"
fi

snakemake $global_profile $workflow_profile \
    --snakefile modules/$module/workflow/Snakefile \
    --configfile config/$module/config.yaml \
    "$@"
