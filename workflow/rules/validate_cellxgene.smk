# workflow/rules/validate_cellxgene.smk
# Validate the CELLxGENE object with the schema validator and save its report.
#
# Ported from the brca_major_atlas CELLxGENE curation step 10_validate.sh. The
# validator reads both matrices in full (a compute job). It exits with status 1
# on an invalid object, and pipefail passes that status on, so the rule fails and
# the invalid object is not consumed downstream. The report holds a line with
# is_valid=True when the object is valid, and one `WARNING:` line per warning.
#
# Runs in its own env (envs/cxg_schema.yaml) because cellxgene-schema pins
# numpy/anndata versions that conflict with the scanpy env. Defined only when
# `reconcile.enabled` is true (the gene lift is what makes the object schema-valid).

if reconcile_enabled():

    if use_subsets:
        rule validate_cellxgene:
            """Validate a named subset's CELLxGENE object against the schema."""
            input:
                cellxgene = out_dir / "{subset}/cellxgene.h5ad"
            output:
                report = out_dir / "{subset}/validation_report.txt"
            log:
                log_dir / "{subset}/validate_cellxgene.log"
            conda:
                "../envs/cxg_schema.yaml"
            shell:
                "set -o pipefail; "
                "cellxgene-schema validate {input.cellxgene} 2>&1 | tee {output.report}"
    else:
        rule validate_cellxgene:
            """Validate the CELLxGENE object against the schema."""
            input:
                cellxgene = out_dir / "cellxgene.h5ad"
            output:
                report = out_dir / "validation_report.txt"
            log:
                log_dir / "validate_cellxgene.log"
            conda:
                "../envs/cxg_schema.yaml"
            shell:
                "set -o pipefail; "
                "cellxgene-schema validate {input.cellxgene} 2>&1 | tee {output.report}"
