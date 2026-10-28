# workflow/rules/reconcile_gene_space.smk
# Optional gene-space reconciliation (gene lift).
#
# Lifts the object's gene space from its source GENCODE build (the Cell Ranger
# reference) to a target build (e.g. GENCODE v48, the CELLxGENE 7.1.0 gene set):
# retired identifiers are renamed to their successor or dropped, following the
# brca_major_atlas CELLxGENE curation crosswalk. Counts are carried unchanged.
#
# Early mode (default): reads merged.h5ad, feeds process_merged_anndata the
#   reconciled object -> HVG/PCA/clustering run on the final gene set.
# Late mode: reads processed.h5ad, feeds format_cellxgene the reconciled object
#   -> the published embedding is preserved; var-axis slicing also slices
#   varm['PCs'], layers['counts'] and X.
#
# The rule is defined ONLY when `reconcile.enabled` is true, so a config without
# the `reconcile:` block leaves the DAG byte-identical to upstream.

if reconcile_enabled():

    _reconcile_refs = reference_paths()
    _reconcile_params = config.get('reconcile', {})

    if use_subsets:
        rule reconcile_gene_space:
            """Lift the gene space of a named subset's object to the target build."""
            input:
                object = reconcile_source_object,
                source_gtf = _reconcile_refs['source_gtf'],
                target_gtf = _reconcile_refs['target_gtf'],
                stable_id_history = _reconcile_refs['stable_id_history'],
            output:
                reconciled = out_dir / "{subset}/reconciled.h5ad",
                crosswalk = out_dir / "{subset}/gene_crosswalk.csv",
            log:
                log_dir / "{subset}/reconcile_gene_space.log"
            params:
                reconcile = _reconcile_params
            conda:
                "../envs/scanpy.yaml"
            script:
                "../scripts/reconcile_gene_space.py"
    else:
        rule reconcile_gene_space:
            """Lift the object's gene space to the target build."""
            input:
                object = reconcile_source_object,
                source_gtf = _reconcile_refs['source_gtf'],
                target_gtf = _reconcile_refs['target_gtf'],
                stable_id_history = _reconcile_refs['stable_id_history'],
            output:
                reconciled = out_dir / "reconciled.h5ad",
                crosswalk = out_dir / "gene_crosswalk.csv",
            log:
                log_dir / "reconcile_gene_space.log"
            params:
                reconcile = _reconcile_params
            conda:
                "../envs/scanpy.yaml"
            script:
                "../scripts/reconcile_gene_space.py"
