# workflow/rules/metadata_anndata.smk
# Rules for attaching experimental metadata to AnnData

rule attach_metadata_anndata:
    """
    Attach experimental metadata (e.g., sample information, clinical data)
    to the merged AnnData object. Joins metadata on sample/assignment column.
    """
    input:
        merged = out_dir / "merged.h5ad",
        metadata = config['deps']['metadata']
    output:
        annotated = out_dir / "merged_annotated.h5ad"
    log:
        log_dir / "attach_metadata_anndata.log"
    params:
        join_column = config.get('params', {}).get('metadata_join_column', 'sample_id')
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/attach_metadata_anndata.py"
