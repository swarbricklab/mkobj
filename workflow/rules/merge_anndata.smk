# workflow/rules/merge_anndata.smk
# Rules for merging per-capture AnnData objects

rule merge_anndata_captures:
    """
    Merge all per-capture AnnData objects into a single merged object.
    Uses anndata.concat for proper concatenation.
    """
    input:
        h5ad_files = expand(out_dir / "per_capture/{capture}.h5ad", capture=captures)
    output:
        merged = out_dir / "merged.h5ad"
    log:
        log_dir / "merge_anndata.log"
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/merge_anndata.py"
