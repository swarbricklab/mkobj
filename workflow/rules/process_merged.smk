# workflow/rules/process_merged.smk
# Rules for normalizing, computing HVGs, PCA, UMAP, and clustering
# on the merged AnnData object.
#
# Embeddings are computed on an HVG subset then transferred back to
# the full-gene object — no genes are dropped.

_process_params = dict(
    n_top_genes = config.get('params', {}).get('n_top_genes', 2000),
    hvg_batch_key = config.get('params', {}).get('hvg_batch_key', None),
    n_pcs = config.get('params', {}).get('n_pcs', 50),
    leiden_resolutions = config.get('params', {}).get('leiden_resolutions', [0.1, 0.3, 0.5, 1.0]),
)

if use_subsets:
    rule process_merged_anndata:
        """
        Normalize, compute HVGs, PCA, UMAP, and Leiden clustering
        on the merged AnnData object for a named subset.
        """
        input:
            merged = out_dir / "{subset}/merged.h5ad"
        output:
            processed = out_dir / "{subset}/processed.h5ad"
        log:
            log_dir / "{subset}/process_merged.log"
        params:
            **_process_params
        conda:
            "../envs/scanpy.yaml"
        script:
            "../scripts/process_merged_data.py"
else:
    rule process_merged_anndata:
        """
        Normalize, compute HVGs, PCA, UMAP, and Leiden clustering
        on the merged AnnData object.
        """
        input:
            merged = out_dir / "merged.h5ad"
        output:
            processed = out_dir / "processed.h5ad"
        log:
            log_dir / "process_merged.log"
        params:
            **_process_params
        conda:
            "../envs/scanpy.yaml"
        script:
            "../scripts/process_merged_data.py"
