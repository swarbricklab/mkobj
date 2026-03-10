# workflow/rules/format_cellxgene.smk
# Rules for formatting the processed AnnData object for CELLxGENE upload.
#
# Maps metadata to CELLxGENE schema fields, renames var columns,
# and fills defaults for required fields.

_cellxgene_params = config.get('cellxgene', {})

if use_subsets:
    rule format_cellxgene:
        """
        Format the processed AnnData object for CELLxGENE upload.
        """
        input:
            processed = out_dir / "{subset}/processed.h5ad"
        output:
            cellxgene = out_dir / "{subset}/cellxgene.h5ad"
        log:
            log_dir / "{subset}/format_cellxgene.log"
        params:
            cellxgene = _cellxgene_params
        conda:
            "../envs/scanpy.yaml"
        script:
            "../scripts/format_cellxgene.py"
else:
    rule format_cellxgene:
        """
        Format the processed AnnData object for CELLxGENE upload.
        """
        input:
            processed = out_dir / "processed.h5ad"
        output:
            cellxgene = out_dir / "cellxgene.h5ad"
        log:
            log_dir / "format_cellxgene.log"
        params:
            cellxgene = _cellxgene_params
        conda:
            "../envs/scanpy.yaml"
        script:
            "../scripts/format_cellxgene.py"
