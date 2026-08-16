# workflow/rules/format_cellxgene.smk
# Rules for formatting the processed AnnData object for CELLxGENE upload.
#
# Joins external metadata tables, maps metadata to CELLxGENE schema fields,
# renames var columns, and fills defaults for required fields.

_cellxgene_params = config.get('cellxgene', {})

# Join tables declared under cellxgene.metadata are listed as rule inputs so
# Snakemake re-runs the format step when one of them changes.
_cellxgene_metadata = [
    entry['path']
    for entry in (_cellxgene_params.get('metadata') or [])
    if entry.get('path')
]

if use_subsets:
    rule format_cellxgene:
        """
        Format the processed AnnData object for CELLxGENE upload.
        """
        input:
            processed = out_dir / "{subset}/processed.h5ad",
            metadata = _cellxgene_metadata
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
            processed = out_dir / "processed.h5ad",
            metadata = _cellxgene_metadata
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
