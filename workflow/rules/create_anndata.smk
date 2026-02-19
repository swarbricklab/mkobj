# workflow/rules/create_anndata.smk
# Rules for creating individual AnnData objects per capture

rule create_anndata_object:
    """
    Create an AnnData object from Cell Ranger filtered feature barcode matrix.
    Handles both unimodal (gene expression only) and multimodal (GEX + antibody capture) data.
    Attaches sample assignments, cell annotations, and ambient profiles if available.
    Optionally subsets to samples listed in samples.csv for the current capture.
    Prefixes barcodes with capture ID to ensure uniqueness when merging.
    """
    input:
        matrix_dir = lambda wc: Path(config['deps']['cellranger']) / wc.capture
    output:
        h5ad = temp(out_dir / "per_capture/{capture}.h5ad")
    log:
        log_dir / "create_anndata/{capture}.log"
    params:
        assignment_root = config['deps'].get('demux', ''),
        annotation_root = config['deps'].get('annotation', ''),
        ambient_root = config['deps'].get('ambient', ''),
        samples_file = config['deps'].get('samples', ''),
        modality = config.get('params', {}).get('modality', 'auto')
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/create_anndata.py"
