# workflow/rules/create_anndata.smk
# Rules for creating individual AnnData objects per capture
#
# Per-capture objects are created once (with all cells) and shared across
# all subsets.  Cell-level filtering is applied at merge time.

rule create_anndata_object:
    """
    Create an AnnData object from Cell Ranger filtered feature barcode matrix.
    Handles both unimodal (gene expression only) and multimodal (GEX + antibody capture) data.
    Attaches sample assignments, cell annotations, and ambient profiles if available.
    Prefixes barcodes with capture ID to ensure uniqueness when merging.
    """
    input:
        # Declare the auxiliary dep roots as inputs (not just params) so a change to
        # any of them triggers a rerun and a missing root is caught at DAG build.
        # Only those the config actually sets are declared (annotation is dropped for
        # the raw deposit). The scripts still read the paths from params below; the
        # inputs exist for Snakemake's change-detection.
        unpack(lambda wc: {
            "matrix_dir": str(Path(config['deps']['cellranger']) / wc.capture),
            **{k: config['deps'][k] for k in ("demux", "ambient", "annotation", "samples")
               if config['deps'].get(k)},
        })
    output:
        h5ad = temp(out_dir / "per_capture_raw/{capture}.h5ad")
    log:
        log_dir / "create_anndata/{capture}.log"
    params:
        assignment_root = config['deps'].get('demux', ''),
        annotation_root = config['deps'].get('annotation', ''),
        ambient_root = config['deps'].get('ambient', ''),
        samples_file = config['deps'].get('samples', ''),
        modality = config.get('params', {}).get('modality', 'auto')
    # conda:
        # "../envs/scanpy.yaml"
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/create_anndata.py"
