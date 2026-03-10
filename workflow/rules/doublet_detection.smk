# workflow/rules/doublet_detection.smk
# Rule for detecting doublets in individual AnnData objects per capture
#

rule detect_doublets:
    """
    Detect doublets using Scrublet and optionally cross-lineage markers.
    """
    input:
        h5ad = out_dir / "per_capture_raw/{capture}.h5ad"
    output:
        h5ad = temp(out_dir / "per_capture/{capture}.h5ad")
    log:
        log_dir / "doublet_detection/{capture}.log"
    params:
        lineage_markers = config['deps'].get('lineage_markers', ''),
        scrublet_threshold = config.get('params', {}).get('scrublet_threshold', 0.2),
        cross_lineage_threshold = config.get('params', {}).get('cross_lineage_threshold', 0.7)
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/detect_doublets.py"
