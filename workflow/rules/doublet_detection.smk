# workflow/rules/doublet_detection.smk
# Rules for detecting doublets per capture.
#
# Engines (params.doublets.engines): scrublet (detect_doublets; also adds the
# doublet columns to the per-capture object the merge reads) and scdblfinder
# (detect_doublets_scdblfinder; calls only).
#
# Each engine writes its own long-format calls to
#   doublets/engines/<engine>/<capture>.csv   (temp)
# and combine_doublet_calls row-binds them into the persisted
#   doublets/<capture>.csv
# (one row per called barcode per engine; schema in scripts/doublet_io.py).

rule detect_doublets:
    """
    Detect doublets using Scrublet and optionally cross-lineage markers.
    """
    input:
        h5ad = out_dir / "per_capture_raw/{capture}.h5ad"
    output:
        h5ad = temp(out_dir / "per_capture/{capture}.h5ad"),
        calls = temp(out_dir / "doublets/engines/scrublet/{capture}.csv")
    log:
        log_dir / "doublet_detection/{capture}.log"
    params:
        lineage_markers = config['deps'].get('lineage_markers', ''),
        scrublet_threshold = config.get('params', {}).get('scrublet_threshold', 0.2),
        cross_lineage_threshold = config.get('params', {}).get('cross_lineage_threshold', 0.7)
    # conda:
        # "../envs/scanpy.yaml"
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/detect_doublets.py"


rule detect_doublets_scdblfinder:
    """
    scDblFinder on every called cell of a capture (package defaults). Writes
    calls only: the merged object's predicted_doublet still comes from Scrublet
    while the two engines are compared (#41).
    """
    input:
        matrix_dir = lambda wc: Path(config['deps']['cellranger']) / wc.capture
    output:
        calls = temp(out_dir / "doublets/engines/scdblfinder/{capture}.csv")
    log:
        log_dir / "doublet_detection/scdblfinder_{capture}.log"
    params:
        seed = config.get('params', {}).get('doublets', {}).get('seed', 1)
    conda:
        "../envs/scdblfinder.yaml"
    script:
        "../scripts/detect_doublets_scdblfinder.R"


rule combine_doublet_calls:
    """
    Persist per-capture doublet calls from every configured engine as one
    long-format CSV (not temp), for downstream modules (somatifind,
    ambient_calibration) that need per-barcode calls without the merged object.
    """
    input:
        calls = lambda wc: expand(
            out_dir / "doublets/engines/{engine}/{capture}.csv",
            engine=doublet_engines(), capture=wc.capture
        )
    output:
        csv = out_dir / "doublets/{capture}.csv"
    log:
        log_dir / "doublet_detection/combine_{capture}.log"
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/combine_doublet_calls.py"
