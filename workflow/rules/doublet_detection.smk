# workflow/rules/doublet_detection.smk
# Rules for detecting doublets per capture.
#
# Engines (params.doublets.engines): scrublet (detect_doublets) and scdblfinder
# (detect_doublets_scdblfinder). Both write calls only; the merge restores the
# legacy Scrublet obs columns (doublet_score, predicted_doublet, ...) from
# doublets/<capture>.csv (#49).
#
# Each engine writes its own long-format calls to
#   doublets/engines/<engine>/<capture>.csv   (temp)
# and combine_doublet_calls row-binds them into the persisted
#   doublets/<capture>.csv
# (one row per called barcode per engine; schema in scripts/doublet_io.py).

# Annotation-stage rules read source inputs (Cell Ranger matrix, demux files),
# never the temp per_capture_raw objects: rebuilding those for the merge must
# not invalidate the persisted CSVs (#49).

rule detect_doublets:
    """
    Detect doublets using Scrublet and optionally cross-lineage markers, on
    every called cell of the capture.
    """
    input:
        matrix_dir = lambda wc: annotation_source(Path(config['deps']['cellranger']) / wc.capture)
    output:
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


rule extract_known_doublets:
    """
    Genotype-detected doublets per capture: snp_demux's doublet_posteriors.csv
    when present, otherwise the demux `status` column of the per-capture object.
    """
    input:
        unpack(lambda wc: {
            "matrix_dir": annotation_source(Path(config['deps']['cellranger']) / wc.capture),
            **({"demux": annotation_source(config['deps']['demux'])} if config['deps'].get('demux') else {}),
            **({"samples": config['deps']['samples']} if config['deps'].get('samples') else {}),
        })
    output:
        csv = out_dir / "doublets/known/{capture}.csv"
    log:
        log_dir / "doublet_detection/known_{capture}.log"
    params:
        demux_root = config['deps'].get('demux', ''),
        samples_file = config['deps'].get('samples', '')
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/extract_known_doublets.py"


rule scdblfinder_params:
    """
    Choose scDblFinder's expected doublet rate per capture (params.doublets.dbr_source)
    and record it, with the reason for any fallback.
    """
    input:
        unpack(lambda wc: {
            "matrix_dir": annotation_source(Path(config['deps']['cellranger']) / wc.capture),
            **({"captures_file": config['deps']['captures']}
               if config['deps'].get('captures') else {}),
            **({"lineages": str(out_dir / f"lineages/{wc.capture}.csv")}
               if doublet_params()['clusters'] != 'none' else {}),
        })
    output:
        csv = out_dir / "doublets/scdblfinder_params/{capture}.csv",
        clusters = temp(out_dir / "doublets/scdblfinder_clusters/{capture}.csv")
    log:
        log_dir / "doublet_detection/scdblfinder_params_{capture}.log"
    params:
        clusters = doublet_params()['clusters'],
        min_cluster_size = doublet_params()['min_cluster_size'],
        dbr_source = doublet_params()['dbr_source'],
        per1k = doublet_params()['dbr_per1k'],
        captures_file = config['deps'].get('captures', ''),
        demux_root = config['deps'].get('demux', '')
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/scdblfinder_params.py"


rule detect_doublets_scdblfinder:
    """
    scDblFinder on every called cell of a capture (package defaults). Writes
    calls only: the merged object's predicted_doublet still comes from Scrublet
    while the two engines are compared (#41).
    """
    input:
        unpack(lambda wc: {
            "matrix_dir": annotation_source(Path(config['deps']['cellranger']) / wc.capture),
            "params_csv": str(out_dir / f"doublets/scdblfinder_params/{wc.capture}.csv"),
            "clusters": str(out_dir / f"doublets/scdblfinder_clusters/{wc.capture}.csv"),
            **({"known": str(out_dir / f"doublets/known/{wc.capture}.csv")}
               if known_doublets_enabled() else {}),
        })
    output:
        calls = temp(out_dir / "doublets/engines/scdblfinder/{capture}.csv")
    log:
        log_dir / "doublet_detection/scdblfinder_{capture}.log"
    params:
        seed = doublet_params()['seed'],
        known_use = doublet_params()['known_use'],
        exclude_genes = doublet_params()['positive_exclude_genes'],
        exclude_patterns = doublet_params()['positive_exclude_patterns']
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
        unpack(lambda wc: {
            "calls": expand(out_dir / "doublets/engines/{engine}/{capture}.csv",
                            engine=doublet_engines(), capture=wc.capture),
            **({"known": str(out_dir / f"doublets/known/{wc.capture}.csv")}
               if known_doublets_enabled() else {}),
        })
    output:
        csv = out_dir / "doublets/{capture}.csv"
    log:
        log_dir / "doublet_detection/combine_{capture}.log"
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/combine_doublet_calls.py"


rule doublet_concordance:
    """
    Per-capture concordance of every engine with the genotype doublets, by
    lineage pair (#48). Enabled by params.doublets.report.
    """
    input:
        unpack(lambda wc: {
            "calls": str(out_dir / f"doublets/{wc.capture}.csv"),
            "known": str(out_dir / f"doublets/known/{wc.capture}.csv"),
            **({"lineages": str(out_dir / f"lineages/{wc.capture}.csv")}
               if lineages_enabled() else {}),
        })
    output:
        csv = out_dir / "doublets/report/concordance_{capture}.csv"
    log:
        log_dir / "doublet_detection/concordance_{capture}.log"
    params:
        demux_root = config['deps'].get('demux', ''),
        min_second_prob = doublet_params()['report_min_second_prob'],
        pair_source = doublet_params()['report_pair_source'],
        ratio_flag = doublet_params()['report_ratio_flag'],
        known_use = doublet_params()['known_use']
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/doublet_concordance.py"


rule doublet_concordance_cohort:
    """Cohort roll-up of the per-capture concordance reports (#48)."""
    input:
        reports = expand(out_dir / "doublets/report/concordance_{capture}.csv", capture=captures)
    output:
        csv = out_dir / "doublets/report/concordance_cohort.csv"
    log:
        log_dir / "doublet_detection/concordance_cohort.log"
    params:
        ratio_flag = doublet_params()['report_ratio_flag']
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/doublet_concordance_cohort.py"


rule doublet_decomposition:
    """
    Cohort-level decomposition of doublets into random co-encapsulation vs
    aggregates, with dissociation protocol as a covariate (#51; exploratory).
    Enabled by params.decomposition.enabled.
    """
    input:
        unpack(lambda wc: {
            "calls": [str(out_dir / f"doublets/{c}.csv") for c in captures],
            "known": [str(out_dir / f"doublets/known/{c}.csv") for c in captures],
            "lineages": [str(out_dir / f"lineages/{c}.csv") for c in captures],
            **({"capture_metadata": config['deps']['capture_metadata']}
               if config['deps'].get('capture_metadata') else {}),
        })
    output:
        per_capture = out_dir / "doublets/decomposition/per_capture.csv",
        pair_enrichment = out_dir / "doublets/decomposition/pair_enrichment.csv",
        regression = out_dir / "doublets/decomposition/regression.csv",
        html = out_dir / "doublets/decomposition/report.html"
    log:
        log_dir / "doublet_detection/decomposition.log"
    params:
        engine = lambda wc: decomposition_cfg()['engine'],
        adjacent_pairs = lambda wc: decomposition_cfg()['adjacent_pairs'],
        f_x_column = lambda wc: decomposition_cfg()['f_x_column'],
        min_second_prob = doublet_params()['report_min_second_prob'],
        demux_root = config['deps'].get('demux', '')
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/doublet_decomposition_cohort.py"
