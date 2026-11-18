# workflow/rules/annotate_lineages.smk
# Per-capture lineage labels (#45), written as persisted CSVs.
#
# Runs on the full per-capture object (every called cell, doublets included:
# downstream consumers filter on doublet calls), before doublet detection and
# cohort subsetting. Gated by params.lineages.enabled (default false).

rule annotate_lineages:
    """
    Apply the fixed reference lineage classifier (params.lineages.model; see
    docs/lineage_classifier.md) to every called cell of a capture.
    """
    input:
        h5ad = out_dir / "per_capture_raw/{capture}.h5ad",
        model_json = lambda wc: str(Path(lineage_cfg()['model']) / "model.json"),
        model_npz = lambda wc: str(Path(lineage_cfg()['model']) / "model.npz")
    output:
        csv = out_dir / "lineages/{capture}.csv"
    log:
        log_dir / "annotate_lineages/{capture}.log"
    params:
        model_dir = lambda wc: lineage_cfg()['model']
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/annotate_lineages.py"
