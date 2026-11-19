# workflow/rules/annotate_lineages.smk
# Per-capture lineage labels (#45), written as persisted CSVs.
#
# Runs on the full per-capture object (every called cell, doublets included:
# downstream consumers filter on doublet calls), before doublet detection and
# cohort subsetting. Gated by params.lineages.enabled (default false).
#
# The classifier is either pre-trained (params.lineages.model) or trained in this
# run from a reference atlas (params.lineages.train -> train_lineage_classifier).

rule train_lineage_classifier:
    """
    Train the lineage classifier from params.lineages.train.reference with
    workflow/scripts/lineage_train.py (the code behind tools/train_lineage_classifier.py).
    Reruns only when the reference, target gene space or training settings change.
    """
    input:
        reference = lambda wc: lineage_cfg()['train']['reference'],
        target_features = lambda wc: lineage_target_features()
    output:
        model_json = out_dir / "lineage_classifier/{version}/model.json",
        model_npz = out_dir / "lineage_classifier/{version}/model.npz",
        metrics = out_dir / "lineage_classifier/{version}/metrics.csv",
        gene_space = out_dir / "lineage_classifier/{version}/gene_space.csv"
    log:
        log_dir / "train_lineage_classifier/{version}.log"
    params:
        train = lambda wc: lineage_cfg()['train'],
        out_dir = lambda wc, output: str(Path(output.model_json).parent)
    threads: 1
    resources:
        mem_mb = lambda wc: lineage_cfg()['train'].get('mem_mb', 32000),
        runtime = lambda wc: lineage_cfg()['train'].get('runtime', 120)
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/train_lineage_classifier.py"


rule annotate_lineages:
    """
    Apply the lineage classifier (pre-trained or trained in this run; see
    docs/lineage_classifier.md) to every called cell of a capture.
    """
    input:
        h5ad = out_dir / "per_capture_raw/{capture}.h5ad",
        model_json = lambda wc: str(lineage_model_dir() / "model.json"),
        model_npz = lambda wc: str(lineage_model_dir() / "model.npz")
    output:
        csv = out_dir / "lineages/{capture}.csv"
    log:
        log_dir / "annotate_lineages/{capture}.log"
    params:
        model_dir = lambda wc: str(lineage_model_dir())
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/annotate_lineages.py"
