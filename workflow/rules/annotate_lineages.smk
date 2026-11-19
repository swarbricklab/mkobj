# workflow/rules/annotate_lineages.smk
# Per-capture lineage labels (#45/#46), written as persisted CSVs.
#
# Runs on the full per-capture object (every called cell, doublets included:
# downstream consumers filter on doublet calls), before doublet detection and
# cohort subsetting. Gated by params.lineages.enabled (default false).
#
# params.lineages.source: model | crosswalk | both. The classifier is either
# pre-trained (params.lineages.model, usually from the reference profile) or trained
# in this run from a reference atlas (params.lineages.train -> train_lineage_classifier).

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


def _annotate_lineages_inputs(wc):
    cfg = lineage_cfg()
    # Source inputs only, not the temp per_capture_raw object (#49).
    inputs = {"matrix_dir": str(Path(config['deps']['cellranger']) / wc.capture)}
    if cfg['source'] in ('model', 'both'):
        inputs["model_json"] = str(lineage_model_dir() / "model.json")
        inputs["model_npz"] = str(lineage_model_dir() / "model.npz")
    if cfg['source'] in ('crosswalk', 'both'):
        inputs["crosswalk"] = cfg['crosswalk']
        ann = lineage_annotation_file(wc.capture)
        if ann is not None and ann.exists():
            inputs["annotation"] = str(ann)
    return inputs


rule annotate_lineages:
    """
    Label every called cell of a capture with the reference lineage classifier
    (docs/lineage_classifier.md) and/or the cellranger-annotate crosswalk.
    """
    input:
        unpack(_annotate_lineages_inputs)
    output:
        csv = out_dir / "lineages/{capture}.csv",
        **({"agreement": out_dir / "lineages/agreement/{capture}.csv"}
           if lineage_cfg()['source'] == 'both' else {})
    log:
        log_dir / "annotate_lineages/{capture}.log"
    params:
        source = lambda wc: lineage_cfg()['source'],
        model_dir = lambda wc: str(lineage_model_dir()) if lineage_cfg()['source'] != 'crosswalk' else '',
        annotation = lambda wc: str(lineage_annotation_file(wc.capture) or ''),
        label_column = lambda wc: lineage_cfg()['crosswalk_label_column'],
        collapse = lambda wc: lineage_cfg()['agreement_collapse']
    container:
        "docker://ghcr.io/swarbricklab/mkobj-scanpy@sha256:f53e2dd5b4a75a4230a95fcd9f3373ba312faffbd35f65d984d9c69e9c4e795d"
    script:
        "../scripts/annotate_lineages.py"
