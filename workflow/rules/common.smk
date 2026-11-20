# workflow/rules/common.smk
# Common functions and utilities for the mkseurat workflow

def get_cellranger_path(wildcards):
    """Return the path to Cell Ranger output for a capture."""
    return Path(config['deps']['cellranger']) / wildcards.capture

def annotation_source(path):
    """A directory input of an annotation-stage rule (Cell Ranger matrix, demux root),
    marked ancient() (#49).

    Snakemake 8 treats an input newer than the outputs as updated unless its recorded
    checksum matches, and it never checksums directories. After `dvc pull`/`checkout`
    every source file gets a fresh mtime, so without ancient() the assemble stage
    (`snakemake all`) would recompute -- and overwrite -- the annotate stage's DVC
    outputs. When the source data really change, DVC re-runs the annotate stage and
    removes its outputs first; in single-stage use, delete the affected
    lineages/doublets outputs (or use -R) after changing Cell Ranger or demux inputs.
    """
    return ancient(str(path))

def get_assignment_path(wildcards):
    """Return the path to sample assignment file for a capture."""
    if 'demux' in config['deps']:
        return Path(config['deps']['demux']) / wildcards.capture / "cell_assignment.tsv"
    return None

def get_annotation_path(wildcards):
    """Return the path to cell annotation file for a capture."""
    if 'annotation' in config['deps']:
        return Path(config['deps']['annotation']) / wildcards.capture / "cell_types.csv"
    return None

def get_ambient_path(wildcards):
    """Return the path to ambient profile file for a capture."""
    if 'ambient' in config['deps']:
        return Path(config['deps']['ambient']) / wildcards.capture / "ambient_summary.csv"
    return None

def has_optional_input(config_key):
    """Check if an optional input is configured."""
    return config_key in config.get('deps', {}) and config['deps'][config_key]


# ── Gene-space reconciliation (reconcile_gene_space stage) ───────────────────────
# All helpers default to the upstream behaviour when the optional `reconcile:` block
# is absent, so a config without it yields a byte-identical DAG (no-op guarantee).

def reconcile_enabled():
    """True when the optional gene-lift stage is configured and switched on."""
    return bool(config.get('reconcile', {}).get('enabled', False))

def reconcile_mode():
    """'early' (default; lift the merged object before processing) or 'late'
    (lift the processed object, preserving the published embedding)."""
    return config.get('reconcile', {}).get('apply_mode', 'early')

def reference_paths():
    """Paths of the three pinned gene references from `reconcile.references`."""
    refs = config.get('reconcile', {}).get('references', {})
    return {
        'source_gtf': refs.get('source_gtf', {}).get('path'),
        'target_gtf': refs.get('target_gtf', {}).get('path'),
        'stable_id_history': refs.get('stable_id_history', {}).get('path'),
    }

def _obj_path(wildcards, filename):
    """Resolve an object path under the right subset directory (or the root)."""
    if use_subsets:
        return out_dir / f"{wildcards.subset}/{filename}"
    return out_dir / filename

def reconcile_source_object(wildcards):
    """Input object of the reconcile rule: the merged object in early mode, the
    processed object in late mode."""
    name = "processed.h5ad" if reconcile_mode() == "late" else "merged.h5ad"
    return _obj_path(wildcards, name)

def process_merged_input(wildcards):
    """Input of process_merged_anndata: the reconciled object when the lift runs
    early, otherwise the stock merged object."""
    if reconcile_enabled() and reconcile_mode() == "early":
        return _obj_path(wildcards, "reconciled.h5ad")
    return _obj_path(wildcards, "merged.h5ad")

def format_cellxgene_input(wildcards):
    """Input of format_cellxgene: the reconciled object when the lift runs late,
    otherwise the stock processed object."""
    if reconcile_enabled() and reconcile_mode() == "late":
        return _obj_path(wildcards, "reconciled.h5ad")
    return _obj_path(wildcards, "processed.h5ad")


# ── Config-gated keep-list filter (filter_deposit stage) ─────────────────────────
# Absent `filter:` block -> helpers are no-ops and the filter targets/rules drop
# from the DAG, so a config without it is byte-identical to the stock pipeline.

def filter_cfg():
    """The optional top-level `filter:` block (empty dict when absent)."""
    return config.get('filter', {})

def filter_enabled():
    """True when the keep-list filter is configured and switched on."""
    return bool(filter_cfg().get('enabled', False))

def filter_output_name():
    """File name of the filtered deposit object (under the subset/root dir)."""
    return filter_cfg().get('output', 'filtered.raw.h5ad')


# ── Doublet detection (#40/#41) ──────────────────────────────────────────────────
# Doublet detection runs unless the raw-deposit switches skip it. Each engine in
# params.doublets.engines writes its own calls; they are combined into the
# persisted doublets/<capture>.csv.

def doublets_enabled():
    """False when params.skip_scrublet / params.raw_only drop doublet detection."""
    p = config.get('params', {})
    return not (p.get('skip_scrublet', False) or p.get('raw_only', False))

DOUBLET_ENGINES = ("scrublet", "scdblfinder")

def doublet_engines():
    """Engines from params.doublets.engines (default: scrublet only)."""
    engines = list(config.get('params', {}).get('doublets', {}).get('engines', ['scrublet']))
    unknown = [e for e in engines if e not in DOUBLET_ENGINES]
    if unknown or not engines:
        raise ValueError(f"params.doublets.engines={engines!r}: allowed {DOUBLET_ENGINES}")
    return engines

def scrublet_in_merge():
    """True when the merge reads the Scrublet-annotated per-capture objects
    (per_capture/) rather than the raw ones (per_capture_raw/)."""
    return doublets_enabled() and "scrublet" in doublet_engines()

def doublet_call_files():
    """The persisted per-capture doublet CSVs (empty when detection is off)."""
    if not doublets_enabled():
        return []
    return [out_dir / f"doublets/{c}.csv" for c in captures]

def known_doublets_enabled():
    """Genotype-detected doublets are extracted whenever SNP demux is configured."""
    return doublets_enabled() and bool(config['deps'].get('demux'))

def doublet_params():
    """params.doublets with defaults filled in."""
    d = dict(config.get('params', {}).get('doublets', {}))
    d.setdefault('seed', 1)
    d.setdefault('known_use', 'discard')
    if d['known_use'] not in ('discard', 'positive'):
        raise ValueError(f"params.doublets.known_use={d['known_use']!r}: use discard|positive")
    # Genes dropped from scDblFinder's expression features when known doublets are
    # training positives, so it learns "two cells" rather than "two donors"
    # (sex-chromosome genes, HLA).
    d.setdefault('positive_exclude_genes', [
        'XIST', 'TSIX', 'RPS4Y1', 'DDX3Y', 'UTY', 'KDM5D', 'EIF1AY', 'ZFY',
        'USP9Y', 'NLGN4Y', 'TXLNGY', 'PRKY'])
    d.setdefault('positive_exclude_patterns', ['^HLA-'])
    d.setdefault('dbr_source', 'loading')
    if d['dbr_source'] not in ('genotype', 'loading', 'none'):
        raise ValueError(f"params.doublets.dbr_source={d['dbr_source']!r}: use genotype|loading|none")
    d['dbr_per1k'] = {'standard': 0.008, 'ht': 0.004, **d.get('dbr_per1k', {})}
    # Lineage labels as scDblFinder clusters (#47); needs params.lineages.enabled.
    d.setdefault('clusters', 'none')
    d.setdefault('min_cluster_size', 50)
    if d['clusters'] not in ('none', 'lineage_fine', 'lineage_major'):
        raise ValueError(f"params.doublets.clusters={d['clusters']!r}: use none|lineage_fine|lineage_major")
    if d['clusters'] != 'none' and not config.get('params', {}).get('lineages', {}).get('enabled'):
        raise ValueError("params.doublets.clusters needs params.lineages.enabled: true")
    # Concordance report (#48)
    d.setdefault('report', False)
    d.setdefault('report_min_second_prob', 0.2)
    d.setdefault('report_ratio_flag', 2.0)
    d.setdefault('report_pair_source', 'scdblfinder_origin')
    if d['report'] and not config['deps'].get('demux'):
        raise ValueError("params.doublets.report needs deps.demux (genotype doublets)")
    return d

def doublet_report_files():
    """Per-capture and cohort concordance CSVs (#48), when enabled."""
    if not (doublets_enabled() and doublet_params()['report']):
        return []
    return ([out_dir / f"doublets/report/concordance_{c}.csv" for c in captures]
            + [out_dir / "doublets/report/concordance_cohort.csv"])

def scdblfinder_param_files():
    """Persisted per-capture scDblFinder parameter records (#43)."""
    if not (doublets_enabled() and 'scdblfinder' in doublet_engines()):
        return []
    return [out_dir / f"doublets/scdblfinder_params/{c}.csv" for c in captures]


# ── Per-capture lineage labels (#45) ─────────────────────────────────────────────
_LINEAGE_TRAIN_REQUIRED = ('reference', 'major_key', 'fine_key')
# Settings a reference profile may supply. They name the reference's classes, so they
# belong to the reference, not the workflow.
_LINEAGE_PROFILE_KEYS = ('model', 'crosswalk', 'agreement_collapse', 'adjacent_pairs')

def _load_lineage_profile(path):
    """Reference profile YAML (e.g. config/lineage_references/wu2021.yaml). Relative
    paths inside it (model, crosswalk) are relative to the profile file."""
    import os
    import yaml
    with open(path) as fh:
        prof = yaml.safe_load(fh) or {}
    unknown = set(prof) - set(_LINEAGE_PROFILE_KEYS) - {'name', 'description'}
    if unknown:
        raise ValueError(f"lineage reference profile {path}: unknown keys {sorted(unknown)}")
    base = Path(path).parent
    for k in ('model', 'crosswalk'):
        if prof.get(k) and not Path(prof[k]).is_absolute():
            prof[k] = os.path.normpath(base / prof[k])
    return prof

def lineage_cfg():
    """params.lineages with defaults (enabled: false), merged with the reference profile.

    The classifier comes from exactly one of:
      model: a pre-trained model directory (e.g. resources/lineage_classifier/wu2021-v1)
      train: a reference atlas to train from in this run (train_lineage_classifier)
    params.lineages.reference names a profile that supplies model, crosswalk and
    agreement_collapse for one reference; keys set in params.lineages win, and a
    profile's model is ignored when train is set. deps.lineage_crosswalk, if set,
    overrides the crosswalk.
    """
    d = dict(config.get('params', {}).get('lineages', {}))
    d.setdefault('enabled', False)
    if d.get('reference'):
        prof = _load_lineage_profile(d['reference'])
        for k in _LINEAGE_PROFILE_KEYS:
            if k == 'model' and d.get('train'):
                continue
            if k not in d and prof.get(k) is not None:
                d[k] = prof[k]
    if config.get('deps', {}).get('lineage_crosswalk'):
        d['crosswalk'] = config['deps']['lineage_crosswalk']
    d.setdefault('source', 'model')
    d.setdefault('crosswalk_label_column', 'coarse_cell_type')
    # Model classes folded onto crosswalk names before computing agreement (e.g. Wu
    # 2021's Cancer/Normal Epithelial -> Epithelial); comes with the reference profile.
    d['agreement_collapse'] = dict(d.get('agreement_collapse') or {})
    if d['source'] not in ('model', 'crosswalk', 'both'):
        raise ValueError(f"params.lineages.source={d['source']!r}: use model|crosswalk|both")
    if d.get('train'):
        t = dict(d['train'])
        missing = [k for k in _LINEAGE_TRAIN_REQUIRED if not t.get(k)]
        if missing:
            raise ValueError(f"params.lineages.train is missing {missing}")
        t.setdefault('version', Path(t['reference']).stem)
        d['train'] = t
    if d['enabled']:
        if d['source'] in ('model', 'both'):
            if d.get('model') and d.get('train'):
                raise ValueError("params.lineages: set either model (pre-trained) or train "
                                 "(train from a reference), not both")
            if not d.get('model') and not d.get('train'):
                raise ValueError(f"params.lineages.source={d['source']} needs a classifier: "
                                 "set params.lineages.reference, model or train")
        if d['source'] in ('crosswalk', 'both'):
            if not d.get('crosswalk'):
                raise ValueError(f"params.lineages.source={d['source']} needs a crosswalk: "
                                 "params.lineages.reference (profile), crosswalk, or "
                                 "deps.lineage_crosswalk")
            if not config.get('deps', {}).get('annotation'):
                raise ValueError(f"params.lineages.source={d['source']} needs deps.annotation")
    return d

def lineage_annotation_file(capture):
    """cellranger annotate labels for a capture, or None without deps.annotation. The
    file may not exist; then the capture's crosswalk labels are NA."""
    if not config.get('deps', {}).get('annotation'):
        return None
    return Path(config['deps']['annotation']) / capture / "cell_types.csv"

def lineage_model_dir():
    """Directory of the classifier used by annotate_lineages: the configured
    pre-trained model, or the one trained in this run under the results directory."""
    d = lineage_cfg()
    if d.get('train'):
        return out_dir / "lineage_classifier" / d['train']['version']
    return Path(d['model'])

def lineage_target_features():
    """Gene space the trained panel is restricted to: params.lineages.train.target_features,
    or the features.tsv.gz of the first capture's Cell Ranger matrix."""
    t = lineage_cfg()['train']
    if t.get('target_features'):
        return t['target_features']
    return str(Path(config['deps']['cellranger']) / captures[0] / "features.tsv.gz")

def lineages_enabled():
    return bool(lineage_cfg()['enabled'])

def lineage_files():
    """The persisted per-capture lineage CSVs (empty when disabled)."""
    if not lineages_enabled():
        return []
    files = [out_dir / f"lineages/{c}.csv" for c in captures]
    if lineage_cfg()['source'] == 'both':
        files += [out_dir / f"lineages/agreement/{c}.csv" for c in captures]
    return files


# ── Attach per-cell annotations to the merged AnnData (#50) ──────────────────────
# Each source attaches only when it exists (lineages enabled / doublet detection
# running / deps.somatifind_calls set). params.attach.<source>: false switches
# one off; with all three off the merged object is unchanged.
def attach_cfg():
    d = dict(config.get('params', {}).get('attach', {}))
    return {k: bool(d.get(k, True)) for k in ('lineages', 'doublets', 'somatifind')}

def merge_doublet_csvs_needed():
    return doublets_enabled() and (scrublet_in_merge() or attach_cfg()['doublets'])

def merge_lineage_csvs_needed():
    return lineages_enabled() and attach_cfg()['lineages']

def somatifind_root():
    root = config['deps'].get('somatifind_calls', '')
    return root if (root and attach_cfg()['somatifind']) else ''


# ── Doublet decomposition (#51; exploratory) ─────────────────────────────────────
def decomposition_cfg():
    d = dict(config.get('params', {}).get('decomposition', {}))
    d.setdefault('enabled', False)
    d.setdefault('engine', 'scdblfinder')
    d.setdefault('f_x_column', 'f_x_posterior')
    # Spatially adjacent lineage pairs, 'A|B' in any order. They name the
    # reference's major classes, so they come from the lineage reference profile
    # (params.lineages.reference -> adjacent_pairs); set here to override.
    if 'adjacent_pairs' not in d:
        d['adjacent_pairs'] = list(lineage_cfg().get('adjacent_pairs') or [])
    if d['enabled']:
        if not (doublets_enabled() and known_doublets_enabled() and lineages_enabled()):
            raise ValueError("params.decomposition needs doublet detection, deps.demux "
                             "and params.lineages.enabled")
        if d['engine'] not in doublet_engines():
            raise ValueError(f"params.decomposition.engine={d['engine']!r} is not in "
                             f"params.doublets.engines={doublet_engines()}")
    return d

def decomposition_files():
    if not decomposition_cfg()['enabled']:
        return []
    return [out_dir / f"doublets/decomposition/{f}"
            for f in ("per_capture.csv", "pair_enrichment.csv", "regression.csv", "report.html")]
