# workflow/rules/common.smk
# Common functions and utilities for the mkseurat workflow

def get_cellranger_path(wildcards):
    """Return the path to Cell Ranger output for a capture."""
    return Path(config['deps']['cellranger']) / wildcards.capture

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
    return d

def scdblfinder_param_files():
    """Persisted per-capture scDblFinder parameter records (#43)."""
    if not (doublets_enabled() and 'scdblfinder' in doublet_engines()):
        return []
    return [out_dir / f"doublets/scdblfinder_params/{c}.csv" for c in captures]


# ── Per-capture lineage labels (#45) ─────────────────────────────────────────────
_LINEAGE_TRAIN_REQUIRED = ('reference', 'major_key', 'fine_key')

def lineage_cfg():
    """params.lineages with defaults (enabled: false).

    The classifier comes from exactly one of:
      model: a pre-trained model directory (e.g. resources/lineage_classifier/wu2021-v1)
      train: a reference atlas to train from in this run (train_lineage_classifier)
    """
    d = dict(config.get('params', {}).get('lineages', {}))
    d.setdefault('enabled', False)
    if d.get('train'):
        t = dict(d['train'])
        missing = [k for k in _LINEAGE_TRAIN_REQUIRED if not t.get(k)]
        if missing:
            raise ValueError(f"params.lineages.train is missing {missing}")
        t.setdefault('version', Path(t['reference']).stem)
        d['train'] = t
    if d['enabled']:
        if d.get('model') and d.get('train'):
            raise ValueError("params.lineages: set either model (pre-trained) or train "
                             "(train from a reference), not both")
        if not d.get('model') and not d.get('train'):
            raise ValueError("params.lineages.enabled is true but neither "
                             "params.lineages.model nor params.lineages.train is set")
    return d

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
    return [out_dir / f"lineages/{c}.csv" for c in captures]
