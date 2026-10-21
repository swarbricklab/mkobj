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
