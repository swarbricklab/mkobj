# workflow/rules/merge_anndata.smk
# Rules for merging per-capture AnnData objects
#
# When subsets are configured, the merge step also applies cell-level
# filtering via the subset's samples.csv.  Per-capture objects are shared
# across subsets and built only once.
#
# "Do even less" raw-deposit path: when params.skip_scrublet (or params.raw_only)
# is set, the merge reads the PRE-scrublet per-capture objects written by
# create_anndata (per_capture_raw/) directly, instead of the detect_doublets
# output (per_capture/).  detect_doublets then never enters the DAG and no
# scrublet columns are added — the merged object is a true raw aggregate.

_skip_scrublet = (
    config.get('params', {}).get('skip_scrublet', False)
    or config.get('params', {}).get('raw_only', False)
)
_per_capture_dir = "per_capture_raw" if _skip_scrublet else "per_capture"

if use_subsets:
    rule merge_anndata_captures:
        """
        Merge per-capture AnnData objects for a named subset.
        Filters cells to the subset's samples.csv before merging.
        """
        input:
            h5ad_files = lambda wc: expand(
                out_dir / (_per_capture_dir + "/{capture}.h5ad"),
                capture=subset_captures[wc.subset]
            )
        output:
            merged = out_dir / "{subset}/merged.h5ad"
        log:
            log_dir / "{subset}/merge_anndata.log"
        params:
            samples_file = lambda wc: subset_samples[wc.subset],
            min_counts = config.get('params', {}).get('min_counts'),
            max_counts = config.get('params', {}).get('max_counts'),
            min_features = config.get('params', {}).get('min_features'),
            max_features = config.get('params', {}).get('max_features'),
            max_pct_mt = config.get('params', {}).get('max_pct_mt'),
            max_pct_ribo = config.get('params', {}).get('max_pct_ribo'),
            filter_doublets = config.get('params', {}).get('filter_doublets', True),
            min_cells_per_gene = config.get('params', {}).get('min_cells_per_gene', 3),
            expected_n_vars = config.get('params', {}).get('expected_n_vars')
        conda:
            "../envs/scanpy.yaml"
        script:
            "../scripts/merge_anndata.py"
else:
    rule merge_anndata_captures:
        """
        Merge all per-capture AnnData objects into a single merged object.
        Uses anndata.concat for proper concatenation.
        """
        input:
            h5ad_files = expand(out_dir / (_per_capture_dir + "/{capture}.h5ad"), capture=captures)
        output:
            merged = out_dir / "merged.h5ad"
        log:
            log_dir / "merge_anndata.log"
        params:
            samples_file = config['deps'].get('samples', ''),
            min_counts = config.get('params', {}).get('min_counts'),
            max_counts = config.get('params', {}).get('max_counts'),
            min_features = config.get('params', {}).get('min_features'),
            max_features = config.get('params', {}).get('max_features'),
            max_pct_mt = config.get('params', {}).get('max_pct_mt'),
            max_pct_ribo = config.get('params', {}).get('max_pct_ribo'),
            filter_doublets = config.get('params', {}).get('filter_doublets', True),
            min_cells_per_gene = config.get('params', {}).get('min_cells_per_gene', 3),
            expected_n_vars = config.get('params', {}).get('expected_n_vars')
        conda:
            "../envs/scanpy.yaml"
        script:
            "../scripts/merge_anndata.py"
