# workflow/rules/filter_deposit.smk
# Config-gated keep-list filter: turn the raw-deposit merged.h5ad into a
# curated subset (e.g. the analysed+control singlets for a public deposit),
# declaratively and reproducibly.
#
# Enabled by a top-level `filter:` block (filter.enabled: true). Three stages:
#   extract_obs   merged.h5ad -> merged_obs.csv          (obs only; cheap)
#   resolve_keep  obs CSV x disposition CSV -> keep/dropped/breakdown  (DuckDB)
#   subset_deposit  merged.h5ad[keep] -> <filter.output>.h5ad   (classic var enc)
#
# The keep_cells / dropped_cells / status_usage_breakdown CSVs are tracked
# outputs, so the deposit's inclusion criteria are versioned alongside it.
# The predicate is generic SQL over (obs columns x the joined disposition value
# column); see config/template.yaml and workflow/scripts/resolve_keep.py.


rule extract_obs:
    """Dump .obs of the merged object to a flat CSV (obs only, backed read)."""
    input:
        h5ad=lambda wc: _obj_path(wc, filter_cfg().get("input", "merged.h5ad")),
    output:
        obs_csv=out_dir / ("{subset}/merged_obs.csv" if use_subsets else "merged_obs.csv"),
    log:
        log_dir / ("{subset}/extract_obs.log" if use_subsets else "extract_obs.log"),
    params:
        obs_columns=filter_cfg().get("obs_columns"),
    # conda:
        # "../envs/filter.yaml"
    container:
        "docker://ghcr.io/swarbricklab/mkobj-filter@sha256:faaceb40950d0fe4877fce8e69067d53aacf8fbafc6909fcefb7430bc7675cb9"
    script:
        "../scripts/extract_obs.py"


rule resolve_keep:
    """LEFT JOIN obs x disposition in DuckDB; evaluate the predicate; emit the
    keep-list, the dropped-cells reasons and the status x value breakdown."""
    input:
        obs_csv=out_dir / ("{subset}/merged_obs.csv" if use_subsets else "merged_obs.csv"),
        disposition=filter_cfg().get("disposition_table", ""),
    output:
        keep=out_dir / ("{subset}/keep_cells.csv" if use_subsets else "keep_cells.csv"),
        dropped=out_dir / ("{subset}/dropped_cells.csv" if use_subsets else "dropped_cells.csv"),
        breakdown=out_dir / (
            "{subset}/status_usage_breakdown.csv" if use_subsets else "status_usage_breakdown.csv"
        ),
    log:
        log_dir / ("{subset}/resolve_keep.log" if use_subsets else "resolve_keep.log"),
    params:
        disposition_value_col=filter_cfg().get("disposition_value_col", "usage"),
        joins=filter_cfg().get("joins", [["sample_id", "sample_id"]]),
        aliases=filter_cfg().get("aliases", {}),
        disposition_filter=filter_cfg().get("disposition_filter", {}),
        predicate=filter_cfg().get("predicate"),
        exclude_capture_regex=filter_cfg().get("exclude_capture_regex"),
        exclude_column=filter_cfg().get("exclude_column", "capture_id"),
    # conda:
        # "../envs/filter.yaml"
    container:
        "docker://ghcr.io/swarbricklab/mkobj-filter@sha256:faaceb40950d0fe4877fce8e69067d53aacf8fbafc6909fcefb7430bc7675cb9"
    script:
        "../scripts/resolve_keep.py"


rule subset_deposit:
    """Subset the merged object to the keep-list; write the deposit object with
    classic var/obs index encoding + a gene-axis assertion. Optionally casts X
    to a compact integer dtype and gzip+shuffle-compresses it (x_dtype /
    compression_level)."""
    input:
        h5ad=lambda wc: _obj_path(wc, filter_cfg().get("input", "merged.h5ad")),
        keep=out_dir / ("{subset}/keep_cells.csv" if use_subsets else "keep_cells.csv"),
    output:
        deposit=out_dir / ("{subset}/" + filter_output_name() if use_subsets else filter_output_name()),
    log:
        log_dir / ("{subset}/subset_deposit.log" if use_subsets else "subset_deposit.log"),
    params:
        expected_n_vars=filter_cfg().get("expected_n_vars"),
        flag_low_expression_genes=filter_cfg().get("flag_low_expression_genes", False),
        min_cells_per_gene=filter_cfg().get("min_cells_per_gene", 3),
        x_dtype=filter_cfg().get("x_dtype"),
        compression_level=filter_cfg().get("compression_level"),
    # conda:
        # "../envs/filter.yaml"
    container:
        "docker://ghcr.io/swarbricklab/mkobj-filter@sha256:faaceb40950d0fe4877fce8e69067d53aacf8fbafc6909fcefb7430bc7675cb9"
    script:
        "../scripts/subset_deposit.py"
