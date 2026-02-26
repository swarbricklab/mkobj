# workflow/rules/merge.smk
# Rules for merging per-capture Seurat objects
#
# When subsets are configured, the merge step also applies cell-level
# filtering via the subset's samples.csv.  Per-capture objects are shared
# across subsets and built only once.

if use_subsets:
    rule merge_captures:
        """
        Merge per-capture Seurat objects for a named subset.
        Filters cells to the subset's samples.csv before merging.
        """
        input:
            rds_files = lambda wc: expand(
                out_dir / "per_capture/{capture}.rds",
                capture=subset_captures[wc.subset]
            )
        output:
            merged = out_dir / "{subset}/merged.qs"
        log:
            log_dir / "{subset}/merge_captures.log"
        params:
            samples_file = lambda wc: subset_samples[wc.subset]
        conda:
            "../envs/seurat.yaml"
        script:
            "../scripts/merge_captures.R"
else:
    rule merge_captures:
        """
        Merge all per-capture Seurat objects into a single merged object.
        Uses Seurat's merge function and joins layers for proper integration.
        """
        input:
            rds_files = expand(out_dir / "per_capture/{capture}.rds", capture=captures)
        output:
            merged = out_dir / "merged.qs"
        log:
            log_dir / "merge_captures.log"
        params:
            samples_file = config['deps'].get('samples', '')
        conda:
            "../envs/seurat.yaml"
        script:
            "../scripts/merge_captures.R"
