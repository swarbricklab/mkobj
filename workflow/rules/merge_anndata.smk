# workflow/rules/merge_anndata.smk
# Rules for merging per-capture AnnData objects
#
# When subsets are configured, the merge step also applies cell-level
# filtering via the subset's samples.csv.  Per-capture objects are shared
# across subsets and built only once.

if use_subsets:
    rule merge_anndata_captures:
        """
        Merge per-capture AnnData objects for a named subset.
        Filters cells to the subset's samples.csv before merging.
        """
        input:
            h5ad_files = lambda wc: expand(
                out_dir / "per_capture/{capture}.h5ad",
                capture=subset_captures[wc.subset]
            )
        output:
            merged = out_dir / "{subset}/merged.h5ad"
        log:
            log_dir / "{subset}/merge_anndata.log"
        params:
            samples_file = lambda wc: subset_samples[wc.subset]
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
            h5ad_files = expand(out_dir / "per_capture/{capture}.h5ad", capture=captures)
        output:
            merged = out_dir / "merged.h5ad"
        log:
            log_dir / "merge_anndata.log"
        params:
            samples_file = config['deps'].get('samples', '')
        conda:
            "../envs/scanpy.yaml"
        script:
            "../scripts/merge_anndata.py"
