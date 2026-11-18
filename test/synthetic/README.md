# Synthetic test fixture

A small, fully synthetic input set for running the AnnData track locally (no
cluster, no DVC pull). Ground truth is known for every barcode, which makes it
useful for doublet and lineage work.

```bash
python test/synthetic/make_fixture.py            # writes test/synthetic/data/
snakemake --snakefile workflow/Snakefile --configfile config/test_synthetic.yaml \
    --cores 2 --use-conda test/synthetic/results/merged.h5ad
```

Without containers or conda, run inside any env that has the packages listed in
`workflow/envs/scanpy.yaml` plus `snakemake` (drop `--use-conda`).

| Capture | Samples        | In cohort (`samples.csv`) |
| ------- | -------------- | ------------------------- |
| CAP_A   | S1 (F), S2 (M) | S1, S2                    |
| CAP_B   | S3 (F), S4 (M) | S3 (S4 is a passenger)    |

Each capture has 700 cells across four fine types in two major lineages
(Epithelial: LumHR, Basal; T-cells: CD4, CD8), plus 56 artificial doublets (summed
counts). Cross-donor doublets are called `doublet` by the fake demux; within-donor
doublets look like singlets to it. `truth/<capture>.csv` records donor(s), fine and
major type and `is_doublet` for each barcode.

Unit tests for pure helpers live in `test/unit/` and run with `pytest test/unit`.
