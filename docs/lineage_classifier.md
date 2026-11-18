# Per-capture lineage classifier

`annotate_lineages` (#45) labels every called cell of each capture, before doublet
detection and cohort subsetting, with a classifier trained **once** on a reference
atlas. This page records the choice of model and how to (re)build it (#44).

## Requirements

- Per cell and per capture: no integration, no clustering across captures.
- Runs in the existing scanpy container: CPU only, no per-run training, minutes per capture.
- Two levels: **major** (Wu 2021 `celltype_major`) and **fine** (`celltype_minor`). Fine
  labels feed scDblFinder's `clusters` (#47), and major labels feed SComatic pseudobulks.
- Versioned, checksummed, DVC-tracked.

## Choice: multinomial logistic regression on a fixed gene panel

The approach is CellTypist's, implemented in `workflow/scripts/lineage_model.py`:

1. Panel = highly variable genes (Cell Ranger flavour) ∪ top t-test markers per fine
   class, restricted to genes present in the target Cell Ranger reference. HVGs alone can
   miss strong, high-mean markers.
2. Normalise each cell to 10k counts, `log1p`, z-scale with the reference mean/SD, clip at ±10.
3. One L2-regularised multinomial logistic regression per level (scikit-learn,
   training only).
4. The model is saved **framework-free**: coefficients, intercepts, means and SDs go in
   `model.npz`, and genes, classes, the fine→major map, metadata and the sha256 of
   `model.npz` go in `model.json`. Prediction is a numpy matrix product and softmax, so it
   doesn't depend on the scikit-learn version, and the checksum is checked on every load.
   Confidence is the softmax probability of the chosen class.

### Alternatives considered

| Option | Why not (for this step) |
| --- | --- |
| CellTypist package | Same model family, but it adds a dependency to the container and its pickled models are tied to the package version. We keep the method and drop the dependency. |
| scANVI / scArches | Needs a GPU for reasonable speed and per-query fine-tuning, which is per-run training. |
| SingleR / Azimuth | Need the full reference at run time (memory, I/O) and R. Azimuth is built around specific references. |
| kNN label transfer (scanpy `ingest`) | Needs the reference embedding at run time, and is sensitive to batch effects without integration. |
| Marker scoring (existing cross-lineage markers) | No fine level, thresholds have to be tuned per lineage, and there's no calibrated confidence. |
| Cohort Leiden + curated cluster map (`ebca` `sc/major_annotation`) | Cohort-level, so it can't run per capture before integration. |

## Building the Wu 2021 model

```bash
python tools/train_lineage_classifier.py \
    --reference <wu2021_miniatlas>.h5ad \
    --major-key celltype_major --fine-key celltype_minor \
    --target-features <refdata-gex-GRCh38-2020-A>/features.tsv.gz \
    --version wu2021-v1 --out resources/lineage_classifier/wu2021-v1
dvc add resources/lineage_classifier/wu2021-v1
git add resources/lineage_classifier/wu2021-v1.dvc
```

The reference must hold raw counts (`--counts-layer` if they are not in `X`). Fine
classes with fewer than `--min-cells-per-class` cells are folded into their major class.

Outputs next to the model:

- `metrics.csv`: held-out accuracy, balanced accuracy and median confidence per level
  (20% split stratified by fine label; the shipped model is then refit on all cells).
- `confusion_<level>.csv`: held-out confusion matrices.
- `gene_space.csv`: reference genes, target genes, genes matched by Ensembl id and by
  symbol only, unmatched genes, and the panel size.

### Gene space

Query matrices from `refdata-gex-GRCh38-2020-A` carry Ensembl ids and symbols. At
prediction time panel genes are matched by version-stripped Ensembl id first, then by
symbol, and a matrix column is never used for two panel genes. Missing panel genes are
zero-filled, and the count of matched genes is logged. If the reference only has gene
symbols, as the public Wu 2021 matrices do, every match is by symbol. Check the
`unmatched` count in `gene_space.csv`, since symbols that changed between annotation
releases fall out of the panel.

## Status

- [x] Model code, training CLI, unit tests (`test/unit/test_lineage_model.py`).
- [x] End-to-end check on the synthetic fixture (`test/synthetic/make_lineage_model.py`).
- [ ] Train on the Wu 2021 mini-atlas (cluster), record held-out accuracy, `dvc add`.
- [ ] Time prediction on one real capture.
