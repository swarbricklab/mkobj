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

## References and models

The classifier is not tied to one atlas. Each trained model is a directory under
`resources/lineage_classifier/<reference>-v<N>/`, tracked by DVC (remote `mkobj`,
`/g/data/a56/dvc/analysis/mkobj`), and a run picks one with `params.lineages.model`.
Swapping references is a config change; nothing in the model code knows about Wu 2021.

A run can also start from scratch: set `params.lineages.train` (a reference atlas and
its label columns; see `config/README.md`) with `model: null`, and the
`train_lineage_classifier` rule trains a model into `<results>/lineage_classifier/<version>/`
before `annotate_lineages` uses it. It is the same code as the CLI below
(`workflow/scripts/lineage_train.py`), so the same reference and settings give the same
model. Snakemake retrains only when the reference, gene space or settings change. A model
worth sharing can be copied to `resources/lineage_classifier/` and DVC-tracked.

| Model | Reference | Labels (major / fine) | Status |
| --- | --- | --- | --- |
| `wu2021-v1` | Wu et al. 2021 mini-atlas (`miniatlas_cellxgene_7.1.0.h5ad`, 100k cells, 26 donors) | `celltype_major` (9) / `celltype_minor` (29) | current default |
| `ebca-v1` | `ebca` atlas | to be decided | planned |

`model.json` records the reference file, the label keys, the training split and the
held-out metrics, so a labels CSV can always be traced back to its reference.

The class vocabulary changes with the reference. Settings that name classes (the
cellranger crosswalk, `agreement_collapse`, `adjacent_pairs`) belong to the reference,
not the workflow, and travel with the model as defaults.

To fetch the models after cloning: `dvc pull resources/lineage_classifier/`.

### Building a model

```bash
python tools/train_lineage_classifier.py \
    --reference <atlas>.h5ad \
    --major-key <major label column> --fine-key <fine label column> \
    --group-key donor_id \
    --target-features <refdata-gex-GRCh38-2020-A>/features.tsv.gz \
    --version <reference>-v<N> --out resources/lineage_classifier/<reference>-v<N>
dvc add resources/lineage_classifier/<reference>-v<N>
dvc push resources/lineage_classifier/<reference>-v<N>.dvc
git add resources/lineage_classifier/<reference>-v<N>.dvc
```

Add `--use-raw` when counts are in `.raw` (as in the cellxgene Wu 2021 object) and
`--gene-id-key` when Ensembl ids are not in `var.index`. Report the donor-held-out
metrics (`--group-key`) for a new cohort; a per-cell split shares donors between train and
test and flatters accuracy. Never overwrite a published version: train `-v<N+1>`.

The reference must hold raw counts (`--counts-layer` if they are not in `X`). Fine
classes with fewer than `--min-cells-per-class` cells are folded into their major class.

Outputs next to the model:

- `metrics.csv`: held-out accuracy, balanced accuracy and median confidence per level
  (20% of cells stratified by fine label, or 20% of `--group-key` groups such as donors;
  the shipped model is then refit on all cells).
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
- [x] `wu2021-v1` trained on the Wu 2021 mini-atlas and DVC-tracked. Held out
      (cells / donors): major 0.992 / 0.994, fine 0.882 / 0.861 accuracy. Fine-level
      confidence is overconfident; major-level confidence is calibrated.
- [x] Prediction on one real capture: under 1 s.
- [ ] `ebca` model.
