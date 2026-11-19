# Within-lineage doublet score (exploratory, #52)

Same-lineage within-donor doublets (for example T + T from one donor) can't be seen by
genotype demux and are nearly invisible to expression engines. Cross-donor
same-lineage doublets are confirmed examples of the same physical event, so they show
what such a doublet looks like in QC space relative to singlets of that lineage.

`rule within_lineage_score` (`params.within_lineage.enabled`) writes:

- `doublets/qc/<capture>.csv`: per-cell `log_total_counts`, `log_n_genes`,
  `pct_counts_mt`, `pct_counts_ribo`. These come from the Cell Ranger matrix in the
  annotation stage, and are z-scored **within capture** before fitting.
- `doublets/within_lineage/separability.csv`, per major lineage: `n_labelled`,
  `n_unlabelled`, cross-validated `auc_cv`, the PU label frequency `c_hat`,
  `est_prevalence_unlabelled` and the logistic coefficients.
- `doublets/within_lineage/per_capture.csv`, per capture × lineage: `n_singlets`,
  `predicted_within_same_lineage`, and with snp_demux#55 also
  `expected_within_same_lineage = n_singlets × d_within × share of L|L among cross-donor
  known doublets` and the ratio between the two.

## Method: positive-unlabelled learning

The genotype singlets of lineage L are **not** clean negatives, because they contain
exactly the within-donor doublets being counted. So this is PU learning
(Elkan & Noto, 2008):

- Labelled positives are cross-donor known doublets whose lineage pair is `L|L`
  (pair from the classifier's top two classes, as in #48). Unlabelled cells are the
  genotype singlets labelled L.
- Fit `g(x) = P(labelled | x)` by logistic regression with out-of-fold predictions.
- `c = mean g` over held-out labelled cells.
- Each unlabelled cell's doublet probability is `w = (1 − c)/c · g/(1 − g)`, clipped
  to [0, 1]. The per-capture estimate is the sum of `w`.

A threshold-based count corrected with the false-positive rate on the same singlets
(Rogan–Gladen) was tried first and rejected. The hidden doublets inflate the measured
false-positive rate, so that estimate is about 0 by construction. On the synthetic
fixture it gave 0–1.8 against a truth of 6–9 per capture and lineage.

## Assumptions and known bias

- **SCAR**: genotype detection picks doublets independently of their QC features.
  This is plausible, since cross- and within-donor doublets are the same event, but it
  breaks if donors differ systematically in RNA content.
- `auc_cv` (labelled vs unlabelled) is **conservative**, because the unlabelled set
  contains positives.
- **Upward floor.** When the classes overlap, `c` is underestimated and the estimate
  has a positive bias that grows as separability drops. In simulation (1,500 singlets
  per capture, 5 seeds), the estimate tracks the hidden count about 1:1, but from a
  floor:

  | separability (AUC) | hidden = 0 | 30 | 60 | 120 |
  |---|---|---|---|---|
  | 0.96 | 35 | 64 | 92 | 154 |
  | 0.81 | 315 | 350 | 389 | 459 |

  So **differences across captures** are informative, and **absolute levels are not**
  unless separability is very high. Compare against `expected_within_same_lineage`
  with that in mind.

## Decision rule

Promote the score to an engine output only if a lineage shows clearly useful
separation (say `auc_cv ≥ 0.9` on held-out known doublets) **and** the per-capture
estimates follow `expected_within_same_lineage` across captures. Otherwise record the
separability table and stop.
