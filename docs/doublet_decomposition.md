# Doublet decomposition (exploratory, #51)

`rule doublet_decomposition` (`params.decomposition.enabled`) splits doublets into
**random co-encapsulation** and **cell aggregates** using the whole cohort. It writes
`doublets/decomposition/{per_capture,pair_enrichment,regression}.csv` and `report.html`.

## Model

Doublets come from two sources:

- **Random co-encapsulation.** This scales with the number of cells loaded. Its donor
  pairing follows the genotype formula: the cross-donor share is `1 − Σpᵢ²`, so the
  expected total rate is `d ≈ f_x / (1 − Σpᵢ²)` and the within-donor rate is
  `d_within = d − f_x` (swarbricklab/snp_demux#55).
- **Aggregates.** Incomplete dissociation gives aggregates that are always
  **same-sample**. Re-aggregation after dissociation can be **cross-sample**, more so
  when samples were dissociated together (`co_dissociated`) than when they were
  dissociated separately and pooled (`separate_then_pooled`). Aggregates should be
  enriched for **spatially adjacent** lineage pairs and for larger multiplets.

## Quantities

Per capture (`per_capture.csv`, engine `params.decomposition.engine`):

- `within_donor_calls`: expression-engine calls among **genotype singlets**.
- `expected_within_donor = n_singlets × d_within`: random co-encapsulation that genotypes can't see.
- `excess_within_donor = within_donor_calls − expected_within_donor`, and `excess_rate`
  per singlet. A positive excess points to same-sample aggregates, or to engine false
  positives (see the concordance report, #48).
- `f_x` comes from `doublet_summary.csv` (`f_x_posterior`, else `f_x_hard`), or
  otherwise falls back to `n_known / n_cells`.
- `cells_loaded` comes from `deps.capture_metadata` if present. Otherwise the
  called-cell count is used as a proxy, and `cells_loaded_source` records which.

Lineage-pair composition (`pair_enrichment.csv`, overall and per protocol):

- The **template** is the lineage pairs of cross-donor known doublets. These are
  random by construction, at least for cross-donor events.
- The **test** set is the lineage pairs of within-donor calls.
- `enrichment = share_within / share_cross`, with a pseudocount of 0.5. Pairs listed
  in `adjacent_pairs` are flagged `adjacent`. The list names the reference's major
  classes, so it comes from the lineage reference profile
  (`config/lineage_references/<name>.yaml`); `params.decomposition.adjacent_pairs`
  overrides it.
- Pairs are labelled as in the concordance report (`params.doublets.report_pair_source`):
  scDblFinder's most likely origin where available (needs `params.doublets.clusters`),
  else the classifier's top two major classes. The top-two fallback badly undercounts
  cross-lineage pairs on real data (#62: 6/53 genotype doublets vs ~82% expected), so
  run the decomposition with `clusters: lineage_fine`. `per_capture.csv` records
  `pair_source_used`.

Regression (`regression.csv`): OLS of `f_x ~ cells_loaded/1000 + dissociation_protocol`
across captures, with treatment coding and the alphabetically first protocol as
reference. Once loading is accounted for, a positive `co_dissociated` effect on the
cross-donor fraction indicates cross-sample re-aggregation.

## Assumptions and caveats

1. `d_within` is a **lower bound** (snp_demux#55): genotype methods miss some
   cross-donor doublets. So the excess is an upper bound on aggregate-driven calls.
2. Expression engines rarely detect same-lineage doublets (see #48 recall by
   category), so within-donor calls are biased towards cross-lineage pairs **whatever
   the mechanism**. That is why the comparison is against the cross-donor template
   *as seen through the same engine*. Even so, recall differs between pair types, so
   read enrichment as relative, not absolute.
3. Lineage pairs are estimates (scDblFinder's origin, or the classifier's runner-up
   class). When one parent dominates the UMIs the second parent is underestimated.
4. `cells_loaded` from the called-cell count is circular for loading effects. Supply
   real loading numbers in `deps.capture_metadata` (`capture_id, dissociation_protocol, cells_loaded`).
5. Engine false positives and real aggregates both raise `excess_within_donor`. Check
   the `singlet_flag` in the concordance report before reading excess as aggregates.
6. The cross-donor template is only as good as the genotype doublet calls. On the
   NeoMet_Pool_012 demo about half of the consensus doublets had only 1-2 of 5
   demuxer votes and singlet-like or lower UMI, i.e. look spurious
   (swarbricklab/snp_demux#65, swarbricklab/sc_pileup#29). Spurious "doublets" are
   singlets, so they pull the template towards same-lineage pairs, and inflate `f_x`
   and `d_within`. Resolve those issues before reading the enrichment or the
   regression.

## Inputs

- `deps.capture_metadata`: CSV with `capture_id`, `dissociation_protocol`
  (`separate_then_pooled` | `co_dissociated`, from wetlab notes) and optional `cells_loaded`.
- `doublets/<capture>.csv`, `doublets/known/<capture>.csv`, `lineages/<capture>.csv`.
- `<demux>/<capture>/doublet_summary.csv` (snp_demux#54/#55) for `f_x`, `d_within`.

## Not yet done

- Per-barcode sex calls (`sex_calls.csv`, swarbricklab/bcarc_chromium#34) as an extra
  signal for sex-discordant doublets.
- Multiplet size (triplets and larger), which would need a UMI/feature-count model.
