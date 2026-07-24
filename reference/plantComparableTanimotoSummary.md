# Summarize plant-pair Tanimoto by comparable chemistry

\`plantComparableTanimotoSummary()\` summarizes plant-compound pair
Tanimoto rows after requiring both compounds in a pair to share the same
comparable chemistry scope or group. This keeps like-with-like
comparisons separate from the unfiltered whole-chemistry plant-pair
summary.

## Usage

``` r
plantComparableTanimotoSummary(
  plant_pair_tanimoto = NULL,
  plant_compound_pair_tanimoto = NULL,
  membership,
  include_unknown = FALSE,
  thresholds = c(0.5, 0.7, 0.85, 0.95),
  max_pair_rows = Inf
)
```

## Arguments

- plant_pair_tanimoto:

  Optional unfiltered plant-pair summary.

- plant_compound_pair_tanimoto:

  Data frame or CSV/CSV.GZ path containing plant-compound pair Tanimoto
  rows.

- membership:

  Enriched plant-compound membership table with \`comparison_scope\`,
  \`comparison_group\`, and \`comparable_for_matrix\`.

- include_unknown:

  Logical. If \`FALSE\`, unknown/non-comparable chemistry is excluded
  from filtered summaries.

- thresholds:

  Numeric Tanimoto thresholds to count.

- max_pair_rows:

  Maximum plant-compound pair rows to read from a file.

## Value

List with \`Overall\`, \`ScopeFiltered\`, \`GroupFiltered\`, and
\`Diagnostics\`.
