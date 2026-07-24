# Summarize normalized chemical traits

`chemicalTraitSummary()` produces compact per-compound summaries from a
`ChemicalTraits` table or a categorate result. Use it to see how many
traits were extracted by domain or source, which sources contributed
evidence, and which trait groups/values dominate each compound.

## Usage

``` r
chemicalTraitSummary(
  traits,
  by = c("type", "source", "compound"),
  min_confidence = 0,
  top_n = 8
)
```

## Arguments

- traits:

  A `ChemicalTraits` data frame or a categorate result list that
  contains `ChemicalTraits`.

- by:

  Summary level: `"type"` summarizes by `TraitType`, `"source"`
  summarizes by `SourceDatabase`, and `"compound"` summarizes all traits
  per compound.

- min_confidence:

  Minimum confidence score for included traits. Accepts a numeric score
  or `"low"`, `"medium"`, or `"high"`.

- top_n:

  Number of top trait groups and values to include in collapsed summary
  columns.

## Value

A data frame with counts, confidence summaries, source coverage, and top
trait groups/values.

## Examples

``` r
if (FALSE) { # \dontrun{
result = categorate(compounds, library_data, detail = "full")
chemicalTraitSummary(result)
chemicalTraitSummary(result, by = "source", min_confidence = "high")
} # }
```
