# Compare chemicals by normalized traits

\`chemicalTraitSimilarity()\` computes pairwise similarity from
\`ChemicalTraits\`, a categorate result, or an already-built trait
matrix. It reports shared and distinct trait counts plus collapsed
shared/distinct matrix keys so researchers can quickly see why compounds
group together.

## Usage

``` r
chemicalTraitSimilarity(
  traits,
  profile = c("core", "full", "bioactivity", "safety", "ecology", "kegg", "sensory",
    "biomedical"),
  min_confidence = 0,
  max_traits = Inf,
  top_n = 12
)
```

## Arguments

- traits:

  A \`ChemicalTraits\` data frame, a categorate result containing
  \`ChemicalTraits\`, or a trait matrix with \`Query\`/\`CID\` columns.

- profile:

  Trait matrix profile used when \`traits\` is not already a matrix. See
  \`chemicalTraitMatrix()\`.

- min_confidence:

  Minimum confidence score for included traits. Accepts a numeric score
  or \`"low"\`, \`"medium"\`, or \`"high"\`.

- max_traits:

  Maximum number of trait columns to include before computing
  similarity.

- top_n:

  Number of shared/distinct matrix keys to include in collapsed
  explanatory columns.

## Value

A pairwise data frame with Jaccard and overlap similarities.

## Examples

``` r
if (FALSE) { # \dontrun{
result = categorate(compounds, library_data, detail = "full")
chemicalTraitSimilarity(result)
chemicalTraitSimilarity(result, profile = "bioactivity")
} # }
```
