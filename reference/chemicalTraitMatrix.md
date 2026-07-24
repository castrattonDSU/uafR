# Build a matrix from normalized chemical traits

`chemicalTraitMatrix()` converts a `ChemicalTraits` table, or a
`categorate(detail = "research")`/`categorate(detail = "full")` result,
into a wide matrix for filtering, clustering, ordination, heatmaps, and
model inputs. The long `ChemicalTraits` table remains the complete
evidence table; this helper lets users choose compact or source-specific
matrix views without rerunning web queries.

## Usage

``` r
chemicalTraitMatrix(
  traits,
  profile = c("core", "full", "bioactivity", "safety", "ecology", "kegg", "sensory",
    "biomedical"),
  mode = c("binary", "count", "confidence"),
  min_confidence = 0,
  max_traits = Inf
)
```

## Arguments

- traits:

  A `ChemicalTraits` data frame or a categorate result list that
  contains `ChemicalTraits`.

- profile:

  Matrix profile. `"core"` keeps compact cross-domain grouping traits.
  `"full"` keeps every eligible trait. Other profiles keep traits for
  one analysis domain: `"bioactivity"`, `"safety"`, `"ecology"`,
  `"kegg"`, `"sensory"`, or `"biomedical"`.

- mode:

  Matrix value mode. `"binary"` stores 0/1 presence, `"count"` stores
  trait-row counts, and `"confidence"` stores the maximum confidence
  score for each compound-trait pair.

- min_confidence:

  Minimum confidence score for included traits. Accepts a numeric score
  or `"low"`, `"medium"`, or `"high"`.

- max_traits:

  Maximum number of trait columns to include, ranked by prevalence and
  confidence.

## Value

A data frame with `Query`, `CID`, and one column per selected trait.

## Examples

``` r
if (FALSE) { # \dontrun{
result = categorate(compounds, library_data, detail = "full")
core_matrix = result$ChemicalTraitMatrix
bioactivity_matrix = chemicalTraitMatrix(result, profile = "bioactivity")
confidence_matrix = chemicalTraitMatrix(result$ChemicalTraits,
                                       profile = "full",
                                       mode = "confidence",
                                       min_confidence = "high")
} # }
```
