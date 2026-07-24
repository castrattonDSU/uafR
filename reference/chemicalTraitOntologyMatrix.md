# Build a matrix from controlled chemical trait ontology terms

`chemicalTraitOntologyMatrix()` converts `ChemicalTraitOntology`, a
categorate result, or a `ChemicalTraits` table into a wide matrix whose
columns are controlled ontology keys. This gives researchers a more
stable analysis surface than source-specific trait names while
preserving source identifiers, raw evidence, and confidence metadata in
`ChemicalTraits` and `ChemicalTraitOntology`.

## Usage

``` r
chemicalTraitOntologyMatrix(
  ontology,
  mode = c("binary", "count", "confidence"),
  min_confidence = 0,
  max_terms = Inf
)
```

## Arguments

- ontology:

  A `ChemicalTraitOntology` data frame, a categorate result, or a
  `ChemicalTraits` data frame.

- mode:

  Matrix value mode. `"binary"` stores 0/1 presence, `"count"` stores
  ontology-row counts, and `"confidence"` stores the maximum confidence
  score for each compound-ontology pair.

- min_confidence:

  Minimum confidence score for included ontology rows. Accepts a numeric
  score or `"low"`, `"medium"`, or `"high"`.

- max_terms:

  Maximum number of ontology columns to include, ranked by prevalence
  and confidence.

## Value

A data frame with `Query`, `CID`, and one column per selected ontology
key.

## Examples

``` r
if (FALSE) { # \dontrun{
result = categorate(compounds, library_data, detail = "research")
ontology_matrix = chemicalTraitOntologyMatrix(result,
                                              min_confidence = "medium")
} # }
```
