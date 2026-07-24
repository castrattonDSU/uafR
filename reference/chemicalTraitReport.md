# Build a researcher-facing chemical trait report

\`chemicalTraitReport()\` summarizes the normalized trait system into
one compact row per compound. It combines trait counts, ontology
domains, evidence coverage, source databases, high-confidence terms,
external IDs, and nearest ontology neighbors so results are easier to
inspect before filtering, clustering, modeling, or exporting.

## Usage

``` r
chemicalTraitReport(x, min_confidence = 0, top_n = 8, neighbor_count = 3)
```

## Arguments

- x:

  A categorate result, a \`ChemicalTraitOntology\` data frame, or a
  \`ChemicalTraits\` data frame.

- min_confidence:

  Minimum confidence score for included rows. Accepts a numeric score or
  \`"low"\`, \`"medium"\`, or \`"high"\`.

- top_n:

  Maximum number of terms, IDs, URLs, and shared keys to collapse into
  each report cell.

- neighbor_count:

  Maximum number of nearest ontology neighbors to include per compound.
  Set to \`0\` to skip neighbor summaries.

## Value

A data frame with one row per compound and compact analysis-ready
summaries.

## Examples

``` r
if (FALSE) { # \dontrun{
result = categorate(compounds, library_data, detail = "research")
report = chemicalTraitReport(result, min_confidence = "medium")
} # }
```
