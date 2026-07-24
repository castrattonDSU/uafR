# Trace trait and ontology keys back to source evidence

`chemicalTraitEvidence()` returns an audit table for `ChemicalTraits`,
`ChemicalTraitOntology`, or a categorate result. It is designed to
explain where matrix or ontology keys came from: source trait keys,
source databases, evidence IDs, evidence text/URLs, extraction rules,
confidence scores, and external identifiers when available.

## Usage

``` r
chemicalTraitEvidence(
  x,
  keys = NULL,
  type = c("auto", "ontology", "trait"),
  min_confidence = 0
)
```

## Arguments

- x:

  A categorate result, a `ChemicalTraitOntology` data frame, or a
  `ChemicalTraits` data frame.

- keys:

  Optional vector of ontology keys, trait matrix keys, source trait
  keys, external IDs, or terms to keep. Use values such as
  `"safety__ghs_hazard_code__h319"`, `"hazard__hazard_code__h319"`,
  `"H319"`, or `"map00590"`.

- type:

  Evidence source. `"auto"` uses `ChemicalTraitOntology` when available
  and otherwise falls back to `ChemicalTraits`. `"ontology"` returns
  ontology-key evidence. `"trait"` returns source-trait evidence.

- min_confidence:

  Minimum confidence score for included evidence rows. Accepts a numeric
  score or `"low"`, `"medium"`, or `"high"`.

## Value

A data frame with analysis keys, source traits, source evidence,
external identifiers, extraction rules, and confidence metadata.

## Examples

``` r
if (FALSE) { # \dontrun{
result = categorate(compounds, library_data, detail = "research")
chemicalTraitEvidence(result, keys = "safety__ghs_hazard_code__h319")
chemicalTraitEvidence(result,
                      keys = "hazard__hazard_code__h319",
                      type = "trait")
} # }
```
