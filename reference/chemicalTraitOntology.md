# Map normalized chemical traits to a controlled ontology

`chemicalTraitOntology()` converts extracted `ChemicalTraits` into a
conservative, source-backed ontology table. The function only maps
whitelisted `TraitType`/`TraitGroup` combinations whose values are
already discrete in `ChemicalTraits`; it does not infer new biology or
parse long descriptive text.

## Usage

``` r
chemicalTraitOntology(traits, min_confidence = 0, include_unmapped = FALSE)
```

## Arguments

- traits:

  A `ChemicalTraits` data frame or a categorate result list that
  contains `ChemicalTraits`.

- min_confidence:

  Minimum confidence score for included traits. Accepts a numeric score
  or `"low"`, `"medium"`, or `"high"`.

- include_unmapped:

  Logical. If `TRUE`, unmapped discrete traits are kept under
  `OntologyDomain = "unmapped"` for auditing. Defaults to `FALSE`.

## Value

A data frame with controlled ontology domains, groups, terms, external
identifiers when available, source trait keys, source evidence, and
confidence metadata.

## Examples

``` r
if (FALSE) { # \dontrun{
result = categorate(compounds, library_data, detail = "research")
ontology = chemicalTraitOntology(result, min_confidence = "medium")
} # }
```
