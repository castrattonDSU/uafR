# Summarize resumable categorate batches

\`summarizeCategorateBatches()\` reports which batch files are usable,
which contain saved errors, and whether PubChem property coverage meets
the requested threshold. Use this before combining large enrichment
outputs.

## Usage

``` r
summarizeCategorateBatches(categorate_batches, min_property_ratio = 0.9)
```

## Arguments

- categorate_batches:

  A batch directory path, object returned by
  \`readCategorateBatchDirectory()\`, list of categorate result objects,
  or character vector of \`.rds\` files.

- min_property_ratio:

  Minimum acceptable ratio of \`PubChemProperties\` rows to resolved
  PubChem CIDs.

## Value

A data frame with one row per batch.

## Examples

``` r
if (FALSE) { # \dontrun{
summarizeCategorateBatches("categorate_batches")
} # }
```
