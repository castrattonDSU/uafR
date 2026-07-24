# Combine data frames from categorate batch results

\`combineCategorateTables()\` binds selected result tables across
successful categorate batches. Columns are aligned before binding, and
batch metadata is added by default so every row remains traceable to the
source \`.rds\` file.

## Usage

``` r
combineCategorateTables(
  categorate_batches,
  tables = NULL,
  include_batch_metadata = TRUE,
  include_empty = FALSE,
  min_property_ratio = 0.9,
  max_cell_chars = 30000
)
```

## Arguments

- categorate_batches:

  A batch directory path, object returned by
  \`readCategorateBatchDirectory()\`, list of categorate result objects,
  or character vector of \`.rds\` files.

- tables:

  Character vector of categorate table names to combine. If \`NULL\`, a
  curated analysis-ready set is used.

- include_batch_metadata:

  Logical. If \`TRUE\`, add \`CategorateBatchIndex\`,
  \`CategorateBatchName\`, and \`CategorateBatchFile\`.

- include_empty:

  Logical. If \`TRUE\`, include empty tables when present.

- min_property_ratio:

  Minimum acceptable ratio of \`PubChemProperties\` rows to resolved
  PubChem CIDs.

- max_cell_chars:

  Maximum characters retained in a single character cell.

## Value

A named list of combined data frames.

## Examples

``` r
if (FALSE) { # \dontrun{
combined = combineCategorateTables(
  "categorate_batches",
  tables = c("ChemicalTraits", "PubChemProperties")
)
combined$ChemicalTraits
} # }
```
