# Read a directory of resumable categorate batch files

\`readCategorateBatchDirectory()\` loads \`.rds\` files produced by a
resumable \`categorate()\` or categorate-enrichment batch run. It keeps
successful result objects and saved error objects together, then adds a
batch summary so downstream analysis can detect failed or incomplete
batches before combining tables.

## Usage

``` r
readCategorateBatchDirectory(
  path,
  pattern = "categorate_batch_[0-9]+[.]rds$",
  min_property_ratio = 0.9
)
```

## Arguments

- path:

  Directory containing batch \`.rds\` files.

- pattern:

  Regular expression used to identify batch files.

- min_property_ratio:

  Minimum acceptable ratio of \`PubChemProperties\` rows to resolved
  PubChem CIDs for a successful batch.

## Value

A list with class \`"uaf_categorate_batch_directory"\` containing
\`Batches\`, \`BatchSummary\`, and \`Source\`.

## Examples

``` r
if (FALSE) { # \dontrun{
batches = readCategorateBatchDirectory("categorate_batches")
batches$BatchSummary
} # }
```
