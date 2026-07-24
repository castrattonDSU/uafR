# Export plant phytochemistry result tables

Export plant phytochemistry result tables

## Usage

``` r
exportPlantPhytochemistryWorkbook(
  x,
  path,
  format = c("auto", "xlsx", "csv"),
  tables = NULL,
  preset = c("all", "analysis_ready"),
  include_empty = TRUE,
  overwrite = FALSE,
  max_cell_chars = 30000
)
```

## Arguments

- x:

  Plant phytochemistry result.

- path:

  Output directory for \`format = "csv"\` or \`.xlsx\` path for \`format
  = "xlsx"\`.

- format:

  One of \`"csv"\`, \`"xlsx"\`, or \`"auto"\`.

- tables:

  Optional table names to export.

- preset:

  One of \`"all"\` or \`"analysis_ready"\`. The latter filters to
  direct/curated analysis-ready occurrence evidence before export while
  preserving validation and provenance tables.

- include_empty:

  Logical. If \`TRUE\`, include empty schema tables.

- overwrite:

  Logical. If \`TRUE\`, replace existing output.

- max_cell_chars:

  Maximum characters retained in any character cell.

## Value

Export manifest data frame.
