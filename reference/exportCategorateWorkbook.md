# Export enriched categorate results for review and sharing

\`exportCategorateWorkbook()\` writes the most useful \`categorate()\`
result tables to a researcher-friendly export bundle. It always supports
a directory of CSV files with an \`ExportManifest\` table. If
\`openxlsx\` or \`writexl\` is installed, it can also write a
multi-sheet \`.xlsx\` workbook.

## Usage

``` r
exportCategorateWorkbook(
  x,
  path,
  tables = NULL,
  format = c("auto", "xlsx", "csv"),
  include_raw = FALSE,
  include_empty = TRUE,
  overwrite = FALSE,
  max_cell_chars = 30000
)
```

## Arguments

- x:

  A list returned by \`categorate()\`, preferably with \`detail =
  "research"\` or \`detail = "full"\`.

- path:

  Output path. For \`format = "csv"\`, this is a directory. For \`format
  = "xlsx"\`, this is an \`.xlsx\` file.

- tables:

  Optional character vector of result table names to export. If
  \`NULL\`, a curated analysis-ready set is exported.

- format:

  Export format: \`"csv"\`, \`"xlsx"\`, or \`"auto"\`. \`"auto"\` writes
  \`.xlsx\` when \`path\` ends in \`.xlsx\` and an Excel writer is
  installed; otherwise it writes a CSV bundle.

- include_raw:

  Logical. If \`TRUE\` and \`tables = NULL\`, export every data frame in
  \`x\`, including raw source tables. If \`FALSE\`, export the curated
  analysis and diagnostics tables.

- include_empty:

  Logical. If \`TRUE\`, include empty data frames so the bundle
  preserves expected schema. If \`FALSE\`, omit empty tables.

- overwrite:

  Logical. If \`TRUE\`, replace an existing output file or directory.

- max_cell_chars:

  Maximum characters retained in a single cell. Longer values are
  truncated to keep spreadsheets responsive.

## Value

A manifest data frame describing exported tables, row/column counts,
sheet names, file names, and output path.

## Examples

``` r
if (FALSE) { # \dontrun{
result = categorate(compounds, library_data, detail = "research")
manifest = exportCategorateWorkbook(result, "categorate_export",
                                    format = "csv", overwrite = TRUE)
manifest
} # }
```
