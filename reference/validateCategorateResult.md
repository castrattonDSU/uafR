# Validate enriched categorate results

\`validateCategorateResult()\` audits a \`categorate()\` result against
the uafR data dictionary. It checks expected tables, required columns,
basic column types, controlled values, duplicate keys, table
completeness, and source coverage diagnostics. The return value is
designed to be stored with results and inspected before downstream
analyses.

## Usage

``` r
validateCategorateResult(x, tables = NULL, strict = FALSE)
```

## Arguments

- x:

  A list returned by \`categorate()\`, preferably with \`detail =
  "research"\` or \`detail = "full"\`.

- tables:

  Optional character vector of table names to validate. If \`NULL\`, all
  dictionary tables relevant to the result are checked.

- strict:

  Logical. If \`TRUE\`, missing optional documented columns are reported
  as warnings.

## Value

A list with \`Summary\`, \`TableQuality\`, \`SourceDiagnostics\`,
\`Issues\`, and \`DataDictionary\` data frames.

## Examples

``` r
if (FALSE) { # \dontrun{
result = categorate(compounds, library_data, detail = "research")
audit = validateCategorateResult(result)
audit$Summary
audit$Issues
} # }
```
