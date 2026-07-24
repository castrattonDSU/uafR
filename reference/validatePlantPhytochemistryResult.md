# Validate plant phytochemistry results

Validate plant phytochemistry results

## Usage

``` r
validatePlantPhytochemistryResult(x, strict = FALSE)
```

## Arguments

- x:

  A result returned by \`resolvePlantPhytochemistry()\` or a named list
  containing plant phytochemistry tables.

- strict:

  Logical. If \`TRUE\`, missing optional documented columns are reported
  as warnings.

## Value

A list with \`Summary\`, \`TableQuality\`, \`ProviderDiagnostics\`,
\`Issues\`, and \`DataDictionary\`.
