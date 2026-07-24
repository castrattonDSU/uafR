# Summarize standardized chemical measurements

\`chemicalMeasurementSummary()\` condenses \`ChemicalMeasurements\` into
one row per compound and property using standardized units where
possible. It is designed for filtering and plotting measurement behavior
without manually parsing source units or evidence text.

## Usage

``` r
chemicalMeasurementSummary(x)
```

## Arguments

- x:

  A categorate result containing \`ChemicalMeasurements\`, or a
  \`ChemicalMeasurements\` data frame.

## Value

A data frame with one row per compound-property combination.

## Examples

``` r
if (FALSE) { # \dontrun{
result = categorate(compounds, library_data, detail = "research")
chemicalMeasurementSummary(result)
} # }
```
