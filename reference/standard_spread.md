# Custom Output from \`spreadOut()\`

A list containing all of the output from \`spreadOut(standard_data)\`.
Has all necessary elements for downstream functions, specifically
\`mzExacto()\`

## Usage

``` r
standard_spread
```

## Format

\## \`standard_spread\` A list with 7 data frames and 1 list.

- Area:

  Spread Out Area/Quantity of Chemical

- Compounds:

  Spread Out Name of Tentative Chemicals

- MZ:

  Spread Out Captured M/Z Values

- MatchFactor:

  Spread Out Accuracy of Tentative Matches

- RT:

  Spread Out Retention Times

- Mass:

  Spread Out Compound Mass

- rtBYmass:

  Spread Out ID Codes (Retention Time Pasted to Mass)

- webInfo:

  A list containing all published names, top m/z peaks, exact masses,
  and optimal retention times for every tentative chemical in
  \`standard_data\`
