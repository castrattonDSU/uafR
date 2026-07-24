# Join plant chemistry summaries to metadata

Join plant chemistry summaries to metadata

## Usage

``` r
joinPlantChemistryMetadata(
  metadata,
  plant_chemistry,
  species_col = "species",
  join_level = "species"
)
```

## Arguments

- metadata:

  User metadata data frame.

- plant_chemistry:

  Plant phytochemistry result, summary table, or matrix.

- species_col:

  Species column in \`metadata\`.

- join_level:

  One of \`"species"\`, \`"genus"\`, or \`"family"\`.

## Value

Joined data frame.
