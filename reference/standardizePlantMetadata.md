# Standardize plant metadata for chemistry projects

\`standardizePlantMetadata()\` normalizes user-supplied plant metadata
for joins into plant chemistry bundles. It derives genus from the
supplied species name when genus is absent, but it does not infer family
or accepted names. Missing family and accepted-name coverage are
reported through status columns so downstream bundles stay transparent.

## Usage

``` r
standardizePlantMetadata(
  x,
  species_col = NULL,
  accepted_species_col = NULL,
  genus_col = NULL,
  family_col = NULL,
  group_col = NULL,
  role_col = NULL,
  metadata_source = "user_supplied"
)
```

## Arguments

- x:

  Data frame, CSV path, or character vector of species names.

- species_col:

  Optional species column name.

- accepted_species_col:

  Optional accepted species column name.

- genus_col:

  Optional genus column name.

- family_col:

  Optional family column name.

- group_col:

  Optional project group column name.

- role_col:

  Optional project role/label column name.

- metadata_source:

  Source label written to \`metadata_source\`.

## Value

Standardized plant metadata data frame.
