# Build plant phytochemistry matrices

Build plant phytochemistry matrices

## Usage

``` r
plantPhytochemistryMatrix(
  x,
  level = c("species", "genus", "family"),
  profile = c("core", "full", "ecology", "metabolism", "safety", "bioactivity"),
  mode = c("binary", "count", "confidence"),
  min_confidence = "medium",
  max_traits = Inf
)
```

## Arguments

- x:

  Plant phytochemistry result or normalized occurrence table.

- level:

  One of \`"species"\`, \`"genus"\`, or \`"family"\`.

- profile:

  One of \`"core"\`, \`"full"\`, \`"ecology"\`, \`"metabolism"\`,
  \`"safety"\`, or \`"bioactivity"\`.

- mode:

  One of \`"binary"\`, \`"count"\`, or \`"confidence"\`.

- min_confidence:

  Minimum confidence to include.

- max_traits:

  Maximum number of columns to retain.

## Value

Wide matrix-like data frame.
