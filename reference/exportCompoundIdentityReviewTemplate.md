# Export a compound identity review template

Export a compound identity review template

## Usage

``` r
exportCompoundIdentityReviewTemplate(
  x,
  path = NULL,
  include_resolved = FALSE,
  overwrite = FALSE
)
```

## Arguments

- x:

  Plant phytochemistry result, \`CompoundResolution\`, or identity audit
  table.

- path:

  Optional CSV path to write.

- include_resolved:

  Logical. If \`FALSE\`, resolved/no-review rows are omitted.

- overwrite:

  Logical. If \`FALSE\`, an existing \`path\` is not replaced.

## Value

Review template data frame.
