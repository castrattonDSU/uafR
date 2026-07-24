# Apply compound identity review decisions

User-facing wrapper around \`applyPlantCompoundIdentityReview()\` for
the production identity-review workflow.

## Usage

``` r
applyCompoundIdentityReview(
  x,
  review_table,
  reviewer = NA_character_,
  require_identity = TRUE
)
```

## Arguments

- x:

  Plant phytochemistry result, named list with \`CompoundResolution\`,
  or \`CompoundResolution\` table.

- review_table:

  Completed review table.

- reviewer:

  Optional reviewer name.

- require_identity:

  Logical. Require explicit proposed identity fields for update/replace
  decisions.

## Value

Updated object in the same shape as \`x\`.
