# Apply reviewed plant compound identity decisions

Applies a completed \`plantCompoundIdentityReviewTable()\` worksheet to
a plant phytochemistry result or \`CompoundResolution\` table. This
function is intentionally conservative: it only updates identity fields
when a reviewer supplies explicit proposed CID, InChIKey, SMILES,
formula, or source fields. It does not infer spelling corrections, merge
ambiguous source structures, or promote class/product labels
automatically.

Supported \`review_decision\` values are: \`"needs_review"\`/\`"keep"\`
(leave unchanged), \`"accept_resolved"\` (append a review note without
changing identity fields), \`"update_identity"\` or
\`"replace_identity"\` (copy proposed identity fields and mark the row
resolved when at least one identity field is supplied), and \`"reject"\`
or \`"exclude"\` (mark the row unresolved with \`resolution_source =
"review_excluded"\`).

## Usage

``` r
applyPlantCompoundIdentityReview(
  x,
  review_table,
  reviewer = NA_character_,
  require_identity = TRUE
)
```

## Arguments

- x:

  Plant phytochemistry result, named list with \`CompoundResolution\`,
  or a \`CompoundResolution\` data frame.

- review_table:

  Completed table from \`plantCompoundIdentityReviewTable()\`.

- reviewer:

  Optional reviewer name used when \`reviewed_by\` is empty.

- require_identity:

  Logical. If \`TRUE\`, update/replace decisions require at least one
  proposed CID, InChIKey, SMILES, molecular formula, or compound name.

## Value

Updated plant phytochemistry result or \`CompoundResolution\` table.
