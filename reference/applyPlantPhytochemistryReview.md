# Apply reviewed plant phytochemistry evidence decisions

Applies a completed review table created by
\`plantPhytochemistryReviewTable()\`. Supported decisions are
\`"needs_review"\`/\`"keep"\` (leave unchanged), \`"keep_candidate"\`
(mark as reviewed but keep candidate/fallback status),
\`"update_context"\` (copy proposed plant
part/tissue/method/confidence/source fields), \`"promote_curated"\`
(promote to \`manual_curated\` evidence with review provenance), and
\`"reject"\`/\`"exclude"\` (remove the occurrence row).

## Usage

``` r
applyPlantPhytochemistryReview(
  x,
  review_table,
  reviewer = NA_character_,
  require_citation = TRUE
)
```

## Arguments

- x:

  Plant phytochemistry result or normalized occurrence table.

- review_table:

  Completed review table.

- reviewer:

  Optional reviewer name used when \`reviewed_by\` is empty.

- require_citation:

  Logical. If \`TRUE\`, \`promote_curated\` requires either an existing
  evidence URL or \`proposed_citation_or_url\`.

## Value

Updated plant phytochemistry result or occurrence table.
