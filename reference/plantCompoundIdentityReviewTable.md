# Create a plant compound identity review table

Builds an auditable worksheet for unresolved or ambiguous compound
identity rows from species-first plant phytochemistry workflows. The
table classifies likely failure modes, ranks high-impact rows by
occurrence count and species coverage, and recommends conservative next
actions. It does not change compound identities and does not infer
spelling corrections.

## Usage

``` r
plantCompoundIdentityReviewTable(
  x,
  occurrences = NULL,
  include_resolved = FALSE,
  min_query_count = 1
)
```

## Arguments

- x:

  Plant phytochemistry result, \`CompoundResolution\` data frame, or a
  list containing \`CompoundResolution\` and optionally
  \`PlantCompoundOccurrences\`.

- occurrences:

  Optional normalized occurrence table used to add source, species, and
  record-count context when \`x\` is only a resolution table.

- include_resolved:

  Logical. If \`FALSE\`, already resolved rows are omitted.

- min_query_count:

  Minimum \`query_count\` to include.

## Value

A data frame suitable for CSV export, review, and re-import into a
project-specific curation workflow.
