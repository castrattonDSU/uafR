# Create a plant phytochemistry evidence review table

Builds a review worksheet for candidate, fallback, literature, or
unresolved plant-compound occurrence evidence. The table is intended for
human review: rows are not promoted automatically. Fill
\`review_decision\` and optional proposed fields, then pass the
completed table to \`applyPlantPhytochemistryReview()\`.

## Usage

``` r
plantPhytochemistryReviewTable(
  x,
  occurrence_status = c("candidate", "taxon_fallback", "literature_reported",
    "unresolved"),
  include_analysis_ready = FALSE
)
```

## Arguments

- x:

  Plant phytochemistry result or normalized occurrence table.

- occurrence_status:

  Occurrence statuses to include. Use \`"all"\` to include every
  occurrence row.

- include_analysis_ready:

  Logical. If \`FALSE\`, rows already marked analysis-ready are omitted
  unless \`occurrence_status = "all"\`.

## Value

A data frame suitable for CSV export, editing, and re-import.
