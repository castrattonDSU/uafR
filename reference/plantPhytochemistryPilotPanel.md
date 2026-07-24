# Recommended plant panels for phytochemistry pilot runs

Returns small species panels useful for calibrating species-first plant
phytochemistry discovery before scaling. The panels are not claims that
each plant has complete public chemistry coverage; they are designed to
include a mix of expected data depth, project relevance, and likely edge
cases.

## Usage

``` r
plantPhytochemistryPilotPanel(profile = c("general", "remediation"))
```

## Arguments

- profile:

  One of \`"general"\` or \`"remediation"\`.

## Value

Data frame with \`species\`, \`panel_role\`, \`expected_data_depth\`,
\`review_focus\`, and \`rationale\` columns.
