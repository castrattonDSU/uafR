# Summarize pilot readiness checks

Builds a compact QA report for a plant phytochemistry result or pilot
result. The report is intended for deciding whether a species panel is
ready for larger-scale runs, manual review, or provider-specific
hardening. It does not certify complete chemistry coverage.

## Usage

``` r
plantPhytochemistryQAReport(
  x,
  min_species_with_records_fraction = 0.5,
  min_analysis_ready_species_fraction = 0.25,
  min_context_known_fraction = 0.25,
  min_resolved_fraction = 0.5,
  max_review_occurrence_fraction = 1
)
```

## Arguments

- x:

  Plant phytochemistry result.

- min_species_with_records_fraction:

  Minimum fraction of species that should have at least one public or
  curated occurrence row.

- min_analysis_ready_species_fraction:

  Minimum fraction of species that should have at least one
  analysis-ready compound.

- min_context_known_fraction:

  Minimum overall context-known occurrence fraction.

- min_resolved_fraction:

  Minimum fraction of attempted compound identities that should resolve.

- max_review_occurrence_fraction:

  Maximum acceptable review-required occurrence fraction before the
  pilot is considered review-heavy.

## Value

QA report data frame with check, status, value, threshold, details, and
recommendation columns.
