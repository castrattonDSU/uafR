# uafR 0.2.0

## Major improvements

- Added research-grade PubChem and KEGG enrichment outputs to `categorate()`
  with `detail = "research"` and `detail = "full"`.
- Added normalized `Chemical*` analysis tables for traits, ontology mappings,
  matrices, evidence, reports, measurements, hazards, uses, taxonomy,
  occurrences, bioactivity, targets, potencies, pathway roles, and KEGG reaction
  participants.
- Added standardized measurement values, canonical units, relation fields,
  behavior bins, and `ChemicalMeasurementSummary`.
- Added `pubchemProfile()` and `keggProfile()` for source-specific enrichment
  workflows.
- Added data-dictionary, validation, table-quality, source-diagnostic, and
  export helpers for reproducible downstream analysis.

## Reliability and packaging

- Restored compatibility between `categorate()` output, bundled categorate data,
  `exactoThese()`, examples, and tests.
- Fixed external-standard calibration coefficient handling in `standardifyIt()`.
- Added explicit `spreadOut()` input validation.
- Converted live web-service tests to opt-in integration checks.
- Added cross-platform GitHub Actions checks for Linux, macOS, and Windows.
- Added GitHub/Bioconductor installation guidance and Bioconductor dependency
  metadata for `ChemmineR` and `fmcsR`.
