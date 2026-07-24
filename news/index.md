# Changelog

## uafR 0.2.0

### Major improvements

- Added research-grade PubChem and KEGG enrichment outputs to
  [`categorate()`](https://castrattonDSU.github.io/uafR/reference/categorate.md)
  with `detail = "research"` and `detail = "full"`.
- Added normalized `Chemical*` analysis tables for traits, ontology
  mappings, matrices, evidence, reports, measurements, hazards, uses,
  taxonomy, occurrences, bioactivity, targets, potencies, pathway roles,
  and KEGG reaction participants.
- Added standardized measurement values, canonical units, relation
  fields, behavior bins, and `ChemicalMeasurementSummary`.
- Added
  [`pubchemProfile()`](https://castrattonDSU.github.io/uafR/reference/pubchemProfile.md)
  and
  [`keggProfile()`](https://castrattonDSU.github.io/uafR/reference/keggProfile.md)
  for source-specific enrichment workflows.
- Added data-dictionary, validation, table-quality, source-diagnostic,
  and export helpers for reproducible downstream analysis.

### Reliability and packaging

- Restored compatibility between
  [`categorate()`](https://castrattonDSU.github.io/uafR/reference/categorate.md)
  output, bundled categorate data,
  [`exactoThese()`](https://castrattonDSU.github.io/uafR/reference/exactoThese.md),
  examples, and tests.
- Fixed external-standard calibration coefficient handling in
  [`standardifyIt()`](https://castrattonDSU.github.io/uafR/reference/standardifyIt.md).
- Added explicit
  [`spreadOut()`](https://castrattonDSU.github.io/uafR/reference/spreadOut.md)
  input validation.
- Converted live web-service tests to opt-in integration checks.
- Added cross-platform GitHub Actions checks for Linux, macOS, and
  Windows.
- Added GitHub/Bioconductor installation guidance and Bioconductor
  dependency metadata for `ChemmineR` and `fmcsR`.
