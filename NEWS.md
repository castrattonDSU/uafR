# uafR 0.2.0

## Major improvements

- Added a phase-one species-first plant phytochemistry resolver with
  `resolvePlantPhytochemistry()`, curated intake normalization, provider
  diagnostics, conservative PubMed/PubTator candidate handling, species
  summaries, matrices, validation, scoring, and CSV/XLSX export helpers.
- Added plant occurrence evidence classification, analysis-ready occurrence
  flags, biological-context grouping for plant part/tissue/method fields,
  evidence quality scores, `filterPlantPhytochemistryEvidence()`, and an
  `exportPlantPhytochemistryWorkbook(preset = "analysis_ready")` export mode.
- Added live-capable plant provider adapters for KNApSAcK organism lookup,
  conservative LOTUS taxon-evidence parsing, PubChem taxonomy annotations, and
  PubTator candidate chemical co-mentions, all covered by mocked no-network
  tests.
- Added PubChem-only compound enrichment fallback for plant phytochemistry
  workflows when no `chemical_library` is supplied, plus plant-level matrices
  that include normalized `ChemicalTraits` features.
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
