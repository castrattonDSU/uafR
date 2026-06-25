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
- Added duplicate occurrence evidence collapsing so repeated provider rows for
  the same plant, compound, source, and record are merged without losing useful
  plant-part, method, evidence-text, or curation details.
- Added `plantPhytochemistryReviewTable()` and
  `applyPlantPhytochemistryReview()` for human review of candidate,
  fallback, unresolved, or literature-derived plant-compound evidence before
  promotion into curated occurrence records.
- Added `resolvePlantCompoundIdentities()` for fast PubChem identity-only
  compound resolution and `runPlantPhytochemistryBatch()` for staged,
  resumable species-first discovery across larger plant lists with checkpoint
  files, analysis-ready exports, review-required exports, run manifests, and
  optional richer compound enrichment after filtering.
- Hardened plant compound identity resolution so normalized compound keys
  preserve chemically meaningful Greek-letter and plus/minus prefixes, while
  PubChem lookups retry conservative deterministic aliases and report
  alias-derived matches as `resolved_alias`.
- Added source-backed plant compound identity recovery from local LOTUS indexes
  so LOTUS SMILES, InChIKeys, formulas, and CIDs are used before PubChem name
  lookup. Ambiguous source structures and resolved rows with class/product-like
  labels are retained in `plantCompoundIdentityReviewTable()` instead of being
  silently accepted. Completed review worksheets can be replayed with
  `applyPlantCompoundIdentityReview()` for reproducible identity curation.
- Added core PubChem Fingerprint2D Tanimoto workflows with
  `chemicalTanimotoSimilarity()` and `plantChemicalTanimotoSimilarity()`.
  These functions compute source-labeled compound-compound similarity, plant or
  group pair summaries, and optional streamed cross-group compound-pair files
  for downstream analyses such as plant chemistry and phylogenetic similarity
  comparisons.
- Hardened PubChem enrichment for larger DSI-style plant chemistry runs by
  allowing explicit `cid:<PubChem CID>` query tokens in `pubchemProfile()`,
  adding adaptive live-request spacing, respecting PubChem throttling headers
  where available, using exponential backoff for `429`/`503` service-busy
  responses, falling back from failed bulk property requests to single-CID
  property requests, and adding quality-gated enrichment-only categorate
  batching with cooldowns to `tools/run_dsi_categorate_tanimoto.R`.
- Added `readCategorateBatchDirectory()`, `summarizeCategorateBatches()`,
  `combineCategorateTables()`, and `exportPlantChemistryAnalysisBundle()` to
  turn resumable categorate batch directories into auditable, analysis-ready
  CSV/XLSX bundles for plant chemistry, Tanimoto, and downstream modeling
  workflows.
- Hardened plant chemistry analysis bundles for publication handoff with
  RFC4180-style CSV serialization, `finalizePlantChemistryAnalysisBundle()`,
  `validatePlantChemistryAnalysisBundle()`, enriched plant-compound membership
  exports, source/validation summaries, species-level chemistry summaries,
  missing-species coverage tables, data dictionaries, README/methods text, and
  analysis-ready plant-pair support keys.
- Added reusable plant chemistry project handoff tools:
  `standardizePlantMetadata()`, `plantOccurrenceEvidenceGrade()`,
  `plantComparableTanimotoSummary()`, `exportPlantChemistryFeatureSet()`, and
  `runPlantChemistryProject()`. Finalized bundles now include evidence-grade
  summaries, review-required occurrence tables, species feature matrices, and
  optional comparable scope/group Tanimoto summaries for downstream modeling
  and phylogeny/chemistry analyses.
- Added a comparable-chemistry layer for plant workflows with
  `plantChemistryComparability()` and `plantComparableChemistryMatrix()` so
  primary metabolites, specialized metabolites, volatile-specialized chemistry,
  lipids/fatty acids, hormone signals, and unknown chemistry can be separated
  before downstream comparison.
- Expanded plant comparable-chemistry classification to prioritize
  source-backed class fields from normalized LOTUS, PubChem classification,
  KEGG, ontology, term, and derived-group tables before lower-confidence
  compound-name rules.
- Added `plantContextEvidence()` and context-aware comparable matrices so
  source-backed plant-part, tissue, and method evidence can be extracted,
  audited, and used to compare biologically comparable chemistry subsets.
- Added provider-specific context extraction and `plantProviderContextAudit()`
  so LOTUS, KNApSAcK, PubChem Taxonomy, PubMed, and PubTator context coverage,
  confidence, review burden, and parser-hardening needs can be audited by
  source before large plant runs.
- Added `enrichPlantContextEvidence()` for source-backed biological-context
  enrichment from provided literature/source text or optional cached PubMed
  PMID/DOI fetches. It fills plant-part, tissue, and method context only from
  source text that mentions the species/genus or compound in a chemical-source
  context, preserving extraction rules and review flags. PubMed context fetches
  now prioritize informative source text and references shared across more
  species/compound occurrences when `max_sources` limits live requests.
- Added `runPlantPhytochemistryPilot()` and
  `tools/run_plant_phytochemistry_pilot.R` to run small live-provider pilot
  panels, write compact species summaries, QA reports, review-needed tables,
  and context-aware comparable matrices before scaling to large plant lists.
- Added `plantPhytochemistryPilotPanel()` and
  `plantPhytochemistryQAReport()` for repeatable 15-species pilot panels and
  explicit readiness checks before larger production runs.
- Added live-capable plant provider adapters for KNApSAcK organism lookup,
  conservative LOTUS taxon-evidence parsing, PubChem taxonomy annotations, and
  PubTator candidate chemical co-mentions, all covered by mocked no-network
  tests.
- Expanded the PubChem Taxonomy plant adapter to read source-backed
  `consolidatedcompoundtaxonomy` external tables for metabolites, natural
  products, and food compounds, and added provider elapsed-time, timeout, and
  error-message diagnostics for large plant runs.
- Hardened live LOTUS species queries by streaming only a bounded prefix of the
  unpaged simple-search response, using optional `curl` transport when
  available, and capping live LOTUS reads to a conservative number of records
  per species by default so large LOTUS payloads cannot stall full plant runs.
- Added local LOTUS index support with `standardizeLotusIndex()`,
  `queryLotusIndex()`, and `lotus_index` arguments for plant resolver, batch,
  pilot, and CLI workflows. This is the recommended LOTUS path for medium and
  large plant panels because it avoids repeated unpaged live API payloads while
  preserving species, genus, and family evidence labels.
- Added `buildLotusIndex()` and `tools/build_lotus_index.R` so flat LOTUS
  exports can be converted into a compact, reusable, manifest-backed local
  index before large species-first plant runs.
- Added `tools/flatten_lotus_mongo_dump.py`, a Python standard-library
  streaming flattener for the official LOTUS MongoDB ZIP download. It writes
  an auditable flat compound-taxon-reference export, a compact local uafR LOTUS
  index, and an optional manifest-backed lookup directory without requiring
  MongoDB command-line tools.
- Added PubChem-only compound enrichment fallback for plant phytochemistry
  workflows when no `chemical_library` is supplied, resumable PubChem-only
  enrichment batching for longer species runs, and plant-level matrices that
  include normalized `ChemicalTraits` features.
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
