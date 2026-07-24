# uafR 0.4.0.9000 (development)

## Major improvements

- PubChem record-mode enrichment now filters requested headings and sources
  while traversing PUG-View records rather than materializing every unrelated
  annotation first. Large LOTUS classification responses are parsed with
  preallocated, vectorized structures, eliminating quadratic post-processing
  while retaining every source-backed hierarchy returned for the compound.
- Resumable enrichment health checks now accept source-backed InChIKey/SMILES
  identities that have no PubChem CID. Source-provided CIDs must still resolve,
  every query must retain an explicit PubChem outcome, and properties remain
  required for every resolved CID.
- KEGG FIND queries now transliterate Greek characters, remove leading
  optical-rotation markers, and convert unsupported punctuation into safe
  keyword terms before requesting the KEGG REST API. Any remaining per-name
  HTTP 400 response is retained as an auditable rejected search candidate
  instead of failing an otherwise valid multi-compound enrichment batch.
- Production plant panels now bound rich `detail = "research"` enrichment to
  a configurable, deterministic evidence-priority sample (1,000 compounds by
  default). Selection first attempts one eligible structure-resolved compound
  per represented species, then fills remaining capacity by direct-species
  and source-support strength. Every deferred identity remains in the
  discovery result and exclusion audit, with a selection summary that reports
  species coverage; `Inf` remains an explicit, review-required opt-in for
  smaller runs.
- Hardened large PubChem identity and categorate enrichment runs. Identity
  batches now omit unused synonym requests, require complete properties for
  resolved CIDs, write manifest/retry files, and pause on retry-exhausted or
  service-busy responses instead of caching temporary failures as no-hits.
  Resumable categorate batches now default to one complete PUG-View record per
  CID with local heading/source filtering, avoiding the prior per-heading and
  per-source request multiplication while preserving the legacy filtered mode.
- Plant-panel discovery now defers fallback context and analysis assembly to
  the final merge. Mixed deferred and precomputed provider results derive
  context from the complete occurrence table before reconciling existing
  evidence, preventing direct-species context omissions and repeated work.
- Biological-context extraction now prefilters rows with actual context
  signals and applies best exact or conservative relaxed evidence through
  vectorized key matching. This avoids per-evidence data-frame copies during
  large plant-panel merges without changing source-backed context rules.
- Large plant-panel assembly now indexes query/species matches, occurrence
  counts, compound-property records, species summaries, and identity-review
  context instead of repeatedly scanning complete occurrence tables. The
  automatic discovery matrix is limited to the 5,000 most frequent features
  to prevent memory-heavy wide-table expansion; complete long-form occurrence
  and comparability tables remain available, and users can still request an
  explicit uncapped matrix with `plantPhytochemistryMatrix(max_traits = Inf)`.
- Full-panel validation now reports representative unresolved-identity
  examples instead of concatenating every unresolved name, and merged results
  perform one final validation pass rather than validating the same assembled
  object twice.
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
- Hardened large plant runs with versioned atomic discovery checkpoints,
  incremental chunk manifests, failed-query and retry queues, exact run
  signatures, corrupt-checkpoint rebuilding, incomplete-discovery enrichment
  deferral, live-run safety gates for 100+ species, and provider circuit
  breakers for repeated 429/503 responses. Offline acceptance tests now cover
  interruption and exact resume behavior across 705 synthetic species.
  `planPlantChemistryRun()` now also returns an input-name audit and blocks a
  species-level ready verdict when genus-only, `sp.`/`spp.`, blank, or
  duplicate names require curation. Provider and batch error diagnostics now
  redact credential-like URL/query values before printing or export.
- Reduced large-panel discovery recomputation with checkpoint schema `3.3.0`.
  Species chunks now checkpoint only normalized provider outputs; literature
  context, comparability, summaries, matrices, identity review, and validation
  are derived once after all chunks are combined. Exact duplicate literature
  records are collapsed independently of retrieval timestamps, while distinct
  provider records, citations, evidence sentences, and chemical mentions are
  preserved.
- Corrected production panel lifecycle and merge contracts. A successful or
  validated-reused stage now clears stale root failure/service-pause markers
  while retaining failure history in the progress ledger. Merging direct and
  fallback provider stages now emits one canonical `PlantQueries` row per
  submitted query with stable global query IDs instead of appending
  provider-local fallback query IDs.
- Strengthened production input provenance so preflight manifests distinguish
  the canonical plant input, distinct plant metadata, exclusion ledgers, and
  supporting taxonomy/normalization ledgers instead of labeling every
  secondary file as generic supporting material.
- Bound production panel runs to checked software artifacts. The panel runner
  can now require a clean-commit release-manifest JSON and matching source
  tarball; package version, Git commit, filename, byte count, MD5, and SHA-256
  are validated and included in the run signature. Both artifacts and portable
  source manifests are copied into the Tanimoto server handoff, whose commands
  install the exact tarball and validate the release manifest before live work.
- Hardened live NCBI behavior for large plant panels. PubChem occurrence
  discovery now resolves exact NCBI Taxonomy IDs through the official Datasets
  v2 suggestion endpoint and rejects non-exact suggestions. Live JSON/XML
  requests use bounded exponential retries, do not cache provider error
  payloads, and classify retriable HTTP 5xx, timeout, and backend failures as
  resumable service-busy conditions. The production pilot now pauses with exit
  status 75 before cache replay when discovery is incomplete.
- Hardened PUG-View annotation parsing for information records with no usable
  value or markup payload. These source records now retain a stable explicit
  missing-value row instead of aborting multi-compound research enrichment.
- Added deferred analysis assembly to `runPlantPhytochemistryBatch()` and the
  production panel provider stages. Large provider stages now write normalized
  discovery, diagnostics, identity, and resume artifacts without repeatedly
  rebuilding context, matrices, comparability, review, and filtered products;
  those products are built after provider results are merged.
- Hardened plant biological-context extraction by parsing each unique source
  record once and retaining sentence-local method distinctions. PubChem
  taxonomy rows no longer interpret CIDs as PMIDs, copy broad collection
  citation lists into occurrence evidence, or infer plant context from
  unrelated citation prose. Aggregate PubChem taxonomy associations remain
  review-required; exact local LOTUS, official NPASS, and exact-taxon
  KNApSAcK records can enter direct-evidence filters without being silently
  promoted beyond their source record.
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
- Corrected a critical legacy PubChem cache-key collision risk by replacing the
  former weighted-text hash with an MD5 key bound to the exact request URL.
  Legacy cache files are ignored automatically. PubChem fingerprint results are
  now checked against input InChIKeys, mismatches are excluded, verified aliases
  collapse to one canonical structure, and identity audit tables remain in the
  result. PubChem-derived identity and Tanimoto outputs created with the legacy
  cache implementation should be regenerated rather than reused.
- Added `preparePlantTanimotoInput()` and source-only behavior in
  `resolvePlantCompoundIdentities()`. The new offline preparation stage joins
  exact LOTUS source-record structures, keeps repeated evidence separate from
  unique plant-structure membership, removes inherited fingerprint payloads,
  excludes unresolved or review-required identities by default, and writes a
  checksummed, portable server handoff without computing any pairwise values.
- Added production server execution gates for plant Tanimoto handoffs.
  `tools/run_plant_tanimoto_server.R` now supports explicit `preflight`,
  `smoke`, `summary`, and `full` modes; deterministic plant-balanced smoke
  samples; durable status/progress/completion artifacts; release-manifest and
  disk checks; cache-preserving PubChem service-busy pauses; atomic staging;
  row-count validation; and compressed pair-output sharding.
  `chemicalTanimotoSimilarity()` now exposes optional progress callbacks,
  service-busy circuit breaking, and shard row limits while preserving prior
  defaults. Private release manifests can be generated only from a clean Git
  worktree after package tests, source checks, and clean-library installation.
- Added the canonical resumable `runPlantChemistryPanel()` workflow and
  installed `run_plant_chemistry_panel.R` CLI for staged preflight, local-index
  preparation, all-provider pilot/discovery, identity resolution, research and
  priority-full enrichment, Tanimoto handoff, and final bundle reconciliation.
  Production handoffs retain the complete plant universe, and server summaries
  now include explicit insufficient-support plant pairs plus exact comparable
  scope/group summaries.
- Hardened the all-provider pilot so discovery is record-bounded, the research
  gate accepts only source-backed CID/full-InChIKey identities, and its replay
  fails if any live request escapes the validated checkpoint/cache path.
  Interrupted pilot runs now regenerate disposable exports without deleting
  successful checkpoints.
- Hardened KEGG name resolution for production enrichment. Chemical names with
  commas use KEGG's documented `+` keyword form; only exact normalized KEGG
  synonyms enter biochemical tables; broad substring hits remain auditable in
  `KEGGSearchCandidates`; stale linked-entry 404s are cached as no-record; and
  429/503 or other transport failures propagate to resumable batch handling.
- Added `buildNpassIndex()`, `queryNpassIndex()`, provider availability
  manifests, and a lossless sharded NPASS 3.0/NPASS-2026 local adapter. The
  resource downloader validates remote byte counts, preserves interrupted
  partial transfers, resumes byte-range downloads, and records source and
  shard checksums. Temporary index shards use RDS rather than quoted text so
  provider fields cannot be truncated by delimiter/quote parsing.
- Added `tools/run_representative_lotus_pilot.R` for deterministic large-panel
  scale-up checks. It selects exact-species-key, genus-key/no-exact-key, and
  no-key strata; runs species-level discovery against a supplied local LOTUS
  index only; preserves explicit no-hit accounting; validates manifests and
  source-record identifiers; and writes a checksummed pilot bundle without
  starting the full panel or any live enrichment provider.
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
- Added production-hardening helpers for internal release readiness:
  `uafRWorkflowGuide()`, `uafRApiStability()`, `uafRSchemaMetadata()`,
  `uafRProviderContracts()`, `standardizeProviderDiagnostics()`,
  `inspectUafRCache()`, `summarizeUafRCache()`,
  `planPlantChemistryRun()`, `plantOccurrenceEvidenceDictionary()`,
  `chemistryComparisonDictionary()`, chemistry-classification override
  helpers, plant-evidence filter helpers, `uafRClaimGuidance()`, and
  `estimateTanimotoOutput()`. Plant chemistry bundle manifests now include
  schema/package metadata, and model-ready feature exports support count,
  binary, fraction, and confidence-weighted matrices.
- Added production release and recovery hardening with
  `validatePlantChemistryRunManifest()`, `writePlantChemistryRetryQueue()`,
  `rerunFailedPlantQueries()`, `standardizeCompoundIdentityAudit()`,
  `validateCompoundIdentityAudit()`,
  `exportCompoundIdentityReviewTemplate()`, and
  `applyCompoundIdentityReview()`. The release-check wrapper now installs the
  built source tarball into a clean temporary library by default and can run a
  offline student-bundle acceptance test when supplied a bundle path.
- Added a simulated offline plant chemistry example under
  `inst/extdata/offline_plant_chemistry` plus
  `tools/build_offline_plant_chemistry_example.R` so users and CI can generate
  a finalized plant chemistry bundle without live provider access.
- Hardened the production plant-chemistry path so manifest-backed LOTUS lookup
  directories work through `runPlantPhytochemistryBatch()`, feature matrices
  retain one ordered project species universe, nested feature-manifest paths
  are validated, and generic database provenance is no longer counted as
  biological plant-part/tissue or analytical-method context. Review-required
  exports now report separate evidence, identity, structure, context,
  comparability, and citation reasons with actionable next steps.
- Unified `chemistryComparisonDictionary()` and classification overrides with
  the comparison scopes, metabolism domains, biosynthetic families, behaviors,
  and groups emitted by the plant comparability engine. Common legacy aliases
  are normalized, while invalid or inconsistent override values now fail
  explicitly instead of entering comparable matrices.
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

- Replaced the legacy package vignette with a verified offline end-to-end
  workflow and expanded the pkgdown reference index to cover the full public
  API.
- Made finalized plant chemistry manifests portable across machines with
  bundle-relative artifact paths, file sizes, MD5 checksums, and checksum
  validation. External file references now omit source-machine absolute paths
  and document when large-file hashing is intentionally skipped.
- Established schema version `1.0.0`, development package version metadata,
  package citation guidance, and a conservative third-party provider-source
  notice for release preparation.
- Restored compatibility between `categorate()` output, bundled categorate data,
  `exactoThese()`, examples, and tests.
- Fixed external-standard calibration coefficient handling in `standardifyIt()`.
- Added explicit `spreadOut()` input validation.
- Converted live web-service tests to opt-in integration checks.
- Added cross-platform GitHub Actions checks for Linux, macOS, and Windows.
- Added GitHub/Bioconductor installation guidance and Bioconductor dependency
  metadata for `ChemmineR` and `fmcsR`.
