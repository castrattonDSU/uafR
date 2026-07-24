# Run a staged plant phytochemistry batch workflow

Runs species-first discovery in resumable species chunks, then
optionally resolves compounds through a staged identity or enrichment
pass. The default compound stage is \`"identity"\`, but panels of 100 or
more species should first be run with \`compound_resolution_profile =
"none"\` and a local LOTUS index. Identity and richer PubChem enrichment
should follow only after the discovery manifest and occurrence evidence
have been reviewed.

## Usage

``` r
runPlantPhytochemistryBatch(
  plants,
  sources = c("lotus", "knapsack", "npass", "pubchem", "pubmed", "pubtator"),
  taxon_fallback = c("species", "genus"),
  out_dir = NULL,
  cache_dir = NULL,
  species_chunk_size = 25,
  compound_resolution_profile = c("identity", "none", "research", "full"),
  occurrence_status = c("direct_reported", "curated_reported"),
  analysis_ready = TRUE,
  min_confidence = "medium",
  max_compounds_per_species = Inf,
  max_unique_compounds = Inf,
  chemical_library = NULL,
  cache = TRUE,
  throttle = 0.5,
  ncbi_email = Sys.getenv("NCBI_EMAIL", ""),
  ncbi_tool = Sys.getenv("NCBI_TOOL", "uafR"),
  ncbi_api_key = Sys.getenv("NCBI_API_KEY", ""),
  max_pubmed_records = 50,
  max_provider_records = max_pubmed_records,
  lotus_index = Sys.getenv("UAFR_LOTUS_INDEX", ""),
  provider_indexes = NULL,
  require_all_providers = FALSE,
  request_timeout = 30,
  provider_results = NULL,
  discovery_fun = NULL,
  enrichment_fun = NULL,
  pubchem_fun = NULL,
  request_fun = NULL,
  pubtator_request_fun = NULL,
  compound_request_fun = NULL,
  compound_batch_size = 50,
  resume = TRUE,
  progress = interactive(),
  overwrite = FALSE,
  defer_derived = FALSE,
  strict = FALSE,
  stop_on_error = FALSE,
  allow_large_live_run = FALSE,
  service_busy_pause_threshold = 1,
  refresh = FALSE,
  ...
)
```

## Arguments

- plants:

  Character vector or data frame of plant names.

- sources:

  Public provider names passed to \`resolvePlantPhytochemistry()\`.

- taxon_fallback:

  Fallback ranks passed to \`resolvePlantPhytochemistry()\`.

- out_dir:

  Optional output directory for checkpoints and CSV/JSON products.

- cache_dir:

  Cache directory. If omitted and \`out_dir\` is supplied, a \`cache\`
  subdirectory under \`out_dir\` is used.

- species_chunk_size:

  Number of species per discovery chunk.

- compound_resolution_profile:

  One of \`"identity"\`, \`"none"\`, \`"research"\`, or \`"full"\`.
  \`"identity"\` performs a fast PubChem identity pass; \`"research"\`
  and \`"full"\` run richer uafR compound enrichment on the filtered
  resolution set.

- occurrence_status:

  Occurrence statuses used to select compounds for resolution. Defaults
  to direct and curated reported records.

- analysis_ready:

  Optional logical used to select compounds for resolution.

- min_confidence:

  Minimum confidence used to select compounds for resolution.

- max_compounds_per_species:

  Optional cap on unique compounds selected per species for compound
  resolution.

- max_unique_compounds:

  Optional cap on total unique compounds selected for compound
  resolution.

- chemical_library:

  Optional library passed to rich compound enrichment.

- cache:

  Logical. If \`TRUE\`, provider requests and batch checkpoints are
  reused where possible.

- throttle:

  Seconds to wait between uncached provider requests.

- ncbi_email:

  Optional NCBI email.

- ncbi_tool:

  Optional NCBI tool name.

- ncbi_api_key:

  Optional NCBI API key.

- max_pubmed_records:

  Maximum PubMed records per plant.

- max_provider_records:

  Maximum occurrence rows retained per plant from non-PubMed providers.

- lotus_index:

  Optional local LOTUS index as a data frame, flat file path, or
  manifest-backed lookup directory. Passed to
  \`resolvePlantPhytochemistry()\` and used instead of the live LOTUS
  simple API when \`"lotus"\` is enabled.

- provider_indexes:

  Optional named list of local provider indexes.
  \`provider_indexes\$lotus\` and \`provider_indexes\$npass\` are
  supported. \`lotus_index\` remains a backward-compatible alias.

- require_all_providers:

  Logical. If \`TRUE\`, each requested provider must be configured and
  each completed discovery chunk must report a successful provider
  stage. Intended for audited production runs after pilot testing.

- request_timeout:

  Maximum seconds allowed for an uncached provider request.

- provider_results:

  Optional mocked or pre-fetched provider results.

- discovery_fun:

  Optional replacement for \`resolvePlantPhytochemistry()\`. This is
  primarily intended for deterministic tests and fully local provider
  workflows.

- enrichment_fun:

  Optional rich enrichment function for tests or cached workflows.

- pubchem_fun:

  Optional replacement for \`pubchemProfile()\`.

- request_fun:

  Optional provider request function.

- pubtator_request_fun:

  Optional PubTator request function.

- compound_request_fun:

  Optional request function passed to \`pubchemProfile()\` during
  compound resolution.

- compound_batch_size:

  Maximum number of unique compounds per compound resolution batch.

- resume:

  Logical. If \`TRUE\`, reuse discovery and compound checkpoints.

- progress:

  Logical. If \`TRUE\`, print simple progress messages.

- overwrite:

  Logical. If \`TRUE\`, replace existing files in \`out_dir\`.

- defer_derived:

  Logical. If \`TRUE\`, retain normalized discovery, diagnostics,
  identity, and operational tables while deferring context linking,
  summaries, matrices, comparability, identity review, and their
  filtered exports. This is intended for provider stages that will be
  merged before analysis; the default preserves the complete standalone
  result.

- strict:

  Logical passed to validation.

- stop_on_error:

  Logical. If \`TRUE\`, stop on the first failed discovery chunk;
  otherwise record the failed chunk and continue.

- allow_large_live_run:

  Logical. Large runs (100 or more species) require persistent
  checkpoints and, by default, local/injected discovery. Set this to
  \`TRUE\` only after a small live pilot and review of
  \`planPlantChemistryRun()\`.

- service_busy_pause_threshold:

  Number of consecutive discovery chunks reporting rate-limit or
  service-busy errors before the run stops scheduling new chunks. Use
  \`Inf\` to disable this automatic pause.

- refresh:

  Logical passed to provider adapters.

- ...:

  Additional arguments passed to rich compound enrichment.

## Value

A \`uaf_plant_phytochemistry\` result with additional
\`BatchRunManifest\`, \`BatchChunkManifest\`, \`FailedQueries\`, and
\`RetryQueue\` tables.
