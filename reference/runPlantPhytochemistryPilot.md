# Run a plant phytochemistry pilot workflow

\`runPlantPhytochemistryPilot()\` is a higher-level workflow for small
live provider pilots before scaling to hundreds of species. It runs
\`runPlantPhytochemistryBatch()\`, writes the full audit bundle, and
adds pilot-specific files that make provider coverage, review needs,
biological context, and comparable-chemistry matrices easy to inspect.

The pilot does not claim complete plant metabolomes. It separates direct
species evidence, fallback/candidate evidence, source-backed biological
context, and chemistry-comparison scope so users can decide which
records are suitable for downstream analysis.

## Usage

``` r
runPlantPhytochemistryPilot(
  plants,
  out_dir,
  sources = c("lotus", "knapsack", "npass", "pubchem", "pubmed", "pubtator"),
  taxon_fallback = c("species", "genus"),
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
  max_pubmed_records = 25,
  max_provider_records = 100,
  lotus_index = Sys.getenv("UAFR_LOTUS_INDEX", ""),
  request_timeout = 30,
  provider_results = NULL,
  enrichment_fun = NULL,
  pubchem_fun = NULL,
  request_fun = NULL,
  pubtator_request_fun = NULL,
  compound_request_fun = NULL,
  compound_batch_size = 50,
  matrix_feature = c("comparison_group", "biosynthetic_family", "chemical_behavior",
    "compound", "comparison_scope"),
  matrix_mode = c("binary", "count", "confidence"),
  max_matrix_features = Inf,
  resume = TRUE,
  progress = interactive(),
  overwrite = FALSE,
  strict = FALSE,
  stop_on_error = FALSE,
  refresh = FALSE,
  ...
)
```

## Arguments

- plants:

  Character vector or data frame of plant names.

- out_dir:

  Output directory for audit tables, pilot reports, and matrices.

- sources:

  Provider sources passed to \`runPlantPhytochemistryBatch()\`.

- taxon_fallback:

  Taxon fallback ranks passed to \`runPlantPhytochemistryBatch()\`.

- cache_dir:

  Cache directory. Defaults to \`file.path(out_dir, "cache")\`.

- species_chunk_size:

  Number of species per discovery chunk.

- compound_resolution_profile:

  One of \`"identity"\`, \`"none"\`, \`"research"\`, or \`"full"\`.

- occurrence_status:

  Occurrence statuses retained for compound resolution.

- analysis_ready:

  Logical filter for compound-resolution input.

- min_confidence:

  Minimum confidence for resolution and matrices.

- max_compounds_per_species:

  Maximum analysis-ready compounds resolved per species.

- max_unique_compounds:

  Maximum unique compounds resolved across the pilot.

- chemical_library:

  Optional library passed to rich enrichment.

- cache:

  Logical. If \`TRUE\`, use cache/checkpoint files.

- throttle:

  Seconds to wait between uncached provider requests.

- ncbi_email:

  Optional NCBI email.

- ncbi_tool:

  Optional NCBI tool name.

- ncbi_api_key:

  Optional NCBI API key. It is not written to output files.

- max_pubmed_records:

  Maximum PubMed records per plant.

- max_provider_records:

  Maximum non-PubMed provider records per plant.

- lotus_index:

  Optional local LOTUS index as a data frame, flat file path, or
  manifest-backed lookup directory. Passed to
  \`runPlantPhytochemistryBatch()\` and used instead of the live LOTUS
  simple API when \`"lotus"\` is enabled.

- request_timeout:

  Maximum seconds allowed for an uncached provider request.

- provider_results:

  Optional mocked or pre-fetched provider results.

- enrichment_fun:

  Optional rich enrichment function for tests or cached workflows.

- pubchem_fun:

  Optional replacement for \`pubchemProfile()\`.

- request_fun:

  Optional provider request function.

- pubtator_request_fun:

  Optional PubTator request function.

- compound_request_fun:

  Optional request function for compound resolution.

- compound_batch_size:

  Maximum compounds per compound-resolution batch.

- matrix_feature:

  Matrix feature family passed to \`plantComparableChemistryMatrix()\`.

- matrix_mode:

  Matrix mode passed to \`plantComparableChemistryMatrix()\`.

- max_matrix_features:

  Maximum feature columns retained per pilot matrix.

- resume:

  Logical. If \`TRUE\`, reuse checkpoints where available.

- progress:

  Logical. If \`TRUE\`, print progress messages.

- overwrite:

  Logical. If \`TRUE\`, replace existing output files.

- strict:

  Logical passed to validation.

- stop_on_error:

  Logical. If \`TRUE\`, stop on first failed chunk.

- refresh:

  Logical passed to provider adapters.

- ...:

  Additional arguments passed to rich enrichment.

## Value

A \`uaf_plant_phytochemistry\` result with \`PilotSummary\`,
\`PilotReviewNeeded\`, \`PilotMatrices\`, and \`PilotExportManifest\`
tables.
