# Run resumable categorate enrichment batches

Runs the PubChem, KEGG, source-annotation, normalized trait, and
validation portion of \`categorate()\` in restartable batches. It
intentionally bypasses the legacy structural-library/FMCS stage, which
is unnecessary when a plant workflow already has source-backed compound
identities. Successful batches are published atomically; failed or
service-busy batches remain retryable.

## Usage

``` r
runCategorateEnrichmentBatches(
  compounds,
  out_dir,
  cache_dir = NULL,
  detail = c("research", "full"),
  batch_size = 25,
  pubchem_throttle = 1.1,
  kegg_throttle = 0.5,
  assay_detail_limit = ifelse(detail[[1]] == "full", 10, 0),
  pubchem_annotation_mode = c("record", "filtered"),
  cache = TRUE,
  resume = TRUE,
  overwrite = FALSE,
  require_source_identity = FALSE,
  trait_matrix_profile = "core",
  trait_matrix_mode = "binary",
  trait_matrix_min_confidence = 0,
  trait_matrix_max_traits = Inf,
  request_fun = NULL,
  kegg_request_fun = NULL,
  enrichment_fun = NULL,
  periodic_cooldown_batches = 10,
  periodic_cooldown_seconds = 60,
  stop_on_error = FALSE,
  progress = interactive()
)
```

## Arguments

- compounds:

  Character vector or data frame of compounds. Data frames may include
  \`compound_name\`, \`CID\`, \`InChIKey\`, \`SMILES\`, and
  \`compound_id\`. Source-backed CID is preferred, followed by a full
  InChIKey, then the exact compound name.

- out_dir:

  Persistent batch output directory.

- cache_dir:

  Persistent provider cache directory. Defaults to \`file.path(out_dir,
  "cache")\`.

- detail:

  Enrichment level, \`"research"\` or \`"full"\`.

- batch_size:

  Number of unique compounds per batch.

- pubchem_throttle:

  Minimum delay between uncached PubChem requests.

- kegg_throttle:

  Minimum delay between uncached KEGG requests.

- assay_detail_limit:

  Maximum full-detail PubChem assay descriptions per compound. Research
  runs should normally use zero.

- pubchem_annotation_mode:

  PubChem PUG-View request strategy. \`"record"\` retrieves one full
  record per CID and filters locally, which is the production default
  for large batches. \`"filtered"\` retains the legacy per-heading and
  per-source request strategy.

- cache:

  Logical. Provider response caching should remain enabled for
  production runs.

- resume:

  Logical. Reuse validated successful batch files.

- overwrite:

  Logical. Replace incompatible or completed batch artifacts. Provider
  response caches are never deleted.

- require_source_identity:

  Logical. If \`TRUE\`, only rows with a source-backed CID or full
  InChIKey are eligible.

- trait_matrix_profile, trait_matrix_mode, trait_matrix_min_confidence,
  trait_matrix_max_traits:

  Trait-matrix settings passed to the categorate research enrichment
  layer.

- request_fun, kegg_request_fun:

  Optional mocked request functions.

- enrichment_fun:

  Optional complete batch enrichment function used for deterministic
  tests. It receives \`compounds\`, \`detail\`, and batch settings.

- periodic_cooldown_batches:

  Number of newly completed batches between cooldowns. Use \`Inf\` to
  disable.

- periodic_cooldown_seconds:

  Cooldown duration.

- stop_on_error:

  Logical. Stop immediately after a non-service error.

- progress:

  Logical. Print batch progress.

## Value

A list with class \`"uaf_categorate_enrichment_batches"\` containing
\`Batches\`, \`BatchSummary\`, \`BatchManifest\`, \`QueryMap\`,
\`FailedQueries\`, \`RetryQueue\`, \`RunSummary\`, and \`OutputPath\`.
