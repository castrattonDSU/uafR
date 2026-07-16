categorate_enrichment_fixture = function(compounds, detail,
                                         query_overrides, ...) {
  cids = vapply(compounds, function(query) {
    selector = query_overrides[[query]]
    if (grepl("^cid:", selector)) {
      return(as.integer(sub("^cid:", "", selector)))
    }
    if (query == "salicin") return(439503L)
    if (query == "unresolved candidate") return(999001L)
    900000L + match(query, compounds)
  }, integer(1))
  list(
    PubChemIdentity = data.frame(
      Query = compounds,
      CID = cids,
      MatchStatus = "resolved_fixture",
      stringsAsFactors = FALSE
    ),
    PubChemProperties = data.frame(
      Query = compounds,
      CID = cids,
      MolecularFormula = paste0("C", seq_along(compounds), "H2"),
      stringsAsFactors = FALSE
    ),
    KEGGMatches = data.frame(
      Query = compounds,
      KEGG_ID = NA_character_,
      MatchStatus = "no_records",
      stringsAsFactors = FALSE
    ),
    ValidationSummary = data.frame(
      Status = "pass", ErrorCount = 0L, WarningCount = 0L,
      stringsAsFactors = FALSE
    )
  )
}

test_that("categorate enrichment batches preserve source query precedence", {
  out_dir = tempfile("categorate_enrichment_")
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  compounds = data.frame(
    compound_id = c("caffeine", "salicin", "candidate"),
    compound_name = c("caffeine", "salicin", "unresolved candidate"),
    CID = c(2519L, NA, NA),
    InChIKey = c(NA, "NGFMICBWJRZIBI-UHFFFAOYSA-N", NA),
    stringsAsFactors = FALSE
  )

  result = runCategorateEnrichmentBatches(
    compounds = compounds,
    out_dir = out_dir,
    cache = TRUE,
    batch_size = 2,
    periodic_cooldown_batches = Inf,
    periodic_cooldown_seconds = 0,
    enrichment_fun = categorate_enrichment_fixture,
    progress = FALSE
  )

  expect_s3_class(result, "uaf_categorate_enrichment_batches")
  expect_equal(result$QueryMap$QueryType,
               c("source_cid", "source_inchikey", "exact_name"))
  expect_equal(result$QueryMap$PubChemQuery,
               c("cid:2519", "NGFMICBWJRZIBI-UHFFFAOYSA-N",
                 "unresolved candidate"))
  expect_true(all(result$BatchManifest$status == "completed"))
  expect_equal(nrow(result$BatchSummary), 2)
  expect_equal(attr(result, "exit_status"), 0L)
  expect_true(all(file.exists(result$BatchManifest$batch_file)))

  reused = runCategorateEnrichmentBatches(
    compounds = compounds,
    out_dir = out_dir,
    cache = TRUE,
    batch_size = 2,
    periodic_cooldown_batches = Inf,
    periodic_cooldown_seconds = 0,
    enrichment_fun = function(...) stop("completed batch should be reused"),
    progress = FALSE
  )
  expect_true(all(reused$BatchManifest$status == "completed"))
  expect_true(all(reused$BatchManifest$checkpoint_status == "reused_valid"))
})

test_that("incomplete enrichment is not published as a successful checkpoint", {
  out_dir = tempfile("categorate_incomplete_")
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  incomplete = function(compounds, ...) {
    list(
      PubChemIdentity = data.frame(
        Query = compounds, CID = 2519L, MatchStatus = "resolved",
        stringsAsFactors = FALSE
      ),
      PubChemProperties = data.frame(),
      ValidationSummary = data.frame(Status = "warning", ErrorCount = 0L,
                                     stringsAsFactors = FALSE)
    )
  }
  result = runCategorateEnrichmentBatches(
    "caffeine", out_dir = out_dir, cache = TRUE,
    periodic_cooldown_batches = Inf, periodic_cooldown_seconds = 0,
    enrichment_fun = incomplete, progress = FALSE
  )
  expect_equal(result$BatchManifest$status, "failed")
  expect_equal(result$BatchManifest$checkpoint_status, "not_published")
  expect_false(file.exists(result$BatchManifest$batch_file))
  expect_equal(nrow(result$RetryQueue), 1)
  expect_equal(attr(result, "exit_status"), 1L)
})

test_that("service-busy enrichment pauses with retryable state", {
  out_dir = tempfile("categorate_busy_")
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  busy = function(...) {
    stop(structure(
      list(message = "PubChem service-busy circuit breaker opened after HTTP 503",
           call = NULL, status_code = 503L),
      class = c("uaf_pubchem_service_busy", "error", "condition")
    ))
  }
  paused = runCategorateEnrichmentBatches(
    c("caffeine", "salicin"), out_dir = out_dir, cache = TRUE,
    batch_size = 1, periodic_cooldown_batches = Inf,
    periodic_cooldown_seconds = 0, enrichment_fun = busy, progress = FALSE
  )
  expect_equal(paused$BatchManifest$status,
               c("paused_service_busy", "pending"))
  expect_equal(nrow(paused$RetryQueue), 2)
  expect_equal(attr(paused, "exit_status"), 75L)
  expect_equal(paused$RunSummary$Paused, "Yes")

  resumed = runCategorateEnrichmentBatches(
    c("caffeine", "salicin"), out_dir = out_dir, cache = TRUE,
    batch_size = 1, periodic_cooldown_batches = Inf,
    periodic_cooldown_seconds = 0,
    enrichment_fun = categorate_enrichment_fixture, progress = FALSE
  )
  expect_true(all(resumed$BatchManifest$status == "completed"))
  expect_equal(nrow(resumed$RetryQueue), 0)
})

test_that("PubChem source query overrides retain displayed compound names", {
  fetch = function(url) stop("CID overrides must not call the identity endpoint")
  overrides = .pubchem_query_override_map(c(caffeine = "cid:2519"),
                                           "caffeine")
  identity = .pubchem_resolve_cids("caffeine", fetch, overrides)
  expect_equal(identity$Query, "caffeine")
  expect_equal(identity$CID, 2519L)
  expect_equal(identity$QueriedName, "cid:2519")
  expect_equal(identity$MatchStatus, "resolved_source_cid")
})

test_that("conflicting structures for one compound name require review", {
  query_map = .categorate_enrichment_query_map(data.frame(
    compound_name = c("ambiguous compound", "ambiguous compound"),
    CID = c(1L, 2L),
    stringsAsFactors = FALSE
  ), require_source_identity = TRUE)
  expect_true(all(query_map$EnrichmentEligible == "No"))
  expect_true(all(query_map$ExclusionReason == "conflicting_identity_for_name"))
})
