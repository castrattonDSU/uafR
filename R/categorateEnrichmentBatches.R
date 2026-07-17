# Resumable research/full enrichment without the legacy FMCS library stage.

#' Run resumable categorate enrichment batches
#'
#' @description
#' Runs the PubChem, KEGG, source-annotation, normalized trait, and validation
#' portion of `categorate()` in restartable batches. It intentionally bypasses
#' the legacy structural-library/FMCS stage, which is unnecessary when a plant
#' workflow already has source-backed compound identities. Successful batches
#' are published atomically; failed or service-busy batches remain retryable.
#'
#' @param compounds Character vector or data frame of compounds. Data frames may
#' include `compound_name`, `CID`, `InChIKey`, `SMILES`, and `compound_id`.
#' Source-backed CID is preferred, followed by a full InChIKey, then the exact
#' compound name.
#' @param out_dir Persistent batch output directory.
#' @param cache_dir Persistent provider cache directory. Defaults to
#' `file.path(out_dir, "cache")`.
#' @param detail Enrichment level, `"research"` or `"full"`.
#' @param batch_size Number of unique compounds per batch.
#' @param pubchem_throttle Minimum delay between uncached PubChem requests.
#' @param kegg_throttle Minimum delay between uncached KEGG requests.
#' @param assay_detail_limit Maximum full-detail PubChem assay descriptions per
#' compound. Research runs should normally use zero.
#' @param pubchem_annotation_mode PubChem PUG-View request strategy. `"record"`
#' retrieves one full record per CID and filters locally, which is the
#' production default for large batches. `"filtered"` retains the legacy
#' per-heading and per-source request strategy.
#' @param cache Logical. Provider response caching should remain enabled for
#' production runs.
#' @param resume Logical. Reuse validated successful batch files.
#' @param overwrite Logical. Replace incompatible or completed batch artifacts.
#' Provider response caches are never deleted.
#' @param require_source_identity Logical. If `TRUE`, only rows with a
#' source-backed CID or full InChIKey are eligible.
#' @param trait_matrix_profile,trait_matrix_mode,trait_matrix_min_confidence,trait_matrix_max_traits
#' Trait-matrix settings passed to the categorate research enrichment layer.
#' @param request_fun,kegg_request_fun Optional mocked request functions.
#' @param enrichment_fun Optional complete batch enrichment function used for
#' deterministic tests. It receives `compounds`, `detail`, and batch settings.
#' @param periodic_cooldown_batches Number of newly completed batches between
#' cooldowns. Use `Inf` to disable.
#' @param periodic_cooldown_seconds Cooldown duration.
#' @param stop_on_error Logical. Stop immediately after a non-service error.
#' @param progress Logical. Print batch progress.
#'
#' @return A list with class `"uaf_categorate_enrichment_batches"` containing
#' `Batches`, `BatchSummary`, `BatchManifest`, `QueryMap`, `FailedQueries`,
#' `RetryQueue`, `RunSummary`, and `OutputPath`.
#'
#' @export
runCategorateEnrichmentBatches = function(
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
    progress = interactive()) {
  detail = match.arg(detail)
  pubchem_annotation_mode = match.arg(pubchem_annotation_mode)
  out_dir = .categorate_enrichment_output_dir(out_dir)
  cache_dir = .uaf_first_non_empty_text(cache_dir,
                                        file.path(out_dir, "cache"))
  cache_dir = normalizePath(cache_dir, winslash = "/", mustWork = FALSE)
  if (!isTRUE(cache)) {
    warning("Production enrichment should use `cache = TRUE`; uncached runs are not restart-safe.",
            call. = FALSE)
  }
  batch_size = .categorate_enrichment_positive_integer(batch_size,
                                                        "batch_size")
  pubchem_throttle = .categorate_enrichment_delay(pubchem_throttle,
                                                   "pubchem_throttle")
  kegg_throttle = .categorate_enrichment_delay(kegg_throttle,
                                                "kegg_throttle")
  cooldown_seconds = .categorate_enrichment_delay(periodic_cooldown_seconds,
                                                   "periodic_cooldown_seconds")
  cooldown_batches = suppressWarnings(as.numeric(periodic_cooldown_batches[[1]]))
  if ((!is.finite(cooldown_batches) && !is.infinite(cooldown_batches)) ||
      cooldown_batches < 1) {
    stop("`periodic_cooldown_batches` must be a positive number or `Inf`.",
         call. = FALSE)
  }

  query_map = .categorate_enrichment_query_map(
    compounds, require_source_identity = require_source_identity
  )
  eligible = query_map[query_map$EnrichmentEligible == "Yes", , drop = FALSE]
  if (nrow(eligible) < 1) {
    stop("No compounds are eligible for categorate enrichment. Inspect the query map for missing or conflicting identities.",
         call. = FALSE)
  }
  queries = unique(eligible$Query)
  groups = split(queries, ceiling(seq_along(queries) / batch_size))
  batch_dir = file.path(out_dir, paste0(detail, "_batches"))
  dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  run_signature = .categorate_enrichment_run_signature(
    query_map, detail, batch_size, pubchem_throttle, kegg_throttle,
    assay_detail_limit, pubchem_annotation_mode,
    trait_matrix_profile, trait_matrix_mode,
    trait_matrix_min_confidence, trait_matrix_max_traits
  )
  manifest_path = file.path(out_dir, paste0(detail, "_batch_manifest.csv"))
  prior_manifest = .categorate_enrichment_prior_manifest(
    manifest_path, run_signature, resume, overwrite
  )
  batch_files = file.path(
    batch_dir,
    paste0("categorate_batch_", sprintf("%04d", seq_along(groups)), ".rds")
  )
  manifest = .categorate_enrichment_manifest(groups, batch_files,
                                              run_signature)
  manifest = .categorate_enrichment_merge_manifest(manifest, prior_manifest)
  .plant_atomic_write_csv(query_map,
                          file.path(out_dir, "categorate_query_map.csv"))
  .plant_atomic_write_csv(manifest, manifest_path)

  completed_this_run = 0L
  paused = FALSE
  pause_reason = NA_character_
  batches = vector("list", length(groups))
  for (i in seq_along(groups)) {
    batch_queries = groups[[i]]
    batch_map = eligible[eligible$Query %in% batch_queries, , drop = FALSE]
    checkpoint = .categorate_enrichment_read_checkpoint(
      batch_files[[i]], run_signature, i, batch_queries
    )
    if (isTRUE(resume) && !isTRUE(overwrite) && checkpoint$ok) {
      batches[[i]] = checkpoint$result
      manifest$status[[i]] = "completed"
      manifest$checkpoint_status[[i]] = "reused_valid"
      manifest$completed_at[[i]] = .uaf_first_non_empty_text(
        manifest$completed_at[[i]], .plant_timestamp()
      )
      .plant_atomic_write_csv(manifest, manifest_path)
      next
    }
    if (file.exists(batch_files[[i]]) && !isTRUE(overwrite) &&
        checkpoint$reason %in% c("signature_mismatch", "batch_mismatch")) {
      stop("Existing enrichment checkpoint is incompatible with this run: ",
           batch_files[[i]], ". Use a new output directory or `overwrite = TRUE`.",
           call. = FALSE)
    }

    started = Sys.time()
    manifest$status[[i]] = "running"
    manifest$started_at[[i]] = .plant_timestamp()
    manifest$attempt_count[[i]] = manifest$attempt_count[[i]] + 1L
    .plant_atomic_write_csv(manifest, manifest_path)
    if (isTRUE(progress)) {
      message("uafR categorate ", detail, " batch ", i, "/",
              length(groups), ": ", length(batch_queries), " compound(s)")
    }

    override_map = stats::setNames(batch_map$PubChemQuery, batch_map$Query)
    condition = NULL
    result = tryCatch({
      if (is.function(enrichment_fun)) {
        enrichment_fun(
          compounds = batch_queries,
          detail = detail,
          query_overrides = override_map,
          cache = cache,
          cache_dir = cache_dir,
          pubchem_throttle = pubchem_throttle,
          kegg_throttle = kegg_throttle,
          assay_detail_limit = assay_detail_limit,
          pubchem_annotation_mode = pubchem_annotation_mode
        )
      } else {
        .categorate_research_enrichment(
          compounds = batch_queries,
          data_list = .categorate_enrichment_empty_standard_tables(),
          detail = detail,
          cache = cache,
          cache_dir = cache_dir,
          throttle = pubchem_throttle,
          assay_detail_limit = assay_detail_limit,
          trait_matrix_profile = trait_matrix_profile,
          trait_matrix_mode = trait_matrix_mode,
          trait_matrix_min_confidence = trait_matrix_min_confidence,
          trait_matrix_max_traits = trait_matrix_max_traits,
          request_fun = request_fun,
          kegg_request_fun = kegg_request_fun,
          pubchem_query_overrides = override_map,
          pubchem_annotation_mode = pubchem_annotation_mode,
          kegg_throttle = kegg_throttle,
          strict_sources = TRUE
        )
      }
    }, error = function(error) {
      condition <<- error
      NULL
    })
    health = if (is.null(condition)) {
      .categorate_enrichment_batch_health(result, batch_queries, batch_map)
    } else {
      list(ok = FALSE, message = .plant_redact_secrets(conditionMessage(condition)))
    }
    if (!isTRUE(health$ok)) {
      error_message = .uaf_first_non_empty_text(
        health$message, "Enrichment batch failed validation."
      )
      service_busy = .plant_service_busy_message(error_message) ||
        inherits(condition, "uaf_pubchem_service_busy") ||
        inherits(condition, "uaf_pubchem_request_failed")
      manifest$status[[i]] = ifelse(service_busy, "paused_service_busy",
                                    "failed")
      manifest$error_message[[i]] = error_message
      manifest$finished_at[[i]] = .plant_timestamp()
      manifest$elapsed_seconds[[i]] = round(
        as.numeric(difftime(Sys.time(), started, units = "secs")), 3
      )
      manifest$checkpoint_status[[i]] = "not_published"
      .plant_atomic_write_csv(manifest, manifest_path)
      if (service_busy) {
        paused = TRUE
        pause_reason = paste(
          "PubChem reported rate limiting or service unavailability.",
          "Validated completed batches and provider caches were preserved;",
          "resume later with the same paths."
        )
        break
      }
      if (isTRUE(stop_on_error)) stop(error_message, call. = FALSE)
      next
    }

    result$CategorateQueryMap = batch_map
    result$BatchMetadata = .categorate_enrichment_batch_metadata(
      run_signature, i, batch_queries, detail
    )
    .plant_atomic_save_rds(result, batch_files[[i]])
    batches[[i]] = result
    completed_this_run = completed_this_run + 1L
    manifest$status[[i]] = "completed"
    manifest$checkpoint_status[[i]] = "published"
    manifest$finished_at[[i]] = .plant_timestamp()
    manifest$completed_at[[i]] = manifest$finished_at[[i]]
    manifest$elapsed_seconds[[i]] = round(
      as.numeric(difftime(Sys.time(), started, units = "secs")), 3
    )
    manifest$resolved_count[[i]] = health$resolved_count
    manifest$property_count[[i]] = health$property_count
    manifest$error_message[[i]] = NA_character_
    .plant_atomic_write_csv(manifest, manifest_path)

    if (is.finite(cooldown_batches) && cooldown_seconds > 0 &&
        completed_this_run %% as.integer(cooldown_batches) == 0L &&
        i < length(groups)) {
      if (isTRUE(progress)) {
        message("uafR enrichment cooldown: ", cooldown_seconds, " seconds")
      }
      Sys.sleep(cooldown_seconds)
    }
  }

  for (i in seq_along(groups)) {
    if (is.null(batches[[i]]) && file.exists(batch_files[[i]])) {
      checkpoint = .categorate_enrichment_read_checkpoint(
        batch_files[[i]], run_signature, i, groups[[i]]
      )
      if (checkpoint$ok) batches[[i]] = checkpoint$result
    }
  }
  failed = manifest[manifest$status %in% c("failed", "paused_service_busy"),
                    , drop = FALSE]
  retry = manifest[manifest$status != "completed", , drop = FALSE]
  .plant_atomic_write_csv(failed,
                          file.path(out_dir, paste0(detail, "_failed_queries.csv")))
  .plant_atomic_write_csv(retry,
                          file.path(out_dir, paste0(detail, "_retry_queue.csv")))
  summary = .categorate_enrichment_summary(manifest, query_map, detail,
                                            paused, pause_reason)
  .plant_atomic_write_json(
    list(
      Format = "uafR_categorate_enrichment_run",
      SchemaVersion = "1.0.0",
      RunSignature = run_signature,
      Summary = summary,
      Paused = paused,
      PauseReason = pause_reason
    ),
    file.path(out_dir, paste0(detail, "_run_status.json")),
    dataframe = "rows", pretty = TRUE, auto_unbox = TRUE, na = "null"
  )
  completed_batches = batches[vapply(batches, is.list, logical(1))]
  names(completed_batches) = sprintf(
    "batch_%04d", which(vapply(batches, is.list, logical(1)))
  )
  attr(completed_batches, "batch_files") = batch_files[
    vapply(batches, is.list, logical(1))
  ]
  attr(completed_batches, "batch_indices") = which(
    vapply(batches, is.list, logical(1))
  )
  batch_summary = if (length(completed_batches) > 0) {
    summarizeCategorateBatches(completed_batches)
  } else .uaf_empty_table(.categorate_enrichment_batch_summary_cols())
  out = list(
    Batches = completed_batches,
    BatchSummary = batch_summary,
    BatchManifest = manifest,
    QueryMap = query_map,
    FailedQueries = failed,
    RetryQueue = retry,
    RunSummary = summary,
    OutputPath = out_dir
  )
  attr(out, "exit_status") = ifelse(paused, 75L,
                                    ifelse(nrow(retry) > 0, 1L, 0L))
  class(out) = c("uaf_categorate_enrichment_batches", "list")
  out
}

.categorate_enrichment_output_dir = function(path) {
  path = .uaf_first_non_empty_text(path)
  if (is.na(path)) stop("`out_dir` is required.", call. = FALSE)
  path = normalizePath(path, winslash = "/", mustWork = FALSE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Could not create `out_dir`: ", path,
                              call. = FALSE)
  path
}

.categorate_enrichment_positive_integer = function(x, name) {
  value = suppressWarnings(as.integer(x[[1]]))
  if (!is.finite(value) || value < 1) {
    stop("`", name, "` must be a positive integer.", call. = FALSE)
  }
  value
}

.categorate_enrichment_delay = function(x, name) {
  value = suppressWarnings(as.numeric(x[[1]]))
  if (!is.finite(value) || value < 0) {
    stop("`", name, "` must be a non-negative number.", call. = FALSE)
  }
  value
}

.categorate_enrichment_query_map = function(x, require_source_identity) {
  if (is.character(x)) {
    input = data.frame(compound_name = x, stringsAsFactors = FALSE)
  } else if (is.data.frame(x)) {
    input = x
  } else {
    stop("`compounds` must be a character vector or data frame.", call. = FALSE)
  }
  names(input) = .plant_normalize_column_names(names(input))
  name_col = intersect(c("compound_name", "resolved_name", "query", "name"),
                       names(input))
  if (length(name_col) < 1) name_col = names(input)[[1]] else
    name_col = name_col[[1]]
  n = nrow(input)
  values = function(candidates, default = NA_character_) {
    column = intersect(candidates, names(input))
    if (length(column) < 1) return(rep(default, n))
    as.character(input[[column[[1]]]])
  }
  compound_name = .uaf_squish_text(input[[name_col]])
  cid_raw = values(c("cid", "pubchem_cid"))
  cid = suppressWarnings(as.integer(sub("^cid[:= ]*", "", cid_raw,
                                         ignore.case = TRUE)))
  inchikey = toupper(.uaf_squish_text(values(c("inchikey", "inchi_key"))))
  valid_key = vapply(inchikey, .uaf_is_inchikey, logical(1))
  inchikey[!valid_key] = NA_character_
  smiles = .uaf_squish_text(values(c("smiles", "isomeric_smiles",
                                     "canonical_smiles")))
  compound_id = .uaf_squish_text(values(c("compound_id", "stable_compound_id")))
  fallback_id = ifelse(!is.na(inchikey), inchikey,
                       ifelse(!is.na(cid), paste0("cid_", cid),
                              .plant_slug(compound_name)))
  compound_id[is.na(compound_id) | compound_id == ""] =
    fallback_id[is.na(compound_id) | compound_id == ""]
  pubchem_query = ifelse(!is.na(cid), paste0("cid:", cid),
                         ifelse(!is.na(inchikey), inchikey, compound_name))
  query_type = ifelse(!is.na(cid), "source_cid",
                      ifelse(!is.na(inchikey), "source_inchikey", "exact_name"))
  source_identity = !is.na(cid) | !is.na(inchikey)
  eligible = !is.na(compound_name) & compound_name != "" &
    !is.na(pubchem_query) & pubchem_query != ""
  if (isTRUE(require_source_identity)) eligible = eligible & source_identity
  reason = ifelse(is.na(compound_name) | compound_name == "", "missing_name",
                  ifelse(isTRUE(require_source_identity) & !source_identity,
                         "source_identity_required", NA_character_))
  out = data.frame(
    compound_id = compound_id,
    compound_name = compound_name,
    compound_name_clean = .plant_clean_compound(compound_name),
    Query = compound_name,
    PubChemQuery = pubchem_query,
    QueryType = query_type,
    CID = cid,
    InChIKey = inchikey,
    SMILES = smiles,
    SourceIdentityAvailable = .uaf_yes_no(source_identity),
    EnrichmentEligible = .uaf_yes_no(eligible),
    ExclusionReason = reason,
    stringsAsFactors = FALSE
  )
  conflict_key = split(seq_len(nrow(out)), out$compound_name_clean)
  for (indices in conflict_key) {
    selectors = unique(.uaf_non_empty(out$PubChemQuery[indices]))
    if (length(selectors) > 1) {
      out$EnrichmentEligible[indices] = "No"
      out$ExclusionReason[indices] = "conflicting_identity_for_name"
    }
  }
  exact_key = paste(out$compound_name_clean, out$PubChemQuery, sep = "\r")
  out = out[!duplicated(exact_key), , drop = FALSE]
  row.names(out) = NULL
  out
}

.categorate_enrichment_run_signature = function(query_map, detail, batch_size,
                                                pubchem_throttle,
                                                kegg_throttle,
                                                assay_detail_limit,
                                                pubchem_annotation_mode,
                                                trait_profile, trait_mode,
                                                trait_confidence, max_traits) {
  version = tryCatch(as.character(utils::packageVersion("uafR")),
                     error = function(error) "development")
  payload = c(
    "uafR_categorate_enrichment_v2", version, detail, batch_size,
    pubchem_throttle, kegg_throttle, assay_detail_limit,
    pubchem_annotation_mode, trait_profile,
    trait_mode, trait_confidence, max_traits,
    apply(query_map[, c("compound_id", "Query", "PubChemQuery",
                        "EnrichmentEligible"), drop = FALSE], 1, paste,
          collapse = "|")
  )
  .pubchem_url_hash(paste(payload, collapse = "\r"))
}

.categorate_enrichment_manifest_cols = function() {
  c("run_signature", "batch_index", "query_start", "query_end",
    "query_count", "queries", "status", "attempt_count", "started_at",
    "finished_at", "completed_at", "elapsed_seconds", "resolved_count",
    "property_count", "checkpoint_status", "batch_file", "error_message")
}

.categorate_enrichment_manifest = function(groups, files, signature) {
  starts = cumsum(c(1L, vapply(groups, length, integer(1))))
  starts = starts[seq_along(groups)]
  ends = starts + vapply(groups, length, integer(1)) - 1L
  data.frame(
    run_signature = signature,
    batch_index = seq_along(groups),
    query_start = starts,
    query_end = ends,
    query_count = vapply(groups, length, integer(1)),
    queries = vapply(groups, paste, collapse = "; ", FUN.VALUE = character(1)),
    status = "pending",
    attempt_count = 0L,
    started_at = NA_character_,
    finished_at = NA_character_,
    completed_at = NA_character_,
    elapsed_seconds = NA_real_,
    resolved_count = 0L,
    property_count = 0L,
    checkpoint_status = "not_checked",
    batch_file = normalizePath(files, winslash = "/", mustWork = FALSE),
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

.categorate_enrichment_prior_manifest = function(path, signature, resume,
                                                 overwrite) {
  if (!isTRUE(resume) || !file.exists(path)) return(NULL)
  prior = tryCatch(utils::read.csv(path, stringsAsFactors = FALSE,
                                   check.names = FALSE),
                   error = function(error) NULL)
  if (!is.data.frame(prior)) return(NULL)
  signatures = unique(.uaf_non_empty(prior$run_signature))
  if (length(signatures) > 0 && !identical(signatures, signature) &&
      !isTRUE(overwrite)) {
    stop("Existing categorate enrichment manifest has a different run signature. Use a new output directory or `overwrite = TRUE`.",
         call. = FALSE)
  }
  if (!identical(signatures, signature)) return(NULL)
  prior
}

.categorate_enrichment_merge_manifest = function(current, prior) {
  if (!is.data.frame(prior) || nrow(prior) < 1) return(current)
  match_index = match(current$batch_index, prior$batch_index)
  hit = !is.na(match_index)
  for (column in intersect(c("attempt_count", "started_at", "finished_at",
                             "completed_at", "elapsed_seconds",
                             "resolved_count", "property_count",
                             "checkpoint_status", "error_message"),
                           names(prior))) {
    current[[column]][hit] = prior[[column]][match_index[hit]]
  }
  current
}

.categorate_enrichment_batch_metadata = function(signature, index, queries,
                                                  detail) {
  data.frame(
    Format = "uafR_categorate_enrichment_batch",
    SchemaVersion = "1.0.0",
    RunSignature = signature,
    BatchIndex = as.integer(index),
    Detail = detail,
    QueryCount = length(queries),
    QueryHash = .pubchem_url_hash(paste(queries, collapse = "\r")),
    PackageVersion = tryCatch(as.character(utils::packageVersion("uafR")),
                              error = function(error) "development"),
    CreatedAt = .plant_timestamp(),
    stringsAsFactors = FALSE
  )
}

.categorate_enrichment_read_checkpoint = function(path, signature, index,
                                                  queries) {
  if (!file.exists(path)) return(list(ok = FALSE, reason = "missing"))
  result = tryCatch(readRDS(path), error = function(error) NULL)
  if (!is.list(result) || is.data.frame(result)) {
    return(list(ok = FALSE, reason = "not_success_result"))
  }
  metadata = result$BatchMetadata
  if (!is.data.frame(metadata) || nrow(metadata) != 1) {
    return(list(ok = FALSE, reason = "metadata_missing"))
  }
  if (!identical(as.character(metadata$RunSignature[[1]]), signature)) {
    return(list(ok = FALSE, reason = "signature_mismatch"))
  }
  expected_hash = .pubchem_url_hash(paste(queries, collapse = "\r"))
  if (!identical(as.integer(metadata$BatchIndex[[1]]), as.integer(index)) ||
      !identical(as.character(metadata$QueryHash[[1]]), expected_hash)) {
    return(list(ok = FALSE, reason = "batch_mismatch"))
  }
  map = result$CategorateQueryMap
  health = .categorate_enrichment_batch_health(result, queries, map)
  if (!health$ok) return(list(ok = FALSE, reason = "incomplete"))
  list(ok = TRUE, reason = "valid", result = result)
}

.categorate_enrichment_batch_health = function(result, queries, query_map) {
  if (!is.list(result) || is.data.frame(result)) {
    return(list(ok = FALSE, message = "Batch result is not a categorate list."))
  }
  identity = result$PubChemIdentity
  properties = result$PubChemProperties
  if (!is.data.frame(identity) ||
      !all(c("Query", "CID", "MatchStatus") %in% names(identity))) {
    return(list(ok = FALSE, message = "PubChemIdentity is missing or malformed."))
  }
  missing_queries = setdiff(queries, .uaf_non_empty(identity$Query))
  if (length(missing_queries) > 0) {
    return(list(ok = FALSE, message = paste(
      "PubChem identity did not account for:", paste(missing_queries,
                                                      collapse = ", ")
    )))
  }
  expected_source = query_map$Query[
    query_map$SourceIdentityAvailable == "Yes" &
      query_map$EnrichmentEligible == "Yes"
  ]
  unresolved_source = expected_source[!expected_source %in%
    identity$Query[!is.na(identity$CID)]]
  if (length(unresolved_source) > 0) {
    return(list(ok = FALSE, message = paste(
      "Source-resolved compounds did not resolve in PubChem:",
      paste(unresolved_source, collapse = ", ")
    )))
  }
  resolved_cids = unique(identity$CID[!is.na(identity$CID)])
  property_cids = if (is.data.frame(properties) && "CID" %in% names(properties)) {
    unique(properties$CID[!is.na(properties$CID)])
  } else integer()
  missing_properties = setdiff(resolved_cids, property_cids)
  if (length(missing_properties) > 0) {
    return(list(ok = FALSE, message = paste(
      "PubChem properties are incomplete for CID(s):",
      paste(missing_properties, collapse = ", ")
    )))
  }
  validation = result$ValidationSummary
  if (is.data.frame(validation) && "ErrorCount" %in% names(validation) &&
      any(suppressWarnings(as.integer(validation$ErrorCount)) > 0,
          na.rm = TRUE)) {
    return(list(ok = FALSE,
                message = "Categorate validation reported one or more errors."))
  }
  list(ok = TRUE, message = NA_character_,
       resolved_count = length(resolved_cids),
       property_count = length(property_cids))
}

.categorate_enrichment_empty_standard_tables = function() {
  list(
    reactives = .uaf_empty_table(c("reactives", "Chemical")),
    LOTUS = .uaf_empty_table(c("LOTUS", "Chemical")),
    KEGG = .uaf_empty_table(c("KEGG", "Chemical")),
    FEMA = .uaf_empty_table(c("FEMA", "Chemical")),
    FDA_SPL = .uaf_empty_table(c("FDA_SPL", "Chemical")),
    FMCS = data.frame(),
    FunctionalGroups = data.frame(),
    BestChemMatch = data.frame()
  )
}

.categorate_enrichment_summary = function(manifest, query_map, detail, paused,
                                          pause_reason) {
  data.frame(
    Detail = detail,
    QueryCount = sum(query_map$EnrichmentEligible == "Yes"),
    ExcludedQueryCount = sum(query_map$EnrichmentEligible != "Yes"),
    BatchCount = nrow(manifest),
    CompletedBatchCount = sum(manifest$status == "completed"),
    FailedBatchCount = sum(manifest$status == "failed"),
    PendingBatchCount = sum(!manifest$status %in% c("completed", "failed")),
    Paused = .uaf_yes_no(paused),
    PauseReason = pause_reason,
    Completed = .uaf_yes_no(all(manifest$status == "completed")),
    UpdatedAt = .plant_timestamp(),
    stringsAsFactors = FALSE
  )
}

.categorate_enrichment_batch_summary_cols = function() {
  c("BatchIndex", "BatchName", "BatchFile", "Status", "ResolvedCIDCount",
    "PubChemPropertyRows", "PubChemPropertyRatio", "ValidationStatus",
    "QualityIssue", "Error")
}
