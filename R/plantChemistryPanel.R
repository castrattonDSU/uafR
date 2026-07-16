# Resumable production orchestration for large species-first plant panels.

.plant_panel_modes = function() {
  c("preflight", "build-indexes", "pilot", "discovery", "identity",
    "research-enrichment", "full-enrichment", "tanimoto-handoff",
    "finalize", "all")
}

#' Run a resumable plant chemistry production panel
#'
#' @description
#' Coordinates the production stages required to turn a curated plant list
#' into source-backed plant-compound evidence, conservative compound identity
#' resolution, research/full categorate enrichment, a Tanimoto server handoff,
#' and a validated analysis bundle. Each stage writes atomic lifecycle files
#' and can be resumed with the same input, cache, and output paths.
#'
#' Live provider work is never implied by `preflight` or `build-indexes`.
#' `mode = "all"` executes stages in order and stops with exit status 75 when a
#' provider circuit breaker requests a pause. Re-running the same command
#' reuses validated checkpoints and provider caches. Official NPASS resource
#' downloads are exact-byte checked and retain resumable `.partial` files;
#' local index construction uses lossless RDS staging and source/shard
#' checksums.
#'
#' @param mode One of `preflight`, `build-indexes`, `pilot`, `discovery`,
#' `identity`, `research-enrichment`, `full-enrichment`,
#' `tanimoto-handoff`, `finalize`, or `all`.
#' @param plant_csv Curated CSV containing one launch-safe plant query per row.
#' @param out_dir Persistent production output directory.
#' @param species_col Column containing canonical plant queries.
#' @param metadata_csv Optional CSV used as plant metadata. Defaults to
#' `plant_csv` so supplied taxonomy and project labels remain available.
#' @param exclusion_files Optional exclusion ledgers copied into the immutable
#' input snapshot. Counts are deduplicated across ledgers.
#' @param supporting_files Optional taxonomy, normalization, or other reviewed
#' ledgers snapshotted for provenance but not counted as exclusions.
#' @param cache_dir Persistent provider cache directory.
#' @param lotus_index Manifest-backed LOTUS lookup index or compatible local
#' index file/data frame.
#' @param npass_index Manifest-backed NPASS lookup index.
#' @param npass_raw_dir Directory containing the four official NPASS 3.0 text
#' files used by `buildNpassIndex()`.
#' @param download_npass Logical. Download missing NPASS files from the official
#' NPASS download host during `build-indexes`.
#' @param sources Plant discovery providers. Production defaults to all six
#' supported species-first providers.
#' @param taxon_fallback Optional fallback ranks run only after direct species
#' discovery. Fallback rows remain separate from direct evidence.
#' @param expected_species_count,expected_exclusion_count Optional contract
#' counts. Production runs should set these to their reviewed values.
#' @param require_all_providers Logical. If `TRUE`, unavailable or unqueried
#' configured providers fail the relevant quality gate.
#' @param species_chunk_size,compound_batch_size Resumable batch sizes.
#' @param pilot_count,pilot_species Deterministic pilot size or explicit pilot
#' species.
#' @param pilot_enrichment_compounds Maximum source-identified pilot compounds
#' submitted to the research enrichment gate. Pilot discovery is bounded to at
#' most `max(25, pilot_enrichment_compounds)` records per provider and species
#' (and never exceeds an explicitly smaller `max_provider_records` value), so a
#' pilot cannot silently become a full production identity run.
#' @param max_pubmed_records Maximum PubMed records fetched per species.
#' @param max_provider_records Maximum records per non-PubMed provider/species.
#' @param provider_throttle General live discovery delay in seconds.
#' @param knapsack_throttle KNApSAcK delay in seconds.
#' @param pubchem_throttle PubChem identity/enrichment delay in seconds.
#' @param kegg_throttle KEGG enrichment delay in seconds.
#' @param request_timeout Live request timeout in seconds.
#' @param full_enrichment_limit Maximum deterministic priority compounds sent
#' to `detail = "full"` enrichment.
#' @param release_manifest Optional checked uafR release-manifest JSON produced
#' by `tools/check_package_release.R --release-manifest`. Production runs should
#' supply this together with `source_tarball`.
#' @param source_tarball Optional checked uafR source tarball referenced by
#' `release_manifest`. The manifest and tarball are checksum-validated, included
#' in the run signature, snapshotted, and copied into the Tanimoto handoff.
#' @param require_release_artifact Logical. If `TRUE`, preflight fails unless a
#' clean-commit release manifest and matching source tarball are supplied.
#' @param server_results_dir Optional directory containing completed Tanimoto
#' server outputs used by `finalize`.
#' @param project_id Project label written to manifests and bundle methods.
#' @param resume Logical. Reuse validated stage markers and checkpoints.
#' @param overwrite Logical. Rebuild selected stage outputs while preserving
#' provider caches.
#' @param progress Logical. Print concise stage progress.
#' @param provider_results Optional provider fixture results for offline tests.
#' @param request_fun,pubtator_request_fun,compound_request_fun Optional request
#' functions used for deterministic tests.
#' @param enrichment_fun Optional complete enrichment function used in tests.
#'
#' @return A list with `Status`, `Progress`, `StageResults`, `OutputPath`, and
#' `ExitStatus`. Exit status 75 means a service-busy pause; 2 means finalization
#' is awaiting externally generated server results.
#'
#' @export
runPlantChemistryPanel = function(
    mode = .plant_panel_modes(),
    plant_csv,
    out_dir,
    species_col = "species",
    metadata_csv = NULL,
    exclusion_files = NULL,
    supporting_files = NULL,
    cache_dir = NULL,
    lotus_index = Sys.getenv("UAFR_LOTUS_INDEX", ""),
    npass_index = Sys.getenv("UAFR_NPASS_INDEX", ""),
    npass_raw_dir = NULL,
    download_npass = FALSE,
    sources = c("lotus", "npass", "knapsack", "pubchem", "pubmed",
                "pubtator"),
    taxon_fallback = c("genus", "family"),
    expected_species_count = NULL,
    expected_exclusion_count = NULL,
    require_all_providers = TRUE,
    species_chunk_size = 25,
    compound_batch_size = 25,
    pilot_count = 12,
    pilot_species = NULL,
    pilot_enrichment_compounds = 25,
    max_pubmed_records = 100,
    max_provider_records = Inf,
    provider_throttle = 0.5,
    knapsack_throttle = 1,
    pubchem_throttle = 1.1,
    kegg_throttle = 0.5,
    request_timeout = 60,
    full_enrichment_limit = 250,
    release_manifest = NULL,
    source_tarball = NULL,
    require_release_artifact = FALSE,
    server_results_dir = NULL,
    project_id = NULL,
    resume = TRUE,
    overwrite = FALSE,
    progress = interactive(),
    provider_results = NULL,
    request_fun = NULL,
    pubtator_request_fun = NULL,
    compound_request_fun = NULL,
    enrichment_fun = NULL) {
  mode = match.arg(mode)
  context = .plant_panel_context(
    mode = mode, plant_csv = plant_csv, out_dir = out_dir,
    species_col = species_col, metadata_csv = metadata_csv,
    exclusion_files = exclusion_files, supporting_files = supporting_files,
    cache_dir = cache_dir,
    lotus_index = lotus_index, npass_index = npass_index,
    npass_raw_dir = npass_raw_dir, download_npass = download_npass,
    sources = sources, taxon_fallback = taxon_fallback,
    expected_species_count = expected_species_count,
    expected_exclusion_count = expected_exclusion_count,
    require_all_providers = require_all_providers,
    species_chunk_size = species_chunk_size,
    compound_batch_size = compound_batch_size, pilot_count = pilot_count,
    pilot_species = pilot_species,
    pilot_enrichment_compounds = pilot_enrichment_compounds,
    max_pubmed_records = max_pubmed_records,
    max_provider_records = max_provider_records,
    provider_throttle = provider_throttle,
    knapsack_throttle = knapsack_throttle,
    pubchem_throttle = pubchem_throttle, kegg_throttle = kegg_throttle,
    request_timeout = request_timeout,
    full_enrichment_limit = full_enrichment_limit,
    release_manifest = release_manifest, source_tarball = source_tarball,
    require_release_artifact = require_release_artifact,
    server_results_dir = server_results_dir, project_id = project_id,
    resume = resume, overwrite = overwrite, progress = progress,
    provider_results = provider_results, request_fun = request_fun,
    pubtator_request_fun = pubtator_request_fun,
    compound_request_fun = compound_request_fun,
    enrichment_fun = enrichment_fun
  )
  stages = if (identical(mode, "all")) {
    setdiff(.plant_panel_modes(), "all")
  } else {
    mode
  }
  stage_results = list()
  exit_status = 0L
  .plant_panel_write_status(context, "running", stage = stages[[1]],
                            message = "Plant chemistry panel started.")

  for (stage in stages) {
    reused = isTRUE(context$config$resume) &&
      !isTRUE(context$config$overwrite) &&
      .plant_panel_stage_marker_valid(context, stage)
    if (reused) {
      .plant_panel_clear_transient_markers(context)
      .plant_panel_progress(context, stage, "reused", Sys.time(), Sys.time(),
                            "Validated completed stage reused.")
      stage_results[[stage]] = list(status = "reused", exit_status = 0L)
      next
    }
    started = Sys.time()
    .plant_panel_progress(context, stage, "running", started, NA,
                          "Stage started.")
    .plant_panel_write_status(context, "running", stage = stage,
                              message = "Stage started.")
    result = tryCatch(
      .plant_panel_run_stage(stage, context),
      interrupt = function(condition) {
        .plant_panel_stage_failure(context, stage, condition, "interrupted")
        stop(condition)
      },
      error = function(condition) {
        .plant_panel_stage_failure(context, stage, condition, "failed")
        stop(condition)
      }
    )
    if (is.list(result$config_updates)) {
      for (name in names(result$config_updates)) {
        context$config[[name]] = result$config_updates[[name]]
      }
    }
    status = .uaf_first_non_empty_text(result$status, "completed")
    exit_status = suppressWarnings(as.integer(
      .uaf_first_non_empty_text(result$exit_status, 0L)
    ))
    if (identical(status, "completed")) {
      .plant_panel_write_stage_marker(
        context, stage, result$artifacts, started, result$message
      )
      .plant_panel_clear_transient_markers(context)
    }
    .plant_panel_progress(context, stage, status, started, Sys.time(),
                          result$message)
    stage_results[[stage]] = result
    if (exit_status != 0L || !identical(status, "completed")) {
      state = if (exit_status == 75L) "paused_service_busy" else
        if (exit_status == 2L) "awaiting_server" else status
      .plant_panel_write_status(context, state, stage = stage,
                                message = result$message)
      if (exit_status == 75L) {
        .plant_atomic_write_json(
          list(stage = stage, state = state, message = result$message,
               paused_at = .plant_timestamp()),
          file.path(context$paths$root, "PAUSED_SERVICE_BUSY.json")
        )
      }
      break
    }
  }

  if (exit_status == 0L && all(vapply(stage_results, function(x) {
    .uaf_first_non_empty_text(x$status) %in% c("completed", "reused")
  }, logical(1)))) {
    .plant_panel_write_status(context, "completed",
                              stage = utils::tail(stages, 1),
                              message = "Requested panel stages completed.")
    if (identical(mode, "all")) {
      .plant_atomic_write_text(
        c("workflow: runPlantChemistryPanel", "status: completed",
          paste0("run_signature: ", context$run_signature),
          paste0("completed_at: ", .plant_timestamp())),
        file.path(context$paths$root, "PIPELINE_COMPLETED.txt")
      )
    }
  }
  status = .plant_panel_read_json(context$paths$status)
  progress_table = .plant_panel_read_csv(context$paths$progress)
  out = list(Status = status, Progress = progress_table,
             StageResults = stage_results, OutputPath = context$paths$root,
             ExitStatus = exit_status)
  class(out) = c("uaf_plant_chemistry_panel", "list")
  out
}

.plant_panel_context = function(mode, plant_csv, out_dir, species_col,
                                 metadata_csv, exclusion_files,
                                 supporting_files, cache_dir,
                                 lotus_index, npass_index, npass_raw_dir,
                                 download_npass, sources, taxon_fallback,
                                 expected_species_count,
                                 expected_exclusion_count,
                                 require_all_providers, species_chunk_size,
                                 compound_batch_size, pilot_count,
                                 pilot_species, pilot_enrichment_compounds,
                                 max_pubmed_records, max_provider_records,
                                 provider_throttle, knapsack_throttle,
                                 pubchem_throttle, kegg_throttle,
                                 request_timeout, full_enrichment_limit,
                                 release_manifest, source_tarball,
                                 require_release_artifact,
                                 server_results_dir, project_id, resume,
                                 overwrite, progress, provider_results,
                                 request_fun, pubtator_request_fun,
                                 compound_request_fun, enrichment_fun) {
  plant_csv = .plant_panel_existing_file(plant_csv, "plant_csv")
  root = .plant_panel_path(out_dir, "out_dir")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  cache_dir = .uaf_first_non_empty_text(cache_dir,
                                        file.path(root, "provider_cache"))
  cache_dir = normalizePath(cache_dir, winslash = "/", mustWork = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  paths = list(
    root = root,
    inputs = file.path(root, "inputs"),
    resources = file.path(root, "resources"),
    stages = file.path(root, "stages"),
    markers = file.path(root, "stage_markers"),
    handoff = file.path(root, "tanimoto_handoff"),
    bundle = file.path(root, "plant_chemistry_analysis_bundle"),
    status = file.path(root, "pipeline_status.json"),
    progress = file.path(root, "pipeline_progress.csv")
  )
  invisible(lapply(paths[c("inputs", "resources", "stages", "markers")],
                   dir.create, recursive = TRUE, showWarnings = FALSE))
  sources = .plant_normalize_sources(sources)
  if (length(sources) < 1) stop("At least one provider is required.",
                                call. = FALSE)
  existing_indexes = .plant_panel_index_config(paths$resources)
  lotus_index = .uaf_first_non_empty_text(lotus_index,
                                           existing_indexes$lotus)
  npass_index = .uaf_first_non_empty_text(npass_index,
                                           existing_indexes$npass)
  package_version = tryCatch(as.character(utils::packageVersion("uafR")),
                             error = function(error) "development")
  config = list(
    mode = mode, plant_csv = plant_csv, species_col = species_col,
    metadata_csv = .uaf_first_non_empty_text(metadata_csv, plant_csv),
    exclusion_files = .plant_panel_existing_files(exclusion_files),
    supporting_files = .plant_panel_existing_files(supporting_files),
    cache_dir = cache_dir, lotus_index = lotus_index,
    npass_index = npass_index, npass_raw_dir = npass_raw_dir,
    download_npass = isTRUE(download_npass), sources = sources,
    taxon_fallback = unique(tolower(.uaf_non_empty(taxon_fallback))),
    expected_species_count = .plant_panel_optional_integer(
      expected_species_count, "expected_species_count"
    ),
    expected_exclusion_count = .plant_panel_optional_integer(
      expected_exclusion_count, "expected_exclusion_count"
    ),
    require_all_providers = isTRUE(require_all_providers),
    species_chunk_size = .plant_panel_positive_integer(
      species_chunk_size, "species_chunk_size"
    ),
    compound_batch_size = .plant_panel_positive_integer(
      compound_batch_size, "compound_batch_size"
    ),
    pilot_count = .plant_panel_positive_integer(pilot_count, "pilot_count"),
    pilot_species = unique(.uaf_non_empty(pilot_species)),
    pilot_enrichment_compounds = .plant_panel_positive_integer(
      pilot_enrichment_compounds, "pilot_enrichment_compounds"
    ),
    max_pubmed_records = .plant_panel_positive_integer(
      max_pubmed_records, "max_pubmed_records"
    ),
    max_provider_records = .plant_panel_positive_or_inf(
      max_provider_records, "max_provider_records"
    ),
    provider_throttle = .plant_panel_nonnegative(provider_throttle,
                                                  "provider_throttle"),
    knapsack_throttle = .plant_panel_nonnegative(knapsack_throttle,
                                                  "knapsack_throttle"),
    pubchem_throttle = .plant_panel_nonnegative(pubchem_throttle,
                                                 "pubchem_throttle"),
    kegg_throttle = .plant_panel_nonnegative(kegg_throttle,
                                              "kegg_throttle"),
    request_timeout = .plant_panel_nonnegative(request_timeout,
                                                "request_timeout"),
    full_enrichment_limit = .plant_panel_positive_integer(
      full_enrichment_limit, "full_enrichment_limit"
    ),
    release_manifest = .plant_panel_optional_existing_file(
      release_manifest, "release_manifest"
    ),
    source_tarball = .plant_panel_optional_existing_file(
      source_tarball, "source_tarball"
    ),
    require_release_artifact = isTRUE(require_release_artifact),
    server_results_dir = .uaf_first_non_empty_text(server_results_dir),
    project_id = .uaf_first_non_empty_text(project_id,
                                            "plant_chemistry_panel"),
    resume = isTRUE(resume), overwrite = isTRUE(overwrite),
    progress = isTRUE(progress), package_version = package_version,
    provider_results = provider_results, request_fun = request_fun,
    pubtator_request_fun = pubtator_request_fun,
    compound_request_fun = compound_request_fun,
    enrichment_fun = enrichment_fun
  )
  signature_values = c(
    "uafR_plant_chemistry_panel_v3", package_version,
    unname(tools::md5sum(plant_csv)[[1]]), species_col,
    .plant_panel_input_signature(c(
      config$metadata_csv, config$exclusion_files, config$supporting_files
    )),
    paste(sources, collapse = ","),
    paste(config$taxon_fallback, collapse = ","),
    config$require_all_providers,
    config$species_chunk_size, config$compound_batch_size,
    config$pilot_count, paste(config$pilot_species, collapse = ","),
    config$pilot_enrichment_compounds,
    config$max_pubmed_records, config$max_provider_records,
    config$provider_throttle, config$knapsack_throttle,
    config$pubchem_throttle, config$kegg_throttle, config$request_timeout,
    config$full_enrichment_limit,
    .plant_panel_input_signature(c(config$release_manifest,
                                   config$source_tarball)),
    config$require_release_artifact,
    .pubchem_url_hash(paste(
      Sys.getenv("NCBI_EMAIL", ""), Sys.getenv("NCBI_TOOL", "uafR"),
      sep = "\r"
    )),
    .plant_batch_provider_results_signature(provider_results)
  )
  list(config = config, paths = paths,
       run_signature = .pubchem_url_hash(paste(signature_values,
                                                collapse = "\r")))
}

.plant_panel_input_signature = function(paths) {
  paths = unique(.uaf_non_empty(paths))
  paths = paths[file.exists(paths) & !dir.exists(paths)]
  if (length(paths) < 1L) return("supporting_inputs:none")
  paths = normalizePath(paths, winslash = "/", mustWork = TRUE)
  paste(paste(paths, unname(tools::md5sum(paths)), sep = "="),
        collapse = "|")
}

.plant_panel_release_check = function(check, status, observed, expected) {
  data.frame(
    check = check, status = status, observed = as.character(observed),
    expected = as.character(expected), stringsAsFactors = FALSE
  )
}

.plant_panel_validate_release_artifacts = function(config) {
  manifest_path = config$release_manifest
  tarball_path = config$source_tarball
  manifest_present = !is.na(manifest_path) && file.exists(manifest_path)
  tarball_present = !is.na(tarball_path) && file.exists(tarball_path)
  required = isTRUE(config$require_release_artifact)
  if (!manifest_present && !tarball_present) {
    return(list(
      Checks = .plant_panel_release_check(
        "release_artifact_contract", ifelse(required, "fail", "pass"),
        ifelse(required, "not supplied", "optional; not supplied"),
        ifelse(required, "checked manifest and matching source tarball",
               "optional for non-production runs")
      ),
      Manifest = NULL, GitCommit = NA_character_
    ))
  }
  presence = rbind(
    .plant_panel_release_check(
      "release_manifest_present", ifelse(manifest_present, "pass", "fail"),
      ifelse(manifest_present, basename(manifest_path), "missing"),
      "existing checked release-manifest JSON"
    ),
    .plant_panel_release_check(
      "source_tarball_present", ifelse(tarball_present, "pass", "fail"),
      ifelse(tarball_present, basename(tarball_path), "missing"),
      "existing checked uafR source tarball"
    )
  )
  if (!manifest_present || !tarball_present) {
    return(list(Checks = presence, Manifest = NULL,
                GitCommit = NA_character_))
  }
  manifest = tryCatch(
    jsonlite::read_json(manifest_path, simplifyVector = TRUE),
    error = function(error) error
  )
  if (inherits(manifest, "error") || !is.list(manifest)) {
    message = if (inherits(manifest, "error")) {
      .plant_redact_secrets(conditionMessage(manifest))
    } else "JSON root is not an object"
    return(list(
      Checks = rbind(
        presence,
        .plant_panel_release_check("release_manifest_readable", "fail",
                                   message, "valid JSON object")
      ),
      Manifest = NULL, GitCommit = NA_character_
    ))
  }
  value = function(name) .uaf_first_non_empty_text(manifest[[name]])
  package_version = value("package_version")
  git_commit = value("git_commit")
  dirty_state = tolower(value("dirty_state_check"))
  expected_name = value("source_tarball")
  expected_bytes = suppressWarnings(as.numeric(manifest$source_tarball_bytes))
  expected_md5 = tolower(value("source_tarball_md5"))
  expected_sha256 = tolower(value("source_tarball_sha256"))
  observed_bytes = as.numeric(file.info(tarball_path)$size)
  observed_md5 = tolower(unname(tools::md5sum(tarball_path)[[1L]]))
  observed_sha256 = tolower(.plant_sha256_file(tarball_path))
  checks = rbind(
    presence,
    .plant_panel_release_check("release_manifest_readable", "pass",
                               basename(manifest_path), "valid JSON object"),
    .plant_panel_release_check(
      "release_manifest_package_version",
      ifelse(identical(package_version, config$package_version),
             "pass", "fail"),
      package_version, config$package_version
    ),
    .plant_panel_release_check(
      "release_manifest_clean_git_commit",
      ifelse(identical(dirty_state, "clean") &&
               !is.na(git_commit) &&
               grepl("^[0-9a-fA-F]{40}$", git_commit), "pass", "fail"),
      paste("dirty_state", dirty_state, "commit", git_commit),
      "dirty_state clean and a 40-character Git commit"
    ),
    .plant_panel_release_check(
      "release_manifest_tarball_name",
      ifelse(!is.na(expected_name) &&
               identical(basename(tarball_path), basename(expected_name)),
             "pass", "fail"),
      basename(tarball_path), expected_name
    ),
    .plant_panel_release_check(
      "release_manifest_tarball_bytes",
      ifelse(is.finite(expected_bytes) &&
               identical(observed_bytes, expected_bytes), "pass", "fail"),
      observed_bytes, expected_bytes
    ),
    .plant_panel_release_check(
      "release_manifest_tarball_md5",
      ifelse(!is.na(expected_md5) && identical(observed_md5, expected_md5),
             "pass", "fail"),
      observed_md5, expected_md5
    ),
    .plant_panel_release_check(
      "release_manifest_tarball_sha256",
      ifelse(!is.na(expected_sha256) &&
               identical(observed_sha256, expected_sha256), "pass", "fail"),
      observed_sha256, expected_sha256
    )
  )
  list(Checks = checks, Manifest = manifest, GitCommit = git_commit)
}

.plant_panel_run_stage = function(stage, context) {
  switch(
    stage,
    preflight = .plant_panel_stage_preflight(context),
    `build-indexes` = .plant_panel_stage_build_indexes(context),
    pilot = .plant_panel_stage_pilot(context),
    discovery = .plant_panel_stage_discovery(context),
    identity = .plant_panel_stage_identity(context),
    `research-enrichment` = .plant_panel_stage_research(context),
    `full-enrichment` = .plant_panel_stage_full(context),
    `tanimoto-handoff` = .plant_panel_stage_tanimoto_handoff(context),
    finalize = .plant_panel_stage_finalize(context),
    stop("Unsupported plant panel stage: ", stage, call. = FALSE)
  )
}

.plant_panel_stage_preflight = function(context) {
  config = context$config
  input = .plant_panel_read_input(config$plant_csv, config$species_col)
  species = input$species
  queries = .plant_queries(input, taxon_fallback = "species")
  aliases = .plant_query_aliases(input, queries)
  missing_authority = is.na(.uaf_squish_text(aliases$authority))
  aliases$authority[missing_authority &
                      aliases$alias_type == "submitted_name"] =
    paste0("curated_input:", basename(config$plant_csv))
  aliases$authority[missing_authority & grepl(
    "wfo", aliases$source_field, ignore.case = TRUE
  )] = paste0("WFO:", .uaf_first_non_empty_text(input$wfo_classification_version,
                                                  "version_not_supplied"))
  aliases$authority[
    is.na(.uaf_squish_text(aliases$authority)) & aliases$verified == "No"
  ] = "unverified_input_not_query_authorized"

  release_validation = .plant_panel_validate_release_artifacts(config)

  snapshot_spec = data.frame(
    path = c(config$plant_csv, config$metadata_csv,
             config$exclusion_files, config$supporting_files,
             config$release_manifest, config$source_tarball),
    role = c("plant_input", "plant_metadata",
             rep("exclusion_ledger", length(config$exclusion_files)),
             rep("supporting_ledger", length(config$supporting_files)),
             "release_manifest", "source_tarball"),
    stringsAsFactors = FALSE
  )
  snapshot_spec = snapshot_spec[
    !is.na(snapshot_spec$path) & nzchar(snapshot_spec$path) &
      file.exists(snapshot_spec$path), , drop = FALSE
  ]
  snapshot_spec = snapshot_spec[
    !duplicated(normalizePath(snapshot_spec$path, winslash = "/",
                              mustWork = TRUE)), , drop = FALSE
  ]
  snapshots = .plant_panel_snapshot_files(snapshot_spec$path,
                                           context$paths$inputs,
                                           context$run_signature,
                                           config$overwrite,
                                           roles = snapshot_spec$role)
  query_file = file.path(context$paths$inputs, "PlantQueries.csv")
  alias_file = file.path(context$paths$inputs, "PlantQueryAliases.csv")
  .plant_atomic_write_csv(queries, query_file)
  .plant_atomic_write_csv(aliases, alias_file)
  exclusion_rows = .plant_panel_exclusion_count(config$exclusion_files)
  launch_safe = if ("taxonomy_name_launch_safe" %in% names(input)) {
    .plant_panel_truthy(input$taxonomy_name_launch_safe)
  } else {
    rep(TRUE, nrow(input))
  }
  checks = data.frame(
    check = c("input_readable", "species_nonblank", "species_unique",
              "species_count_contract", "exclusion_count_contract",
              "launch_safe_queries", "alias_authority_present",
              "ncbi_email_configured"),
    status = c(
      "pass",
      ifelse(all(!is.na(species) & species != ""), "pass", "fail"),
      ifelse(anyDuplicated(.plant_clean_name(species)) == 0L,
             "pass", "fail"),
      ifelse(is.null(config$expected_species_count) ||
               nrow(input) == config$expected_species_count,
             "pass", "fail"),
      ifelse(is.null(config$expected_exclusion_count) ||
               exclusion_rows == config$expected_exclusion_count,
             "pass", "fail"),
      ifelse(all(launch_safe), "pass", "fail"),
      ifelse(all(!is.na(.uaf_squish_text(aliases$authority))),
             "pass", "fail"),
      ifelse(length(.uaf_non_empty(Sys.getenv("NCBI_EMAIL", ""))) > 0,
             "pass", "warning")
    ),
    observed = c(
      config$plant_csv,
      sum(is.na(species) | species == ""),
      anyDuplicated(.plant_clean_name(species)),
      nrow(input), exclusion_rows, sum(launch_safe),
      sum(is.na(.uaf_squish_text(aliases$authority))),
      ifelse(nzchar(Sys.getenv("NCBI_EMAIL", "")), "configured",
             "not_configured")
    ),
    expected = c("readable CSV", "0 blank", "0 duplicates",
                 .uaf_first_non_empty_text(config$expected_species_count,
                                           "not_fixed"),
                 .uaf_first_non_empty_text(config$expected_exclusion_count,
                                           "not_fixed"),
                 nrow(input), "0 missing", "configured before live NCBI work"),
    stringsAsFactors = FALSE
  )
  checks = rbind(checks, release_validation$Checks)
  checks_file = file.path(context$paths$stages, "preflight",
                          "preflight_checks.csv")
  dir.create(dirname(checks_file), recursive = TRUE, showWarnings = FALSE)
  .plant_atomic_write_csv(checks, checks_file)
  manifest_file = file.path(context$paths$inputs, "input_manifest.csv")
  .plant_atomic_write_csv(snapshots, manifest_file)
  plan = planPlantChemistryRun(
    plants = species, sources = config$sources,
    cache_dir = config$cache_dir, lotus_index = config$lotus_index,
    species_chunk_size = config$species_chunk_size,
    compound_batch_size = config$compound_batch_size,
    max_pubmed_records = config$max_pubmed_records
  )
  plan_dir = file.path(context$paths$stages, "preflight", "run_plan")
  dir.create(plan_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(plan)) {
    if (is.data.frame(plan[[name]])) {
      .plant_atomic_write_csv(
        plan[[name]], file.path(plan_dir, paste0(name, ".csv"))
      )
    }
  }
  contract = list(
    workflow = "runPlantChemistryPanel",
    schema_version = "1.0.0", run_signature = context$run_signature,
    package_version = config$package_version,
    project_id = config$project_id,
    created_at = .plant_timestamp(),
    plant_csv = config$plant_csv, plant_csv_md5 = snapshots$md5[
      snapshots$source_path == config$plant_csv
    ],
    species_count = nrow(input), exclusion_count = exclusion_rows,
    sources = config$sources, taxon_fallback = config$taxon_fallback,
    release_artifact_required = config$require_release_artifact,
    release_manifest = if (!is.na(config$release_manifest)) {
      basename(config$release_manifest)
    } else NULL,
    release_manifest_md5 = if (!is.na(config$release_manifest)) {
      unname(tools::md5sum(config$release_manifest)[[1L]])
    } else NULL,
    source_tarball = if (!is.na(config$source_tarball)) {
      basename(config$source_tarball)
    } else NULL,
    source_tarball_md5 = if (!is.na(config$source_tarball)) {
      unname(tools::md5sum(config$source_tarball)[[1L]])
    } else NULL,
    source_tarball_sha256 = if (!is.na(config$source_tarball)) {
      .plant_sha256_file(config$source_tarball)
    } else NULL,
    git_commit = release_validation$GitCommit,
    scientific_limit = paste(
      "Records are reported public-source evidence, not a complete metabolome",
      "and not confirmation in a project sample."
    )
  )
  contract_file = file.path(context$paths$inputs, "run_contract.json")
  .plant_atomic_write_json(contract, contract_file)
  if (any(checks$status == "fail")) {
    stop("Panel preflight failed: ",
         paste(checks$check[checks$status == "fail"], collapse = ", "),
         call. = FALSE)
  }
  list(status = "completed", exit_status = 0L,
       message = paste(nrow(input), "plant queries passed preflight."),
       artifacts = c(checks_file, manifest_file, query_file, alias_file,
                     contract_file,
                     file.path(context$paths$inputs,
                               snapshots$snapshot_file)))
}

.plant_panel_stage_build_indexes = function(context) {
  config = context$config
  lotus = .plant_lotus_index_or_null(config$lotus_index)
  if ("lotus" %in% config$sources && is.null(lotus)) {
    stop("A valid local LOTUS index is required for production.",
         call. = FALSE)
  }
  npass = .plant_npass_index_or_null(config$npass_index)
  if ("npass" %in% config$sources && is.null(.plant_npass_lookup_info(npass))) {
    raw = .plant_panel_npass_raw_paths(config$npass_raw_dir)
    if (isTRUE(config$download_npass)) {
      raw = .plant_panel_download_npass(config$npass_raw_dir,
                                        progress = config$progress,
                                        request_timeout = max(
                                          3600, config$request_timeout
                                        ))
    }
    missing = names(raw)[!file.exists(raw)]
    if (length(missing) > 0) {
      stop("NPASS index is unavailable and source files are missing: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    npass = file.path(context$paths$resources, "NPASS_3.0_2026_index")
    buildNpassIndex(
      general_info = raw[["general_info"]],
      structures = raw[["structures"]],
      species_pairs = raw[["species_pairs"]],
      species_info = raw[["species_info"]],
      out_dir = npass, overwrite = config$overwrite,
      progress = config$progress
    )
  }
  indexes = list(lotus = .uaf_first_non_empty_text(config$lotus_index),
                 npass = .uaf_first_non_empty_text(npass))
  index_file = file.path(context$paths$resources, "provider_indexes.json")
  .plant_atomic_write_json(indexes, index_file)
  availability = plantProviderAvailability(
    sources = config$sources,
    provider_indexes = list(lotus = indexes$lotus, npass = indexes$npass),
    lotus_index = indexes$lotus, probe_live = FALSE
  )
  availability_file = file.path(context$paths$resources,
                                "ProviderResourceManifest.csv")
  .plant_atomic_write_csv(availability, availability_file)
  local_required = availability$provider %in% c("lotus", "npass")
  if (any(local_required & availability$availability_status != "available")) {
    stop("One or more required local provider indexes failed validation.",
         call. = FALSE)
  }
  list(
    status = "completed", exit_status = 0L,
    message = "Local provider resource indexes are manifest-backed and readable.",
    artifacts = c(index_file, availability_file,
                  .plant_panel_index_manifest(indexes$lotus),
                  .plant_panel_index_manifest(indexes$npass)),
    config_updates = list(lotus_index = indexes$lotus,
                          npass_index = indexes$npass)
  )
}

.plant_panel_stage_pilot = function(context) {
  .plant_panel_require_resource_stage(context)
  config = context$config
  input = .plant_panel_read_input(config$plant_csv, config$species_col)
  pilot = .plant_panel_select_pilot(input, config$pilot_species,
                                    config$pilot_count)
  stage_dir = file.path(context$paths$stages, "pilot")
  discovery_dir = file.path(stage_dir, "discovery")
  indexes = .plant_panel_provider_indexes(context)
  pilot_file = file.path(stage_dir, "pilot_species.csv")
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  .plant_atomic_write_csv(pilot, pilot_file)
  pilot_record_limit = .plant_panel_pilot_record_limit(
    config$max_provider_records, config$pilot_enrichment_compounds
  )
  pilot_pubmed_limit = min(config$max_pubmed_records, pilot_record_limit)
  .plant_panel_clear_derived_exports(discovery_dir)

  result = runPlantPhytochemistryBatch(
    plants = pilot,
    sources = config$sources,
    taxon_fallback = "species",
    out_dir = discovery_dir,
    cache_dir = config$cache_dir,
    species_chunk_size = min(config$species_chunk_size, nrow(pilot)),
    compound_resolution_profile = "none",
    occurrence_status = c("direct_reported", "curated_reported"),
    analysis_ready = TRUE, min_confidence = "medium",
    cache = TRUE,
    throttle = max(config$provider_throttle, config$knapsack_throttle),
    ncbi_email = Sys.getenv("NCBI_EMAIL", ""),
    ncbi_tool = Sys.getenv("NCBI_TOOL", "uafR"),
    ncbi_api_key = Sys.getenv("NCBI_API_KEY", ""),
    max_pubmed_records = pilot_pubmed_limit,
    max_provider_records = pilot_record_limit,
    lotus_index = indexes$lotus, provider_indexes = indexes,
    require_all_providers = config$require_all_providers,
    request_timeout = config$request_timeout,
    provider_results = config$provider_results,
    request_fun = config$request_fun,
    pubtator_request_fun = config$pubtator_request_fun,
    compound_request_fun = config$compound_request_fun,
    compound_batch_size = config$compound_batch_size,
    resume = config$resume, progress = config$progress,
    overwrite = config$overwrite, stop_on_error = FALSE,
    allow_large_live_run = TRUE
  )
  result_file = file.path(stage_dir, "pilot_result.rds")
  .plant_atomic_save_rds(result, result_file)
  cache_fingerprint = .plant_panel_result_fingerprint(result)
  # The replay verifies checkpoint reuse. Remove only derived top-level exports
  # so `overwrite = FALSE` can protect checkpoints while regenerating the same
  # auditable tables from them.
  .plant_panel_clear_derived_exports(discovery_dir)
  rerun = runPlantPhytochemistryBatch(
    plants = pilot,
    sources = config$sources,
    taxon_fallback = "species",
    out_dir = discovery_dir,
    cache_dir = config$cache_dir,
    species_chunk_size = min(config$species_chunk_size, nrow(pilot)),
    compound_resolution_profile = "none",
    cache = TRUE,
    throttle = max(config$provider_throttle, config$knapsack_throttle),
    ncbi_email = Sys.getenv("NCBI_EMAIL", ""),
    ncbi_tool = Sys.getenv("NCBI_TOOL", "uafR"),
    ncbi_api_key = Sys.getenv("NCBI_API_KEY", ""),
    max_pubmed_records = pilot_pubmed_limit,
    max_provider_records = pilot_record_limit,
    lotus_index = indexes$lotus, provider_indexes = indexes,
    require_all_providers = config$require_all_providers,
    request_timeout = config$request_timeout,
    provider_results = config$provider_results,
    request_fun = .plant_panel_cache_only_request,
    pubtator_request_fun = .plant_panel_cache_only_request,
    compound_request_fun = .plant_panel_cache_only_request,
    compound_batch_size = config$compound_batch_size,
    resume = TRUE, progress = FALSE, overwrite = FALSE,
    stop_on_error = FALSE, allow_large_live_run = TRUE
  )
  rerun_fingerprint = .plant_panel_result_fingerprint(rerun)
  provider_gate = .plant_panel_provider_gate(result, config$sources,
                                               nrow(pilot))
  cache_gate = data.frame(
    check = "cache_only_rerun_identical",
    status = ifelse(identical(cache_fingerprint, rerun_fingerprint),
                    "pass", "fail"),
    observed = rerun_fingerprint, expected = cache_fingerprint,
    stringsAsFactors = FALSE
  )
  gate = rbind(provider_gate, cache_gate)

  enrichment_artifacts = character()
  eligible = .plant_panel_pilot_enrichment_input(
    result$SourceCompoundIdentity, config$pilot_enrichment_compounds
  )
  selection_file = file.path(stage_dir, "pilot_enrichment_selection.csv")
  .plant_atomic_write_csv(eligible, selection_file)
  enrichment_artifacts = c(enrichment_artifacts, selection_file)
  gate = rbind(gate, data.frame(
    check = c("pilot_discovery_record_limit",
              "pilot_enrichment_source_identity"),
    status = c(
      ifelse(is.finite(pilot_record_limit) && pilot_record_limit >= 1L,
             "pass", "fail"),
      ifelse(nrow(eligible) > 0L &&
               all(eligible$SourceIdentityAvailable == "Yes"),
             "pass", "fail")
    ),
    observed = c(as.character(pilot_record_limit), as.character(nrow(eligible))),
    expected = c("finite positive per-provider/species cap",
                 paste0("1-", config$pilot_enrichment_compounds,
                        " source-identified compounds")),
    stringsAsFactors = FALSE
  ))
  if (nrow(eligible) > 0) {
    enrichment = runCategorateEnrichmentBatches(
      eligible, out_dir = file.path(stage_dir, "research_enrichment"),
      cache_dir = file.path(config$cache_dir, "pilot_enrichment"),
      detail = "research", batch_size = config$compound_batch_size,
      pubchem_throttle = config$pubchem_throttle,
      kegg_throttle = config$kegg_throttle, assay_detail_limit = 0,
      cache = TRUE, resume = config$resume,
      overwrite = config$overwrite,
      require_source_identity = TRUE,
      request_fun = config$compound_request_fun,
      enrichment_fun = config$enrichment_fun,
      periodic_cooldown_batches = Inf,
      progress = config$progress
    )
    enrichment_file = file.path(stage_dir, "pilot_enrichment_result.rds")
    .plant_atomic_save_rds(enrichment, enrichment_file)
    enrichment_artifacts = enrichment_file
    enrichment_exit = attr(enrichment, "exit_status", exact = TRUE)
    if (identical(enrichment_exit, 75L)) {
      return(list(status = "paused_service_busy", exit_status = 75L,
                  message = "Pilot enrichment paused after a service-busy response.",
                  artifacts = c(pilot_file, result_file, enrichment_file)))
    }
    if (nrow(enrichment$RetryQueue) > 0) {
      gate = rbind(gate, data.frame(
        check = "pilot_enrichment_complete", status = "fail",
        observed = nrow(enrichment$RetryQueue), expected = 0,
        stringsAsFactors = FALSE
      ))
    }
  }
  gate_file = file.path(stage_dir, "pilot_quality_gate.csv")
  .plant_atomic_write_csv(gate, gate_file)
  if (any(gate$status == "fail")) {
    stop("All-provider pilot failed: ",
         paste(gate$check[gate$status == "fail"], collapse = ", "),
         call. = FALSE)
  }
  list(status = "completed", exit_status = 0L,
       message = paste("All-provider pilot passed for", nrow(pilot),
                       "species."),
       artifacts = c(pilot_file, result_file, gate_file,
                     enrichment_artifacts))
}

.plant_panel_pilot_record_limit = function(max_provider_records,
                                            enrichment_compounds) {
  target = max(25L, as.integer(enrichment_compounds[[1L]]))
  configured = suppressWarnings(as.numeric(max_provider_records[[1L]]))
  if (is.finite(configured)) {
    return(as.integer(min(configured, target)))
  }
  as.integer(target)
}

.plant_panel_cache_only_request = function(...) {
  stop("Pilot cache-only replay attempted a live provider request.",
       call. = FALSE)
}

.plant_panel_clear_derived_exports = function(path) {
  if (!dir.exists(path)) return(invisible(character()))
  exports = list.files(path, full.names = TRUE, recursive = FALSE,
                       all.files = TRUE, no.. = TRUE)
  exports = exports[!dir.exists(exports)]
  if (length(exports) > 0L) unlink(exports, force = TRUE)
  invisible(exports)
}

.plant_panel_pilot_enrichment_input = function(source_identity, limit) {
  empty = .categorate_enrichment_query_map(
    data.frame(compound_name = character(), stringsAsFactors = FALSE),
    require_source_identity = TRUE
  )
  if (!is.data.frame(source_identity) || nrow(source_identity) < 1L) {
    return(empty)
  }
  query_map = .categorate_enrichment_query_map(
    source_identity, require_source_identity = TRUE
  )
  query_map = query_map[
    query_map$EnrichmentEligible == "Yes" &
      query_map$SourceIdentityAvailable == "Yes", , drop = FALSE
  ]
  if (nrow(query_map) < 1L) return(query_map)
  identity_priority = match(
    query_map$QueryType, c("source_cid", "source_inchikey")
  )
  identity_priority[is.na(identity_priority)] = 99L
  query_map = query_map[order(
    identity_priority, query_map$compound_name_clean,
    query_map$PubChemQuery, query_map$compound_id
  ), , drop = FALSE]
  query_map = query_map[!duplicated(query_map$PubChemQuery), , drop = FALSE]
  query_map = utils::head(query_map, as.integer(limit[[1L]]))
  row.names(query_map) = NULL
  query_map
}

.plant_panel_stage_discovery = function(context) {
  .plant_panel_require_resource_stage(context)
  config = context$config
  input = .plant_panel_read_input(config$plant_csv, config$species_col)
  indexes = .plant_panel_provider_indexes(context)
  root = file.path(context$paths$stages, "discovery")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  results = list()
  artifacts = character()
  for (provider in config$sources) {
    stage_sources = provider
    supplied = config$provider_results
    if (provider == "pubtator" && "pubmed" %in% config$sources) {
      pubmed_result = results[["pubmed"]]
      if (!is.null(pubmed_result)) {
        supplied = supplied %||% list()
        supplied$pubmed = list(
          PlantCompoundOccurrences = pubmed_result$PlantCompoundOccurrences[
            pubmed_result$PlantCompoundOccurrences$source_database ==
              "PubMed", , drop = FALSE
          ],
          LiteratureCandidates = pubmed_result$LiteratureCandidates[
            pubmed_result$LiteratureCandidates$source_database ==
              "PubMed", , drop = FALSE
          ]
        )
        stage_sources = c("pubmed", "pubtator")
      }
    }
    provider_dir = file.path(root, provider, "direct_species")
    .plant_panel_clear_derived_exports(provider_dir)
    result = runPlantPhytochemistryBatch(
      plants = input, sources = stage_sources,
      taxon_fallback = "species", out_dir = provider_dir,
      cache_dir = config$cache_dir,
      species_chunk_size = config$species_chunk_size,
      compound_resolution_profile = "none", cache = TRUE,
      throttle = ifelse(provider == "knapsack",
                        config$knapsack_throttle,
                        config$provider_throttle),
      ncbi_email = Sys.getenv("NCBI_EMAIL", ""),
      ncbi_tool = Sys.getenv("NCBI_TOOL", "uafR"),
      ncbi_api_key = Sys.getenv("NCBI_API_KEY", ""),
      max_pubmed_records = config$max_pubmed_records,
      max_provider_records = config$max_provider_records,
      lotus_index = indexes$lotus, provider_indexes = indexes,
      require_all_providers = config$require_all_providers,
      request_timeout = config$request_timeout,
      provider_results = supplied,
      request_fun = config$request_fun,
      pubtator_request_fun = config$pubtator_request_fun,
      resume = config$resume, progress = config$progress,
      overwrite = config$overwrite, stop_on_error = FALSE,
      allow_large_live_run = TRUE, service_busy_pause_threshold = 1
    )
    result_file = file.path(root, provider,
                            paste0(provider, "_direct_result.rds"))
    .plant_atomic_save_rds(result, result_file)
    artifacts = c(artifacts, result_file,
                  file.path(provider_dir, "discovery_chunk_manifest.csv"),
                  file.path(provider_dir, "retry_queue.csv"))
    results[[provider]] = result
    run_status = .uaf_first_non_empty_text(
      result$BatchRunManifest$run_status, "incomplete"
    )
    if (!identical(run_status, "completed")) {
      rate_limited = any(result$BatchChunkManifest$status == "rate_limited")
      return(list(
        status = ifelse(rate_limited, "paused_service_busy",
                        "incomplete"),
        exit_status = ifelse(rate_limited, 75L, 1L),
        message = paste(provider, "discovery is incomplete; retry rows were preserved."),
        artifacts = artifacts
      ))
    }
  }

  fallback_results = .plant_panel_run_fallbacks(
    context, input, indexes, results
  )
  results = c(results, fallback_results$results)
  artifacts = c(artifacts, fallback_results$artifacts)
  merged = mergePlantPhytochemistryResults(results)
  merged_file = file.path(root, "plant_phytochemistry_discovery_merged.rds")
  .plant_atomic_save_rds(merged, merged_file)
  artifacts = c(artifacts, merged_file)
  table_map = list(
    PlantQueries = merged$PlantQueries,
    PlantQueryAliases = merged$PlantQueryAliases,
    ProviderDiagnostics = merged$ProviderDiagnostics,
    ProviderQueryAccounting = merged$ProviderQueryAccounting,
    ProviderResourceManifest = merged$ProviderResourceManifest,
    PlantCompoundOccurrences = merged$PlantCompoundOccurrences,
    PlantContextEvidence = merged$PlantContextEvidence,
    ProviderContextAudit = merged$ProviderContextAudit,
    LiteratureCandidates = merged$LiteratureCandidates,
    SourceCompoundIdentity = merged$SourceCompoundIdentity
  )
  for (name in names(table_map)) {
    path = file.path(root, paste0(name, ".csv"))
    .plant_atomic_write_csv(table_map[[name]], path)
    artifacts = c(artifacts, path)
  }
  gate = .plant_panel_provider_gate(merged, config$sources, nrow(input))
  gate_file = file.path(root, "discovery_quality_gate.csv")
  .plant_atomic_write_csv(gate, gate_file)
  artifacts = c(artifacts, gate_file)
  if (any(gate$status == "fail")) {
    stop("Full discovery failed provider accounting: ",
         paste(gate$check[gate$status == "fail"], collapse = ", "),
         call. = FALSE)
  }
  list(status = "completed", exit_status = 0L,
       message = paste("All configured providers accounted for",
                       nrow(input), "species."), artifacts = artifacts)
}

.plant_panel_stage_identity = function(context) {
  config = context$config
  discovery_file = .plant_panel_required_file(
    file.path(context$paths$stages, "discovery",
              "plant_phytochemistry_discovery_merged.rds"),
    "Run discovery before identity resolution."
  )
  result = readRDS(discovery_file)
  identity_selection = .plant_panel_identity_occurrences(result)
  identity_input = result
  identity_input$PlantCompoundOccurrences = identity_selection$included
  resolved_identity = resolvePlantCompoundIdentities(
    identity_input, cache = TRUE,
    cache_dir = file.path(config$cache_dir, "compound_identity"),
    throttle = config$pubchem_throttle,
    batch_size = config$compound_batch_size, resume = config$resume,
    progress = config$progress,
    compound_request_fun = config$compound_request_fun,
    lotus_index = .plant_panel_provider_indexes(context)$lotus
  )
  identity_categorate = attr(resolved_identity, "CategorateResult",
                             exact = TRUE)
  source_identity = attr(resolved_identity, "SourceCompoundIdentity",
                         exact = TRUE)
  resolution = .plant_merge_compound_resolution(
    resolved_identity, identity_selection$excluded_resolution
  )
  result$CompoundResolution = resolution
  result$CategorateResult = identity_categorate
  if (is.data.frame(source_identity)) {
    result$SourceCompoundIdentity = .plant_merge_source_identity_tables(
      list(result$SourceCompoundIdentity, source_identity)
    )
  }
  result$CompoundIdentityReview = plantCompoundIdentityReviewTable(result)
  result = .plant_panel_refresh_derived(result)
  stage_dir = file.path(context$paths$stages, "identity")
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  result_file = file.path(stage_dir, "plant_phytochemistry_identity.rds")
  resolution_file = file.path(stage_dir, "CompoundResolution.csv")
  review_file = file.path(stage_dir, "CompoundIdentityReview.csv")
  selection_file = file.path(stage_dir, "IdentityOccurrenceSelection.csv")
  exclusion_file = file.path(stage_dir, "IdentityCandidateExclusions.csv")
  .plant_atomic_save_rds(result, result_file)
  .plant_atomic_write_csv(resolution, resolution_file)
  .plant_atomic_write_csv(result$CompoundIdentityReview, review_file)
  .plant_atomic_write_csv(identity_selection$included, selection_file)
  .plant_atomic_write_csv(identity_selection$excluded, exclusion_file)
  unique_compounds = unique(.uaf_non_empty(
    result$PlantCompoundOccurrences$compound_name_clean
  ))
  represented = unique(.uaf_non_empty(resolution$compound_name_clean))
  summary = data.frame(
    unique_discovered_compounds = length(unique_compounds),
    identity_rows = nrow(resolution),
    resolved_compounds = sum(resolution$resolved %in% TRUE),
    unresolved_compounds = sum(!(resolution$resolved %in% TRUE)),
    review_required_rows = nrow(result$CompoundIdentityReview),
    all_compounds_accounted = .uaf_yes_no(
      length(setdiff(unique_compounds, represented)) == 0
    ),
    stringsAsFactors = FALSE
  )
  summary_file = file.path(stage_dir, "identity_summary.csv")
  .plant_atomic_write_csv(summary, summary_file)
  if (summary$all_compounds_accounted != "Yes") {
    stop("Identity resolution omitted one or more discovered compound keys.",
         call. = FALSE)
  }
  busy = any(vapply(.uaf_non_empty(resolution$notes),
                    .plant_service_busy_message, logical(1)))
  if (busy && sum(resolution$resolved %in% TRUE) == 0) {
    return(list(status = "paused_service_busy", exit_status = 75L,
                message = "Compound identity stage paused after service-busy responses.",
                artifacts = c(result_file, resolution_file, review_file,
                              summary_file, selection_file, exclusion_file)))
  }
  list(status = "completed", exit_status = 0L,
       message = paste(sum(resolution$resolved %in% TRUE), "of",
                       nrow(resolution), "compound identities resolved; unresolved rows retained."),
       artifacts = c(result_file, resolution_file, review_file, summary_file,
                     selection_file, exclusion_file))
}

.plant_panel_identity_occurrences = function(result) {
  occurrences = .plant_normalize_occurrences(
    result$PlantCompoundOccurrences
  )
  if (nrow(occurrences) < 1L) {
    return(list(
      included = occurrences,
      excluded = occurrences,
      excluded_resolution = .uaf_empty_table(
        .plant_compound_resolution_cols()
      )
    ))
  }
  status = tolower(.uaf_squish_text(occurrences$occurrence_status))
  tier = tolower(.uaf_squish_text(occurrences$evidence_tier))
  candidate = status == "candidate" |
    grepl("candidate", tier, fixed = TRUE) |
    grepl("co.?mention", tier, perl = TRUE)
  has_name = !is.na(.uaf_squish_text(occurrences$compound_name)) &
    .uaf_squish_text(occurrences$compound_name) != ""
  include = has_name & !candidate & status %in%
    c("direct_reported", "curated_reported", "taxon_fallback")
  included = occurrences[include, , drop = FALSE]
  excluded = occurrences[!include, , drop = FALSE]
  excluded_resolution = .plant_panel_explicit_unresolved(excluded)
  list(included = included, excluded = excluded,
       excluded_resolution = excluded_resolution)
}

.plant_panel_explicit_unresolved = function(occurrences) {
  occurrences = .plant_normalize_occurrences(occurrences)
  keys = unique(.uaf_non_empty(occurrences$compound_name_clean))
  if (length(keys) < 1L) {
    return(.uaf_empty_table(.plant_compound_resolution_cols()))
  }
  rows = lapply(keys, function(key) {
    part = occurrences[
      occurrences$compound_name_clean == key, , drop = FALSE
    ]
    names_observed = .uaf_non_empty(part$compound_name)
    name = if (length(names_observed) > 0L) names_observed[[1L]] else key
    data.frame(
      compound_name = name,
      compound_name_clean = key,
      query_count = nrow(part),
      resolved = FALSE,
      CID = NA_character_,
      InChIKey = NA_character_,
      SMILES = NA_character_,
      MolecularFormula = NA_character_,
      resolution_source = "not_submitted_candidate_evidence",
      notes = paste(
        "Retained unresolved without PubChem lookup because available rows",
        "were candidate/co-mention evidence or lacked an eligible occurrence",
        "status. Review source evidence before identity promotion."
      ),
      stringsAsFactors = FALSE
    )
  })
  .plant_bind_tables(rows, .plant_compound_resolution_cols())
}

.plant_panel_stage_research = function(context) {
  config = context$config
  identity_file = .plant_panel_required_file(
    file.path(context$paths$stages, "identity",
              "plant_phytochemistry_identity.rds"),
    "Run identity before research enrichment."
  )
  result = readRDS(identity_file)
  selection = .plant_panel_enrichment_selection(result)
  stage_dir = file.path(context$paths$stages, "research_enrichment")
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  selection_file = file.path(stage_dir, "research_enrichment_selection.csv")
  exclusion_file = file.path(stage_dir, "research_enrichment_exclusions.csv")
  .plant_atomic_write_csv(selection$included, selection_file)
  .plant_atomic_write_csv(selection$excluded, exclusion_file)
  if (nrow(selection$included) < 1) {
    stop("No conservatively resolved compounds are eligible for research enrichment.",
         call. = FALSE)
  }
  enrichment = runCategorateEnrichmentBatches(
    compounds = selection$included,
    out_dir = stage_dir,
    cache_dir = file.path(config$cache_dir, "categorate_research"),
    detail = "research", batch_size = config$compound_batch_size,
    pubchem_throttle = config$pubchem_throttle,
    kegg_throttle = config$kegg_throttle, assay_detail_limit = 0,
    cache = TRUE, resume = config$resume, overwrite = config$overwrite,
    request_fun = config$compound_request_fun,
    enrichment_fun = config$enrichment_fun,
    periodic_cooldown_batches = 10, periodic_cooldown_seconds = 60,
    stop_on_error = FALSE, progress = config$progress
  )
  enrichment_file = file.path(stage_dir, "research_enrichment_result.rds")
  .plant_atomic_save_rds(enrichment, enrichment_file)
  exit_status = attr(enrichment, "exit_status", exact = TRUE)
  if (identical(exit_status, 75L)) {
    return(list(status = "paused_service_busy", exit_status = 75L,
                message = "Research enrichment paused; completed batches and caches were preserved.",
                artifacts = c(selection_file, exclusion_file,
                              enrichment_file)))
  }
  if (nrow(enrichment$RetryQueue) > 0) {
    return(list(status = "incomplete", exit_status = 1L,
                message = "Research enrichment has retryable incomplete batches.",
                artifacts = c(selection_file, exclusion_file,
                              enrichment_file)))
  }
  combined = combineCategorateTables(
    enrichment$Batches, include_batch_metadata = FALSE,
    include_empty = TRUE, min_property_ratio = 0
  )
  result$CategorateResult = combined
  result = .plant_panel_refresh_derived(result)
  result_file = file.path(stage_dir,
                          "plant_phytochemistry_research_enriched.rds")
  combined_file = file.path(stage_dir, "research_combined_tables.rds")
  .plant_atomic_save_rds(result, result_file)
  .plant_atomic_save_rds(combined, combined_file)
  list(status = "completed", exit_status = 0L,
       message = paste(nrow(selection$included),
                       "resolved compounds completed research enrichment."),
       artifacts = c(selection_file, exclusion_file, enrichment_file,
                     result_file, combined_file,
                     file.path(stage_dir, "research_batch_manifest.csv"),
                     file.path(stage_dir, "research_retry_queue.csv")))
}

.plant_panel_stage_full = function(context) {
  config = context$config
  result_file = .plant_panel_required_file(
    file.path(context$paths$stages, "research_enrichment",
              "plant_phytochemistry_research_enriched.rds"),
    "Run research enrichment before full enrichment."
  )
  result = readRDS(result_file)
  priority = .plant_panel_full_priority(result,
                                        config$full_enrichment_limit)
  stage_dir = file.path(context$paths$stages, "full_enrichment")
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  priority_file = file.path(stage_dir, "full_enrichment_priority.csv")
  .plant_atomic_write_csv(priority, priority_file)
  if (nrow(priority) < 1) {
    note = file.path(stage_dir, "NO_PRIORITY_COMPOUNDS.txt")
    .plant_atomic_write_text(
      "No compounds met the deterministic direct-evidence priority rule.",
      note
    )
    return(list(status = "completed", exit_status = 0L,
                message = "No compounds met the full-detail priority rule; research records were preserved.",
                artifacts = c(priority_file, note)))
  }
  enrichment = runCategorateEnrichmentBatches(
    compounds = priority,
    out_dir = stage_dir,
    cache_dir = file.path(config$cache_dir, "categorate_full"),
    detail = "full", batch_size = config$compound_batch_size,
    pubchem_throttle = config$pubchem_throttle,
    kegg_throttle = config$kegg_throttle, assay_detail_limit = 10,
    cache = TRUE, resume = config$resume, overwrite = config$overwrite,
    request_fun = config$compound_request_fun,
    enrichment_fun = config$enrichment_fun,
    periodic_cooldown_batches = 10, periodic_cooldown_seconds = 60,
    stop_on_error = FALSE, progress = config$progress
  )
  enrichment_file = file.path(stage_dir, "full_enrichment_result.rds")
  .plant_atomic_save_rds(enrichment, enrichment_file)
  exit_status = attr(enrichment, "exit_status", exact = TRUE)
  if (identical(exit_status, 75L)) {
    return(list(status = "paused_service_busy", exit_status = 75L,
                message = "Full enrichment paused; completed batches and caches were preserved.",
                artifacts = c(priority_file, enrichment_file)))
  }
  if (nrow(enrichment$RetryQueue) > 0) {
    return(list(status = "incomplete", exit_status = 1L,
                message = "Full enrichment has retryable incomplete batches.",
                artifacts = c(priority_file, enrichment_file)))
  }
  list(status = "completed", exit_status = 0L,
       message = paste(nrow(priority), "priority compounds completed full enrichment."),
       artifacts = c(priority_file, enrichment_file,
                     file.path(stage_dir, "full_batch_manifest.csv"),
                     file.path(stage_dir, "full_retry_queue.csv")))
}

.plant_panel_stage_tanimoto_handoff = function(context) {
  release_validation = .plant_panel_validate_release_artifacts(context$config)
  if (any(release_validation$Checks$status == "fail")) {
    stop("Tanimoto handoff release-artifact validation failed: ",
         paste(release_validation$Checks$check[
           release_validation$Checks$status == "fail"
         ], collapse = ", "), call. = FALSE)
  }
  result_file = .plant_panel_required_file(
    file.path(context$paths$stages, "research_enrichment",
              "plant_phytochemistry_research_enriched.rds"),
    "Run research enrichment before preparing the Tanimoto handoff."
  )
  result = .plant_panel_refresh_derived(readRDS(result_file))
  handoff = preparePlantTanimotoInput(
    result, level = "species",
    occurrence_status = c("direct_reported", "curated_reported"),
    analysis_ready = TRUE, min_confidence = "medium",
    include_review_required = FALSE,
    out_dir = context$paths$handoff,
    overwrite = context$config$overwrite,
    strict = TRUE
  )
  all_species = .plant_panel_read_input(
    context$config$plant_csv, context$config$species_col
  )$species
  species_universe_file = file.path(context$paths$handoff,
                                    "plant_species_universe.csv")
  .plant_atomic_write_csv(
    data.frame(species = unique(all_species), stringsAsFactors = FALSE),
    species_universe_file
  )
  expected_pairs = choose(length(unique(all_species)), 2)
  estimate = data.frame(
    plant_count = length(unique(all_species)),
    plants_with_tanimoto_membership =
      length(unique(handoff$PlantCompoundMembership$species)),
    membership_row_count = nrow(handoff$PlantCompoundMembership),
    unique_structure_count = nrow(handoff$CompoundInput),
    expected_unordered_plant_pair_count = format(expected_pairs,
                                                  scientific = FALSE),
    expected_compound_pair_count = format(
      choose(nrow(handoff$CompoundInput), 2), scientific = FALSE
    ),
    minimum_server_free_disk_multiplier = 2.5,
    stringsAsFactors = FALSE
  )
  estimate_file = file.path(context$paths$handoff,
                            "tanimoto_production_estimate.csv")
  .plant_atomic_write_csv(estimate, estimate_file)
  runner = .plant_panel_tanimoto_runner_path()
  runner_copy = file.path(context$paths$handoff,
                          "run_plant_tanimoto_server.R")
  .plant_panel_copy_file(runner, runner_copy, overwrite = TRUE)
  handoff_support_sources = c(
    context$config$plant_csv,
    file.path(context$paths$inputs, "run_contract.json"),
    file.path(context$paths$inputs, "input_manifest.csv")
  )
  handoff_support_sources = handoff_support_sources[
    file.exists(handoff_support_sources)
  ]
  handoff_support_artifacts = file.path(
    context$paths$handoff, basename(handoff_support_sources)
  )
  for (i in seq_along(handoff_support_sources)) {
    .plant_panel_copy_file(handoff_support_sources[[i]],
                           handoff_support_artifacts[[i]], overwrite = TRUE)
  }
  release_artifacts = character()
  release_manifest_copy = NA_character_
  source_tarball_copy = NA_character_
  if (!is.na(context$config$release_manifest)) {
    release_manifest_copy = file.path(
      context$paths$handoff, basename(context$config$release_manifest)
    )
    .plant_panel_copy_file(context$config$release_manifest,
                           release_manifest_copy, overwrite = TRUE)
    release_artifacts = c(release_artifacts, release_manifest_copy)
  }
  if (!is.na(context$config$source_tarball)) {
    source_tarball_copy = file.path(
      context$paths$handoff, basename(context$config$source_tarball)
    )
    .plant_panel_copy_file(context$config$source_tarball,
                           source_tarball_copy, overwrite = TRUE)
    release_artifacts = c(release_artifacts, source_tarball_copy)
  }
  source_manifest = .plant_panel_source_snapshot(
    context, runner_copy, species_universe_file, release_artifacts,
    handoff_support_artifacts
  )
  source_manifest_file = file.path(context$paths$handoff,
                                   "source_release_manifest.csv")
  .plant_atomic_write_csv(source_manifest, source_manifest_file)
  input_file = file.path(
    context$paths$handoff,
    "plant_compound_membership_tanimoto_ready.csv"
  )
  input_manifest = file.path(context$paths$handoff,
                             "tanimoto_input_export_manifest.csv")
  server_workdir = Sys.getenv("UAFR_SERVER_WORKDIR", "")
  remote_input = if (nzchar(server_workdir)) {
    file.path(server_workdir, basename(input_file))
  } else {
    "<UAFR_SERVER_WORKDIR>/plant_compound_membership_tanimoto_ready.csv"
  }
  release_arg = if (!is.na(release_manifest_copy)) {
    paste("--release-manifest", shQuote(basename(release_manifest_copy)))
  } else ""
  commands = c(
    "# Run these commands on the configured compute server after transferring the entire handoff directory.",
    if (!is.na(source_tarball_copy)) {
      paste("R CMD INSTALL", shQuote(basename(source_tarball_copy)))
    } else {
      "# Install the exact uafR production package before running live modes."
    },
    paste("Rscript run_plant_tanimoto_server.R --input",
          shQuote(remote_input), "--out-dir",
          shQuote(file.path(ifelse(nzchar(server_workdir), server_workdir,
                                   "<UAFR_SERVER_WORKDIR>"),
                            "tanimoto_output")),
          "--cache-dir",
          shQuote(file.path(ifelse(nzchar(server_workdir), server_workdir,
                                   "<UAFR_SERVER_WORKDIR>"),
                            "pubchem_tanimoto_cache")),
          "--manifest", shQuote(basename(input_manifest)),
          "--species-universe", shQuote(basename(species_universe_file)),
          release_arg,
          "--mode preflight"),
    "export UAFR_CONFIRM_SERVER_TANIMOTO=YES",
    paste("Rscript run_plant_tanimoto_server.R --input",
          shQuote(remote_input), "--out-dir tanimoto_output",
          "--cache-dir pubchem_tanimoto_cache --manifest",
          shQuote(basename(input_manifest)),
          "--species-universe", shQuote(basename(species_universe_file)),
          release_arg,
          "--mode smoke --smoke-structures 25 --throttle 1.1"),
    paste("Rscript run_plant_tanimoto_server.R --input",
          shQuote(remote_input), "--out-dir tanimoto_output",
          "--cache-dir pubchem_tanimoto_cache --manifest",
          shQuote(basename(input_manifest)),
          "--species-universe", shQuote(basename(species_universe_file)),
          release_arg,
          "--mode summary --throttle 1.1"),
    paste("Rscript run_plant_tanimoto_server.R --input",
          shQuote(remote_input), "--out-dir tanimoto_output",
          "--cache-dir pubchem_tanimoto_cache --manifest",
          shQuote(basename(input_manifest)),
          "--species-universe", shQuote(basename(species_universe_file)),
          release_arg,
          "--mode full --full-pairs true --pair-shard-rows 5000000",
          "--throttle 1.1")
  )
  command_file = file.path(context$paths$handoff,
                           "RUN_TANIMOTO_SERVER.txt")
  .plant_atomic_write_text(commands, command_file)
  list(status = "completed", exit_status = 0L,
       message = paste(nrow(handoff$CompoundInput),
                       "verified structures prepared for server Tanimoto analysis."),
       artifacts = c(input_file, input_manifest, species_universe_file,
                     estimate_file, runner_copy,
                     source_manifest_file, command_file,
                     release_artifacts, handoff_support_artifacts))
}

.plant_panel_stage_finalize = function(context) {
  config = context$config
  server_root = .uaf_first_non_empty_text(
    config$server_results_dir,
    file.path(context$paths$root, "tanimoto_server_output")
  )
  full_dir = if (dir.exists(file.path(server_root, "full"))) {
    file.path(server_root, "full")
  } else server_root
  summary_file = file.path(full_dir, "plant_pair_tanimoto_summary.csv")
  if (!file.exists(summary_file)) {
    note = paste(
      "Finalization is waiting for completed server Tanimoto outputs.",
      "Expected:", summary_file,
      "Run the commands in tanimoto_handoff/RUN_TANIMOTO_SERVER.txt, then",
      "rerun mode = 'finalize' with server_results_dir set to the downloaded",
      "server output directory."
    )
    awaiting_file = file.path(context$paths$root,
                              "AWAITING_TANIMOTO_SERVER.txt")
    .plant_atomic_write_text(note, awaiting_file)
    return(list(status = "awaiting_server", exit_status = 2L,
                message = note, artifacts = awaiting_file))
  }
  research_dir = file.path(context$paths$stages, "research_enrichment")
  enrichment_file = .plant_panel_required_file(
    file.path(research_dir, "research_enrichment_result.rds"),
    "Research enrichment output is missing."
  )
  research = readRDS(enrichment_file)
  handoff_membership = .plant_panel_required_file(
    file.path(context$paths$handoff,
              "plant_compound_membership_tanimoto_ready.csv"),
    "Tanimoto handoff membership is missing."
  )
  identity_result = readRDS(.plant_panel_required_file(
    file.path(context$paths$stages, "identity",
              "plant_phytochemistry_identity.rds"),
    "Identity result is missing."
  ))
  fingerprints = file.path(full_dir, "pubchem_fingerprint_audit.csv")
  if (!file.exists(fingerprints)) fingerprints = NULL
  pair_files = list.files(
    full_dir,
    pattern = "^group_compound_pair_tanimoto.*[.]csv([.]gz)?$",
    full.names = TRUE
  )
  scope_summary_file = file.path(
    full_dir, "comparable_scope_tanimoto_summary.csv"
  )
  group_summary_file = file.path(
    full_dir, "comparable_group_tanimoto_summary.csv"
  )
  if (!file.exists(scope_summary_file) || !file.exists(group_summary_file)) {
    stop(paste(
      "Completed server output is missing comparable chemistry summaries.",
      "Use the bundled production Tanimoto runner and rerun full mode."
    ), call. = FALSE)
  }
  metadata = .plant_panel_read_input(config$metadata_csv,
                                     config$species_col)
  plant_list = .plant_panel_read_input(config$plant_csv,
                                       config$species_col)
  file_references = .plant_panel_file_references(context, server_root)
  manifest = exportPlantChemistryAnalysisBundle(
    categorate_batches = research$Batches,
    path = context$paths$bundle,
    plant_membership = handoff_membership,
    species_pair_tanimoto = summary_file,
    resolved_compounds = identity_result$CompoundResolution,
    pubchem_fingerprints = fingerprints,
    plant_compound_pair_tanimoto = NULL,
    file_references = file_references,
    plant_list = plant_list,
    metadata = metadata,
    project_id = config$project_id,
    format = "csv", overwrite = TRUE, finalize = TRUE,
    include_feature_exports = TRUE,
    include_comparable_tanimoto = FALSE,
    validate_export = TRUE
  )
  .plant_panel_copy_file(
    scope_summary_file,
    file.path(context$paths$bundle,
              "22_ComparableScopeTanimotoSummary.csv"),
    overwrite = TRUE
  )
  .plant_panel_copy_file(
    group_summary_file,
    file.path(context$paths$bundle,
              "23_ComparableGroupTanimotoSummary.csv"),
    overwrite = TRUE
  )
  discovery = readRDS(.plant_panel_required_file(
    file.path(context$paths$stages, "discovery",
              "plant_phytochemistry_discovery_merged.rds"),
    "Merged discovery result is missing."
  ))
  .plant_panel_add_discovery_bundle_tables(context$paths$bundle, discovery)
  manifest = finalizePlantChemistryAnalysisBundle(
    context$paths$bundle, plant_list = plant_list, metadata = metadata,
    plant_compound_pair_tanimoto = NULL,
    project_id = config$project_id, overwrite = TRUE,
    include_feature_exports = TRUE,
    include_comparable_tanimoto = FALSE,
    validate_export = TRUE
  )
  validation = validatePlantChemistryAnalysisBundle(context$paths$bundle)
  validation_file = file.path(context$paths$bundle,
                              "99_FinalBundleValidation.csv")
  .plant_atomic_write_csv(validation$Summary, validation_file)
  plant_pairs = utils::read.csv(summary_file, stringsAsFactors = FALSE,
                                check.names = FALSE)
  scope_pairs = utils::read.csv(scope_summary_file,
                                stringsAsFactors = FALSE,
                                check.names = FALSE)
  group_pairs = utils::read.csv(group_summary_file,
                                stringsAsFactors = FALSE,
                                check.names = FALSE)
  server_validation_file = file.path(full_dir,
                                     "server_tanimoto_validation.csv")
  server_validation = if (file.exists(server_validation_file)) {
    utils::read.csv(server_validation_file, stringsAsFactors = FALSE,
                    check.names = FALSE)
  } else data.frame()
  pair_manifest_file = file.path(server_root,
                                 "server_tanimoto_output_manifest.csv")
  pair_manifest = if (file.exists(pair_manifest_file)) {
    utils::read.csv(pair_manifest_file, stringsAsFactors = FALSE,
                    check.names = FALSE)
  } else data.frame()
  shard_rows = if (nrow(pair_manifest) > 0L &&
      all(c("file", "row_count") %in% names(pair_manifest))) {
    sum(suppressWarnings(as.numeric(pair_manifest$row_count[
      grepl("group_compound_pair_tanimoto", pair_manifest$file,
            fixed = TRUE)
    ])), na.rm = TRUE)
  } else 0
  expected_pairs = choose(nrow(plant_list), 2)
  reconciliation = data.frame(
    check = c("species_accounted", "provider_retry_queue_empty",
              "plant_pair_universe_complete",
              "unsupported_pairs_explicit", "server_validation_passed",
              "pairwise_shards_manifested",
              "comparable_scope_summary_available",
              "comparable_group_summary_available",
              "bundle_export_ready"),
    status = c(
      ifelse(nrow(plant_list) == length(unique(plant_list$species)),
             "pass", "fail"),
      ifelse(all(discovery$ProviderQueryAccounting$retry_required == "No"),
             "pass", "fail"),
      ifelse(nrow(plant_pairs) == expected_pairs, "pass", "fail"),
      ifelse("support_status" %in% names(plant_pairs) &&
               all(plant_pairs$support_status %in%
                     c("computed", "insufficient_support")),
             "pass", "fail"),
      ifelse(nrow(server_validation) > 0L &&
               all(server_validation$status == "pass"), "pass", "fail"),
      ifelse(length(pair_files) > 0L && shard_rows > 0,
             "pass", "fail"),
      ifelse(nrow(scope_pairs) > 0L, "pass", "fail"),
      ifelse(nrow(group_pairs) > 0L, "pass", "fail"),
      ifelse(.plant_panel_bundle_ready(validation), "pass", "fail")
    ),
    observed = c(nrow(plant_list),
                 sum(discovery$ProviderQueryAccounting$retry_required == "Yes"),
                 nrow(plant_pairs),
                 if ("support_status" %in% names(plant_pairs)) {
                   paste(table(plant_pairs$support_status), collapse = "; ")
                 } else "missing support_status",
                 if (nrow(server_validation) > 0L) {
                   paste(unique(server_validation$status), collapse = "; ")
                 } else "missing",
                 paste(length(pair_files), "shards;", shard_rows, "rows"),
                 nrow(scope_pairs), nrow(group_pairs),
                 .plant_panel_validation_status(validation)),
    expected = c(nrow(plant_list), 0, expected_pairs,
                 "computed or insufficient_support", "all pass",
                 ">=1 shard and >0 rows", ">0 rows", ">0 rows",
                 "export ready"),
    stringsAsFactors = FALSE
  )
  reconciliation_file = file.path(context$paths$bundle,
                                  "99b_FinalReconciliation.csv")
  .plant_atomic_write_csv(reconciliation, reconciliation_file)
  if (any(reconciliation$status == "fail")) {
    stop("Final production reconciliation failed: ",
         paste(reconciliation$check[reconciliation$status == "fail"],
               collapse = ", "), call. = FALSE)
  }
  list(status = "completed", exit_status = 0L,
       message = "Final plant chemistry analysis bundle passed reconciliation.",
       artifacts = c(file.path(context$paths$bundle,
                               "01_ExportManifest.csv"),
                     validation_file, reconciliation_file))
}

.plant_panel_run_fallbacks = function(context, input, indexes,
                                       direct_results) {
  config = context$config
  ranks = intersect(config$taxon_fallback, c("genus", "family"))
  providers = intersect(config$sources, c("lotus", "npass", "knapsack"))
  if (length(ranks) < 1 || length(providers) < 1) {
    return(list(results = list(), artifacts = character()))
  }
  out_results = list()
  artifacts = character()
  for (provider in providers) {
    no_hit_species = .plant_panel_provider_no_hit_species(
      direct_results[[provider]], provider
    )
    fallback_input = input[
      .plant_clean_name(input$species) %in%
        .plant_clean_name(no_hit_species), , drop = FALSE
    ]
    if (nrow(fallback_input) < 1L) next
    out_dir = file.path(context$paths$stages, "discovery", provider,
                        "fallback")
    .plant_panel_clear_derived_exports(out_dir)
    selection_file = file.path(
      context$paths$stages, "discovery", provider,
      paste0(provider, "_fallback_species.csv")
    )
    .plant_atomic_write_csv(fallback_input, selection_file)
    result = runPlantPhytochemistryBatch(
      plants = fallback_input, sources = provider,
      taxon_fallback = c("species", ranks), out_dir = out_dir,
      cache_dir = config$cache_dir,
      species_chunk_size = config$species_chunk_size,
      compound_resolution_profile = "none", cache = TRUE,
      throttle = ifelse(provider == "knapsack",
                        config$knapsack_throttle,
                        config$provider_throttle),
      max_pubmed_records = config$max_pubmed_records,
      max_provider_records = config$max_provider_records,
      lotus_index = indexes$lotus, provider_indexes = indexes,
      require_all_providers = config$require_all_providers,
      request_timeout = config$request_timeout,
      provider_results = config$provider_results,
      request_fun = config$request_fun,
      resume = config$resume, progress = config$progress,
      overwrite = config$overwrite, stop_on_error = FALSE,
      allow_large_live_run = TRUE
    )
    fallback = result
    keep = fallback$PlantCompoundOccurrences$matched_rank %in% ranks
    fallback$PlantCompoundOccurrences =
      fallback$PlantCompoundOccurrences[keep, , drop = FALSE]
    fallback$PlantContextEvidence = plantContextEvidence(
      fallback$PlantCompoundOccurrences
    )
    fallback = .plant_panel_refresh_derived(fallback)
    name = paste0(provider, "_fallback")
    file = file.path(context$paths$stages, "discovery", provider,
                     paste0(name, "_result.rds"))
    .plant_atomic_save_rds(fallback, file)
    fallback_csv = file.path(context$paths$stages, "discovery", provider,
                             paste0(name, "_occurrences.csv"))
    .plant_atomic_write_csv(fallback$PlantCompoundOccurrences, fallback_csv)
    out_results[[name]] = fallback
    artifacts = c(artifacts, selection_file, file, fallback_csv)
  }
  list(results = out_results, artifacts = artifacts)
}

.plant_panel_provider_no_hit_species = function(result, provider) {
  if (!is.list(result) || !is.data.frame(result$ProviderQueryAccounting)) {
    return(character())
  }
  accounting = result$ProviderQueryAccounting
  rows = tolower(.uaf_squish_text(accounting$provider)) ==
    tolower(provider) &
    tolower(.uaf_squish_text(accounting$query_status)) == "no_records"
  species = .uaf_squish_text(accounting$species)
  if ("query_plant" %in% names(accounting)) {
    missing = is.na(species) | species == ""
    species[missing] = .uaf_squish_text(accounting$query_plant[missing])
  }
  unique(.uaf_non_empty(species[rows]))
}

.plant_panel_refresh_derived = function(result) {
  result$SpeciesChemistrySummary = summarizePlantPhytochemistry(
    plant_compounds = result$PlantCompoundOccurrences,
    categorate_result = result$CategorateResult,
    compound_resolution = result$CompoundResolution,
    plant_queries = result$PlantQueries,
    provider_diagnostics = result$ProviderDiagnostics
  )
  result$SpeciesChemistryMatrix = plantPhytochemistryMatrix(
    result, level = "species", profile = "core", mode = "binary",
    min_confidence = "medium"
  )
  result$ChemistryComparability = plantChemistryComparability(
    result, min_confidence = "low"
  )
  result$ComparableChemistryMatrix = plantComparableChemistryMatrix(
    result, level = "species", comparison_scope = "all_classified",
    mode = "binary", min_comparability_confidence = "medium"
  )
  result$CompoundIdentityReview = plantCompoundIdentityReviewTable(result)
  result$Validation = validatePlantPhytochemistryResult(result)
  result
}

.plant_panel_enrichment_selection = function(result) {
  resolution = result$CompoundResolution
  review = plantCompoundIdentityReviewTable(result)
  review_keys = unique(.uaf_non_empty(review$compound_name_clean))
  resolved = resolution$resolved %in% TRUE
  identity = !is.na(.uaf_squish_text(resolution$CID)) |
    vapply(resolution$InChIKey, .uaf_is_inchikey, logical(1))
  review_ok = !resolution$compound_name_clean %in% review_keys
  keep = resolved & identity & review_ok
  included = resolution[keep, , drop = FALSE]
  excluded = resolution[!keep, , drop = FALSE]
  excluded$enrichment_exclusion_reason = vapply(
    which(!keep), function(i) {
      reason = character()
      if (!resolved[[i]]) reason = c(reason, "unresolved")
      if (!identity[[i]]) reason = c(reason, "no_verified_cid_or_inchikey")
      if (!review_ok[[i]]) reason = c(reason, "identity_review_required")
      paste(reason, collapse = "; ")
    }, character(1)
  )
  list(included = included, excluded = excluded)
}

.plant_panel_full_priority = function(result, limit) {
  resolution = result$CompoundResolution
  occurrences = result$PlantCompoundOccurrences
  direct = occurrences[
    occurrences$matched_rank == "species" &
      occurrences$occurrence_status %in% c("direct_reported",
                                            "curated_reported") &
      occurrences$analysis_ready == "Yes", , drop = FALSE
  ]
  keys = unique(.uaf_non_empty(resolution$compound_name_clean))
  rows = lapply(keys, function(key) {
    id = resolution[resolution$compound_name_clean == key, , drop = FALSE]
    occ = direct[direct$compound_name_clean == key, , drop = FALSE]
    if (nrow(id) < 1) return(NULL)
    source_count = length(unique(.uaf_non_empty(occ$source_database)))
    species_count = length(unique(.uaf_non_empty(occ$species)))
    structure = id$resolved[[1]] %in% TRUE &&
      (!is.na(.uaf_squish_text(id$CID[[1]])) ||
         .uaf_is_inchikey(id$InChIKey[[1]]))
    eligible = structure && nrow(occ) > 0 &&
      (source_count >= 2 || species_count >= 5)
    data.frame(
      compound_id = if (.uaf_is_inchikey(id$InChIKey[[1]])) {
        id$InChIKey[[1]]
      } else if (!is.na(.uaf_squish_text(id$CID[[1]]))) {
        paste0("cid_", id$CID[[1]])
      } else key,
      compound_name = id$compound_name[[1]],
      compound_name_clean = key,
      CID = id$CID[[1]], InChIKey = id$InChIKey[[1]],
      SMILES = id$SMILES[[1]], MolecularFormula = id$MolecularFormula[[1]],
      source_database_count = source_count,
      species_count = species_count,
      direct_evidence_count = nrow(occ), priority_eligible = eligible,
      stringsAsFactors = FALSE
    )
  })
  out = .plant_bind_tables(rows, c(
    "compound_id", "compound_name", "compound_name_clean", "CID",
    "InChIKey", "SMILES", "MolecularFormula", "source_database_count",
    "species_count", "direct_evidence_count", "priority_eligible"
  ))
  out = out[out$priority_eligible %in% TRUE, , drop = FALSE]
  out = out[order(-out$source_database_count, -out$species_count,
                  -out$direct_evidence_count, out$compound_id), , drop = FALSE]
  utils::head(out, limit)
}

.plant_panel_provider_indexes = function(context) {
  configured = .plant_panel_index_config(context$paths$resources)
  list(
    lotus = .uaf_first_non_empty_text(context$config$lotus_index,
                                      configured$lotus),
    npass = .uaf_first_non_empty_text(context$config$npass_index,
                                      configured$npass)
  )
}

.plant_panel_provider_gate = function(result, providers, species_count) {
  accounting = result$ProviderQueryAccounting
  diagnostics = result$ProviderDiagnostics
  rows = list()
  for (provider in providers) {
    x = accounting[accounting$provider == provider, , drop = FALSE]
    d = diagnostics[diagnostics$provider == provider, , drop = FALSE]
    rows[[length(rows) + 1L]] = data.frame(
      check = paste0(provider, "_species_accounting"),
      status = ifelse(nrow(x) == species_count &&
                        all(x$query_status %in% c("records", "no_records")) &&
                        all(x$retry_required == "No"), "pass", "fail"),
      observed = paste(nrow(x), "rows;",
                       sum(x$retry_required == "Yes"), "retry"),
      expected = paste(species_count, "rows; 0 retry"),
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] = data.frame(
      check = paste0(provider, "_diagnostic_status"),
      status = ifelse(nrow(d) > 0 &&
                        all(d$status %in% c("ok", "no_records")) &&
                        all(d$queried == "Yes") && all(d$available == "Yes"),
                      "pass", "fail"),
      observed = paste(unique(d$status), collapse = "; "),
      expected = "ok or no_records; queried and available",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.plant_panel_select_pilot = function(input, requested, count) {
  species = input$species
  requested = intersect(unique(.uaf_non_empty(requested)), species)
  selected = requested
  if (length(selected) < count) {
    preferred = c("Camellia sinensis", "Zea mays", "Salix nigra",
                  "Arabidopsis thaliana", "Acer negundo",
                  "Lotus corniculatus", "Coffea arabica",
                  "Nicotiana tabacum", "Quercus robur", "Glycine max")
    selected = unique(c(selected, intersect(preferred, species)))
  }
  if (length(selected) < count) {
    remaining = setdiff(species, selected)
    if (length(remaining) > 0) {
      positions = unique(round(seq(1, length(remaining),
                                   length.out = min(count - length(selected),
                                                    length(remaining)))))
      selected = unique(c(selected, remaining[positions]))
    }
  }
  input[match(utils::head(selected, count), species), , drop = FALSE]
}

.plant_panel_run_manifest_cols = function() {
  c("stage", "status", "started_at", "finished_at", "elapsed_seconds",
    "message", "run_signature")
}

.plant_panel_progress = function(context, stage, status, started, finished,
                                  message) {
  row = data.frame(
    stage = stage, status = status,
    started_at = if (is.na(started[[1]])) NA_character_ else
      format(started, "%Y-%m-%dT%H:%M:%S%z"),
    finished_at = if (length(finished) < 1 || is.na(finished[[1]])) {
      NA_character_
    } else format(finished, "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_seconds = if (length(finished) < 1 || is.na(finished[[1]])) {
      NA_real_
    } else round(as.numeric(difftime(finished, started, units = "secs")), 3),
    message = .plant_redact_secrets(message),
    run_signature = context$run_signature,
    stringsAsFactors = FALSE
  )
  old = .plant_panel_read_csv(context$paths$progress)
  old = if (is.data.frame(old) && nrow(old) > 0) {
    .plant_bind_tables(list(old), .plant_panel_run_manifest_cols())
  } else {
    .uaf_empty_table(.plant_panel_run_manifest_cols())
  }
  .plant_atomic_write_csv(rbind(old, row), context$paths$progress)
  invisible(row)
}

.plant_panel_write_status = function(context, state, stage, message) {
  status = list(
    workflow = "runPlantChemistryPanel", schema_version = "1.0.0",
    run_signature = context$run_signature,
    package_version = context$config$package_version,
    project_id = context$config$project_id,
    mode = context$config$mode, state = state, stage = stage,
    updated_at = .plant_timestamp(),
    plant_csv = context$config$plant_csv,
    output_dir = context$paths$root,
    cache_dir = context$config$cache_dir,
    message = .plant_redact_secrets(message)
  )
  .plant_atomic_write_json(status, context$paths$status)
  invisible(status)
}

.plant_panel_stage_failure = function(context, stage, condition, state) {
  message = .uaf_first_non_empty_text(
    .plant_redact_secrets(conditionMessage(condition)),
    ifelse(identical(state, "interrupted"),
           "Stage interrupted by user.", "Stage failed without a message.")
  )
  .plant_panel_progress(context, stage, state, Sys.time(), Sys.time(), message)
  .plant_panel_write_status(context, state, stage, message)
  .plant_atomic_write_json(
    list(stage = stage, state = state, message = message,
         failed_at = .plant_timestamp(), run_signature = context$run_signature),
    file.path(context$paths$root, "FAILED.json")
  )
  invisible(NULL)
}

.plant_panel_clear_transient_markers = function(context) {
  marker_files = file.path(
    context$paths$root,
    c("FAILED.json", "PAUSED_SERVICE_BUSY.json")
  )
  existing = marker_files[file.exists(marker_files)]
  if (length(existing) > 0) unlink(existing, force = TRUE)
  invisible(existing)
}

.plant_panel_stage_marker_file = function(context, stage) {
  file.path(context$paths$markers,
            paste0(gsub("-", "_", stage, fixed = TRUE),
                   "_COMPLETED.json"))
}

.plant_panel_write_stage_marker = function(context, stage, artifacts,
                                            started, message) {
  artifacts = unique(.uaf_non_empty(artifacts))
  artifacts = artifacts[file.exists(artifacts) & !dir.exists(artifacts)]
  manifest = .plant_panel_artifact_manifest(artifacts,
                                             context$paths$root)
  marker = list(
    workflow = "runPlantChemistryPanel", stage = stage,
    status = "completed", run_signature = context$run_signature,
    package_version = context$config$package_version,
    started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    completed_at = .plant_timestamp(),
    message = .plant_redact_secrets(message), artifacts = manifest
  )
  .plant_atomic_write_json(marker,
                           .plant_panel_stage_marker_file(context, stage))
  invisible(marker)
}

.plant_panel_stage_marker_valid = function(context, stage) {
  marker_file = .plant_panel_stage_marker_file(context, stage)
  marker = .plant_panel_read_json(marker_file)
  if (!is.list(marker) ||
      !identical(.uaf_first_non_empty_text(marker$status), "completed") ||
      !identical(.uaf_first_non_empty_text(marker$run_signature),
                 context$run_signature)) return(FALSE)
  manifest = marker$artifacts
  if (is.null(manifest)) return(FALSE)
  if (is.list(manifest) && !is.data.frame(manifest)) {
    manifest = tryCatch(as.data.frame(manifest, stringsAsFactors = FALSE),
                        error = function(error) data.frame())
  }
  if (!is.data.frame(manifest) || nrow(manifest) < 1) return(FALSE)
  artifact_paths = .uaf_squish_text(manifest$path)
  absolute = grepl("^/", artifact_paths) |
    grepl("^[A-Za-z]:[/\\\\]", artifact_paths)
  paths = artifact_paths
  paths[!absolute] = file.path(context$paths$root,
                               artifact_paths[!absolute])
  if (!all(file.exists(paths))) return(FALSE)
  observed = unname(tools::md5sum(paths))
  all(observed == manifest$md5)
}

.plant_panel_artifact_manifest = function(paths, root) {
  if (length(paths) < 1) {
    return(data.frame(path = character(), bytes = numeric(),
                      md5 = character(), sha256 = character(),
                      stringsAsFactors = FALSE))
  }
  normalized = normalizePath(paths, winslash = "/", mustWork = TRUE)
  root_norm = normalizePath(root, winslash = "/", mustWork = TRUE)
  relative = ifelse(startsWith(normalized, paste0(root_norm, "/")),
                    substring(normalized, nchar(root_norm) + 2L), normalized)
  data.frame(
    path = relative, bytes = as.numeric(file.info(normalized)$size),
    md5 = unname(tools::md5sum(normalized)),
    sha256 = vapply(normalized, .plant_sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
}

.plant_panel_snapshot_files = function(paths, destination, run_signature,
                                        overwrite, roles = NULL) {
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  if (is.null(roles)) {
    roles = c("plant_input", rep("supporting_ledger",
                                  max(0L, length(paths) - 1L)))
  }
  if (length(roles) != length(paths)) {
    stop("`roles` must contain one value per input snapshot.", call. = FALSE)
  }
  rows = lapply(seq_along(paths), function(i) {
    source = normalizePath(paths[[i]], winslash = "/", mustWork = TRUE)
    target = file.path(destination, basename(source))
    if (file.exists(target)) {
      same = identical(unname(tools::md5sum(source)[[1]]),
                       unname(tools::md5sum(target)[[1]]))
      if (!same && !isTRUE(overwrite)) {
        stop("Input snapshot differs from the current source: ", target,
             call. = FALSE)
      }
    }
    .plant_panel_copy_file(source, target, overwrite = TRUE)
    row_count = if (grepl("[.]csv$", target, ignore.case = TRUE)) {
      tryCatch(nrow(utils::read.csv(
        target, stringsAsFactors = FALSE, check.names = FALSE
      )), error = function(error) NA_integer_)
    } else if (grepl("[.](tsv|txt)$", target, ignore.case = TRUE)) {
      tryCatch(nrow(utils::read.delim(
        target, stringsAsFactors = FALSE, check.names = FALSE,
        quote = "", comment.char = ""
      )), error = function(error) NA_integer_)
    } else {
      NA_integer_
    }
    data.frame(
      role = as.character(roles[[i]]),
      source_path = source,
      snapshot_file = basename(target),
      bytes = as.numeric(file.info(target)$size),
      row_count = row_count,
      md5 = unname(tools::md5sum(target)[[1]]),
      sha256 = .plant_sha256_file(target),
      run_signature = run_signature, snapshotted_at = .plant_timestamp(),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.plant_panel_copy_file = function(source, target, overwrite = FALSE) {
  if (!file.exists(source) || dir.exists(source)) {
    stop("Source file does not exist: ", source, call. = FALSE)
  }
  if (file.exists(target) && !isTRUE(overwrite)) {
    stop("Target file exists: ", target, call. = FALSE)
  }
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  temp = tempfile(paste0(".", basename(target), "_"),
                  tmpdir = dirname(target))
  on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
  if (!file.copy(source, temp, overwrite = TRUE, copy.mode = TRUE,
                 copy.date = TRUE)) {
    stop("Could not copy source file: ", source, call. = FALSE)
  }
  .plant_atomic_replace(temp, target)
}

.plant_panel_read_input = function(path, species_col) {
  path = .plant_panel_existing_file(path, "plant input")
  x = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                      na.strings = c("", "NA"))
  if (nrow(x) < 1 || ncol(x) < 1) stop("Plant input contains no rows.",
                                       call. = FALSE)
  normalized = .plant_normalize_column_names(names(x))
  requested = .plant_normalize_column_names(species_col)
  index = match(requested, normalized)
  if (is.na(index)) stop("Plant input does not contain species column `",
                         species_col, "`.", call. = FALSE)
  if (names(x)[[index]] != "species") names(x)[[index]] = "species"
  x$species = .plant_canonical_taxon_name(x$species)
  x
}

.plant_panel_existing_file = function(path, name) {
  path = .uaf_first_non_empty_text(path)
  if (is.na(path) || !file.exists(path) || dir.exists(path)) {
    stop("`", name, "` must be an existing file.", call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.plant_panel_optional_existing_file = function(path, name) {
  path = .uaf_first_non_empty_text(path)
  if (is.na(path)) return(NA_character_)
  .plant_panel_existing_file(path, name)
}

.plant_panel_existing_files = function(paths) {
  paths = unique(.uaf_non_empty(paths))
  if (length(paths) < 1) return(character())
  missing = paths[!file.exists(paths) | dir.exists(paths)]
  if (length(missing) > 0) {
    stop("Supporting input file(s) missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  normalizePath(paths, winslash = "/", mustWork = TRUE)
}

.plant_panel_exclusion_count = function(paths) {
  paths = paths[file.exists(paths)]
  if (length(paths) < 1) return(0L)
  tables = lapply(paths, function(path) {
    x = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    names(x) = .plant_normalize_column_names(names(x))
    x
  })
  keys = unlist(lapply(tables, function(x) {
    order = if ("input_order" %in% names(x)) x$input_order else
      seq_len(nrow(x))
    name = if ("original_species" %in% names(x)) x$original_species else if (
      "species" %in% names(x)
    ) x$species else rep(NA_character_, nrow(x))
    paste(.uaf_squish_text(order), .plant_clean_name(name), sep = "\r")
  }), use.names = FALSE)
  length(unique(keys[!is.na(keys) & keys != ""]))
}

.plant_panel_path = function(path, name) {
  path = .uaf_first_non_empty_text(path)
  if (is.na(path)) stop("`", name, "` is required.", call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.plant_panel_required_file = function(path, message) {
  if (!file.exists(path) || dir.exists(path)) stop(message, call. = FALSE)
  path
}

.plant_panel_read_json = function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE),
           error = function(error) NULL)
}

.plant_panel_read_csv = function(path) {
  if (!file.exists(path)) return(data.frame())
  tryCatch(utils::read.csv(path, stringsAsFactors = FALSE,
                           check.names = FALSE),
           error = function(error) data.frame())
}

.plant_panel_optional_integer = function(x, name) {
  if (is.null(x) || length(x) < 1 || is.na(x[[1]])) return(NULL)
  .plant_panel_positive_integer(x, name)
}

.plant_panel_positive_integer = function(x, name) {
  value = suppressWarnings(as.integer(x[[1]]))
  if (is.na(value) || value < 1) stop("`", name,
                                      "` must be a positive integer.",
                                      call. = FALSE)
  value
}

.plant_panel_positive_or_inf = function(x, name) {
  value = suppressWarnings(as.numeric(x[[1]]))
  if (is.na(value) || value < 1) stop("`", name,
                                      "` must be positive or Inf.",
                                      call. = FALSE)
  value
}

.plant_panel_nonnegative = function(x, name) {
  value = suppressWarnings(as.numeric(x[[1]]))
  if (!is.finite(value) || value < 0) stop("`", name,
                                           "` must be non-negative.",
                                           call. = FALSE)
  value
}

.plant_panel_truthy = function(x) {
  tolower(.uaf_squish_text(x)) %in% c("true", "t", "yes", "y", "1")
}

.plant_panel_index_config = function(resource_dir) {
  path = file.path(resource_dir, "provider_indexes.json")
  x = .plant_panel_read_json(path)
  if (!is.list(x)) return(list(lotus = NA_character_, npass = NA_character_))
  list(lotus = .uaf_first_non_empty_text(x$lotus),
       npass = .uaf_first_non_empty_text(x$npass))
}

.plant_panel_index_manifest = function(index) {
  index = .uaf_first_non_empty_text(index)
  if (is.na(index)) return(NA_character_)
  if (dir.exists(index)) {
    candidates = c(file.path(index, "manifest.json"),
                   file.path(index, "lookup_manifest.json"))
    hit = candidates[file.exists(candidates)]
    if (length(hit) > 0) return(hit[[1]])
  }
  if (file.exists(index)) return(index)
  NA_character_
}

.plant_panel_require_resource_stage = function(context) {
  marker = .plant_panel_stage_marker_file(context, "build-indexes")
  indexes = .plant_panel_provider_indexes(context)
  availability = plantProviderAvailability(
    sources = context$config$sources, provider_indexes = indexes,
    lotus_index = indexes$lotus, probe_live = FALSE
  )
  local = availability$provider %in% c("lotus", "npass")
  if (!file.exists(marker) ||
      any(local & availability$availability_status != "available")) {
    stop("Run mode = 'build-indexes' before this stage.", call. = FALSE)
  }
  invisible(TRUE)
}

.plant_panel_npass_raw_paths = function(root) {
  root = .uaf_first_non_empty_text(root)
  if (is.na(root)) root = ""
  files = c(
    general_info = "NPASS3.0_naturalproducts_generalinfo.txt",
    structures = "NPASS3.0_naturalproducts_structure.txt",
    species_pairs = "NPASS3.0_naturalproducts_species_pair.txt",
    species_info = "NPASS3.0_species_info.txt"
  )
  stats::setNames(file.path(root, files), names(files))
}

.plant_panel_npass_urls = function() {
  base = "https://bidd.group/NPASS/downloadFiles/"
  paths = .plant_panel_npass_raw_paths(".")
  files = basename(paths)
  stats::setNames(paste0(base, files), names(paths))
}

.plant_panel_download_npass = function(root, progress = TRUE,
                                        request_timeout = 3600,
                                        max_attempts = 5L) {
  root = .uaf_first_non_empty_text(root)
  if (is.na(root)) stop("`npass_raw_dir` is required for NPASS download.",
                        call. = FALSE)
  root = normalizePath(root, winslash = "/", mustWork = FALSE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  paths = .plant_panel_npass_raw_paths(root)
  urls = .plant_panel_npass_urls()
  for (name in names(paths)) {
    if (isTRUE(progress)) message("uafR NPASS download: ", name)
    .plant_panel_download_resource(
      urls[[name]], paths[[name]], request_timeout = request_timeout,
      max_attempts = max_attempts, progress = progress
    )
  }
  paths
}

.plant_panel_download_resource = function(
    url, path, request_timeout = 3600, max_attempts = 5L,
    retry_wait = 5, progress = TRUE,
    metadata_fun = .plant_panel_remote_metadata,
    transfer_fun = .plant_panel_transfer_resource,
    sleep_fun = Sys.sleep) {
  request_timeout = suppressWarnings(as.numeric(request_timeout[[1]]))
  if (!is.finite(request_timeout) || request_timeout < 1) {
    stop("`request_timeout` must be a positive number.", call. = FALSE)
  }
  max_attempts = suppressWarnings(as.integer(max_attempts[[1]]))
  if (is.na(max_attempts) || max_attempts < 1L) {
    stop("`max_attempts` must be a positive integer.", call. = FALSE)
  }
  retry_wait = suppressWarnings(as.numeric(retry_wait[[1]]))
  if (!is.finite(retry_wait) || retry_wait < 0) {
    stop("`retry_wait` must be non-negative.", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial = paste0(path, ".partial")
  metadata = tryCatch(
    metadata_fun(url, request_timeout),
    error = function(error) list(expected_bytes = NA_real_,
                                 error = conditionMessage(error))
  )
  expected = suppressWarnings(as.numeric(metadata$expected_bytes[[1]]))
  if (!is.finite(expected) || expected < 1) expected = NA_real_
  if (file.exists(path) && !dir.exists(path)) {
    existing = as.numeric(file.info(path)$size)
    if ((!is.na(expected) && existing == expected) ||
        (is.na(expected) && existing >= 100)) {
      return(invisible(path))
    }
    if (!is.na(expected) && existing > 0 && existing < expected) {
      partial_size = if (file.exists(partial)) {
        as.numeric(file.info(partial)$size)
      } else 0
      if (existing > partial_size) {
        if (file.exists(partial)) unlink(partial, force = TRUE)
        if (!file.rename(path, partial)) {
          if (!file.copy(path, partial, overwrite = TRUE)) {
            stop("Could not stage incomplete resource for resume: ", path,
                 call. = FALSE)
          }
          unlink(path, force = TRUE)
        }
      } else {
        unlink(path, force = TRUE)
      }
    } else {
      unlink(path, force = TRUE)
    }
  }
  if (file.exists(partial) && !is.na(expected) &&
      file.info(partial)$size > expected) {
    unlink(partial, force = TRUE)
  }
  last_error = NULL
  for (attempt in seq_len(max_attempts)) {
    offset = if (file.exists(partial)) as.numeric(file.info(partial)$size) else
      0
    if (!is.na(expected) && identical(offset, expected)) {
      .plant_atomic_replace(partial, path)
      return(invisible(path))
    }
    last_error = tryCatch({
      transfer_fun(url, partial, offset, request_timeout, expected)
      NULL
    }, error = function(error) error)
    observed = if (file.exists(partial)) {
      as.numeric(file.info(partial)$size)
    } else 0
    complete = (!is.na(expected) && observed > 0 && observed == expected) ||
      (is.na(expected) && observed >= 100 && is.null(last_error))
    if (isTRUE(complete)) {
      .plant_atomic_replace(partial, path)
      return(invisible(path))
    }
    if (!is.na(expected) && observed > expected) {
      unlink(partial, force = TRUE)
      observed = 0
    }
    if (attempt < max_attempts) {
      wait = min(60, retry_wait * 2^(attempt - 1L))
      if (isTRUE(progress)) {
        message(
          "uafR resource download incomplete (attempt ", attempt, "/",
          max_attempts, "; ", format(observed, big.mark = ","),
          if (!is.na(expected)) paste0("/", format(expected, big.mark = ","))
          else "", " bytes). Resuming after ", wait, " second(s)."
        )
      }
      sleep_fun(wait)
    }
  }
  stop(
    "Resource download failed after ", max_attempts, " attempt(s): ", url,
    if (!is.null(last_error)) paste0("; ", conditionMessage(last_error)) else
      "; byte count did not match the remote resource",
    call. = FALSE
  )
}

.plant_panel_remote_metadata = function(url, request_timeout) {
  if (!requireNamespace("curl", quietly = TRUE)) {
    return(list(expected_bytes = NA_real_, etag = NA_character_))
  }
  handle = curl::new_handle(
    nobody = TRUE, followlocation = TRUE, failonerror = TRUE,
    timeout = request_timeout,
    connecttimeout = min(60, request_timeout)
  )
  curl::handle_setheaders(
    handle,
    .list = list(
      "Accept-Encoding" = "identity",
      "User-Agent" = "uafR NPASS resource downloader"
    )
  )
  response = curl::curl_fetch_memory(url, handle = handle)
  headers = curl::parse_headers_list(response$headers)
  expected = suppressWarnings(as.numeric(headers[["content-length"]]))
  list(
    expected_bytes = if (length(expected) && is.finite(expected)) expected else
      NA_real_,
    etag = .uaf_first_non_empty_text(headers[["etag"]]),
    last_modified = .uaf_first_non_empty_text(headers[["last-modified"]]),
    accept_ranges = .uaf_first_non_empty_text(headers[["accept-ranges"]])
  )
}

.plant_panel_transfer_resource = function(url, partial, offset,
                                           request_timeout, expected) {
  if (!requireNamespace("curl", quietly = TRUE)) {
    old_timeout = getOption("timeout")
    options(timeout = max(request_timeout, old_timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)
    if (file.exists(partial)) unlink(partial, force = TRUE)
    utils::download.file(url, partial, mode = "wb", quiet = TRUE,
                         method = "libcurl")
    return(invisible(partial))
  }
  handle = curl::new_handle(
    followlocation = TRUE, failonerror = TRUE,
    timeout = request_timeout,
    connecttimeout = min(60, request_timeout)
  )
  curl::handle_setheaders(
    handle,
    .list = list(
      "Accept-Encoding" = "identity",
      "User-Agent" = "uafR NPASS resource downloader"
    )
  )
  if (offset > 0) {
    curl::handle_setopt(handle, resume_from_large = offset)
  }
  input = NULL
  output = NULL
  on.exit({
    if (!is.null(input)) try(close(input), silent = TRUE)
    if (!is.null(output)) try(close(output), silent = TRUE)
  }, add = TRUE)
  input = curl::curl(url, open = "rb", handle = handle)
  output = file(partial, open = if (offset > 0) "ab" else "wb")
  repeat {
    block = readBin(input, what = "raw", n = 1024L * 1024L)
    if (length(block) < 1L) break
    writeBin(block, output)
  }
  close(input)
  input = NULL
  close(output)
  output = NULL
  invisible(partial)
}

.plant_panel_result_fingerprint = function(result) {
  tables = c("PlantQueries", "PlantQueryAliases", "ProviderDiagnostics",
             "ProviderQueryAccounting", "PlantCompoundOccurrences",
             "LiteratureCandidates", "SourceCompoundIdentity",
             "CompoundResolution")
  path = tempfile("uafR_panel_fingerprint_", fileext = ".rds")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  saveRDS(result[intersect(tables, names(result))], path, version = 3)
  unname(tools::md5sum(path)[[1]])
}

.plant_panel_tanimoto_runner_path = function() {
  installed = system.file("scripts", "run_plant_tanimoto_server.R",
                          package = "uafR")
  if (nzchar(installed) && file.exists(installed)) return(installed)
  candidates = c(
    file.path(getwd(), "inst", "scripts", "run_plant_tanimoto_server.R"),
    file.path(getwd(), "tools", "run_plant_tanimoto_server.R")
  )
  hit = candidates[file.exists(candidates)]
  if (length(hit) < 1) {
    stop("The plant Tanimoto server runner is not installed.", call. = FALSE)
  }
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

.plant_panel_source_snapshot = function(context, runner_copy,
                                         species_universe_file = NULL,
                                         release_artifacts = NULL,
                                         support_artifacts = NULL) {
  files = c(species_universe_file, runner_copy, support_artifacts,
            release_artifacts)
  files = files[file.exists(files)]
  role_by_name = c(
    stats::setNames("plant_input", basename(context$config$plant_csv)),
    stats::setNames("species_universe", basename(species_universe_file)),
    stats::setNames("server_runner", basename(runner_copy)),
    run_contract.json = "run_contract",
    input_manifest.csv = "input_manifest",
    stats::setNames(
      ifelse(grepl("[.]tar[.]gz$", basename(release_artifacts),
                   ignore.case = TRUE), "source_tarball", "release_manifest"),
      basename(release_artifacts)
    )
  )
  data.frame(
    role = unname(role_by_name[basename(files)]),
    file = basename(files), source_path = basename(files),
    bytes = as.numeric(file.info(files)$size),
    md5 = unname(tools::md5sum(files)),
    sha256 = vapply(files, .plant_sha256_file, character(1)),
    package_version = context$config$package_version,
    run_signature = context$run_signature,
    created_at = .plant_timestamp(), stringsAsFactors = FALSE
  )
}

.plant_panel_file_references = function(context, server_root) {
  paths = c(
    release_manifest = context$config$release_manifest,
    source_tarball = context$config$source_tarball,
    input_manifest = file.path(context$paths$inputs, "input_manifest.csv"),
    run_contract = file.path(context$paths$inputs, "run_contract.json"),
    provider_resources = file.path(context$paths$resources,
                                   "ProviderResourceManifest.csv"),
    discovery_gate = file.path(context$paths$stages, "discovery",
                               "discovery_quality_gate.csv"),
    server_manifest = file.path(server_root,
                                "server_tanimoto_output_manifest.csv")
  )
  paths[file.exists(paths)]
}

.plant_panel_add_discovery_bundle_tables = function(path, result) {
  tables = list(
    "30_PlantQueries.csv" = result$PlantQueries,
    "31_PlantQueryAliases.csv" = result$PlantQueryAliases,
    "32_PlantNameResolution.csv" = result$PlantNameResolution,
    "33_ProviderDiagnostics.csv" = result$ProviderDiagnostics,
    "34_ProviderQueryAccounting.csv" = result$ProviderQueryAccounting,
    "35_ProviderResourceManifest.csv" = result$ProviderResourceManifest,
    "36_AllPlantCompoundOccurrences.csv" = result$PlantCompoundOccurrences,
    "37_PlantContextEvidence.csv" = result$PlantContextEvidence,
    "38_ProviderContextAudit.csv" = result$ProviderContextAudit,
    "39_LiteratureCandidates.csv" = result$LiteratureCandidates,
    "40_SourceCompoundIdentity.csv" = result$SourceCompoundIdentity,
    "41_CompoundIdentityReview.csv" = result$CompoundIdentityReview
  )
  for (name in names(tables)) {
    if (is.data.frame(tables[[name]])) {
      .plant_atomic_write_csv(tables[[name]], file.path(path, name))
    }
  }
  .bundle_refresh_manifest(path)
  invisible(path)
}

.plant_panel_validation_status = function(validation) {
  if (!is.list(validation) || !is.data.frame(validation$Summary) ||
      nrow(validation$Summary) < 1) return("unknown")
  fields = c("ExportReadyStatus", "Status", "status")
  field = fields[fields %in% names(validation$Summary)]
  if (length(field) < 1) return("unknown")
  .uaf_first_non_empty_text(validation$Summary[[field[[1]]]])
}

.plant_panel_bundle_ready = function(validation) {
  status = tolower(.plant_panel_validation_status(validation))
  status %in% c("ready", "export_ready", "pass", "passed", "valid") ||
    grepl("ready", status, fixed = TRUE)
}
