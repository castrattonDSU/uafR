#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

usage = function() {
  cat(paste(
    "Run an identity-verified plant Tanimoto handoff on a compute server.",
    "",
    "Usage:",
    "  Rscript tools/run_plant_tanimoto_server.R \\",
    "    --input plant_compound_membership_tanimoto_ready.csv \\",
    "    --out-dir plant_tanimoto_server_output \\",
    "    --cache-dir pubchem_tanimoto_cache \\",
    "    [--manifest tanimoto_input_export_manifest.csv] \\",
    "    [--species-universe plant_species_universe.csv] \\",
    "    [--release-manifest uafR_release_manifest.json] \\",
    "    [--mode preflight|smoke|summary|full] \\",
    "    [--smoke-structures 25] [--progress-every 25] \\",
    "    [--service-busy-limit 2] [--pair-shard-rows 5000000] \\",
    "    [--full-pairs true] [--throttle 1.1] \\",
    "    [--pair-block-size 250] [--overwrite false] \\",
    "    [--overwrite-incomplete false]",
    "",
    "Modes:",
    "  preflight  Validate inputs, dependencies, estimates, and disk only.",
    "  smoke      Resolve a deterministic plant-balanced structure subset.",
    "  summary    Resolve all fingerprints and write plant-pair summaries.",
    "  full       Write compound pairs and optional cross-plant pair rows.",
    "",
    "Live modes require UAFR_CONFIRM_SERVER_TANIMOTO=YES. If --mode is",
    "omitted, the historical behavior is preserved: preflight exits with",
    "status 2 unless that confirmation is set, then the full mode runs.",
    sep = "\n"
  ))
}

parse_args = function(args) {
  out = list()
  i = 1L
  while (i <= length(args)) {
    token = args[[i]]
    if (token %in% c("-h", "--help")) {
      usage()
      quit(save = "no", status = 0L)
    }
    if (!startsWith(token, "--")) {
      stop("Unexpected argument: ", token, call. = FALSE)
    }
    key = gsub("-", "_", substring(token, 3), fixed = TRUE)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] = "true"
      i = i + 1L
    } else {
      out[[key]] = args[[i + 1L]]
      i = i + 2L
    }
  }
  out
}

required_arg = function(args, name) {
  value = args[[name]]
  if (is.null(value) || is.na(value) || !nzchar(value)) {
    stop("Missing required --", gsub("_", "-", name, fixed = TRUE),
         " argument.", call. = FALSE)
  }
  value
}

as_flag = function(x, default = FALSE) {
  if (is.null(x)) return(default)
  value = tolower(trimws(as.character(x[[1]])))
  if (!value %in% c("true", "false", "yes", "no", "1", "0")) {
    stop("Logical arguments must be true or false; received: ", value,
         call. = FALSE)
  }
  value %in% c("true", "yes", "1")
}

as_positive_number = function(x, default, name, integer = FALSE) {
  value = suppressWarnings(as.numeric(if (is.null(x)) default else x))
  if (length(value) != 1L || is.na(value) || value < 1) {
    stop("`--", gsub("_", "-", name, fixed = TRUE),
         "` must be a positive number.", call. = FALSE)
  }
  if (integer) as.integer(value) else value
}

atomic_csv = function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp = tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temp, force = TRUE), add = TRUE)
  utils::write.csv(x, temp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(temp, path)) stop("Could not write: ", path, call. = FALSE)
  invisible(path)
}

atomic_json = function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp = tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temp, force = TRUE), add = TRUE)
  jsonlite::write_json(x, temp, pretty = TRUE, auto_unbox = TRUE,
                       na = "null", dataframe = "rows")
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(temp, path)) stop("Could not write: ", path, call. = FALSE)
  invisible(path)
}

atomic_text = function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp = tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temp, force = TRUE), add = TRUE)
  writeLines(enc2utf8(as.character(x)), temp, useBytes = TRUE)
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(temp, path)) stop("Could not write: ", path, call. = FALSE)
  invisible(path)
}

clean_text = function(x) {
  x = trimws(gsub("[[:space:]]+", " ", as.character(x), perl = TRUE))
  x[x == "" | toupper(x) %in% c("NA", "N/A", "NULL")] = NA_character_
  x
}

normalize_cid = function(x) {
  x = clean_text(x)
  numeric = suppressWarnings(as.numeric(x))
  valid = !is.na(x) & grepl("^[0-9]+$", x) & is.finite(numeric) &
    numeric > 0 & numeric == floor(numeric)
  out = rep(NA_character_, length(x))
  out[valid] = format(numeric[valid], scientific = FALSE, trim = TRUE)
  out
}

read_species_universe = function(path, membership_species) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    return(sort(unique(clean_text(membership_species))))
  }
  if (!file.exists(path) || dir.exists(path)) {
    stop("Species-universe file does not exist: ", path, call. = FALSE)
  }
  x = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                      na.strings = c("", "NA", "N/A", "NULL"),
                      fileEncoding = "UTF-8")
  normalized = tolower(gsub("[^a-z0-9]+", "_", names(x)))
  candidates = match(c("species", "accepted_species_name", "plant"),
                     normalized, nomatch = 0L)
  candidates = candidates[candidates > 0L]
  if (length(candidates) < 1L) {
    stop("Species-universe CSV must contain a species column.", call. = FALSE)
  }
  species = clean_text(x[[candidates[[1L]]]])
  species = species[!is.na(species)]
  if (length(species) < 1L) {
    stop("Species-universe CSV contains no usable species names.",
         call. = FALSE)
  }
  if (anyDuplicated(species) > 0L) {
    stop("Species-universe CSV contains duplicate species names.",
         call. = FALSE)
  }
  species
}

complete_plant_pair_summary = function(summary, species_universe,
                                        membership) {
  species_universe = unique(clean_text(species_universe))
  species_universe = species_universe[!is.na(species_universe)]
  if (length(species_universe) < 2L) {
    stop("At least two species are required for plant-pair output.",
         call. = FALSE)
  }
  pairs = utils::combn(sort(species_universe), 2L)
  out = data.frame(
    group_pair_id = paste(pairs[1L, ], pairs[2L, ], sep = "__"),
    group_a = pairs[1L, ],
    group_b = pairs[2L, ],
    group_level = "species",
    compounds_a = 0L,
    compounds_b = 0L,
    compound_pair_count = 0,
    shared_compound_count = 0L,
    mean_tanimoto = NA_real_,
    median_tanimoto = NA_real_,
    p95_tanimoto = NA_real_,
    max_tanimoto = NA_real_,
    top_compound_pair_ids = NA_character_,
    top_compound_pairs = NA_character_,
    fingerprint_source = "PubChem_Fingerprint2D",
    stringsAsFactors = FALSE
  )
  threshold_cols = grep("^compound_pair_count_ge_", names(summary),
                        value = TRUE)
  for (col in threshold_cols) out[[col]] = 0L
  compounds_by_species = split(membership$compound_id, membership$group_id)
  compounds_by_species = lapply(compounds_by_species, unique)
  compound_set = function(species) {
    value = compounds_by_species[[species]]
    if (is.null(value)) character() else value
  }
  count_for = function(species) length(compound_set(species))
  out$compounds_a = vapply(out$group_a, count_for, integer(1))
  out$compounds_b = vapply(out$group_b, count_for, integer(1))
  out$compound_pair_count = out$compounds_a * out$compounds_b
  out$shared_compound_count = vapply(seq_len(nrow(out)), function(i) {
    length(intersect(compound_set(out$group_a[[i]]),
                     compound_set(out$group_b[[i]])))
  }, integer(1))

  summary = as.data.frame(summary, stringsAsFactors = FALSE)
  if (nrow(summary) > 0L) {
    summary_key = paste(pmin(summary$group_a, summary$group_b),
                        pmax(summary$group_a, summary$group_b), sep = "\r")
    output_key = paste(out$group_a, out$group_b, sep = "\r")
    idx = match(output_key, summary_key)
    matched = !is.na(idx)
    copy_cols = intersect(names(summary), names(out))
    for (col in copy_cols) out[[col]][matched] = summary[[col]][idx[matched]]
  } else {
    matched = rep(FALSE, nrow(out))
  }
  out$support_status = ifelse(matched, "computed", "insufficient_support")
  out$support_note = ifelse(
    matched,
    "Tanimoto metrics were computed from at least one cross-plant structure pair.",
    "One or both plants had no usable resolved structure, so Tanimoto metrics are unavailable."
  )
  out
}

valid_inchikey = function(x) {
  x = toupper(clean_text(x))
  !is.na(x) & grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", x)
}

test_writable = function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  probe = tempfile("uafr_write_probe_", tmpdir = path)
  ok = tryCatch({
    writeLines("uafR", probe, useBytes = TRUE)
    file.exists(probe)
  }, error = function(error) FALSE)
  unlink(probe, force = TRUE)
  isTRUE(ok)
}

disk_free_bytes = function(path) {
  command = Sys.which("df")
  if (!nzchar(command)) return(NA_real_)
  output = tryCatch(
    system2(command, c("-Pk", shQuote(normalizePath(path, mustWork = TRUE))),
            stdout = TRUE, stderr = FALSE),
    error = function(error) character()
  )
  if (length(output) < 2L) return(NA_real_)
  fields = strsplit(trimws(tail(output, 1L)), "[[:space:]]+")[[1]]
  if (length(fields) < 4L) return(NA_real_)
  available_kb = suppressWarnings(as.numeric(fields[[4L]]))
  if (!is.finite(available_kb)) return(NA_real_)
  available_kb * 1024
}

verify_input_manifest = function(input, manifest_file) {
  if (is.null(manifest_file) || is.na(manifest_file) || !nzchar(manifest_file)) {
    return(data.frame(
      check = "input_manifest_available", status = "warning",
      detail = "No handoff manifest was supplied; input checksum was not compared.",
      stringsAsFactors = FALSE
    ))
  }
  if (!file.exists(manifest_file)) {
    return(data.frame(
      check = "input_manifest_available", status = "fail",
      detail = paste("Manifest does not exist:", manifest_file),
      stringsAsFactors = FALSE
    ))
  }
  manifest = utils::read.csv(manifest_file, stringsAsFactors = FALSE,
                             check.names = FALSE, na.strings = c("", "NA"))
  file_col = intersect(c("file", "path"), names(manifest))
  if (length(file_col) < 1L || !"md5" %in% names(manifest)) {
    return(data.frame(
      check = "input_manifest_schema", status = "fail",
      detail = "Manifest must contain file/path and md5 columns.",
      stringsAsFactors = FALSE
    ))
  }
  listed = basename(clean_text(manifest[[file_col[[1L]]]]))
  idx = which(listed == basename(input))
  if (length(idx) != 1L) {
    return(data.frame(
      check = "input_manifest_membership", status = "fail",
      detail = paste("Expected exactly one manifest row for", basename(input)),
      stringsAsFactors = FALSE
    ))
  }
  observed = unname(tools::md5sum(input)[[1L]])
  expected = clean_text(manifest$md5[[idx]])
  data.frame(
    check = "input_checksum_matches_manifest",
    status = ifelse(identical(observed, expected), "pass", "fail"),
    detail = paste("expected", expected, "observed", observed),
    stringsAsFactors = FALSE
  )
}

release_manifest_check = function(path, package_version) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    return(data.frame(
      check = "release_manifest_available", status = "warning",
      detail = paste("No release manifest supplied; installed uafR version is",
                     package_version), stringsAsFactors = FALSE
    ))
  }
  if (!file.exists(path)) {
    return(data.frame(
      check = "release_manifest_available", status = "fail",
      detail = paste("Release manifest does not exist:", path),
      stringsAsFactors = FALSE
    ))
  }
  manifest = tryCatch({
    if (grepl("[.]json$", path, ignore.case = TRUE)) {
      jsonlite::read_json(path, simplifyVector = TRUE)
    } else {
      x = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
      if (nrow(x) < 1L) stop("Release manifest contains no rows.")
      as.list(x[1L, , drop = FALSE])
    }
  }, error = function(error) error)
  if (inherits(manifest, "error")) {
    return(data.frame(
      check = "release_manifest_readable", status = "fail",
      detail = conditionMessage(manifest), stringsAsFactors = FALSE
    ))
  }
  version_fields = c("package_version", "uafR_package_version", "version")
  field = version_fields[version_fields %in% names(manifest)]
  if (length(field) < 1L) {
    return(data.frame(
      check = "release_manifest_package_version", status = "fail",
      detail = "Release manifest has no package_version field.",
      stringsAsFactors = FALSE
    ))
  }
  expected = clean_text(manifest[[field[[1L]]]])[[1L]]
  data.frame(
    check = "installed_package_matches_release_manifest",
    status = ifelse(identical(package_version, expected), "pass", "fail"),
    detail = paste("manifest", expected, "installed", package_version),
    stringsAsFactors = FALSE
  )
}

progress_recorder = function(path, run_id, mode) {
  metrics = new.env(parent = emptyenv())
  metrics$request_count = 0L
  metrics$cache_hit_count = 0L
  metrics$error_count = 0L
  metrics$rate_limit_count = 0L
  callback = function(event) {
    event_name = if (is.null(event$event)) "progress" else event$event
    if (identical(event_name, "request_attempt")) {
      metrics$request_count = metrics$request_count + 1L
    }
    if (identical(event_name, "cache_hit")) {
      metrics$cache_hit_count = metrics$cache_hit_count + 1L
    }
    status_code = suppressWarnings(as.integer(event$status_code))
    if (length(status_code) != 1L || is.na(status_code)) {
      status_code = NA_integer_
    }
    if (!is.na(status_code) && status_code %in% c(429L, 503L)) {
      metrics$rate_limit_count = metrics$rate_limit_count + 1L
    }
    if (event_name %in% c("request_exhausted", "service_busy_limit_reached")) {
      metrics$error_count = metrics$error_count + 1L
    }
    row = data.frame(
      run_id = run_id,
      mode = mode,
      event_at = if (is.null(event$event_at)) {
        format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      } else event$event_at,
      phase = if (is.null(event$phase)) NA_character_ else event$phase,
      event = event_name,
      completed = if (is.null(event$completed)) NA_real_ else event$completed,
      total = if (is.null(event$total)) NA_real_ else event$total,
      status_code = status_code,
      cache_hit = if (is.null(event$cache_hit)) FALSE else event$cache_hit,
      detail = if (is.null(event$detail)) NA_character_ else event$detail,
      url = if (is.null(event$url)) NA_character_ else
        substr(event$url, 1L, 500L),
      stringsAsFactors = FALSE
    )
    write_header = !file.exists(path)
    utils::write.table(row, path, sep = ",", row.names = FALSE,
                       col.names = write_header, append = !write_header,
                       quote = TRUE, na = "")
    invisible(TRUE)
  }
  list(callback = callback, metrics = metrics)
}

write_status = function(path, run_id, mode, state, input, out_dir,
                        package_version, message = NA_character_,
                        counts = list()) {
  atomic_json(c(list(
    workflow = "run_plant_tanimoto_server",
    schema_version = "1.0.0",
    run_id = run_id,
    mode = mode,
    state = state,
    updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    package_version = package_version,
    input_file = input,
    output_directory = normalizePath(out_dir, mustWork = FALSE),
    message = message
  ), counts), path)
}

archive_staging = function(staging_dir, incomplete_root, run_id, suffix) {
  if (!dir.exists(staging_dir)) return(NA_character_)
  archive_dir = file.path(incomplete_root, "archive",
                          paste0(basename(staging_dir), "_", run_id, "_", suffix))
  dir.create(dirname(archive_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(staging_dir, archive_dir)) return(staging_dir)
  archive_dir
}

publish_staging = function(staging_dir, final_dir, legacy_layout, overwrite) {
  if (!legacy_layout) {
    if (dir.exists(final_dir)) {
      if (!overwrite) stop("Completed mode output already exists: ", final_dir,
                           call. = FALSE)
      unlink(final_dir, recursive = TRUE, force = TRUE)
    }
    if (!file.rename(staging_dir, final_dir)) {
      stop("Could not atomically publish completed mode directory: ",
           final_dir, call. = FALSE)
    }
    return(invisible(final_dir))
  }
  files = list.files(staging_dir, full.names = TRUE, recursive = TRUE,
                     include.dirs = FALSE)
  for (source in files) {
    relative = substring(source, nchar(staging_dir) + 2L)
    target = file.path(final_dir, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(target)) {
      if (!overwrite) stop("Output already exists: ", target, call. = FALSE)
      unlink(target, force = TRUE)
    }
    if (!file.rename(source, target)) {
      stop("Could not publish completed output: ", target, call. = FALSE)
    }
  }
  unlink(staging_dir, recursive = TRUE, force = TRUE)
  invisible(final_dir)
}

build_output_manifest = function(final_dir, out_dir, mode, run_id,
                                  result, table_rows) {
  files = list.files(final_dir, full.names = TRUE, recursive = TRUE,
                     include.dirs = FALSE)
  files = files[!grepl("[.]partial$", files)]
  relative = substring(files, nchar(normalizePath(out_dir, mustWork = TRUE)) + 2L)
  rows_by_path = numeric()
  if (is.data.frame(result$ExportManifest) && nrow(result$ExportManifest) > 0L) {
    rows_by_path = stats::setNames(
      as.numeric(result$ExportManifest$Rows),
      basename(result$ExportManifest$Path)
    )
  }
  row_count = vapply(files, function(path) {
    name = basename(path)
    if (name %in% names(table_rows)) return(as.numeric(table_rows[[name]]))
    if (name %in% names(rows_by_path)) return(as.numeric(rows_by_path[[name]]))
    NA_real_
  }, numeric(1))
  data.frame(
    mode = mode,
    run_id = run_id,
    file = relative,
    row_count = row_count,
    bytes = unname(file.info(files)$size),
    md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )[order(relative), , drop = FALSE]
}

main = function(command_args = commandArgs(trailingOnly = TRUE),
                request_fun = NULL) {
  args = parse_args(command_args)
  input = normalizePath(required_arg(args, "input"), mustWork = TRUE)
  out_dir = required_arg(args, "out_dir")
  cache_dir = required_arg(args, "cache_dir")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir = normalizePath(out_dir, mustWork = TRUE)
  cache_dir = normalizePath(cache_dir, mustWork = TRUE)

  mode_supplied = !is.null(args$mode)
  requested_mode = if (mode_supplied) tolower(args$mode) else "legacy"
  if (!requested_mode %in% c("legacy", "preflight", "smoke", "summary", "full")) {
    stop("`--mode` must be preflight, smoke, summary, or full.",
         call. = FALSE)
  }
  confirmed = identical(Sys.getenv("UAFR_CONFIRM_SERVER_TANIMOTO", ""), "YES")
  mode = if (identical(requested_mode, "legacy")) {
    if (confirmed) "full" else "preflight"
  } else requested_mode
  legacy_layout = !mode_supplied && identical(mode, "full")

  manifest_file = args$manifest
  if (is.null(manifest_file)) {
    candidate = file.path(dirname(input), "tanimoto_input_export_manifest.csv")
    if (file.exists(candidate)) manifest_file = candidate
  }
  if (!is.null(manifest_file)) {
    manifest_file = normalizePath(manifest_file, mustWork = FALSE)
  }
  release_manifest = args$release_manifest
  if (!is.null(release_manifest)) {
    release_manifest = normalizePath(release_manifest, mustWork = FALSE)
  }

  full_pairs = as_flag(args$full_pairs, default = FALSE)
  overwrite = as_flag(args$overwrite, default = FALSE)
  overwrite_incomplete = as_flag(args$overwrite_incomplete, default = FALSE)
  extreme_output_confirm = as_flag(args$extreme_output_confirm, default = FALSE)
  throttle = suppressWarnings(as.numeric(if (is.null(args$throttle)) 1.1 else
    args$throttle))
  if (!is.finite(throttle) || throttle < 1) {
    stop("`--throttle` must be at least 1 second for a production server run.",
         call. = FALSE)
  }
  pair_block_size = as_positive_number(args$pair_block_size, 250L,
                                        "pair_block_size", integer = TRUE)
  pair_shard_rows = as_positive_number(args$pair_shard_rows, 5000000,
                                        "pair_shard_rows")
  smoke_structures = as_positive_number(args$smoke_structures, 25L,
                                         "smoke_structures", integer = TRUE)
  progress_every = as_positive_number(args$progress_every, 25L,
                                       "progress_every", integer = TRUE)
  service_busy_limit = as_positive_number(args$service_busy_limit, 2L,
                                           "service_busy_limit", integer = TRUE)
  min_fingerprint_usable = suppressWarnings(as.numeric(
    if (is.null(args$min_fingerprint_usable)) 0.8 else args$min_fingerprint_usable
  ))
  if (!is.finite(min_fingerprint_usable) || min_fingerprint_usable < 0 ||
      min_fingerprint_usable > 1) {
    stop("`--min-fingerprint-usable` must be between 0 and 1.",
         call. = FALSE)
  }

  if (!requireNamespace("uafR", quietly = TRUE)) {
    stop("Install the same validated uafR source version on the server first.",
         call. = FALSE)
  }
  package_version = as.character(utils::packageVersion("uafR"))
  run_id = paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid())
  status_path = file.path(out_dir, "server_run_status.json")
  progress_path = file.path(out_dir, "server_progress.csv")

  membership = utils::read.csv(
    input, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = c("", "NA", "N/A", "NULL"), fileEncoding = "UTF-8"
  )
  required = c("species", "compound_id", "compound_name", "InChIKey", "CID")
  missing = setdiff(required, names(membership))
  if (length(missing) > 0L) {
    stop("Input is missing required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  membership$species = clean_text(membership$species)
  membership$compound_id = clean_text(membership$compound_id)
  membership$InChIKey = toupper(clean_text(membership$InChIKey))
  membership$CID = normalize_cid(membership$CID)
  key_ok = valid_inchikey(membership$InChIKey)
  cid_ok = !is.na(membership$CID)
  membership_key = paste(membership$species, membership$compound_id, sep = "\r")
  identity_key = ifelse(key_ok, membership$InChIKey,
                        ifelse(cid_ok, paste0("cid:", membership$CID), NA))
  identity_count = ave(identity_key, membership$compound_id,
                       FUN = function(x) length(unique(x[!is.na(x)])))
  species_universe_file = args$species_universe
  if (is.null(species_universe_file)) {
    candidate = file.path(dirname(input), "plant_species_universe.csv")
    if (file.exists(candidate)) species_universe_file = candidate
  }
  if (!is.null(species_universe_file)) {
    species_universe_file = normalizePath(species_universe_file,
                                          mustWork = FALSE)
  }
  species_universe = read_species_universe(species_universe_file,
                                           membership$species)
  membership_outside_universe = setdiff(unique(membership$species),
                                        species_universe)
  review_count = if ("identity_review_required" %in% names(membership)) {
    sum(tolower(clean_text(membership$identity_review_required)) %in%
          c("true", "yes", "1", "required"), na.rm = TRUE)
  } else 0L

  species_counts = table(membership$species)
  unique_compounds = length(unique(membership$compound_id))
  compound_pair_count = choose(unique_compounds, 2)
  plant_pair_count = choose(length(species_universe), 2)
  cross_plant_pair_count = choose(nrow(membership), 2) -
    sum(vapply(as.numeric(species_counts), function(n) choose(n, 2), numeric(1)))
  requested_pair_rows = if (identical(mode, "full")) {
    compound_pair_count + if (full_pairs) cross_plant_pair_count else 0
  } else 0
  planned_pair_rows = compound_pair_count + if (full_pairs) {
    cross_plant_pair_count
  } else 0
  estimated_output_bytes = max(1024^3, planned_pair_rows * 180 * 2)
  free_bytes = disk_free_bytes(out_dir)
  disk_status = if (is.na(free_bytes)) {
    "fail"
  } else if (free_bytes >= estimated_output_bytes) {
    "pass"
  } else {
    "fail"
  }

  checks = rbind(
    verify_input_manifest(input, manifest_file),
    release_manifest_check(release_manifest, package_version),
    data.frame(
      check = c(
        "membership_nonempty", "species_nonmissing", "compound_id_nonmissing",
        "membership_key_unique", "identity_exactly_one_inchikey_or_cid",
        "compound_id_maps_one_input_identity", "review_required_rows_absent",
        "membership_species_in_universe", "species_universe_unique",
        "fingerprint_payload_absent", "at_least_two_structures",
        "dependency_ChemmineR", "dependency_jsonlite", "output_directory_writable",
        "cache_directory_writable", "disk_space_adequate"
      ),
      status = c(
        ifelse(nrow(membership) > 0L, "pass", "fail"),
        ifelse(all(!is.na(membership$species)), "pass", "fail"),
        ifelse(all(!is.na(membership$compound_id)), "pass", "fail"),
        ifelse(anyDuplicated(membership_key) == 0L, "pass", "fail"),
        ifelse(all(xor(key_ok, cid_ok)), "pass", "fail"),
        ifelse(all(identity_count == 1L), "pass", "fail"),
        ifelse(review_count == 0L, "pass", "fail"),
        ifelse(length(membership_outside_universe) == 0L, "pass", "fail"),
        ifelse(anyDuplicated(species_universe) == 0L, "pass", "fail"),
        ifelse(!any(grepl("fingerprint", names(membership), ignore.case = TRUE)),
               "pass", "fail"),
        ifelse(unique_compounds >= 2L, "pass", "fail"),
        ifelse(requireNamespace("ChemmineR", quietly = TRUE), "pass", "fail"),
        ifelse(requireNamespace("jsonlite", quietly = TRUE), "pass", "fail"),
        ifelse(test_writable(out_dir), "pass", "fail"),
        ifelse(test_writable(cache_dir), "pass", "fail"),
        disk_status
      ),
      detail = c(
        paste(nrow(membership), "membership rows"),
        paste(sum(is.na(membership$species)), "missing species labels"),
        paste(sum(is.na(membership$compound_id)), "missing compound IDs"),
        paste(anyDuplicated(membership_key), "first duplicated-key index"),
        paste(sum(!xor(key_ok, cid_ok)), "identity-format failures"),
        paste(sum(identity_count != 1L), "rows have an inconsistent compound identity"),
        paste(review_count, "review-required rows"),
        paste(length(membership_outside_universe),
              "membership species absent from the supplied universe"),
        paste(length(species_universe), "unique species in the universe"),
        paste(grep("fingerprint", names(membership), ignore.case = TRUE,
                   value = TRUE), collapse = "; "),
        paste(unique_compounds, "unique structures"),
        "ChemmineR is required to decode PubChem Fingerprint2D.",
        "jsonlite is required for lifecycle and PubChem JSON.",
        out_dir,
        cache_dir,
        paste("available", ifelse(is.na(free_bytes), "unknown",
                                  format(free_bytes, scientific = FALSE)),
              "bytes; conservative requirement",
              format(estimated_output_bytes, scientific = FALSE), "bytes")
      ),
      stringsAsFactors = FALSE
    )
  )
  if (planned_pair_rows > 250000000 && identical(mode, "full")) {
    extreme_checks = data.frame(
      check = c("summary_completed_before_extreme_full",
                "extreme_output_explicitly_confirmed"),
      status = c(
        ifelse(file.exists(file.path(out_dir, "SUMMARY_COMPLETED.txt")),
               "pass", "fail"),
        ifelse(extreme_output_confirm, "pass", "fail")
      ),
      detail = c(
        "Outputs above 250 million rows require a completed summary gate.",
        "Pass --extreme-output-confirm true only after reviewing estimates."
      ),
      stringsAsFactors = FALSE
    )
    checks = rbind(checks, extreme_checks)
  }

  validation_status = ifelse(any(checks$status == "fail"), "fail", "pass")
  preflight = data.frame(
    mode = mode,
    input_file = input,
    input_md5 = unname(tools::md5sum(input)[[1L]]),
    membership_row_count = nrow(membership),
    plant_count = length(species_universe),
    plants_with_input_membership = length(species_counts),
    species_universe_file = ifelse(is.null(species_universe_file),
                                   NA_character_, species_universe_file),
    unique_structure_count = unique_compounds,
    compound_pair_count = format(compound_pair_count, scientific = FALSE),
    plant_pair_count = format(plant_pair_count, scientific = FALSE),
    cross_plant_membership_pair_count =
      format(cross_plant_pair_count, scientific = FALSE),
    full_cross_plant_pair_output_requested = full_pairs,
    requested_pair_rows = format(requested_pair_rows, scientific = FALSE),
    estimated_required_disk_bytes =
      format(estimated_output_bytes, scientific = FALSE),
    available_disk_bytes = ifelse(is.na(free_bytes), NA_character_,
                                  format(free_bytes, scientific = FALSE)),
    package_version = package_version,
    throttle_seconds = throttle,
    pair_block_size = pair_block_size,
    pair_shard_rows = pair_shard_rows,
    validation_status = validation_status,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  atomic_csv(checks, file.path(out_dir, "server_tanimoto_preflight_checks.csv"))
  atomic_csv(preflight, file.path(out_dir, "server_tanimoto_preflight_summary.csv"))
  counts = as.list(preflight[1L, c("membership_row_count", "plant_count",
                                   "plants_with_input_membership",
                                   "unique_structure_count", "compound_pair_count",
                                   "plant_pair_count",
                                   "cross_plant_membership_pair_count")])
  write_status(status_path, run_id, mode,
               ifelse(validation_status == "pass", "preflight_passed",
                      "preflight_failed"),
               input, out_dir, package_version,
               message = paste("Preflight", validation_status), counts = counts)
  print(preflight, row.names = FALSE)
  if (validation_status == "fail") {
    print(checks[checks$status == "fail", ], row.names = FALSE)
    stop("Server Tanimoto preflight failed; no PubChem or pairwise work started.",
         call. = FALSE)
  }
  if (identical(mode, "preflight")) {
    if (!mode_supplied) {
      cat("\nPreflight passed. No PubChem or pairwise work was started.\n")
      cat("Set UAFR_CONFIRM_SERVER_TANIMOTO=YES and rerun this exact command.\n")
      return(2L)
    }
    cat("Server Tanimoto preflight passed; no network work was requested.\n")
    return(0L)
  }
  if (!confirmed) {
    stop("Live modes require UAFR_CONFIRM_SERVER_TANIMOTO=YES.",
         call. = FALSE)
  }

  run_membership = membership
  if (identical(mode, "smoke")) {
    selector = getFromNamespace(".tanimoto_smoke_membership", "uafR")
    run_membership = selector(membership, n = smoke_structures,
                              group_col = "species",
                              compound_id_col = "compound_id")
    run_membership$.smoke_selection_rank = NULL
  }
  mode_structure_count = length(unique(run_membership$compound_id))
  if (mode_structure_count < 2L) {
    stop("The selected mode has fewer than two structures.", call. = FALSE)
  }

  incomplete_root = file.path(out_dir, ".incomplete")
  staging_dir = file.path(incomplete_root, mode)
  final_dir = if (legacy_layout) out_dir else file.path(out_dir, mode)
  prior_markers = file.path(out_dir, c("PAUSED_SERVICE_BUSY.json", "FAILED.json"))
  prior_markers = prior_markers[file.exists(prior_markers)]
  if (length(prior_markers) > 0L) {
    marker_archive = file.path(incomplete_root, "archive", "markers")
    dir.create(marker_archive, recursive = TRUE, showWarnings = FALSE)
    for (marker in prior_markers) {
      archived = file.path(
        marker_archive,
        paste0(run_id, "_prior_", basename(marker))
      )
      if (!file.copy(marker, archived, overwrite = FALSE)) {
        stop("Could not archive prior lifecycle marker: ", marker,
             call. = FALSE)
      }
      unlink(marker, force = TRUE)
    }
  }
  if (dir.exists(staging_dir)) {
    if (!overwrite_incomplete) {
      stop("Incomplete staging exists: ", staging_dir,
           ". Inspect it, then pass --overwrite-incomplete true to restart ",
           "from the preserved PubChem cache.", call. = FALSE)
    }
    unlink(staging_dir, recursive = TRUE, force = TRUE)
  }
  if (!legacy_layout && dir.exists(final_dir) && !overwrite) {
    stop("Completed mode output already exists: ", final_dir,
         call. = FALSE)
  }
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  completed_successfully = FALSE
  on.exit({
    pause_marker = file.path(out_dir, "PAUSED_SERVICE_BUSY.json")
    failure_marker = file.path(out_dir, "FAILED.json")
    if (!completed_successfully && !file.exists(pause_marker) &&
        !file.exists(failure_marker)) {
      marker = list(
        workflow = "run_plant_tanimoto_server",
        run_id = run_id,
        mode = mode,
        state = "failed",
        message = paste(
          "Run exited before atomic publication. Inspect server progress and",
          "the archived incomplete staging directory."
        ),
        failed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      )
      atomic_json(marker, failure_marker)
      archive = if (dir.exists(staging_dir)) {
        archive_staging(staging_dir, incomplete_root, run_id, "failed")
      } else {
        NA_character_
      }
      write_status(status_path, run_id, mode, "failed", input, out_dir,
                   package_version,
                   paste(marker$message, "Incomplete staging:", archive),
                   counts)
    }
  }, add = TRUE)
  recorder = progress_recorder(progress_path, run_id, mode)
  write_status(status_path, run_id, mode, "running", input, out_dir,
               package_version, "Live Tanimoto mode started.", counts)
  recorder$callback(list(phase = "runner", event = "mode_started",
                         completed = 0, total = mode_structure_count,
                         detail = paste(mode_structure_count,
                                        "structures submitted.")))

  effective_shard_rows = Inf
  if (identical(mode, "full") &&
      max(compound_pair_count, if (full_pairs) cross_plant_pair_count else 0) >
      25000000) {
    effective_shard_rows = pair_shard_rows
  }
  started = Sys.time()
  result = tryCatch(
    uafR::chemicalTanimotoSimilarity(
      run_membership,
      compound_col = "compound_name",
      cid_col = "CID",
      inchikey_col = "InChIKey",
      smiles_col = if ("SMILES" %in% names(run_membership)) "SMILES" else NULL,
      compound_id_col = "compound_id",
      group_cols = "species",
      return_compound_pairs = mode %in% c("smoke", "full"),
      return_group_compound_pairs = identical(mode, "smoke") ||
        (identical(mode, "full") && full_pairs),
      out_dir = staging_dir,
      pair_block_size = pair_block_size,
      pair_shard_rows = effective_shard_rows,
      cache = TRUE,
      cache_dir = cache_dir,
      throttle = throttle,
      refresh = FALSE,
      request_fun = request_fun,
      name_fallback = FALSE,
      progress_fun = recorder$callback,
      progress_every = progress_every,
      service_busy_limit = service_busy_limit
    ),
    uaf_pubchem_service_busy = function(error) {
      marker = list(
        workflow = "run_plant_tanimoto_server",
        run_id = run_id,
        mode = mode,
        state = "paused_service_busy",
        status_code = error$status_code,
        url = error$url,
        retry_after_seconds = error$retry_after_seconds,
        cache_dir = cache_dir,
        message = conditionMessage(error),
        paused_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      )
      atomic_json(marker, file.path(out_dir, "PAUSED_SERVICE_BUSY.json"))
      archive = archive_staging(staging_dir, incomplete_root, run_id, "paused")
      write_status(status_path, run_id, mode, "paused_service_busy", input,
                   out_dir, package_version,
                   paste(conditionMessage(error), "Incomplete staging:", archive),
                   counts)
      structure(list(status = 75L), class = "uaf_server_exit")
    },
    error = function(error) {
      marker = list(
        workflow = "run_plant_tanimoto_server",
        run_id = run_id,
        mode = mode,
        state = "failed",
        message = conditionMessage(error),
        failed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      )
      atomic_json(marker, file.path(out_dir, "FAILED.json"))
      archive = archive_staging(staging_dir, incomplete_root, run_id, "failed")
      write_status(status_path, run_id, mode, "failed", input, out_dir,
                   package_version,
                   paste(conditionMessage(error), "Incomplete staging:", archive),
                   counts)
      stop(error)
    }
  )
  if (inherits(result, "uaf_server_exit")) return(result$status)
  finished = Sys.time()

  identity = result$CompoundIdentityMap
  canonical = result$CanonicalCompounds
  canonical_membership = result$GroupCompoundMembership
  submitted = mode_structure_count
  usable = nrow(canonical)
  mismatched = sum(identity$identity_match_status ==
                     "input_pubchem_inchikey_mismatch", na.rm = TRUE)
  mismatch_ids = unique(identity$canonical_compound_id[
    identity$identity_match_status == "input_pubchem_inchikey_mismatch"
  ])
  mismatch_ids = mismatch_ids[!is.na(mismatch_ids) & mismatch_ids != ""]
  mismatch_entered = intersect(mismatch_ids, canonical_membership$compound_id)
  usable_fraction = usable / submitted
  actual_compound_pairs = choose(usable, 2)
  actual_cross_pairs = if (nrow(canonical_membership) >= 2L) {
    total = choose(nrow(canonical_membership), 2)
    within = sum(vapply(
      split(canonical_membership$compound_id, canonical_membership$group_id),
      function(x) choose(length(x), 2), numeric(1)
    ))
    total - within
  } else 0
  written_compound_pairs = sum(suppressWarnings(as.numeric(
    result$ExportManifest$Rows[
      result$ExportManifest$Table == "CompoundTanimoto"]
  )), na.rm = TRUE)
  written_group_pairs = sum(suppressWarnings(as.numeric(
    result$ExportManifest$Rows[
      result$ExportManifest$Table == "GroupCompoundTanimoto"]
  )), na.rm = TRUE)
  complete_summary = complete_plant_pair_summary(
    result$GroupPairTanimotoSummary, species_universe,
    canonical_membership
  )
  scope_summary = result$ComparableScopeTanimotoSummary
  group_summary = result$ComparableGroupTanimotoSummary
  summary_keys = complete_summary$group_pair_id
  output_checks = data.frame(
    check = c(
      "at_least_two_canonical_structures", "canonical_compound_ids_unique",
      "fingerprint_usable_fraction", "canonical_membership_keys_unique",
      "mismatches_absent_from_canonical_membership",
      "plant_pair_summary_keys_unique", "plant_pair_universe_complete",
      "comparable_scope_keys_unique", "comparable_group_keys_unique",
      "compound_pair_rows_complete",
      "cross_plant_pair_rows_complete", "smoke_identity_mismatches_absent"
    ),
    status = c(
      ifelse(usable >= 2L, "pass", "fail"),
      ifelse(anyDuplicated(canonical$compound_id) == 0L, "pass", "fail"),
      ifelse(usable_fraction >= min_fingerprint_usable, "pass", "fail"),
      ifelse(anyDuplicated(paste(canonical_membership$group_id,
                                 canonical_membership$compound_id,
                                 sep = "\r")) == 0L, "pass", "fail"),
      ifelse(all(canonical_membership$compound_id %in% canonical$compound_id) &&
               length(mismatch_entered) == 0L,
             "pass", "fail"),
      ifelse(anyDuplicated(summary_keys) == 0L, "pass", "fail"),
      ifelse(nrow(complete_summary) == choose(length(species_universe), 2L),
             "pass", "fail"),
      ifelse(anyDuplicated(paste(scope_summary$species_a,
                                 scope_summary$species_b,
                                 scope_summary$comparison_value,
                                 sep = "\r")) == 0L, "pass", "fail"),
      ifelse(anyDuplicated(paste(group_summary$species_a,
                                 group_summary$species_b,
                                 group_summary$comparison_value,
                                 sep = "\r")) == 0L, "pass", "fail"),
      ifelse(!mode %in% c("smoke", "full") ||
               written_compound_pairs == actual_compound_pairs,
             "pass", "fail"),
      ifelse(!(identical(mode, "smoke") ||
                 (identical(mode, "full") && full_pairs)) ||
               written_group_pairs == actual_cross_pairs,
             "pass", "fail"),
      ifelse(!identical(mode, "smoke") || mismatched == 0L,
             "pass", "fail")
    ),
    detail = c(
      paste(usable, "canonical structures"),
      paste(anyDuplicated(canonical$compound_id), "first duplicate index"),
      paste(round(usable_fraction, 4), "usable; minimum",
            min_fingerprint_usable),
      paste(nrow(canonical_membership), "canonical membership rows"),
      paste(mismatched, "mismatches detected;", length(mismatch_entered),
            "entered canonical membership"),
      paste(length(summary_keys), "plant-pair summary rows"),
      paste(nrow(complete_summary), "written;",
            choose(length(species_universe), 2L), "expected"),
      paste(nrow(scope_summary), "scope-filtered summary rows"),
      paste(nrow(group_summary), "group-filtered summary rows"),
      paste(written_compound_pairs, "written;", actual_compound_pairs,
            "expected"),
      paste(written_group_pairs, "written;", actual_cross_pairs, "expected"),
      paste(mismatched, "smoke mismatches")
    ),
    stringsAsFactors = FALSE
  )
  atomic_csv(output_checks, file.path(staging_dir,
                                      "server_tanimoto_validation.csv"))
  if (any(output_checks$status == "fail")) {
    stop("Completed calculations failed server output validation: ",
         paste(output_checks$check[output_checks$status == "fail"],
               collapse = ", "), call. = FALSE)
  }

  final_prefix = normalizePath(final_dir, mustWork = FALSE)
  staging_prefix = normalizePath(staging_dir, mustWork = TRUE)
  if (is.data.frame(result$ExportManifest) && nrow(result$ExportManifest) > 0L) {
    result$ExportManifest$Path = sub(staging_prefix, final_prefix,
                                     result$ExportManifest$Path, fixed = TRUE)
  }
  tables = list(
    plant_pair_tanimoto_summary = complete_summary,
    comparable_scope_tanimoto_summary = scope_summary,
    comparable_group_tanimoto_summary = group_summary,
    compound_identity_map = result$CompoundIdentityMap,
    excluded_compound_identities = result$ExcludedCompoundIdentities,
    canonical_compounds = result$CanonicalCompounds,
    plant_compound_membership_evidence = result$GroupCompoundMembershipEvidence,
    plant_compound_membership_canonical = result$GroupCompoundMembership,
    pubchem_fingerprint_audit = result$PubChemFingerprints
  )
  table_rows = list()
  for (name in names(tables)) {
    file_name = paste0(name, ".csv")
    atomic_csv(tables[[name]], file.path(staging_dir, file_name))
    table_rows[[file_name]] = nrow(tables[[name]])
  }
  provenance = data.frame(
    workflow = "run_plant_tanimoto_server",
    mode = mode,
    run_id = run_id,
    package_version = package_version,
    input_file = input,
    input_md5 = unname(tools::md5sum(input)[[1L]]),
    cache_dir = cache_dir,
    started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    finished_at = format(finished, "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_seconds = round(as.numeric(difftime(finished, started,
                                                units = "secs")), 3),
    requested_structures = submitted,
    resolved_structures = sum(!is.na(result$CompoundResolution$pubchem_cid)),
    fingerprint_usable_structures = usable,
    excluded_identity_rows = nrow(result$ExcludedCompoundIdentities),
    mismatched_identity_rows = mismatched,
    canonical_structure_count = usable,
    request_count = recorder$metrics$request_count,
    cache_hit_count = recorder$metrics$cache_hit_count,
    error_count = recorder$metrics$error_count,
    rate_limit_count = recorder$metrics$rate_limit_count,
    full_cross_plant_pairs = full_pairs,
    pair_shard_rows = effective_shard_rows,
    name_fallback = FALSE,
    fingerprint_source = "PubChem Fingerprint2D",
    stringsAsFactors = FALSE
  )
  atomic_csv(provenance, file.path(staging_dir, "server_tanimoto_provenance.csv"))
  table_rows[["server_tanimoto_provenance.csv"]] = nrow(provenance)
  table_rows[["server_tanimoto_validation.csv"]] = nrow(output_checks)
  saveRDS(result, file.path(staging_dir, "server_tanimoto_result.rds"),
          compress = "xz")
  publish_staging(staging_dir, final_dir, legacy_layout, overwrite)

  mode_manifest = build_output_manifest(final_dir, out_dir, mode, run_id,
                                        result, table_rows)
  manifest_path = file.path(out_dir, "server_tanimoto_output_manifest.csv")
  existing_manifest = if (file.exists(manifest_path)) {
    utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                    check.names = FALSE)
  } else data.frame()
  if (nrow(existing_manifest) > 0L && "mode" %in% names(existing_manifest)) {
    existing_manifest = existing_manifest[existing_manifest$mode != mode,
                                          , drop = FALSE]
  }
  atomic_csv(rbind(existing_manifest, mode_manifest), manifest_path)

  completion = paste("run_id:", run_id,
                     "completed_at:", format(Sys.time(),
                                               "%Y-%m-%dT%H:%M:%S%z"),
                     "mode:", mode, sep = "\n")
  atomic_text(completion, file.path(out_dir,
                                    "FINGERPRINT_RESOLUTION_COMPLETED.txt"))
  if (mode %in% c("summary", "full")) {
    atomic_text(completion, file.path(out_dir, "SUMMARY_COMPLETED.txt"))
  }
  if (identical(mode, "full")) {
    atomic_text(completion, file.path(out_dir, "FULL_PAIRWISE_COMPLETED.txt"))
  }
  unlink(file.path(out_dir, c("PAUSED_SERVICE_BUSY.json", "FAILED.json")),
         force = TRUE)
  recorder$callback(list(phase = "runner", event = "mode_completed",
                         completed = submitted, total = submitted,
                         detail = paste("Completed", mode, "mode.")))
  write_status(status_path, run_id, mode, "completed", input, out_dir,
               package_version, paste("Completed", mode, "mode."),
               c(counts, list(
                 requested_structures = submitted,
                 fingerprint_usable_structures = usable,
                 output_directory = normalizePath(final_dir, mustWork = TRUE)
               )))
  completed_successfully = TRUE
  cat("Server Tanimoto", mode, "mode completed.\n")
  cat("Output:", normalizePath(final_dir, mustWork = TRUE), "\n")
  0L
}

if (!identical(Sys.getenv("UAFR_SERVER_RUNNER_SOURCE_ONLY", ""), "true")) {
  status = tryCatch(main(), error = function(error) {
    message("Error: ", conditionMessage(error))
    1L
  })
  quit(save = "no", status = as.integer(status))
}
