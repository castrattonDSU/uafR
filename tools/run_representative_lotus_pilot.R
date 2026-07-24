#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path = if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  file.path(getwd(), "tools", "run_representative_lotus_pilot.R")
}
repo_root = normalizePath(file.path(dirname(script_path), ".."),
                          winslash = "/", mustWork = FALSE)

parse_args = function(args) {
  out = list()
  i = 1L
  while (i <= length(args)) {
    token = args[[i]]
    if (token %in% c("-h", "--help")) {
      usage()
      quit(save = "no", status = 0L)
    }
    if (!startsWith(token, "--") || i == length(args)) {
      stop("Arguments must use --name value pairs.", call. = FALSE)
    }
    out[[gsub("-", "_", substring(token, 3), fixed = TRUE)]] = args[[i + 1L]]
    i = i + 2L
  }
  out
}

usage = function() {
  cat(paste(
    "Run a deterministic, source-local LOTUS panel pilot.",
    "",
    "Usage:",
    "  Rscript tools/run_representative_lotus_pilot.R \\",
    "    --plant-csv plant_species_run_input.csv \\",
    "    --lotus-index LOTUS_lookup_index \\",
    "    --out-dir representative_lotus_pilot \\",
    "    --cache-dir representative_lotus_pilot_cache \\",
    "    [--species-col species] \\",
    "    [--exact-col lotus_exact_species_key_available] \\",
    "    [--genus-col lotus_genus_key_available] \\",
    "    [--exact-count 10] [--genus-only-count 10] [--no-key-count 5] \\",
    "    [--overwrite true]",
    "",
    "The pilot uses only the supplied local LOTUS index, species-level",
    "matching, and no compound enrichment. It does not contact PubChem,",
    "PubMed, PubTator, KNApSAcK, NPASS, or live LOTUS services.",
    sep = "\n"
  ))
}

required_arg = function(args, name) {
  value = args[[name]]
  if (is.null(value) || is.na(value) || !nzchar(value)) {
    stop("Missing required --", gsub("_", "-", name, fixed = TRUE),
         " argument.", call. = FALSE)
  }
  value
}

integer_arg = function(args, name, default) {
  value = suppressWarnings(as.integer(if (is.null(args[[name]])) {
    default
  } else args[[name]]))
  if (length(value) != 1L || is.na(value) || value < 0L) {
    stop("`--", gsub("_", "-", name, fixed = TRUE),
         "` must be a non-negative integer.", call. = FALSE)
  }
  value
}

flag_arg = function(args, name, default = FALSE) {
  if (is.null(args[[name]])) return(default)
  value = tolower(trimws(args[[name]]))
  if (!value %in% c("true", "false", "yes", "no", "1", "0")) {
    stop("`--", gsub("_", "-", name, fixed = TRUE),
         "` must be true or false.", call. = FALSE)
  }
  value %in% c("true", "yes", "1")
}

as_logical_column = function(x, name) {
  if (is.logical(x)) return(x %in% TRUE)
  value = tolower(trimws(as.character(x)))
  invalid = !is.na(value) & nzchar(value) &
    !value %in% c("true", "false", "yes", "no", "1", "0")
  if (any(invalid)) {
    stop("Column `", name, "` contains non-logical values: ",
         paste(utils::head(unique(value[invalid]), 10L), collapse = "; "),
         call. = FALSE)
  }
  value %in% c("true", "yes", "1")
}

spread_indices = function(n, size) {
  if (size == 0L) return(integer())
  if (n < size) stop("Requested ", size, " rows from a stratum with only ", n,
                     " candidates.", call. = FALSE)
  indices = unique(as.integer(round(seq(1, n, length.out = size))))
  if (length(indices) < size) {
    indices = c(indices, setdiff(seq_len(n), indices))
  }
  sort(indices[seq_len(size)])
}

select_panel = function(plants, species_col, exact_col, genus_col,
                        exact_count, genus_only_count, no_key_count) {
  required = c(species_col, exact_col, genus_col)
  missing = setdiff(required, names(plants))
  if (length(missing) > 0L) {
    stop("Plant input is missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  species = trimws(as.character(plants[[species_col]]))
  exact = as_logical_column(plants[[exact_col]], exact_col)
  genus = as_logical_column(plants[[genus_col]], genus_col)
  if (any(is.na(species) | !nzchar(species))) {
    stop("Plant input contains blank species names.", call. = FALSE)
  }
  if (anyDuplicated(species) > 0L) {
    stop("Plant input species must be unique before pilot selection.",
         call. = FALSE)
  }

  strata = list(
    exact_species_key = which(exact),
    genus_key_no_exact_species_key = which(!exact & genus),
    no_species_or_genus_key = which(!exact & !genus)
  )
  requested = c(exact_count, genus_only_count, no_key_count)
  rows = lapply(seq_along(strata), function(i) {
    candidates = plants[strata[[i]], , drop = FALSE]
    candidates$.pilot_species = species[strata[[i]]]
    candidates = candidates[order(candidates$.pilot_species), , drop = FALSE]
    chosen = candidates[spread_indices(nrow(candidates), requested[[i]]),
                        , drop = FALSE]
    chosen$pilot_stratum = names(strata)[[i]]
    chosen$pilot_stratum_rank = seq_len(nrow(chosen))
    chosen$.pilot_species = NULL
    chosen
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out$pilot_selection_order = seq_len(nrow(out))
  out
}

main = function(command_args = commandArgs(trailingOnly = TRUE)) {
  if (length(command_args) < 1L) {
    usage()
    return(1L)
  }
  args = parse_args(command_args)
  plant_csv = normalizePath(required_arg(args, "plant_csv"), mustWork = TRUE)
  lotus_index = normalizePath(required_arg(args, "lotus_index"),
                              mustWork = TRUE)
  out_dir = required_arg(args, "out_dir")
  cache_dir = required_arg(args, "cache_dir")
  species_col = if (is.null(args$species_col)) "species" else args$species_col
  exact_col = if (is.null(args$exact_col)) {
    "lotus_exact_species_key_available"
  } else args$exact_col
  genus_col = if (is.null(args$genus_col)) {
    "lotus_genus_key_available"
  } else args$genus_col
  exact_count = integer_arg(args, "exact_count", 10L)
  genus_only_count = integer_arg(args, "genus_only_count", 10L)
  no_key_count = integer_arg(args, "no_key_count", 5L)
  overwrite = flag_arg(args, "overwrite", FALSE)

  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required to load the current uafR source tree.",
         call. = FALSE)
  }
  devtools::load_all(repo_root, quiet = TRUE)
  plants = utils::read.csv(plant_csv, stringsAsFactors = FALSE,
                           check.names = FALSE, fileEncoding = "UTF-8")
  panel = select_panel(
    plants, species_col, exact_col, genus_col,
    exact_count, genus_only_count, no_key_count
  )
  expected_size = exact_count + genus_only_count + no_key_count
  if (nrow(panel) != expected_size || anyDuplicated(panel[[species_col]]) > 0L) {
    stop("Representative panel selection did not produce the requested unique",
         " species count.", call. = FALSE)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  protected = file.path(out_dir, c("representative_pilot_species.csv",
                                   "representative_pilot_manifest.csv"))
  if (!overwrite && any(file.exists(protected))) {
    stop("Representative pilot output exists; preserve it or pass",
         " --overwrite true deliberately.", call. = FALSE)
  }
  uafR:::.plant_atomic_write_csv(
    panel, file.path(out_dir, "representative_pilot_species.csv")
  )
  selection_summary = data.frame(
    pilot_stratum = names(table(panel$pilot_stratum)),
    selected_species_count = as.integer(table(panel$pilot_stratum)),
    stringsAsFactors = FALSE
  )
  uafR:::.plant_atomic_write_csv(
    selection_summary,
    file.path(out_dir, "representative_pilot_selection_summary.csv")
  )

  started = Sys.time()
  result = runPlantPhytochemistryPilot(
    plants = panel[[species_col]],
    out_dir = out_dir,
    sources = "lotus",
    taxon_fallback = "species",
    cache_dir = cache_dir,
    species_chunk_size = expected_size,
    compound_resolution_profile = "none",
    max_compounds_per_species = Inf,
    max_unique_compounds = Inf,
    cache = TRUE,
    throttle = 0,
    max_pubmed_records = 0,
    max_provider_records = Inf,
    lotus_index = lotus_index,
    compound_batch_size = 25,
    resume = TRUE,
    progress = TRUE,
    overwrite = overwrite,
    strict = FALSE,
    stop_on_error = FALSE,
    refresh = FALSE
  )
  finished = Sys.time()

  occurrences = result$PlantCompoundOccurrences
  occurrence_count = table(factor(occurrences$species,
                                  levels = panel[[species_col]]))
  accounting = panel[, c(species_col, "pilot_stratum",
                         "pilot_selection_order"), drop = FALSE]
  names(accounting)[names(accounting) == species_col] = "species"
  accounting$occurrence_count = as.integer(occurrence_count)
  accounting$query_outcome = ifelse(accounting$occurrence_count > 0L,
                                     "direct_species_records", "explicit_no_hit")
  uafR:::.plant_atomic_write_csv(
    accounting, file.path(out_dir, "representative_pilot_query_accounting.csv")
  )

  manifest_check = validatePlantChemistryRunManifest(
    result$BatchChunkManifest, base_dir = out_dir
  )
  diagnostics = result$ProviderDiagnostics
  provider_errors = if (is.data.frame(diagnostics) &&
                         "error_count" %in% names(diagnostics)) {
    sum(suppressWarnings(as.numeric(diagnostics$error_count)), na.rm = TRUE)
  } else 0
  source_ids_complete = nrow(occurrences) == 0L ||
    ("source_record_id" %in% names(occurrences) &&
       all(!is.na(occurrences$source_record_id) &
             nzchar(trimws(occurrences$source_record_id))))
  source_names = if (nrow(occurrences) > 0L) {
    unique(tolower(trimws(occurrences$source_database)))
  } else character()
  resolution = result$CompoundResolution
  resolution_disabled = !is.data.frame(resolution) || nrow(resolution) == 0L ||
    (all(!(resolution$resolved %in% TRUE)) &&
       all(tolower(trimws(resolution$resolution_source)) == "not_attempted") &&
       all(is.na(resolution$CID) |
             !nzchar(trimws(as.character(resolution$CID)))) &&
       all(is.na(resolution$InChIKey) |
             !nzchar(trimws(as.character(resolution$InChIKey)))) &&
       all(is.na(resolution$SMILES) |
             !nzchar(trimws(as.character(resolution$SMILES)))))
  validation = data.frame(
    check = c(
      "requested_panel_size", "species_unique", "strata_counts_exact",
      "all_queries_accounted_for", "batch_manifest_valid",
      "provider_errors_absent", "occurrence_species_within_panel",
      "source_record_ids_complete", "local_lotus_only",
      "compound_resolution_disabled", "full_run_marker_not_managed_here"
    ),
    status = c(
      ifelse(nrow(panel) == expected_size, "pass", "fail"),
      ifelse(anyDuplicated(panel[[species_col]]) == 0L, "pass", "fail"),
      ifelse(identical(
        as.integer(table(factor(panel$pilot_stratum,
                                levels = c("exact_species_key",
                                  "genus_key_no_exact_species_key",
                                  "no_species_or_genus_key")))),
        c(exact_count, genus_only_count, no_key_count)
      ), "pass", "fail"),
      ifelse(nrow(accounting) == expected_size &&
               sum(accounting$occurrence_count) == nrow(occurrences),
             "pass", "fail"),
      ifelse(manifest_check$Summary$validation_status[[1L]] == "pass",
             "pass", "fail"),
      ifelse(provider_errors == 0, "pass", "fail"),
      ifelse(all(occurrences$species %in% panel[[species_col]]),
             "pass", "fail"),
      ifelse(source_ids_complete, "pass", "fail"),
      ifelse(length(source_names) == 0L ||
               all(source_names %in% c("lotus", "lotus_local_index")),
             "pass", "fail"),
      ifelse(resolution_disabled, "pass", "fail"),
      "pass"
    ),
    detail = c(
      paste(nrow(panel), "selected species"),
      paste(length(unique(panel[[species_col]])), "unique species"),
      paste(names(table(panel$pilot_stratum)),
            as.integer(table(panel$pilot_stratum)), collapse = "; "),
      paste(nrow(accounting), "queries;", nrow(occurrences),
            "occurrence rows"),
      manifest_check$Summary$validation_status[[1L]],
      paste(provider_errors, "provider errors"),
      paste(length(unique(occurrences$species)), "species with records"),
      paste(sum(!is.na(occurrences$source_record_id)),
            "rows with source record IDs"),
      paste(source_names, collapse = "; "),
      paste(nrow(resolution), "explicit not-attempted audit rows;",
            sum(resolution$resolved %in% TRUE), "resolved rows"),
      "This tool writes only its separate pilot output and cache directories."
    ),
    stringsAsFactors = FALSE
  )
  uafR:::.plant_atomic_write_csv(
    validation, file.path(out_dir, "representative_pilot_validation.csv")
  )
  provenance = data.frame(
    workflow = "run_representative_lotus_pilot",
    input_file = plant_csv,
    input_md5 = unname(tools::md5sum(plant_csv)[[1L]]),
    lotus_index = lotus_index,
    lotus_manifest_md5 = if (file.exists(file.path(lotus_index,
                                                   "manifest.json"))) {
      unname(tools::md5sum(file.path(lotus_index, "manifest.json"))[[1L]])
    } else NA_character_,
    uafR_version = as.character(utils::packageVersion("uafR")),
    sources = "lotus_local_index_only",
    taxon_fallback = "species",
    compound_resolution_profile = "none",
    selected_species_count = nrow(panel),
    occurrence_row_count = nrow(occurrences),
    species_with_records = length(unique(occurrences$species)),
    started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    finished_at = format(finished, "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_seconds = round(as.numeric(difftime(finished, started,
                                                units = "secs")), 3),
    full_panel_run_started = FALSE,
    interpretation_limit = paste(
      "LOTUS records are reported source evidence, not measurements in project",
      "samples; no-hit results are not evidence of biological absence."
    ),
    stringsAsFactors = FALSE
  )
  uafR:::.plant_atomic_write_csv(
    provenance, file.path(out_dir, "representative_pilot_provenance.csv")
  )
  if (any(validation$status == "fail")) {
    stop("Representative LOTUS pilot failed validation: ",
         paste(validation$check[validation$status == "fail"], collapse = ", "),
         call. = FALSE)
  }

  writeLines(
    c("Representative local LOTUS pilot completed and validated.",
      paste("Completed at:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))),
    file.path(out_dir, "REPRESENTATIVE_PILOT_COMPLETED.txt"), useBytes = TRUE
  )
  files = list.files(out_dir, full.names = TRUE, recursive = TRUE)
  files = files[file.info(files)$isdir %in% FALSE]
  files = files[basename(files) != "representative_pilot_manifest.csv"]
  relative = substring(files, nchar(normalizePath(out_dir, mustWork = TRUE)) + 2L)
  manifest = data.frame(
    file = relative,
    bytes = unname(file.info(files)$size),
    md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
  manifest = manifest[order(manifest$file), , drop = FALSE]
  uafR:::.plant_atomic_write_csv(
    manifest, file.path(out_dir, "representative_pilot_manifest.csv")
  )
  cat("Representative local LOTUS pilot completed.\n")
  cat("Output:", normalizePath(out_dir, mustWork = TRUE), "\n")
  cat("Species:", nrow(panel), " Occurrence rows:", nrow(occurrences), "\n")
  0L
}

if (!identical(Sys.getenv("UAFR_REPRESENTATIVE_PILOT_SOURCE_ONLY", ""),
               "true")) {
  status = tryCatch(main(), error = function(error) {
    message("Error: ", conditionMessage(error))
    1L
  })
  quit(save = "no", status = status)
}
