#!/usr/bin/env Rscript

file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path = if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[[1]])
} else {
  file.path(getwd(), "tools", "run_plant_phytochemistry_pilot.R")
}
repo_root = normalizePath(file.path(dirname(script_path), ".."),
                          winslash = "/", mustWork = FALSE)

load_uafr = function(repo_root) {
  if (file.exists(file.path(repo_root, "DESCRIPTION")) &&
      requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(repo_root, quiet = TRUE)
    return(invisible(TRUE))
  }
  if (requireNamespace("uafR", quietly = TRUE) &&
      exists("runPlantPhytochemistryPilot", asNamespace("uafR"),
             inherits = FALSE)) {
    suppressPackageStartupMessages(library(uafR))
    return(invisible(TRUE))
  }
  stop("Could not load uafR with runPlantPhytochemistryPilot(). Install the ",
       "current package or run this script from the uafR repository with ",
       "devtools installed.", call. = FALSE)
}

parse_args = function(args) {
  out = list()
  i = 1L
  while (i <= length(args)) {
    arg = args[[i]]
    if (!grepl("^--", arg)) {
      stop("Unexpected argument: ", arg, call. = FALSE)
    }
    key = sub("^--", "", arg)
    value = TRUE
    if (grepl("=", key, fixed = TRUE)) {
      parts = strsplit(key, "=", fixed = TRUE)[[1]]
      key = parts[[1]]
      value = paste(parts[-1], collapse = "=")
    } else if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
      value = args[[i + 1L]]
      i = i + 1L
    }
    key = gsub("-", "_", key)
    out[[key]] = value
    i = i + 1L
  }
  out
}

usage = function() {
  cat(
    "Usage:\n",
    "Rscript tools/run_plant_phytochemistry_pilot.R \\\n",
    "  --plant-csv path/to/plants.csv \\\n",
    "  --species-col species \\\n",
    "  --out-dir path/to/plant_pilot \\\n",
    "  --cache-dir path/to/uafR_plant_cache \\\n",
    "  --lotus-index path/to/lotus_compact_index.csv \\\n",
    "  --sources lotus,knapsack,npass,pubchem,pubmed,pubtator \\\n",
    "  --compound-resolution-profile identity \\\n",
    "  --max-pubmed-records 25 \\\n",
    "  --max-provider-records 100 \\\n",
    "  --request-timeout 30 \\\n",
    "  --overwrite\n\n",
    "Or use the built-in 15-species pilot panel:\n",
    "Rscript tools/run_plant_phytochemistry_pilot.R \\\n",
    "  --default-panel \\\n",
    "  --out-dir plant_phytochemistry_pilot \\\n",
    "  --compound-resolution-profile identity \\\n",
    "  --overwrite\n\n",
    "Notes:\n",
    "- The default run is a small pilot workflow with cached, resumable ",
    "provider discovery.\n",
    "- Use --compound-resolution-profile research only after the identity ",
    "pilot behaves as expected.\n",
    "- Use --lotus-index, or set UAFR_LOTUS_INDEX, for scalable LOTUS ",
    "species-compound discovery.\n",
    "- NCBI_EMAIL, NCBI_TOOL, and NCBI_API_KEY can be supplied through the ",
    "environment or explicit arguments.\n",
    sep = ""
  )
}

required_arg = function(args, name) {
  value = args[[name]]
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    stop("Missing required argument `--", gsub("_", "-", name), "`.",
         call. = FALSE)
  }
  value
}

optional_arg = function(args, name, default = NULL) {
  value = args[[name]]
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    return(default)
  }
  value
}

flag_arg = function(args, name, default = FALSE) {
  value = args[[name]]
  if (is.null(value)) return(default)
  if (identical(value, TRUE)) return(TRUE)
  tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

split_arg = function(value, default) {
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    return(default)
  }
  pieces = unlist(strsplit(value, "\\s*[,;]\\s*", perl = TRUE),
                  use.names = FALSE)
  pieces = trimws(pieces)
  pieces[nzchar(pieces)]
}

numeric_arg = function(args, name, default) {
  value = optional_arg(args, name, NULL)
  if (is.null(value)) return(default)
  if (tolower(value) %in% c("inf", "infinite")) return(Inf)
  out = suppressWarnings(as.numeric(value))
  if (length(out) != 1 || is.na(out) ||
      (!is.finite(out) && !is.infinite(out))) {
    stop("Argument `--", gsub("_", "-", name), "` must be numeric.",
         call. = FALSE)
  }
  out
}

integer_arg = function(args, name, default) {
  as.integer(numeric_arg(args, name, default))
}

resolve_species_col = function(table, species_col) {
  if (species_col %in% names(table)) return(species_col)
  hit = which(tolower(names(table)) == tolower(species_col))
  if (length(hit) > 0) return(names(table)[[hit[[1]]]])
  if (identical(species_col, "species") && ncol(table) >= 1) {
    message("Column `species` not found; using first CSV column `",
            names(table)[[1]], "`.")
    return(names(table)[[1]])
  }
  stop("Plant CSV does not contain species column `", species_col, "`.",
       call. = FALSE)
}

main = function() {
  args = commandArgs(trailingOnly = TRUE)
  if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
    usage()
    quit(status = ifelse(length(args) == 0, 1L, 0L))
  }
  parsed = parse_args(args)
  load_uafr(repo_root)

  out_dir = required_arg(parsed, "out_dir")
  if (flag_arg(parsed, "default_panel", FALSE)) {
    panel = plantPhytochemistryPilotPanel(
      optional_arg(parsed, "panel_profile", "general")
    )
    plants = panel$species
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(panel, file.path(out_dir, "pilot_species_panel.csv"),
                     row.names = FALSE, na = "")
  } else {
    plant_csv = required_arg(parsed, "plant_csv")
    plant_table = utils::read.csv(plant_csv, stringsAsFactors = FALSE,
                                  check.names = FALSE)
    species_col = resolve_species_col(
      plant_table,
      optional_arg(parsed, "species_col", "species")
    )
    plants = plant_table[[species_col]]
  }

  result = runPlantPhytochemistryPilot(
    plants = plants,
    out_dir = out_dir,
    sources = split_arg(optional_arg(parsed, "sources", NULL),
                        c("lotus", "knapsack", "npass", "pubchem",
                          "pubmed", "pubtator")),
    taxon_fallback = split_arg(optional_arg(parsed, "taxon_fallback", NULL),
                               c("species", "genus")),
    cache_dir = optional_arg(parsed, "cache_dir", file.path(out_dir, "cache")),
    species_chunk_size = integer_arg(parsed, "species_chunk_size", 25),
    compound_resolution_profile =
      optional_arg(parsed, "compound_resolution_profile", "identity"),
    min_confidence = optional_arg(parsed, "min_confidence", "medium"),
    max_compounds_per_species =
      numeric_arg(parsed, "max_compounds_per_species", Inf),
    max_unique_compounds = numeric_arg(parsed, "max_unique_compounds", Inf),
    cache = !flag_arg(parsed, "no_cache", FALSE),
    throttle = numeric_arg(parsed, "throttle", 0.5),
    ncbi_email = optional_arg(parsed, "ncbi_email",
                              Sys.getenv("NCBI_EMAIL", "")),
    ncbi_tool = optional_arg(parsed, "ncbi_tool",
                             Sys.getenv("NCBI_TOOL", "uafR")),
    ncbi_api_key = optional_arg(parsed, "ncbi_api_key",
                                Sys.getenv("NCBI_API_KEY", "")),
    max_pubmed_records = integer_arg(parsed, "max_pubmed_records", 25),
    max_provider_records = integer_arg(parsed, "max_provider_records", 100),
    lotus_index = optional_arg(parsed, "lotus_index",
                               Sys.getenv("UAFR_LOTUS_INDEX", "")),
    request_timeout = numeric_arg(parsed, "request_timeout", 30),
    compound_batch_size = integer_arg(parsed, "compound_batch_size", 50),
    matrix_feature = optional_arg(parsed, "matrix_feature",
                                  "comparison_group"),
    matrix_mode = optional_arg(parsed, "matrix_mode", "binary"),
    max_matrix_features = numeric_arg(parsed, "max_matrix_features", Inf),
    resume = !flag_arg(parsed, "no_resume", FALSE),
    progress = flag_arg(parsed, "progress", TRUE),
    overwrite = flag_arg(parsed, "overwrite", FALSE),
    strict = flag_arg(parsed, "strict", FALSE),
    stop_on_error = flag_arg(parsed, "stop_on_error", FALSE),
    refresh = flag_arg(parsed, "refresh", FALSE)
  )

  cat("Plant phytochemistry pilot complete.\n")
  cat("Output directory: ", normalizePath(out_dir, winslash = "/",
                                         mustWork = FALSE), "\n", sep = "")
  cat("Pilot summary: ", file.path(out_dir, "pilot_summary.csv"), "\n",
      sep = "")
  cat("Review table: ", file.path(out_dir, "pilot_review_needed.csv"), "\n",
      sep = "")
  cat("Pilot manifest: ", file.path(out_dir, "pilot_export_manifest.csv"),
      "\n", sep = "")
  cat("QA report: ", file.path(out_dir, "pilot_qa_report.csv"), "\n",
      sep = "")
  cat("Species processed: ", nrow(result$PilotSummary), "\n", sep = "")
}

main()
