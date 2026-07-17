#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

usage = function() {
  cat(paste(
    "Run the resumable uafR plant chemistry production panel.",
    "",
    "Usage:",
    "  Rscript run_plant_chemistry_panel.R \\",
    "    --mode preflight|build-indexes|pilot|discovery|identity|research-enrichment|full-enrichment|tanimoto-handoff|finalize|all \\",
    "    --plant-csv plant_species.csv \\",
    "    --out-dir plant_phytochemistry_full_panel_v2 \\",
    "    [--config run_config.json] \\",
    "    [--species-col species] [--metadata-csv plant_species.csv] \\",
    "    [--exclusion-files file1.csv,file2.csv] \\",
    "    [--supporting-files taxonomy.csv,normalization.csv] \\",
    "    [--cache-dir persistent_cache] \\",
    "    [--lotus-index LOTUS_lookup_index] \\",
    "    [--npass-index NPASS_index] [--npass-raw-dir NPASS_raw] \\",
    "    [--download-npass false] \\",
    "    [--sources lotus,npass,knapsack,pubchem,pubmed,pubtator] \\",
    "    [--taxon-fallback genus,family] \\",
    "    [--expected-species-count 701] [--expected-exclusion-count 25] \\",
    "    [--species-chunk-size 25] [--compound-batch-size 25] \\",
    "    [--pilot-count 12] [--pilot-species Species one,Species two] \\",
    "    [--pilot-enrichment-compounds 25] \\",
    "    [--max-pubmed-records 100] [--max-provider-records Inf] \\",
    "    [--provider-throttle 0.5] [--knapsack-throttle 1] \\",
    "    [--pubchem-throttle 1.1] [--kegg-throttle 0.5] \\",
    "    [--request-timeout 60] [--research-enrichment-limit 1000] \\",
    "    [--full-enrichment-limit 250] \\",
    "    [--release-manifest uafR_release_manifest.json] \\",
    "    [--source-tarball uafR_<version>.tar.gz] \\",
    "    [--require-release-artifact true] \\",
    "    [--project-id plant_chemistry_panel] \\",
    "    [--server-results-dir downloaded_server_output] \\",
    "    [--resume true] [--overwrite false] [--progress true]",
    "",
    "Credentials:",
    "  Set NCBI_EMAIL, NCBI_TOOL, and optional NCBI_API_KEY in the",
    "  environment. API keys are intentionally not accepted as CLI arguments",
    "  and are never written to manifests.",
    "",
    "Exit status 75 means a public service requested a pause. Re-run the same",
    "command later to resume from validated checkpoints and caches.",
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
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token,
                                        call. = FALSE)
    raw = substring(token, 3L)
    if (grepl("=", raw, fixed = TRUE)) {
      parts = strsplit(raw, "=", fixed = TRUE)[[1]]
      key = parts[[1]]
      value = paste(parts[-1], collapse = "=")
      i = i + 1L
    } else {
      key = raw
      if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
        value = "true"
        i = i + 1L
      } else {
        value = args[[i + 1L]]
        i = i + 2L
      }
    }
    key = gsub("-", "_", key, fixed = TRUE)
    out[[key]] = value
  }
  out
}

load_config = function(path) {
  if (is.null(path)) return(list())
  if (!file.exists(path)) stop("Config file does not exist: ", path,
                               call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite package is required to read --config.", call. = FALSE)
  }
  x = jsonlite::read_json(path, simplifyVector = TRUE)
  if (!is.list(x)) stop("Config JSON must contain one object.", call. = FALSE)
  x
}

as_flag = function(x, name) {
  value = tolower(trimws(as.character(x[[1]])))
  if (!value %in% c("true", "false", "yes", "no", "1", "0")) {
    stop("--", gsub("_", "-", name, fixed = TRUE),
         " must be true or false.", call. = FALSE)
  }
  value %in% c("true", "yes", "1")
}

as_number = function(x, name, integer = FALSE) {
  value = suppressWarnings(as.numeric(x[[1]]))
  if (is.na(value) && tolower(as.character(x[[1]])) %in% c("inf", "infinity")) {
    value = Inf
  }
  if (is.na(value)) stop("--", gsub("_", "-", name, fixed = TRUE),
                         " must be numeric.", call. = FALSE)
  if (isTRUE(integer) && is.finite(value)) as.integer(value) else value
}

as_vector = function(x) {
  if (length(x) > 1L) return(as.character(x))
  value = trimws(as.character(x[[1]]))
  if (!nzchar(value)) return(character())
  trimws(unlist(strsplit(value, "[,;]", perl = TRUE), use.names = FALSE))
}

load_uafr = function() {
  dev_root = Sys.getenv("UAFR_DEV_ROOT", "")
  if (nzchar(dev_root) && file.exists(file.path(dev_root, "DESCRIPTION")) &&
      requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(dev_root, quiet = TRUE)
    return(invisible(TRUE))
  }
  if (!requireNamespace("uafR", quietly = TRUE)) {
    stop("uafR is not installed. Install the production source package first.",
         call. = FALSE)
  }
  suppressPackageStartupMessages(library(uafR))
  invisible(TRUE)
}

normalize_values = function(x) {
  vector_fields = c("sources", "taxon_fallback", "exclusion_files",
                    "supporting_files", "pilot_species")
  logical_fields = c("download_npass", "require_all_providers", "resume",
                     "overwrite", "progress", "require_release_artifact")
  integer_fields = c("expected_species_count", "expected_exclusion_count",
                     "species_chunk_size", "compound_batch_size",
                     "pilot_count", "pilot_enrichment_compounds",
                     "max_pubmed_records", "full_enrichment_limit")
  numeric_fields = c("max_provider_records", "provider_throttle",
                     "knapsack_throttle", "pubchem_throttle",
                     "kegg_throttle", "request_timeout",
                     "research_enrichment_limit")
  for (name in intersect(names(x), vector_fields)) {
    x[[name]] = as_vector(x[[name]])
  }
  for (name in intersect(names(x), logical_fields)) {
    if (!is.logical(x[[name]])) x[[name]] = as_flag(x[[name]], name)
  }
  for (name in intersect(names(x), integer_fields)) {
    if (!is.numeric(x[[name]])) x[[name]] = as_number(x[[name]], name, TRUE)
  }
  for (name in intersect(names(x), numeric_fields)) {
    if (!is.numeric(x[[name]])) x[[name]] = as_number(x[[name]], name)
  }
  x
}

main = function() {
  cli = parse_args(commandArgs(trailingOnly = TRUE))
  config = load_config(cli$config)
  cli$config = NULL
  args = modifyList(config, cli)
  required = c("mode", "plant_csv", "out_dir")
  missing = required[!required %in% names(args) |
                       vapply(required, function(name) {
                         is.null(args[[name]]) || !nzchar(as.character(args[[name]][[1]]))
                       }, logical(1))]
  if (length(missing) > 0) {
    usage()
    stop("Missing required argument(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if ("ncbi_api_key" %in% names(args)) {
    stop("NCBI API keys must be supplied through NCBI_API_KEY, not CLI/config.",
         call. = FALSE)
  }
  args = normalize_values(args)
  allowed = names(formals(uafR::runPlantChemistryPanel))
  unknown = setdiff(names(args), allowed)
  if (length(unknown) > 0) {
    stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  result = do.call(uafR::runPlantChemistryPanel, args)
  cat("uafR plant chemistry panel state:", result$Status$state, "\n")
  cat("Stage:", result$Status$stage, "\n")
  cat("Output:", result$OutputPath, "\n")
  as.integer(result$ExitStatus)
}

load_uafr()
status = tryCatch(main(), error = function(error) {
  message("Error: ", conditionMessage(error))
  1L
})
quit(save = "no", status = status)
