#!/usr/bin/env Rscript

file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path = if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[[1]])
} else {
  file.path(getwd(), "tools", "run_plant_chemistry_project.R")
}
repo_root = normalizePath(file.path(dirname(script_path), ".."),
                          winslash = "/", mustWork = FALSE)

load_uafr = function(repo_root) {
  if (file.exists(file.path(repo_root, "DESCRIPTION")) &&
      requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(repo_root, quiet = TRUE)
    return(invisible(TRUE))
  }
  if (requireNamespace("uafR", quietly = TRUE)) {
    suppressPackageStartupMessages(library(uafR))
    return(invisible(TRUE))
  }
  stop("Could not load uafR. Install the package or run this script from the ",
       "uafR repository with devtools installed.", call. = FALSE)
}

parse_args = function(args) {
  out = list()
  i = 1L
  while (i <= length(args)) {
    arg = args[[i]]
    if (!grepl("^--", arg)) stop("Unexpected argument: ", arg,
                                  call. = FALSE)
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
    out[[gsub("-", "_", key)]] = value
    i = i + 1L
  }
  out
}

optional_arg = function(args, name, default = NULL) {
  value = args[[name]]
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) default else
    value
}

flag_arg = function(args, name, default = FALSE) {
  value = args[[name]]
  if (is.null(value)) return(default)
  if (identical(value, TRUE)) return(TRUE)
  tolower(value) %in% c("true", "t", "yes", "y", "1")
}

usage = function() {
  paste(
    "Usage:",
    "Rscript tools/run_plant_chemistry_project.R \\",
    "  --plant-list plants.csv \\",
    "  --out-dir results/plant_chemistry_project \\",
    "  [--metadata plant_metadata.csv] \\",
    "  [--plant-compounds curated_plant_compounds.csv] \\",
    "  [--plant-result plant_phytochemistry_result.rds] \\",
    "  [--categorate-batches categorate_batches] \\",
    "  [--species-pair-tanimoto species_pair_tanimoto.csv] \\",
    "  [--resolved-compounds resolved_compounds.csv] \\",
    "  [--pubchem-fingerprints pubchem_fingerprints.csv] \\",
    "  [--plant-compound-pair-tanimoto plant_compound_pair_tanimoto.csv.gz] \\",
    "  [--project-id project_label] [--run-discovery true] [--overwrite true]",
    sep = "\n"
  )
}

args = parse_args(commandArgs(TRUE))
if (isTRUE(args$help) || is.null(args$plant_list) || is.null(args$out_dir)) {
  message(usage())
  if (isTRUE(args$help)) quit(status = 0)
  quit(status = 1)
}

load_uafr(repo_root)

manifest = runPlantChemistryProject(
  plant_list = args$plant_list,
  output_dir = args$out_dir,
  metadata = optional_arg(args, "metadata"),
  plant_compounds = optional_arg(args, "plant_compounds"),
  plant_result = optional_arg(args, "plant_result"),
  categorate_batches = optional_arg(args, "categorate_batches"),
  species_pair_tanimoto = optional_arg(args, "species_pair_tanimoto"),
  resolved_compounds = optional_arg(args, "resolved_compounds"),
  pubchem_fingerprints = optional_arg(args, "pubchem_fingerprints"),
  plant_compound_pair_tanimoto =
    optional_arg(args, "plant_compound_pair_tanimoto"),
  project_id = optional_arg(args, "project_id"),
  run_discovery = flag_arg(args, "run_discovery", FALSE),
  overwrite = flag_arg(args, "overwrite", FALSE)
)

message("Plant chemistry project complete.")
message("  output: ", manifest$Project$output_dir[[1]])
message("  bundle: ", manifest$Project$bundle_dir[[1]])
message("  validation: ", manifest$Project$bundle_validation_status[[1]])
