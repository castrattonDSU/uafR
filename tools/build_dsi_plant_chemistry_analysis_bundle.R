#!/usr/bin/env Rscript

file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path = if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[[1]])
} else {
  file.path(getwd(), "tools", "build_dsi_plant_chemistry_analysis_bundle.R")
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
  stop("Could not load uafR. Install the current package or run this script ",
       "from the uafR repository with devtools installed.", call. = FALSE)
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

optional_arg = function(args, name, default) {
  value = args[[name]]
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) default else
    value
}

args = parse_args(commandArgs(TRUE))
load_uafr(repo_root)

categorate_batch_dir = optional_arg(
  args,
  "categorate_batch_dir",
  file.path(repo_root, "lotus_cache", "exports",
            "dsi_categorate_research_pubchem_cid_20260611",
            "categorate_batches")
)
tanimoto_dir = optional_arg(
  args,
  "tanimoto_dir",
  file.path(repo_root, "lotus_cache", "exports",
            "dsi_categorate_tanimoto_20260611")
)
plant_list = optional_arg(
  args,
  "plant_list",
  "/Users/chasestratton/Desktop/DSU_Hornets/Projects/DSI_Dataset/submissions/final_species_round2_2026.csv"
)
out_dir = optional_arg(
  args,
  "out_dir",
  file.path(repo_root, "lotus_cache", "exports",
            "dsi_plant_chemistry_analysis_bundle_20260624_publication_ready")
)

message("Building DSI plant chemistry analysis bundle:")
message("  categorate batches: ", categorate_batch_dir)
message("  Tanimoto directory: ", tanimoto_dir)
message("  plant list: ", plant_list)
message("  output: ", out_dir)

manifest = exportPlantChemistryAnalysisBundle(
  categorate_batches = categorate_batch_dir,
  path = out_dir,
  plant_membership = file.path(tanimoto_dir,
                               "dsi_species_compound_membership.csv"),
  species_pair_tanimoto = file.path(tanimoto_dir,
                                    "dsi_species_pair_tanimoto_summary.csv"),
  resolved_compounds = file.path(tanimoto_dir,
                                 "dsi_resolved_compounds_for_categorate.csv"),
  pubchem_fingerprints = file.path(tanimoto_dir,
                                   "dsi_pubchem_fingerprints.csv"),
  file_references = c(
    plant_compound_pairs = file.path(tanimoto_dir,
                                     "dsi_plant_compound_pair_tanimoto.csv.gz"),
    compound_pairs = file.path(tanimoto_dir,
                               "dsi_compound_pair_tanimoto.csv.gz"),
    categorate_batches = categorate_batch_dir
  ),
  plant_list = plant_list,
  project_id = "DSI_Dataset",
  tables = c("ChemicalTraitSummary", "DerivedGroups", "PubChemProperties",
             "SourceCoverage", "ValidationSummary", "ValidationIssues"),
  format = "csv",
  overwrite = TRUE,
  finalize = TRUE,
  validate_export = TRUE
)

validation = validatePlantChemistryAnalysisBundle(
  out_dir,
  use_python = TRUE,
  use_pandas = TRUE
)

message("Bundle written: ", out_dir)
message("Export readiness: ", validation$Summary$ExportReadyStatus[[1]])
print(manifest[, c("Table", "FileName", "RowCount", "ColumnCount")],
      row.names = FALSE)

if (!identical(validation$Summary$ExportReadyStatus[[1]], "pass")) {
  print(validation$CSVValidation[validation$CSVValidation$Status != "pass", ],
        row.names = FALSE)
  print(validation$RequiredColumns[validation$RequiredColumns$Status == "fail", ],
        row.names = FALSE)
  stop("Bundle validation did not pass.", call. = FALSE)
}
