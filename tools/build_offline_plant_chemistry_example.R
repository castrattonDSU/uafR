#!/usr/bin/env Rscript

message("Building offline uafR plant chemistry example bundle...")

load_uafr = function() {
  if (requireNamespace("uafR", quietly = TRUE)) {
    suppressPackageStartupMessages(library(uafR))
    return(invisible(TRUE))
  }
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
    return(invisible(TRUE))
  }
  stop("Install uafR or run this script from the uafR repository with devtools available.",
       call. = FALSE)
}

load_uafr()

repo_example_dir = file.path("inst", "extdata", "offline_plant_chemistry")
installed_example_dir = system.file("extdata", "offline_plant_chemistry",
                                    package = "uafR")
example_dir = if (dir.exists(repo_example_dir)) {
  repo_example_dir
} else {
  installed_example_dir
}
if (!dir.exists(example_dir)) {
  stop("Offline example data directory was not found.", call. = FALSE)
}

parse_args = function(args) {
  out = list()
  i = 1L
  while (i <= length(args)) {
    token = args[[i]]
    if (token %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript tools/build_offline_plant_chemistry_example.R ",
        "[--out-dir PATH]\n",
        sep = ""
      )
      quit(save = "no", status = 0L)
    }
    if (identical(token, "--out-dir")) {
      if (i == length(args)) {
        stop("`--out-dir` requires a path.", call. = FALSE)
      }
      out$out_dir = args[[i + 1L]]
      i = i + 2L
      next
    }
    if (startsWith(token, "--out-dir=")) {
      out$out_dir = substring(token, nchar("--out-dir=") + 1L)
      i = i + 1L
      next
    }
    if (!startsWith(token, "-") && is.null(out$out_dir)) {
      out$out_dir = token
      i = i + 1L
      next
    }
    stop("Unknown argument: ", token, call. = FALSE)
  }
  out
}

args = parse_args(commandArgs(trailingOnly = TRUE))
out_dir = if (!is.null(args$out_dir) && nzchar(args$out_dir)) {
  args$out_dir
} else {
  file.path(tempdir(), "uafR_offline_plant_chemistry_example")
}

plants = file.path(example_dir, "plant_list.csv")
metadata = file.path(example_dir, "plant_metadata.csv")
plant_compounds = file.path(example_dir, "plant_compounds.csv")
pair_tanimoto = file.path(example_dir, "plant_compound_pair_tanimoto.csv")

manifest = runPlantChemistryProject(
  plant_list = plants,
  metadata = metadata,
  plant_compounds = plant_compounds,
  plant_compound_pair_tanimoto = pair_tanimoto,
  output_dir = out_dir,
  project_id = "offline_simulated_plant_chemistry_example",
  overwrite = TRUE
)

validation = validatePlantChemistryAnalysisBundle(
  manifest$Project$bundle_dir[[1]]
)
if (!identical(validation$Summary$ExportReadyStatus[[1]], "pass")) {
  stop("Offline example bundle validation did not pass.", call. = FALSE)
}

message("Offline example bundle written to: ", out_dir)
message("Bundle validation: ", validation$Summary$ExportReadyStatus[[1]])
