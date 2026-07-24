# Offline smoke tests for the public training support scripts.

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop("Run this script from the uafR repository root.", call. = FALSE)
}

if (!requireNamespace("uafR", quietly = TRUE)) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    message("Loading uafR from the source tree for repository smoke tests")
    pkgload::load_all(repo_root, quiet = TRUE)
  } else {
    stop(
      "uafR is not installed and pkgload is unavailable. Install uafR or ",
      "install the pkgload package before running repository smoke tests.",
      call. = FALSE
    )
  }
}

script_files <- c(
  "00_install_check.R",
  "01_project_setup.R",
  "04_core_workflow.R",
  "05_categorate_research.R",
  "06_trait_matrix_analysis.R",
  "07_visualization.R",
  "09_reproducibility_check.R",
  "10_species_phytochemistry.R"
)

message("Sourcing training scripts in an isolated environment")
training_env <- new.env(parent = globalenv())
for (script in script_files) {
  sys.source(
    file.path(repo_root, "training", "scripts", script),
    envir = training_env
  )
}

message("Checking project setup in a temporary directory")
tmp_project <- file.path(tempdir(), paste0("uafR-training-", Sys.getpid()))
training_env$create_training_project(tmp_project)
project_report <- training_env$check_training_project(tmp_project)
stopifnot(all(project_report$Directories$Exists))
stopifnot(all(project_report$Files$Exists))

run_live <- identical(Sys.getenv("UAFR_TRAINING_LIVE"), "1")

if (requireNamespace("uafR", quietly = TRUE)) {
  message("Running offline core workflow with package example data")
  core <- training_env$run_core_workflow(live_lookup = FALSE)
  stopifnot(is.list(core))
  stopifnot("exact" %in% names(core))

  message("Running the categorate teaching wrapper with an offline fixture")
  mock_result <- list(
    SourceCoverage = data.frame(
      Query = c("aspirin", "caffeine"),
      PubChem = c(TRUE, TRUE),
      stringsAsFactors = FALSE
    ),
    ChemicalTraitReport = data.frame(
      Query = c("aspirin", "caffeine"),
      TraitCount = c(2L, 3L),
      stringsAsFactors = FALSE
    ),
    ChemicalTraitMatrix = data.frame(
      Query = c("aspirin", "caffeine"),
      trait_example = c(1L, 0L),
      stringsAsFactors = FALSE
    )
  )
  mock_categorate <- function(...) mock_result
  mock_validate <- function(x) {
    list(Summary = data.frame(Status = "PASS", stringsAsFactors = FALSE))
  }
  categorate_output <- file.path(tempdir(), "uafR-training-categorate")
  categorate_smoke <- training_env$run_categorate_smoke(
    compounds = c("aspirin", "caffeine"),
    output_dir = categorate_output,
    categorate_fun = mock_categorate,
    validation_fun = mock_validate,
    chemical_library = data.frame(
      Chemical = c("aspirin", "caffeine"),
      stringsAsFactors = FALSE
    )
  )
  stopifnot(is.list(categorate_smoke))
  stopifnot(file.exists(file.path(
    categorate_output,
    "chemical_trait_report.csv"
  )))

  message("Running the offline species-first phytochemistry extension")
  plant_output <- file.path(tempdir(), "uafR-training-plant")
  plant_smoke <- training_env$run_species_phytochemistry_smoke(
    output_dir = plant_output,
    verbose = FALSE
  )
  stopifnot(inherits(plant_smoke, "uaf_plant_phytochemistry"))
  plant_manifest_files <- list.files(
    file.path(plant_output, "plant_phytochemistry_export"),
    pattern = "ExportManifest[.]csv$",
    full.names = TRUE
  )
  stopifnot(length(plant_manifest_files) == 1L)
} else {
  stop(
    "uafR is not installed; the public training smoke test requires the ",
    "current package.",
    call. = FALSE
  )
}

if (run_live && requireNamespace("uafR", quietly = TRUE)) {
  message("Running service-dependent core workflow with package example data")
  core <- training_env$run_core_workflow(live_lookup = TRUE)
  stopifnot(is.list(core))
  stopifnot("exact" %in% names(core))

  message("Running the two-compound live categorate teaching workflow")
  live_categorate <- training_env$run_categorate_smoke(
    compounds = c("aspirin", "caffeine"),
    cache = TRUE,
    throttle = 0.2
  )
  stopifnot(is.list(live_categorate))
} else if (run_live) {
  message("uafR is not installed; skipping service-dependent workflow smoke test")
} else {
  message("Skipping service-dependent uafR workflow smoke test")
  message("Set UAFR_TRAINING_LIVE=1 to include it intentionally")
}

message("Training smoke tests completed")
