# Small research-grade categorate smoke test.

run_categorate_smoke <- function(compounds = c("aspirin", "caffeine"),
                                 cache = TRUE,
                                 throttle = 0.2,
                                 output_dir = NULL,
                                 categorate_fun = NULL,
                                 validation_fun = NULL,
                                 chemical_library = NULL) {
  if (is.null(categorate_fun) || is.null(validation_fun) ||
      is.null(chemical_library)) {
    if (!requireNamespace("uafR", quietly = TRUE)) {
      stop("uafR is not installed.", call. = FALSE)
    }
  }

  if (is.null(categorate_fun)) {
    categorate_fun <- uafR::categorate
  }
  if (is.null(validation_fun)) {
    validation_fun <- uafR::validateCategorateResult
  }
  if (is.null(chemical_library)) {
    data_env <- new.env(parent = emptyenv())
    utils::data("library_data", package = "uafR", envir = data_env)
    chemical_library <- data_env$library_data
  }

  result <- categorate_fun(
    compounds = compounds,
    chemical_library = chemical_library,
    input_format = "wide",
    detail = "research",
    cache = cache,
    throttle = throttle,
    assay_detail_limit = 0,
    trait_matrix_profile = "core",
    trait_matrix_min_confidence = "medium",
    trait_matrix_max_traits = 100
  )

  required_tables <- c(
    "SourceCoverage",
    "ChemicalTraitReport",
    "ChemicalTraitMatrix"
  )
  missing_tables <- setdiff(required_tables, names(result))
  if (length(missing_tables) > 0) {
    stop(
      "categorate smoke result is missing required table(s): ",
      paste(missing_tables, collapse = ", "),
      call. = FALSE
    )
  }

  validation <- validation_fun(result)
  if (!is.list(validation) || is.null(validation$Summary)) {
    stop("Validation did not return a Summary table.", call. = FALSE)
  }
  print(validation$Summary)
  print(result$SourceCoverage)
  print(result$ChemicalTraitReport)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, file.path(output_dir, "categorate_research_smoke.rds"))
    utils::write.csv(result$ChemicalTraitReport,
                     file.path(output_dir, "chemical_trait_report.csv"),
                     row.names = FALSE)
    utils::write.csv(result$ChemicalTraitMatrix,
                     file.path(output_dir, "chemical_trait_matrix.csv"),
                     row.names = FALSE)
  }

  result
}
