# Small research-grade categorate smoke test.

run_categorate_smoke <- function(compounds = c("aspirin", "caffeine"),
                                 cache = TRUE,
                                 throttle = 0.2,
                                 output_dir = NULL) {
  if (!requireNamespace("uafR", quietly = TRUE)) {
    stop("uafR is not installed.", call. = FALSE)
  }

  library(uafR)
  data("library_data", package = "uafR")

  result <- categorate(
    compounds = compounds,
    chemical_library = library_data,
    input_format = "wide",
    detail = "research",
    cache = cache,
    throttle = throttle,
    assay_detail_limit = 0,
    trait_matrix_profile = "core",
    trait_matrix_min_confidence = "medium",
    trait_matrix_max_traits = 100
  )

  validation <- validateCategorateResult(result)
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
