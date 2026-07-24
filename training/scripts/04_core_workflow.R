# Core uafR workflow using bundled example data.

run_core_workflow <- function(query_chemicals = c(
                                "Linalool",
                                "Methyl Salicylate",
                                "Limonene",
                                "alpha-Thujene"
                              ),
                              live_lookup = FALSE,
                              output_dir = NULL) {
  if (!requireNamespace("uafR", quietly = TRUE)) {
    stop("uafR is not installed.", call. = FALSE)
  }

  library(uafR)
  data("standard_data", package = "uafR")
  data("standard_spread", package = "uafR")
  data("standard_exacto", package = "uafR")

  if (isTRUE(live_lookup)) {
    spread <- spreadOut(standard_data)
    exact <- mzExacto(spread, query_chemicals)
  } else {
    spread <- standard_spread
    exact <- standard_exacto
  }

  result <- list(
    query_chemicals = query_chemicals,
    live_lookup = isTRUE(live_lookup),
    spread = spread,
    exact = exact
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(exact, file.path(output_dir, "core_exact_matches.csv"),
                     row.names = FALSE)
    saveRDS(result, file.path(output_dir, "core_workflow_result.rds"))
  }

  result
}
