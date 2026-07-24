# Build analysis-ready trait matrices and evidence tables from a categorate result.

build_trait_analysis <- function(result,
                                 output_dir = NULL,
                                 min_confidence = "medium",
                                 max_traits = 200) {
  if (!requireNamespace("uafR", quietly = TRUE)) {
    stop("uafR is not installed.", call. = FALSE)
  }

  library(uafR)

  validation <- validateCategorateResult(result)
  print(validation$Summary)

  core_matrix <- chemicalTraitMatrix(
    result,
    profile = "core",
    mode = "binary",
    min_confidence = min_confidence,
    max_traits = max_traits
  )

  confidence_matrix <- chemicalTraitMatrix(
    result,
    profile = "full",
    mode = "confidence",
    min_confidence = min_confidence,
    max_traits = max_traits
  )

  ontology_matrix <- chemicalTraitOntologyMatrix(
    result,
    mode = "binary",
    min_confidence = min_confidence,
    max_terms = max_traits
  )

  evidence <- chemicalTraitEvidence(result)
  report <- chemicalTraitReport(result)

  out <- list(
    validation = validation,
    core_matrix = core_matrix,
    confidence_matrix = confidence_matrix,
    ontology_matrix = ontology_matrix,
    evidence = evidence,
    report = report
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(core_matrix, file.path(output_dir, "core_trait_matrix.csv"),
                     row.names = FALSE)
    utils::write.csv(confidence_matrix,
                     file.path(output_dir, "confidence_trait_matrix.csv"),
                     row.names = FALSE)
    utils::write.csv(ontology_matrix,
                     file.path(output_dir, "ontology_trait_matrix.csv"),
                     row.names = FALSE)
    utils::write.csv(evidence, file.path(output_dir, "trait_evidence.csv"),
                     row.names = FALSE)
    utils::write.csv(report, file.path(output_dir, "trait_report.csv"),
                     row.names = FALSE)
  }

  out
}
