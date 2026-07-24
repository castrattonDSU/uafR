# Base R visualizations for uafR training outputs.

as_numeric_matrix <- function(data, id_columns = c("Query", "Chemical", "Compound")) {
  if (!is.data.frame(data) || nrow(data) < 1) {
    return(matrix(numeric(), nrow = 0, ncol = 0))
  }

  value_cols <- setdiff(names(data), id_columns)
  if (length(value_cols) < 1) {
    return(matrix(numeric(), nrow = 0, ncol = 0))
  }

  mat <- data[value_cols]
  mat[] <- lapply(mat, function(x) {
    if (is.logical(x)) return(as.integer(x))
    if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
    suppressWarnings(as.numeric(x))
  })
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  mat
}

plot_matrix_heatmap <- function(mat, title, file) {
  if (!is.matrix(mat) || nrow(mat) < 1 || ncol(mat) < 1) {
    return(NA_character_)
  }

  grDevices::png(file, width = 1400, height = 900, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(mar = c(8, 8, 4, 2))
  graphics::image(
    x = seq_len(nrow(mat)),
    y = seq_len(ncol(mat)),
    z = mat,
    col = grDevices::colorRampPalette(c("white", "#C69214", "#7A0019"))(20),
    axes = FALSE,
    main = title,
    xlab = "Chemical",
    ylab = "Trait or source"
  )
  graphics::axis(1, at = seq_len(nrow(mat)), labels = rownames(mat), las = 2,
                 cex.axis = 0.7)
  graphics::axis(2, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2,
                 cex.axis = 0.7)
  file
}

build_training_figures <- function(result, output_dir = "figures") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- list()

  if (is.data.frame(result$SourceCoverage) && nrow(result$SourceCoverage) > 0) {
    coverage <- result$SourceCoverage
    rownames(coverage) <- if ("Query" %in% names(coverage)) {
      coverage$Query
    } else {
      seq_len(nrow(coverage))
    }
    coverage_mat <- as_numeric_matrix(coverage)
    paths$source_coverage <- plot_matrix_heatmap(
      coverage_mat,
      "Source coverage by chemical",
      file.path(output_dir, "source_coverage_heatmap.png")
    )
  }

  if (is.data.frame(result$ChemicalTraitMatrix) &&
      nrow(result$ChemicalTraitMatrix) > 0) {
    trait_matrix <- result$ChemicalTraitMatrix
    rownames(trait_matrix) <- if ("Query" %in% names(trait_matrix)) {
      trait_matrix$Query
    } else {
      seq_len(nrow(trait_matrix))
    }
    trait_mat <- as_numeric_matrix(trait_matrix)
    paths$trait_matrix <- plot_matrix_heatmap(
      trait_mat,
      "Chemical trait matrix",
      file.path(output_dir, "chemical_trait_matrix_heatmap.png")
    )

    trait_counts <- rowSums(trait_mat > 0)
    count_file <- file.path(output_dir, "trait_counts_by_chemical.png")
    grDevices::png(count_file, width = 1200, height = 800, res = 160)
    graphics::par(mar = c(8, 5, 4, 1))
    graphics::barplot(
      trait_counts,
      las = 2,
      col = "#7A0019",
      border = NA,
      main = "Trait count by chemical",
      ylab = "Number of traits"
    )
    grDevices::dev.off()
    paths$trait_counts <- count_file
  }

  invisible(paths)
}
