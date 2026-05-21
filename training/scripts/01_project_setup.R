# Create a reproducible project folder for DSU dsDNA Core uafR training.

create_training_project <- function(path = "student-project", overwrite = FALSE) {
  if (dir.exists(path) && !isTRUE(overwrite)) {
    stop("Project directory already exists: ", path, call. = FALSE)
  }

  dirs <- c(
    "analysis",
    "data_raw",
    "data_processed",
    "figures",
    "results",
    "scripts",
    "reports",
    "logs"
  )

  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  for (dir in dirs) {
    dir.create(file.path(path, dir), recursive = TRUE, showWarnings = FALSE)
  }

  readme <- c(
    "# Student Research Project",
    "",
    "## Research question",
    "",
    "Write the current research question here.",
    "",
    "## Data",
    "",
    "- `data_raw/`: unchanged input files or secure-access instructions.",
    "- `data_processed/`: derived tables created by scripts.",
    "",
    "## Scripts",
    "",
    "Number scripts in the order they should be run.",
    "",
    "## Reproducibility",
    "",
    "Record R version, uafR version, and package versions in `logs/`."
  )

  log <- c(
    "# Research Log",
    "",
    paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "## Decisions",
    "",
    "## Problems and fixes",
    "",
    "## Interpretation notes"
  )

  writeLines(readme, file.path(path, "README.md"))
  writeLines(log, file.path(path, "logs", "research_log.md"))

  message("Created project at: ", normalizePath(path, mustWork = FALSE))
  invisible(normalizePath(path, mustWork = FALSE))
}
