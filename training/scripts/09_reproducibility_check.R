# Lightweight reproducibility audit for a student training project.

check_training_project <- function(path = ".") {
  required_dirs <- c(
    "analysis",
    "data_raw",
    "data_processed",
    "figures",
    "results",
    "scripts",
    "reports",
    "logs"
  )

  dir_status <- data.frame(
    Directory = required_dirs,
    Exists = dir.exists(file.path(path, required_dirs)),
    stringsAsFactors = FALSE
  )

  file_status <- data.frame(
    File = c("README.md", file.path("logs", "research_log.md")),
    Exists = file.exists(file.path(path, c("README.md", file.path("logs", "research_log.md")))),
    stringsAsFactors = FALSE
  )

  script_files <- list.files(file.path(path, "scripts"), pattern = "\\.R$",
                             full.names = TRUE)
  result_files <- list.files(file.path(path, "results"), full.names = TRUE)
  figure_files <- list.files(file.path(path, "figures"), full.names = TRUE)

  package_status <- data.frame(
    Package = c("uafR", "ChemmineR", "fmcsR", "jsonlite", "webchem"),
    Installed = vapply(c("uafR", "ChemmineR", "fmcsR", "jsonlite", "webchem"),
                       requireNamespace, logical(1), quietly = TRUE),
    stringsAsFactors = FALSE
  )

  report <- list(
    Path = normalizePath(path, mustWork = FALSE),
    Directories = dir_status,
    Files = file_status,
    ScriptCount = length(script_files),
    ResultFileCount = length(result_files),
    FigureFileCount = length(figure_files),
    Packages = package_status,
    SessionInfo = utils::sessionInfo()
  )

  print(dir_status, row.names = FALSE)
  print(file_status, row.names = FALSE)
  print(package_status, row.names = FALSE)
  message("Scripts: ", length(script_files))
  message("Result files: ", length(result_files))
  message("Figure files: ", length(figure_files))

  invisible(report)
}
