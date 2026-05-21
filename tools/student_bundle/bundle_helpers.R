# Shared helper functions for the uafR private student bundle.

uafr_bundle_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    path <- sub("^--file=", "", file_arg[[length(file_arg)]])
    return(dirname(normalizePath(path, winslash = "/", mustWork = TRUE)))
  }

  frames <- sys.frames()
  for (frame in rev(frames)) {
    if (!is.null(frame$ofile)) {
      return(dirname(normalizePath(frame$ofile, winslash = "/", mustWork = TRUE)))
    }
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(path)) {
      return(dirname(normalizePath(path, winslash = "/", mustWork = TRUE)))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

uafr_bundle_root <- function(start = uafr_bundle_script_dir()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "README_STUDENT_INSTALL.md")) &&
        dir.exists(file.path(current, "packages"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find the uafR bundle root. Make sure the whole bundle folder was unzipped.",
           call. = FALSE)
    }
    current <- parent
  }
}

uafr_log_file <- function(bundle_dir, name) {
  file.path(bundle_dir, paste0(name, ".txt"))
}

uafr_note <- function(..., log_file = NULL) {
  text <- paste0(...)
  cat(text, "\n", sep = "")
  if (!is.null(log_file)) {
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", text, "\n",
        file = log_file, append = TRUE, sep = "")
  }
  invisible(text)
}

uafr_header <- function(title, log_file = NULL) {
  if (!is.null(log_file)) {
    cat("", file = log_file)
  }
  line <- paste(rep("=", nchar(title)), collapse = "")
  uafr_note(line, log_file = log_file)
  uafr_note(title, log_file = log_file)
  uafr_note(line, log_file = log_file)
  uafr_note("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
            log_file = log_file)
  uafr_note("R: ", R.version.string, log_file = log_file)
  uafr_note("Platform: ", R.version$platform, log_file = log_file)
  invisible(TRUE)
}

uafr_step <- function(number, total, title, log_file = NULL) {
  uafr_note("", log_file = log_file)
  uafr_note("Step ", number, " of ", total, ": ", title, log_file = log_file)
  invisible(TRUE)
}

uafr_find_package_archive <- function(bundle_dir) {
  packages_dir <- file.path(bundle_dir, "packages")
  tarballs <- list.files(packages_dir, pattern = "^uafR_[0-9].*[.]tar[.]gz$",
                         full.names = TRUE)
  if (length(tarballs) == 0) {
    stop("No uafR package archive was found in ", packages_dir,
         ". Make sure the whole bundle folder was unzipped.", call. = FALSE)
  }
  if (length(tarballs) > 1) {
    tarballs <- tarballs[order(file.info(tarballs)$mtime, decreasing = TRUE)]
  }
  normalizePath(tarballs[[1]], winslash = "/", mustWork = TRUE)
}

uafr_package_status <- function(packages) {
  data.frame(
    Package = packages,
    Installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    Version = vapply(packages, function(pkg) {
      if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
      as.character(utils::packageVersion(pkg))
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

uafr_install_missing_cran <- function(packages, log_file = NULL) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) {
    uafr_note("CRAN dependencies already installed: ",
              paste(packages, collapse = ", "), log_file = log_file)
    return(invisible(TRUE))
  }

  uafr_note("Installing CRAN package(s): ", paste(missing, collapse = ", "),
            log_file = log_file)
  utils::install.packages(missing, repos = "https://cloud.r-project.org")
  invisible(TRUE)
}

uafr_require_installed <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Required package(s) are still missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

uafr_library_writable <- function() {
  lib <- .libPaths()[1]
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  test_file <- tempfile("uafr-write-test-", tmpdir = lib)
  ok <- tryCatch({
    file.create(test_file)
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (file.exists(test_file)) {
    unlink(test_file)
  }
  isTRUE(ok)
}

uafr_download_ok <- function(url, timeout = 10) {
  dest <- tempfile("uafr-download-test-")
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  on.exit(if (file.exists(dest)) unlink(dest), add = TRUE)
  options(timeout = timeout)
  result <- tryCatch(
    utils::download.file(url, destfile = dest, quiet = TRUE, mode = "wb"),
    error = function(e) 1L,
    warning = function(w) 1L
  )
  identical(result, 0L) && file.exists(dest) && file.info(dest)$size > 0
}

uafr_status_row <- function(check, status, details) {
  data.frame(Check = check, Status = status, Details = details,
             stringsAsFactors = FALSE)
}

uafr_print_status <- function(status_table, log_file = NULL) {
  print(status_table, row.names = FALSE)
  if (!is.null(log_file)) {
    capture.output(print(status_table, row.names = FALSE),
                   file = log_file, append = TRUE)
  }
  invisible(status_table)
}

uafr_run_verification <- function(bundle_dir, log_file = NULL) {
  verify_script <- file.path(bundle_dir, "verify_uafR_install.R")
  if (!file.exists(verify_script)) {
    warning("Verification script was not found: ", verify_script)
    return(invisible(FALSE))
  }
  uafr_note("Running verification script: ", basename(verify_script),
            log_file = log_file)
  source(verify_script, local = new.env(parent = globalenv()))
  invisible(TRUE)
}

uafr_install_from_bundle <- function(bundle_dir,
                                     log_file = NULL,
                                     action = c("install", "update")) {
  action <- match.arg(action)
  options(repos = c(CRAN = "https://cloud.r-project.org"))

  archive <- uafr_find_package_archive(bundle_dir)

  uafr_step(1, 5, "Checking the bundle and R library", log_file = log_file)
  uafr_note("Bundle folder: ", bundle_dir, log_file = log_file)
  uafr_note("Package archive: ", basename(archive), log_file = log_file)
  uafr_note("R library: ", .libPaths()[1], log_file = log_file)
  if (!uafr_library_writable()) {
    stop("The first R library is not writable: ", .libPaths()[1],
         ". In RStudio, restart R and try again. If the problem remains, ask for help setting a user library.",
         call. = FALSE)
  }

  uafr_step(2, 5, "Installing CRAN dependencies", log_file = log_file)
  uafr_install_missing_cran(c("jsonlite", "webchem"), log_file = log_file)

  uafr_step(3, 5, "Installing Bioconductor dependencies", log_file = log_file)
  bioc_packages <- c("ChemmineR", "fmcsR")
  missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace,
                                        logical(1), quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    uafr_install_missing_cran("BiocManager", log_file = log_file)
    uafr_require_installed("BiocManager")
    uafr_note("Installing Bioconductor package(s): ",
              paste(missing_bioc, collapse = ", "), log_file = log_file)
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  } else {
    uafr_note("Bioconductor dependencies already installed: ",
              paste(bioc_packages, collapse = ", "), log_file = log_file)
  }
  uafr_require_installed(c("ChemmineR", "fmcsR", "jsonlite", "webchem"))

  uafr_step(4, 5, paste("Installing uafR from the local bundle archive"),
            log_file = log_file)
  if (identical(action, "update") && requireNamespace("uafR", quietly = TRUE)) {
    uafr_note("Existing uafR version: ",
              as.character(utils::packageVersion("uafR")), log_file = log_file)
  }
  utils::install.packages(archive, repos = NULL, type = "source")

  uafr_step(5, 5, "Verifying the installation", log_file = log_file)
  uafr_run_verification(bundle_dir, log_file = log_file)
  invisible(TRUE)
}
