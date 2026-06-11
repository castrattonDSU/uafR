#!/usr/bin/env Rscript

file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path = if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[[1]])
} else {
  file.path(getwd(), "tools", "check_package_release.R")
}
repo_root = normalizePath(file.path(dirname(script_path), ".."),
                          winslash = "/", mustWork = FALSE)

parse_args = function(args) {
  out = list()
  i = 1L
  while (i <= length(args)) {
    arg = args[[i]]
    if (!grepl("^--", arg)) {
      stop("Unexpected argument: ", arg, call. = FALSE)
    }
    key = sub("^--", "", arg)
    value = TRUE
    if (grepl("=", key, fixed = TRUE)) {
      parts = strsplit(key, "=", fixed = TRUE)[[1]]
      key = parts[[1]]
      value = paste(parts[-1], collapse = "=")
    } else if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
      value = args[[i + 1L]]
      i = i + 1L
    }
    out[[gsub("-", "_", key)]] = value
    i = i + 1L
  }
  out
}

usage = function() {
  cat(
    "Usage:\n",
    "Rscript tools/check_package_release.R [options]\n\n",
    "Options:\n",
    "  --out-dir DIR        Directory for build/check artifacts.\n",
    "  --as-cran            Add --as-cran to R CMD check.\n",
    "  --run-live-tests     Set UAFR_RUN_LIVE_TESTS=true for integration tests.\n",
    "  --help, -h           Show this message.\n\n",
    "This script builds a source tarball outside the repository and checks the\n",
    "built artifact. That is the production package-check path; checking the\n",
    "live source directory can report local hidden files and generated check\n",
    "directories that are not included in source builds.\n",
    sep = ""
  )
}

flag_arg = function(args, name, default = FALSE) {
  value = args[[name]]
  if (is.null(value)) return(default)
  if (identical(value, TRUE)) return(TRUE)
  tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

optional_arg = function(args, name, default = NULL) {
  value = args[[name]]
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    return(default)
  }
  value
}

run_cmd = function(args, workdir, log_file) {
  old = setwd(workdir)
  on.exit(setwd(old), add = TRUE)
  status = system2(file.path(R.home("bin"), "R"),
                   args = args,
                   stdout = log_file,
                   stderr = log_file)
  if (!identical(status, 0L)) {
    cat(readLines(log_file, warn = FALSE), sep = "\n")
    stop("Command failed: R ", paste(args, collapse = " "),
         "\nLog file: ", log_file, call. = FALSE)
  }
  invisible(TRUE)
}

main = function() {
  args = commandArgs(trailingOnly = TRUE)
  if (any(args %in% c("--help", "-h"))) {
    usage()
    quit(status = 0L)
  }
  parsed = parse_args(args)
  if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
    stop("Could not find DESCRIPTION at repository root: ", repo_root,
         call. = FALSE)
  }

  desc = read.dcf(file.path(repo_root, "DESCRIPTION"))
  package = desc[1, "Package"]
  version = desc[1, "Version"]
  stamp = format(Sys.time(), "%Y%m%d_%H%M%S")
  out_dir = optional_arg(
    parsed,
    "out_dir",
    file.path(tempdir(), paste0(package, "_release_check_", stamp))
  )
  out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  build_dir = file.path(out_dir, "build")
  check_dir = file.path(out_dir, "check")
  dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

  Sys.setenv(UAFR_RUN_LIVE_TESTS = ifelse(flag_arg(parsed, "run_live_tests"),
                                          "true", "false"))

  build_log = file.path(out_dir, "R-CMD-build.log")
  check_log = file.path(out_dir, "R-CMD-check.log")
  tarball = file.path(build_dir, paste0(package, "_", version, ".tar.gz"))

  cat("Repository: ", repo_root, "\n", sep = "")
  cat("Output:     ", out_dir, "\n", sep = "")
  cat("Live tests: ", Sys.getenv("UAFR_RUN_LIVE_TESTS"), "\n", sep = "")
  cat("Building source tarball...\n")
  run_cmd(c("CMD", "build", repo_root), build_dir, build_log)
  if (!file.exists(tarball)) {
    built = list.files(build_dir, pattern = "[.]tar[.]gz$", full.names = TRUE)
    tarball = if (length(built) > 0) built[[1]] else tarball
  }
  if (!file.exists(tarball)) {
    stop("Build completed, but source tarball was not found in ", build_dir,
         call. = FALSE)
  }

  check_args = c("CMD", "check", "--no-manual")
  if (flag_arg(parsed, "as_cran")) check_args = c(check_args, "--as-cran")
  check_args = c(check_args, tarball)
  cat("Checking source tarball...\n")
  run_cmd(check_args, check_dir, check_log)

  cat("Release check complete.\n")
  cat("Tarball:   ", tarball, "\n", sep = "")
  cat("Build log: ", build_log, "\n", sep = "")
  cat("Check log: ", check_log, "\n", sep = "")
}

main()
