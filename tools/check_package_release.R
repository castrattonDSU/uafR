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
    "  --skip-install-check Skip clean-library installation of the built tarball.\n",
    "  --student-bundle DIR Run the student bundle acceptance test in DIR.\n",
    "  --release-manifest FILE\n",
    "                       Write a checked JSON/CSV release manifest. This\n",
    "                       requires a clean Git worktree and runs devtools tests.\n",
    "  --help, -h           Show this message.\n\n",
    "This script builds a source tarball outside the repository and checks the\n",
    "built artifact. That is the production package-check path; checking the\n",
    "live source directory can report local hidden files and generated check\n",
    "directories that are not included in source builds.\n",
    sep = ""
  )
}

git_value = function(args) {
  git = Sys.which("git")
  if (!nzchar(git)) stop("Git is required for a release manifest.",
                         call. = FALSE)
  output = system2(git, c("-C", shQuote(repo_root), args),
                   stdout = TRUE, stderr = TRUE)
  status = attr(output, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop("Git command failed: git ", paste(args, collapse = " "),
         "\n", paste(output, collapse = "\n"), call. = FALSE)
  }
  output
}

require_clean_git = function() {
  status = git_value(c("status", "--porcelain=v1", "--untracked-files=all"))
  status = status[nzchar(status)]
  if (length(status) > 0L) {
    preview = paste(utils::head(status, 25L), collapse = "\n")
    stop(
      "Refusing to generate a release manifest from a dirty worktree.\n",
      preview,
      if (length(status) > 25L) "\n... additional paths omitted" else "",
      call. = FALSE
    )
  }
  trimws(git_value(c("rev-parse", "HEAD"))[[1L]])
}

sha256_file = function(path) {
  sha256sum = Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    output = system2(sha256sum, shQuote(path), stdout = TRUE, stderr = TRUE)
    return(strsplit(trimws(output[[1L]]), "[[:space:]]+")[[1L]][[1L]])
  }
  shasum = Sys.which("shasum")
  if (nzchar(shasum)) {
    output = system2(shasum, c("-a", "256", shQuote(path)),
                     stdout = TRUE, stderr = TRUE)
    return(strsplit(trimws(output[[1L]]), "[[:space:]]+")[[1L]][[1L]])
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    con = file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    return(gsub(":", "", tolower(as.character(openssl::sha256(con))),
                fixed = TRUE))
  }
  stop("No SHA-256 implementation is available. Install `openssl` or provide",
       " sha256sum/shasum before generating a release manifest.",
       call. = FALSE)
}

atomic_release_manifest = function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp = tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temp, force = TRUE), add = TRUE)
  if (grepl("[.]json$", path, ignore.case = TRUE)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("jsonlite is required for a JSON release manifest.", call. = FALSE)
    }
    jsonlite::write_json(x, temp, pretty = TRUE, auto_unbox = TRUE,
                         na = "null")
  } else {
    utils::write.csv(as.data.frame(x, stringsAsFactors = FALSE), temp,
                     row.names = FALSE, na = "")
  }
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(temp, path)) {
    stop("Could not atomically write release manifest: ", path,
         call. = FALSE)
  }
  invisible(path)
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

run_r_expr = function(expr, workdir, log_file, env = character()) {
  old = setwd(workdir)
  on.exit(setwd(old), add = TRUE)
  status = system2(file.path(R.home("bin"), "R"),
                   args = c("--vanilla", "-s", "-e", shQuote(expr)),
                   stdout = log_file,
                   stderr = log_file,
                   env = env)
  if (!identical(status, 0L)) {
    cat(readLines(log_file, warn = FALSE), sep = "\n")
    stop("R expression failed. Log file: ", log_file, call. = FALSE)
  }
  invisible(TRUE)
}

run_rscript = function(script, args, workdir, log_file, env = character()) {
  old = setwd(workdir)
  on.exit(setwd(old), add = TRUE)
  status = system2(file.path(R.home("bin"), "Rscript"),
                   args = c(script, args),
                   stdout = log_file,
                   stderr = log_file,
                   env = env)
  if (!identical(status, 0L)) {
    cat(readLines(log_file, warn = FALSE), sep = "\n")
    stop("Rscript failed: ", script, "\nLog file: ", log_file,
         call. = FALSE)
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
  release_manifest = optional_arg(parsed, "release_manifest", NULL)
  git_commit = NA_character_
  if (!is.null(release_manifest)) {
    release_manifest = normalizePath(release_manifest, winslash = "/",
                                     mustWork = FALSE)
    git_commit = require_clean_git()
  }
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
  install_log = file.path(out_dir, "clean-install.log")
  student_log = file.path(out_dir, "student-bundle-acceptance.log")
  test_log = file.path(out_dir, "devtools-test.log")
  tarball = file.path(build_dir, paste0(package, "_", version, ".tar.gz"))

  if (!is.null(release_manifest)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      stop("devtools is required when --release-manifest is requested.",
           call. = FALSE)
    }
    cat("Running package tests for the release manifest...\n")
    run_r_expr(
      paste0("devtools::test(", shQuote(repo_root),
             ", reporter = 'summary', stop_on_failure = TRUE)"),
      repo_root, test_log
    )
  }

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

  if (!flag_arg(parsed, "skip_install_check")) {
    install_lib = file.path(out_dir, "clean_library")
    dir.create(install_lib, recursive = TRUE, showWarnings = FALSE)
    install_expr = paste(
      "install.packages(",
      shQuote(tarball),
      ", repos = NULL, type = 'source', lib = ",
      shQuote(install_lib),
      ");",
      ".libPaths(c(",
      shQuote(install_lib),
      ", .libPaths()));",
      "library(",
      package,
      ");",
      "stopifnot(is.data.frame(uafRWorkflowGuide()));",
      "stopifnot(is.data.frame(uafRApiStability()));",
      "stopifnot(is.data.frame(uafRClaimGuidance()));",
      sep = ""
    )
    cat("Installing built tarball into a clean temporary library...\n")
    run_r_expr(install_expr, out_dir, install_log,
               env = paste0("R_LIBS_USER=", install_lib))
  }

  student_bundle = optional_arg(parsed, "student_bundle", NULL)
  if (!is.null(student_bundle)) {
    student_bundle = normalizePath(student_bundle, winslash = "/",
                                   mustWork = FALSE)
    acceptance = file.path(student_bundle, "run_student_acceptance_test.R")
    if (!file.exists(acceptance)) {
      stop("Student bundle acceptance script was not found: ", acceptance,
           call. = FALSE)
    }
    cat("Running student bundle offline acceptance test...\n")
    run_rscript("run_student_acceptance_test.R", character(), student_bundle,
                student_log,
                env = c("UAFR_STUDENT_LIVE=0",
                        "UAFR_RUN_LIVE_TESTS=false"))
  }

  cat("Release check complete.\n")
  cat("Tarball:   ", tarball, "\n", sep = "")
  cat("Build log: ", build_log, "\n", sep = "")
  cat("Check log: ", check_log, "\n", sep = "")
  if (!flag_arg(parsed, "skip_install_check")) {
    cat("Install log:", install_log, "\n")
  }
  if (!is.null(optional_arg(parsed, "student_bundle", NULL))) {
    cat("Student log:", student_log, "\n")
  }
  if (!is.null(release_manifest)) {
    release = list(
      workflow = "uafR_private_release_candidate",
      package = package,
      package_version = version,
      git_commit = git_commit,
      dirty_state_check = "clean",
      built_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      r_version = R.version.string,
      platform = R.version$platform,
      source_tarball = basename(tarball),
      source_tarball_path = basename(tarball),
      source_tarball_bytes = unname(file.info(tarball)$size),
      source_tarball_md5 = unname(tools::md5sum(tarball)[[1L]]),
      source_tarball_sha256 = sha256_file(tarball),
      devtools_test_result = "pass",
      r_cmd_check_result = "pass",
      clean_library_install_result = if (
        flag_arg(parsed, "skip_install_check")
      ) "skipped" else "pass",
      live_tests_enabled = Sys.getenv("UAFR_RUN_LIVE_TESTS")
    )
    atomic_release_manifest(release, release_manifest)
    cat("Release manifest:", release_manifest, "\n")
  }
}

main()
