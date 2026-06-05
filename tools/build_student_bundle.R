#!/usr/bin/env Rscript

find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    desc <- file.path(current, "DESCRIPTION")
    if (file.exists(desc)) {
      dcf <- tryCatch(read.dcf(desc), error = function(e) NULL)
      if (!is.null(dcf) && identical(unname(dcf[1, "Package"]), "uafR")) {
        return(current)
      }
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find the uafR repository root.", call. = FALSE)
    }
    current <- parent
  }
}

copy_required <- function(from, to) {
  if (!file.exists(from)) {
    stop("Required file is missing: ", from, call. = FALSE)
  }
  ok <- file.copy(from, to, overwrite = TRUE, copy.date = TRUE)
  if (!ok) {
    stop("Could not copy ", from, " to ", to, call. = FALSE)
  }
}

copy_dir_files <- function(from, to, pattern = NULL) {
  if (!dir.exists(from)) return(invisible(FALSE))
  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(from, pattern = pattern, full.names = TRUE)
  for (file in files) {
    copy_required(file, file.path(to, basename(file)))
  }
  invisible(TRUE)
}

run_command <- function(command, args, wd = getwd()) {
  old <- setwd(wd)
  on.exit(setwd(old), add = TRUE)
  output <- system2(command, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(status = status, output = output)
}

pdf_page_count <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  pdfinfo <- Sys.which("pdfinfo")
  if (!nzchar(pdfinfo)) return(NA_integer_)
  info <- tryCatch(
    system2(pdfinfo, path, stdout = TRUE, stderr = TRUE),
    error = function(e) character()
  )
  page_line <- grep("^Pages:", info, value = TRUE)
  if (length(page_line) == 0) return(NA_integer_)
  suppressWarnings(as.integer(sub("^Pages:[[:space:]]*", "", page_line[[1]])))
}

format_md5 <- function(path) {
  if (!file.exists(path)) return("missing")
  unname(tools::md5sum(path))
}

bundle_file_table <- function(bundle_dir, exclude = "CHECKSUMS.csv") {
  files <- list.files(bundle_dir, recursive = TRUE, full.names = TRUE,
                      all.files = FALSE, no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  relative <- substring(files, nchar(bundle_dir) + 2L)
  keep <- !(relative %in% exclude)
  files <- files[keep]
  relative <- relative[keep]
  data.frame(
    Path = relative,
    Bytes = unname(file.info(files)$size),
    MD5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
}

repo_root <- find_repo_root()
args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else "student_bundle"
output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)

desc <- read.dcf(file.path(repo_root, "DESCRIPTION"))
package <- unname(desc[1, "Package"])
version <- unname(desc[1, "Version"])
bundle_name <- paste0(package, "_student_bundle_", version)
bundle_dir <- file.path(output_root, bundle_name)
zip_path <- file.path(output_root, paste0(bundle_name, ".zip"))

message("Building ", package, " ", version, " student bundle")
message("Repository: ", repo_root)

build_dir <- tempfile("uafr-package-build-")
dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

r_binary <- file.path(R.home("bin"), "R")
build <- run_command(
  r_binary,
  c("CMD", "build", "--no-build-vignettes", "--no-manual", repo_root),
  wd = build_dir
)
cat(paste(build$output, collapse = "\n"), "\n")
if (!identical(build$status, 0L)) {
  stop("R CMD build failed.", call. = FALSE)
}

tarball <- list.files(
  build_dir,
  pattern = paste0("^", package, "_", gsub(".", "[.]", version, fixed = TRUE), "[.]tar[.]gz$"),
  full.names = TRUE
)
if (length(tarball) != 1) {
  stop("Could not find the built package archive in ", build_dir, call. = FALSE)
}

if (dir.exists(bundle_dir)) {
  unlink(bundle_dir, recursive = TRUE, force = TRUE)
}
if (file.exists(zip_path)) {
  unlink(zip_path, force = TRUE)
}

dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(bundle_dir, "packages"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(bundle_dir, "training"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(bundle_dir, "examples"), recursive = TRUE, showWarnings = FALSE)

template_dir <- file.path(repo_root, "tools", "student_bundle")
copy_required(file.path(template_dir, "START_HERE.md"),
              file.path(bundle_dir, "START_HERE.md"))
copy_required(file.path(template_dir, "README_STUDENT_INSTALL.md"),
              file.path(bundle_dir, "README_STUDENT_INSTALL.md"))
copy_required(file.path(template_dir, "uafR_QUICK_REFERENCE.md"),
              file.path(bundle_dir, "uafR_QUICK_REFERENCE.md"))
copy_required(file.path(template_dir, "bundle_helpers.R"),
              file.path(bundle_dir, "bundle_helpers.R"))
copy_required(file.path(template_dir, "preflight_check.R"),
              file.path(bundle_dir, "preflight_check.R"))
copy_required(file.path(template_dir, "verify_bundle_integrity.R"),
              file.path(bundle_dir, "verify_bundle_integrity.R"))
copy_required(file.path(template_dir, "install_uafR_from_bundle.R"),
              file.path(bundle_dir, "install_uafR_from_bundle.R"))
copy_required(file.path(template_dir, "update_uafR_from_bundle.R"),
              file.path(bundle_dir, "update_uafR_from_bundle.R"))
copy_required(file.path(template_dir, "verify_uafR_install.R"),
              file.path(bundle_dir, "verify_uafR_install.R"))
copy_required(file.path(template_dir, "run_student_acceptance_test.R"),
              file.path(bundle_dir, "run_student_acceptance_test.R"))
copy_required(file.path(template_dir, "examples", "test_install.R"),
              file.path(bundle_dir, "examples", "test_install.R"))
copy_required(tarball, file.path(bundle_dir, "packages", basename(tarball)))

training_pdf <- file.path(repo_root, "training", "main.pdf")
if (file.exists(training_pdf)) {
  copy_required(training_pdf, file.path(bundle_dir, "training", "uafR_training_manual.pdf"))
} else {
  warning("training/main.pdf was not found. Render the training manual before building a student bundle.")
}

copy_dir_files(file.path(repo_root, "training", "scripts"),
               file.path(bundle_dir, "training", "scripts"),
               pattern = "[.]R$")
copy_dir_files(file.path(repo_root, "training", "data"),
               file.path(bundle_dir, "training", "data"),
               pattern = "[.](csv|md|txt)$")

git_commit <- tryCatch(
  run_command("git", c("rev-parse", "--short", "HEAD"), wd = repo_root)$output[[1]],
  error = function(e) NA_character_
)
git_status <- tryCatch(
  run_command("git", c("status", "--short"), wd = repo_root)$output,
  error = function(e) character()
)
manual_bundle_pdf <- file.path(bundle_dir, "training", "uafR_training_manual.pdf")
manual_pages <- pdf_page_count(manual_bundle_pdf)
manual_pages_label <- ifelse(is.na(manual_pages), "unknown", as.character(manual_pages))
archive_in_bundle <- file.path(bundle_dir, "packages", basename(tarball))
manifest <- c(
  paste0("Package: ", package),
  paste0("Version: ", version),
  paste0("Built: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R: ", R.version.string),
  paste0("Source commit: ", ifelse(is.na(git_commit), "unknown", git_commit)),
  paste0("Source tree clean: ", ifelse(length(git_status) == 0, "yes", "no")),
  paste0("Manual pages: ", manual_pages_label),
  paste0("Package archive MD5: ", format_md5(archive_in_bundle)),
  paste0("Training manual MD5: ", format_md5(manual_bundle_pdf)),
  paste0("Full checksum table: CHECKSUMS.csv"),
  "",
  "Bundle contents:",
  "- START_HERE.md: first student checklist",
  "- README_STUDENT_INSTALL.md: detailed student installation directions",
  "- uafR_QUICK_REFERENCE.md: one-page student command reference",
  "- CHECKSUMS.csv: file sizes and MD5 checksums for this generated bundle",
  "- packages/: local uafR source package archive",
  "- preflight_check.R: machine readiness check",
  "- verify_bundle_integrity.R: confirms copied bundle files match CHECKSUMS.csv",
  "- install_uafR_from_bundle.R: student installer",
  "- update_uafR_from_bundle.R: reinstall/update script for newer bundles",
  "- verify_uafR_install.R: post-install verification",
  "- run_student_acceptance_test.R: offline student readiness test",
  "- training/uafR_training_manual.pdf: training manual",
  "- training/scripts/: classroom support scripts",
  "- training/data/: classroom data templates",
  "- examples/test_install.R: rerun the acceptance test from the examples folder"
)
writeLines(manifest, file.path(bundle_dir, "MANIFEST.txt"))

checksum_table <- bundle_file_table(bundle_dir)
checksum_table <- checksum_table[order(checksum_table$Path), , drop = FALSE]
utils::write.csv(checksum_table, file.path(bundle_dir, "CHECKSUMS.csv"),
                 row.names = FALSE)

old <- setwd(output_root)
on.exit(setwd(old), add = TRUE)
zip_result <- tryCatch(
  utils::zip(zipfile = basename(zip_path), files = bundle_name, flags = "-r9X"),
  error = function(e) e
)
if (inherits(zip_result, "error")) {
  warning("Could not create zip archive: ", conditionMessage(zip_result))
} else {
  message("Created zip archive: ", zip_path)
}

message("Created bundle folder: ", bundle_dir)
message("Done.")
