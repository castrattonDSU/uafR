#!/usr/bin/env Rscript

find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    description <- file.path(current, "DESCRIPTION")
    if (file.exists(description)) {
      dcf <- tryCatch(read.dcf(description), error = function(e) NULL)
      if (!is.null(dcf) &&
          identical(unname(dcf[1, "Package"]), "uafR")) {
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

run_command <- function(command, args, wd) {
  old <- setwd(wd)
  on.exit(setwd(old), add = TRUE)
  output <- system2(command, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(status = status, output = output)
}

pdf_page_count <- function(path) {
  pdfinfo <- Sys.which("pdfinfo")
  if (!nzchar(pdfinfo)) return(NA_integer_)
  output <- system2(pdfinfo, path, stdout = TRUE, stderr = TRUE)
  line <- grep("^Pages:", output, value = TRUE)
  if (length(line) == 0L) return(NA_integer_)
  suppressWarnings(as.integer(sub("^Pages:[[:space:]]*", "", line[[1]])))
}

build_manual <- function(latexmk, training_dir, source, output) {
  source_pdf <- file.path(
    training_dir,
    paste0(tools::file_path_sans_ext(source), ".pdf")
  )
  stable_pdf <- file.path(training_dir, output)
  log_path <- file.path(
    training_dir,
    paste0(tools::file_path_sans_ext(source), ".log")
  )

  message("Building ", source)
  build <- run_command(
    latexmk,
    c(
      "-silent",
      "-pdf",
      "-interaction=nonstopmode",
      "-halt-on-error",
      source
    ),
    wd = training_dir
  )
  verbose <- identical(
    tolower(Sys.getenv("UAFR_LATEX_VERBOSE", "false")),
    "true"
  )
  if (!identical(build$status, 0L) || verbose) {
    cat(paste(build$output, collapse = "\n"), "\n")
  }
  if (!identical(build$status, 0L)) {
    stop("LaTeX build failed for ", source, ".", call. = FALSE)
  }
  if (!file.exists(source_pdf) || file.info(source_pdf)$size < 50000) {
    stop("LaTeX did not create a usable PDF for ", source, ".", call. = FALSE)
  }

  log_text <- if (file.exists(log_path)) {
    paste(readLines(log_path, warn = FALSE), collapse = "\n")
  } else {
    ""
  }
  unresolved <- c(
    "There were undefined references",
    "Citation [`'][^\\n]+[`'] on page [^\\n]+ undefined",
    "Reference [`'][^\\n]+[`'] on page [^\\n]+ undefined"
  )
  unresolved_found <- vapply(
    unresolved,
    grepl,
    logical(1),
    x = log_text,
    perl = TRUE
  )
  if (any(unresolved_found)) {
    stop(
      "Unresolved references or citations remain in ",
      basename(log_path),
      ".",
      call. = FALSE
    )
  }

  copied <- file.copy(source_pdf, stable_pdf, overwrite = TRUE)
  if (!copied) {
    stop("Could not write stable PDF: ", stable_pdf, call. = FALSE)
  }

  pages <- pdf_page_count(stable_pdf)
  message(
    "Created ",
    stable_pdf,
    if (is.na(pages)) "" else paste0(" (", pages, " pages)")
  )

  clean <- run_command(latexmk, c("-c", source), wd = training_dir)
  if (!identical(clean$status, 0L)) {
    warning("latexmk cleanup returned a nonzero status for ", source, ".")
  }
  if (!identical(normalizePath(source_pdf, mustWork = FALSE),
                 normalizePath(stable_pdf, mustWork = FALSE))) {
    unlink(source_pdf, force = TRUE)
  }
  invisible(stable_pdf)
}

args <- commandArgs(trailingOnly = TRUE)
allowed <- c("--student-only", "--instructor-only")
unknown <- setdiff(args, allowed)
if (length(unknown) > 0L) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       call. = FALSE)
}
if (all(allowed %in% args)) {
  stop("Choose only one of --student-only or --instructor-only.",
       call. = FALSE)
}

repo_root <- find_repo_root()
training_dir <- file.path(repo_root, "training")
latexmk <- Sys.which("latexmk")
if (!nzchar(latexmk)) {
  stop(
    "latexmk is required to build the public manuals. Install a LaTeX ",
    "distribution with latexmk and rerun this command.",
    call. = FALSE
  )
}

targets <- data.frame(
  source = c("main.tex", "instructor_guide.tex"),
  output = c("uafR_training_manual.pdf", "uafR_instructor_guide.pdf"),
  stringsAsFactors = FALSE
)
if ("--student-only" %in% args) targets <- targets[1, , drop = FALSE]
if ("--instructor-only" %in% args) targets <- targets[2, , drop = FALSE]

for (index in seq_len(nrow(targets))) {
  build_manual(
    latexmk = latexmk,
    training_dir = training_dir,
    source = targets$source[[index]],
    output = targets$output[[index]]
  )
}

message("Training manual build completed.")
