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

read_text <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"),
        collapse = "\n")
}

png_dimensions <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  header <- readBin(con, what = "raw", n = 24L)
  if (length(header) != 24L ||
      !identical(as.integer(header[1:8]),
                 c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))) {
    stop("Not a valid PNG header: ", path, call. = FALSE)
  }
  decode <- function(bytes) {
    sum(as.integer(bytes) * c(256^3, 256^2, 256, 1))
  }
  c(width = decode(header[17:20]), height = decode(header[21:24]))
}

pdf_page_count <- function(path) {
  pdfinfo <- Sys.which("pdfinfo")
  if (!nzchar(pdfinfo)) return(NA_integer_)
  output <- system2(pdfinfo, path, stdout = TRUE, stderr = TRUE)
  line <- grep("^Pages:", output, value = TRUE)
  if (length(line) == 0L) return(NA_integer_)
  suppressWarnings(as.integer(sub("^Pages:[[:space:]]*", "", line[[1]])))
}

pdf_text <- function(path) {
  pdftotext <- Sys.which("pdftotext")
  if (!nzchar(pdftotext)) return(NA_character_)
  paste(
    system2(pdftotext, c(path, "-"), stdout = TRUE, stderr = TRUE),
    collapse = "\n"
  )
}

repo_root <- find_repo_root()
training_dir <- file.path(repo_root, "training")
issues <- character()
pass <- function(message) cat("PASS:", message, "\n")
fail <- function(message) issues <<- c(issues, message)

required_files <- c(
  "training/main.tex",
  "training/instructor_guide.tex",
  "training/preamble.tex",
  "training/PUBLICATION_NOTICE.md",
  "training/assets/dsdna_core_logo.png",
  "training/uafR_training_manual.pdf",
  "training/uafR_instructor_guide.pdf",
  "training/program_notes/instructor_preface.tex",
  "training/program_notes/program_implementation_notes.tex",
  "training/program_notes/dsdna_core_review_notes.tex",
  "training/program_notes/bundle_maintenance.tex"
)
missing <- required_files[!file.exists(file.path(repo_root, required_files))]
if (length(missing) > 0L) {
  fail(paste("Missing required publication files:",
             paste(missing, collapse = ", ")))
} else {
  pass("required manual sources, PDFs, logo, and publication notice exist")
}

main_text <- read_text(file.path(training_dir, "main.tex"))
instructor_text <- read_text(file.path(training_dir, "instructor_guide.tex"))
if (grepl("program_notes/", main_text, fixed = TRUE)) {
  fail("Student main.tex includes instructor-only program_notes content.")
} else {
  pass("student edition excludes instructor-only source files")
}
required_instructor_includes <- c(
  "program_notes/program_implementation_notes",
  "program_notes/dsdna_core_review_notes",
  "program_notes/bundle_maintenance"
)
missing_instructor <- required_instructor_includes[
  !vapply(
    required_instructor_includes,
    grepl,
    logical(1),
    x = instructor_text,
    fixed = TRUE
  )
]
if (length(missing_instructor) > 0L) {
  fail(paste("Instructor guide is missing:",
             paste(missing_instructor, collapse = ", ")))
} else {
  pass("instructor edition includes implementation, review, and release guidance")
}

if (grepl("\\\\today", main_text) ||
    grepl("\\\\today", instructor_text)) {
  fail("Public manuals use nondeterministic \\today metadata.")
} else {
  pass("manuals use fixed version and revision metadata")
}

daily_text <- read_text(file.path(training_dir, "chapters", "daily_plans.tex"))
day_matches <- gregexpr(
  "\\\\begin\\{dailybox\\}\\{Day ([0-9]+):",
  daily_text,
  perl = TRUE
)
day_strings <- regmatches(daily_text, day_matches)[[1]]
days <- as.integer(sub(
  ".*Day ([0-9]+):.*",
  "\\1",
  day_strings,
  perl = TRUE
))
if (!identical(days, 1:70)) {
  fail(
    paste0(
      "Daily plan must contain Day 1 through Day 70 exactly once; found: ",
      paste(days, collapse = ", ")
    )
  )
} else {
  pass("student manual contains the complete ordered 70-day sequence")
}

required_scripts <- c(
  "00_install_check.R",
  "01_project_setup.R",
  "04_core_workflow.R",
  "05_categorate_research.R",
  "06_trait_matrix_analysis.R",
  "07_visualization.R",
  "09_reproducibility_check.R",
  "10_species_phytochemistry.R",
  "run_training_smoke_tests.R"
)
missing_scripts <- required_scripts[
  !file.exists(file.path(training_dir, "scripts", required_scripts))
]
if (length(missing_scripts) > 0L) {
  fail(paste("Missing training scripts:",
             paste(missing_scripts, collapse = ", ")))
} else {
  pass("all documented training scripts exist")
}

namespace <- readLines(file.path(repo_root, "NAMESPACE"), warn = FALSE)
export_lines <- grep("^export\\(", namespace, value = TRUE)
exports <- sub("^export\\((.*)\\)$", "\\1", export_lines)
exports <- gsub('^"|"$', "", exports)
required_functions <- c(
  "spreadOut",
  "mzExacto",
  "exactoThese",
  "categorate",
  "validateCategorateResult",
  "chemicalTraitMatrix",
  "chemicalTraitOntologyMatrix",
  "chemicalTraitEvidence",
  "resolvePlantPhytochemistry",
  "filterPlantPhytochemistryEvidence",
  "plantPhytochemistryMatrix",
  "scorePlantChemistryCandidates",
  "exportPlantPhytochemistryWorkbook"
)
missing_functions <- setdiff(required_functions, exports)
if (length(missing_functions) > 0L) {
  fail(paste("Documented uafR functions are not exported:",
             paste(missing_functions, collapse = ", ")))
} else {
  pass("core and species-first functions used by the curriculum are exported")
}

logo_path <- file.path(training_dir, "assets", "dsdna_core_logo.png")
if (file.exists(logo_path)) {
  dimensions <- png_dimensions(logo_path)
  aspect <- unname(dimensions[["width"]] / dimensions[["height"]])
  if (dimensions[["width"]] < 1200L ||
      aspect < 2.0 || aspect > 2.5) {
    fail(
      paste0(
        "Logo dimensions/aspect are unsuitable: ",
        dimensions[["width"]],
        "x",
        dimensions[["height"]],
        ", aspect ",
        round(aspect, 3)
      )
    )
  } else {
    pass(
      paste0(
        "official logo is print-sized and aspect-preserving (",
        dimensions[["width"]],
        "x",
        dimensions[["height"]],
        ")"
      )
    )
  }
}

scan_roots <- c(
  file.path(training_dir),
  file.path(repo_root, "tools", "student_bundle")
)
scan_files <- unlist(lapply(
  scan_roots,
  list.files,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "[.](tex|md|R|csv|txt)$"
))
forbidden <- c(
  "current private uafR",
  "package remains private",
  "keeps the package repository private",
  "private student bundle",
  "/Users/[A-Za-z0-9._-]+/",
  "[A-Za-z]:\\\\\\\\Users\\\\\\\\",
  "\\bTODO\\b",
  "\\bFIXME\\b",
  "\\bplaceholder\\b",
  "\\bCodex\\b"
)
for (path in scan_files) {
  text <- read_text(path)
  hits <- forbidden[vapply(
    forbidden,
    grepl,
    logical(1),
    x = text,
    ignore.case = TRUE,
    perl = TRUE
  )]
  if (length(hits) > 0L) {
    fail(
      paste0(
        "Public source scan failed for ",
        substring(path, nchar(repo_root) + 2L),
        ": ",
        paste(hits, collapse = ", ")
      )
    )
  }
}
if (!any(grepl("^Public source scan failed", issues))) {
  pass("public materials contain no stale private-repository claims, local paths, placeholders, or development markers")
}

plant_script <- read_text(file.path(
  training_dir,
  "scripts",
  "10_species_phytochemistry.R"
))
if (!grepl("Simulated training row", plant_script, fixed = TRUE)) {
  fail("Species-first teaching rows are not explicitly labeled simulated.")
} else {
  pass("species-first teaching rows are explicitly labeled simulated")
}

training_data_files <- list.files(
  file.path(training_dir, "data"),
  full.names = TRUE,
  pattern = "[.](md|csv)$",
  ignore.case = TRUE
)
training_data_text <- paste(
  vapply(training_data_files, read_text, character(1)),
  collapse = "\n"
)
if (!grepl("simulated", training_data_text, ignore.case = TRUE)) {
  fail("Training data/templates are not explicitly labeled simulated.")
} else {
  pass("training data templates explicitly distinguish simulated rows from evidence")
}

student_pdf <- file.path(training_dir, "uafR_training_manual.pdf")
instructor_pdf <- file.path(training_dir, "uafR_instructor_guide.pdf")
if (file.exists(student_pdf) && file.exists(instructor_pdf)) {
  student_pages <- pdf_page_count(student_pdf)
  instructor_pages <- pdf_page_count(instructor_pdf)
  if (!is.na(student_pages) && student_pages < 100L) {
    fail(paste("Student manual is unexpectedly short:", student_pages,
               "pages."))
  }
  if (!is.na(instructor_pages) && instructor_pages < 10L) {
    fail(paste("Instructor guide is unexpectedly short:", instructor_pages,
               "pages."))
  }
  if (is.na(student_pages) || is.na(instructor_pages)) {
    pass("manual PDFs exist; pdfinfo is unavailable for page-count checks")
  } else {
    pass(
      paste0(
        "manual page counts are plausible (student ",
        student_pages,
        "; instructor ",
        instructor_pages,
        ")"
      )
    )
  }

  student_pdf_text <- pdf_text(student_pdf)
  instructor_pdf_text <- pdf_text(instructor_pdf)
  instructor_only_titles <- c(
    "Program Implementation Guide",
    "Instructor Review Notes and Answer Guidance",
    "Bundle Maintenance and Curriculum Release"
  )
  if (!is.na(student_pdf_text) &&
      any(vapply(
        instructor_only_titles,
        grepl,
        logical(1),
        x = student_pdf_text,
        fixed = TRUE
      ))) {
    fail("Student PDF contains instructor-only chapters.")
  } else if (!is.na(student_pdf_text)) {
    pass("student PDF text excludes instructor-only chapters")
  }
  if (!is.na(instructor_pdf_text) &&
      !all(vapply(
        instructor_only_titles,
        grepl,
        logical(1),
        x = instructor_pdf_text,
        fixed = TRUE
      ))) {
    fail("Instructor PDF is missing one or more instructor-only chapters.")
  } else if (!is.na(instructor_pdf_text)) {
    pass("instructor PDF contains all intended instructor-only chapters")
  }
}

if (length(issues) > 0L) {
  cat("\nTraining publication audit failed:\n")
  cat(paste0("- ", unique(issues), "\n"))
  stop(length(unique(issues)), " publication issue(s) require attention.",
       call. = FALSE)
}

cat("\nTraining publication audit: PASS\n")
