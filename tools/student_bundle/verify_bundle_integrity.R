# uafR offline student bundle integrity verification.
# This checks generated bundle file sizes and MD5 checksums.

script_dir <- (function() {
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
})()

helper_path <- file.path(script_dir, "bundle_helpers.R")
if (!file.exists(helper_path)) {
  helper_path <- file.path(script_dir, "..", "bundle_helpers.R")
}
source(helper_path)

bundle_dir <- uafr_bundle_root(script_dir)
log_file <- uafr_log_file(bundle_dir, "uafR_integrity_log")
checksum_file <- file.path(bundle_dir, "CHECKSUMS.csv")

uafr_header("uafR Student Bundle Integrity Check", log_file = log_file)

if (!file.exists(checksum_file)) {
  stop("CHECKSUMS.csv was not found in the bundle folder: ", bundle_dir,
       call. = FALSE)
}

expected <- utils::read.csv(checksum_file, stringsAsFactors = FALSE)
required_columns <- c("Path", "Bytes", "MD5")
missing_columns <- setdiff(required_columns, names(expected))
if (length(missing_columns) > 0) {
  stop("CHECKSUMS.csv is missing required column(s): ",
       paste(missing_columns, collapse = ", "), call. = FALSE)
}

check_one <- function(path, expected_bytes, expected_md5) {
  full_path <- file.path(bundle_dir, path)
  if (!file.exists(full_path)) {
    return(data.frame(
      Path = path,
      Status = "FAIL",
      Details = "File is missing.",
      stringsAsFactors = FALSE
    ))
  }

  actual_bytes <- unname(file.info(full_path)$size)
  actual_md5 <- unname(tools::md5sum(full_path))
  size_ok <- identical(as.numeric(actual_bytes), as.numeric(expected_bytes))
  md5_ok <- identical(tolower(actual_md5), tolower(expected_md5))

  status <- if (size_ok && md5_ok) "PASS" else "FAIL"
  details <- if (identical(status, "PASS")) {
    paste0(actual_bytes, " bytes; MD5 matched")
  } else {
    paste0(
      "Expected ", expected_bytes, " bytes / ", expected_md5,
      "; found ", actual_bytes, " bytes / ", actual_md5
    )
  }

  data.frame(Path = path, Status = status, Details = details,
             stringsAsFactors = FALSE)
}

results <- do.call(
  rbind,
  Map(check_one, expected$Path, expected$Bytes, expected$MD5)
)

summary_table <- data.frame(
  Status = c("PASS", "FAIL"),
  Files = c(sum(results$Status == "PASS"), sum(results$Status == "FAIL")),
  stringsAsFactors = FALSE
)

uafr_note("Checksum file: ", checksum_file, log_file = log_file)
uafr_note("Files checked: ", nrow(results), log_file = log_file)
uafr_note("", log_file = log_file)
print(summary_table, row.names = FALSE)
capture.output(print(summary_table, row.names = FALSE),
               file = log_file, append = TRUE)

failures <- results[results$Status == "FAIL", , drop = FALSE]
if (nrow(failures) > 0) {
  uafr_note("", log_file = log_file)
  uafr_note("Failed file check(s):", log_file = log_file)
  print(failures, row.names = FALSE)
  capture.output(print(failures, row.names = FALSE),
                 file = log_file, append = TRUE)
  uafr_note("", log_file = log_file)
  uafr_note("Bundle integrity result: FAIL.", log_file = log_file)
  stop("Bundle integrity check failed. See ", log_file, call. = FALSE)
}

uafr_note("", log_file = log_file)
uafr_note("Bundle integrity result: PASS.", log_file = log_file)
uafr_note("Log file: ", log_file, log_file = log_file)
