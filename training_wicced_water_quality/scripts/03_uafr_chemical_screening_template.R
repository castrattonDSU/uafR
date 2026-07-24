message("WiCCED simulated GC-MS uafR-compatible pipeline")

locate_manual_dir <- function() {
  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  candidates <- character()
  if (length(cmd_file) > 0) {
    script_path <- normalizePath(sub("^--file=", "", cmd_file[[1]]), mustWork = FALSE)
    candidates <- c(candidates, dirname(dirname(script_path)), dirname(script_path))
  }
  candidates <- c(
    candidates,
    getwd(),
    dirname(getwd()),
    file.path(getwd(), "training_wicced_water_quality")
  )
  for (candidate in unique(candidates)) {
    if (dir.exists(file.path(candidate, "scripts")) &&
        dir.exists(file.path(candidate, "data")) &&
        file.exists(file.path(candidate, "README.md"))) {
      return(normalizePath(candidate, mustWork = FALSE))
    }
  }
  normalizePath(getwd(), mustWork = FALSE)
}

manual_dir <- locate_manual_dir()
data_dir <- file.path(manual_dir, "data")

args <- commandArgs(trailingOnly = TRUE)
hits_csv <- if (length(args) >= 1) args[[1]] else file.path(data_dir, "wicced_simulated_gcms_hits.csv")
queries_csv <- if (length(args) >= 2) args[[2]] else file.path(data_dir, "wicced_simulated_gcms_queries.csv")
reference_csv <- if (length(args) >= 3) args[[3]] else file.path(data_dir, "wicced_simulated_gcms_compound_reference.csv")
metadata_csv <- if (length(args) >= 4) args[[4]] else file.path(data_dir, "wicced_simulated_gcms_sample_metadata.csv")
out_dir <- if (length(args) >= 5) args[[5]] else file.path(manual_dir, "results")
mode <- if (length(args) >= 6) args[[6]] else Sys.getenv("WICCED_GCMS_MODE", "offline")
mode <- tolower(mode)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

read_required_csv <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, " file: ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

hits <- read_required_csv(hits_csv, "simulated hits")
queries <- read_required_csv(queries_csv, "query")
reference <- read_required_csv(reference_csv, "compound reference")
metadata <- read_required_csv(metadata_csv, "sample metadata")

required_cols <- c("Component.RT", "Component.Area", "Base.Peak.MZ",
                   "File.Name", "Compound.Name", "Match.Factor")
missing_cols <- setdiff(required_cols, names(hits))
if (length(missing_cols) > 0) {
  stop("Simulated GC-MS hits are missing required uafR columns: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

numeric_cols <- c("Component.RT", "Component.Area", "Base.Peak.MZ", "Match.Factor")
for (col in numeric_cols) {
  hits[[col]] <- suppressWarnings(as.numeric(hits[[col]]))
  if (any(is.na(hits[[col]]))) {
    stop("Column ", col, " contains missing or non-numeric values.", call. = FALSE)
  }
}

peak_keys <- paste(hits$Component.RT, hits$Component.Area, hits$Base.Peak.MZ, sep = " | ")
if (any(duplicated(peak_keys))) {
  stop("Simulated GC-MS hits contain duplicate uafR peak keys.", call. = FALSE)
}

if (!all(c("Compound.Name", "Exact.Mass", "Primary.MZ") %in% names(reference))) {
  stop("Compound reference must contain Compound.Name, Exact.Mass, and Primary.MZ.", call. = FALSE)
}
reference$Exact.Mass <- as.numeric(reference$Exact.Mass)
reference$Primary.MZ <- as.character(reference$Primary.MZ)

missing_ref <- setdiff(unique(hits$Compound.Name), reference$Compound.Name)
if (length(missing_ref) > 0) {
  stop("Missing simulated reference rows for: ", paste(missing_ref, collapse = ", "),
       call. = FALSE)
}

query_col <- if ("Compound.Name" %in% names(queries)) "Compound.Name" else names(queries)[[1]]
query_chemicals <- unique(queries[[query_col]])
missing_queries <- setdiff(query_chemicals, reference$Compound.Name)
if (length(missing_queries) > 0) {
  stop("Missing simulated query reference rows for: ",
       paste(missing_queries, collapse = ", "), call. = FALSE)
}

make_empty_matrix <- function(n_rows, sample_names) {
  data.frame(matrix(NA, nrow = n_rows, ncol = length(sample_names),
                    dimnames = list(NULL, sample_names)),
             check.names = FALSE)
}

build_simulated_spread <- function(hits, reference) {
  sample_names <- unique(hits$File.Name)
  spread <- list(
    Area = make_empty_matrix(nrow(hits), sample_names),
    Compounds = make_empty_matrix(nrow(hits), sample_names),
    MZ = make_empty_matrix(nrow(hits), sample_names),
    MatchFactor = make_empty_matrix(nrow(hits), sample_names),
    RT = make_empty_matrix(nrow(hits), sample_names),
    Mass = make_empty_matrix(nrow(hits), sample_names),
    rtBYmass = make_empty_matrix(nrow(hits), sample_names)
  )
  web_info <- vector("list", nrow(hits))

  for (i in seq_len(nrow(hits))) {
    sample <- hits$File.Name[[i]]
    compound <- hits$Compound.Name[[i]]
    ref_row <- reference[match(compound, reference$Compound.Name), ]
    mass <- ref_row$Exact.Mass[[1]]

    spread$Area[i, sample] <- hits$Component.Area[[i]]
    spread$Compounds[i, sample] <- compound
    spread$MZ[i, sample] <- hits$Base.Peak.MZ[[i]]
    spread$MatchFactor[i, sample] <- hits$Match.Factor[[i]]
    spread$RT[i, sample] <- hits$Component.RT[[i]]
    spread$Mass[i, sample] <- mass
    spread$rtBYmass[i, sample] <- paste0(hits$Component.RT[[i]], " | ", mass)
    web_info[[i]] <- list(
      names = compound,
      primary_mz = ref_row$Primary.MZ[[1]],
      exact_mass = mass,
      retention_time = hits$Component.RT[[i]]
    )
  }

  names(web_info) <- hits$Compound.Name
  spread$webInfo <- web_info
  spread
}

build_simulated_exact <- function(hits, reference, query_chemicals) {
  sample_names <- unique(hits$File.Name)
  rows <- lapply(query_chemicals, function(compound) {
    compound_hits <- hits[hits$Compound.Name == compound, , drop = FALSE]
    ref_row <- reference[match(compound, reference$Compound.Name), ]
    sample_areas <- vapply(sample_names, function(sample) {
      sum(compound_hits$Component.Area[compound_hits$File.Name == sample], na.rm = TRUE)
    }, numeric(1))
    detected <- compound_hits$Component.Area > 0
    weighted_rt <- if (any(detected)) {
      stats::weighted.mean(compound_hits$Component.RT[detected],
                           compound_hits$Component.Area[detected])
    } else {
      NA_real_
    }
    best_match <- if (nrow(compound_hits) > 0) max(compound_hits$Match.Factor) else NA_real_
    data.frame(
      Compound = compound,
      Mass = ref_row$Exact.Mass[[1]],
      RT = weighted_rt,
      `Best Match` = best_match,
      as.list(sample_areas),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  exact <- do.call(rbind, rows)
  rownames(exact) <- NULL
  exact
}

load_uafr_for_live_mode <- function() {
  if (requireNamespace("uafR", quietly = TRUE)) return(TRUE)
  repo_root <- dirname(dirname(manual_dir))
  if (file.exists(file.path(repo_root, "DESCRIPTION")) &&
      requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(repo_root, quiet = TRUE)
    return(exists("spreadOut", mode = "function") && exists("mzExacto", mode = "function"))
  }
  FALSE
}

if (mode == "live") {
  if (!load_uafr_for_live_mode()) {
    stop("Live mode requires uafR to be installed or loadable from the repository.",
         call. = FALSE)
  }
  message("Running live uafR spreadOut() and mzExacto(). This may query web services.")
  spread <- spreadOut(hits)
  exact <- mzExacto(spread, query_chemicals)
} else {
  message("Running offline simulated uafR-compatible spread/exact workflow.")
  spread <- build_simulated_spread(hits, reference)
  exact <- build_simulated_exact(hits, reference, query_chemicals)
}

sample_cols <- setdiff(names(exact), c("Compound", "Mass", "RT", "Best Match"))
long <- do.call(rbind, lapply(sample_cols, function(sample) {
  data.frame(
    File.Name = sample,
    Compound = exact$Compound,
    Absolute.Area = as.numeric(exact[[sample]]),
    stringsAsFactors = FALSE
  )
}))
long <- long[long$Absolute.Area > 0, , drop = FALSE]
sample_totals <- stats::aggregate(Absolute.Area ~ File.Name, long, sum)
long$Total.Sample.Area <- sample_totals$Absolute.Area[match(long$File.Name, sample_totals$File.Name)]
long$Relative.Area <- long$Absolute.Area / long$Total.Sample.Area

if ("File.Name" %in% names(metadata)) {
  long <- cbind(long, metadata[match(long$File.Name, metadata$File.Name),
                              setdiff(names(metadata), c("File.Name", "Sample.Name")),
                              drop = FALSE])
}

summary_by_context <- stats::aggregate(
  Absolute.Area ~ Site.Class + Salinity.Scenario + Water.Quality.Context + Compound,
  long,
  sum
)
summary_by_context <- summary_by_context[order(summary_by_context$Site.Class,
                                               -summary_by_context$Absolute.Area), ]

observed <- stats::aggregate(Absolute.Area ~ Compound, long, sum)
observed$Observed.Samples <- stats::aggregate(
  File.Name ~ Compound,
  long,
  function(x) length(unique(x))
)$File.Name
evidence_notes <- merge(reference, observed, by.x = "Compound.Name", by.y = "Compound",
                        all.x = TRUE)
evidence_notes$Absolute.Area[is.na(evidence_notes$Absolute.Area)] <- 0
evidence_notes$Observed.Samples[is.na(evidence_notes$Observed.Samples)] <- 0
evidence_notes$Data.Status <- "simulated_training_data_not_environmental_evidence"

relative_check <- stats::aggregate(Relative.Area ~ File.Name, long, sum)
if (max(abs(relative_check$Relative.Area - 1)) > 1e-8) {
  stop("Relative areas do not sum to 1 by sample.", call. = FALSE)
}
if (!all(query_chemicals %in% exact$Compound)) {
  stop("Not all query chemicals were represented in the exact-match output.",
       call. = FALSE)
}

saveRDS(spread, file.path(out_dir, "wicced_simulated_spread.rds"))
utils::write.csv(exact, file.path(out_dir, "wicced_simulated_exact_matches.csv"),
                 row.names = FALSE)
utils::write.csv(long, file.path(out_dir, "wicced_simulated_long_abundance.csv"),
                 row.names = FALSE)
utils::write.csv(summary_by_context,
                 file.path(out_dir, "wicced_simulated_sample_group_summary.csv"),
                 row.names = FALSE)
utils::write.csv(evidence_notes,
                 file.path(out_dir, "wicced_simulated_compound_evidence_notes.csv"),
                 row.names = FALSE)

summary_lines <- c(
  "WiCCED simulated GC-MS pipeline summary",
  paste("Mode:", mode),
  paste("Rows in simulated hit table:", nrow(hits)),
  paste("Samples:", length(unique(hits$File.Name))),
  paste("Unique compounds in hits:", length(unique(hits$Compound.Name))),
  paste("Query compounds:", length(query_chemicals)),
  paste("Output directory:", normalizePath(out_dir)),
  "Relative areas sum to 1 within each sample after exact-match aggregation.",
  "All data are simulated teaching data and must not be interpreted as environmental evidence."
)
writeLines(summary_lines, file.path(out_dir, "wicced_simulated_gcms_pipeline_summary.txt"))

message(paste(summary_lines, collapse = "\n"))
