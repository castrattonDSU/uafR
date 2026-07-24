#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

repo_root = normalizePath(getwd(), mustWork = TRUE)
default_result = file.path(
  repo_root, "lotus_cache", "exports",
  "dsi_lotus_identity_source_resolved_20260611",
  "dsi_lotus_phyto_identity_source_resolved_result.rds"
)
default_fingerprints = file.path(
  repo_root, "lotus_cache", "exports",
  "dsi_plant_chemistry_analysis_bundle_20260624_publication_ready",
  "06_PubChemFingerprints.csv"
)
default_out = file.path(
  repo_root, "lotus_cache", "exports",
  "dsi_tanimoto_server_handoff_source_only_20260715"
)

parse_args = function(args) {
  out = list()
  i = 1L
  while (i <= length(args)) {
    key = args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must use --name value pairs.", call. = FALSE)
    }
    out[[gsub("-", "_", substring(key, 3), fixed = TRUE)]] = args[[i + 1L]]
    i = i + 2L
  }
  out
}

as_flag = function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(trimws(x)) %in% c("true", "yes", "1")
}

old_url_hash = function(url) {
  vapply(url, function(value) {
    if (is.na(value) || !nzchar(value)) return(NA_character_)
    ints = utf8ToInt(value)
    weighted = sum((ints * seq_along(ints)) %% .Machine$integer.max)
    paste0(nchar(value), "_",
           sprintf("%08x", as.integer(weighted %% .Machine$integer.max)))
  }, character(1))
}

clean_text = function(x) {
  x = trimws(as.character(x))
  x[x == "" | toupper(x) %in% c("NA", "N/A", "NULL")] = NA_character_
  x
}

args = parse_args(commandArgs(trailingOnly = TRUE))
plant_result_file = if (is.null(args$plant_result)) default_result else
  args$plant_result
legacy_fingerprint_file = if (is.null(args$legacy_fingerprints)) {
  default_fingerprints
} else args$legacy_fingerprints
out_dir = if (is.null(args$out_dir)) default_out else args$out_dir
overwrite = as_flag(args$overwrite, FALSE)
for (path in c(plant_result_file, legacy_fingerprint_file)) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("This audit requires devtools for the current uafR source tree.")
}
devtools::load_all(repo_root, quiet = TRUE)

result = readRDS(plant_result_file)
required_tables = c("PlantCompoundOccurrences", "CompoundResolution",
                    "CategorateResult", "CompoundIdentityReview")
missing_tables = setdiff(required_tables, names(result))
if (length(missing_tables) > 0) {
  stop("Plant result is missing: ", paste(missing_tables, collapse = ", "))
}
source_identity = result$CategorateResult$SourceCompoundIdentity
if (!is.data.frame(source_identity) || nrow(source_identity) < 1) {
  stop("Plant result does not contain source-backed LOTUS identities.")
}

# Quarantine every name- or PubChem-resolved structure. Exact LOTUS record
# identities are joined separately by preparePlantTanimotoInput().
source_resolution = result$CompoundResolution
source_keep = source_resolution$resolution_source == "LOTUS_source_identity" &
  source_resolution$resolved %in% TRUE
source_resolution$resolved[!source_keep] = FALSE
for (col in c("CID", "InChIKey", "SMILES", "MolecularFormula")) {
  source_resolution[[col]][!source_keep] = NA
}
source_resolution$resolution_source[!source_keep] =
  "excluded_legacy_pubchem_resolution"
source_resolution$notes[!source_keep] = paste(
  "Legacy PubChem-derived identity quarantined; exact source-record identity",
  "is required for this recovery handoff."
)

safe_input = list(
  PlantCompoundOccurrences = result$PlantCompoundOccurrences,
  CompoundResolution = source_resolution,
  SourceCompoundIdentity = source_identity,
  CompoundIdentityReview = result$CompoundIdentityReview
)
prepared = preparePlantTanimotoInput(
  safe_input,
  occurrence_status = c("direct_reported", "curated_reported"),
  analysis_ready = TRUE,
  min_confidence = "medium",
  include_review_required = FALSE,
  out_dir = out_dir,
  overwrite = overwrite,
  strict = TRUE
)

legacy = utils::read.csv(
  legacy_fingerprint_file, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA", "N/A", "NULL"), fileEncoding = "UTF-8"
)
required_legacy = c("compound_id", "compound_name", "pubchem_cid",
                    "InChIKey_source", "InChIKey", "cid_resolution_url")
missing_legacy = setdiff(required_legacy, names(legacy))
if (length(missing_legacy) > 0) {
  stop("Legacy fingerprint table is missing: ",
       paste(missing_legacy, collapse = ", "))
}
source_key = toupper(clean_text(legacy$InChIKey_source))
fetched_key = toupper(clean_text(legacy$InChIKey))
legacy$legacy_cache_hash = old_url_hash(clean_text(legacy$cid_resolution_url))
hash_url_count = ave(
  clean_text(legacy$cid_resolution_url), legacy$legacy_cache_hash,
  FUN = function(x) length(unique(x[!is.na(x)]))
)
legacy$distinct_urls_for_legacy_hash = suppressWarnings(as.integer(hash_url_count))
legacy$legacy_hash_collision = legacy$distinct_urls_for_legacy_hash > 1L
legacy$source_pubchem_inchikey_match = !is.na(source_key) &
  !is.na(fetched_key) & source_key == fetched_key
legacy$identity_integrity_status = ifelse(
  is.na(source_key), "source_inchikey_missing",
  ifelse(is.na(fetched_key), "pubchem_inchikey_missing",
         ifelse(source_key == fetched_key, "verified_exact_inchikey",
                "source_pubchem_inchikey_mismatch"))
)
legacy_audit = legacy[
  legacy$legacy_hash_collision %in% TRUE |
    legacy$identity_integrity_status != "verified_exact_inchikey",
  c("compound_id", "compound_name", "pubchem_cid", "InChIKey_source",
    "InChIKey", "identity_integrity_status", "legacy_cache_hash",
    "distinct_urls_for_legacy_hash", "legacy_hash_collision",
    "cid_resolution_status", "cid_resolution_source", "cid_resolution_url")
]
legacy_audit = legacy_audit[order(
  legacy_audit$legacy_cache_hash,
  legacy_audit$identity_integrity_status,
  legacy_audit$compound_name
), , drop = FALSE]

membership = prepared$PlantCompoundMembership
membership_by_plant = split(membership$compound_id, membership$species)
cross_plant_membership_pair_count = choose(nrow(membership), 2) - sum(
  vapply(membership_by_plant, function(x) choose(length(x), 2), numeric(1))
)

legacy_summary = data.frame(
  legacy_fingerprint_rows = nrow(legacy),
  source_inchikey_count = sum(!is.na(source_key)),
  exact_source_pubchem_inchikey_match_count =
    sum(legacy$identity_integrity_status == "verified_exact_inchikey"),
  source_pubchem_inchikey_mismatch_count =
    sum(legacy$identity_integrity_status ==
          "source_pubchem_inchikey_mismatch"),
  legacy_hash_collision_row_count =
    sum(legacy$legacy_hash_collision %in% TRUE),
  legacy_hash_collision_group_count = length(unique(
    legacy$legacy_cache_hash[legacy$legacy_hash_collision %in% TRUE]
  )),
  safe_input_evidence_rows = prepared$Summary$input_evidence_rows,
  safe_plant_structure_memberships =
    prepared$Summary$plant_compound_membership_count,
  safe_unique_structures = prepared$Summary$unique_structure_count,
  safe_name_structure_ambiguity_rows =
    prepared$Summary$name_structure_ambiguity_count,
  safe_excluded_review_rows = prepared$Summary$excluded_evidence_rows,
  estimated_unique_compound_pair_count =
    choose(prepared$Summary$unique_structure_count, 2),
  estimated_plant_pair_count = choose(
    prepared$Summary$plants_with_tanimoto_ready_compounds, 2
  ),
  estimated_cross_plant_membership_pair_count =
    cross_plant_membership_pair_count,
  safe_validation_status = prepared$Summary$validation_status,
  pairwise_comparisons_run = FALSE,
  stringsAsFactors = FALSE
)

write_csv = function(x, name) {
  uafR:::.plant_atomic_write_csv(x, file.path(out_dir, name))
}
write_csv(source_identity, "source_compound_identity_audit.csv")
write_csv(source_resolution, "compound_resolution_source_only.csv")
write_csv(legacy_audit, "legacy_pubchem_cache_collision_audit.csv")
write_csv(legacy_summary, "legacy_pubchem_cache_collision_summary.csv")

notice = c(
  "# Legacy PubChem/Tanimoto Outputs Are Quarantined",
  "",
  "The legacy DSI PubChem CID and Fingerprint2D assignments in the June 2026",
  "analysis bundle must not be used for scientific analysis. The former cache",
  "key was not collision resistant: distinct request URLs could share one file",
  "name, allowing a valid response for one InChIKey to be reused for another.",
  "",
  paste("This audit found",
        legacy_summary$source_pubchem_inchikey_mismatch_count,
        "source-versus-PubChem InChIKey mismatches across",
        legacy_summary$legacy_hash_collision_group_count,
        "colliding legacy hash groups."),
  "",
  "The source-backed LOTUS occurrence records were not discarded. The safe",
  "handoff rejoined exact LOTUS source-record InChIKeys, removed every inherited",
  "fingerprint, blanked CIDs whenever an InChIKey is available, excluded review",
  "rows, and deduplicated only species-by-structure membership.",
  "Compound labels associated with multiple exact structures remain separate",
  "and are listed in `plant_compound_name_structure_audit.csv`.",
  "",
  "Use `plant_compound_membership_tanimoto_ready.csv` for a new server run with",
  "the collision-safe cache. Do not use the legacy resolved-compound, PubChem",
  "fingerprint, compound-pair, or plant-pair files. No pairwise comparisons were",
  "run while creating this recovery bundle."
)
writeLines(notice, file.path(out_dir, "DO_NOT_REUSE_LEGACY_TANIMOTO.md"),
           useBytes = TRUE)
invisible(file.copy(
  file.path(repo_root, "tools", "run_plant_tanimoto_server.R"),
  file.path(out_dir, "run_plant_tanimoto_server.R"), overwrite = TRUE
))

server_readme = c(
  "# DSI Source-Backed Plant Tanimoto Server Handoff",
  "",
  "This directory is the validated input handoff. It contains no newly",
  "calculated fingerprints or pairwise Tanimoto values.",
  "",
  "Copy the complete directory and install the same validated uafR source",
  "version on the server. Copy its checked release manifest into this",
  "directory as `uafR_release_manifest.json`. Run the offline preflight:",
  "",
  "```sh",
  "Rscript run_plant_tanimoto_server.R \\",
  "  --input plant_compound_membership_tanimoto_ready.csv \\",
  "  --manifest tanimoto_input_export_manifest.csv \\",
  "  --release-manifest uafR_release_manifest.json \\",
  "  --out-dir server_tanimoto_output \\",
  "  --cache-dir pubchem_tanimoto_cache \\",
  "  --mode preflight --full-pairs true --throttle 1.1",
  "```",
  "",
  paste("Expected unique structures:",
        legacy_summary$safe_unique_structures),
  paste("Expected unique compound pairs:",
        legacy_summary$estimated_unique_compound_pair_count),
  paste("Expected plant pairs:",
        legacy_summary$estimated_plant_pair_count),
  paste("Expected cross-plant membership pairs:",
        legacy_summary$estimated_cross_plant_membership_pair_count),
  "",
  "The preflight must exit successfully without contacting PubChem. Inspect",
  "`server_tanimoto_preflight_checks.csv` and",
  "`server_tanimoto_preflight_summary.csv`. Only after every check passes,",
  "run the deterministic smoke gate:",
  "",
  "```sh",
  "UAFR_CONFIRM_SERVER_TANIMOTO=YES Rscript run_plant_tanimoto_server.R \\",
  "  --input plant_compound_membership_tanimoto_ready.csv \\",
  "  --manifest tanimoto_input_export_manifest.csv \\",
  "  --release-manifest uafR_release_manifest.json \\",
  "  --out-dir server_tanimoto_output \\",
  "  --cache-dir pubchem_tanimoto_cache \\",
  "  --mode smoke --smoke-structures 25 --throttle 1.1",
  "```",
  "",
  "Require zero smoke identity mismatches, at least 80% fingerprint usability,",
  "valid pair rows, and a clean cache-hit rerun. Then run `--mode summary` for",
  "all structures. Run `--mode full --full-pairs true` only after the summary",
  "identity and exclusion audits pass. Reuse the same input, release manifest,",
  "output directory, cache directory, and throttle in every command.",
  "",
  "The server runner disables compound-name fallback, verifies source",
  "InChIKeys against PubChem, uses the collision-safe cache, pauses safely after",
  "repeated service-busy responses, and publishes validated pair files from",
  "atomic staging. Preserve this input bundle and its manifests unchanged.",
  "Review `plant_compound_name_structure_audit.csv` before interpreting names.",
  "Do not reuse the quarantined legacy CID, fingerprint, or Tanimoto outputs."
)
writeLines(server_readme, file.path(out_dir, "SERVER_HANDOFF_README.md"),
           useBytes = TRUE)

files = list.files(out_dir, full.names = TRUE, recursive = FALSE)
files = files[file.info(files)$isdir %in% FALSE]
files = files[basename(files) != "recovery_handoff_manifest.csv"]
manifest = data.frame(
  file = basename(files),
  bytes = unname(file.info(files)$size),
  md5 = unname(tools::md5sum(files)),
  stringsAsFactors = FALSE
)
manifest = manifest[order(manifest$file), , drop = FALSE]
write_csv(manifest, "recovery_handoff_manifest.csv")

print(legacy_summary, row.names = FALSE)
cat("Source-only recovery handoff:", normalizePath(out_dir), "\n")
