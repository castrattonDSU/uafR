#' uafR workflow guide
#'
#' @description
#' Returns a compact decision table that maps common research situations to the
#' recommended uafR workflow. Use this as the first stop when deciding which
#' function family to run.
#'
#' @return A data frame with workflow recommendations and stability labels.
#'
#' @export
uafRWorkflowGuide = function() {
  data.frame(
    user_has = c(
      "GC-MS hit table",
      "compound names or PubChem CIDs",
      "plant species list",
      "curated plant-compound table",
      "large plant panel",
      "compound or plant groups for Tanimoto",
      "analysis bundle for downstream modeling",
      "student training package"
    ),
    recommended_workflow = c(
      "Core GC-MS workflow",
      "Research-grade chemical enrichment",
      "Species-first plant phytochemistry",
      "Curated plant chemistry fallback",
      "Large cached plant chemistry run",
      "PubChem Fingerprint2D Tanimoto",
      "Finalized plant chemistry bundle",
      "Student source bundle workflow"
    ),
    primary_functions = c(
      "spreadOut(); mzExacto(); exactoThese()",
      "categorate(detail = 'research'); pubchemProfile(); keggProfile()",
      "resolvePlantPhytochemistry(); validatePlantPhytochemistryResult()",
      "standardizePlantCompoundIntake(); runPlantChemistryProject()",
      "planPlantChemistryRun(); runPlantPhytochemistryBatch(); runPlantChemistryProject()",
      "chemicalTanimotoSimilarity(); plantChemicalTanimotoSimilarity(); plantComparableTanimotoSummary()",
      "finalizePlantChemistryAnalysisBundle(); validatePlantChemistryAnalysisBundle(); exportPlantChemistryFeatureSet()",
      "tools/build_student_bundle.R; tools/student_bundle/*"
    ),
    stable_outputs = c(
      "standard_spread-like tables and exact mass subsets",
      "categorate result tables, validation, workbook exports",
      "PlantCompoundOccurrences, ProviderDiagnostics, SpeciesChemistrySummary",
      "PlantCompoundMembership and finalized analysis bundle tables",
      "batch manifests, review tables, cache-backed result files",
      "CompoundResolution, PubChemFingerprints, plant/group pair summaries",
      "ExportManifest, DataDictionary, enriched membership, feature matrices",
      "source ZIP, install scripts, acceptance tests"
    ),
    stability = c("stable", "stable", "experimental", "stable",
                  "experimental", "stable", "stable", "project_specific"),
    interpretation_caution = c(
      "GC-MS library hits are tentative unless confirmed with standards.",
      "Database annotations are source evidence, not experimental confirmation.",
      "Public plant chemistry records are incomplete and uneven by species.",
      "Curated rows are only as strong as supplied citations and review notes.",
      "Large live-provider runs require throttling, caching, and resume plans.",
      "Structural similarity is not functional equivalence.",
      "Model-ready features are analysis inputs, not efficacy claims.",
      "Training examples may include simulated data and must be labeled."
    ),
    stringsAsFactors = FALSE
  )
}

#' uafR API stability table
#'
#' @description
#' Returns the current public API stability classification for core user-facing
#' functions. Stability labels are documentation commitments for downstream
#' scripts and training material.
#'
#' @return A data frame with function names, workflow areas, and stability.
#'
#' @export
uafRApiStability = function() {
  rows = list(
    .api_row("spreadOut", "Core GC-MS", "stable",
             "Input column requirements and returned spread table structure are supported."),
    .api_row("mzExacto", "Core GC-MS", "stable",
             "Exact-mass matching workflow is supported."),
    .api_row("exactoThese", "Core GC-MS", "stable",
             "Subsetting categorate/exacto candidates is supported."),
    .api_row("categorate", "Chemical enrichment", "stable",
             "`detail = 'research'` result table names are supported."),
    .api_row("pubchemProfile", "Chemical enrichment", "stable",
             "PubChem profile table names are supported; live source content can change."),
    .api_row("keggProfile", "Chemical enrichment", "stable",
             "KEGG profile table names are supported; live source content can change."),
    .api_row("validateCategorateResult", "Validation", "stable",
             "Validation summary/table-quality outputs are supported."),
    .api_row("exportCategorateWorkbook", "Export", "stable",
             "CSV/XLSX export behavior is supported."),
    .api_row("resolvePlantPhytochemistry", "Plant chemistry", "experimental",
             "Species-first discovery is supported but provider coverage may evolve."),
    .api_row("runPlantPhytochemistryBatch", "Plant chemistry", "experimental",
             "Large-run batch workflow is active and cache/resume behavior may expand."),
    .api_row("runPlantChemistryProject", "Plant chemistry", "stable",
             "Curated/cached project-bundle handoff workflow is supported."),
    .api_row("chemicalTanimotoSimilarity", "Tanimoto", "stable",
             "Compound/group Tanimoto output names are supported."),
    .api_row("plantChemicalTanimotoSimilarity", "Tanimoto", "stable",
             "Plant-labeled Tanimoto aliases are supported."),
    .api_row("plantComparableTanimotoSummary", "Tanimoto", "stable",
             "Comparable scope/group filtered summaries are supported."),
    .api_row("exportPlantChemistryFeatureSet", "Model-ready export", "stable",
             "Species feature matrices and manifest outputs are supported."),
    .api_row("tools/run_dsi_categorate_tanimoto.R", "Project wrapper",
             "project_specific", "DSI-specific script; not a general package API.")
  )
  do.call(rbind, rows)
}

#' uafR schema metadata
#'
#' @description
#' Creates package/schema metadata rows for result objects, manifests, and
#' downstream bundles.
#'
#' @param workflow_name Name of the workflow creating the object.
#' @param workflow_parameters Optional named list of important parameters.
#' @param schema_version uafR schema version. Defaults to the current internal
#' schema version.
#'
#' @return One-row data frame with schema and package metadata.
#'
#' @export
uafRSchemaMetadata = function(workflow_name = "unspecified",
                              workflow_parameters = NULL,
                              schema_version = .uaf_schema_version()) {
  data.frame(
    uafR_schema_version = schema_version,
    uafR_package_version = as.character(utils::packageVersion("uafR")),
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    workflow_name = .uaf_first_non_empty_text(workflow_name, "unspecified"),
    workflow_parameters = .uaf_parameter_json(workflow_parameters),
    stringsAsFactors = FALSE
  )
}

#' uafR provider contracts
#'
#' @description
#' Returns a source-contract table describing how uafR interprets public data
#' providers. The table is documentation and audit metadata; it does not query
#' any live services.
#'
#' @param providers Optional providers to return.
#'
#' @return Data frame with provider contract records.
#'
#' @export
uafRProviderContracts = function(providers = NULL) {
  out = data.frame(
    provider = c("PubChem", "KEGG", "LOTUS local index", "LOTUS live search",
                 "KNApSAcK", "NPASS", "PubMed", "PubTator", "FDA/SPL",
                 "FEMA", "MeSH"),
    source_type = c("public API", "public API", "local public export",
                    "public API", "public web/API", "public database",
                    "public API", "public API", "PubChem source annotation",
                    "PubChem source annotation", "PubChem source annotation"),
    query_method = c(
      "PUG-REST and PUG-View by name/CID/InChIKey",
      "KEGG REST find/get/link",
      "standardizeLotusIndex()/queryLotusIndex() over local flat export",
      "bounded simple-search response parsing",
      "organism/metabolite lookup where available",
      "species-source records where available",
      "NCBI E-utilities ESearch/ESummary/EFetch",
      "PubTator entity annotation by PMID/PMC where available",
      "source-filtered PubChem PUG-View annotations",
      "source-filtered PubChem PUG-View annotations",
      "source-filtered PubChem PUG-View annotations"
    ),
    evidence_interpretation = c(
      "compound identity, properties, annotations, source links, and literature context",
      "biochemical pathway/reaction/enzyme context, not pathway activity",
      "source-backed taxon-compound occurrence from local LOTUS records",
      "bounded candidate/direct occurrence depending on returned taxon fields",
      "source-backed species-metabolite evidence when record fields support it",
      "source-backed natural-product/source evidence when record fields support it",
      "literature record discovery; not compound occurrence unless curated",
      "candidate chemical/species co-mentions; review required",
      "regulatory/drug-label source evidence, not exposure or risk assessment",
      "flavor/fragrance source evidence, not measured sample abundance",
      "biomedical vocabulary/source evidence, not clinical mechanism"
    ),
    default_cache = c("Yes", "Yes", "local index", "Yes", "Yes", "Yes",
                      "Yes", "Yes", "via PubChem cache", "via PubChem cache",
                      "via PubChem cache"),
    no_hit_interpretation = rep(
      "No returned record is a coverage limitation, not evidence of biological absence.",
      11
    ),
    diagnostic_expectation = c(
      "request/cache/record/error counts through PubChem/categorate diagnostics",
      "request/cache/record/error counts through KEGG profile diagnostics",
      "local index manifest and row counts",
      "bounded response count and timeout/error status",
      "provider diagnostics row with record/error counts",
      "provider diagnostics row with record/error counts",
      "provider diagnostics row with PMID counts and errors",
      "provider diagnostics row with candidate counts and errors",
      "SourceCoverage and PubChem source annotation rows",
      "SourceCoverage and PubChem source annotation rows",
      "SourceCoverage and PubChem source annotation rows"
    ),
    stability = c("stable", "stable", "stable", "experimental",
                  "experimental", "experimental", "stable", "experimental",
                  "stable", "stable", "stable"),
    stringsAsFactors = FALSE
  )
  if (!is.null(providers)) {
    keep = tolower(out$provider) %in% tolower(.uaf_non_empty(providers))
    out = out[keep, , drop = FALSE]
  }
  row.names(out) = NULL
  out
}

#' Standardize provider diagnostics
#'
#' @description
#' Converts provider diagnostic tables into a stable cross-provider schema.
#'
#' @param x Provider diagnostics data frame.
#'
#' @return Standardized provider diagnostics data frame.
#'
#' @export
standardizeProviderDiagnostics = function(x) {
  cols = .uaf_provider_diagnostic_cols()
  if (!is.data.frame(x) || nrow(x) < 1) return(.uaf_empty_df(cols))
  x = as.data.frame(x, stringsAsFactors = FALSE)
  names(x) = .plant_normalize_column_names(names(x))
  out = data.frame(
    provider = .uaf_diag_col(x, "provider"),
    enabled = .uaf_diag_yes_no(.uaf_diag_col(x, "enabled")),
    queried = .uaf_diag_yes_no(.uaf_diag_col(x, "queried")),
    query_count = .uaf_diag_num(x, c("query_count", "request_count")),
    cache_hit_count = .uaf_diag_num(x, "cache_hit_count"),
    record_count = .uaf_diag_num(x, "record_count"),
    error_count = .uaf_diag_num(x, "error_count"),
    timeout_count = .uaf_diag_num(x, "timeout_count"),
    rate_limit_count = .uaf_diag_num(x, "rate_limit_count"),
    warning_message = .uaf_diag_first(x, c("warning_message", "message",
                                           "error_messages")),
    retrieved_at = .uaf_diag_first(x, c("retrieved_at", "RetrievedAt")),
    status = .uaf_diag_first(x, c("status", "QualityStatus")),
    stringsAsFactors = FALSE
  )
  out$no_hit_reason = ifelse(
    out$record_count == 0 & out$error_count == 0,
    "no_records_found_or_source_not_applicable",
    ifelse(out$error_count > 0, "provider_error_or_timeout", NA_character_)
  )
  out[, cols, drop = FALSE]
}

#' Inspect a uafR cache directory
#'
#' @description
#' Lists cache files and summarizes size, extension, and likely provider from
#' path names. This function never queries live services.
#'
#' @param cache_dir Cache directory.
#' @param recursive Logical. If `TRUE`, inspect nested files.
#'
#' @return Cache inventory data frame.
#'
#' @export
inspectUafRCache = function(cache_dir, recursive = TRUE) {
  if (missing(cache_dir) || length(.uaf_non_empty(cache_dir)) != 1) {
    stop("`cache_dir` must be one directory path.", call. = FALSE)
  }
  cache_dir = normalizePath(cache_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(cache_dir)) {
    return(.uaf_empty_df(c("cache_dir", "file", "relative_path", "provider",
                           "extension", "size_bytes", "modified_at")))
  }
  files = list.files(cache_dir, recursive = recursive, full.names = TRUE,
                     all.files = FALSE, no.. = TRUE)
  files = files[file.exists(files) & !dir.exists(files)]
  if (length(files) < 1) {
    return(.uaf_empty_df(c("cache_dir", "file", "relative_path", "provider",
                           "extension", "size_bytes", "modified_at")))
  }
  info = file.info(files)
  rel = sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", cache_dir),
                   "/?"), "", normalizePath(files, winslash = "/",
                                            mustWork = FALSE))
  data.frame(
    cache_dir = cache_dir,
    file = normalizePath(files, winslash = "/", mustWork = FALSE),
    relative_path = rel,
    provider = .uaf_cache_provider(rel),
    extension = tolower(sub("^.*[.]", "", basename(files))),
    size_bytes = as.numeric(info$size),
    modified_at = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
}

#' Summarize a uafR cache directory
#'
#' @param cache_dir Cache directory.
#' @param recursive Logical. If `TRUE`, inspect nested files.
#'
#' @return Provider-level cache summary.
#'
#' @export
summarizeUafRCache = function(cache_dir, recursive = TRUE) {
  inv = inspectUafRCache(cache_dir, recursive = recursive)
  cols = c("provider", "file_count", "total_size_bytes", "total_size_mb",
           "json_count", "rds_count", "csv_count", "cache_dir")
  if (nrow(inv) < 1) return(.uaf_empty_df(cols))
  rows = lapply(split(inv, inv$provider), function(x) {
    data.frame(
      provider = x$provider[[1]],
      file_count = nrow(x),
      total_size_bytes = sum(x$size_bytes, na.rm = TRUE),
      total_size_mb = round(sum(x$size_bytes, na.rm = TRUE) / 1024^2, 4),
      json_count = sum(x$extension == "json", na.rm = TRUE),
      rds_count = sum(x$extension == "rds", na.rm = TRUE),
      csv_count = sum(x$extension %in% c("csv", "gz"), na.rm = TRUE),
      cache_dir = x$cache_dir[[1]],
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out[, cols, drop = FALSE]
}

#' Plan a large plant chemistry run
#'
#' @description
#' Estimates run size, request burden, cache availability, pairwise output
#' sizes, and recommended conservative defaults before running a large plant
#' chemistry workflow. It does not query public services.
#'
#' @param plants Character vector, data frame, or CSV path of plant names.
#' @param compounds Optional compound table/vector for categorate/Tanimoto
#' estimates.
#' @param sources Provider sources expected for plant discovery.
#' @param cache_dir Optional cache directory to inspect.
#' @param lotus_index Optional local LOTUS index path.
#' @param expected_compounds_per_plant Estimated compounds per plant when no
#' compound table is supplied.
#' @param write_full_pairwise Logical. If `TRUE`, estimate full pairwise file
#' burden.
#'
#' @return List with `Summary`, `ProviderPlan`, `CacheSummary`,
#' `OutputEstimates`, and `Recommendations`.
#'
#' @export
planPlantChemistryRun = function(plants,
                                 compounds = NULL,
                                 sources = c("lotus", "pubmed", "pubtator"),
                                 cache_dir = NULL,
                                 lotus_index = NULL,
                                 expected_compounds_per_plant = 25,
                                 write_full_pairwise = FALSE) {
  plant_table = .uaf_plan_plants(plants)
  compound_count = .uaf_plan_compound_count(compounds)
  if (compound_count < 1) {
    compound_count = max(1L, nrow(plant_table) *
                           as.integer(expected_compounds_per_plant))
  }
  sources = tolower(.uaf_non_empty(sources))
  provider_plan = data.frame(
    provider = sources,
    enabled = "Yes",
    estimated_queries = nrow(plant_table),
    recommended_mode = ifelse(sources == "lotus" &
                                !is.null(lotus_index) &
                                nzchar(.uaf_first_non_empty_text(lotus_index)),
                              "local_index", "cached_live_or_mocked"),
    recommended_throttle_seconds = ifelse(sources %in%
                                            c("pubmed", "pubtator"), 0.34, 0.5),
    caveat = "Provider coverage is incomplete; no-hit results are not biological absence.",
    stringsAsFactors = FALSE
  )
  cache_summary = if (!is.null(cache_dir)) {
    summarizeUafRCache(cache_dir)
  } else {
    .uaf_empty_df(c("provider", "file_count", "total_size_bytes",
                    "total_size_mb", "json_count", "rds_count", "csv_count",
                    "cache_dir"))
  }
  compound_pairs = if (compound_count >= 2) choose(compound_count, 2) else 0
  membership_rows = nrow(plant_table) * expected_compounds_per_plant
  plant_compound_pairs = if (membership_rows >= 2) choose(membership_rows, 2) else 0
  estimates = data.frame(
    species_count = nrow(plant_table),
    estimated_unique_compounds = compound_count,
    estimated_membership_rows = membership_rows,
    estimated_compound_pair_rows = compound_pairs,
    estimated_plant_compound_pair_rows = plant_compound_pairs,
    write_full_pairwise = .uaf_yes_no(write_full_pairwise),
    approximate_compound_pair_csv_gb =
      round((compound_pairs * 220) / 1024^3, 3),
    approximate_plant_compound_pair_csv_gb =
      round((plant_compound_pairs * 260) / 1024^3, 3),
    stringsAsFactors = FALSE
  )
  rec = .uaf_run_recommendations(nrow(plant_table), compound_count,
                                 plant_compound_pairs, write_full_pairwise)
  summary = data.frame(
    species_count = nrow(plant_table),
    source_count = length(sources),
    estimated_provider_queries = nrow(plant_table) * length(sources),
    cache_dir_supplied = .uaf_yes_no(!is.null(cache_dir)),
    lotus_index_supplied = .uaf_yes_no(!is.null(lotus_index) &&
                                         nzchar(.uaf_first_non_empty_text(lotus_index))),
    runtime_tier = rec$runtime_tier,
    risk_level = rec$risk_level,
    stringsAsFactors = FALSE
  )
  out = list(Summary = summary,
             PlantQueries = plant_table,
             ProviderPlan = provider_plan,
             CacheSummary = cache_summary,
             OutputEstimates = estimates,
             Recommendations = rec$table)
  class(out) = c("uaf_plant_run_plan", "list")
  out
}

#' Evidence grade dictionary
#'
#' @return Data frame defining canonical plant occurrence evidence grades.
#'
#' @export
plantOccurrenceEvidenceDictionary = function() {
  data.frame(
    evidence_grade = c("direct_species_database_record",
                       "direct_species_literature_supported_record",
                       "source_backed_genus_family_fallback",
                       "pubtator_pubmed_candidate_only",
                       "unresolved_or_review_required",
                       "excluded_by_review"),
    evidence_grade_rank = c(1L, 2L, 3L, 4L, 5L, 99L),
    definition = c(
      "A source database record reports the compound at the queried species level.",
      "A literature-supported record reports the compound at the queried species level.",
      "A source-backed record applies at genus or family level rather than direct species level.",
      "A PubMed or PubTator co-mention/candidate record that requires human review.",
      "The compound/evidence is unresolved or too weak for analysis without review.",
      "A reviewer excluded the row from downstream use."
    ),
    allowed_downstream_use = c(
      "Use for source-backed species-level exploratory analysis when identity and comparability filters pass.",
      "Use for source-backed species-level exploratory analysis; cite and inspect the source.",
      "Use as broader taxonomic context only; keep separate from direct species evidence.",
      "Use for literature triage only until curated.",
      "Do not use for analysis before review.",
      "Do not use unless review decision changes."
    ),
    manuscript_safe_wording = c(
      "Public database records report this compound for the species.",
      "Literature-linked records report this compound for the species.",
      "Related taxonomic records report this compound; species-level evidence was not established here.",
      "Literature mining identified a candidate co-mention requiring curation.",
      "The row was retained for audit but not interpreted as resolved occurrence evidence.",
      "The row was excluded by review."
    ),
    reviewer_challenge_question = c(
      "Does the source record clearly resolve the species and compound?",
      "Does the cited paper/source directly report occurrence rather than mention?",
      "Is fallback evidence being separated from direct species evidence?",
      "Was the co-mention manually verified before use?",
      "Why was this row included in any analysis?",
      "Is there a documented reason to reverse the exclusion?"
    ),
    stringsAsFactors = FALSE
  )
}

#' Filter direct plant evidence
#'
#' @param x Occurrence or enriched membership table.
#' @param require_structure Logical. If `TRUE`, require SMILES, InChIKey, or CID.
#'
#' @return Filtered data frame.
#'
#' @export
filterPlantEvidenceDirect = function(x, require_structure = TRUE) {
  .filter_plant_evidence(x,
                         grades = c("direct_species_database_record",
                                    "direct_species_literature_supported_record"),
                         require_structure = require_structure,
                         require_comparable = FALSE,
                         review_required = FALSE)
}

#' Filter comparable plant evidence
#'
#' @param x Occurrence or enriched membership table.
#' @param require_direct Logical. If `TRUE`, retain only direct species grades.
#'
#' @return Filtered data frame.
#'
#' @export
filterPlantEvidenceComparable = function(x, require_direct = TRUE) {
  grades = if (isTRUE(require_direct)) {
    c("direct_species_database_record",
      "direct_species_literature_supported_record")
  } else {
    NULL
  }
  .filter_plant_evidence(x, grades = grades, require_structure = TRUE,
                         require_comparable = TRUE, review_required = FALSE)
}

#' Filter plant evidence requiring review
#'
#' @param x Occurrence or enriched membership table.
#'
#' @return Filtered data frame.
#'
#' @export
filterPlantEvidenceReviewRequired = function(x) {
  .filter_plant_evidence(x, grades = NULL, require_structure = FALSE,
                         require_comparable = FALSE, review_required = TRUE)
}

#' Chemistry comparison dictionary
#'
#' @return Data frame defining default comparable chemistry scopes.
#'
#' @export
chemistryComparisonDictionary = function() {
  data.frame(
    comparison_scope = c("primary_metabolites", "specialized_metabolites",
                         "volatile_specialized", "lipids_fatty_acids",
                         "hormones_signaling", "xenobiotic_or_contaminant",
                         "unknown"),
    comparison_group = c("primary_metabolism", "specialized_metabolism",
                         "volatile_fraction", "lipid_metabolism",
                         "plant_signaling", "external_or_contaminant",
                         "unknown"),
    comparison_subgroup = c("sugars_amino_acids_organic_acids",
                            "phenolics_terpenoids_alkaloids_and_related",
                            "volatile_terpenoids_aromatics_and_aliphatics",
                            "fatty_acids_lipids_waxes",
                            "hormones_and_signaling_molecules",
                            "xenobiotic_or_environmental_compounds",
                            "unclassified"),
    metabolism_domain = c("primary", "specialized", "specialized",
                          "primary_or_storage", "signaling", "external",
                          "unknown"),
    biosynthetic_family = c("central_metabolism", "mixed_specialized",
                            "volatile_specialized", "lipid",
                            "hormone_signal", "xenobiotic", "unknown"),
    chemical_behavior = c("polar_metabolite", "specialized_metabolite",
                          "volatile_or_semivolatile", "lipophilic",
                          "bioactive_signal", "external_context",
                          "unknown"),
    classification_source = c(rep("uafR_default_dictionary", 7)),
    classification_confidence = c("medium", "medium", "medium", "medium",
                                  "medium", "low", "low"),
    review_required = c("No", "No", "No", "No", "No", "Yes", "Yes"),
    comparable_for_matrix = c("Yes", "Yes", "Yes", "Yes", "Yes", "No", "No"),
    stringsAsFactors = FALSE
  )
}

#' Standardize chemistry classification overrides
#'
#' @param x Data frame or CSV path with override rows.
#'
#' @return Standardized override table.
#'
#' @export
standardizeChemistryClassificationOverrides = function(x) {
  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    x = utils::read.csv(x, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame or CSV path.", call. = FALSE)
  }
  names(x) = .plant_normalize_column_names(names(x))
  if (!any(c("compound_id", "compound_name_clean", "compound_name") %in%
           names(x))) {
    stop("Overrides require `compound_id`, `compound_name_clean`, or ",
         "`compound_name`.", call. = FALSE)
  }
  cols = c("compound_id", "compound_name", "compound_name_clean",
           "comparison_scope", "comparison_group", "comparison_subgroup",
           "metabolism_domain", "biosynthetic_family", "chemical_behavior",
           "classification_source", "classification_confidence",
           "review_required", "comparable_for_matrix", "override_note")
  for (col in cols) if (!col %in% names(x)) x[[col]] = NA_character_
  x$compound_name_clean = .uaf_first_non_empty_vec(
    x$compound_name_clean,
    if ("compound_name" %in% names(x)) .plant_clean_compound(x$compound_name) else
      NA_character_
  )
  x$classification_source = .uaf_first_non_empty_vec(
    x$classification_source, "user_override"
  )
  x$classification_confidence = .uaf_first_non_empty_vec(
    x$classification_confidence, "high"
  )
  x$review_required = .uaf_first_non_empty_vec(x$review_required, "No")
  x$comparable_for_matrix = .uaf_first_non_empty_vec(
    x$comparable_for_matrix,
    ifelse(tolower(x$comparison_scope) == "unknown", "No", "Yes")
  )
  x[, cols, drop = FALSE]
}

#' Apply chemistry classification overrides
#'
#' @param comparability Chemistry comparability or enriched membership table.
#' @param overrides Override data frame or CSV path.
#'
#' @return Updated comparability table.
#'
#' @export
applyChemistryClassificationOverrides = function(comparability, overrides) {
  if (!is.data.frame(comparability)) {
    stop("`comparability` must be a data frame.", call. = FALSE)
  }
  out = as.data.frame(comparability, stringsAsFactors = FALSE)
  overrides = standardizeChemistryClassificationOverrides(overrides)
  if (nrow(out) < 1 || nrow(overrides) < 1) return(out)
  for (col in names(overrides)) if (!col %in% names(out)) out[[col]] = NA_character_
  keys = list(
    compound_id = match(out$compound_id, overrides$compound_id),
    compound_name_clean = match(.uaf_first_non_empty_vec(
      out$compound_name_clean,
      if ("compound_name" %in% names(out)) .plant_clean_compound(out$compound_name) else
        NA_character_
    ), overrides$compound_name_clean)
  )
  idx = keys$compound_id
  missing = is.na(idx)
  idx[missing] = keys$compound_name_clean[missing]
  override_cols = c("comparison_scope", "comparison_group",
                    "comparison_subgroup", "metabolism_domain",
                    "biosynthetic_family", "chemical_behavior",
                    "classification_source", "classification_confidence",
                    "review_required", "comparable_for_matrix")
  hit = !is.na(idx)
  for (col in override_cols) {
    out[[col]][hit] = .uaf_first_non_empty_vec(overrides[[col]][idx[hit]],
                                                out[[col]][hit])
  }
  if ("override_note" %in% names(out)) {
    out$override_note[hit] = .uaf_first_non_empty_vec(
      overrides$override_note[idx[hit]], out$override_note[hit]
    )
  }
  out
}

#' uafR claim guidance
#'
#' @description
#' Provides manuscript-safe wording and claim limits for common uafR evidence
#' types.
#'
#' @param evidence_type Optional evidence type(s) to return.
#'
#' @return Data frame with scientific guardrail language.
#'
#' @export
uafRClaimGuidance = function(evidence_type = NULL) {
  out = data.frame(
    evidence_type = c("tentative_gcms_match", "database_occurrence",
                      "source_annotation", "genus_family_fallback",
                      "pubtator_pubmed_candidate", "pathway_context",
                      "hazard_annotation", "tanimoto_similarity",
                      "model_ready_feature"),
    can_support = c(
      "A tentative annotation or candidate identity for follow-up.",
      "A public source reports a compound/taxon or compound/source association.",
      "A source reports a descriptor, classification, or annotation.",
      "Broader taxonomic context when species evidence is absent.",
      "Literature triage candidates requiring human review.",
      "Potential biochemical context in source databases.",
      "Screening-level hazard or safety annotation.",
      "Structural similarity based on PubChem Fingerprint2D.",
      "A reproducible analysis input feature."
    ),
    cannot_support = c(
      "Confirmed chemical identity without authentic standards or stronger evidence.",
      "Presence in the user's sample or complete metabolome coverage.",
      "Experimental confirmation or mechanism by itself.",
      "Direct species occurrence unless separately curated.",
      "Confirmed occurrence or mechanism without source verification.",
      "Pathway activity in the organism or sample.",
      "Risk assessment, exposure, or toxicity in the study system.",
      "Functional equivalence, shared bioactivity, or shared ecological role.",
      "Causality, efficacy, remediation performance, or biological mechanism."
    ),
    recommended_safe_wording = c(
      "The feature was tentatively annotated as...",
      "Public database records report...",
      "Source annotations describe...",
      "Genus/family-level records suggest broader context...",
      "Literature mining identified candidate mentions of...",
      "KEGG/source records provide pathway context for...",
      "Hazard annotations indicate screening-level safety considerations...",
      "Compounds/plants showed structural similarity by PubChem Fingerprint2D...",
      "Features were used as source-backed covariates for downstream analysis..."
    ),
    review_required = c("Yes", "No", "No", "Yes", "Yes", "Yes", "Yes",
                        "No", "No"),
    stringsAsFactors = FALSE
  )
  if (!is.null(evidence_type)) {
    out = out[tolower(out$evidence_type) %in%
                tolower(.uaf_non_empty(evidence_type)), , drop = FALSE]
  }
  row.names(out) = NULL
  out
}

#' Estimate Tanimoto output size
#'
#' @param compound_count Number of unique structure-resolved compounds.
#' @param membership_count Optional plant/group-compound membership rows.
#'
#' @return Data frame with pair-count and approximate CSV-size estimates.
#'
#' @export
estimateTanimotoOutput = function(compound_count, membership_count = NULL) {
  compound_count = suppressWarnings(as.numeric(compound_count[[1]]))
  if (is.na(compound_count) || compound_count < 0) compound_count = 0
  membership_count = if (is.null(membership_count)) NA_real_ else
    suppressWarnings(as.numeric(membership_count[[1]]))
  compound_pairs = if (compound_count >= 2) choose(compound_count, 2) else 0
  membership_pairs = if (!is.na(membership_count) && membership_count >= 2) {
    choose(membership_count, 2)
  } else {
    NA_real_
  }
  data.frame(
    compound_count = compound_count,
    compound_pair_rows = compound_pairs,
    approximate_compound_pair_csv_gb = round((compound_pairs * 220) / 1024^3, 3),
    membership_count = membership_count,
    plant_compound_pair_rows = membership_pairs,
    approximate_plant_compound_pair_csv_gb =
      round((membership_pairs * 260) / 1024^3, 3),
    recommendation = if (compound_pairs > 1000000 ||
                         (!is.na(membership_pairs) &&
                          membership_pairs > 1000000)) {
      "Stream full pairwise outputs to compressed CSV and keep summaries for routine analysis."
    } else {
      "In-memory pairwise output is likely manageable, but summaries remain preferred for handoff."
    },
    stringsAsFactors = FALSE
  )
}

.api_row = function(function_name, workflow, stability, contract) {
  data.frame(function_name = function_name,
             workflow = workflow,
             stability = stability,
             contract = contract,
             stringsAsFactors = FALSE)
}

.uaf_schema_version = function() "2026-06-25"

.uaf_parameter_json = function(x) {
  if (is.null(x)) return(NA_character_)
  tryCatch(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"),
           error = function(e) NA_character_)
}

.uaf_provider_diagnostic_cols = function() {
  c("provider", "enabled", "queried", "query_count", "cache_hit_count",
    "record_count", "error_count", "timeout_count", "rate_limit_count",
    "warning_message", "retrieved_at", "status", "no_hit_reason")
}

.uaf_empty_df = function(cols) {
  as.data.frame(stats::setNames(rep(list(character()), length(cols)), cols),
                stringsAsFactors = FALSE)
}

.uaf_diag_col = function(x, col) {
  if (col %in% names(x)) as.character(x[[col]]) else rep(NA_character_, nrow(x))
}

.uaf_diag_first = function(x, cols) {
  out = rep(NA_character_, nrow(x))
  for (col in cols) {
    col_norm = .plant_normalize_column_names(col)
    if (col_norm %in% names(x)) {
      out = .uaf_first_non_empty_vec(out, as.character(x[[col_norm]]))
    } else if (col %in% names(x)) {
      out = .uaf_first_non_empty_vec(out, as.character(x[[col]]))
    }
  }
  out
}

.uaf_diag_num = function(x, cols) {
  val = suppressWarnings(as.numeric(.uaf_diag_first(x, cols)))
  val[is.na(val)] = 0
  val
}

.uaf_diag_yes_no = function(x) {
  x = tolower(.uaf_first_non_empty_vec(x, "No"))
  ifelse(x %in% c("yes", "true", "t", "1", "enabled"), "Yes", "No")
}

.uaf_first_non_empty_vec = function(...) {
  values = list(...)
  if (length(values) < 1) return(character())
  n = max(vapply(values, length, integer(1)))
  values = lapply(values, function(x) {
    x = as.character(x)
    if (length(x) == 1 && n > 1) x = rep(x, n)
    if (length(x) < n) x = c(x, rep(NA_character_, n - length(x)))
    x = .uaf_squish_text(x)
    x
  })
  out = rep(NA_character_, n)
  for (value in values) {
    hit = (is.na(out) | out == "") & !is.na(value) & value != ""
    out[hit] = value[hit]
  }
  out
}

.uaf_cache_provider = function(path) {
  path = tolower(path)
  out = rep("unknown", length(path))
  out[grepl("pubchem", path)] = "PubChem"
  out[grepl("kegg", path)] = "KEGG"
  out[grepl("lotus", path)] = "LOTUS"
  out[grepl("pubmed|ncbi", path)] = "PubMed"
  out[grepl("pubtator", path)] = "PubTator"
  out[grepl("categorate", path)] = "categorate"
  out[grepl("tanimoto", path)] = "Tanimoto"
  out
}

.uaf_plan_plants = function(plants) {
  if (is.character(plants) && length(plants) == 1 && file.exists(plants)) {
    plants = utils::read.csv(plants, stringsAsFactors = FALSE,
                             check.names = FALSE)
  }
  if (is.data.frame(plants)) {
    nm = .plant_normalize_column_names(names(plants))
    names(plants) = nm
    col = if ("species" %in% names(plants)) "species" else names(plants)[[1]]
    species = plants[[col]]
  } else {
    species = plants
  }
  species = unique(.uaf_non_empty(as.character(species)))
  data.frame(species = species,
             species_id = tolower(gsub("[^a-z0-9]+", "_", species)),
             stringsAsFactors = FALSE)
}

.uaf_plan_compound_count = function(compounds) {
  if (is.null(compounds)) return(0L)
  if (is.character(compounds) && length(compounds) == 1 && file.exists(compounds)) {
    compounds = utils::read.csv(compounds, stringsAsFactors = FALSE,
                                check.names = FALSE)
  }
  if (is.data.frame(compounds)) {
    cols = .plant_normalize_column_names(names(compounds))
    names(compounds) = cols
    key_col = c("compound_id", "cid", "pubchem_cid", "compound_name")
    key_col = key_col[key_col %in% names(compounds)]
    if (length(key_col) > 0) {
      return(length(unique(.uaf_non_empty(compounds[[key_col[[1]]]]))))
    }
    return(nrow(compounds))
  }
  length(unique(.uaf_non_empty(as.character(compounds))))
}

.uaf_run_recommendations = function(species_count, compound_count,
                                    plant_compound_pairs,
                                    write_full_pairwise) {
  runtime_tier = if (species_count >= 500 || compound_count >= 10000) {
    "very_large"
  } else if (species_count >= 100 || compound_count >= 1000) {
    "large"
  } else {
    "small_or_moderate"
  }
  risk = if (isTRUE(write_full_pairwise) && plant_compound_pairs > 1000000) {
    "high_pairwise_output_risk"
  } else if (runtime_tier == "very_large") {
    "high_provider_runtime_risk"
  } else {
    "manageable_with_cache"
  }
  recs = c(
    "Use local LOTUS lookup indexes for medium and large plant panels.",
    "Run live providers only with cache enabled and conservative throttling.",
    "Run categorate/PubChem enrichment in resumable batches.",
    "Write full plant-compound pair Tanimoto files only when explicitly needed.",
    "Validate finalized bundles before downstream modeling."
  )
  list(runtime_tier = runtime_tier,
       risk_level = risk,
       table = data.frame(priority = seq_along(recs),
                          recommendation = recs,
                          stringsAsFactors = FALSE))
}

.filter_plant_evidence = function(x, grades = NULL, require_structure = FALSE,
                                  require_comparable = FALSE,
                                  review_required = FALSE) {
  if (!is.data.frame(x)) stop("`x` must be a data frame.", call. = FALSE)
  dat = as.data.frame(x, stringsAsFactors = FALSE)
  grades_tbl = plantOccurrenceEvidenceGrade(dat)
  keep = rep(TRUE, nrow(dat))
  if (!is.null(grades)) keep = keep & grades_tbl$evidence_grade %in% grades
  if (isTRUE(require_structure)) {
    keep = keep & .bundle_truthy(grades_tbl$structure_resolved)
  }
  if (isTRUE(require_comparable)) {
    keep = keep & .bundle_truthy(grades_tbl$comparable_for_analysis)
  }
  if (isTRUE(review_required)) {
    keep = .bundle_truthy(grades_tbl$review_required)
  } else {
    keep = keep & !.bundle_truthy(grades_tbl$review_required)
  }
  out = dat[keep, , drop = FALSE]
  if (!"evidence_grade" %in% names(out)) {
    out$evidence_grade = grades_tbl$evidence_grade[keep]
  }
  row.names(out) = NULL
  out
}
