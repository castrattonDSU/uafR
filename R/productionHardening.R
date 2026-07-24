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
      "planPlantChemistryRun(); runPlantChemistryPanel(); runPlantPhytochemistryBatch()",
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
    .api_row("runPlantChemistryPanel", "Plant chemistry", "experimental",
             "Staged production-panel modes and checkpoint contracts are supported; provider coverage may evolve."),
    .api_row("runPlantChemistryProject", "Plant chemistry", "stable",
             "Curated/cached project-bundle handoff workflow is supported."),
    .api_row("buildLotusIndex", "Plant providers", "stable",
             "Flat-export indexing and build-manifest outputs are supported."),
    .api_row("queryLotusIndex", "Plant providers", "stable",
             "Direct species and labeled fallback occurrence output is supported."),
    .api_row("buildNpassIndex", "Plant providers", "experimental",
             "Official-resource local indexing is supported; upstream release layouts may evolve."),
    .api_row("queryNpassIndex", "Plant providers", "experimental",
             "Local source-record lookup is supported; provider-specific fields may expand."),
    .api_row("writePlantChemistryRetryQueue", "Large-run recovery", "stable",
             "Retry queue schema for failed plant/project batches is supported."),
    .api_row("validatePlantChemistryRunManifest", "Large-run recovery", "stable",
             "Run manifest validation outputs are supported."),
    .api_row("rerunFailedPlantQueries", "Large-run recovery", "experimental",
             "Generic retry runner contract is supported; project-specific runners may evolve."),
    .api_row("standardizeCompoundIdentityAudit", "Identity review", "stable",
             "Compound identity audit table schema is supported."),
    .api_row("validateCompoundIdentityAudit", "Identity review", "stable",
             "Identity audit validation outputs are supported."),
    .api_row("exportCompoundIdentityReviewTemplate", "Identity review", "stable",
             "CSV review template schema is supported."),
    .api_row("applyCompoundIdentityReview", "Identity review", "stable",
             "Review replay wrapper preserves original plant identity review behavior."),
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
    source_url = c(
      "https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest",
      "https://www.kegg.jp/kegg/rest/keggapi.html",
      "https://lotus.naturalproducts.net/",
      "https://lotus.naturalproducts.net/",
      "https://www.knapsackfamily.com/KNApSAcK/",
      "https://bidd.group/NPASS/",
      "https://www.ncbi.nlm.nih.gov/books/NBK25501/",
      "https://www.ncbi.nlm.nih.gov/research/pubtator3/api",
      "https://www.fda.gov/industry/fda-data-standards-advisory-board/structured-product-labeling-resources",
      "https://www.femaflavor.org/",
      "https://www.nlm.nih.gov/mesh/meshhome.html"
    ),
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
    fields_extracted = c(
      "CID; names; formula; SMILES; InChIKey; properties; annotations; references",
      "KEGG IDs; names; definitions; equations; pathways; reactions; enzymes; modules; links",
      "taxon names; compound names/IDs; structures; source records; references; available context",
      "matched taxon fields; compound names/IDs; source records; available context",
      "species/metabolite names and IDs; source record fields when available",
      "species source; natural-product IDs; composition/activity context when available",
      "PMID; DOI; title; abstract; publication metadata",
      "PMID/PMC; species mentions; chemical mentions; normalized entity IDs when available",
      "source annotation headings; text; record links through PubChem",
      "source annotation headings; text; record links through PubChem",
      "source annotation headings; text; record links through PubChem"
    ),
    terms_review_status = rep(
      "not_asserted_check_current_provider_terms", 11
    ),
    redistribution_default = rep(
      "do_not_redistribute_provider_database_copies_without_terms_review", 11
    ),
    contract_reviewed_on = rep("2026-07-14", 11),
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
#' @param species_chunk_size Planned species discovery chunk size.
#' @param compound_batch_size Planned compound identity/enrichment batch size.
#' @param max_pubmed_records Planned maximum PubMed records per species.
#'
#' @return List with `Summary`, `PlantQueries`, `InputNameAudit`,
#' `ProviderPlan`, `CacheSummary`, `OutputEstimates`, `ReadinessChecks`,
#' `RunConfiguration`, and `Recommendations`.
#'
#' @export
planPlantChemistryRun = function(plants,
                                 compounds = NULL,
                                 sources = c("lotus", "pubmed", "pubtator"),
                                 cache_dir = NULL,
                                 lotus_index = NULL,
                                 expected_compounds_per_plant = 25,
                                 write_full_pairwise = FALSE,
                                 species_chunk_size = 25,
                                 compound_batch_size = 25,
                                 max_pubmed_records = 25) {
  name_audit = .uaf_plan_plant_name_audit(plants)
  plant_table = .uaf_plan_plants(name_audit)
  species_count = nrow(plant_table)
  parsed_species_count = sum(plant_table$query_status == "parsed_species")
  review_required_name_count = sum(plant_table$review_required)
  blank_input_count = sum(name_audit$query_status == "blank_input")
  duplicate_input_count = sum(name_audit$duplicate_input)
  species_chunk_size = .uaf_plan_positive_integer(species_chunk_size, 25L)
  compound_batch_size = .uaf_plan_positive_integer(compound_batch_size, 25L)
  max_pubmed_records = .uaf_plan_positive_integer(max_pubmed_records, 25L)
  expected_compounds_per_plant = .uaf_plan_positive_integer(
    expected_compounds_per_plant, 25L
  )
  compound_count = .uaf_plan_compound_count(compounds)
  compounds_supplied = compound_count > 0
  if (compound_count < 1) {
    compound_count = max(1L, species_count * expected_compounds_per_plant)
  }
  sources = tolower(.uaf_non_empty(sources))
  lotus_available = .uaf_plan_lotus_available(lotus_index)
  network_requests = vapply(sources, function(provider) {
    if (provider == "lotus" && lotus_available) return(0L)
    if (provider == "npass") return(0L)
    if (provider == "pubmed") return(as.integer(species_count * 2L))
    as.integer(species_count)
  }, integer(1))
  provider_plan = data.frame(
    provider = sources,
    enabled = "Yes",
    estimated_species_queries = species_count,
    estimated_network_request_lower_bound = network_requests,
    recommended_stage = .uaf_plan_provider_stage(sources, lotus_available),
    recommended_mode = .uaf_plan_provider_mode(
      sources, lotus_available, species_count
    ),
    recommended_throttle_seconds = .uaf_plan_provider_throttle(
      sources, lotus_available
    ),
    max_records_per_species = ifelse(sources == "pubmed",
                                     max_pubmed_records, NA_integer_),
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
  membership_rows = species_count * expected_compounds_per_plant
  plant_pairs = if (species_count >= 2) choose(species_count, 2) else 0
  plant_compound_pairs = plant_pairs * expected_compounds_per_plant^2
  discovery_chunks = if (species_count > 0) {
    ceiling(species_count / species_chunk_size)
  } else 0L
  compound_batches = if (compound_count > 0) {
    ceiling(compound_count / compound_batch_size)
  } else 0L
  estimates = data.frame(
    species_count = species_count,
    planned_discovery_chunk_count = discovery_chunks,
    planned_species_chunk_size = species_chunk_size,
    estimated_unique_compounds = compound_count,
    compound_count_basis = ifelse(
      compounds_supplied, "supplied_compound_input",
      "upper_bound_from_species_times_expected_compounds"
    ),
    planned_compound_batch_count = compound_batches,
    planned_compound_batch_size = compound_batch_size,
    estimated_membership_rows = membership_rows,
    estimated_species_pair_rows = plant_pairs,
    estimated_compound_pair_rows = compound_pairs,
    estimated_cross_plant_compound_pair_rows = plant_compound_pairs,
    estimated_plant_compound_pair_rows = plant_compound_pairs,
    estimated_network_request_lower_bound = sum(network_requests),
    write_full_pairwise = .uaf_yes_no(write_full_pairwise),
    approximate_compound_pair_csv_gb =
      round((compound_pairs * 220) / 1024^3, 3),
    approximate_cross_plant_compound_pair_csv_gb =
      round((plant_compound_pairs * 260) / 1024^3, 3),
    approximate_plant_compound_pair_csv_gb =
      round((plant_compound_pairs * 260) / 1024^3, 3),
    stringsAsFactors = FALSE
  )
  readiness = .uaf_plan_readiness_checks(
    species_count = species_count,
    sources = sources,
    cache_dir = cache_dir,
    lotus_available = lotus_available,
    write_full_pairwise = write_full_pairwise,
    species_chunk_size = species_chunk_size,
    compound_batch_size = compound_batch_size,
    parsed_species_count = parsed_species_count,
    review_required_name_count = review_required_name_count,
    blank_input_count = blank_input_count,
    duplicate_input_count = duplicate_input_count
  )
  readiness_status = if (any(readiness$status == "fail")) {
    "not_ready"
  } else if (any(readiness$status == "warn")) {
    "staged_run_required"
  } else {
    "ready"
  }
  run_configuration = .uaf_plan_run_configuration(
    species_chunk_size, compound_batch_size, max_pubmed_records,
    species_count
  )
  rec = .uaf_run_recommendations(species_count, compound_count,
                                 plant_compound_pairs, write_full_pairwise,
                                 lotus_available, sources)
  summary = data.frame(
    input_name_count = nrow(name_audit),
    species_count = species_count,
    parsed_species_count = parsed_species_count,
    review_required_name_count = review_required_name_count,
    blank_input_count = blank_input_count,
    duplicate_input_count = duplicate_input_count,
    source_count = length(sources),
    planned_discovery_chunk_count = discovery_chunks,
    planned_compound_batch_count = compound_batches,
    estimated_provider_queries = species_count * length(sources),
    estimated_network_request_lower_bound = sum(network_requests),
    cache_dir_supplied = .uaf_yes_no(!is.null(cache_dir)),
    lotus_index_supplied = .uaf_yes_no(lotus_available),
    readiness_status = readiness_status,
    runtime_tier = rec$runtime_tier,
    risk_level = rec$risk_level,
    stringsAsFactors = FALSE
  )
  out = list(Summary = summary,
             PlantQueries = plant_table,
             InputNameAudit = name_audit,
             ProviderPlan = provider_plan,
             CacheSummary = cache_summary,
             OutputEstimates = estimates,
             ReadinessChecks = readiness,
             RunConfiguration = run_configuration,
             Recommendations = rec$table)
  class(out) = c("uaf_plant_run_plan", "list")
  out
}

#' Validate a plant chemistry run manifest
#'
#' @description
#' Standardizes and validates a batch/run manifest from large plant chemistry
#' projects. The function is offline-only: it checks table shape, status values,
#' and output-file existence when paths are supplied, but it never reruns failed
#' queries or contacts providers.
#'
#' @param manifest Data frame, CSV path, JSON path, or list containing a manifest
#' table.
#' @param base_dir Optional directory used to resolve relative output paths.
#'
#' @return A list with `Summary`, `TableQuality`, `Issues`, `RetryQueue`, and
#' standardized `Manifest` tables.
#'
#' @export
validatePlantChemistryRunManifest = function(manifest, base_dir = NULL) {
  dat = .uaf_standardize_run_manifest(manifest)
  if (nrow(dat) < 1) {
    issues = data.frame(
      issue_type = "empty_manifest",
      severity = "fail",
      batch_index = NA_integer_,
      message = "The run manifest has no rows.",
      stringsAsFactors = FALSE
    )
    return(list(
      Summary = data.frame(batch_count = 0L, completed_count = 0L,
                           failed_count = 0L, retry_count = 0L,
                           missing_output_count = 0L,
                           validation_status = "fail",
                           stringsAsFactors = FALSE),
      TableQuality = .uaf_run_manifest_quality(dat),
      Issues = issues,
      RetryQueue = .uaf_empty_df(.uaf_retry_queue_cols()),
      Manifest = dat
    ))
  }
  status = tolower(.uaf_first_non_empty_vec(dat$status, "unknown"))
  completed = status %in% c("complete", "completed", "success", "succeeded",
                            "ok", "pass")
  failed = .uaf_manifest_retry_status(status)
  output_path = .uaf_manifest_output_path(dat, base_dir)
  has_output_col = length(.uaf_non_empty(dat$output_file)) > 0 ||
    length(.uaf_non_empty(dat$output_path)) > 0
  output_missing = completed & has_output_col & !file.exists(output_path)

  issue_rows = list()
  if (any(!status %in% .uaf_allowed_manifest_statuses())) {
    bad = which(!status %in% .uaf_allowed_manifest_statuses())
    issue_rows[[length(issue_rows) + 1L]] = data.frame(
      issue_type = "unknown_status",
      severity = "warn",
      batch_index = dat$batch_index[bad],
      message = paste0("Unknown batch status: ", dat$status[bad]),
      stringsAsFactors = FALSE
    )
  }
  if (any(output_missing)) {
    bad = which(output_missing)
    issue_rows[[length(issue_rows) + 1L]] = data.frame(
      issue_type = "completed_output_missing",
      severity = "fail",
      batch_index = dat$batch_index[bad],
      message = paste0("Completed batch output was not found: ",
                       output_path[bad]),
      stringsAsFactors = FALSE
    )
  }
  issues = if (length(issue_rows) > 0) {
    do.call(rbind, issue_rows)
  } else {
    data.frame(issue_type = character(), severity = character(),
               batch_index = integer(), message = character(),
               stringsAsFactors = FALSE)
  }
  retry_queue = writePlantChemistryRetryQueue(dat)
  validation_status = if (any(issues$severity == "fail")) {
    "fail"
  } else if (any(issues$severity == "warn") || nrow(retry_queue) > 0) {
    "warn"
  } else {
    "pass"
  }
  list(
    Summary = data.frame(
      batch_count = nrow(dat),
      completed_count = sum(completed),
      failed_count = sum(failed),
      retry_count = nrow(retry_queue),
      missing_output_count = sum(output_missing),
      validation_status = validation_status,
      stringsAsFactors = FALSE
    ),
    TableQuality = .uaf_run_manifest_quality(dat),
    Issues = issues,
    RetryQueue = retry_queue,
    Manifest = dat
  )
}

#' Write a plant chemistry retry queue
#'
#' @description
#' Extracts failed, incomplete, timed-out, rate-limited, or not-started batches
#' from a large-run manifest into a reproducible retry queue. The queue can be
#' inspected, edited, archived, or passed to `rerunFailedPlantQueries()`.
#'
#' @param manifest Data frame, CSV path, JSON path, or list containing a manifest
#' table.
#' @param path Optional CSV path to write.
#' @param include_status Batch statuses to include in the retry queue.
#' @param overwrite Logical. If `FALSE`, an existing `path` is not replaced.
#'
#' @return Retry queue data frame.
#'
#' @export
writePlantChemistryRetryQueue = function(manifest,
                                         path = NULL,
                                         include_status = c(
                                           "failed", "error", "timeout",
                                           "timed_out", "rate_limited",
                                           "incomplete", "planned",
                                           "not_started", "running",
                                           "started", "stopped", "retry"
                                         ),
                                         overwrite = FALSE) {
  dat = .uaf_standardize_run_manifest(manifest)
  include_status = tolower(.uaf_non_empty(include_status))
  if (nrow(dat) < 1) {
    queue = .uaf_empty_df(.uaf_retry_queue_cols())
  } else {
    status = tolower(.uaf_first_non_empty_vec(dat$status, "unknown"))
    keep = status %in% include_status | .uaf_manifest_retry_status(status)
    dat = dat[keep, , drop = FALSE]
    if (nrow(dat) < 1) {
      queue = .uaf_empty_df(.uaf_retry_queue_cols())
    } else {
      queue = data.frame(
        retry_id = sprintf("retry_%04d", seq_len(nrow(dat))),
        batch_index = dat$batch_index,
        query_start = dat$query_start,
        query_end = dat$query_end,
        query_count = dat$query_count,
        query_label = dat$query_label,
        original_status = dat$status,
        retry_reason = .uaf_retry_reason(dat$status, dat$error_message),
        retry_attempt = suppressWarnings(as.integer(dat$retry_count)) + 1L,
        output_file = dat$output_file,
        output_path = dat$output_path,
        error_message = dat$error_message,
        recommended_action = .uaf_retry_action(dat$status,
                                               dat$error_message),
        created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        stringsAsFactors = FALSE
      )
      queue$retry_attempt[is.na(queue$retry_attempt)] = 1L
    }
  }
  if (!is.null(path)) {
    path = .uaf_first_non_empty_text(path)
    if (is.na(path) || path == "") {
      stop("`path` must be a non-empty CSV path.", call. = FALSE)
    }
    if (file.exists(path) && !isTRUE(overwrite)) {
      stop("Retry queue already exists. Use `overwrite = TRUE`: ", path,
           call. = FALSE)
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(queue, path, row.names = FALSE, na = "")
  }
  queue
}

#' Rerun failed plant chemistry query batches
#'
#' @description
#' Provides a generic retry runner for queue rows produced by
#' `writePlantChemistryRetryQueue()`. By default `dry_run = TRUE`, so the
#' function only returns an execution plan. To actually rerun work, provide a
#' `runner_fun` that accepts one retry-queue row as its first argument and set
#' `dry_run = FALSE`.
#'
#' @param retry_queue Retry queue data frame, CSV path, or manifest-like object.
#' @param runner_fun Optional function called for each retry row when
#' `dry_run = FALSE`.
#' @param dry_run Logical. If `TRUE`, do not execute retries.
#' @param out_dir Optional directory for `retry_run_manifest.csv`.
#' @param overwrite Logical. If `FALSE`, existing retry manifests are preserved.
#' @param ... Additional arguments passed to `runner_fun`.
#'
#' @return List with `Plan`, `RetryQueue`, and `RetryRunManifest`.
#'
#' @export
rerunFailedPlantQueries = function(retry_queue,
                                   runner_fun = NULL,
                                   dry_run = TRUE,
                                   out_dir = NULL,
                                   overwrite = FALSE,
                                   ...) {
  queue = .uaf_standardize_retry_queue(retry_queue)
  plan = data.frame(
    retry_row_count = nrow(queue),
    dry_run = .uaf_yes_no(dry_run),
    runner_supplied = .uaf_yes_no(is.function(runner_fun)),
    status = if (nrow(queue) < 1) "nothing_to_retry" else
      if (isTRUE(dry_run)) "planned" else "executed",
    note = if (isTRUE(dry_run) || !is.function(runner_fun)) {
      "No provider queries were run. Provide `runner_fun` and set `dry_run = FALSE` to execute retries."
    } else {
      "Retry rows were passed to the supplied runner function."
    },
    stringsAsFactors = FALSE
  )
  if (nrow(queue) < 1 || isTRUE(dry_run) || !is.function(runner_fun)) {
    run_manifest = .uaf_empty_df(.uaf_retry_run_manifest_cols())
  } else {
    rows = lapply(seq_len(nrow(queue)), function(i) {
      row = queue[i, , drop = FALSE]
      started = Sys.time()
      result_path = NA_character_
      message = NA_character_
      status = "completed"
      result = tryCatch(
        runner_fun(row, ...),
        error = function(e) {
          status <<- "failed"
          message <<- conditionMessage(e)
          NULL
        }
      )
      if (!is.null(result)) {
        if (is.character(result) && length(result) > 0) {
          result_path = result[[1]]
        } else if (is.list(result) && !is.null(result$output_path)) {
          result_path = as.character(result$output_path[[1]])
        }
      }
      finished = Sys.time()
      data.frame(
        retry_id = row$retry_id,
        batch_index = row$batch_index,
        status = status,
        started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
        finished_at = format(finished, "%Y-%m-%dT%H:%M:%S%z"),
        elapsed_seconds = round(as.numeric(difftime(finished, started,
                                                    units = "secs")), 3),
        output_path = result_path,
        error_message = message,
        stringsAsFactors = FALSE
      )
    })
    run_manifest = do.call(rbind, rows)
  }
  if (!is.null(out_dir)) {
    out_dir = .uaf_first_non_empty_text(out_dir)
    if (is.na(out_dir) || out_dir == "") {
      stop("`out_dir` must be a non-empty directory.", call. = FALSE)
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    manifest_file = file.path(out_dir, "retry_run_manifest.csv")
    if (file.exists(manifest_file) && !isTRUE(overwrite)) {
      stop("Retry run manifest exists. Use `overwrite = TRUE`: ",
           manifest_file, call. = FALSE)
    }
    utils::write.csv(run_manifest, manifest_file, row.names = FALSE, na = "")
  }
  list(Plan = plan, RetryQueue = queue, RetryRunManifest = run_manifest)
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
  primary = c("carbohydrate", "amino_acid", "organic_acid",
              "nucleoside_nucleotide")
  specialized = c(
    "terpenoid", "phenolic_phenylpropanoid", "flavonoid",
    "alkaloid_nitrogenous", "organosulfur", "glucosinolate",
    "benzoxazinoid", "cyanogenic_glycoside", "saponin",
    "steroid_triterpenoid", "polyketide"
  )
  rows = list()
  add_rows = function(scope, family, domain, behavior,
                      comparable = "Yes", review = "No",
                      confidence = "medium") {
    group = if (identical(scope, "volatile_specialized_metabolites")) {
      paste0("volatile_", family)
    } else {
      family
    }
    rows[[length(rows) + 1L]] <<- data.frame(
      comparison_scope = scope,
      comparison_group = group,
      comparison_subgroup = "source_specific_or_family",
      metabolism_domain = domain,
      biosynthetic_family = family,
      chemical_behavior = behavior,
      classification_source = "uafR_comparability_schema",
      classification_confidence = confidence,
      review_required = review,
      comparable_for_matrix = comparable,
      stringsAsFactors = FALSE
    )
  }
  add_rows("primary_metabolites", primary, "primary_metabolism",
           "nonvolatile_or_unspecified")
  add_rows("specialized_metabolites", specialized,
           "specialized_metabolism", "nonvolatile_or_unspecified")
  add_rows("volatile_specialized_metabolites", specialized,
           "specialized_metabolism", "volatile_semivolatile")
  add_rows("lipids_fatty_acids", "fatty_acid_lipid", "lipid_metabolism",
           "nonvolatile_or_unspecified")
  add_rows("plant_hormone_signaling", "plant_hormone_signal",
           "plant_hormone_signaling", "nonvolatile_or_unspecified")
  add_rows("broad_or_uncertain", "broad_or_uncertain",
           "broad_or_uncertain", "nonvolatile_or_unspecified",
           comparable = "No", review = "Yes", confidence = "low")
  add_rows("xenobiotic_or_contaminant", "unknown",
           "xenobiotic_or_contaminant", "unknown",
           comparable = "No", review = "Yes", confidence = "low")
  add_rows("unknown", "unknown", "unknown", "unknown",
           comparable = "No", review = "Yes", confidence = "low")
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
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
  for (col in c("comparison_scope", "comparison_group",
                "metabolism_domain", "biosynthetic_family",
                "chemical_behavior")) {
    x[[col]] = .uaf_classification_alias(x[[col]], col)
  }
  x$classification_confidence = .uaf_classification_alias(
    x$classification_confidence, "classification_confidence"
  )
  x$review_required = .uaf_classification_alias(x$review_required,
                                                 "yes_no")
  x$comparable_for_matrix = .uaf_classification_alias(
    x$comparable_for_matrix, "yes_no"
  )
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
  dictionary = chemistryComparisonDictionary()
  scope_idx = match(x$comparison_scope, dictionary$comparison_scope)
  default_review = dictionary$review_required[scope_idx]
  default_comparable = dictionary$comparable_for_matrix[scope_idx]
  x$review_required = .uaf_first_non_empty_vec(
    x$review_required, default_review, "No"
  )
  x$comparable_for_matrix = .uaf_first_non_empty_vec(
    x$comparable_for_matrix, default_comparable,
    ifelse(tolower(x$comparison_scope) %in%
             c("unknown", "broad_or_uncertain",
               "xenobiotic_or_contaminant"), "No", "Yes")
  )
  .uaf_validate_classification_overrides(x, dictionary)
  x[, cols, drop = FALSE]
}

.uaf_classification_alias = function(x, field) {
  out = tolower(.uaf_squish_text(x))
  out = gsub("[^a-z0-9]+", "_", out)
  out = gsub("^_+|_+$", "", out)
  out[out == ""] = NA_character_
  aliases = switch(
    field,
    comparison_scope = c(
      volatile_specialized = "volatile_specialized_metabolites",
      hormones_signaling = "plant_hormone_signaling",
      hormone_signaling = "plant_hormone_signaling"
    ),
    comparison_group = c(
      terpenoids = "terpenoid",
      phenolics = "phenolic_phenylpropanoid",
      flavonoids = "flavonoid",
      alkaloids = "alkaloid_nitrogenous",
      amino_acids = "amino_acid",
      organic_acids = "organic_acid",
      carbohydrates = "carbohydrate",
      fatty_acids = "fatty_acid_lipid",
      volatile_terpenoids = "volatile_terpenoid"
    ),
    metabolism_domain = c(
      primary = "primary_metabolism",
      specialized = "specialized_metabolism",
      signaling = "plant_hormone_signaling",
      primary_or_storage = "lipid_metabolism",
      external = "xenobiotic_or_contaminant"
    ),
    biosynthetic_family = c(
      terpenoids = "terpenoid",
      phenolics = "phenolic_phenylpropanoid",
      flavonoids = "flavonoid",
      alkaloids = "alkaloid_nitrogenous",
      fatty_acids = "fatty_acid_lipid",
      hormones = "plant_hormone_signal"
    ),
    chemical_behavior = c(
      volatile_or_semivolatile = "volatile_semivolatile",
      volatile = "volatile_semivolatile",
      nonvolatile = "nonvolatile_or_unspecified"
    ),
    classification_confidence = c(),
    yes_no = c(true = "yes", false = "no", y = "yes", n = "no",
               `1` = "yes", `0` = "no"),
    c()
  )
  hit = !is.na(out) & out %in% names(aliases)
  out[hit] = unname(aliases[out[hit]])
  if (identical(field, "yes_no")) {
    out[out == "yes"] = "Yes"
    out[out == "no"] = "No"
  }
  out
}

.uaf_validate_classification_overrides = function(x, dictionary) {
  allowed = list(
    comparison_scope = unique(dictionary$comparison_scope),
    comparison_group = unique(dictionary$comparison_group),
    metabolism_domain = unique(dictionary$metabolism_domain),
    biosynthetic_family = unique(dictionary$biosynthetic_family),
    chemical_behavior = unique(dictionary$chemical_behavior),
    classification_confidence = .plant_confidence_values(),
    review_required = .plant_yes_no_values(),
    comparable_for_matrix = .plant_yes_no_values()
  )
  for (field in names(allowed)) {
    value = x[[field]]
    bad = .uaf_non_empty(value[!is.na(value) & !value %in% allowed[[field]]])
    if (length(bad) > 0) {
      stop("Invalid `", field, "` override value(s): ",
           paste(unique(bad), collapse = ", "), ". Allowed values: ",
           paste(allowed[[field]], collapse = ", "), ".", call. = FALSE)
    }
  }
  complete = !is.na(x$comparison_scope) & !is.na(x$comparison_group)
  if (any(complete)) {
    supplied = paste(x$comparison_scope[complete],
                     x$comparison_group[complete], sep = "\r")
    valid = paste(dictionary$comparison_scope, dictionary$comparison_group,
                  sep = "\r")
    if (any(!supplied %in% valid)) {
      bad = unique(gsub("\r", " / ", supplied[!supplied %in% valid],
                        fixed = TRUE))
      stop("Inconsistent comparison scope/group override pair(s): ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
    }
  }
  noncomparable = x$comparison_scope %in%
    c("unknown", "broad_or_uncertain", "xenobiotic_or_contaminant")
  if (any(noncomparable & x$comparable_for_matrix == "Yes", na.rm = TRUE)) {
    stop("Unknown, broad/uncertain, and xenobiotic/contaminant scopes cannot ",
         "be marked comparable for a default matrix.", call. = FALSE)
  }
  invisible(TRUE)
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

#' Standardize a compound identity audit table
#'
#' @description
#' Converts uafR compound-resolution outputs into a stable audit schema for
#' review, handoff, and reproducible downstream filtering. The audit table does
#' not invent structures or identifiers; unresolved fields remain missing.
#'
#' @param x Plant phytochemistry result, `CompoundResolution` data frame, or
#' compatible identity table.
#' @param occurrences Optional occurrence table used to add query counts when
#' the input lacks them.
#'
#' @return Compound identity audit data frame.
#'
#' @export
standardizeCompoundIdentityAudit = function(x, occurrences = NULL) {
  resolution = .uaf_identity_resolution_from_input(x)
  if (nrow(resolution) < 1) {
    return(.uaf_empty_df(.uaf_identity_audit_cols()))
  }
  if (!is.null(occurrences) && is.data.frame(occurrences) &&
      !"query_count" %in% names(resolution)) {
    occurrences = .plant_normalize_occurrences(occurrences)
    counts = table(occurrences$compound_name_clean)
    resolution$query_count = as.integer(counts[resolution$compound_name_clean])
    resolution$query_count[is.na(resolution$query_count)] = 0L
  }
  rows = lapply(seq_len(nrow(resolution)), function(i) {
    row = resolution[i, , drop = FALSE]
    issue = .plant_identity_issue(row)
    identity = .uaf_identity_values(row)
    flags = .uaf_identity_risk_flags(row, issue)
    review_required = issue$type != "resolved" ||
      issue$type %in% .plant_identity_review_required_types()
    data.frame(
      query_name = row$compound_name,
      query_name_clean = row$compound_name_clean,
      resolved_name = row$compound_name,
      cid = identity$cid,
      inchikey = identity$inchikey,
      inchikey_first_block = .uaf_inchikey_first_block(identity$inchikey),
      smiles = identity$smiles,
      molecular_formula = identity$formula,
      resolution_method = .uaf_identity_resolution_method(row),
      resolution_source = row$resolution_source,
      resolution_confidence = .uaf_identity_confidence(row, issue),
      alias_used = .uaf_identity_alias_used(row),
      ambiguity_flag = .uaf_yes_no(flags$ambiguity_flag),
      review_required = .uaf_yes_no(review_required),
      recommended_action = issue$decision,
      salt_hydrate_flag = .uaf_yes_no(flags$salt_hydrate),
      stereochemistry_unspecified_flag =
        .uaf_yes_no(flags$stereochemistry_unspecified),
      mixture_common_name_flag = .uaf_yes_no(flags$mixture_common_name),
      class_like_name_flag = .uaf_yes_no(flags$class_like_name),
      plant_source_name_flag = .uaf_yes_no(flags$plant_source_name),
      synonym_only_match_flag = .uaf_yes_no(flags$synonym_only_match),
      alias_derived_match_flag = .uaf_yes_no(flags$alias_derived_match),
      multiple_candidate_flag = .uaf_yes_no(flags$multiple_candidate),
      identity_issue_type = issue$type,
      review_reason = issue$reason,
      notes = row$notes,
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  out = out[, .uaf_identity_audit_cols(), drop = FALSE]
  row.names(out) = NULL
  out
}

#' Validate a compound identity audit table
#'
#' @param x Identity audit table or object accepted by
#' `standardizeCompoundIdentityAudit()`.
#'
#' @return List with `Summary`, `TableQuality`, `Issues`, and `Audit`.
#'
#' @export
validateCompoundIdentityAudit = function(x) {
  audit = standardizeCompoundIdentityAudit(x)
  required = .uaf_identity_audit_cols()
  present = required %in% names(audit)
  quality = data.frame(
    required_column = required,
    present = .uaf_yes_no(present),
    stringsAsFactors = FALSE
  )
  issues = list()
  if (any(!present)) {
    issues[[length(issues) + 1L]] = data.frame(
      issue_type = "missing_required_column",
      severity = "fail",
      row_id = NA_integer_,
      message = paste("Missing columns:",
                      paste(required[!present], collapse = ", ")),
      stringsAsFactors = FALSE
    )
  }
  dup = duplicated(audit$query_name_clean) & .bundle_known(audit$query_name_clean)
  if (any(dup)) {
    issues[[length(issues) + 1L]] = data.frame(
      issue_type = "duplicate_query_name_clean",
      severity = "warn",
      row_id = which(dup),
      message = paste0("Duplicate cleaned query name: ",
                       audit$query_name_clean[dup]),
      stringsAsFactors = FALSE
    )
  }
  review = .bundle_truthy(audit$review_required)
  structure_resolved = .bundle_known(audit$smiles) |
    .bundle_known(audit$inchikey) | .bundle_known(audit$cid)
  issues_df = if (length(issues) > 0) {
    do.call(rbind, issues)
  } else {
    data.frame(issue_type = character(), severity = character(),
               row_id = integer(), message = character(),
               stringsAsFactors = FALSE)
  }
  status = if (any(issues_df$severity == "fail")) "fail" else
    if (any(issues_df$severity == "warn") || any(review)) "warn" else "pass"
  list(
    Summary = data.frame(
      identity_count = nrow(audit),
      structure_resolved_count = sum(structure_resolved),
      review_required_count = sum(review),
      ambiguity_flag_count = sum(.bundle_truthy(audit$ambiguity_flag)),
      validation_status = status,
      stringsAsFactors = FALSE
    ),
    TableQuality = quality,
    Issues = issues_df,
    Audit = audit
  )
}

#' Export a compound identity review template
#'
#' @param x Plant phytochemistry result, `CompoundResolution`, or identity
#' audit table.
#' @param path Optional CSV path to write.
#' @param include_resolved Logical. If `FALSE`, resolved/no-review rows are
#' omitted.
#' @param overwrite Logical. If `FALSE`, an existing `path` is not replaced.
#'
#' @return Review template data frame.
#'
#' @export
exportCompoundIdentityReviewTemplate = function(x,
                                                path = NULL,
                                                include_resolved = FALSE,
                                                overwrite = FALSE) {
  if (inherits(x, "uaf_plant_phytochemistry") ||
      (is.list(x) && is.data.frame(x$CompoundResolution)) ||
      .uaf_is_compound_resolution_table(x)) {
    template = plantCompoundIdentityReviewTable(x,
                                                include_resolved =
                                                  include_resolved)
  } else {
    audit = standardizeCompoundIdentityAudit(x)
    keep = rep(TRUE, nrow(audit))
    if (!isTRUE(include_resolved)) {
      keep = .bundle_truthy(audit$review_required) |
        !(.bundle_known(audit$smiles) | .bundle_known(audit$inchikey) |
            .bundle_known(audit$cid))
    }
    audit = audit[keep, , drop = FALSE]
    template = data.frame(
      review_id = sprintf("compound_identity_%04d", seq_len(nrow(audit))),
      review_decision = "needs_review",
      reviewed_by = NA_character_,
      reviewed_at = NA_character_,
      review_note = NA_character_,
      query_name = audit$query_name,
      query_name_clean = audit$query_name_clean,
      cid = audit$cid,
      inchikey = audit$inchikey,
      smiles = audit$smiles,
      molecular_formula = audit$molecular_formula,
      identity_issue_type = audit$identity_issue_type,
      recommended_action = audit$recommended_action,
      proposed_compound_name = NA_character_,
      proposed_cid = NA_character_,
      proposed_inchikey = NA_character_,
      proposed_smiles = NA_character_,
      proposed_molecular_formula = NA_character_,
      proposed_resolution_source = NA_character_,
      review_reason = audit$review_reason,
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(path)) {
    path = .uaf_first_non_empty_text(path)
    if (is.na(path) || path == "") {
      stop("`path` must be a non-empty CSV path.", call. = FALSE)
    }
    if (file.exists(path) && !isTRUE(overwrite)) {
      stop("Review template already exists. Use `overwrite = TRUE`: ", path,
           call. = FALSE)
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(template, path, row.names = FALSE, na = "")
  }
  template
}

#' Apply compound identity review decisions
#'
#' @description
#' User-facing wrapper around `applyPlantCompoundIdentityReview()` for the
#' production identity-review workflow.
#'
#' @param x Plant phytochemistry result, named list with `CompoundResolution`,
#' or `CompoundResolution` table.
#' @param review_table Completed review table.
#' @param reviewer Optional reviewer name.
#' @param require_identity Logical. Require explicit proposed identity fields
#' for update/replace decisions.
#'
#' @return Updated object in the same shape as `x`.
#'
#' @export
applyCompoundIdentityReview = function(x,
                                       review_table,
                                       reviewer = NA_character_,
                                       require_identity = TRUE) {
  applyPlantCompoundIdentityReview(
    x = x,
    review_table = review_table,
    reviewer = reviewer,
    require_identity = require_identity
  )
}

.uaf_is_compound_resolution_table = function(x) {
  is.data.frame(x) &&
    all(.plant_compound_resolution_cols() %in% names(x))
}

.api_row = function(function_name, workflow, stability, contract) {
  data.frame(function_name = function_name,
             workflow = workflow,
             stability = stability,
             contract = contract,
             stringsAsFactors = FALSE)
}

.uaf_schema_version = function() "1.0.0"

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

.uaf_plan_raw_species = function(plants) {
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
  as.character(species)
}

.uaf_plan_plant_name_audit = function(plants) {
  input_name = .uaf_plan_raw_species(plants)
  input_name_clean = .uaf_squish_text(input_name)
  canonical_species = .plant_canonical_taxon_name(input_name_clean)
  blank = is.na(canonical_species) | canonical_species == ""
  species_level = rep(FALSE, length(canonical_species))
  species_level[!blank] = .plant_is_binomial(canonical_species[!blank])
  duplicate_input = rep(FALSE, length(canonical_species))
  duplicate_input[!blank] = duplicated(canonical_species[!blank])
  query_status = ifelse(
    blank, "blank_input",
    ifelse(!species_level, "name_needs_review",
           ifelse(duplicate_input, "duplicate_input", "parsed_species"))
  )
  review_reason = ifelse(
    blank, "Blank plant name; remove or replace before running.",
    ifelse(
      !species_level,
      paste(
        "Not a clear species-level binomial; supply an accepted species name",
        "or document deliberate lower-rank use."
      ),
      ifelse(duplicate_input,
             "Duplicate canonical plant name; only the first query is retained.",
             NA_character_)
    )
  )
  data.frame(
    input_order = seq_along(input_name),
    input_name = input_name,
    canonical_species = canonical_species,
    species_id = .plant_slug(canonical_species),
    species_level_name = species_level,
    duplicate_input = duplicate_input,
    review_required = !blank & !species_level,
    query_status = query_status,
    review_reason = review_reason,
    stringsAsFactors = FALSE
  )
}

.uaf_plan_plants = function(name_audit) {
  keep = name_audit$query_status != "blank_input" &
    !duplicated(name_audit$canonical_species)
  out = name_audit[keep, c("canonical_species", "species_id",
                           "species_level_name", "review_required"),
                   drop = FALSE]
  names(out)[names(out) == "canonical_species"] = "species"
  out$query_status = ifelse(out$species_level_name, "parsed_species",
                            "name_needs_review")
  out$species_level_name = NULL
  row.names(out) = NULL
  out
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

.uaf_plan_positive_integer = function(value, default) {
  if (is.null(value) || length(value) < 1) return(as.integer(default))
  value = suppressWarnings(as.numeric(value[[1]]))
  if (!is.finite(value) || value < 1) return(as.integer(default))
  as.integer(value)
}

.uaf_plan_lotus_available = function(lotus_index) {
  if (is.data.frame(lotus_index)) return(nrow(lotus_index) > 0)
  if (is.list(lotus_index) && !is.data.frame(lotus_index)) {
    return(length(lotus_index) > 0)
  }
  path = .uaf_first_non_empty_text(lotus_index)
  !is.na(path) && path != "" && (file.exists(path) || dir.exists(path))
}

.uaf_plan_provider_stage = function(providers, lotus_available) {
  vapply(providers, function(provider) {
    if (provider == "lotus" && lotus_available) return("1_local_discovery")
    if (provider == "lotus") return("0_pilot_before_large_run")
    if (provider %in% c("pubchem")) return("2_targeted_occurrence_review")
    if (provider %in% c("pubmed", "pubtator")) {
      return("3_targeted_literature_candidates")
    }
    if (provider %in% c("knapsack", "npass")) {
      return("3_provider_specific_optional")
    }
    "3_optional_provider"
  }, character(1))
}

.uaf_plan_provider_mode = function(providers, lotus_available,
                                    species_count) {
  vapply(providers, function(provider) {
    if (provider == "lotus" && lotus_available) return("local_index")
    if (provider == "lotus" && species_count >= 100) {
      return("live_not_recommended_at_scale")
    }
    if (provider == "npass") return("prefetched_or_mocked_only")
    if (provider %in% c("pubmed", "pubtator")) {
      return("cached_targeted_live_stage")
    }
    "cached_pilot_then_stage"
  }, character(1))
}

.uaf_plan_provider_throttle = function(providers, lotus_available) {
  vapply(providers, function(provider) {
    if (provider == "lotus" && lotus_available) return(0)
    if (provider == "pubmed") return(0.34)
    if (provider %in% c("pubtator", "pubchem")) return(0.5)
    1
  }, numeric(1))
}

.uaf_plan_readiness_checks = function(species_count, sources, cache_dir,
                                       lotus_available, write_full_pairwise,
                                       species_chunk_size,
                                       compound_batch_size,
                                       parsed_species_count,
                                       review_required_name_count,
                                       blank_input_count,
                                       duplicate_input_count) {
  large = species_count >= 100
  very_large = species_count >= 500
  live_sources = intersect(sources,
                           c("lotus", "knapsack", "pubchem", "pubmed",
                             "pubtator"))
  lotus_required = "lotus" %in% sources && large
  rows = list(
    data.frame(
      check = "non_empty_unique_species",
      status = ifelse(
        species_count < 1, "fail",
        ifelse(blank_input_count > 0 || duplicate_input_count > 0,
               "warn", "pass")
      ),
      observed = paste0(
        species_count, " unique; ", blank_input_count, " blank; ",
        duplicate_input_count, " duplicate"
      ),
      recommendation = "Resolve blank and duplicate plant names before running.",
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "species_level_name_resolution",
      status = ifelse(review_required_name_count > 0, "fail", "pass"),
      observed = paste0(
        parsed_species_count, " parsed species; ",
        review_required_name_count, " review required"
      ),
      recommendation = paste(
        "Resolve genus-only, sp./spp., and otherwise ambiguous names before",
        "a species-level run. Preserve deliberate lower-rank queries in a",
        "separate, explicitly labeled analysis."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "local_lotus_for_large_panel",
      status = ifelse(!lotus_required || lotus_available, "pass", "fail"),
      observed = .uaf_yes_no(lotus_available),
      recommendation = paste(
        "Use a manifest-backed local LOTUS lookup index for 100+ species;",
        "do not launch a large live LOTUS simple-search run."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "persistent_cache_directory",
      status = ifelse(!large || !is.null(cache_dir), "pass", "fail"),
      observed = .uaf_yes_no(!is.null(cache_dir)),
      recommendation = "Supply a persistent cache_dir outside temporary storage.",
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "live_provider_staging",
      status = ifelse(very_large && length(setdiff(live_sources, "lotus")) > 0,
                      "warn", "pass"),
      observed = paste(live_sources, collapse = "; "),
      recommendation = paste(
        "Run local LOTUS discovery first, then target literature and other live",
        "providers to no-hit, priority, or review-required species."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "npass_live_availability",
      status = ifelse("npass" %in% sources, "warn", "pass"),
      observed = .uaf_yes_no("npass" %in% sources),
      recommendation = paste(
        "NPASS species discovery requires a documented local/prefetched source;",
        "do not count it as a live provider."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "full_pairwise_output_disabled",
      status = ifelse(large && isTRUE(write_full_pairwise), "fail", "pass"),
      observed = .uaf_yes_no(write_full_pairwise),
      recommendation = paste(
        "Generate plant-pair summaries by default; write full compound-pair",
        "files only after inspecting the size estimate."
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "conservative_species_chunk_size",
      status = ifelse(species_chunk_size <= 25, "pass", "warn"),
      observed = as.character(species_chunk_size),
      recommendation = "Use 10-25 species per discovery chunk for large runs.",
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "conservative_compound_batch_size",
      status = ifelse(compound_batch_size <= 50, "pass", "warn"),
      observed = as.character(compound_batch_size),
      recommendation = "Use 20-50 compounds per PubChem identity batch.",
      stringsAsFactors = FALSE
    )
  )
  do.call(rbind, rows)
}

.uaf_plan_run_configuration = function(species_chunk_size,
                                        compound_batch_size,
                                        max_pubmed_records,
                                        species_count) {
  data.frame(
    stage = c(
      "0_preflight", "1_local_discovery", "2_evidence_review",
      "3_identity_resolution", "4_targeted_literature",
      "5_rich_enrichment", "6_similarity_and_export"
    ),
    setting = c(
      "plan_and_validate_names", "species_chunk_size",
      "direct_source_backed_filter", "compound_batch_size",
      "max_pubmed_records", "reviewed_subset_only",
      "summary_first"
    ),
    recommended_value = c(
      paste0(species_count, " unique species"),
      as.character(species_chunk_size),
      "direct species + source-backed + review decisions",
      as.character(compound_batch_size),
      as.character(max_pubmed_records),
      "research/full enrichment after identity audit",
      "plant-pair and comparable-scope summaries"
    ),
    rationale = c(
      "Freeze the input universe and preserve a run signature.",
      "Limits failure scope and makes discovery restartable.",
      "Prevents candidate co-mentions and fallback records from silently entering analysis.",
      "Reduces PubChem service-busy risk and limits retry cost.",
      "Controls literature volume; expand only for priority species.",
      "Avoids expensive enrichment of unresolved or excluded compounds.",
      "Avoids unnecessary full pairwise output explosion."
    ),
    stringsAsFactors = FALSE
  )
}

.uaf_run_recommendations = function(species_count, compound_count,
                                    plant_compound_pairs,
                                    write_full_pairwise,
                                    lotus_available = FALSE,
                                    sources = character()) {
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
    if (species_count >= 100 && "lotus" %in% sources && !lotus_available) {
      "Build or supply a local LOTUS lookup index before discovery."
    },
    "Freeze and validate the plant-name input before creating checkpoints.",
    "Run local/direct species discovery first with compound resolution disabled.",
    "Review fallback, candidate, context, and identity flags before enrichment.",
    "Resolve structures in PubChem batches of 20-50 with persistent cache and resume enabled.",
    "Run PubMed/PubTator and other live providers as targeted follow-up stages.",
    "Write plant-pair summaries by default; full compound-pair files are opt-in.",
    "Validate the completed run manifest and final analysis bundle before modeling."
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

.uaf_run_manifest_cols = function() {
  c("batch_index", "query_start", "query_end", "query_count", "query_label",
    "status", "started_at", "finished_at", "elapsed_seconds",
    "cache_hit_count", "request_count", "retry_count", "error_message",
    "output_file", "output_path")
}

.uaf_retry_queue_cols = function() {
  c("retry_id", "batch_index", "query_start", "query_end", "query_count",
    "query_label", "original_status", "retry_reason", "retry_attempt",
    "output_file", "output_path", "error_message", "recommended_action",
    "created_at")
}

.uaf_retry_run_manifest_cols = function() {
  c("retry_id", "batch_index", "status", "started_at", "finished_at",
    "elapsed_seconds", "output_path", "error_message")
}

.uaf_standardize_run_manifest = function(manifest) {
  dat = .uaf_read_manifest_like(manifest)
  cols = .uaf_run_manifest_cols()
  if (!is.data.frame(dat) || nrow(dat) < 1) return(.uaf_empty_df(cols))
  raw_names = names(dat)
  norm_names = .plant_normalize_column_names(raw_names)
  names(dat) = norm_names
  alias = list(
    batch_index = c("batch_index", "batch", "batch_id", "batch_number",
                    "BatchIndex"),
    query_start = c("query_start", "start_index", "first_query", "range_start"),
    query_end = c("query_end", "end_index", "last_query", "range_end"),
    query_count = c("query_count", "queries", "n_queries", "record_count"),
    query_label = c("query_label", "query", "plant", "species", "batch_label"),
    status = c("status", "batch_status", "run_status"),
    started_at = c("started_at", "start_time", "started", "StartedAt"),
    finished_at = c("finished_at", "end_time", "finished", "FinishedAt"),
    elapsed_seconds = c("elapsed_seconds", "elapsed", "duration_seconds"),
    cache_hit_count = c("cache_hit_count", "cache_hits"),
    request_count = c("request_count", "requests", "query_request_count"),
    retry_count = c("retry_count", "retries", "attempt", "attempt_count"),
    error_message = c("error_message", "error", "message", "warning_message"),
    output_file = c("output_file", "file", "result_file"),
    output_path = c("output_path", "path", "result_path", "rds_path")
  )
  out = as.data.frame(
    stats::setNames(rep(list(rep(NA_character_, nrow(dat))), length(cols)),
                    cols),
    stringsAsFactors = FALSE
  )
  for (col in cols) {
    out[[col]] = .uaf_manifest_col(dat, alias[[col]], nrow(dat))
  }
  out$batch_index = .uaf_integer_or_sequence(out$batch_index, nrow(out))
  for (col in c("query_start", "query_end", "query_count",
                "cache_hit_count", "request_count", "retry_count")) {
    out[[col]] = suppressWarnings(as.integer(out[[col]]))
    out[[col]][is.na(out[[col]])] = 0L
  }
  out$elapsed_seconds = suppressWarnings(as.numeric(out$elapsed_seconds))
  out$elapsed_seconds[is.na(out$elapsed_seconds)] = 0
  out$status = tolower(.uaf_first_non_empty_vec(out$status, "unknown"))
  out$.__raw_column_names = paste(raw_names, collapse = "; ")
  out[, cols, drop = FALSE]
}

.uaf_read_manifest_like = function(x) {
  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    if (grepl("[.]json$", x, ignore.case = TRUE)) {
      obj = jsonlite::fromJSON(x, simplifyDataFrame = TRUE)
      return(.uaf_read_manifest_like(obj))
    }
    return(utils::read.csv(x, stringsAsFactors = FALSE,
                           check.names = FALSE))
  }
  if (is.data.frame(x)) return(x)
  if (is.list(x)) {
    candidates = c("RunManifest", "BatchManifest", "Manifest", "manifest",
                   "Project", "Batches")
    for (name in candidates) {
      if (!is.null(x[[name]]) && is.data.frame(x[[name]])) return(x[[name]])
    }
    df_idx = which(vapply(x, is.data.frame, logical(1)))
    if (length(df_idx) > 0) return(x[[df_idx[[1]]]])
  }
  .uaf_empty_df(character())
}

.uaf_manifest_col = function(dat, candidates, n) {
  candidates = unique(.plant_normalize_column_names(candidates))
  hit = candidates[candidates %in% names(dat)]
  if (length(hit) > 0) return(as.character(dat[[hit[[1]]]]))
  rep(NA_character_, n)
}

.uaf_integer_or_sequence = function(x, n) {
  y = suppressWarnings(as.integer(x))
  if (length(y) != n) y = rep(NA_integer_, n)
  missing = is.na(y) | y < 1
  y[missing] = seq_len(n)[missing]
  y
}

.uaf_allowed_manifest_statuses = function() {
  c("unknown", "planned", "not_started", "running", "started",
    "complete", "completed", "success", "succeeded", "ok", "pass",
    "failed", "error", "timeout", "timed_out", "rate_limited",
    "incomplete", "stopped", "retry", "skipped")
}

.uaf_manifest_retry_status = function(status) {
  status = tolower(.uaf_first_non_empty_vec(status, "unknown"))
  status %in% c("failed", "error", "timeout", "timed_out", "rate_limited",
                "incomplete", "planned", "not_started", "running",
                "started", "stopped", "retry")
}

.uaf_manifest_output_path = function(dat, base_dir = NULL) {
  path = .uaf_first_non_empty_vec(dat$output_path, dat$output_file)
  if (!is.null(base_dir)) {
    rel = !is.na(path) & path != "" & !grepl("^(/|[A-Za-z]:)", path)
    path[rel] = file.path(base_dir, path[rel])
  }
  path
}

.uaf_run_manifest_quality = function(dat) {
  cols = .uaf_run_manifest_cols()
  data.frame(
    required_column = cols,
    present = .uaf_yes_no(cols %in% names(dat)),
    non_empty_count = vapply(cols, function(col) {
      if (!col %in% names(dat)) return(0L)
      length(.uaf_non_empty(dat[[col]]))
    }, integer(1)),
    stringsAsFactors = FALSE
  )
}

.uaf_retry_reason = function(status, error_message) {
  status = tolower(.uaf_first_non_empty_vec(status, "unknown"))
  error_message = .uaf_first_non_empty_vec(error_message, "")
  ifelse(status %in% c("timeout", "timed_out"), "timeout",
         ifelse(status == "rate_limited" |
                  grepl("429|503|busy|rate", error_message,
                        ignore.case = TRUE),
                "rate_limit_or_service_busy",
                ifelse(status %in% c("planned", "not_started", "running",
                                     "started", "incomplete", "stopped"),
                       "incomplete_or_not_started", "provider_or_batch_error")))
}

.uaf_retry_action = function(status, error_message) {
  reason = .uaf_retry_reason(status, error_message)
  ifelse(reason == "rate_limit_or_service_busy",
         "Resume later with a longer throttle/cooldown and reuse the existing cache.",
         ifelse(reason == "timeout",
                "Retry the batch with cache enabled and a smaller batch size if it repeats.",
                ifelse(reason == "incomplete_or_not_started",
                       "Run only this queued batch and preserve successful prior outputs.",
                       "Inspect the error message, fix inputs if needed, then rerun only this batch.")))
}

.uaf_standardize_retry_queue = function(x) {
  dat = .uaf_read_manifest_like(x)
  if (!is.data.frame(dat) || nrow(dat) < 1) {
    return(.uaf_empty_df(.uaf_retry_queue_cols()))
  }
  names(dat) = .plant_normalize_column_names(names(dat))
  if (!"retry_id" %in% names(dat)) {
    dat = writePlantChemistryRetryQueue(dat)
  }
  cols = .uaf_retry_queue_cols()
  for (col in cols) if (!col %in% names(dat)) dat[[col]] = NA_character_
  dat[, cols, drop = FALSE]
}

.uaf_identity_audit_cols = function() {
  c("query_name", "query_name_clean", "resolved_name", "cid", "inchikey",
    "inchikey_first_block", "smiles", "molecular_formula",
    "resolution_method", "resolution_source", "resolution_confidence",
    "alias_used", "ambiguity_flag", "review_required",
    "recommended_action", "salt_hydrate_flag",
    "stereochemistry_unspecified_flag", "mixture_common_name_flag",
    "class_like_name_flag", "plant_source_name_flag",
    "synonym_only_match_flag", "alias_derived_match_flag",
    "multiple_candidate_flag", "identity_issue_type", "review_reason",
    "notes")
}

.uaf_identity_resolution_from_input = function(x) {
  if (inherits(x, "uaf_plant_phytochemistry") &&
      is.data.frame(x$CompoundResolution)) {
    x = x$CompoundResolution
  } else if (is.list(x) && !is.data.frame(x) &&
             is.data.frame(x$CompoundResolution)) {
    x = x$CompoundResolution
  }
  if (!is.data.frame(x)) return(.uaf_empty_table(.plant_compound_resolution_cols()))
  raw = as.data.frame(x, stringsAsFactors = FALSE)
  names(raw) = .uaf_restore_identity_names(.plant_normalize_column_names(names(raw)))
  for (col in .plant_compound_resolution_cols()) {
    if (!col %in% names(raw)) raw[[col]] = NA
  }
  raw$compound_name = .uaf_first_non_empty_vec(raw$compound_name,
                                               raw$query_name,
                                               raw$resolved_name)
  raw$compound_name_clean = .uaf_first_non_empty_vec(
    raw$compound_name_clean,
    raw$query_name_clean,
    .plant_clean_compound(raw$compound_name)
  )
  raw$query_count = suppressWarnings(as.integer(raw$query_count))
  raw$query_count[is.na(raw$query_count)] = 1L
  raw$resolved = .bundle_truthy(raw$resolved) |
    .bundle_known(raw$CID) | .bundle_known(raw$InChIKey) |
    .bundle_known(raw$SMILES) | .bundle_known(raw$MolecularFormula)
  .plant_bind_tables(list(raw), .plant_compound_resolution_cols())
}

.uaf_restore_identity_names = function(x) {
  map = c(cid = "CID", pubchem_cid = "CID", inchikey = "InChIKey",
          inchi_key = "InChIKey", smiles = "SMILES",
          canonical_smiles = "SMILES", isomeric_smiles = "SMILES",
          molecularformula = "MolecularFormula",
          molecular_formula = "MolecularFormula")
  out = x
  hit = out %in% names(map)
  out[hit] = unname(map[out[hit]])
  out
}

.uaf_identity_values = function(row) {
  list(
    cid = .uaf_first_non_empty_text(row$CID),
    inchikey = .uaf_first_non_empty_text(row$InChIKey),
    smiles = .uaf_first_non_empty_text(row$SMILES),
    formula = .uaf_first_non_empty_text(row$MolecularFormula)
  )
}

.uaf_inchikey_first_block = function(x) {
  x = .uaf_first_non_empty_vec(x)
  ifelse(is.na(x) | x == "", NA_character_, sub("-.*$", "", x))
}

.uaf_identity_resolution_method = function(row) {
  source = tolower(.uaf_first_non_empty_text(row$resolution_source, ""))
  if (source == "" || is.na(source)) return("unresolved")
  if (grepl("lotus|source", source)) return("source_backed_structure")
  if (grepl("pubchem|cid", source)) return("pubchem_lookup")
  if (grepl("manual|review", source)) return("manual_review")
  "other_source"
}

.uaf_identity_confidence = function(row, issue) {
  if (!row$resolved %in% TRUE) return("low")
  if (issue$type %in% .plant_identity_review_required_types()) return("low")
  method = .uaf_identity_resolution_method(row)
  if (method %in% c("source_backed_structure", "manual_review")) return("high")
  if (method == "pubchem_lookup") return("medium")
  "medium"
}

.uaf_identity_alias_used = function(row) {
  notes = tolower(.uaf_first_non_empty_text(row$notes, ""))
  .uaf_yes_no(grepl("alias|synonym|transliter", notes))
}

.uaf_identity_risk_flags = function(row, issue) {
  text = tolower(paste(.uaf_first_non_empty_text(row$compound_name, ""),
                       .uaf_first_non_empty_text(row$notes, "")))
  issue_type = tolower(.uaf_first_non_empty_text(issue$type, ""))
  list(
    ambiguity_flag = issue$type != "resolved" ||
      grepl("ambiguous|multiple|mixture|isomer", text),
    salt_hydrate = grepl("salt|hydrate|hydrochloride|sodium|potassium",
                         text),
    stereochemistry_unspecified = .bundle_known(row$SMILES) &&
      !grepl("@|/|\\\\", .uaf_first_non_empty_text(row$SMILES, "")) &&
      grepl("stereo|isomer|alpha|beta|gamma|delta", text),
    mixture_common_name = grepl("mixture|extract|fraction|isomers|derivatives",
                                text) ||
      grepl("mixture|isomer", issue_type),
    class_like_name = grepl("class|family|broad_class", issue_type),
    plant_source_name = grepl("source|product|material", issue_type),
    synonym_only_match = grepl("synonym", text),
    alias_derived_match = grepl("alias|transliter", text),
    multiple_candidate = grepl("multiple|ambiguous", text) ||
      grepl("ambiguous", issue_type)
  )
}
