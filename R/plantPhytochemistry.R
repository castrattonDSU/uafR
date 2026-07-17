#' Species-first plant phytochemistry resolver
#'
#' @description
#' `resolvePlantPhytochemistry()` builds a species-first phytochemistry result
#' from plant names, public-provider adapters, literature candidates, and
#' optional curated species-compound rows. The resolver is intentionally
#' conservative: no compound is fabricated, literature co-mentions remain
#' candidate evidence, and provider no-hit states are reported through
#' diagnostics.
#'
#' @param plants Character vector or data frame of plant names. If a data frame
#' is supplied, a `species` column is preferred; otherwise the first column is
#' treated as the plant name. Optional `genus` and `family` columns are used
#' when present.
#' @param sources Character vector of provider names. Supported names are
#' `"lotus"`, `"knapsack"`, `"npass"`, `"pubchem"`, `"pubmed"`, and
#' `"pubtator"`.
#' @param taxon_fallback Character vector of fallback ranks to record in query
#' planning. Species-level evidence is queried first. Providers that support
#' broader lookup may also use genus or family terms; fallback evidence is
#' labeled separately from direct species evidence.
#' @param enrich_compounds Logical. If `TRUE`, attempt compound enrichment after
#' occurrence discovery.
#' @param chemical_library Optional chemical library passed to `categorate()`
#' when compound enrichment is requested. If omitted, compound enrichment falls
#' back to `pubchemProfile()` and returns PubChem-derived categorate-like
#' analysis tables without FMCS library matching.
#' @param detail One of `"research"`, `"full"`, or `"none"`.
#' @param cache Logical. If `TRUE`, provider requests may use cache files.
#' @param cache_dir Cache directory. Defaults to a uafR cache directory when
#' needed.
#' @param lotus_index Optional local LOTUS index as a data frame, path to a
#' CSV, TSV, JSON, JSONL, NDJSON, or RDS file, or manifest-backed lookup
#' directory produced by `tools/flatten_lotus_mongo_dump.py --lookup-dir`. When
#' supplied and `"lotus"` is enabled, uafR queries this local index instead of
#' the live LOTUS simple API. Lookup directories are the recommended path for
#' hundreds of plant names.
#' @param provider_indexes Optional named list of local provider indexes.
#' `provider_indexes$lotus` and `provider_indexes$npass` are currently
#' supported. `lotus_index` remains a backward-compatible alias.
#' @param require_all_providers Logical. If `TRUE`, fail when a requested local
#' resource is unavailable or a provider finishes with an incomplete/error
#' status. Use this for audited production panels after a live pilot.
#' @param throttle Seconds to wait between uncached provider requests.
#' @param ncbi_email Optional NCBI email. Defaults to `Sys.getenv("NCBI_EMAIL")`.
#' @param ncbi_tool Optional NCBI tool name. Defaults to `Sys.getenv("NCBI_TOOL",
#' "uafR")`.
#' @param ncbi_api_key Optional NCBI API key. Defaults to
#' `Sys.getenv("NCBI_API_KEY")`. It is used in request URLs but is not persisted
#' in result tables.
#' @param max_pubmed_records Maximum PubMed records per plant to request.
#' @param max_provider_records Maximum occurrence rows to retain per plant from
#' non-PubMed provider adapters.
#' @param request_timeout Maximum seconds allowed for an uncached provider
#' request before the underlying R connection times out. Cached and mocked
#' requests are not delayed by this setting.
#' @param enrich_context Logical. If `TRUE`, attempt source-backed PubMed
#' context enrichment for occurrence rows that have PMID or DOI provenance.
#' This can add plant-part, tissue, or method context only when source text
#' contains the relevant terms; it does not fabricate missing context.
#' @param context_sources Optional source text table, such as
#' `LiteratureCandidates` or a data frame with `pmid`, `doi`, `title`,
#' `abstract`, and `evidence_text`, used for source-backed context enrichment
#' without live requests.
#' @param max_context_sources Maximum unique PMID/DOI records to fetch when
#' `enrich_context = TRUE`.
#' @param min_confidence Minimum confidence for summary/matrix features.
#' @param enrichment_batch_size Maximum number of unique compounds per
#' PubChem-only enrichment batch. Use `Inf` for one batch.
#' @param resume_enrichment Logical. If `TRUE`, reuse saved PubChem-only
#' enrichment batch files when available.
#' @param progress Logical. If `TRUE`, print simple enrichment progress
#' messages.
#' @param defer_derived Logical. If `TRUE`, return provider-normalized discovery
#' tables while deferring context linking, summaries, matrices, identity review,
#' and validation. This staging option is used by the resumable batch runner so
#' expensive derived products are built once after chunks are combined.
#' @param strict Logical. If `TRUE`, validation treats optional missing fields
#' more strictly.
#' @param curated_data Optional local species-compound table to standardize and
#' include.
#' @param provider_results Optional named list of mocked or pre-fetched provider
#' results. This is intended for tests, offline workflows, and future cached
#' provider download parsers.
#' @param enrichment_fun Optional function used instead of `categorate()` for
#' compound enrichment. It should accept a character vector of compounds and
#' return a categorate-like list.
#' @param request_fun Optional URL request function used by live-capable
#' provider adapters in tests or advanced use.
#' @param pubtator_request_fun Optional URL request function used by the
#' PubTator adapter in tests or advanced use.
#' @param refresh Logical. Reserved for future provider cache refresh support.
#' @param ... Additional arguments passed to `categorate()`, `enrichment_fun`,
#' or the PubChem-only enrichment fallback. Use `compound_request_fun` in `...`
#' to mock or customize compound-centered PubChem requests.
#'
#' @return A list with class `"uaf_plant_phytochemistry"` containing normalized
#' query, occurrence, literature, enrichment, summary, matrix, validation,
#' dictionary, and provenance tables.
#'
#' @examples
#' \dontrun{
#' plants = c("Salix nigra", "Camellia sinensis", "Zea mays")
#' phyto = resolvePlantPhytochemistry(
#'   plants = plants,
#'   sources = c("lotus", "pubmed", "pubtator"),
#'   taxon_fallback = c("species", "genus"),
#'   enrich_compounds = FALSE
#' )
#' phyto$SpeciesChemistrySummary
#' validatePlantPhytochemistryResult(phyto)$Summary
#' }
#'
#' @export
resolvePlantPhytochemistry = function(plants,
                                      sources = c("lotus", "knapsack",
                                                  "npass", "pubchem",
                                                  "pubmed", "pubtator"),
                                      taxon_fallback = c("species", "genus"),
                                      enrich_compounds = TRUE,
                                      chemical_library = NULL,
                                      detail = c("research", "full", "none"),
                                      cache = TRUE,
                                      cache_dir = NULL,
                                      lotus_index = Sys.getenv("UAFR_LOTUS_INDEX", ""),
                                      provider_indexes = NULL,
                                      require_all_providers = FALSE,
                                      throttle = 0.2,
                                      ncbi_email = Sys.getenv("NCBI_EMAIL", ""),
                                      ncbi_tool = Sys.getenv("NCBI_TOOL", "uafR"),
                                      ncbi_api_key = Sys.getenv("NCBI_API_KEY", ""),
                                      max_pubmed_records = 50,
                                      max_provider_records = max_pubmed_records,
                                      request_timeout = 30,
                                      enrich_context = FALSE,
                                      context_sources = NULL,
                                      max_context_sources = 100,
                                      min_confidence = "medium",
                                      enrichment_batch_size = Inf,
                                      resume_enrichment = TRUE,
                                      progress = interactive(),
                                      defer_derived = FALSE,
                                      strict = FALSE,
                                      curated_data = NULL,
                                      provider_results = NULL,
                                      enrichment_fun = NULL,
                                      request_fun = NULL,
                                      pubtator_request_fun = NULL,
                                      refresh = FALSE,
                                      ...) {
  detail = match.arg(detail)
  if (isTRUE(defer_derived) && isTRUE(enrich_compounds) && detail != "none") {
    stop("`defer_derived = TRUE` requires compound enrichment to be disabled.",
         call. = FALSE)
  }
  plant_queries = .plant_queries(plants, taxon_fallback)
  plant_resolution = .plant_name_resolution(plant_queries)
  sources = .plant_normalize_sources(sources)
  cache_dir = .plant_cache_dir(cache_dir)
  provider_indexes = .plant_provider_indexes(provider_indexes, lotus_index)
  lotus_index = provider_indexes$lotus
  aliases = .plant_query_aliases(plants, plant_queries)
  resources = plantProviderAvailability(
    sources = sources,
    provider_indexes = provider_indexes,
    lotus_index = lotus_index
  )
  if (isTRUE(require_all_providers)) {
    unavailable = resources$provider[
      resources$availability_status != "available"
    ]
    if (length(unavailable) > 0) {
      stop("Required plant provider resource(s) unavailable: ",
           paste(unavailable, collapse = ", "), call. = FALSE)
    }
  }

  dispatch = .plant_provider_dispatch(
    plant_queries = plant_queries,
    plant_aliases = aliases,
    sources = sources,
    cache = cache,
    cache_dir = cache_dir,
    throttle = throttle,
    ncbi_email = ncbi_email,
    ncbi_tool = ncbi_tool,
    ncbi_api_key = ncbi_api_key,
    max_pubmed_records = max_pubmed_records,
    max_provider_records = max_provider_records,
    lotus_index = lotus_index,
    provider_indexes = provider_indexes,
    request_timeout = request_timeout,
    provider_results = provider_results,
    request_fun = request_fun,
    pubtator_request_fun = pubtator_request_fun,
    refresh = refresh,
    progress = progress
  )
  if (isTRUE(require_all_providers)) {
    incomplete = dispatch$ProviderDiagnostics$provider[
      tolower(dispatch$ProviderDiagnostics$status) %in%
        c("warning", "error", "failed", "timeout", "timed_out",
          "rate_limited", "service_unavailable", "not_queried",
          "not_implemented", "unavailable")
    ]
    if (length(incomplete) > 0) {
      stop("Required plant provider stage(s) incomplete: ",
           paste(unique(incomplete), collapse = ", "), call. = FALSE)
    }
  }

  occurrence_parts = list(dispatch$PlantCompoundOccurrences)
  if (!is.null(curated_data)) {
    occurrence_parts[[length(occurrence_parts) + 1]] =
      standardizePlantCompoundIntake(curated_data)
  }
  occurrences = .plant_bind_occurrences(occurrence_parts)
  occurrences = .plant_match_occurrences_to_queries(occurrences, plant_queries)
  if (isTRUE(defer_derived)) {
    context_evidence = .uaf_empty_table(.plant_context_evidence_cols())
    provider_context_audit =
      .uaf_empty_table(.plant_provider_context_audit_cols())
  } else {
    context_evidence = plantContextEvidence(occurrences)
    context_source_rows = .plant_bind_tables(
      list(context_sources, dispatch$LiteratureCandidates),
      .plant_literature_cols()
    )
    if (isTRUE(enrich_context) || nrow(context_source_rows) > 0) {
      source_context = enrichPlantContextEvidence(
        occurrences,
        context_sources = context_source_rows,
        fetch_pubmed = isTRUE(enrich_context),
        cache = cache,
        cache_dir = file.path(cache_dir, "context_enrichment"),
        throttle = throttle,
        ncbi_email = ncbi_email,
        ncbi_tool = ncbi_tool,
        ncbi_api_key = ncbi_api_key,
        max_sources = max_context_sources,
        request_fun = request_fun,
        request_timeout = request_timeout,
        min_confidence = "low",
        apply = FALSE
      )
      context_evidence = .plant_bind_tables(
        list(context_evidence, source_context),
        .plant_context_evidence_cols()
      )
    }
    if (nrow(context_evidence) > 0) {
      occurrences = .plant_apply_context_evidence(occurrences,
                                                  context_evidence)
    }
    occurrences = .plant_clean_context_conflicts(occurrences)
    provider_context_audit = plantProviderContextAudit(
      list(PlantCompoundOccurrences = occurrences,
           PlantContextEvidence = context_evidence)
    )
  }

  enrichment = if (isTRUE(enrich_compounds) && detail != "none") {
    enrichPlantCompounds(occurrences,
                         chemical_library = chemical_library,
                         detail = detail,
                         cache = cache,
                         cache_dir = file.path(cache_dir, "compound_enrichment"),
                         throttle = throttle,
                         enrichment_fun = enrichment_fun,
                         batch_size = enrichment_batch_size,
                         resume = resume_enrichment,
                         progress = progress,
                         ...)
  } else {
    list(CategorateResult = NULL,
         CompoundResolution = .plant_compound_resolution(occurrences, NULL),
         TraitEvidence = .uaf_empty_table(.plant_trait_evidence_cols()),
         Provenance = .plant_provenance("compound_enrichment", "not_run",
                                        NA_character_, NA_character_, 0,
                                        "Compound enrichment was not requested."))
  }

  if (isTRUE(defer_derived)) {
    comparability = .uaf_empty_table(.plant_comparability_cols())
    summary = .uaf_empty_table(.plant_summary_cols())
    matrix = .uaf_empty_table(c("species"))
    comparable_matrix = .uaf_empty_table(c("species"))
    compound_identity_review = .plant_empty_compound_identity_review()
  } else {
    comparability = plantChemistryComparability(
      list(PlantCompoundOccurrences = occurrences,
           CategorateResult = enrichment$CategorateResult),
      min_confidence = "low"
    )
    summary = summarizePlantPhytochemistry(
      plant_compounds = occurrences,
      categorate_result = enrichment$CategorateResult,
      compound_resolution = enrichment$CompoundResolution,
      plant_queries = plant_queries,
      provider_diagnostics = dispatch$ProviderDiagnostics,
      comparability = comparability
    )
    matrix = plantPhytochemistryMatrix(
      list(PlantCompoundOccurrences = occurrences,
           CategorateResult = enrichment$CategorateResult),
      level = "species",
      profile = "core",
      mode = "binary",
      min_confidence = min_confidence,
      max_traits = .plant_automatic_matrix_max_traits()
    )
    comparable_matrix = plantComparableChemistryMatrix(
      list(ChemistryComparability = comparability),
      comparison_scope = "specialized_metabolites",
      level = "species",
      mode = "binary",
      min_comparability_confidence = min_confidence
    )
    compound_identity_review = plantCompoundIdentityReviewTable(
      list(PlantCompoundOccurrences = occurrences,
           CompoundResolution = enrichment$CompoundResolution)
    )
  }

  provenance = .plant_bind_tables(list(
    dispatch$Provenance,
    enrichment$Provenance,
    .plant_provenance("resolvePlantPhytochemistry", "uafR",
                      paste(plant_queries$query_plant, collapse = "; "),
                      NA_character_, nrow(occurrences),
                      "Species-first plant phytochemistry result assembled.")
  ), .plant_provenance_cols())

  out = list(
    PlantQueries = plant_queries,
    PlantQueryAliases = aliases,
    PlantNameResolution = plant_resolution,
    ProviderDiagnostics = dispatch$ProviderDiagnostics,
    ProviderQueryAccounting = dispatch$ProviderQueryAccounting,
    ProviderResourceManifest = resources,
    PlantCompoundOccurrences = occurrences,
    PlantContextEvidence = context_evidence,
    ProviderContextAudit = provider_context_audit,
    LiteratureCandidates = dispatch$LiteratureCandidates,
    SourceCompoundIdentity = dispatch$SourceCompoundIdentity,
    CompoundResolution = enrichment$CompoundResolution,
    CompoundIdentityReview = compound_identity_review,
    CategorateResult = enrichment$CategorateResult,
    SpeciesChemistrySummary = summary,
    SpeciesChemistryMatrix = matrix,
    ChemistryComparability = comparability,
    ComparableChemistryMatrix = comparable_matrix,
    TraitEvidence = enrichment$TraitEvidence,
    Validation = list(Summary = data.frame(),
                      TableQuality = data.frame(),
                      ProviderDiagnostics = data.frame(),
                      Issues = data.frame(),
                      DataDictionary = data.frame()),
    DataDictionary = plantPhytochemistrySchema(),
    Provenance = provenance
  )
  class(out) = c("uaf_plant_phytochemistry", class(out))
  if (!isTRUE(defer_derived)) {
    out$Validation = validatePlantPhytochemistryResult(out, strict = strict)
  }
  out
}

#' Describe plant phytochemistry output schemas
#'
#' @param tables Optional table names to return. If `NULL`, all schemas are
#' returned.
#'
#' @return A data frame describing normalized table, column, type, required,
#' role, allowed-value, and description contracts.
#'
#' @export
plantPhytochemistrySchema = function(tables = NULL) {
  dictionary = .plant_data_dictionary()
  tables = .uaf_non_empty(tables)
  if (length(tables) > 0) {
    dictionary = dictionary[dictionary$Table %in% tables, , drop = FALSE]
  }
  row.names(dictionary) = NULL
  dictionary
}

#' Validate plant phytochemistry results
#'
#' @param x A result returned by `resolvePlantPhytochemistry()` or a named list
#' containing plant phytochemistry tables.
#' @param strict Logical. If `TRUE`, missing optional documented columns are
#' reported as warnings.
#'
#' @return A list with `Summary`, `TableQuality`, `ProviderDiagnostics`,
#' `Issues`, and `DataDictionary`.
#'
#' @export
validatePlantPhytochemistryResult = function(x, strict = FALSE) {
  dictionary = plantPhytochemistrySchema()
  expected_tables = unique(dictionary$Table)
  issue_rows = list()
  quality_rows = list()

  if (!is.list(x) || is.data.frame(x)) {
    issues = .plant_validation_issue("error", NA_character_, NA_character_,
                                     "Input is not a plant phytochemistry result list",
                                     "list", class(x)[[1]], NA_integer_,
                                     NA_character_)
    return(.plant_validation_result(.uaf_empty_table(.plant_table_quality_cols()),
                                    .uaf_empty_table(.plant_provider_diagnostic_cols()),
                                    issues, dictionary))
  }

  for (table_name in expected_tables) {
    table_dictionary = dictionary[dictionary$Table == table_name, ,
                                  drop = FALSE]
    table = x[[table_name]]
    quality_rows[[length(quality_rows) + 1]] =
      .plant_table_quality(table_name, table, table_dictionary, strict)
    issues = .plant_validate_table(table_name, table, table_dictionary, strict)
    if (nrow(issues) > 0) issue_rows[[length(issue_rows) + 1]] = issues
  }

  issues = .plant_bind_tables(issue_rows, .plant_validation_issue_cols())
  issues = .plant_extra_validation_issues(x, issues)
  table_quality = .plant_bind_tables(quality_rows, .plant_table_quality_cols())
  provider_diagnostics = if (is.data.frame(x$ProviderDiagnostics)) {
    x$ProviderDiagnostics
  } else {
    .uaf_empty_table(.plant_provider_diagnostic_cols())
  }

  .plant_validation_result(table_quality, provider_diagnostics, issues,
                           dictionary)
}

#' Standardize curated plant-compound intake
#'
#' @param x Data frame containing at least `species` and `compound_name`.
#' Recommended columns include `source_database`, `citation_or_url`, and
#' `evidence_tier`.
#'
#' @return A normalized `PlantCompoundOccurrences` data frame.
#'
#' @export
standardizePlantCompoundIntake = function(x) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame with plant-compound rows.", call. = FALSE)
  }
  names(x) = .plant_normalize_column_names(names(x))
  if (!"species" %in% names(x) || !"compound_name" %in% names(x)) {
    stop("Curated plant-compound intake requires `species` and ",
         "`compound_name` columns.", call. = FALSE)
  }

  out = .plant_empty_occurrences()
  if (nrow(x) < 1) return(out)
  species = .uaf_squish_text(x$species)
  genus = if ("genus" %in% names(x)) .uaf_squish_text(x$genus) else
    .plant_genus(species)
  family = if ("family" %in% names(x)) .uaf_squish_text(x$family) else
    rep(NA_character_, nrow(x))
  source_database = .plant_col_or_default(x, "source_database", "manual")
  citation = .plant_col_or_default(x, "citation_or_url", NA_character_)
  evidence_tier = .plant_col_or_default(x, "evidence_tier", "manual_curated")
  confidence = .plant_evidence_tier_confidence(evidence_tier)

  out = data.frame(
    query_plant = species,
    query_plant_clean = .plant_clean_name(species),
    matched_taxon = species,
    matched_rank = "species",
    species = species,
    genus = genus,
    family = family,
    compound_name = .uaf_squish_text(x$compound_name),
    compound_name_clean = .plant_clean_compound(x$compound_name),
    compound_id = .plant_col_or_default(x, "compound_id", NA_character_),
    compound_id_type = .plant_col_or_default(x, "compound_id_type", NA_character_),
    source_database = source_database,
    source_record_id = .plant_col_or_default(x, "source_record_id", NA_character_),
    evidence_text = .plant_col_or_default(x, "evidence_note", NA_character_),
    evidence_url = citation,
    reference_id = .plant_col_or_default(x, "reference_id", NA_character_),
    pmid = .plant_col_or_default(x, "pmid", NA_character_),
    doi = .plant_col_or_default(x, "doi", NA_character_),
    plant_part = .plant_col_or_default(x, "plant_part", NA_character_),
    tissue = .plant_col_or_default(x, "tissue", NA_character_),
    method = .plant_col_or_default(x, "method", NA_character_),
    occurrence_type = .plant_col_or_default(x, "occurrence_type", "reported"),
    retrieved_at = .plant_timestamp(),
    confidence = confidence,
    curation_flag = .plant_col_or_default(x, "curation_status", "curated"),
    evidence_tier = evidence_tier,
    stringsAsFactors = FALSE
  )
  out = .plant_normalize_occurrences(out, source_hint = "manual")
  row.names(out) = NULL
  out
}

#' Standardize a local LOTUS species-compound index
#'
#' @description
#' Converts a flat LOTUS export into a compact species-compound index that uafR
#' can query locally. This is the recommended LOTUS route for medium or large
#' plant panels because the live LOTUS simple API is unpaged and can return very
#' large payloads for common species.
#'
#' The input may be a data frame or a path to a CSV, TSV, JSON, JSONL, NDJSON,
#' or RDS file produced from a LOTUS MongoDB, Wikidata, SDF metadata, or other
#' flat export. Column names are matched flexibly; useful fields include
#' `species`, `genus`, `family`, `allTaxa`, `traditional_name`, `iupac_name`,
#' `lotus_id`, `inchikey`, `smiles`, `molecular_formula`, `doi`, `pmid`,
#' `plant_part`, `tissue`, and `method`.
#'
#' @param x Data frame or path to a flat LOTUS export.
#' @param source_file Optional source label recorded in the standardized index.
#'
#' @return A normalized LOTUS index data frame with one row per
#' taxon-compound record where taxon evidence can be extracted.
#'
#' @export
standardizeLotusIndex = function(x, source_file = NULL) {
  if (is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)) {
    source_file = .uaf_first_non_empty_text(source_file, x)
    x = .plant_read_lotus_index_source(x)
  }
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame or a path to a flat LOTUS export.",
         call. = FALSE)
  }
  if (nrow(x) < 1) return(.plant_empty_lotus_index())
  names(x) = .plant_normalize_column_names(names(x))
  if (all(.plant_lotus_index_cols() %in% names(x))) {
    return(.plant_finalize_lotus_index(x[, .plant_lotus_index_cols(),
                                         drop = FALSE]))
  }
  rows = lapply(seq_len(nrow(x)), function(i) {
    .plant_standardize_lotus_index_row(x[i, , drop = FALSE], source_file)
  })
  out = .plant_bind_tables(rows, .plant_lotus_index_cols())
  if (nrow(out) < 1) return(.plant_empty_lotus_index())
  .plant_finalize_lotus_index(out)
}

.plant_finalize_lotus_index = function(out) {
  out = .plant_bind_tables(list(out), .plant_lotus_index_cols())
  if (nrow(out) < 1) return(.plant_empty_lotus_index())
  out$species = .uaf_squish_text(out$species)
  out$species_clean = .plant_clean_name(out$species)
  out$genus = .uaf_squish_text(ifelse(is.na(out$genus) | out$genus == "",
                                      .plant_genus(out$species), out$genus))
  out$family = .uaf_squish_text(out$family)
  out$compound_name = .uaf_squish_text(out$compound_name)
  out$compound_name_clean = .plant_clean_compound(out$compound_name)
  out = out[!is.na(out$compound_name) & out$compound_name != "" &
              (!is.na(out$species) | !is.na(out$genus) |
                 !is.na(out$family)), , drop = FALSE]
  out = unique(out)
  row.names(out) = NULL
  out
}

#' Query a local LOTUS index for plant phytochemistry records
#'
#' @description
#' Searches a standardized local LOTUS index for direct species records and,
#' when requested, genus or family fallback records. The returned table uses the
#' same `PlantCompoundOccurrences` schema as `resolvePlantPhytochemistry()`.
#'
#' @param plants Character vector or data frame of plant names.
#' @param lotus_index Data frame or path accepted by `standardizeLotusIndex()`,
#' or a manifest-backed lookup directory produced by
#' `tools/flatten_lotus_mongo_dump.py --lookup-dir`.
#' @param taxon_fallback Fallback ranks. Species records are always queried;
#' optional `"genus"` and `"family"` records are labeled as fallback evidence.
#' @param max_records Maximum records retained per input plant.
#'
#' @return A normalized `PlantCompoundOccurrences` data frame.
#'
#' @export
queryLotusIndex = function(plants, lotus_index,
                           taxon_fallback = c("species", "genus"),
                           max_records = Inf) {
  plant_queries = .plant_queries(plants, taxon_fallback)
  aliases = .plant_query_aliases(plants, plant_queries)
  .plant_query_lotus_index(
    plant_queries, lotus_index, max_records, plant_aliases = aliases
  )
}

#' Build a compact local LOTUS index
#'
#' @description
#' Builds the compact LOTUS index used by `queryLotusIndex()` and the
#' `lotus_index` argument in plant phytochemistry workflows. Input can be one
#' or more flat LOTUS exports, a directory containing supported flat exports, or
#' an in-memory data frame. Supported file formats are CSV, TSV, JSON, JSONL,
#' NDJSON, and RDS.
#'
#' Raw LOTUS MongoDB downloads are commonly distributed as zipped BSON dumps.
#' Use `tools/flatten_lotus_mongo_dump.py` to stream the official LOTUS MongoDB
#' ZIP into an auditable flat CSV, a compact CSV, and optionally a
#' manifest-backed lookup directory. Then pass the flat/compact file to this
#' function or pass the lookup directory directly to `queryLotusIndex()` or
#' `resolvePlantPhytochemistry(lotus_index = ...)`.
#'
#' @param input Data frame, file path, directory path, or character vector of
#' file/directory paths containing flat LOTUS export rows.
#' @param out_file Optional output path for the compact LOTUS index.
#' @param format Output format. `"auto"` infers from `out_file`; otherwise use
#' `"csv"` or `"rds"`.
#' @param overwrite Logical. If `FALSE`, existing output files are not
#' replaced.
#' @param manifest_file Optional JSON manifest path. If omitted and `out_file`
#' is supplied, a sidecar `*_manifest.json` file is written.
#' @param recursive Logical. If `TRUE`, directory inputs are searched
#' recursively for supported flat export files.
#' @param strict Logical. If `TRUE`, stop on the first unreadable input source.
#' If `FALSE`, unreadable sources are recorded in the build manifest and other
#' sources are still processed.
#' @param progress Logical. If `TRUE`, print source-processing progress.
#'
#' @return A list with `LotusIndex`, `BuildSummary`, and `BuildManifest`.
#'
#' @export
buildLotusIndex = function(input,
                           out_file = NULL,
                           format = c("auto", "csv", "rds"),
                           overwrite = FALSE,
                           manifest_file = NULL,
                           recursive = TRUE,
                           strict = FALSE,
                           progress = interactive()) {
  format = match.arg(format)
  sources = .plant_lotus_build_sources(input, recursive = recursive)
  if (length(sources) < 1) {
    stop("No LOTUS input sources were found.", call. = FALSE)
  }

  index_parts = list()
  manifest_rows = list()
  for (i in seq_along(sources)) {
    source = sources[[i]]
    .plant_progress(progress, "uafR LOTUS index source ", i, "/",
                    length(sources), ": ", source$label)
    started = Sys.time()
    status = "ok"
    error_message = NA_character_
    raw_rows = NA_integer_
    index_rows = 0L
    index = .plant_empty_lotus_index()
    raw = tryCatch(
      .plant_lotus_build_read_source(source),
      error = function(error) {
        if (isTRUE(strict)) stop(error)
        status <<- "error"
        error_message <<- .plant_redact_secrets(conditionMessage(error))
        NULL
      }
    )
    if (is.data.frame(raw)) {
      raw_rows = nrow(raw)
      index = tryCatch(
        standardizeLotusIndex(raw, source_file = source$label),
        error = function(error) {
          if (isTRUE(strict)) stop(error)
          status <<- "error"
          error_message <<- .plant_redact_secrets(conditionMessage(error))
          .plant_empty_lotus_index()
        }
      )
      index_rows = nrow(index)
      if (index_rows < 1 && status == "ok") status = "no_taxon_records"
      if (index_rows > 0) index_parts[[length(index_parts) + 1]] = index
    }
    manifest_rows[[length(manifest_rows) + 1]] = data.frame(
      source_id = i,
      source_label = source$label,
      source_type = source$type,
      source_path = source$path,
      raw_row_count = raw_rows,
      index_row_count = index_rows,
      status = status,
      error_message = error_message,
      elapsed_seconds = round(as.numeric(difftime(Sys.time(), started,
                                                  units = "secs")), 3),
      processed_at = .plant_timestamp(),
      stringsAsFactors = FALSE
    )
  }

  index = .plant_bind_tables(index_parts, .plant_lotus_index_cols())
  pre_dedup_rows = nrow(index)
  index = .plant_deduplicate_lotus_index(index)
  if (nrow(index) < 1) {
    stop("No usable LOTUS taxon-compound rows were extracted. Check that the ",
         "input includes compound names and taxon fields such as `allTaxa`, ",
         "`species`, `organism`, or `taxonomy`.", call. = FALSE)
  }

  output_file = .uaf_first_non_empty_text(out_file)
  resolved_format = .plant_lotus_index_output_format(output_file, format)
  written_file = NA_character_
  if (!is.na(output_file)) {
    written_file = .plant_write_lotus_index_file(index, output_file,
                                                 resolved_format, overwrite)
  }
  manifest_file = .plant_lotus_manifest_file(manifest_file, output_file)
  summary = .plant_lotus_index_build_summary(index, manifest_rows,
                                             pre_dedup_rows, written_file,
                                             manifest_file)
  manifest = .plant_bind_tables(manifest_rows,
                                .plant_lotus_build_manifest_cols())
  if (!is.na(manifest_file)) {
    if (file.exists(manifest_file) && !isTRUE(overwrite)) {
      stop("Manifest file exists and `overwrite = FALSE`: ", manifest_file,
           call. = FALSE)
    }
    dir.create(dirname(manifest_file), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(
      list(BuildSummary = summary, BuildManifest = manifest),
      path = manifest_file,
      dataframe = "rows",
      pretty = TRUE,
      na = "null"
    )
  }

  out = list(LotusIndex = index,
             BuildSummary = summary,
             BuildManifest = manifest)
  class(out) = c("uaf_lotus_index_build", class(out))
  out
}

#' Enrich plant compounds with existing uafR compound workflows
#'
#' @param plant_compounds A normalized or curated plant-compound table.
#' @param chemical_library Optional library passed to `categorate()`.
#' @param detail One of `"research"`, `"full"`, or `"none"`.
#' @param cache Logical.
#' @param cache_dir Cache directory.
#' @param throttle Request throttle.
#' @param enrichment_fun Optional injected enrichment function for tests or
#' cached workflows.
#' @param pubchem_fun Optional replacement for `pubchemProfile()` used by the
#' PubChem-only fallback. Intended for tests and advanced cached workflows.
#' @param compound_request_fun Optional request function passed to
#' `pubchemProfile()` when the PubChem-only fallback is used.
#' @param batch_size Maximum number of unique compounds per PubChem-only
#' enrichment batch. Use `Inf` for one batch.
#' @param resume Logical. If `TRUE`, reuse saved PubChem-only batch `.rds`
#' files in `cache_dir`.
#' @param progress Logical. If `TRUE`, print simple progress messages during
#' PubChem-only enrichment.
#' @param ... Additional arguments passed to `categorate()`, `enrichment_fun`,
#' or the PubChem-only enrichment fallback.
#'
#' @return A list with `CategorateResult`, `CompoundResolution`,
#' `TraitEvidence`, and `Provenance`.
#'
#' @export
enrichPlantCompounds = function(plant_compounds,
                                chemical_library = NULL,
                                detail = c("research", "full", "none"),
                                cache = TRUE,
                                cache_dir = NULL,
                                throttle = 0.2,
                                enrichment_fun = NULL,
                                pubchem_fun = NULL,
                                compound_request_fun = NULL,
                                batch_size = Inf,
                                resume = TRUE,
                                progress = interactive(),
                                ...) {
  detail = match.arg(detail)
  occurrences = if (is.data.frame(plant_compounds) &&
                    all(.plant_occurrence_cols() %in% names(plant_compounds))) {
    .plant_normalize_occurrences(plant_compounds)
  } else {
    standardizePlantCompoundIntake(plant_compounds)
  }
  compounds = unique(.uaf_non_empty(occurrences$compound_name))
  if (length(compounds) < 1 || detail == "none") {
    return(list(
      CategorateResult = NULL,
      CompoundResolution = .plant_compound_resolution(occurrences, NULL),
      TraitEvidence = .uaf_empty_table(.plant_trait_evidence_cols()),
      Provenance = .plant_provenance("compound_enrichment", "not_run",
                                     paste(compounds, collapse = "; "),
                                     NA_character_, 0,
                                     "No compounds were available for enrichment or detail = 'none'.")
    ))
  }

  categorate_result = tryCatch({
    if (!is.null(enrichment_fun)) {
      enrichment_fun(compounds, ...)
    } else if (!is.null(chemical_library)) {
      categorate(compounds = compounds,
                 chemical_library = chemical_library,
                 detail = detail,
                 cache = cache,
                 cache_dir = cache_dir,
                 throttle = throttle,
                 ...)
    } else {
      .plant_pubchem_only_enrichment(
        compounds = compounds,
        detail = detail,
        cache = cache,
        cache_dir = cache_dir,
        throttle = throttle,
        pubchem_fun = pubchem_fun,
        request_fun = compound_request_fun,
        batch_size = batch_size,
        resume = resume,
        progress = progress,
        ...
      )
    }
  }, error = function(error) {
    warning("Plant compound enrichment failed: ", conditionMessage(error),
            call. = FALSE)
    NULL
  })

  list(
    CategorateResult = categorate_result,
    CompoundResolution = .plant_compound_resolution(occurrences,
                                                    categorate_result),
    TraitEvidence = .plant_trait_evidence(occurrences, categorate_result),
    Provenance = .plant_provenance("compound_enrichment",
                                   .plant_enrichment_source(categorate_result),
                                   paste(compounds, collapse = "; "),
                                   NA_character_, length(compounds),
                                   ifelse(is.null(categorate_result),
                                          "Compound enrichment not available.",
                                          .plant_enrichment_note(categorate_result)))
  )
}

#' Resolve plant compound identities with a fast PubChem pass
#'
#' @description
#' Resolves unique compound names from a plant occurrence table or plant
#' phytochemistry result to PubChem identity fields without pulling the richer
#' annotation tables used by `categorate(detail = "research")`. This is the
#' preferred first compound-resolution step for large species-first runs.
#'
#' @param plant_compounds Plant phytochemistry result, normalized occurrence
#' table, or curated plant-compound intake table.
#' @param cache Logical. If `TRUE`, PubChem responses and identity batch results
#' can be reused.
#' @param cache_dir Cache directory.
#' @param throttle Seconds to wait between uncached PubChem requests.
#' @param batch_size Maximum number of unique compounds per identity batch.
#' @param resume Logical. If `TRUE`, reuse saved identity batch `.rds` files.
#' @param progress Logical. If `TRUE`, print simple progress messages.
#' @param pubchem_fun Optional replacement for `pubchemProfile()` used in tests
#' or advanced cached workflows.
#' @param compound_request_fun Optional request function passed to
#' `pubchemProfile()`.
#' @param lotus_index Optional local LOTUS index as a data frame, flat file, or
#' manifest-backed lookup directory. When supplied, source-backed LOTUS SMILES,
#' InChIKeys, formulas, and PubChem CIDs are used before PubChem name lookup.
#' @param source_only Logical. If `TRUE`, resolve only exact identities available
#' from `lotus_index` and do not make PubChem requests. Unresolved and ambiguous
#' names remain explicit for later review or server-side resolution.
#'
#' @return A `CompoundResolution` data frame. The underlying identity-only
#' PubChem tables are attached as the `"CategorateResult"` attribute. Compound
#' keys preserve chemically meaningful Greek-letter and plus/minus prefixes so
#' isomers such as alpha-pinene and beta-pinene are not collapsed together.
#'
#' @export
resolvePlantCompoundIdentities = function(plant_compounds,
                                          cache = TRUE,
                                          cache_dir = NULL,
                                          throttle = 0.2,
                                          batch_size = 100,
                                          resume = TRUE,
                                          progress = interactive(),
                                          pubchem_fun = NULL,
                                          compound_request_fun = NULL,
                                          lotus_index = NULL,
                                          source_only = FALSE) {
  occurrences = .plant_occurrences_from_input(plant_compounds)
  source_identity = .plant_source_identity_from_input(plant_compounds)
  identity = .plant_resolve_compound_identities(
    occurrences = occurrences,
    cache = cache,
    cache_dir = cache_dir,
    throttle = throttle,
    batch_size = batch_size,
    resume = resume,
    progress = progress,
    pubchem_fun = pubchem_fun,
    request_fun = compound_request_fun,
    lotus_index = lotus_index,
    source_identity = source_identity,
    source_only = source_only
  )
  out = identity$CompoundResolution
  attr(out, "CategorateResult") = identity$CategorateResult
  attr(out, "SourceCompoundIdentity") =
    identity$CategorateResult$SourceCompoundIdentity
  attr(out, "Provenance") = identity$Provenance
  out
}

#' Create a plant compound identity review table
#'
#' @description
#' Builds an auditable worksheet for unresolved or ambiguous compound identity
#' rows from species-first plant phytochemistry workflows. The table classifies
#' likely failure modes, ranks high-impact rows by occurrence count and species
#' coverage, and recommends conservative next actions. It does not change
#' compound identities and does not infer spelling corrections.
#'
#' @param x Plant phytochemistry result, `CompoundResolution` data frame, or a
#' list containing `CompoundResolution` and optionally
#' `PlantCompoundOccurrences`.
#' @param occurrences Optional normalized occurrence table used to add source,
#' species, and record-count context when `x` is only a resolution table.
#' @param include_resolved Logical. If `FALSE`, already resolved rows are
#' omitted.
#' @param min_query_count Minimum `query_count` to include.
#'
#' @return A data frame suitable for CSV export, review, and re-import into a
#' project-specific curation workflow.
#'
#' @export
plantCompoundIdentityReviewTable = function(x,
                                            occurrences = NULL,
                                            include_resolved = FALSE,
                                            min_query_count = 1) {
  resolution = if (inherits(x, "uaf_plant_phytochemistry") &&
                   is.data.frame(x$CompoundResolution)) {
    x$CompoundResolution
  } else if (is.list(x) && is.data.frame(x$CompoundResolution)) {
    x$CompoundResolution
  } else if (is.data.frame(x)) {
    x
  } else {
    .uaf_empty_table(.plant_compound_resolution_cols())
  }
  if (is.null(occurrences)) {
    occurrences = if (inherits(x, "uaf_plant_phytochemistry") &&
                      is.data.frame(x$PlantCompoundOccurrences)) {
      x$PlantCompoundOccurrences
    } else if (is.list(x) && is.data.frame(x$PlantCompoundOccurrences)) {
      x$PlantCompoundOccurrences
    } else {
      .plant_empty_occurrences()
    }
  }
  resolution = .plant_bind_tables(list(resolution),
                                  .plant_compound_resolution_cols())
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(resolution) < 1) return(.plant_empty_compound_identity_review())
  issues = lapply(seq_len(nrow(resolution)), function(i) {
    .plant_identity_issue(resolution[i, , drop = FALSE])
  })
  review_required_issue = vapply(issues, function(issue) {
    issue$type %in% .plant_identity_review_required_types()
  }, logical(1))
  keep = rep(TRUE, nrow(resolution))
  if (!isTRUE(include_resolved)) {
    keep = keep & (!(resolution$resolved %in% TRUE) | review_required_issue)
  }
  counts = suppressWarnings(as.integer(resolution$query_count))
  counts[!is.finite(counts)] = 0L
  keep = keep & counts >= min_query_count
  resolution = resolution[keep, , drop = FALSE]
  issues = issues[keep]
  if (nrow(resolution) < 1) return(.plant_empty_compound_identity_review())
  occurrence_index = split(seq_len(nrow(occurrences)),
                           occurrences$compound_name_clean)
  occurrence_info = lapply(resolution$compound_name_clean, function(key) {
    idx = occurrence_index[[key]]
    if (is.null(idx)) idx = integer()
    source_records = unique(.uaf_non_empty(occurrences$source_record_id[idx]))
    source_databases = unique(.uaf_non_empty(occurrences$source_database[idx]))
    species = unique(.uaf_non_empty(occurrences$species[idx]))
    list(
      occurrence_count = length(idx),
      source_record_count = length(source_records),
      source_databases = .pubchem_collapse(source_databases),
      source_record_ids = .pubchem_collapse(utils::head(source_records, 25)),
      species_count = length(species),
      example_species = .pubchem_collapse(utils::head(species, 10))
    )
  })
  issue_type = vapply(issues, `[[`, character(1), "type")
  occurrence_count = vapply(occurrence_info, `[[`, integer(1),
                            "occurrence_count")
  species_count = vapply(occurrence_info, `[[`, integer(1), "species_count")
  priority = .plant_identity_review_priority_vector(
    resolved = resolution$resolved,
    query_count = resolution$query_count,
    occurrence_count = occurrence_count,
    species_count = species_count,
    issue_type = issue_type
  )
  out = data.frame(
    review_id = sprintf("compound_identity_%04d", seq_len(nrow(resolution))),
    review_decision = "needs_review",
    reviewed_by = NA_character_,
    reviewed_at = NA_character_,
    review_note = NA_character_,
    compound_name = resolution$compound_name,
    compound_name_clean = resolution$compound_name_clean,
    query_count = resolution$query_count,
    resolved = resolution$resolved,
    CID = resolution$CID,
    InChIKey = resolution$InChIKey,
    SMILES = resolution$SMILES,
    MolecularFormula = resolution$MolecularFormula,
    resolution_source = resolution$resolution_source,
    identity_issue_type = issue_type,
    review_priority = priority,
    recommended_decision = vapply(issues, `[[`, character(1), "decision"),
    proposed_compound_name = NA_character_,
    proposed_cid = NA_character_,
    proposed_inchikey = NA_character_,
    proposed_smiles = NA_character_,
    proposed_molecular_formula = NA_character_,
    proposed_resolution_source = NA_character_,
    suggested_query = vapply(issues, `[[`, character(1), "suggested_query"),
    suggested_source = vapply(issues, `[[`, character(1), "suggested_source"),
    source_record_count = vapply(occurrence_info, `[[`, integer(1),
                                 "source_record_count"),
    source_databases = vapply(occurrence_info, `[[`, character(1),
                               "source_databases"),
    source_record_ids = vapply(occurrence_info, `[[`, character(1),
                                "source_record_ids"),
    species_count = species_count,
    example_species = vapply(occurrence_info, `[[`, character(1),
                             "example_species"),
    review_reason = vapply(issues, `[[`, character(1), "reason"),
    notes = resolution$notes,
    stringsAsFactors = FALSE
  )
  out = .plant_bind_tables(list(out), .plant_compound_identity_review_cols())
  out$.priority_rank = match(out$review_priority,
                             c("high", "medium", "low", "resolved"))
  out$.priority_rank[is.na(out$.priority_rank)] = 99L
  out = out[order(out$.priority_rank,
                  -suppressWarnings(as.integer(out$query_count)),
                  out$identity_issue_type,
                  out$compound_name_clean), , drop = FALSE]
  out$.priority_rank = NULL
  row.names(out) = NULL
  out
}

#' Apply reviewed plant compound identity decisions
#'
#' @description
#' Applies a completed `plantCompoundIdentityReviewTable()` worksheet to a
#' plant phytochemistry result or `CompoundResolution` table. This function is
#' intentionally conservative: it only updates identity fields when a reviewer
#' supplies explicit proposed CID, InChIKey, SMILES, formula, or source fields.
#' It does not infer spelling corrections, merge ambiguous source structures,
#' or promote class/product labels automatically.
#'
#' Supported `review_decision` values are:
#' `"needs_review"`/`"keep"` (leave unchanged), `"accept_resolved"` (append a
#' review note without changing identity fields), `"update_identity"` or
#' `"replace_identity"` (copy proposed identity fields and mark the row
#' resolved when at least one identity field is supplied), and `"reject"` or
#' `"exclude"` (mark the row unresolved with `resolution_source =
#' "review_excluded"`).
#'
#' @param x Plant phytochemistry result, named list with `CompoundResolution`,
#' or a `CompoundResolution` data frame.
#' @param review_table Completed table from `plantCompoundIdentityReviewTable()`.
#' @param reviewer Optional reviewer name used when `reviewed_by` is empty.
#' @param require_identity Logical. If `TRUE`, update/replace decisions require
#' at least one proposed CID, InChIKey, SMILES, molecular formula, or compound
#' name.
#'
#' @return Updated plant phytochemistry result or `CompoundResolution` table.
#'
#' @export
applyPlantCompoundIdentityReview = function(x,
                                            review_table,
                                            reviewer = NA_character_,
                                            require_identity = TRUE) {
  if (!is.data.frame(review_table)) {
    stop("`review_table` must be a data frame.", call. = FALSE)
  }
  is_result = inherits(x, "uaf_plant_phytochemistry") ||
    (is.list(x) && !is.data.frame(x) && is.data.frame(x$CompoundResolution))
  resolution = if (is_result) x$CompoundResolution else x
  resolution = .plant_bind_tables(list(resolution),
                                  .plant_compound_resolution_cols())
  if (nrow(resolution) < 1 || nrow(review_table) < 1) {
    return(if (is_result) x else resolution)
  }
  names(review_table) = .plant_normalize_column_names(names(review_table))
  for (col in .plant_compound_identity_review_cols()) {
    if (!col %in% names(review_table)) {
      review_table[[col]] = rep(NA_character_, nrow(review_table))
    }
  }
  for (i in seq_len(nrow(review_table))) {
    review = review_table[i, , drop = FALSE]
    key = .uaf_first_non_empty_text(review$compound_name_clean,
                                    .plant_clean_compound(review$compound_name))
    if (is.na(key) || key == "") next
    idx = which(resolution$compound_name_clean == key)
    if (length(idx) < 1) next
    decision = .plant_compound_identity_review_decision(
      review$review_decision
    )
    stamp = .plant_compound_identity_review_stamp(review, reviewer)
    if (decision %in% c("needs_review", "keep")) next
    if (decision %in% c("reject", "exclude")) {
      resolution$resolved[idx] = FALSE
      resolution$resolution_source[idx] = "review_excluded"
      resolution$notes[idx] = vapply(
        resolution$notes[idx],
        .plant_append_note,
        character(1),
        note = paste("Identity excluded by review.", stamp)
      )
      next
    }
    if (decision == "accept_resolved") {
      resolution$notes[idx] = vapply(
        resolution$notes[idx],
        .plant_append_note,
        character(1),
        note = paste("Resolved identity accepted by review.", stamp)
      )
      next
    }
    proposed = .plant_reviewed_identity_values(review)
    has_identity = length(.uaf_non_empty(unlist(proposed, use.names = FALSE))) >
      0
    if (isTRUE(require_identity) && !has_identity) {
      warning("Skipping compound identity review row ", i,
              ": update/replace decision has no proposed identity fields.",
              call. = FALSE)
      next
    }
    if (!is.na(proposed$compound_name)) {
      resolution$compound_name[idx] = proposed$compound_name
    }
    if (!is.na(proposed$CID)) resolution$CID[idx] = proposed$CID
    if (!is.na(proposed$InChIKey)) resolution$InChIKey[idx] = proposed$InChIKey
    if (!is.na(proposed$SMILES)) resolution$SMILES[idx] = proposed$SMILES
    if (!is.na(proposed$MolecularFormula)) {
      resolution$MolecularFormula[idx] = proposed$MolecularFormula
    }
    resolution$resolved[idx] = length(.uaf_non_empty(c(
      resolution$CID[idx], resolution$InChIKey[idx], resolution$SMILES[idx],
      resolution$MolecularFormula[idx]
    ))) > 0
    resolution$resolution_source[idx] = .uaf_first_non_empty_text(
      proposed$resolution_source, "manual_identity_review"
    )
    resolution$notes[idx] = vapply(
      resolution$notes[idx],
      .plant_append_note,
      character(1),
      note = paste("Identity updated by review.", stamp)
    )
  }
  resolution = .plant_collapse_compound_resolution_keys(resolution)
  if (!is_result) return(resolution)
  x$CompoundResolution = resolution
  x$CompoundIdentityReview = plantCompoundIdentityReviewTable(x)
  if (is.data.frame(x$PlantCompoundOccurrences)) {
    x$SpeciesChemistrySummary = summarizePlantPhytochemistry(
      plant_compounds = x$PlantCompoundOccurrences,
      categorate_result = x$CategorateResult,
      compound_resolution = resolution,
      plant_queries = x$PlantQueries,
      provider_diagnostics = x$ProviderDiagnostics,
      comparability = x$ChemistryComparability
    )
  }
  x$Validation = validatePlantPhytochemistryResult(x)
  x
}

#' Run a staged plant phytochemistry batch workflow
#'
#' @description
#' Runs species-first discovery in resumable species chunks, then optionally
#' resolves compounds through a staged identity or enrichment pass. The default
#' compound stage is `"identity"`, but panels of 100 or more species should
#' first be run with `compound_resolution_profile = "none"` and a local LOTUS
#' index. Identity and richer PubChem enrichment should follow only after the
#' discovery manifest and occurrence evidence have been reviewed.
#'
#' @param plants Character vector or data frame of plant names.
#' @param sources Public provider names passed to `resolvePlantPhytochemistry()`.
#' @param taxon_fallback Fallback ranks passed to `resolvePlantPhytochemistry()`.
#' @param out_dir Optional output directory for checkpoints and CSV/JSON
#' products.
#' @param cache_dir Cache directory. If omitted and `out_dir` is supplied, a
#' `cache` subdirectory under `out_dir` is used.
#' @param species_chunk_size Number of species per discovery chunk.
#' @param compound_resolution_profile One of `"identity"`, `"none"`,
#' `"research"`, or `"full"`. `"identity"` performs a fast PubChem identity
#' pass; `"research"` and `"full"` run richer uafR compound enrichment on the
#' filtered resolution set.
#' @param occurrence_status Occurrence statuses used to select compounds for
#' resolution. Defaults to direct and curated reported records.
#' @param analysis_ready Optional logical used to select compounds for
#' resolution.
#' @param min_confidence Minimum confidence used to select compounds for
#' resolution.
#' @param max_compounds_per_species Optional cap on unique compounds selected
#' per species for compound resolution.
#' @param max_unique_compounds Optional cap on total unique compounds selected
#' for compound resolution.
#' @param chemical_library Optional library passed to rich compound enrichment.
#' @param cache Logical. If `TRUE`, provider requests and batch checkpoints are
#' reused where possible.
#' @param throttle Seconds to wait between uncached provider requests.
#' @param ncbi_email Optional NCBI email.
#' @param ncbi_tool Optional NCBI tool name.
#' @param ncbi_api_key Optional NCBI API key.
#' @param max_pubmed_records Maximum PubMed records per plant.
#' @param max_provider_records Maximum occurrence rows retained per plant from
#' non-PubMed providers.
#' @param lotus_index Optional local LOTUS index as a data frame, flat file
#' path, or manifest-backed lookup directory. Passed to
#' `resolvePlantPhytochemistry()` and used instead of the live LOTUS simple API
#' when `"lotus"` is enabled.
#' @param provider_indexes Optional named list of local provider indexes.
#' `provider_indexes$lotus` and `provider_indexes$npass` are supported.
#' `lotus_index` remains a backward-compatible alias.
#' @param require_all_providers Logical. If `TRUE`, each requested provider must
#' be configured and each completed discovery chunk must report a successful
#' provider stage. Intended for audited production runs after pilot testing.
#' @param request_timeout Maximum seconds allowed for an uncached provider
#' request.
#' @param provider_results Optional mocked or pre-fetched provider results.
#' @param discovery_fun Optional replacement for
#' `resolvePlantPhytochemistry()`. This is primarily intended for deterministic
#' tests and fully local provider workflows.
#' @param enrichment_fun Optional rich enrichment function for tests or cached
#' workflows.
#' @param pubchem_fun Optional replacement for `pubchemProfile()`.
#' @param request_fun Optional provider request function.
#' @param pubtator_request_fun Optional PubTator request function.
#' @param compound_request_fun Optional request function passed to
#' `pubchemProfile()` during compound resolution.
#' @param compound_batch_size Maximum number of unique compounds per compound
#' resolution batch.
#' @param resume Logical. If `TRUE`, reuse discovery and compound checkpoints.
#' @param progress Logical. If `TRUE`, print simple progress messages.
#' @param overwrite Logical. If `TRUE`, replace existing files in `out_dir`.
#' @param defer_derived Logical. If `TRUE`, retain normalized discovery,
#' diagnostics, identity, and operational tables while deferring context
#' linking, summaries, matrices, comparability, identity review, and their
#' filtered exports. This is intended for provider stages that will be merged
#' before analysis; the default preserves the complete standalone result.
#' @param strict Logical passed to validation.
#' @param stop_on_error Logical. If `TRUE`, stop on the first failed discovery
#' chunk; otherwise record the failed chunk and continue.
#' @param allow_large_live_run Logical. Large runs (100 or more species) require
#' persistent checkpoints and, by default, local/injected discovery. Set this
#' to `TRUE` only after a small live pilot and review of
#' `planPlantChemistryRun()`.
#' @param service_busy_pause_threshold Number of consecutive discovery chunks
#' reporting rate-limit or service-busy errors before the run stops scheduling
#' new chunks. Use `Inf` to disable this automatic pause.
#' @param refresh Logical passed to provider adapters.
#' @param ... Additional arguments passed to rich compound enrichment.
#'
#' @return A `uaf_plant_phytochemistry` result with additional
#' `BatchRunManifest`, `BatchChunkManifest`, `FailedQueries`, and `RetryQueue`
#' tables.
#'
#' @export
runPlantPhytochemistryBatch = function(
    plants,
    sources = c("lotus", "knapsack", "npass", "pubchem",
                "pubmed", "pubtator"),
    taxon_fallback = c("species", "genus"),
    out_dir = NULL,
    cache_dir = NULL,
    species_chunk_size = 25,
    compound_resolution_profile = c("identity", "none", "research", "full"),
    occurrence_status = c("direct_reported", "curated_reported"),
    analysis_ready = TRUE,
    min_confidence = "medium",
    max_compounds_per_species = Inf,
    max_unique_compounds = Inf,
    chemical_library = NULL,
    cache = TRUE,
    throttle = 0.5,
    ncbi_email = Sys.getenv("NCBI_EMAIL", ""),
    ncbi_tool = Sys.getenv("NCBI_TOOL", "uafR"),
    ncbi_api_key = Sys.getenv("NCBI_API_KEY", ""),
    max_pubmed_records = 50,
    max_provider_records = max_pubmed_records,
    lotus_index = Sys.getenv("UAFR_LOTUS_INDEX", ""),
    provider_indexes = NULL,
    require_all_providers = FALSE,
    request_timeout = 30,
    provider_results = NULL,
    discovery_fun = NULL,
    enrichment_fun = NULL,
    pubchem_fun = NULL,
    request_fun = NULL,
    pubtator_request_fun = NULL,
    compound_request_fun = NULL,
    compound_batch_size = 50,
    resume = TRUE,
    progress = interactive(),
    overwrite = FALSE,
    defer_derived = FALSE,
    strict = FALSE,
    stop_on_error = FALSE,
    allow_large_live_run = FALSE,
    service_busy_pause_threshold = 1,
    refresh = FALSE,
    ...) {
  compound_resolution_profile = match.arg(compound_resolution_profile)
  started = Sys.time()
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (is.null(cache_dir) && !is.null(out_dir)) {
    cache_dir = file.path(out_dir, "cache")
  }
  cache_dir = .plant_cache_dir(cache_dir)
  checkpoint_dir = if (is.null(out_dir)) {
    file.path(cache_dir, "batch_checkpoints")
  } else {
    file.path(out_dir, "checkpoints")
  }
  if (isTRUE(cache) || !is.null(out_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }
  sources_normalized = .plant_normalize_sources(sources)
  provider_indexes = .plant_provider_indexes(provider_indexes, lotus_index)
  lotus_enabled = "lotus" %in% sources_normalized
  lotus_index = if (lotus_enabled) {
    .plant_lotus_index_or_null(provider_indexes$lotus)
  } else NULL
  if (!is.null(lotus_index) && is.null(.plant_lotus_lookup_info(lotus_index))) {
    lotus_index = standardizeLotusIndex(lotus_index)
  }
  if (!is.null(lotus_index)) provider_indexes$lotus = lotus_index
  provider_index_signature = .plant_provider_indexes_signature(
    provider_indexes, sources_normalized
  )
  provider_resources = plantProviderAvailability(
    sources = sources_normalized,
    provider_indexes = provider_indexes,
    lotus_index = lotus_index
  )
  if (isTRUE(require_all_providers)) {
    unavailable = provider_resources$provider[
      provider_resources$availability_status != "available"
    ]
    if (length(unavailable) > 0) {
      stop("Required plant provider resource(s) unavailable: ",
           paste(unique(unavailable), collapse = ", "), call. = FALSE)
    }
  }

  plant_queries = .plant_queries(plants, taxon_fallback)
  plant_aliases = .plant_query_aliases(plants, plant_queries)
  species = unique(.uaf_non_empty(plant_queries$species))
  if (length(species) < 1) {
    stop("`plants` must contain at least one non-empty species name.",
         call. = FALSE)
  }
  injected_discovery = !is.null(provider_results) || !is.null(discovery_fun)
  if (length(species) >= 100 && !injected_discovery) {
    if (!isTRUE(cache)) {
      stop("Runs with 100 or more species require `cache = TRUE`.",
           call. = FALSE)
    }
    if (is.null(out_dir)) {
      stop(paste(
        "Runs with 100 or more species require `out_dir` so incremental",
        "manifests, failed queries, and retry queues are preserved."
      ), call. = FALSE)
    }
    if (lotus_enabled && is.null(lotus_index) &&
        !isTRUE(allow_large_live_run)) {
      stop(paste(
        "Large live LOTUS simple-search runs are disabled by default.",
        "Supply a local `lotus_index`, or complete a small pilot and set",
        "`allow_large_live_run = TRUE` explicitly."
      ), call. = FALSE)
    }
    local_sources = character()
    if (!is.null(lotus_index)) local_sources = c(local_sources, "lotus")
    if (!is.null(.plant_npass_index_or_null(provider_indexes$npass))) {
      local_sources = c(local_sources, "npass")
    }
    live_sources = setdiff(sources_normalized, local_sources)
    if (length(live_sources) > 0 && !isTRUE(allow_large_live_run)) {
      stop(paste0(
        "Large live provider stages are disabled by default (selected: ",
        paste(live_sources, collapse = ", "), "). Run local LOTUS discovery",
        " first, then use targeted follow-up panels; or set",
        " `allow_large_live_run = TRUE` after reviewing the preflight plan."
      ), call. = FALSE)
    }
  }
  chunk_size = .plant_enrichment_batch_size(species_chunk_size,
                                            length(species))
  species_chunks = split(species, ceiling(seq_along(species) / chunk_size))
  species_chunk_inputs = .plant_batch_chunk_inputs(
    plants, plant_queries, species_chunks
  )
  chunk_results = vector("list", length(species_chunks))
  run_signature = .plant_batch_run_signature(
    species = species,
    sources = sources_normalized,
    taxon_fallback = taxon_fallback,
    provider_index_signature = provider_index_signature,
    max_pubmed_records = max_pubmed_records,
    max_provider_records = max_provider_records,
    provider_results = provider_results,
    query_signature = .plant_batch_query_signature(plant_queries,
                                                   plant_aliases)
  )
  checkpoint_files = vapply(seq_along(species_chunks), function(i) {
    chunk_signature = .plant_batch_chunk_signature(
      run_signature, i, species_chunks[[i]]
    )
    file.path(
      checkpoint_dir,
      paste0("discovery_v2_chunk_", sprintf("%04d", i), "_",
             chunk_signature, ".rds")
    )
  }, character(1))
  chunk_manifest = .plant_batch_manifest_plan(
    species_chunks = species_chunks,
    checkpoint_files = checkpoint_files,
    run_signature = run_signature,
    out_dir = out_dir,
    cache = cache
  )
  chunk_manifest = .plant_batch_merge_prior_manifest(
    chunk_manifest, out_dir, run_signature, resume, overwrite
  )
  .plant_write_batch_operational_files(chunk_manifest, out_dir)
  discovery_fun = discovery_fun %||% resolvePlantPhytochemistry
  if (!is.function(discovery_fun)) {
    stop("`discovery_fun` must be a function when supplied.", call. = FALSE)
  }
  busy_threshold = suppressWarnings(as.numeric(
    service_busy_pause_threshold[[1]]
  ))
  if (is.na(busy_threshold) || busy_threshold < 1) {
    stop("`service_busy_pause_threshold` must be at least 1 or `Inf`.",
         call. = FALSE)
  }
  consecutive_busy_chunks = 0L
  pause_reason = NA_character_

  for (i in seq_along(species_chunks)) {
    chunk_species = species_chunks[[i]]
    chunk_input = species_chunk_inputs[[i]]
    chunk_provider_results = .plant_batch_subset_provider_results(
      provider_results, chunk_species
    )
    checkpoint_file = checkpoint_files[[i]]
    chunk_signature = chunk_manifest$chunk_signature[[i]]
    chunk_started = Sys.time()
    status = "completed"
    error_message = NA_character_
    result = NULL
    checkpoint = list(ok = FALSE, reason = "not_checked", result = NULL)
    if (isTRUE(cache) && isTRUE(resume) && !isTRUE(refresh)) {
      checkpoint = .plant_read_batch_checkpoint(
        checkpoint_file,
        expected_run_signature = run_signature,
        expected_chunk_signature = chunk_signature,
        expected_species = chunk_species
      )
    }
    chunk_manifest$checkpoint_read_status[[i]] = checkpoint$reason
    if (isTRUE(checkpoint$ok)) {
      if (isTRUE(progress)) {
        message("uafR plant discovery chunk ", i, "/",
                length(species_chunks), ": using cached result")
      }
      result = checkpoint$result
      counts = .plant_batch_result_counts(result)
      chunk_results[[i]] = result
      chunk_manifest$status[[i]] = "completed"
      chunk_manifest$checkpoint_status[[i]] = "cache_hit"
      chunk_manifest$checkpoint_cache_hit[[i]] = "Yes"
      chunk_manifest$started_at[[i]] = format(
        chunk_started, "%Y-%m-%dT%H:%M:%S%z"
      )
      chunk_manifest$finished_at[[i]] = .plant_timestamp()
      chunk_manifest$completed_at[[i]] = chunk_manifest$finished_at[[i]]
      chunk_manifest$elapsed_seconds[[i]] = round(as.numeric(difftime(
        Sys.time(), chunk_started, units = "secs"
      )), 3)
      chunk_manifest$cache_hit_count[[i]] = counts$cache_hit_count + 1L
      chunk_manifest$request_count[[i]] = counts$request_count
      chunk_manifest$occurrence_count[[i]] = counts$occurrence_count
      chunk_manifest$literature_candidate_count[[i]] =
        counts$literature_candidate_count
      chunk_manifest$provider_diagnostic_count[[i]] =
        counts$provider_diagnostic_count
      chunk_manifest$error_count[[i]] = counts$error_count
      chunk_manifest$warning_count[[i]] = counts$warning_count
      .plant_write_batch_operational_files(chunk_manifest, out_dir)
      consecutive_busy_chunks = 0L
      next
    }

    if (file.exists(checkpoint_file) && checkpoint$reason != "not_checked") {
      chunk_manifest$checkpoint_status[[i]] = paste0(
        "rebuild_", checkpoint$reason
      )
    } else if (!isTRUE(cache)) {
      chunk_manifest$checkpoint_status[[i]] = "cache_disabled"
    } else {
      chunk_manifest$checkpoint_status[[i]] = "cache_miss"
    }
    prior_retry = tolower(.uaf_first_non_empty_text(
      chunk_manifest$previous_status[[i]], ""
    )) %in% c("failed", "error", "timeout", "timed_out", "rate_limited",
              "incomplete", "running", "started", "stopped", "retry")
    invalid_checkpoint = file.exists(checkpoint_file) &&
      checkpoint$reason != "not_checked" && checkpoint$reason != "missing"
    if (prior_retry || invalid_checkpoint) {
      chunk_manifest$retry_count[[i]] =
        suppressWarnings(as.integer(chunk_manifest$retry_count[[i]])) + 1L
    }
    chunk_manifest$status[[i]] = "running"
    chunk_manifest$started_at[[i]] = format(
      chunk_started, "%Y-%m-%dT%H:%M:%S%z"
    )
    .plant_write_batch_operational_files(chunk_manifest, out_dir)

    if (isTRUE(progress)) {
      message("uafR plant discovery chunk ", i, "/",
              length(species_chunks), ": querying ", length(chunk_species),
              " species")
    }
    captured_condition = NULL
    was_interrupted = FALSE
    discovery_args = list(
      plants = chunk_input,
      sources = sources,
      taxon_fallback = taxon_fallback,
      enrich_compounds = FALSE,
      detail = "none",
      cache = cache,
      cache_dir = cache_dir,
      throttle = throttle,
      ncbi_email = ncbi_email,
      ncbi_tool = ncbi_tool,
      ncbi_api_key = ncbi_api_key,
      max_pubmed_records = max_pubmed_records,
      max_provider_records = max_provider_records,
      lotus_index = lotus_index,
      provider_indexes = provider_indexes,
      require_all_providers = require_all_providers,
      request_timeout = request_timeout,
      provider_results = chunk_provider_results,
      request_fun = request_fun,
      pubtator_request_fun = pubtator_request_fun,
      progress = progress,
      refresh = refresh
    )
    discovery_formals = names(formals(discovery_fun))
    if ("defer_derived" %in% discovery_formals ||
        "..." %in% discovery_formals) {
      discovery_args$defer_derived = TRUE
    }
    result = tryCatch(
      do.call(discovery_fun, discovery_args),
      interrupt = function(condition) {
        captured_condition <<- condition
        was_interrupted <<- TRUE
        NULL
      },
      error = function(condition) {
        captured_condition <<- condition
        NULL
      }
    )
    if (is.null(captured_condition) &&
        !.plant_valid_discovery_result(result)) {
      captured_condition = simpleError(paste(
        "Discovery function returned an invalid result; required tables are",
        "PlantQueries, ProviderDiagnostics, PlantCompoundOccurrences, and",
        "LiteratureCandidates."
      ))
    }

    health = NULL
    if (!is.null(captured_condition)) {
      error_message = .plant_redact_secrets(
        conditionMessage(captured_condition)
      )
      status = if (isTRUE(was_interrupted)) {
        "stopped"
      } else if (.plant_service_busy_message(error_message)) {
        "rate_limited"
      } else {
        "failed"
      }
      result = .plant_empty_batch_result(chunk_input, taxon_fallback,
                                         error_message)
      chunk_manifest$checkpoint_status[[i]] = "not_written_failure"
    } else {
      health = .plant_batch_result_health(result)
      if (health$error_count > 0) {
        status = if (isTRUE(health$service_busy)) "rate_limited" else "failed"
        error_message = health$error_message
        chunk_manifest$checkpoint_status[[i]] = "not_written_provider_error"
      } else if (isTRUE(cache)) {
        envelope = .plant_batch_checkpoint_envelope(
          result = result,
          run_signature = run_signature,
          chunk_signature = chunk_signature,
          species = chunk_species
        )
        save_error = tryCatch({
          .plant_atomic_save_rds(envelope, checkpoint_file)
          NULL
        }, error = function(condition) condition)
        if (inherits(save_error, "condition")) {
          status = "failed"
          error_message = .plant_redact_secrets(paste(
            "Checkpoint write failed:", conditionMessage(save_error)
          ))
          chunk_manifest$checkpoint_status[[i]] = "write_failed"
        } else {
          chunk_manifest$checkpoint_status[[i]] = "written"
        }
      } else {
        chunk_manifest$checkpoint_status[[i]] = "cache_disabled"
      }
    }
    chunk_results[[i]] = result
    counts = .plant_batch_result_counts(result)
    elapsed = round(as.numeric(difftime(Sys.time(), chunk_started,
                                        units = "secs")), 3)
    chunk_manifest$status[[i]] = status
    chunk_manifest$error_message[[i]] = error_message
    chunk_manifest$elapsed_seconds[[i]] = elapsed
    chunk_manifest$finished_at[[i]] = .plant_timestamp()
    chunk_manifest$completed_at[[i]] = chunk_manifest$finished_at[[i]]
    chunk_manifest$cache_hit_count[[i]] = counts$cache_hit_count
    chunk_manifest$request_count[[i]] = counts$request_count
    chunk_manifest$occurrence_count[[i]] = counts$occurrence_count
    chunk_manifest$literature_candidate_count[[i]] =
      counts$literature_candidate_count
    chunk_manifest$provider_diagnostic_count[[i]] =
      counts$provider_diagnostic_count
    chunk_manifest$error_count[[i]] = counts$error_count
    chunk_manifest$warning_count[[i]] = counts$warning_count
    .plant_write_batch_operational_files(chunk_manifest, out_dir)

    if (status == "rate_limited") {
      consecutive_busy_chunks = consecutive_busy_chunks + 1L
    } else {
      consecutive_busy_chunks = 0L
    }
    if (isTRUE(was_interrupted) ||
        (isTRUE(stop_on_error) && status %in%
           c("failed", "rate_limited", "stopped"))) {
      if (isTRUE(was_interrupted) && !is.null(captured_condition)) {
        stop(captured_condition)
      }
      stop(simpleError(error_message))
    }
    if (is.finite(busy_threshold) &&
        consecutive_busy_chunks >= busy_threshold) {
      pause_reason = paste0(
        "Paused after ", consecutive_busy_chunks,
        " consecutive rate-limit/service-busy discovery chunks. Resume later",
        " with the same out_dir and cache_dir."
      )
      if (isTRUE(progress)) message(pause_reason)
      break
    }
  }

  combined = .plant_combine_batch_results(
    chunk_results, plant_queries, taxon_fallback,
    link_source_context = !isTRUE(defer_derived),
    defer_derived = defer_derived
  )
  combined$BatchChunkManifest = chunk_manifest
  combined$FailedQueries = .plant_batch_failed_queries(chunk_manifest)
  combined$RetryQueue = writePlantChemistryRetryQueue(chunk_manifest)
  discovery_complete = all(chunk_manifest$status == "completed")

  resolution_occurrences = .plant_batch_resolution_occurrences(
    combined$PlantCompoundOccurrences,
    occurrence_status = occurrence_status,
    analysis_ready = analysis_ready,
    min_confidence = min_confidence,
    max_compounds_per_species = max_compounds_per_species,
    max_unique_compounds = max_unique_compounds
  )

  compound_provenance = .plant_provenance(
    "compound_resolution", compound_resolution_profile,
    paste(unique(.uaf_non_empty(resolution_occurrences$compound_name)),
          collapse = "; "),
    NA_character_, length(unique(.uaf_non_empty(
      resolution_occurrences$compound_name_clean))),
    if (discovery_complete) {
      "Compound resolution stage was not run."
    } else {
      "Compound resolution was deferred because discovery is incomplete."
    }
  )
  categorate_result = NULL
  trait_evidence = .uaf_empty_table(.plant_trait_evidence_cols())
  compound_resolution = .plant_compound_resolution(
    combined$PlantCompoundOccurrences, NULL
  )
  if (!discovery_complete || compound_resolution_profile == "none") {
    compound_resolution =
      .plant_mark_all_compounds_not_attempted(compound_resolution)
    if (!discovery_complete && nrow(compound_resolution) > 0) {
      compound_resolution$notes = paste(
        "Compound resolution was deferred because one or more discovery",
        "chunks are failed, stopped, rate-limited, or not started."
      )
    }
  } else {
    compound_resolution = .plant_mark_unattempted_compounds(
      compound_resolution, resolution_occurrences
    )
  }

  compound_stage_status = if (!discovery_complete) {
    "deferred_incomplete_discovery"
  } else if (compound_resolution_profile == "none") {
    "not_requested"
  } else if (nrow(resolution_occurrences) < 1) {
    "no_compounds_selected"
  } else {
    "completed"
  }

  if (discovery_complete && compound_resolution_profile == "identity" &&
      nrow(resolution_occurrences) > 0) {
    identity = .plant_resolve_compound_identities(
      occurrences = resolution_occurrences,
      cache = cache,
      cache_dir = file.path(cache_dir, "compound_identity"),
      throttle = throttle,
      batch_size = compound_batch_size,
      resume = resume,
      progress = progress,
      pubchem_fun = pubchem_fun,
      request_fun = compound_request_fun,
      lotus_index = lotus_index,
      source_identity = combined$SourceCompoundIdentity
    )
    categorate_result = identity$CategorateResult
    compound_resolution = .plant_merge_compound_resolution(
      compound_resolution, identity$CompoundResolution
    )
    compound_provenance = identity$Provenance
  } else if (discovery_complete &&
             compound_resolution_profile %in% c("research", "full") &&
             nrow(resolution_occurrences) > 0) {
    enrichment = enrichPlantCompounds(
      resolution_occurrences,
      chemical_library = chemical_library,
      detail = compound_resolution_profile,
      cache = cache,
      cache_dir = file.path(cache_dir, "compound_enrichment"),
      throttle = throttle,
      enrichment_fun = enrichment_fun,
      pubchem_fun = pubchem_fun,
      compound_request_fun = compound_request_fun,
      batch_size = compound_batch_size,
      resume = resume,
      progress = progress,
      ...
    )
    categorate_result = enrichment$CategorateResult
    compound_resolution = .plant_merge_compound_resolution(
      compound_resolution, enrichment$CompoundResolution
    )
    trait_evidence = enrichment$TraitEvidence
    compound_provenance = enrichment$Provenance
  }

  combined$CategorateResult = categorate_result
  combined$CompoundResolution = compound_resolution
  combined$TraitEvidence = trait_evidence
  combined$Provenance = .plant_bind_tables(
    list(combined$Provenance, compound_provenance),
    .plant_provenance_cols()
  )
  if (!isTRUE(defer_derived)) {
    combined$ChemistryComparability = plantChemistryComparability(
      list(PlantCompoundOccurrences = combined$PlantCompoundOccurrences,
           CategorateResult = categorate_result),
      min_confidence = "low"
    )
    combined$SpeciesChemistrySummary = summarizePlantPhytochemistry(
      plant_compounds = combined$PlantCompoundOccurrences,
      categorate_result = categorate_result,
      compound_resolution = compound_resolution,
      plant_queries = combined$PlantQueries,
      provider_diagnostics = combined$ProviderDiagnostics,
      comparability = combined$ChemistryComparability
    )
    combined$SpeciesChemistryMatrix = plantPhytochemistryMatrix(
      list(PlantCompoundOccurrences = combined$PlantCompoundOccurrences,
           CategorateResult = categorate_result),
      level = "species",
      profile = "core",
      mode = "binary",
      min_confidence = min_confidence,
      max_traits = .plant_automatic_matrix_max_traits()
    )
    combined$ComparableChemistryMatrix = plantComparableChemistryMatrix(
      list(ChemistryComparability = combined$ChemistryComparability),
      comparison_scope = "specialized_metabolites",
      level = "species",
      mode = "binary",
      min_comparability_confidence = min_confidence
    )
    combined$CompoundIdentityReview = plantCompoundIdentityReviewTable(combined)
  }

  completed = Sys.time()
  run_manifest = .plant_batch_run_manifest(
    result = combined,
    chunk_manifest = chunk_manifest,
    started = started,
    completed = completed,
    sources = sources,
    compound_resolution_profile = compound_resolution_profile,
    resolution_occurrences = resolution_occurrences,
    out_dir = out_dir,
    cache_dir = cache_dir,
    discovery_complete = discovery_complete,
    compound_stage_status = compound_stage_status,
    run_signature = run_signature,
    pause_reason = pause_reason
  )
  combined$BatchRunManifest = run_manifest
  if (isTRUE(defer_derived)) {
    combined$Validation = list(
      Summary = data.frame(), TableQuality = data.frame(),
      ProviderDiagnostics = data.frame(), Issues = data.frame(),
      DataDictionary = data.frame()
    )
    combined$BatchRunManifest$validation_status = "deferred"
  } else {
    combined$Validation = validatePlantPhytochemistryResult(combined,
                                                            strict = strict)
    combined$BatchRunManifest$validation_status =
      .uaf_first_non_empty_text(combined$Validation$Summary$Status)
  }
  class(combined) = unique(c("uaf_plant_phytochemistry", class(combined)))
  if (!is.null(out_dir)) {
    combined$BatchExportManifest = .plant_write_batch_outputs(
      combined, out_dir, overwrite,
      include_analysis = !isTRUE(defer_derived)
    )
  }
  combined
}

#' Recommended plant panels for phytochemistry pilot runs
#'
#' @description
#' Returns small species panels useful for calibrating species-first plant
#' phytochemistry discovery before scaling. The panels are not claims that each
#' plant has complete public chemistry coverage; they are designed to include a
#' mix of expected data depth, project relevance, and likely edge cases.
#'
#' @param profile One of `"general"` or `"remediation"`.
#'
#' @return Data frame with `species`, `panel_role`, `expected_data_depth`,
#' `review_focus`, and `rationale` columns.
#'
#' @export
plantPhytochemistryPilotPanel = function(profile = c("general",
                                                     "remediation")) {
  profile = match.arg(profile)
  general = data.frame(
    species = c("Camellia sinensis", "Zea mays", "Arabidopsis thaliana",
                "Salix nigra", "Mentha piperita",
                "Lavandula angustifolia", "Salvia rosmarinus",
                "Ocimum basilicum", "Glycine max", "Oryza sativa",
                "Populus deltoides", "Panicum virgatum",
                "Phragmites australis", "Typha latifolia",
                "Spartina alterniflora"),
    panel_role = c("well_studied_crop_or_model",
                   "well_studied_crop_or_model",
                   "well_studied_crop_or_model",
                   "natural_products_and_ecology",
                   "volatile_natural_products",
                   "volatile_natural_products",
                   "volatile_natural_products",
                   "volatile_natural_products",
                   "well_studied_crop_or_model",
                   "well_studied_crop_or_model",
                   "environmental_or_remediation_relevance",
                   "environmental_or_remediation_relevance",
                   "environmental_or_remediation_relevance",
                   "environmental_or_remediation_relevance",
                   "environmental_or_remediation_relevance"),
    expected_data_depth = c("high", "high", "high", "moderate", "high",
                            "high", "moderate", "high", "high", "high",
                            "moderate", "moderate", "variable",
                            "variable", "variable"),
    review_focus = c("source coverage and primary/specialized separation",
                     "benzoxazinoid and primary-metabolism separation",
                     "model-plant literature and candidate filtering",
                     "direct species versus genus fallback evidence",
                     "volatile-specialized chemistry and method context",
                     "volatile-specialized chemistry and method context",
                     "synonyms and volatile-specialized chemistry",
                     "volatile-specialized chemistry and leaf context",
                     "primary versus specialized chemistry",
                     "primary versus specialized chemistry",
                     "direct species records and tissue context",
                     "remediation relevance and source sparsity",
                     "rhizosphere/root context and source sparsity",
                     "aquatic/wetland context and source sparsity",
                     "coastal/wetland context and source sparsity"),
    rationale = c(
      "High-publication plant useful for checking rich database/literature behavior.",
      "Well-studied crop useful for testing primary and specialized chemistry separation.",
      "Model plant useful for evaluating literature-candidate behavior.",
      "Natural-products/ecology species useful for fallback and context review.",
      "Essential-oil plant useful for volatile chemistry matrices.",
      "Essential-oil plant useful for volatile chemistry matrices.",
      "Aromatic plant with synonym risk useful for name/provider review.",
      "Aromatic plant useful for leaf and volatile chemistry context.",
      "Crop species useful for primary-metabolism and source coverage checks.",
      "Crop species useful for primary-metabolism and source coverage checks.",
      "Woody plant useful for tissue/context and environmental comparison.",
      "Grass species useful for environmental and remediation-focused pilots.",
      "Wetland grass useful for sparse public-record and rhizosphere checks.",
      "Wetland plant useful for sparse public-record and context checks.",
      "Coastal grass useful for sparse public-record and context checks."
    ),
    stringsAsFactors = FALSE
  )
  if (profile == "general") return(general)
  general[general$panel_role %in%
            c("environmental_or_remediation_relevance",
              "well_studied_crop_or_model",
              "natural_products_and_ecology"), , drop = FALSE]
}

#' Run a plant phytochemistry pilot workflow
#'
#' @description
#' `runPlantPhytochemistryPilot()` is a higher-level workflow for small live
#' provider pilots before scaling to hundreds of species. It runs
#' `runPlantPhytochemistryBatch()`, writes the full audit bundle, and adds
#' pilot-specific files that make provider coverage, review needs, biological
#' context, and comparable-chemistry matrices easy to inspect.
#'
#' The pilot does not claim complete plant metabolomes. It separates direct
#' species evidence, fallback/candidate evidence, source-backed biological
#' context, and chemistry-comparison scope so users can decide which records are
#' suitable for downstream analysis.
#'
#' @param plants Character vector or data frame of plant names.
#' @param out_dir Output directory for audit tables, pilot reports, and
#' matrices.
#' @param sources Provider sources passed to `runPlantPhytochemistryBatch()`.
#' @param taxon_fallback Taxon fallback ranks passed to
#' `runPlantPhytochemistryBatch()`.
#' @param cache_dir Cache directory. Defaults to `file.path(out_dir, "cache")`.
#' @param species_chunk_size Number of species per discovery chunk.
#' @param compound_resolution_profile One of `"identity"`, `"none"`,
#' `"research"`, or `"full"`.
#' @param occurrence_status Occurrence statuses retained for compound
#' resolution.
#' @param analysis_ready Logical filter for compound-resolution input.
#' @param min_confidence Minimum confidence for resolution and matrices.
#' @param max_compounds_per_species Maximum analysis-ready compounds resolved
#' per species.
#' @param max_unique_compounds Maximum unique compounds resolved across the
#' pilot.
#' @param chemical_library Optional library passed to rich enrichment.
#' @param cache Logical. If `TRUE`, use cache/checkpoint files.
#' @param throttle Seconds to wait between uncached provider requests.
#' @param ncbi_email Optional NCBI email.
#' @param ncbi_tool Optional NCBI tool name.
#' @param ncbi_api_key Optional NCBI API key. It is not written to output
#' files.
#' @param max_pubmed_records Maximum PubMed records per plant.
#' @param max_provider_records Maximum non-PubMed provider records per plant.
#' @param lotus_index Optional local LOTUS index as a data frame, flat file
#' path, or manifest-backed lookup directory. Passed to
#' `runPlantPhytochemistryBatch()` and used instead of the live LOTUS simple API
#' when `"lotus"` is enabled.
#' @param request_timeout Maximum seconds allowed for an uncached provider
#' request.
#' @param provider_results Optional mocked or pre-fetched provider results.
#' @param enrichment_fun Optional rich enrichment function for tests or cached
#' workflows.
#' @param pubchem_fun Optional replacement for `pubchemProfile()`.
#' @param request_fun Optional provider request function.
#' @param pubtator_request_fun Optional PubTator request function.
#' @param compound_request_fun Optional request function for compound
#' resolution.
#' @param compound_batch_size Maximum compounds per compound-resolution batch.
#' @param matrix_feature Matrix feature family passed to
#' `plantComparableChemistryMatrix()`.
#' @param matrix_mode Matrix mode passed to `plantComparableChemistryMatrix()`.
#' @param max_matrix_features Maximum feature columns retained per pilot
#' matrix.
#' @param resume Logical. If `TRUE`, reuse checkpoints where available.
#' @param progress Logical. If `TRUE`, print progress messages.
#' @param overwrite Logical. If `TRUE`, replace existing output files.
#' @param strict Logical passed to validation.
#' @param stop_on_error Logical. If `TRUE`, stop on first failed chunk.
#' @param refresh Logical passed to provider adapters.
#' @param ... Additional arguments passed to rich enrichment.
#'
#' @return A `uaf_plant_phytochemistry` result with `PilotSummary`,
#' `PilotReviewNeeded`, `PilotMatrices`, and `PilotExportManifest` tables.
#'
#' @export
runPlantPhytochemistryPilot = function(
    plants,
    out_dir,
    sources = c("lotus", "knapsack", "npass", "pubchem",
                "pubmed", "pubtator"),
    taxon_fallback = c("species", "genus"),
    cache_dir = NULL,
    species_chunk_size = 25,
    compound_resolution_profile = c("identity", "none", "research", "full"),
    occurrence_status = c("direct_reported", "curated_reported"),
    analysis_ready = TRUE,
    min_confidence = "medium",
    max_compounds_per_species = Inf,
    max_unique_compounds = Inf,
    chemical_library = NULL,
    cache = TRUE,
    throttle = 0.5,
    ncbi_email = Sys.getenv("NCBI_EMAIL", ""),
    ncbi_tool = Sys.getenv("NCBI_TOOL", "uafR"),
    ncbi_api_key = Sys.getenv("NCBI_API_KEY", ""),
    max_pubmed_records = 25,
    max_provider_records = 100,
    lotus_index = Sys.getenv("UAFR_LOTUS_INDEX", ""),
    request_timeout = 30,
    provider_results = NULL,
    enrichment_fun = NULL,
    pubchem_fun = NULL,
    request_fun = NULL,
    pubtator_request_fun = NULL,
    compound_request_fun = NULL,
    compound_batch_size = 50,
    matrix_feature = c("comparison_group", "biosynthetic_family",
                       "chemical_behavior", "compound",
                       "comparison_scope"),
    matrix_mode = c("binary", "count", "confidence"),
    max_matrix_features = Inf,
    resume = TRUE,
    progress = interactive(),
    overwrite = FALSE,
    strict = FALSE,
    stop_on_error = FALSE,
    refresh = FALSE,
    ...) {
  compound_resolution_profile = match.arg(compound_resolution_profile)
  matrix_feature = match.arg(matrix_feature)
  matrix_mode = match.arg(matrix_mode)
  if (missing(out_dir) || is.null(out_dir) ||
      length(.uaf_non_empty(out_dir)) < 1) {
    stop("`out_dir` is required.", call. = FALSE)
  }
  out_dir = .uaf_non_empty(out_dir)[[1]]
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(cache_dir)) cache_dir = file.path(out_dir, "cache")

  result = runPlantPhytochemistryBatch(
    plants = plants,
    sources = sources,
    taxon_fallback = taxon_fallback,
    out_dir = out_dir,
    cache_dir = cache_dir,
    species_chunk_size = species_chunk_size,
    compound_resolution_profile = compound_resolution_profile,
    occurrence_status = occurrence_status,
    analysis_ready = analysis_ready,
    min_confidence = min_confidence,
    max_compounds_per_species = max_compounds_per_species,
    max_unique_compounds = max_unique_compounds,
    chemical_library = chemical_library,
    cache = cache,
    throttle = throttle,
    ncbi_email = ncbi_email,
    ncbi_tool = ncbi_tool,
    ncbi_api_key = ncbi_api_key,
    max_pubmed_records = max_pubmed_records,
    max_provider_records = max_provider_records,
    lotus_index = lotus_index,
    request_timeout = request_timeout,
    provider_results = provider_results,
    enrichment_fun = enrichment_fun,
    pubchem_fun = pubchem_fun,
    request_fun = request_fun,
    pubtator_request_fun = pubtator_request_fun,
    compound_request_fun = compound_request_fun,
    compound_batch_size = compound_batch_size,
    resume = resume,
    progress = progress,
    overwrite = overwrite,
    strict = strict,
    stop_on_error = stop_on_error,
    refresh = refresh,
    ...
  )

  pilot_summary = .plant_pilot_summary(result)
  review_needed = plantPhytochemistryReviewTable(result)
  pilot_matrices = .plant_pilot_matrices(
    result,
    min_confidence = min_confidence,
    feature = matrix_feature,
    mode = matrix_mode,
    max_features = max_matrix_features
  )
  result$PilotSummary = pilot_summary
  result$PilotReviewNeeded = review_needed
  result$PilotMatrices = pilot_matrices
  qa_report = plantPhytochemistryQAReport(result)
  result$PilotQAReport = qa_report
  pilot_manifest = .plant_write_pilot_outputs(
    result = result,
    pilot_summary = pilot_summary,
    review_needed = review_needed,
    pilot_matrices = pilot_matrices,
    qa_report = qa_report,
    out_dir = out_dir,
    overwrite = overwrite
  )

  result$PilotExportManifest = pilot_manifest
  class(result) = unique(c("uaf_plant_phytochemistry", class(result)))
  result
}

#' Summarize pilot readiness checks
#'
#' @description
#' Builds a compact QA report for a plant phytochemistry result or pilot result.
#' The report is intended for deciding whether a species panel is ready for
#' larger-scale runs, manual review, or provider-specific hardening. It does not
#' certify complete chemistry coverage.
#'
#' @param x Plant phytochemistry result.
#' @param min_species_with_records_fraction Minimum fraction of species that
#' should have at least one public or curated occurrence row.
#' @param min_analysis_ready_species_fraction Minimum fraction of species that
#' should have at least one analysis-ready compound.
#' @param min_context_known_fraction Minimum overall context-known occurrence
#' fraction.
#' @param min_resolved_fraction Minimum fraction of attempted compound
#' identities that should resolve.
#' @param max_review_occurrence_fraction Maximum acceptable review-required
#' occurrence fraction before the pilot is considered review-heavy.
#'
#' @return QA report data frame with check, status, value, threshold, details,
#' and recommendation columns.
#'
#' @export
plantPhytochemistryQAReport = function(
    x,
    min_species_with_records_fraction = 0.50,
    min_analysis_ready_species_fraction = 0.25,
    min_context_known_fraction = 0.25,
    min_resolved_fraction = 0.50,
    max_review_occurrence_fraction = 1.00) {
  if (!is.list(x) || is.data.frame(x)) {
    stop("`x` must be a plant phytochemistry result list.", call. = FALSE)
  }
  summary = if (is.data.frame(x$PilotSummary)) {
    x$PilotSummary
  } else {
    .plant_pilot_summary(x)
  }
  occurrences = .plant_normalize_occurrences(x$PlantCompoundOccurrences)
  resolution = if (is.data.frame(x$CompoundResolution)) {
    x$CompoundResolution
  } else {
    .uaf_empty_table(.plant_compound_resolution_cols())
  }
  diagnostics = if (is.data.frame(x$ProviderDiagnostics)) {
    x$ProviderDiagnostics
  } else {
    .uaf_empty_table(.plant_provider_diagnostic_cols())
  }
  review = if (is.data.frame(x$PilotReviewNeeded)) {
    x$PilotReviewNeeded
  } else {
    plantPhytochemistryReviewTable(x)
  }
  matrices = if (is.list(x$PilotMatrices)) x$PilotMatrices else list()
  species_n = max(1, nrow(summary))
  species_with_records = sum(summary$compound_count > 0, na.rm = TRUE)
  analysis_ready_species =
    sum(summary$analysis_ready_compound_count > 0, na.rm = TRUE)
  context_known = if (nrow(occurrences) > 0) {
    sum(occurrences$biological_context_status != "context_missing",
        na.rm = TRUE) / nrow(occurrences)
  } else {
    0
  }
  attempted = resolution[resolution$resolution_source != "not_attempted", ,
                         drop = FALSE]
  resolved_fraction = if (nrow(attempted) > 0) {
    sum(attempted$resolved %in% TRUE, na.rm = TRUE) / nrow(attempted)
  } else {
    NA_real_
  }
  review_fraction = if (nrow(occurrences) > 0) {
    nrow(review) / nrow(occurrences)
  } else {
    0
  }
  validation_status = if (is.list(x$Validation) &&
                          is.data.frame(x$Validation$Summary) &&
                          "Status" %in% names(x$Validation$Summary)) {
    .uaf_first_non_empty_text(x$Validation$Summary$Status)
  } else {
    NA_character_
  }
  provider_errors = sum(suppressWarnings(as.numeric(diagnostics$error_count)),
                        na.rm = TRUE)
  provider_records = sum(suppressWarnings(as.numeric(diagnostics$record_count)),
                         na.rm = TRUE)
  matrix_nonempty = vapply(matrices, function(item) {
    is.data.frame(item) && nrow(item) > 0 && ncol(item) > 1
  }, logical(1))
  rows = list(
    .plant_qa_row(
      "validation_status",
      ifelse(identical(validation_status, "pass"), "pass",
             ifelse(identical(validation_status, "warning"),
                    "warning", "fail")),
      validation_status,
      "pass",
      "Package-level validation status for the assembled result.",
      "Inspect Validation$Issues before interpreting matrices."
    ),
    .plant_qa_row(
      "provider_errors",
      ifelse(provider_errors == 0, "pass", "warning"),
      provider_errors,
      "0",
      paste("Provider records:", provider_records),
      "Inspect ProviderDiagnostics for timeouts, unavailable providers, or parser failures."
    ),
    .plant_qa_row(
      "species_with_records_fraction",
      .plant_qa_threshold_status(species_with_records / species_n,
                                 min_species_with_records_fraction),
      round(species_with_records / species_n, 3),
      paste0(">=", min_species_with_records_fraction),
      paste(species_with_records, "of", species_n,
            "species have occurrence records."),
      "Try alternate species names or curated intake for no-record species."
    ),
    .plant_qa_row(
      "analysis_ready_species_fraction",
      .plant_qa_threshold_status(analysis_ready_species / species_n,
                                 min_analysis_ready_species_fraction),
      round(analysis_ready_species / species_n, 3),
      paste0(">=", min_analysis_ready_species_fraction),
      paste(analysis_ready_species, "of", species_n,
            "species have at least one analysis-ready compound."),
      "Review candidate/fallback rows before scaling or running rich enrichment."
    ),
    .plant_qa_row(
      "context_known_fraction",
      .plant_qa_threshold_status(context_known,
                                 min_context_known_fraction),
      round(context_known, 3),
      paste0(">=", min_context_known_fraction),
      paste(round(100 * context_known, 1),
            "percent of occurrence rows have plant-part, tissue, or method context."),
      "Use PlantContextEvidence and source records to improve context-sparse rows."
    ),
    .plant_qa_row(
      "identity_resolved_fraction",
      if (is.na(resolved_fraction)) "warning" else
        .plant_qa_threshold_status(resolved_fraction, min_resolved_fraction),
      ifelse(is.na(resolved_fraction), NA, round(resolved_fraction, 3)),
      paste0(">=", min_resolved_fraction),
      ifelse(is.na(resolved_fraction),
             "No compound identity resolution was attempted.",
             paste(round(100 * resolved_fraction, 1),
                   "percent of attempted compounds resolved.")),
      "Use identity resolution before rich enrichment; inspect unresolved names."
    ),
    .plant_qa_row(
      "review_burden_fraction",
      ifelse(review_fraction <= max_review_occurrence_fraction, "pass",
             "warning"),
      round(review_fraction, 3),
      paste0("<=", max_review_occurrence_fraction),
      paste(nrow(review), "review rows for", nrow(occurrences),
            "occurrence rows."),
      "High review burden means provider output needs curation or parser hardening before scale-up."
    ),
    .plant_qa_row(
      "nonempty_context_matrices",
      ifelse(any(matrix_nonempty), "pass", "warning"),
      sum(matrix_nonempty),
      ">=1",
      paste(sum(matrix_nonempty), "of", length(matrix_nonempty),
            "pilot comparable matrices have feature columns."),
      "If matrices are empty, inspect ChemistryComparability and context filters."
    )
  )
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

#' Summarize plant phytochemistry
#'
#' @param plant_compounds Plant-compound occurrence table or a
#' `uaf_plant_phytochemistry` result.
#' @param categorate_result Optional categorate-like enrichment result.
#' @param compound_resolution Optional compound resolution table.
#' @param plant_queries Optional plant query table.
#' @param provider_diagnostics Optional provider diagnostics table.
#' @param comparability Optional precomputed `ChemistryComparability` table.
#'
#' @return Species-level summary table.
#'
#' @export
summarizePlantPhytochemistry = function(plant_compounds,
                                        categorate_result = NULL,
                                        compound_resolution = NULL,
                                        plant_queries = NULL,
                                        provider_diagnostics = NULL,
                                        comparability = NULL) {
  if (inherits(plant_compounds, "uaf_plant_phytochemistry")) {
    x = plant_compounds
    plant_compounds = x$PlantCompoundOccurrences
    categorate_result = x$CategorateResult
    compound_resolution = x$CompoundResolution
    plant_queries = x$PlantQueries
    provider_diagnostics = x$ProviderDiagnostics
    comparability = x$ChemistryComparability
  }
  occurrences = .plant_normalize_occurrences(plant_compounds)
  if (is.null(plant_queries)) plant_queries = .plant_queries(unique(occurrences$species))
  if (is.null(compound_resolution)) {
    compound_resolution = .plant_compound_resolution(occurrences,
                                                     categorate_result)
  }
  if (is.null(comparability)) {
    comparability = plantChemistryComparability(
      list(PlantCompoundOccurrences = occurrences,
           CategorateResult = categorate_result),
      min_confidence = "low"
    )
  }

  base_species = unique(.uaf_non_empty(c(plant_queries$species,
                                         occurrences$species)))
  if (length(base_species) < 1) return(.uaf_empty_table(.plant_summary_cols()))
  occurrence_index = split(seq_len(nrow(occurrences)), occurrences$species)
  query_index = split(seq_len(nrow(plant_queries)), plant_queries$species)
  comparability_index = if (is.data.frame(comparability) &&
                            nrow(comparability) > 0L) {
    split(seq_len(nrow(comparability)), comparability$species)
  } else {
    list()
  }
  resolution_key = .plant_clean_compound(
    compound_resolution$compound_name_clean
  )
  resolution_unique = !anyDuplicated(resolution_key)
  resolution_index = if (resolution_unique) {
    NULL
  } else {
    split(seq_len(nrow(compound_resolution)), resolution_key)
  }
  rows = lapply(base_species, function(species) {
    group_idx = occurrence_index[[species]]
    if (is.null(group_idx)) group_idx = integer()
    group = occurrences[group_idx, , drop = FALSE]
    query_idx = query_index[[species]]
    if (is.null(query_idx)) query_idx = integer()
    query = plant_queries[query_idx, , drop = FALSE]
    group_compounds = unique(.uaf_non_empty(group$compound_name_clean))
    if (resolution_unique) {
      resolved_idx = match(group_compounds, resolution_key, nomatch = 0L)
      resolved_idx = resolved_idx[resolved_idx > 0L]
    } else {
      resolved_idx = unlist(resolution_index[group_compounds],
                            use.names = FALSE)
      resolved_idx = resolved_idx[!is.na(resolved_idx)]
    }
    resolved = compound_resolution[resolved_idx, , drop = FALSE]
    source_count = length(unique(.uaf_non_empty(group$source_database)))
    direct_species = group[group$matched_rank == "species", , drop = FALSE]
    genus_fallback = group[group$matched_rank == "genus", , drop = FALSE]
    family_fallback = group[group$matched_rank == "family", , drop = FALSE]
    literature = group[grepl("literature|pubtator", group$evidence_tier,
                             ignore.case = TRUE), , drop = FALSE]
    database_rows = group[!grepl("literature|pubtator", group$evidence_tier,
                                 ignore.case = TRUE), , drop = FALSE]
    ready = group[group$analysis_ready == "Yes", , drop = FALSE]
    quality = suppressWarnings(as.numeric(group$evidence_quality_score))
    context_known = group[group$biological_context_status !=
                            "context_missing", , drop = FALSE]
    trait_info = .plant_summary_trait_info(group, categorate_result)
    group_comparability_idx = comparability_index[[species]]
    if (is.null(group_comparability_idx)) {
      group_comparability_idx = integer()
    }
    group_comparability = comparability[group_comparability_idx, ,
                                        drop = FALSE]
    comparability_info = .plant_summary_comparability_info(group,
                                                           categorate_result,
                                                           group_comparability)
    data.frame(
      species = species,
      species_slug = .plant_slug(species),
      genus = .uaf_first_non_empty_text(query$genus, .plant_genus(species)),
      family = .uaf_first_non_empty_text(query$family, group$family),
      query_status = ifelse(nrow(group) > 0, "records_found", "no_public_records"),
      compound_count = length(unique(.uaf_non_empty(group$compound_name_clean))),
      resolved_compound_count = sum(resolved$resolved %in% TRUE),
      unresolved_compound_count = sum(!(resolved$resolved %in% TRUE)),
      analysis_ready_compound_count = length(unique(.uaf_non_empty(ready$compound_name_clean))),
      direct_species_compound_count = length(unique(.uaf_non_empty(direct_species$compound_name_clean))),
      genus_level_compound_count = length(unique(.uaf_non_empty(genus_fallback$compound_name_clean))),
      family_level_compound_count = length(unique(.uaf_non_empty(family_fallback$compound_name_clean))),
      literature_candidate_count = nrow(literature),
      database_occurrence_count = nrow(database_rows),
      direct_reported_occurrence_count = sum(group$occurrence_status ==
                                               "direct_reported"),
      curated_reported_occurrence_count = sum(group$occurrence_status ==
                                                "curated_reported"),
      candidate_occurrence_count = sum(group$occurrence_status == "candidate"),
      taxon_fallback_occurrence_count = sum(group$occurrence_status ==
                                              "taxon_fallback"),
      context_known_occurrence_count = nrow(context_known),
      mean_evidence_quality_score = ifelse(length(quality[is.finite(quality)]) > 0,
                                           round(mean(quality[is.finite(quality)]), 3),
                                           NA_real_),
      analysis_ready_fraction = ifelse(nrow(group) > 0,
                                       round(nrow(ready) / nrow(group), 3), 0),
      plant_part_groups = .pubchem_collapse(group$plant_part_group[
        !group$plant_part_group %in% c("unknown", "extract_unspecified")
      ]),
      tissue_groups = .pubchem_collapse(group$tissue_group[
        !group$tissue_group %in% c("unknown", "extract_unspecified")
      ]),
      method_groups = .pubchem_collapse(group$method_group[
        !group$method_group %in% c("unknown")
      ]),
      evidence_tier_summary = .plant_count_summary(group$evidence_tier),
      source_database_count = source_count,
      source_databases = .pubchem_collapse(group$source_database),
      natural_product_superclasses = trait_info$superclasses,
      natural_product_classes = trait_info$classes,
      natural_product_subclasses = trait_info$subclasses,
      dominant_compound_classes = trait_info$dominant_classes,
      kingdoms_observed = trait_info$kingdoms,
      families_observed = .pubchem_collapse(unique(.uaf_non_empty(group$family))),
      is_plant_occurring = any(group$occurrence_status %in%
                                 c("direct_reported", "curated_reported")),
      volatile_proxy_fraction = trait_info$volatile_fraction,
      lipophilic_fraction = trait_info$lipophilic_fraction,
      oxygenated_fraction = trait_info$oxygenated_fraction,
      nitrogenous_fraction = trait_info$nitrogenous_fraction,
      sulfur_containing_fraction = trait_info$sulfur_fraction,
      halogenated_fraction = trait_info$halogen_fraction,
      kegg_pathway_groups = trait_info$kegg_pathway_groups,
      chemical_trait_count = trait_info$chemical_trait_count,
      high_confidence_trait_count = trait_info$high_confidence_trait_count,
      comparable_specialized_compound_count =
        comparability_info$specialized_count,
      comparable_volatile_compound_count = comparability_info$volatile_count,
      comparable_primary_compound_count = comparability_info$primary_count,
      comparable_lipid_compound_count = comparability_info$lipid_count,
      comparable_unknown_compound_count = comparability_info$unknown_count,
      comparison_scopes = comparability_info$comparison_scopes,
      dominant_comparison_groups = comparability_info$dominant_groups,
      source_coverage_score = .plant_source_coverage_score(source_count,
                                                           provider_diagnostics),
      uafR_validation_status = .plant_categorate_validation_status(categorate_result),
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.plant_automatic_matrix_max_traits = function() 5000L

#' Build plant phytochemistry matrices
#'
#' @param x Plant phytochemistry result or normalized occurrence table.
#' @param level One of `"species"`, `"genus"`, or `"family"`.
#' @param profile One of `"core"`, `"full"`, `"ecology"`, `"metabolism"`,
#' `"safety"`, or `"bioactivity"`.
#' @param mode One of `"binary"`, `"count"`, or `"confidence"`.
#' @param min_confidence Minimum confidence to include.
#' @param max_traits Maximum number of columns to retain.
#'
#' @return Wide matrix-like data frame.
#'
#' @export
plantPhytochemistryMatrix = function(x,
                                     level = c("species", "genus", "family"),
                                     profile = c("core", "full", "ecology",
                                                 "metabolism", "safety",
                                                 "bioactivity"),
                                     mode = c("binary", "count", "confidence"),
                                     min_confidence = "medium",
                                     max_traits = Inf) {
  level = match.arg(level)
  profile = match.arg(profile)
  mode = match.arg(mode)
  occurrences = if (inherits(x, "uaf_plant_phytochemistry")) {
    x$PlantCompoundOccurrences
  } else if (is.list(x) && is.data.frame(x$PlantCompoundOccurrences)) {
    x$PlantCompoundOccurrences
  } else {
    x
  }
  categorate_result = if (inherits(x, "uaf_plant_phytochemistry")) {
    x$CategorateResult
  } else if (is.list(x)) {
    x$CategorateResult
  } else {
    NULL
  }
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) {
    return(.uaf_empty_table(c(level)))
  }
  keep = .plant_confidence_score(occurrences$confidence) >=
    .plant_confidence_score(min_confidence)
  occurrences = occurrences[keep, , drop = FALSE]
  if (nrow(occurrences) < 1) {
    return(.uaf_empty_table(c(level)))
  }
  id = occurrences[[level]]
  id[is.na(id) | id == ""] = occurrences$species[is.na(id) | id == ""]
  long = data.frame(
    id = id,
    trait = paste0("compound__", .plant_matrix_key(occurrences$compound_name_clean)),
    value = .plant_confidence_score(occurrences$confidence),
    stringsAsFactors = FALSE
  )
  if (profile %in% c("core", "full", "ecology")) {
    long = rbind(long, data.frame(
      id = id,
      trait = paste0("source__", .plant_matrix_key(occurrences$source_database)),
      value = .plant_confidence_score(occurrences$confidence),
      stringsAsFactors = FALSE
    ))
    long = rbind(long, data.frame(
      id = id,
      trait = paste0("evidence__", .plant_matrix_key(occurrences$evidence_tier)),
      value = .plant_confidence_score(occurrences$confidence),
      stringsAsFactors = FALSE
    ))
    long = rbind(long, data.frame(
      id = id,
      trait = paste0("occurrence_status__",
                     .plant_matrix_key(occurrences$occurrence_status)),
      value = suppressWarnings(as.numeric(occurrences$evidence_quality_score)),
      stringsAsFactors = FALSE
    ))
    long = rbind(long, .plant_context_matrix_long(
      occurrences, level, "plant_part_group", "plant_part",
      exclude = c("unknown", "extract_unspecified")
    ))
    long = rbind(long, .plant_context_matrix_long(
      occurrences, level, "tissue_group", "tissue",
      exclude = c("unknown", "extract_unspecified")
    ))
    long = rbind(long, .plant_context_matrix_long(
      occurrences, level, "method_group", "method",
      exclude = c("unknown")
    ))
  }
  trait_long = .plant_chemical_trait_matrix_long(
    occurrences = occurrences,
    categorate_result = categorate_result,
    level = level,
    profile = profile,
    min_confidence = min_confidence
  )
  if (nrow(trait_long) > 0) {
    long = rbind(long, trait_long)
  }
  long = long[!is.na(long$id) & long$id != "" &
                !is.na(long$trait) & long$trait != "compound__", ,
              drop = FALSE]
  .plant_wide_matrix(long, id_col = level, mode = mode,
                     max_traits = max_traits)
}

#' Classify plant chemistry into comparable analysis scopes
#'
#' @description
#' `plantChemistryComparability()` adds a conservative comparison layer to
#' plant-compound evidence. The output separates biological role
#' (`metabolism_domain`), biosynthetic or chemical family
#' (`biosynthetic_family`), and analytical/physicochemical behavior
#' (`chemical_behavior`). This prevents downstream analyses from mixing
#' unmatched concepts such as primary metabolites, specialized metabolites,
#' lipids, hormone signals, and volatile fractions.
#'
#' The classifier is deterministic and source-transparent. It prioritizes
#' source-backed uafR class fields from `ChemicalClasses`, `LOTUSProfile`,
#' `PubChemClassifications`, `ChemicalTraitOntology`, `ChemicalTerms`,
#' `KEGGClassifications`, `KEGGPathways`, and `DerivedGroups`, then falls back
#' to `ChemicalTraits`/`ChemicalClasses` aggregates, compound-name patterns,
#' and occurrence evidence text. Unknown compounds remain `unknown` and are
#' excluded from comparable matrices unless explicitly requested.
#'
#' @param x Plant phytochemistry result, normalized occurrence table, or named
#' list containing `PlantCompoundOccurrences` and optional `CategorateResult`.
#' @param categorate_result Optional categorate-like enrichment result when
#' `x` is an occurrence table.
#' @param min_confidence Minimum occurrence confidence to include.
#'
#' @return A `ChemistryComparability` data frame with one row per retained
#' plant-compound occurrence.
#'
#' @export
plantChemistryComparability = function(x,
                                       categorate_result = NULL,
                                       min_confidence = "low") {
  occurrences = if (inherits(x, "uaf_plant_phytochemistry")) {
    x$PlantCompoundOccurrences
  } else if (is.list(x) && is.data.frame(x$PlantCompoundOccurrences)) {
    x$PlantCompoundOccurrences
  } else {
    x
  }
  if (is.null(categorate_result)) {
    categorate_result = if (inherits(x, "uaf_plant_phytochemistry")) {
      x$CategorateResult
    } else if (is.list(x)) {
      x$CategorateResult
    } else {
      NULL
    }
  }
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) {
    return(.uaf_empty_table(.plant_comparability_cols()))
  }
  keep = .plant_confidence_score(occurrences$confidence) >=
    .plant_confidence_score(min_confidence)
  occurrences = occurrences[keep, , drop = FALSE]
  if (nrow(occurrences) < 1) {
    return(.uaf_empty_table(.plant_comparability_cols()))
  }
  traits = .plant_comparability_traits(categorate_result)
  classes = .plant_comparability_classes(categorate_result)
  class_sources = .plant_comparability_class_sources(categorate_result)
  if (nrow(traits) < 1L && nrow(classes) < 1L && nrow(class_sources) < 1L) {
    classification = .plant_comparability_unenriched_table(occurrences)
    out = data.frame(
      species = occurrences$species,
      species_slug = .plant_slug(occurrences$species),
      genus = occurrences$genus,
      family = occurrences$family,
      compound_name = occurrences$compound_name,
      compound_name_clean = occurrences$compound_name_clean,
      source_database = occurrences$source_database,
      source_record_id = occurrences$source_record_id,
      pmid = occurrences$pmid,
      evidence_url = occurrences$evidence_url,
      evidence_tier = occurrences$evidence_tier,
      confidence = occurrences$confidence,
      occurrence_status = occurrences$occurrence_status,
      analysis_ready = occurrences$analysis_ready,
      matched_rank = occurrences$matched_rank,
      plant_part_group = occurrences$plant_part_group,
      tissue_group = occurrences$tissue_group,
      method_group = occurrences$method_group,
      biological_context_status = occurrences$biological_context_status,
      evidence_quality_score = occurrences$evidence_quality_score,
      classification,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    out = out[, .plant_comparability_cols(), drop = FALSE]
    row.names(out) = NULL
    return(out)
  }
  rows = lapply(seq_len(nrow(occurrences)), function(i) {
    occ = occurrences[i, , drop = FALSE]
    classification = .plant_comparability_classify(occ, traits, classes,
                                                   class_sources)
    data.frame(
      species = occ$species,
      species_slug = .plant_slug(occ$species),
      genus = occ$genus,
      family = occ$family,
      compound_name = occ$compound_name,
      compound_name_clean = occ$compound_name_clean,
      source_database = occ$source_database,
      source_record_id = occ$source_record_id,
      pmid = occ$pmid,
      evidence_url = occ$evidence_url,
      evidence_tier = occ$evidence_tier,
      confidence = occ$confidence,
      occurrence_status = occ$occurrence_status,
      analysis_ready = occ$analysis_ready,
      matched_rank = occ$matched_rank,
      plant_part_group = occ$plant_part_group,
      tissue_group = occ$tissue_group,
      method_group = occ$method_group,
      biological_context_status = occ$biological_context_status,
      evidence_quality_score = occ$evidence_quality_score,
      metabolism_domain = classification$metabolism_domain,
      biosynthetic_family = classification$biosynthetic_family,
      chemical_behavior = classification$chemical_behavior,
      comparison_scope = classification$comparison_scope,
      comparison_group = classification$comparison_group,
      comparison_subgroup = classification$comparison_subgroup,
      comparability_confidence = classification$comparability_confidence,
      comparability_basis = classification$comparability_basis,
      classification_source = classification$classification_source,
      classification_source_table =
        classification$classification_source_table,
      classification_source_field =
        classification$classification_source_field,
      classification_source_value =
        classification$classification_source_value,
      classification_source_confidence =
        classification$classification_source_confidence,
      comparable_for_matrix = classification$comparable_for_matrix,
      comparison_caveat = classification$comparison_caveat,
      stringsAsFactors = FALSE
    )
  })
  out = .plant_bind_tables(rows, .plant_comparability_cols())
  row.names(out) = NULL
  out
}

#' Build matrices from comparable plant chemistry scopes
#'
#' @description
#' Builds wide matrices only after selecting a defensible chemistry comparison
#' scope. Use this for clustering, ordination, candidate scoring, and figures
#' when chemistry should be compared like-for-like. For example,
#' `"volatile_specialized_metabolites"` compares volatile specialized chemistry,
#' while `"primary_metabolites"` compares primary-metabolism compounds.
#'
#' @param x Plant phytochemistry result, `ChemistryComparability` table, or
#' named list containing `ChemistryComparability`.
#' @param level One of `"species"`, `"genus"`, or `"family"`.
#' @param comparison_scope Scope(s) to retain. Use `"all_classified"` to retain
#' every classified non-unknown scope.
#' @param feature Feature family to turn into columns: `"comparison_group"`,
#' `"biosynthetic_family"`, `"chemical_behavior"`, `"compound"`, or
#' `"comparison_scope"`.
#' @param mode One of `"binary"`, `"count"`, or `"confidence"`.
#' @param min_comparability_confidence Minimum classification confidence.
#' @param plant_part_group Optional plant-part groups to retain.
#' @param tissue_group Optional tissue groups to retain.
#' @param method_group Optional method groups to retain.
#' @param biological_context_status Optional context-status values to retain.
#' @param require_context Logical. If `TRUE`, retain only rows with known
#' plant part/tissue or analytical method context.
#' @param include_unknown Logical. If `TRUE`, unknown comparison groups can be
#' retained.
#' @param max_features Maximum number of feature columns to retain.
#'
#' @return Wide matrix-like data frame.
#'
#' @export
plantComparableChemistryMatrix = function(x,
                                          level = c("species", "genus",
                                                    "family"),
                                          comparison_scope =
                                            "specialized_metabolites",
                                          feature = c("comparison_group",
                                                      "biosynthetic_family",
                                                      "chemical_behavior",
                                                      "compound",
                                                      "comparison_scope"),
                                          mode = c("binary", "count",
                                                   "confidence"),
                                          min_comparability_confidence =
                                            "medium",
                                          plant_part_group = NULL,
                                          tissue_group = NULL,
                                          method_group = NULL,
                                          biological_context_status = NULL,
                                          require_context = FALSE,
                                          include_unknown = FALSE,
                                          max_features = Inf) {
  level = match.arg(level)
  feature = match.arg(feature)
  mode = match.arg(mode)
  comparability = if (inherits(x, "uaf_plant_phytochemistry") &&
                      is.data.frame(x$ChemistryComparability)) {
    x$ChemistryComparability
  } else if (is.list(x) && is.data.frame(x$ChemistryComparability)) {
    x$ChemistryComparability
  } else if (is.data.frame(x) &&
             all(c("comparison_scope", "comparison_group") %in% names(x))) {
    x
  } else {
    plantChemistryComparability(x, min_confidence = "low")
  }
  comparability = .plant_normalize_comparability(comparability)
  if (nrow(comparability) < 1) return(.uaf_empty_table(c(level)))
  keep = .plant_confidence_score(comparability$comparability_confidence) >=
    .plant_confidence_score(min_comparability_confidence)
  scope = .uaf_non_empty(comparison_scope)
  if (length(scope) > 0 && !"all" %in% tolower(scope)) {
    if ("all_classified" %in% tolower(scope)) {
      keep = keep & comparability$comparison_scope != "unknown"
    } else {
      keep = keep & comparability$comparison_scope %in% scope
    }
  }
  if (!isTRUE(include_unknown)) {
    keep = keep & comparability$comparison_scope != "unknown" &
      comparability$comparison_group != "unknown" &
      comparability$comparable_for_matrix == "Yes"
  }
  keep = .plant_apply_optional_filter(keep, comparability$plant_part_group,
                                      plant_part_group)
  keep = .plant_apply_optional_filter(keep, comparability$tissue_group,
                                      tissue_group)
  keep = .plant_apply_optional_filter(keep, comparability$method_group,
                                      method_group)
  if ("biological_context_status" %in% names(comparability)) {
    keep = .plant_apply_optional_filter(
      keep, comparability$biological_context_status,
      biological_context_status
    )
  }
  if (isTRUE(require_context)) {
    has_part = !comparability$plant_part_group %in%
      c("unknown", "extract_unspecified")
    has_tissue = !comparability$tissue_group %in%
      c("unknown", "extract_unspecified")
    has_method = !comparability$method_group %in%
      c("unknown", "database_record", "literature_curation")
    keep = keep & (has_part | has_tissue | has_method)
  }
  comparability = comparability[keep, , drop = FALSE]
  if (nrow(comparability) < 1) return(.uaf_empty_table(c(level)))

  id = comparability[[level]]
  id[is.na(id) | id == ""] = comparability$species[is.na(id) | id == ""]
  feature_values = switch(
    feature,
    comparison_group = comparability$comparison_group,
    biosynthetic_family = comparability$biosynthetic_family,
    chemical_behavior = comparability$chemical_behavior,
    compound = comparability$compound_name_clean,
    comparison_scope = comparability$comparison_scope
  )
  prefix = switch(
    feature,
    comparison_group = "comparison_group",
    biosynthetic_family = "biosynthetic_family",
    chemical_behavior = "chemical_behavior",
    compound = "compound",
    comparison_scope = "comparison_scope"
  )
  values = pmin(
    .plant_confidence_score(comparability$confidence),
    .plant_confidence_score(comparability$comparability_confidence),
    suppressWarnings(as.numeric(comparability$evidence_quality_score)),
    na.rm = TRUE
  )
  values[!is.finite(values)] = pmin(
    .plant_confidence_score(comparability$confidence),
    .plant_confidence_score(comparability$comparability_confidence),
    na.rm = TRUE
  )
  long_keyed = data.frame(
    id = id,
    trait = paste0(prefix, "__",
                   .plant_matrix_key(comparability$comparison_scope),
                   "__", .plant_matrix_key(feature_values)),
    compound_name_clean = comparability$compound_name_clean,
    value = values,
    stringsAsFactors = FALSE
  )
  long_keyed = long_keyed[!is.na(long_keyed$id) & long_keyed$id != "" &
                            !is.na(long_keyed$trait) &
                            long_keyed$trait != "", ,
                          drop = FALSE]
  long_keyed = unique(long_keyed)
  long = long_keyed[, c("id", "trait", "value"), drop = FALSE]
  .plant_wide_matrix(long, id_col = level, mode = mode,
                     max_traits = max_features)
}

#' Extract biological context evidence from plant occurrence rows
#'
#' @description
#' `plantContextEvidence()` extracts auditable plant-part, tissue, and method
#' context evidence from normalized plant-compound occurrence rows. It preserves
#' the raw matched text, normalized context group, extraction rule, source field,
#' source record identifiers, confidence, and review flag. The result is used by
#' plant phytochemistry workflows to separate context-known rows from
#' context-missing rows without fabricating biological location.
#'
#' @param x Plant phytochemistry result, normalized occurrence table, or named
#' list containing `PlantCompoundOccurrences`.
#' @param min_confidence Minimum context confidence to retain.
#'
#' @return A `PlantContextEvidence` data frame.
#'
#' @export
plantContextEvidence = function(x, min_confidence = "low") {
  occurrences = if (inherits(x, "uaf_plant_phytochemistry")) {
    x$PlantCompoundOccurrences
  } else if (is.list(x) && is.data.frame(x$PlantCompoundOccurrences)) {
    x$PlantCompoundOccurrences
  } else {
    x
  }
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) {
    return(.uaf_empty_table(.plant_context_evidence_cols()))
  }
  signal_rows = .plant_context_signal_rows(occurrences)
  if (length(signal_rows) < 1L) {
    return(.uaf_empty_table(.plant_context_evidence_cols()))
  }
  occurrences = occurrences[signal_rows, , drop = FALSE]
  rows = list(
    .plant_direct_context_table(occurrences),
    .plant_text_context_table(occurrences)
  )
  out = .plant_normalize_context_evidence(
    .plant_bind_tables(rows, .plant_context_evidence_cols())
  )
  if (nrow(out) < 1) return(out)
  keep = .plant_confidence_score(out$context_confidence) >=
    .plant_confidence_score(min_confidence)
  out = out[keep, , drop = FALSE]
  out = unique(out)
  row.names(out) = NULL
  out
}

#' Enrich plant context evidence from source-backed literature text
#'
#' @description
#' `enrichPlantContextEvidence()` extracts plant-part, tissue, and method
#' context from source text linked to occurrence rows by PMID or DOI. The
#' function is intentionally conservative: source text must contain a context
#' term and must mention the occurrence species/genus, the compound, or a
#' chemical-profile source phrase before context is emitted. Rows preserve the
#' source PMID/DOI, source field, matched text, extraction rule, confidence, and
#' review flag.
#'
#' When `x` is a plant phytochemistry result, the returned result has updated
#' `PlantContextEvidence`, `PlantCompoundOccurrences`, and
#' `ProviderContextAudit` tables. When `x` is an occurrence table and
#' `apply = FALSE`, the function returns only the source-backed
#' `PlantContextEvidence` rows.
#'
#' @param x Plant phytochemistry result, normalized occurrence table, or named
#' list containing `PlantCompoundOccurrences`.
#' @param context_sources Optional source text table. Columns are matched to the
#' `LiteratureCandidates` schema; useful fields include `pmid`, `doi`, `title`,
#' `abstract`, `evidence_text`, and `evidence_url`.
#' @param fetch_pubmed Logical. If `TRUE`, fetch PubMed abstracts for occurrence
#' PMIDs and DOI-to-PMID matches. Network access is never attempted when this is
#' `FALSE`.
#' @param cache Logical. If `TRUE`, cache PubMed requests.
#' @param cache_dir Cache directory for PubMed context requests.
#' @param throttle Seconds to wait between uncached PubMed requests.
#' @param ncbi_email Optional NCBI email.
#' @param ncbi_tool Optional NCBI tool name.
#' @param ncbi_api_key Optional NCBI API key. It is not written to outputs.
#' @param max_sources Maximum unique PubMed source records to fetch. Candidate
#' PMID/DOI records are prioritized by available source text, direct occurrence
#' evidence, and shared species/compound coverage before fetching.
#' @param request_fun Optional request function for tests or controlled HTTP.
#' @param request_timeout Maximum seconds allowed for uncached PubMed requests.
#' @param min_confidence Minimum context confidence to retain.
#' @param apply Logical. If `TRUE`, apply the source-backed context evidence to
#' occurrence rows when returning a result/list.
#'
#' @return A source-backed `PlantContextEvidence` table, or an updated plant
#' phytochemistry result/list when `x` is a result/list and `apply = TRUE`.
#'
#' @export
enrichPlantContextEvidence = function(x,
                                      context_sources = NULL,
                                      fetch_pubmed = FALSE,
                                      cache = TRUE,
                                      cache_dir = NULL,
                                      throttle = 0.34,
                                      ncbi_email = Sys.getenv("NCBI_EMAIL", ""),
                                      ncbi_tool = Sys.getenv("NCBI_TOOL",
                                                             "uafR"),
                                      ncbi_api_key = Sys.getenv("NCBI_API_KEY",
                                                               ""),
                                      max_sources = 100,
                                      request_fun = NULL,
                                      request_timeout = 30,
                                      min_confidence = "low",
                                      apply = inherits(x, "uaf_plant_phytochemistry") ||
                                        (is.list(x) &&
                                           is.data.frame(x$PlantCompoundOccurrences))) {
  is_result = inherits(x, "uaf_plant_phytochemistry") ||
    (is.list(x) && is.data.frame(x$PlantCompoundOccurrences))
  occurrences = if (is_result) x$PlantCompoundOccurrences else x
  occurrences = .plant_normalize_occurrences(occurrences)
  existing_context = if (is_result && is.data.frame(x$PlantContextEvidence)) {
    .plant_normalize_context_evidence(x$PlantContextEvidence)
  } else {
    .uaf_empty_table(.plant_context_evidence_cols())
  }
  if (nrow(occurrences) < 1) {
    return(if (is_result && isTRUE(apply)) x else
      .uaf_empty_table(.plant_context_evidence_cols()))
  }
  cache_dir = .plant_cache_dir(cache_dir)
  sources = .plant_normalize_context_sources(context_sources)
  if (isTRUE(fetch_pubmed)) {
    fetched = .plant_fetch_pubmed_context_sources(
      occurrences = occurrences,
      context_sources = sources,
      cache = cache,
      cache_dir = file.path(cache_dir, "pubmed_context"),
      throttle = throttle,
      ncbi_email = ncbi_email,
      ncbi_tool = ncbi_tool,
      ncbi_api_key = ncbi_api_key,
      max_sources = max_sources,
      request_fun = request_fun,
      request_timeout = request_timeout
    )
    sources = .plant_bind_tables(list(sources, fetched),
                                 .plant_literature_cols())
  }
  source_context = .plant_source_context_evidence(
    occurrences = occurrences,
    context_sources = sources,
    min_confidence = min_confidence
  )
  if (!is_result || !isTRUE(apply)) return(source_context)

  combined_context = .plant_normalize_context_evidence(
    .plant_bind_tables(list(existing_context, source_context),
                       .plant_context_evidence_cols())
  )
  x$PlantContextEvidence = combined_context
  if (nrow(combined_context) > 0) {
    occurrences = .plant_apply_context_evidence(occurrences,
                                                combined_context)
  }
  occurrences = .plant_clean_context_conflicts(occurrences)
  x$PlantCompoundOccurrences = occurrences
  x$ProviderContextAudit = plantProviderContextAudit(
    list(PlantCompoundOccurrences = occurrences,
         PlantContextEvidence = combined_context)
  )
  x
}

.plant_has_context_signal = function(occurrences) {
  length(.plant_context_signal_rows(occurrences)) > 0L
}

.plant_context_signal_rows = function(occurrences) {
  if (!is.data.frame(occurrences) || nrow(occurrences) < 1L) {
    return(integer())
  }
  direct = rep(FALSE, nrow(occurrences))
  for (field in c("plant_part", "tissue", "method")) {
    value = .uaf_squish_text(
      .plant_col_or_default(occurrences, field, NA_character_)
    )
    direct = direct | (!is.na(value) & value != "")
  }
  text = paste(
    .plant_col_or_default(occurrences, "evidence_text", NA_character_),
    .plant_col_or_default(occurrences, "occurrence_type", NA_character_),
    sep = " "
  )
  pattern = paste(.plant_context_patterns()$pattern, collapse = "|")
  text_signal = grepl(pattern, text, ignore.case = TRUE, perl = TRUE)
  text_signal[is.na(text_signal)] = FALSE
  which(direct | text_signal)
}

.plant_normalize_context_sources = function(context_sources) {
  if (is.null(context_sources)) return(.uaf_empty_table(.plant_literature_cols()))
  out = .plant_normalize_literature(context_sources,
                                    source_hint = "source_context")
  out$doi = .plant_normalize_doi(out$doi)
  out$pmid = .plant_normalize_pmid(out$pmid)
  text = paste(out$title, out$abstract, out$evidence_text, sep = " ")
  has_text = !is.na(.uaf_squish_text(text)) & .uaf_squish_text(text) != ""
  out = out[has_text |
              !is.na(out$pmid) | !is.na(out$doi), , drop = FALSE]
  out = unique(out)
  row.names(out) = NULL
  out
}

.plant_source_context_evidence = function(occurrences, context_sources,
                                          min_confidence = "low") {
  occurrences = .plant_normalize_occurrences(occurrences)
  sources = .plant_normalize_context_sources(context_sources)
  if (nrow(occurrences) < 1 || nrow(sources) < 1) {
    return(.uaf_empty_table(.plant_context_evidence_cols()))
  }
  source_signatures = .plant_context_source_signatures(sources)
  source_identity = paste(
    .plant_normalize_pmid(sources$pmid),
    .plant_normalize_doi(sources$doi),
    ifelse(is.na(sources$evidence_url), "", sources$evidence_url),
    source_signatures,
    sep = "\r"
  )
  sources = sources[!duplicated(source_identity), , drop = FALSE]
  source_signatures = source_signatures[!duplicated(source_identity)]
  occurrence_keys = .plant_source_context_keys(occurrences)
  source_keys = .plant_source_context_keys(sources)
  source_index = .plant_source_context_index(source_keys)
  patterns = .plant_context_patterns()
  unique_source_index = which(!duplicated(source_signatures))
  prepared_templates = lapply(unique_source_index, function(i) {
    .plant_prepare_context_source(sources[i, , drop = FALSE], patterns)
  })
  template_index = match(source_signatures,
                         source_signatures[unique_source_index])
  prepared_sources = lapply(seq_len(nrow(sources)), function(i) {
    prepared = prepared_templates[[template_index[[i]]]]
    prepared$source = sources[i, , drop = FALSE]
    prepared
  })
  rows = list()
  row_index = 0L
  for (i in seq_len(nrow(occurrences))) {
    keys = occurrence_keys[[i]]
    if (length(keys) < 1) next
    hit_idx = unique(unlist(source_index[keys], use.names = FALSE))
    hit_idx = hit_idx[!is.na(hit_idx)]
    if (length(hit_idx) < 1) next
    occurrence = occurrences[i, , drop = FALSE]
    occurrence_terms = .plant_context_occurrence_terms(occurrence)
    for (j in hit_idx) {
      evidence = .plant_context_rows_from_prepared_source(
        occurrence = occurrence,
        prepared = prepared_sources[[j]],
        occurrence_terms = occurrence_terms
      )
      if (nrow(evidence) < 1L) next
      row_index = row_index + 1L
      rows[[row_index]] = evidence
    }
  }
  out = .plant_normalize_context_evidence(
    .plant_bind_tables(rows, .plant_context_evidence_cols())
  )
  if (nrow(out) < 1) return(out)
  keep = .plant_confidence_score(out$context_confidence) >=
    .plant_confidence_score(min_confidence)
  out = unique(out[keep, , drop = FALSE])
  row.names(out) = NULL
  out
}

.plant_context_source_signatures = function(sources) {
  provider = vapply(sources$source_database, function(value) {
    .plant_provider_key(.uaf_first_non_empty_text(value, "unknown"))
  }, character(1))
  fields = lapply(c("title", "abstract", "evidence_text"), function(col) {
    value = .plant_col_or_default(sources, col, NA_character_)
    value[is.na(value)] = ""
    value
  })
  do.call(paste, c(list(provider), fields, sep = "\r"))
}

.plant_source_context_keys = function(x) {
  if (!is.data.frame(x) || nrow(x) < 1) return(list())
  pmid = .plant_normalize_pmid(.plant_col_or_default(x, "pmid",
                                                     NA_character_))
  doi = .plant_normalize_doi(.plant_col_or_default(x, "doi", NA_character_))
  lapply(seq_len(nrow(x)), function(i) {
    .uaf_non_empty(c(if (!is.na(pmid[[i]])) paste0("pmid:", pmid[[i]]),
                     if (!is.na(doi[[i]])) paste0("doi:", doi[[i]])))
  })
}

.plant_source_context_index = function(source_keys) {
  if (length(source_keys) < 1) return(list())
  rows = list()
  for (i in seq_along(source_keys)) {
    keys = source_keys[[i]]
    if (length(keys) < 1) next
    for (key in keys) {
      rows[[length(rows) + 1]] = data.frame(
        key = key,
        source_index = i,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) < 1) return(list())
  index = do.call(rbind, rows)
  split(index$source_index, index$key)
}

.plant_context_rows_from_source = function(occurrence, source) {
  prepared = .plant_prepare_context_source(
    source = source,
    patterns = .plant_context_patterns()
  )
  .plant_context_rows_from_prepared_source(
    occurrence = occurrence,
    prepared = prepared,
    occurrence_terms = .plant_context_occurrence_terms(occurrence)
  )
}

.plant_prepare_context_source = function(source, patterns) {
  source = .plant_normalize_literature(source)
  if (nrow(source) < 1) {
    return(list(
      source = .uaf_empty_table(.plant_literature_cols()),
      document_lower = "",
      sentences = data.frame(),
      hits = data.frame()
    ))
  }
  source = source[1, , drop = FALSE]
  text_fields = c(
    title = .uaf_first_non_empty_text(source$title),
    abstract = .uaf_first_non_empty_text(source$abstract),
    evidence_text = .uaf_first_non_empty_text(source$evidence_text)
  )
  document_text = paste(.uaf_non_empty(text_fields), collapse = " ")
  sentence_rows = list()
  hit_rows = list()
  seen_text = character()
  sentence_index = 0L
  provider = .plant_provider_key(source$source_database)
  for (field in names(text_fields)) {
    text = text_fields[[field]]
    if (is.na(text) || text == "") next
    text_key = tolower(.uaf_squish_text(text))
    if (!is.na(text_key) && text_key %in% seen_text) next
    seen_text = c(seen_text, text_key)
    sentences = .plant_split_sentences(text)
    if (length(sentences) < 1) sentences = text
    for (sentence in sentences) {
      raw_hits = vapply(patterns$pattern, function(pattern) {
        .plant_regex_match_text(sentence, pattern)
      }, character(1))
      keep_hits = !is.na(raw_hits) & raw_hits != "" &
        patterns$normalized_context != "extract_unspecified"
      if (!any(keep_hits)) next
      sentence_index = sentence_index + 1L
      sentence_rows[[sentence_index]] = data.frame(
        sentence_index = sentence_index,
        source_field = field,
        text = sentence,
        text_lower = tolower(sentence),
        source_rule = paste0(provider, "_", field, "_sentence"),
        has_scope = grepl(.plant_source_context_scope_pattern(), sentence,
                          ignore.case = TRUE, perl = TRUE),
        stringsAsFactors = FALSE
      )
      idx = which(keep_hits)
      hit_rows[[sentence_index]] = data.frame(
        sentence_index = sentence_index,
        context_type = patterns$context_type[idx],
        raw_context_text = raw_hits[idx],
        normalized_context = patterns$normalized_context[idx],
        rule_name = patterns$rule_name[idx],
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    source = source,
    document_lower = tolower(document_text),
    document_has_scope = grepl(.plant_source_context_scope_pattern(),
                               document_text, ignore.case = TRUE,
                               perl = TRUE),
    sentences = if (length(sentence_rows) > 0) {
      do.call(rbind, sentence_rows)
    } else {
      data.frame()
    },
    hits = if (length(hit_rows) > 0) {
      do.call(rbind, hit_rows)
    } else {
      data.frame()
    }
  )
}

.plant_context_occurrence_terms = function(occurrence) {
  species = .uaf_first_non_empty_text(occurrence$species,
                                      occurrence$query_plant)
  genus = .uaf_first_non_empty_text(occurrence$genus,
                                    .plant_genus(species))
  compound = .uaf_first_non_empty_text(occurrence$compound_name)
  list(
    species = tolower(.uaf_first_non_empty_text(species, "")),
    genus = tolower(.uaf_first_non_empty_text(genus, "")),
    compound = tolower(.uaf_first_non_empty_text(compound, ""))
  )
}

.plant_context_fixed_has = function(text, term, min_chars = 3L,
                                     max_chars = Inf) {
  if (length(term) < 1 || is.na(term) || term == "" ||
      nchar(term) < min_chars || nchar(term) > max_chars) {
    return(rep(FALSE, length(text)))
  }
  grepl(term, text, fixed = TRUE)
}

.plant_context_rows_from_prepared_source = function(
    occurrence, prepared,
    occurrence_terms = .plant_context_occurrence_terms(occurrence)) {
  text_rows = prepared$sentences
  hits = prepared$hits
  if (!is.data.frame(text_rows) || nrow(text_rows) < 1 ||
      !is.data.frame(hits) || nrow(hits) < 1) {
    return(.uaf_empty_table(.plant_context_evidence_cols()))
  }
  has_species = .plant_context_fixed_has(
    text_rows$text_lower, occurrence_terms$species
  ) | .plant_context_fixed_has(text_rows$text_lower, occurrence_terms$genus)
  has_compound = .plant_context_fixed_has(
    text_rows$text_lower, occurrence_terms$compound,
    min_chars = 4L, max_chars = 120L
  )
  document_has_species = .plant_context_fixed_has(
    prepared$document_lower, occurrence_terms$species
  ) | .plant_context_fixed_has(prepared$document_lower,
                               occurrence_terms$genus)
  document_has_compound = .plant_context_fixed_has(
    prepared$document_lower, occurrence_terms$compound,
    min_chars = 4L, max_chars = 120L
  )
  keep = rep(FALSE, nrow(text_rows))
  confidence = rep("low", nrow(text_rows))
  basis = rep("source_text_context_not_linked_to_species_or_compound",
              nrow(text_rows))
  direct = has_species & has_compound
  keep[direct] = TRUE
  confidence[direct] = "medium"
  basis[direct] = "source_backed_species_compound_sentence"
  species_scope = !keep & has_species & text_rows$has_scope
  keep[species_scope] = TRUE
  confidence[species_scope] = "medium"
  basis[species_scope] =
    "source_backed_species_chemical_context_sentence"
  compound_scope = !keep & has_compound & text_rows$has_scope
  keep[compound_scope] = TRUE
  basis[compound_scope] = "source_backed_compound_context_sentence"
  document_scope = !keep & text_rows$has_scope &
    isTRUE(prepared$document_has_scope) &
    isTRUE(document_has_species) & isTRUE(document_has_compound)
  keep[document_scope] = TRUE
  basis[document_scope] = "source_backed_document_context_sentence"
  if (!any(keep)) {
    return(.uaf_empty_table(.plant_context_evidence_cols()))
  }
  sentence_metadata = data.frame(
    sentence_index = text_rows$sentence_index,
    source_field = text_rows$source_field,
    source_rule = text_rows$source_rule,
    context_confidence = confidence,
    evidence_basis = basis,
    keep = keep,
    stringsAsFactors = FALSE
  )
  hit_metadata = sentence_metadata[
    match(hits$sentence_index, sentence_metadata$sentence_index), ,
    drop = FALSE
  ]
  hits = cbind(hits, hit_metadata[, setdiff(names(hit_metadata),
                                            "sentence_index"), drop = FALSE])
  hits = hits[hits$keep, , drop = FALSE]
  if (nrow(hits) < 1) {
    return(.uaf_empty_table(.plant_context_evidence_cols()))
  }
  source_occurrence = .plant_occurrence_with_source(
    occurrence, prepared$source
  )
  out = .plant_context_evidence_rows_from_hits(source_occurrence, hits)
  sentence_groups = split(seq_len(nrow(out)), hits$sentence_index)
  rows = lapply(sentence_groups, function(idx) {
    .plant_bind_tables(
      .plant_prune_context_rows(lapply(idx, function(i) {
        out[i, , drop = FALSE]
      })),
      .plant_context_evidence_cols()
    )
  })
  .plant_bind_tables(rows, .plant_context_evidence_cols())
}

.plant_context_evidence_rows_from_hits = function(occurrence, hits) {
  n = nrow(hits)
  data.frame(
    species = rep(.uaf_first_non_empty_text(occurrence$species), n),
    species_slug = rep(.uaf_first_non_empty_text(
      occurrence$species_slug, .plant_slug(occurrence$species)
    ), n),
    genus = rep(.uaf_first_non_empty_text(occurrence$genus), n),
    family = rep(.uaf_first_non_empty_text(occurrence$family), n),
    compound_name = rep(.uaf_first_non_empty_text(
      occurrence$compound_name
    ), n),
    compound_name_clean = rep(.uaf_first_non_empty_text(
      occurrence$compound_name_clean
    ), n),
    source_database = rep(.uaf_first_non_empty_text(
      occurrence$source_database
    ), n),
    source_record_id = rep(.uaf_first_non_empty_text(
      occurrence$source_record_id
    ), n),
    pmid = rep(.uaf_first_non_empty_text(occurrence$pmid), n),
    evidence_url = rep(.uaf_first_non_empty_text(
      occurrence$evidence_url
    ), n),
    context_type = hits$context_type,
    raw_context_text = hits$raw_context_text,
    normalized_context = hits$normalized_context,
    source_field = hits$source_field,
    extraction_rule = paste0("source_context:", hits$source_rule, ":",
                             hits$rule_name),
    context_confidence = hits$context_confidence,
    evidence_basis = hits$evidence_basis,
    requires_review = .uaf_yes_no(
      hits$context_confidence == "low" |
        hits$normalized_context %in% c("other", "extract_unspecified")
    ),
    retrieved_at = rep(.plant_timestamp(), n),
    stringsAsFactors = FALSE
  )
}

.plant_occurrence_with_source = function(occurrence, source) {
  occurrence$pmid = .uaf_first_non_empty_text(source$pmid, occurrence$pmid)
  occurrence$doi = .uaf_first_non_empty_text(source$doi, occurrence$doi)
  occurrence$evidence_url = .uaf_first_non_empty_text(
    source$evidence_url,
    if (!is.na(.uaf_first_non_empty_text(source$pmid))) {
      paste0("https://pubmed.ncbi.nlm.nih.gov/", source$pmid, "/")
    } else {
      NA_character_
    },
    occurrence$evidence_url
  )
  occurrence$evidence_text = .uaf_first_non_empty_text(
    source$evidence_text, source$title, occurrence$evidence_text
  )
  occurrence
}

.plant_source_context_text_rows = function(occurrence, source) {
  source = .plant_normalize_literature(source)
  text_fields = c(title = .uaf_first_non_empty_text(source$title),
                  abstract = .uaf_first_non_empty_text(source$abstract),
                  evidence_text = .uaf_first_non_empty_text(source$evidence_text))
  document_text = paste(.uaf_non_empty(text_fields), collapse = " ")
  rows = list()
  seen_text = character()
  for (field in names(text_fields)) {
    text = text_fields[[field]]
    if (is.na(text) || text == "") next
    text_key = tolower(.uaf_squish_text(text))
    if (!is.na(text_key) && text_key %in% seen_text) next
    seen_text = c(seen_text, text_key)
    sentences = .plant_split_sentences(text)
    if (length(sentences) < 1) sentences = text
    for (sentence in sentences) {
      gate = .plant_source_context_sentence_gate(occurrence, sentence,
                                                 document_text)
      if (!gate$keep) next
      rows[[length(rows) + 1]] = data.frame(
        source_field = field,
        text = sentence,
        source_rule = paste0(.plant_provider_key(source$source_database),
                             "_", field, "_sentence"),
        context_confidence = gate$confidence,
        evidence_basis = gate$basis,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) < 1) {
    return(.uaf_empty_table(c("source_field", "text", "source_rule",
                              "context_confidence", "evidence_basis")))
  }
  do.call(rbind, rows)
}

.plant_source_context_sentence_gate = function(occurrence, sentence,
                                               document_text = sentence) {
  sentence = .uaf_first_non_empty_text(sentence)
  if (is.na(sentence) || sentence == "") {
    return(list(keep = FALSE, confidence = "low",
                basis = "source_text_no_context"))
  }
  has_context = .plant_text_has_any_context_pattern(sentence)
  if (!has_context) {
    return(list(keep = FALSE, confidence = "low",
                basis = "source_text_no_context"))
  }
  species = .uaf_first_non_empty_text(occurrence$species,
                                      occurrence$query_plant)
  genus = .uaf_first_non_empty_text(occurrence$genus,
                                    .plant_genus(species))
  compound = .uaf_first_non_empty_text(occurrence$compound_name)
  has_species = .plant_text_has_term(sentence, species) ||
    .plant_text_has_term(sentence, genus)
  has_compound = .plant_text_has_compound(sentence, compound)
  has_scope = grepl(.plant_source_context_scope_pattern(), sentence,
                    ignore.case = TRUE, perl = TRUE)
  document_has_species = .plant_text_has_term(document_text, species) ||
    .plant_text_has_term(document_text, genus)
  document_has_compound = .plant_text_has_compound(document_text, compound)
  document_has_scope = grepl(.plant_source_context_scope_pattern(),
                             document_text, ignore.case = TRUE, perl = TRUE)
  if (has_species && has_compound) {
    return(list(keep = TRUE, confidence = "medium",
                basis = "source_backed_species_compound_sentence"))
  }
  if (has_species && has_scope) {
    return(list(keep = TRUE, confidence = "medium",
                basis = "source_backed_species_chemical_context_sentence"))
  }
  if (has_compound && has_scope) {
    return(list(keep = TRUE, confidence = "low",
                basis = "source_backed_compound_context_sentence"))
  }
  if (has_scope && document_has_scope && document_has_species &&
      document_has_compound) {
    return(list(keep = TRUE, confidence = "low",
                basis = "source_backed_document_context_sentence"))
  }
  list(keep = FALSE, confidence = "low",
       basis = "source_text_context_not_linked_to_species_or_compound")
}

.plant_text_has_any_context_pattern = function(text) {
  patterns = .plant_context_patterns()$pattern
  any(vapply(patterns, function(pattern) {
    grepl(pattern, text, ignore.case = TRUE, perl = TRUE)
  }, logical(1)), na.rm = TRUE)
}

.plant_source_context_scope_pattern = function() {
  paste(c("phytochemical", "chemical constituents?", "constituents?",
          "metabolites?", "secondary metabolites?", "natural products?",
          "essential oils?", "volatile oils?", "extracts?", "isolat",
          "identified", "profile", "composition", "gc[- ]?ms", "lc[- ]?ms",
          "hplc", "uplc", "nmr"), collapse = "|")
}

.plant_text_has_term = function(text, term) {
  term = .uaf_first_non_empty_text(term)
  if (is.na(term) || term == "" || nchar(term) < 3) return(FALSE)
  grepl(tolower(term), tolower(.uaf_first_non_empty_text(text, "")),
        fixed = TRUE)
}

.plant_text_has_compound = function(text, compound) {
  compound = .uaf_first_non_empty_text(compound)
  if (is.na(compound) || compound == "" || nchar(compound) < 4 ||
      nchar(compound) > 120) {
    return(FALSE)
  }
  .plant_text_has_term(text, compound)
}

.plant_fetch_pubmed_context_sources = function(occurrences, context_sources,
                                               cache, cache_dir, throttle,
                                               ncbi_email, ncbi_tool,
                                               ncbi_api_key, max_sources,
                                               request_fun,
                                               request_timeout) {
  occurrences = .plant_normalize_occurrences(occurrences)
  context_sources = .plant_normalize_context_sources(context_sources)
  max_sources = suppressWarnings(as.integer(max_sources[[1]]))
  if (!is.finite(max_sources) || max_sources < 1) max_sources = 100L
  candidates = .plant_pubmed_context_candidates(occurrences, context_sources)
  pmid_candidates = candidates[candidates$key_type == "pmid", , drop = FALSE]
  pmids = unique(.uaf_non_empty(pmid_candidates$pmid))
  pmids = utils::head(pmids, max_sources)
  if (nrow(candidates) > 0 && length(pmids) < max_sources) {
    doi_candidates = candidates[candidates$key_type == "doi", , drop = FALSE]
    for (doi in unique(.uaf_non_empty(doi_candidates$doi))) {
      if (length(pmids) >= max_sources) break
      doi_pmid = .uaf_first_non_empty_text(
        doi_candidates$pmid[doi_candidates$doi == doi]
      )
      if (!is.na(doi_pmid) && doi_pmid %in% pmids) next
      doi_pmids = .plant_pubmed_ids_for_doi(
        doi, cache, cache_dir, throttle, ncbi_email, ncbi_tool, ncbi_api_key,
        request_fun, request_timeout
      )
      pmids = unique(.uaf_non_empty(c(pmids, doi_pmids)))
    }
  }
  pmids = utils::head(pmids, max_sources)
  if (length(pmids) < 1) return(.uaf_empty_table(.plant_literature_cols()))
  batches = split(pmids, ceiling(seq_along(pmids) / 50))
  rows = lapply(batches, function(ids) {
    url = .plant_ncbi_url(
      endpoint = "efetch.fcgi",
      params = c(db = "pubmed", id = paste(ids, collapse = ","),
                 retmode = "xml", tool = ncbi_tool, email = ncbi_email),
      api_key = ncbi_api_key
    )
    txt = tryCatch(
      .plant_fetch_text(url, cache, cache_dir, throttle, request_fun,
                        timeout = request_timeout),
      error = function(error) NA_character_
    )
    .plant_pubmed_xml_sources(txt)
  })
  .plant_normalize_context_sources(
    .plant_bind_tables(rows, .plant_literature_cols())
  )
}

.plant_pubmed_context_candidates = function(occurrences, context_sources) {
  cols = c("key_type", "key", "pmid", "doi", "source_origin",
           "species", "compound_name", "source_database", "evidence_tier",
           "confidence", "source_text_signal", "evidence_rank",
           "first_seen")
  occurrences = .plant_normalize_occurrences(occurrences)
  context_sources = .plant_normalize_context_sources(context_sources)
  rows = list()
  if (nrow(occurrences) > 0) {
    rows = c(rows, .plant_pubmed_context_candidate_rows(
      table = occurrences,
      source_origin = "occurrence",
      species_col = "species",
      compound_col = "compound_name",
      text_cols = c("evidence_text", "occurrence_type"),
      first_seen_offset = 0L
    ))
  }
  if (nrow(context_sources) > 0) {
    rows = c(rows, .plant_pubmed_context_candidate_rows(
      table = context_sources,
      source_origin = "context_source",
      species_col = "species",
      compound_col = "chemical_mention",
      text_cols = c("title", "abstract", "evidence_text"),
      first_seen_offset = length(rows)
    ))
  }
  if (length(rows) < 1) return(.uaf_empty_table(cols))
  raw = .plant_bind_tables(rows, cols)
  raw = raw[!is.na(raw$key) & raw$key != "", , drop = FALSE]
  if (nrow(raw) < 1) return(.uaf_empty_table(cols))
  groups = split(seq_len(nrow(raw)), paste(raw$key_type, raw$key, sep = "::"))
  out = lapply(groups, function(idx) {
    group = raw[idx, , drop = FALSE]
    species = unique(.uaf_non_empty(.plant_clean_name(group$species)))
    compounds = unique(.uaf_non_empty(.plant_clean_compound(
      group$compound_name
    )))
    text_signal = max(suppressWarnings(as.integer(group$source_text_signal)),
                      na.rm = TRUE)
    evidence_rank = max(suppressWarnings(as.numeric(group$evidence_rank)),
                        na.rm = TRUE)
    if (!is.finite(text_signal)) text_signal = 0L
    if (!is.finite(evidence_rank)) evidence_rank = 0
    source_origin_rank = ifelse(any(group$source_origin == "context_source"),
                                3, 2)
    key_rank = ifelse(group$key_type[[1]] == "pmid", 2, 1)
    confidence_rank = max(.plant_confidence_score(group$confidence),
                          na.rm = TRUE)
    occurrence_count = nrow(group)
    species_count = length(species)
    compound_count = length(compounds)
    priority_score = (100 * text_signal) + (50 * source_origin_rank) +
      (25 * key_rank) + (10 * evidence_rank) +
      (2 * min(occurrence_count, 50)) + species_count + compound_count +
      confidence_rank
    data.frame(
      key_type = group$key_type[[1]],
      key = group$key[[1]],
      pmid = .uaf_first_non_empty_text(group$pmid),
      doi = .uaf_first_non_empty_text(group$doi),
      source_origin = .pubchem_collapse(unique(.uaf_non_empty(
        group$source_origin
      ))),
      species = .pubchem_collapse(species),
      compound_name = .pubchem_collapse(compounds),
      source_database = .pubchem_collapse(unique(.uaf_non_empty(
        group$source_database
      ))),
      evidence_tier = .pubchem_collapse(unique(.uaf_non_empty(
        group$evidence_tier
      ))),
      confidence = .uaf_first_non_empty_text(group$confidence, "unknown"),
      source_text_signal = as.character(text_signal),
      evidence_rank = as.character(evidence_rank),
      first_seen = as.character(min(suppressWarnings(as.integer(
        group$first_seen
      )), na.rm = TRUE)),
      priority_score = priority_score,
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, out)
  key_rank = ifelse(out$key_type == "pmid", 1L, 2L)
  out = out[order(-out$priority_score,
                  -suppressWarnings(as.numeric(out$source_text_signal)),
                  key_rank,
                  -suppressWarnings(as.numeric(out$evidence_rank)),
                  suppressWarnings(as.integer(out$first_seen)),
                  out$key), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_pubmed_context_candidate_rows = function(table, source_origin,
                                                species_col, compound_col,
                                                text_cols,
                                                first_seen_offset = 0L) {
  if (!is.data.frame(table) || nrow(table) < 1) return(list())
  pmids = .plant_normalize_pmid(.plant_col_or_default(table, "pmid",
                                                       NA_character_))
  dois = .plant_normalize_doi(.plant_col_or_default(table, "doi",
                                                    NA_character_))
  species = .plant_col_or_default(table, species_col, NA_character_)
  compounds = .plant_col_or_default(table, compound_col, NA_character_)
  missing_compound = is.na(compounds) | compounds == ""
  if (any(missing_compound)) {
    compounds[missing_compound] =
      .plant_col_or_default(table, "compound_name", NA_character_)[missing_compound]
  }
  text = .plant_pubmed_context_candidate_text(table, text_cols)
  text_signal = .plant_context_source_text_signal(text)
  evidence_rank = .plant_pubmed_context_evidence_rank(
    evidence_tier = .plant_col_or_default(table, "evidence_tier",
                                          NA_character_),
    source_database = .plant_col_or_default(table, "source_database",
                                            NA_character_),
    confidence = .plant_col_or_default(table, "confidence", NA_character_)
  )
  rows = list()
  for (i in seq_len(nrow(table))) {
    base = list(
      pmid = pmids[[i]],
      doi = dois[[i]],
      source_origin = source_origin,
      species = species[[i]],
      compound_name = compounds[[i]],
      source_database = .plant_col_or_default(table, "source_database",
                                              NA_character_)[[i]],
      evidence_tier = .plant_col_or_default(table, "evidence_tier",
                                            NA_character_)[[i]],
      confidence = .plant_col_or_default(table, "confidence",
                                         NA_character_)[[i]],
      source_text_signal = as.character(text_signal[[i]]),
      evidence_rank = as.character(evidence_rank[[i]]),
      first_seen = as.character(first_seen_offset + i)
    )
    if (!is.na(pmids[[i]]) && pmids[[i]] != "") {
      rows[[length(rows) + 1]] = data.frame(
        key_type = "pmid",
        key = pmids[[i]],
        pmid = base$pmid,
        doi = base$doi,
        source_origin = base$source_origin,
        species = base$species,
        compound_name = base$compound_name,
        source_database = base$source_database,
        evidence_tier = base$evidence_tier,
        confidence = base$confidence,
        source_text_signal = base$source_text_signal,
        evidence_rank = base$evidence_rank,
        first_seen = base$first_seen,
        stringsAsFactors = FALSE
      )
    }
    if (!is.na(dois[[i]]) && dois[[i]] != "") {
      rows[[length(rows) + 1]] = data.frame(
        key_type = "doi",
        key = dois[[i]],
        pmid = base$pmid,
        doi = base$doi,
        source_origin = base$source_origin,
        species = base$species,
        compound_name = base$compound_name,
        source_database = base$source_database,
        evidence_tier = base$evidence_tier,
        confidence = base$confidence,
        source_text_signal = base$source_text_signal,
        evidence_rank = base$evidence_rank,
        first_seen = base$first_seen,
        stringsAsFactors = FALSE
      )
    }
  }
  rows
}

.plant_pubmed_context_candidate_text = function(table, text_cols) {
  pieces = lapply(text_cols, function(col) {
    .plant_col_or_default(table, col, NA_character_)
  })
  if (length(pieces) < 1) return(rep(NA_character_, nrow(table)))
  text = do.call(paste, c(pieces, sep = " "))
  .uaf_squish_text(text)
}

.plant_context_source_text_signal = function(text) {
  text = .uaf_squish_text(text)
  context_patterns = .plant_context_patterns()$pattern
  scope_pattern = .plant_source_context_scope_pattern()
  vapply(text, function(item) {
    if (is.na(item) || item == "") return(0L)
    has_context = any(vapply(context_patterns, function(pattern) {
      grepl(pattern, item, ignore.case = TRUE, perl = TRUE)
    }, logical(1)), na.rm = TRUE)
    has_scope = grepl(scope_pattern, item, ignore.case = TRUE, perl = TRUE)
    if (isTRUE(has_context) && isTRUE(has_scope)) return(2L)
    if (isTRUE(has_context) || isTRUE(has_scope)) return(1L)
    0L
  }, integer(1), USE.NAMES = FALSE)
}

.plant_pubmed_context_evidence_rank = function(evidence_tier, source_database,
                                               confidence) {
  tier = .plant_normalize_evidence_tier(evidence_tier, source_database,
                                        "species")
  out = rep(0, length(tier))
  out[tier %in% c("manual_curated", "direct_species_database")] = 5
  out[tier == "direct_species_literature"] = 4
  out[tier %in% c("genus_database_fallback", "genus_literature_fallback")] = 3
  out[tier %in% c("family_database_fallback")] = 2
  out[tier %in% c("direct_species_pubtator_candidate", "unresolved")] = 1
  out + .plant_confidence_score(confidence)
}

.plant_pubmed_ids_for_doi = function(doi, cache, cache_dir, throttle,
                                     ncbi_email, ncbi_tool, ncbi_api_key,
                                     request_fun, request_timeout) {
  doi = .plant_normalize_doi(doi)
  if (is.na(doi) || doi == "") return(character())
  term = paste0('"', doi, '"[AID]')
  url = .plant_ncbi_url(
    endpoint = "esearch.fcgi",
    params = c(db = "pubmed", term = term, retmode = "json", retmax = "1",
               tool = ncbi_tool, email = ncbi_email),
    api_key = ncbi_api_key
  )
  search = tryCatch(
    .plant_fetch_json(url, cache, cache_dir, throttle, request_fun,
                      timeout = request_timeout),
    error = function(error) NULL
  )
  .plant_pubmed_ids(search)
}

.plant_pubmed_xml_sources = function(xml) {
  cols = .plant_literature_cols()
  xml = .uaf_first_non_empty_text(xml)
  if (is.na(xml) || xml == "") return(.uaf_empty_table(cols))
  matches = gregexpr("(?s)<PubmedArticle.*?</PubmedArticle>", xml,
                     perl = TRUE)
  blocks = regmatches(xml, matches)[[1]]
  if (length(blocks) < 1 || identical(blocks, "-1")) {
    blocks = xml
  }
  rows = lapply(blocks, function(block) {
    pmid = .plant_xml_first(block, "<PMID[^>]*>(.*?)</PMID>")
    title = .plant_xml_first(block, "<ArticleTitle[^>]*>(.*?)</ArticleTitle>")
    abstract_bits = .plant_xml_all(block,
                                   "<AbstractText[^>]*>(.*?)</AbstractText>")
    abstract = .pubchem_collapse(abstract_bits)
    doi = .plant_xml_first(
      block,
      "<ArticleId[^>]*IdType=[\"']doi[\"'][^>]*>(.*?)</ArticleId>"
    )
    data.frame(
      query_plant = NA_character_,
      query_plant_clean = NA_character_,
      species = NA_character_,
      genus = NA_character_,
      family = NA_character_,
      source_database = "PubMed",
      source_record_id = pmid,
      pmid = .plant_normalize_pmid(pmid),
      doi = .plant_normalize_doi(doi),
      title = title,
      abstract = abstract,
      chemical_mention = NA_character_,
      species_mention = NA_character_,
      evidence_text = .uaf_first_non_empty_text(title, abstract),
      evidence_url = ifelse(!is.na(.plant_normalize_pmid(pmid)),
                            paste0("https://pubmed.ncbi.nlm.nih.gov/",
                                   .plant_normalize_pmid(pmid), "/"),
                            NA_character_),
      retrieved_at = .plant_timestamp(),
      confidence = "medium",
      curation_flag = "source_context",
      evidence_tier = "direct_species_literature",
      stringsAsFactors = FALSE
    )
  })
  .plant_bind_tables(rows, cols)
}

.plant_xml_first = function(xml, pattern) {
  out = .plant_xml_all(xml, pattern)
  .uaf_first_non_empty_text(out)
}

.plant_xml_all = function(xml, pattern) {
  xml = .uaf_first_non_empty_text(xml)
  if (is.na(xml) || xml == "") return(character())
  matches = gregexpr(pattern, xml, perl = TRUE, ignore.case = TRUE)
  hits = regmatches(xml, matches)[[1]]
  if (length(hits) < 1 || identical(hits, "-1")) return(character())
  sub(pattern, "\\1", hits, perl = TRUE, ignore.case = TRUE) |>
    vapply(.plant_xml_text, character(1), USE.NAMES = FALSE)
}

.plant_xml_text = function(x) {
  x = gsub("<[^>]+>", " ", x, perl = TRUE)
  replacements = c("&amp;" = "&", "&lt;" = "<", "&gt;" = ">",
                   "&quot;" = "\"", "&apos;" = "'", "&#8211;" = "-",
                   "&#8212;" = "-", "&#x2013;" = "-", "&#x2014;" = "-")
  for (entity in names(replacements)) {
    x = gsub(entity, replacements[[entity]], x, fixed = TRUE)
  }
  .uaf_squish_text(x)
}

.plant_normalize_doi = function(x) {
  x = .uaf_squish_text(x)
  hit = regexpr("10\\.\\d{4,9}/[-._;()/:A-Za-z0-9]+", x, perl = TRUE,
                ignore.case = TRUE)
  out = rep(NA_character_, length(x))
  ok = !is.na(hit) & hit > 0
  out[ok] = tolower(regmatches(x, hit))
  out
}

.plant_normalize_pmid = function(x) {
  x = .uaf_squish_text(x)
  hit = regexpr("\\b\\d{1,9}\\b", x, perl = TRUE)
  out = rep(NA_character_, length(x))
  ok = !is.na(hit) & hit > 0
  out[ok] = regmatches(x, hit)
  out
}

#' Audit provider-specific biological-context coverage
#'
#' @description
#' `plantProviderContextAudit()` summarizes how much source-backed plant-part,
#' tissue, and method context each provider contributes. It is designed for
#' pilot runs and scale-up decisions: sparse or low-confidence providers are
#' flagged for source inspection, parser hardening, or manual curation before
#' context-aware matrices are interpreted.
#'
#' @param x Plant phytochemistry result, normalized occurrence table, or named
#' list containing `PlantCompoundOccurrences`.
#' @param context_evidence Optional `PlantContextEvidence` table. If omitted,
#' context evidence is taken from `x` when available or extracted from
#' occurrences.
#'
#' @return A `ProviderContextAudit` data frame with one row per source database.
#'
#' @export
plantProviderContextAudit = function(x, context_evidence = NULL) {
  occurrences = if (inherits(x, "uaf_plant_phytochemistry")) {
    x$PlantCompoundOccurrences
  } else if (is.list(x) && is.data.frame(x$PlantCompoundOccurrences)) {
    x$PlantCompoundOccurrences
  } else {
    x
  }
  if (is.null(context_evidence)) {
    context_evidence = if (is.list(x) && is.data.frame(x$PlantContextEvidence)) {
      x$PlantContextEvidence
    } else {
      plantContextEvidence(occurrences)
    }
  }
  .plant_provider_context_audit(occurrences, context_evidence)
}

#' Filter plant phytochemistry occurrence evidence
#'
#' @description
#' Keeps plant-compound occurrence rows that meet explicit evidence, confidence,
#' and biological-context criteria. When `x` is a plant phytochemistry result,
#' the returned result is rebuilt with updated occurrence, compound-resolution,
#' trait-evidence, summary, matrix, and validation tables. Literature candidate
#' tables are preserved for audit/provenance.
#'
#' @param x Plant phytochemistry result or normalized occurrence table.
#' @param occurrence_status Character vector of occurrence status values to
#' retain. Defaults to direct and curated reported records. Use `"all"` or
#' `NULL` to skip this filter.
#' @param analysis_ready Optional logical. If `TRUE`, retain only rows marked
#' analysis-ready; if `FALSE`, retain only non-analysis-ready rows; if `NULL`,
#' do not filter on the flag.
#' @param min_confidence Minimum confidence to retain. Use `NULL` to skip.
#' @param plant_part_group Optional plant-part groups to retain.
#' @param tissue_group Optional tissue groups to retain.
#' @param method_group Optional method groups to retain.
#' @param source_database Optional source databases to retain.
#' @param evidence_tier Optional evidence tiers to retain.
#'
#' @return A filtered plant phytochemistry result when `x` is a result object;
#' otherwise a filtered `PlantCompoundOccurrences` data frame.
#'
#' @export
filterPlantPhytochemistryEvidence = function(
    x,
    occurrence_status = c("direct_reported", "curated_reported"),
    analysis_ready = TRUE,
    min_confidence = "medium",
    plant_part_group = NULL,
    tissue_group = NULL,
    method_group = NULL,
    source_database = NULL,
    evidence_tier = NULL) {
  is_result = inherits(x, "uaf_plant_phytochemistry")
  occurrences = if (is_result) x$PlantCompoundOccurrences else x
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) {
    return(if (is_result) .plant_rebuild_filtered_result(x, occurrences) else
      occurrences)
  }
  keep = rep(TRUE, nrow(occurrences))
  status_filter = .uaf_non_empty(occurrence_status)
  if (length(status_filter) > 0 && !"all" %in% tolower(status_filter)) {
    keep = keep & occurrences$occurrence_status %in% status_filter
  }
  if (!is.null(analysis_ready)) {
    ready_value = .uaf_yes_no(isTRUE(analysis_ready))
    keep = keep & occurrences$analysis_ready == ready_value
  }
  if (!is.null(min_confidence)) {
    keep = keep & .plant_confidence_score(occurrences$confidence) >=
      .plant_confidence_score(min_confidence)
  }
  keep = .plant_apply_optional_filter(keep, occurrences$plant_part_group,
                                      plant_part_group)
  keep = .plant_apply_optional_filter(keep, occurrences$tissue_group,
                                      tissue_group)
  keep = .plant_apply_optional_filter(keep, occurrences$method_group,
                                      method_group)
  keep = .plant_apply_optional_filter(keep, occurrences$source_database,
                                      source_database)
  keep = .plant_apply_optional_filter(keep, occurrences$evidence_tier,
                                      evidence_tier)
  filtered = occurrences[keep, , drop = FALSE]
  row.names(filtered) = NULL
  if (is_result) .plant_rebuild_filtered_result(x, filtered) else filtered
}

#' Create a plant phytochemistry evidence review table
#'
#' @description
#' Builds a review worksheet for candidate, fallback, literature, or unresolved
#' plant-compound occurrence evidence. The table is intended for human review:
#' rows are not promoted automatically. Fill `review_decision` and optional
#' proposed fields, then pass the completed table to
#' `applyPlantPhytochemistryReview()`.
#'
#' @param x Plant phytochemistry result or normalized occurrence table.
#' @param occurrence_status Occurrence statuses to include. Use `"all"` to
#' include every occurrence row.
#' @param include_analysis_ready Logical. If `FALSE`, rows already marked
#' analysis-ready are omitted unless `occurrence_status = "all"`.
#'
#' @return A data frame suitable for CSV export, editing, and re-import.
#'
#' @export
plantPhytochemistryReviewTable = function(
    x,
    occurrence_status = c("candidate", "taxon_fallback",
                          "literature_reported", "unresolved"),
    include_analysis_ready = FALSE) {
  occurrences = if (inherits(x, "uaf_plant_phytochemistry")) {
    x$PlantCompoundOccurrences
  } else {
    x
  }
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) {
    return(.uaf_empty_table(.plant_review_cols()))
  }
  status_filter = .uaf_non_empty(occurrence_status)
  keep = rep(TRUE, nrow(occurrences))
  if (length(status_filter) > 0 && !"all" %in% tolower(status_filter)) {
    keep = keep & occurrences$occurrence_status %in% status_filter
  }
  if (!isTRUE(include_analysis_ready) && !"all" %in% tolower(status_filter)) {
    keep = keep & occurrences$analysis_ready != "Yes"
  }
  review = occurrences[keep, , drop = FALSE]
  if (nrow(review) < 1) return(.uaf_empty_table(.plant_review_cols()))
  out = data.frame(
    review_id = sprintf("review_%04d", seq_len(nrow(review))),
    review_decision = "needs_review",
    reviewed_by = NA_character_,
    reviewed_at = NA_character_,
    review_note = NA_character_,
    proposed_evidence_tier = NA_character_,
    proposed_confidence = NA_character_,
    proposed_source_database = NA_character_,
    proposed_citation_or_url = NA_character_,
    proposed_plant_part = NA_character_,
    proposed_tissue = NA_character_,
    proposed_method = NA_character_,
    review_key = .plant_occurrence_key(review),
    review[, .plant_occurrence_cols(), drop = FALSE],
    stringsAsFactors = FALSE
  )
  row.names(out) = NULL
  out
}

#' Apply reviewed plant phytochemistry evidence decisions
#'
#' @description
#' Applies a completed review table created by
#' `plantPhytochemistryReviewTable()`. Supported decisions are
#' `"needs_review"`/`"keep"` (leave unchanged), `"keep_candidate"` (mark as
#' reviewed but keep candidate/fallback status), `"update_context"` (copy
#' proposed plant part/tissue/method/confidence/source fields),
#' `"promote_curated"` (promote to `manual_curated` evidence with review
#' provenance), and `"reject"`/`"exclude"` (remove the occurrence row).
#'
#' @param x Plant phytochemistry result or normalized occurrence table.
#' @param review_table Completed review table.
#' @param reviewer Optional reviewer name used when `reviewed_by` is empty.
#' @param require_citation Logical. If `TRUE`, `promote_curated` requires either
#' an existing evidence URL or `proposed_citation_or_url`.
#'
#' @return Updated plant phytochemistry result or occurrence table.
#'
#' @export
applyPlantPhytochemistryReview = function(x,
                                          review_table,
                                          reviewer = NA_character_,
                                          require_citation = TRUE) {
  if (!is.data.frame(review_table)) {
    stop("`review_table` must be a data frame.", call. = FALSE)
  }
  is_result = inherits(x, "uaf_plant_phytochemistry")
  occurrences = if (is_result) x$PlantCompoundOccurrences else x
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1 || nrow(review_table) < 1) {
    return(if (is_result) .plant_rebuild_filtered_result(x, occurrences) else
      occurrences)
  }
  names(review_table) = .plant_normalize_column_names(names(review_table))
  keep = rep(TRUE, nrow(occurrences))
  for (i in seq_len(nrow(review_table))) {
    review = review_table[i, , drop = FALSE]
    idx = .plant_review_match_indices(occurrences, review)
    if (length(idx) < 1) next
    decision = .plant_review_decision(review$review_decision)
    if (decision %in% c("reject", "exclude")) {
      keep[idx] = FALSE
      next
    }
    if (decision %in% c("needs_review", "keep", "")) next
    if (decision == "keep_candidate") {
      occurrences$curation_flag[idx] = "reviewed_candidate_kept"
      next
    }
    occurrences = .plant_apply_review_context(occurrences, idx, review)
    if (decision == "promote_curated") {
      citation = .uaf_first_non_empty_text(
        review$proposed_citation_or_url,
        occurrences$evidence_url[idx][[1]]
      )
      if (isTRUE(require_citation) && is.na(citation)) {
        warning("Skipping promotion for review row ", i,
                ": no citation or URL supplied.", call. = FALSE)
        next
      }
      occurrences$evidence_tier[idx] = "manual_curated"
      occurrences$confidence[idx] = .uaf_first_non_empty_text(
        review$proposed_confidence, "high"
      )
      occurrences$source_database[idx] = .uaf_first_non_empty_text(
        review$proposed_source_database, "manual_review"
      )
      occurrences$evidence_url[idx] = citation
      occurrences$curation_flag[idx] = "reviewed_curated"
      reviewer_value = .uaf_first_non_empty_text(review$reviewed_by,
                                                 reviewer)
      note = .uaf_first_non_empty_text(review$review_note)
      if (!is.na(reviewer_value) || !is.na(note)) {
        occurrences$evidence_text[idx] = .plant_truncate(
          .pubchem_collapse(c(occurrences$evidence_text[idx],
                              paste("reviewed_by", reviewer_value),
                              note)),
          1200
        )
      }
    }
  }
  updated = .plant_normalize_occurrences(occurrences[keep, , drop = FALSE])
  updated = .plant_collapse_occurrence_evidence(updated)
  if (!is_result) return(updated)
  out = .plant_rebuild_filtered_result(x, updated)
  out$Provenance = .plant_bind_tables(list(
    out$Provenance,
    .plant_provenance("evidence_review", "manual_review",
                      paste(unique(.uaf_non_empty(updated$species)),
                            collapse = "; "),
                      NA_character_, nrow(review_table),
                      "Human review table applied to plant occurrence evidence.")
  ), .plant_provenance_cols())
  out$Validation = validatePlantPhytochemistryResult(out)
  out
}

#' Join plant chemistry summaries to metadata
#'
#' @param metadata User metadata data frame.
#' @param plant_chemistry Plant phytochemistry result, summary table, or matrix.
#' @param species_col Species column in `metadata`.
#' @param join_level One of `"species"`, `"genus"`, or `"family"`.
#'
#' @return Joined data frame.
#'
#' @export
joinPlantChemistryMetadata = function(metadata,
                                      plant_chemistry,
                                      species_col = "species",
                                      join_level = "species") {
  if (!is.data.frame(metadata)) stop("`metadata` must be a data frame.",
                                     call. = FALSE)
  if (!species_col %in% names(metadata)) {
    stop("`species_col` was not found in metadata.", call. = FALSE)
  }
  chemistry = if (inherits(plant_chemistry, "uaf_plant_phytochemistry")) {
    plant_chemistry$SpeciesChemistrySummary
  } else {
    plant_chemistry
  }
  if (!is.data.frame(chemistry)) {
    stop("`plant_chemistry` must be a result object or data frame.",
         call. = FALSE)
  }
  if (!join_level %in% names(chemistry)) {
    stop("`join_level` was not found in the chemistry table.", call. = FALSE)
  }
  metadata$.uaf_join_key = .plant_clean_name(metadata[[species_col]])
  chemistry$.uaf_join_key = .plant_clean_name(chemistry[[join_level]])
  out = merge(metadata, chemistry, by = ".uaf_join_key", all.x = TRUE)
  out$.uaf_join_key = NULL
  row.names(out) = NULL
  out
}

#' Score plant chemistry candidates
#'
#' @param x Plant phytochemistry result or species summary table.
#' @param weights Named numeric vector for transparent score components.
#'
#' @return Data frame with component scores and total prioritization score.
#'
#' @export
scorePlantChemistryCandidates = function(x,
                                         weights = c(direct_species = 0.25,
                                                     resolved_compounds = 0.20,
                                                     source_coverage = 0.20,
                                                     literature = 0.15,
                                                     traits = 0.20)) {
  summary = if (inherits(x, "uaf_plant_phytochemistry")) {
    x$SpeciesChemistrySummary
  } else {
    x
  }
  if (!is.data.frame(summary) || nrow(summary) < 1) {
    return(.uaf_empty_table(c("species", "chemistry_priority_score")))
  }
  weights = weights[!is.na(weights)]
  if (length(weights) < 1 || sum(weights) <= 0) {
    weights = c(direct_species = 1)
  }
  weights = weights / sum(weights)
  norm = function(value) {
    value = suppressWarnings(as.numeric(value))
    max_value = max(value, na.rm = TRUE)
    if (!is.finite(max_value) || max_value <= 0) return(rep(0, length(value)))
    pmin(value / max_value, 1)
  }
  components = data.frame(
    species = summary$species,
    direct_species_score = norm(summary$direct_species_compound_count),
    resolved_compound_score = norm(summary$resolved_compound_count),
    source_coverage_score = suppressWarnings(as.numeric(summary$source_coverage_score)),
    literature_score = norm(summary$literature_candidate_count),
    trait_score = norm(summary$chemical_trait_count),
    stringsAsFactors = FALSE
  )
  components$source_coverage_score[is.na(components$source_coverage_score)] = 0
  total = rep(0, nrow(components))
  if ("direct_species" %in% names(weights)) {
    total = total + weights[["direct_species"]] * components$direct_species_score
  }
  if ("resolved_compounds" %in% names(weights)) {
    total = total + weights[["resolved_compounds"]] * components$resolved_compound_score
  }
  if ("source_coverage" %in% names(weights)) {
    total = total + weights[["source_coverage"]] * components$source_coverage_score
  }
  if ("literature" %in% names(weights)) {
    total = total + weights[["literature"]] * components$literature_score
  }
  if ("traits" %in% names(weights)) {
    total = total + weights[["traits"]] * components$trait_score
  }
  components$chemistry_priority_score = round(total, 4)
  components$score_interpretation = "Prioritization aid only; not evidence of efficacy, occurrence completeness, remediation performance, or causal mechanism."
  components
}

#' Export plant phytochemistry result tables
#'
#' @param x Plant phytochemistry result.
#' @param path Output directory for `format = "csv"` or `.xlsx` path for
#' `format = "xlsx"`.
#' @param format One of `"csv"`, `"xlsx"`, or `"auto"`.
#' @param tables Optional table names to export.
#' @param preset One of `"all"` or `"analysis_ready"`. The latter filters to
#' direct/curated analysis-ready occurrence evidence before export while
#' preserving validation and provenance tables.
#' @param include_empty Logical. If `TRUE`, include empty schema tables.
#' @param overwrite Logical. If `TRUE`, replace existing output.
#' @param max_cell_chars Maximum characters retained in any character cell.
#'
#' @return Export manifest data frame.
#'
#' @export
exportPlantPhytochemistryWorkbook = function(x,
                                             path,
                                             format = c("auto", "xlsx", "csv"),
                                             tables = NULL,
                                             preset = c("all", "analysis_ready"),
                                             include_empty = TRUE,
                                             overwrite = FALSE,
                                             max_cell_chars = 30000) {
  format = match.arg(format)
  preset = match.arg(preset)
  if (missing(path) || is.null(path) || length(.uaf_non_empty(path)) < 1) {
    stop("`path` is required.", call. = FALSE)
  }
  if (!is.list(x) || is.data.frame(x)) {
    stop("`x` must be a plant phytochemistry result list.", call. = FALSE)
  }
  if (preset == "analysis_ready") {
    x = filterPlantPhytochemistryEvidence(
      x,
      occurrence_status = c("direct_reported", "curated_reported"),
      analysis_ready = TRUE,
      min_confidence = "medium"
    )
  }
  resolved = .categorate_export_resolve_path(path, format)
  export_tables = .plant_export_tables(x, tables, include_empty, max_cell_chars)
  if (length(export_tables) < 1) stop("No exportable tables found.",
                                      call. = FALSE)

  all_names = c("ExportManifest", names(export_tables))
  sheet_names = .categorate_export_sheet_names(all_names)
  file_names = .categorate_export_file_names(all_names)
  created_at = .plant_timestamp()
  manifest = data.frame(
    Table = all_names,
    SheetName = sheet_names,
    FileName = if (resolved$format == "csv") file_names else NA_character_,
    RowCount = c(NA_integer_, vapply(export_tables, nrow, integer(1))),
    ColumnCount = c(NA_integer_, vapply(export_tables, ncol, integer(1))),
    Format = resolved$format,
    OutputPath = resolved$path,
    CreatedAt = created_at,
    stringsAsFactors = FALSE
  )
  manifest$RowCount[manifest$Table == "ExportManifest"] = nrow(manifest)
  manifest$ColumnCount[manifest$Table == "ExportManifest"] = ncol(manifest)
  full_tables = c(list(ExportManifest = manifest), export_tables)
  names(full_tables) = sheet_names
  if (resolved$format == "xlsx") {
    .categorate_write_xlsx(full_tables, resolved$path, overwrite)
  } else {
    .categorate_write_csv_bundle(full_tables, file_names, resolved$path,
                                 overwrite)
  }
  row.names(manifest) = NULL
  manifest
}

#' @export
print.uaf_plant_phytochemistry = function(x, ...) {
  cat("uafR plant phytochemistry result\n")
  cat("  plant queries: ", nrow(x$PlantQueries), "\n", sep = "")
  cat("  occurrence rows: ", nrow(x$PlantCompoundOccurrences), "\n", sep = "")
  cat("  literature candidates: ", nrow(x$LiteratureCandidates), "\n", sep = "")
  cat("  resolved compounds: ",
      sum(x$CompoundResolution$resolved %in% TRUE), "\n", sep = "")
  if (is.list(x$Validation) && is.data.frame(x$Validation$Summary)) {
    cat("  validation status: ", x$Validation$Summary$Status[[1]], "\n",
        sep = "")
  }
  invisible(x)
}

.plant_data_dictionary = function() {
	  specs = list(
	    .plant_schema("PlantQueries", .plant_query_cols(),
                  required = c("query_id", "query_plant",
                               "query_plant_clean", "species", "species_slug"),
	                  role = "query",
	                  description = "Input plant names and query-normalized fields."),
	    .plant_schema("PlantQueryAliases", .plant_query_alias_cols(),
	                  required = c("query_id", "query_plant", "species",
	                               "alias", "alias_clean", "alias_type"),
	                  role = "query_alias",
	                  description = "Submitted, normalized, and verified plant-name aliases used for provider matching."),
    .plant_schema("PlantNameResolution", .plant_name_resolution_cols(),
                  required = c("query_id", "query_plant",
                               "query_status", "matched_taxon"),
                  role = "taxonomy",
                  description = "Plant name parsing and resolution status."),
	    .plant_schema("ProviderDiagnostics", .plant_provider_diagnostic_cols(),
                  required = c("provider", "enabled", "queried",
                               "record_count", "error_count", "retrieved_at"),
	                  role = "diagnostics",
	                  description = "Provider availability, requests, records, and warnings."),
	    .plant_schema("ProviderQueryAccounting",
	                  .plant_provider_query_accounting_cols(),
	                  required = c("query_id", "species", "provider",
	                               "query_status", "occurrence_count",
	                               "provider_status", "retry_required"),
	                  role = "query_accounting",
	                  allowed = list(
	                    query_status = c("records", "no_records",
	                                     "retry_required"),
	                    retry_required = .plant_yes_no_values()
	                  ),
	                  description = "One row per plant-provider combination distinguishing records, confirmed no-hit results, and retry-required queries."),
	    .plant_schema("ProviderResourceManifest",
	                  .plant_provider_resource_manifest_cols(),
	                  required = c("provider", "resource_type", "resource_id",
	                               "availability_status", "checked_at"),
	                  role = "provider_resource",
	                  allowed = list(
	                    availability_status = c("available", "unavailable")
	                  ),
	                  description = "Versioned local-index and live-service resources configured for the run."),
    .plant_schema("PlantCompoundOccurrences", .plant_occurrence_cols(),
                  required = c("query_plant", "query_plant_clean",
                               "species", "compound_name",
                               "compound_name_clean", "source_database",
                               "confidence", "evidence_tier",
                               "occurrence_status", "occurrence_basis",
                               "analysis_ready", "evidence_quality_score"),
                  role = "occurrence",
                  allowed = list(
                    confidence = .plant_confidence_values(),
                    evidence_tier = .plant_evidence_tiers(),
                    plant_part_group = .plant_context_group_values(),
                    tissue_group = .plant_context_group_values(),
                    method_group = .plant_method_group_values(),
                    biological_context_status = .plant_context_status_values(),
                    occurrence_status = .plant_occurrence_status_values(),
                    occurrence_basis = .plant_occurrence_basis_values(),
                    analysis_ready = .plant_yes_no_values(),
                    matched_rank = c("species", "genus", "family",
                                     "unknown")
                  ),
                  description = "Normalized species-compound occurrence evidence."),
	    .plant_schema("PlantContextEvidence", .plant_context_evidence_cols(),
	                  required = c("species", "compound_name_clean",
	                               "context_type", "normalized_context",
                               "context_confidence", "source_field",
                               "extraction_rule"),
                  role = "context",
                  allowed = list(
                    context_type = .plant_context_type_values(),
                    context_confidence = .plant_confidence_values(),
                    requires_review = .plant_yes_no_values()
	                  ),
	                  description = "Source-backed plant part, tissue, and method context evidence extracted from occurrence rows."),
	    .plant_schema("ProviderContextAudit", .plant_provider_context_audit_cols(),
	                  required = c("source_database", "occurrence_count",
	                               "context_known_fraction",
	                               "audit_status"),
	                  role = "context_audit",
	                  allowed = list(
	                    audit_status = c("context_ready", "context_sparse",
	                                     "review_heavy", "candidate_only",
	                                     "no_occurrences")
	                  ),
	                  description = "Provider-level biological-context coverage, review burden, and recommended parser or curation action."),
	    .plant_schema("LiteratureCandidates", .plant_literature_cols(),
                  required = c("query_plant", "species", "source_database",
                               "evidence_tier", "confidence"),
                  role = "literature",
	                  allowed = list(confidence = .plant_confidence_values(),
	                                 evidence_tier = .plant_evidence_tiers()),
	                  description = "Publication and annotation candidates; not confirmed occurrence by default."),
	    .plant_schema("SourceCompoundIdentity",
	                  .plant_source_compound_identity_cols(),
	                  required = c("compound_name", "compound_name_clean",
	                               "source_database", "identity_status"),
	                  role = "source_compound_identity",
	                  description = "Chemical identifiers and structures supplied directly by occurrence providers before PubChem resolution."),
    .plant_schema("CompoundResolution", .plant_compound_resolution_cols(),
                  required = c("compound_name", "compound_name_clean",
                               "resolved"),
                  role = "compound_identity",
                  description = "Compound identity resolution from uafR enrichment."),
    .plant_schema("CompoundIdentityReview",
                  .plant_compound_identity_review_cols(),
                  required = c("review_id", "compound_name",
                               "compound_name_clean", "identity_issue_type",
                               "recommended_decision"),
                  role = "compound_identity_review",
                  description = "Ranked review worksheet for unresolved or ambiguous compound identities."),
    .plant_schema("SpeciesChemistrySummary", .plant_summary_cols(),
                  required = c("species", "species_slug", "compound_count",
                               "source_coverage_score"),
                  role = "summary",
                  description = "One row per species with counts, coverage, traits, and caveats."),
    .plant_schema("SpeciesChemistryMatrix", c("species"),
                  required = c("species"),
                  role = "matrix",
                  description = "Wide species-feature matrix; columns vary by profile and data."),
    .plant_schema("ChemistryComparability", .plant_comparability_cols(),
                  required = c("species", "compound_name_clean",
                               "metabolism_domain", "biosynthetic_family",
                               "chemical_behavior", "comparison_scope",
                               "comparison_group",
                               "comparability_confidence",
                               "comparable_for_matrix"),
                  role = "comparability",
                  allowed = list(
                    metabolism_domain = .plant_metabolism_domain_values(),
                    biosynthetic_family = .plant_biosynthetic_family_values(),
                    chemical_behavior = .plant_chemical_behavior_values(),
                    comparison_scope = .plant_comparison_scope_values(),
                    biological_context_status = .plant_context_status_values(),
                    comparability_confidence = .plant_confidence_values(),
                    comparable_for_matrix = .plant_yes_no_values()
                  ),
                  description = "Occurrence-level chemistry scope annotations used to compare like with like."),
    .plant_schema("ComparableChemistryMatrix", c("species"),
                  required = c("species"),
                  role = "matrix",
                  description = "Wide species-feature matrix built from a declared comparable chemistry scope."),
    .plant_schema("TraitEvidence", .plant_trait_evidence_cols(),
                  required = c("species", "compound_name_clean",
                               "EvidenceType", "AnalysisKey",
                               "SourceDatabase", "Confidence"),
                  role = "evidence",
                  description = "Species-compound links to uafR chemical trait evidence."),
    .plant_schema("Validation", c("Summary", "TableQuality",
                                  "ProviderDiagnostics", "Issues",
                                  "DataDictionary"),
                  required = c("Summary", "TableQuality", "Issues"),
                  role = "validation",
                  description = "Nested validation result."),
    .plant_schema("DataDictionary", c("Table", "Column", "ExpectedType",
                                      "Required", "Role", "AllowedValues",
                                      "Description"),
                  required = c("Table", "Column", "Required"),
                  role = "dictionary",
                  description = "Schema table."),
    .plant_schema("Provenance", .plant_provenance_cols(),
                  required = c("step", "source", "retrieved_at",
                               "record_count"),
                  role = "provenance",
                  description = "Source, retrieval, and assembly provenance.")
  )
  out = do.call(rbind, specs)
  row.names(out) = NULL
  out
}

.plant_schema = function(table, cols, required = character(), role,
                         description, allowed = list()) {
  data.frame(
    Table = table,
    Column = cols,
    ExpectedType = .plant_expected_type(cols),
    Required = cols %in% required,
    Role = role,
    AllowedValues = vapply(cols, function(col) {
      if (!is.null(allowed[[col]])) paste(allowed[[col]], collapse = "; ") else ""
    }, character(1)),
    Description = description,
    stringsAsFactors = FALSE
  )
}

.plant_expected_type = function(cols) {
  integer_cols = grepl(
    paste(c("^input_order$", "_count$", "^RowCount$", "^ColumnCount$",
            "^request_count$", "^cache_hit_count$", "^record_count$",
            "^error_count$", "^warning_count$", "^timeout_count$",
            "^rate_limit_count$", "^max_pubmed_records$",
            "^max_provider_records$"),
          collapse = "|"),
    cols,
    ignore.case = TRUE
  )
  numeric_cols = grepl("score|fraction|elapsed_seconds|request_timeout",
                       cols, ignore.case = TRUE)
  ifelse(integer_cols,
         "integer",
         ifelse(numeric_cols, "numeric", "character"))
}

.plant_query_cols = function() {
  c("query_id", "input_order", "query_plant", "query_plant_clean",
    "species", "genus", "family", "species_slug", "taxon_fallback")
}

.plant_name_resolution_cols = function() {
  c("query_id", "query_plant", "query_plant_clean", "matched_taxon",
    "matched_rank", "species", "genus", "family", "species_slug",
    "query_status", "resolution_source", "resolution_note")
}

.plant_provider_diagnostic_cols = function() {
  c("provider", "enabled", "queried", "available", "request_count",
    "cache_hit_count", "record_count", "error_count", "warning_count",
    "timeout_count", "rate_limit_count", "status", "no_hit_reason",
    "warning_message", "message", "retrieved_at", "elapsed_seconds",
    "error_messages")
}

.plant_occurrence_cols = function() {
  c("query_plant", "query_plant_clean", "matched_taxon", "matched_rank",
    "species", "genus", "family", "compound_name", "compound_name_clean",
    "compound_id", "compound_id_type", "source_database",
    "source_record_id", "evidence_text", "evidence_url", "reference_id",
    "pmid", "doi", "plant_part", "plant_part_group", "tissue",
    "tissue_group", "method", "method_group", "biological_context_status",
    "occurrence_type", "occurrence_status", "occurrence_basis",
    "analysis_ready", "evidence_quality_score", "retrieved_at",
    "confidence", "curation_flag", "evidence_tier")
}

.plant_context_evidence_cols = function() {
  c("species", "species_slug", "genus", "family", "compound_name",
    "compound_name_clean", "source_database", "source_record_id", "pmid",
    "evidence_url", "context_type", "raw_context_text",
    "normalized_context", "source_field", "extraction_rule",
    "context_confidence", "evidence_basis", "requires_review",
    "retrieved_at")
}

.plant_provider_context_audit_cols = function() {
  c("source_database", "occurrence_count", "species_count",
    "compound_count", "direct_reported_count", "candidate_count",
    "taxon_fallback_count", "analysis_ready_count", "context_known_count",
    "context_known_fraction", "plant_part_known_count", "tissue_known_count",
    "method_known_count", "context_evidence_count",
    "high_confidence_context_count", "medium_confidence_context_count",
    "low_confidence_context_count", "review_required_context_count",
    "review_required_context_fraction", "context_missing_count",
    "context_missing_fraction", "top_context_groups", "top_extraction_rules",
    "audit_status", "recommended_action")
}

.plant_literature_cols = function() {
  c("query_plant", "query_plant_clean", "species", "genus", "family",
    "source_database", "source_record_id", "pmid", "doi", "title",
    "abstract", "chemical_mention", "species_mention", "evidence_text",
    "evidence_url", "literature_query", "literature_total_hit_count",
    "literature_search_truncated", "abstract_retrieval_status",
    "retrieved_at", "confidence", "curation_flag", "evidence_tier")
}

.plant_compound_resolution_cols = function() {
  c("compound_name", "compound_name_clean", "query_count", "resolved", "CID",
    "InChIKey", "SMILES", "MolecularFormula", "resolution_source", "notes")
}

.plant_source_compound_identity_cols = function() {
  c("compound_name", "compound_name_clean", "source_database",
    "source_record_id", "source_compound_id", "source_compound_id_type",
    "CID", "InChIKey", "SMILES", "MolecularFormula", "evidence_url",
    "evidence_text", "identity_status", "identity_note")
}

.plant_empty_source_compound_identity = function() {
  .uaf_empty_table(.plant_source_compound_identity_cols())
}

.plant_compound_identity_review_cols = function() {
  c("review_id", "review_decision", "reviewed_by", "reviewed_at",
    "review_note", "compound_name", "compound_name_clean", "query_count",
    "resolved", "CID", "InChIKey", "SMILES", "MolecularFormula",
    "resolution_source", "identity_issue_type", "review_priority",
    "recommended_decision", "proposed_compound_name", "proposed_cid",
    "proposed_inchikey", "proposed_smiles", "proposed_molecular_formula",
    "proposed_resolution_source", "suggested_query", "suggested_source",
    "source_record_count", "source_databases", "source_record_ids",
    "species_count", "example_species", "review_reason", "notes")
}

.plant_empty_compound_identity_review = function() {
  .uaf_empty_table(.plant_compound_identity_review_cols())
}

.plant_summary_cols = function() {
  c("species", "species_slug", "genus", "family", "query_status",
    "compound_count", "resolved_compound_count", "unresolved_compound_count",
    "analysis_ready_compound_count",
    "direct_species_compound_count", "genus_level_compound_count",
    "family_level_compound_count", "literature_candidate_count",
    "database_occurrence_count", "direct_reported_occurrence_count",
    "curated_reported_occurrence_count", "candidate_occurrence_count",
    "taxon_fallback_occurrence_count", "context_known_occurrence_count",
    "mean_evidence_quality_score", "analysis_ready_fraction",
    "plant_part_groups", "tissue_groups", "method_groups",
    "evidence_tier_summary",
    "source_database_count", "source_databases",
    "natural_product_superclasses", "natural_product_classes",
    "natural_product_subclasses", "dominant_compound_classes",
    "kingdoms_observed", "families_observed", "is_plant_occurring",
    "volatile_proxy_fraction", "lipophilic_fraction", "oxygenated_fraction",
    "nitrogenous_fraction", "sulfur_containing_fraction",
    "halogenated_fraction", "kegg_pathway_groups", "chemical_trait_count",
    "high_confidence_trait_count",
    "comparable_specialized_compound_count",
    "comparable_volatile_compound_count",
    "comparable_primary_compound_count",
    "comparable_lipid_compound_count",
    "comparable_unknown_compound_count", "comparison_scopes",
    "dominant_comparison_groups", "source_coverage_score",
    "uafR_validation_status")
}

.plant_comparability_cols = function() {
  c("species", "species_slug", "genus", "family", "compound_name",
    "compound_name_clean", "source_database", "source_record_id", "pmid",
    "evidence_url", "evidence_tier", "confidence", "occurrence_status",
    "analysis_ready", "matched_rank", "plant_part_group", "tissue_group",
    "method_group", "biological_context_status", "evidence_quality_score",
    "metabolism_domain", "biosynthetic_family", "chemical_behavior",
    "comparison_scope", "comparison_group", "comparison_subgroup",
    "comparability_confidence", "comparability_basis", "classification_source",
    "classification_source_table", "classification_source_field",
    "classification_source_value", "classification_source_confidence",
    "comparable_for_matrix", "comparison_caveat")
}

.plant_trait_evidence_cols = function() {
  c("species", "species_slug", "compound_name", "compound_name_clean",
    "EvidenceType", "AnalysisKey", "SourceDatabase", "EvidenceText",
    "EvidenceURL", "Confidence")
}

.plant_provenance_cols = function() {
  c("step", "source", "query", "source_url", "cache_file", "retrieved_at",
    "record_count", "notes")
}

.plant_validation_issue_cols = function() {
  c("severity", "table", "column", "issue", "expected", "observed",
    "row_count", "examples")
}

.plant_table_quality_cols = function() {
  c("Table", "Present", "RowCount", "ColumnCount", "MissingRequiredColumns",
    "MissingOptionalColumns", "DuplicateKeyColumns", "Status")
}

.plant_confidence_values = function() c("low", "medium", "high", "unknown")

.plant_yes_no_values = function() c("Yes", "No")

.plant_context_type_values = function() c("plant_part", "tissue", "method")

.plant_metabolism_domain_values = function() {
  c("specialized_metabolism", "primary_metabolism", "lipid_metabolism",
    "plant_hormone_signaling", "xenobiotic_or_contaminant",
    "broad_or_uncertain", "unknown")
}

.plant_biosynthetic_family_values = function() {
  c("terpenoid", "phenolic_phenylpropanoid", "flavonoid",
    "alkaloid_nitrogenous", "organosulfur", "glucosinolate",
    "benzoxazinoid", "cyanogenic_glycoside", "saponin",
    "steroid_triterpenoid", "polyketide", "fatty_acid_lipid",
    "carbohydrate", "amino_acid", "organic_acid",
    "nucleoside_nucleotide", "plant_hormone_signal",
    "broad_or_uncertain", "unknown")
}

.plant_chemical_behavior_values = function() {
  c("volatile_semivolatile", "nonvolatile_or_unspecified", "unknown")
}

.plant_comparison_scope_values = function() {
  c("specialized_metabolites", "volatile_specialized_metabolites",
    "primary_metabolites", "lipids_fatty_acids",
    "plant_hormone_signaling", "xenobiotic_or_contaminant",
    "broad_or_uncertain", "unknown")
}

.plant_context_group_values = function() {
  c("unknown", "root_belowground", "leaf", "stem_shoot", "bark_wood",
    "flower", "fruit_seed", "aerial", "whole_plant",
    "exudate_rhizosphere", "vascular", "secretory", "epidermal",
    "microbial", "extract_unspecified", "other")
}

.plant_method_group_values = function() {
  c("unknown", "gc_ms", "lc_ms", "hplc", "nmr", "mass_spectrometry",
    "chromatography", "spectroscopy", "database_record",
    "literature_curation", "other")
}

.plant_occurrence_status_values = function() {
  c("direct_reported", "curated_reported", "literature_reported",
    "candidate", "taxon_fallback", "unresolved", "unknown")
}

.plant_occurrence_basis_values = function() {
  c("species_database_record", "species_curated_record",
    "species_literature_record", "candidate_co_mention",
    "genus_fallback_record", "family_fallback_record", "unresolved",
    "unknown")
}

.plant_context_status_values = function() {
  c("plant_part_and_method_known", "plant_part_known", "method_known",
    "context_missing")
}

.plant_evidence_tiers = function() {
  c("direct_species_database", "direct_species_literature",
    "direct_species_pubtator_candidate", "genus_database_fallback",
    "genus_literature_fallback", "family_database_fallback",
    "manual_curated", "unresolved")
}

.plant_queries = function(plants, taxon_fallback = c("species", "genus")) {
  if (missing(plants) || is.null(plants)) {
    stop("`plants` must contain at least one plant name.", call. = FALSE)
  }
  if (is.data.frame(plants)) {
    names(plants) = .plant_normalize_column_names(names(plants))
    raw_species = if ("species" %in% names(plants)) plants$species else
      plants[[1]]
    species = .plant_canonical_taxon_name(raw_species)
    genus = if ("genus" %in% names(plants)) {
      .plant_canonical_taxon_name(plants$genus)
    } else {
      .plant_genus(species)
    }
    family = if ("family" %in% names(plants)) plants$family else
      rep(NA_character_, length(species))
  } else {
    raw_species = plants
    species = .plant_canonical_taxon_name(raw_species)
    genus = .plant_genus(species)
    family = rep(NA_character_, length(species))
  }
  raw_species = .uaf_squish_text(raw_species)
  species = .uaf_squish_text(species)
  keep = !is.na(species) & species != ""
  raw_species = raw_species[keep]
  species = species[keep]
  genus = .uaf_squish_text(genus)[keep]
  family = .uaf_squish_text(family)[keep]
  if (length(species) < 1) {
    stop("`plants` must contain at least one non-empty plant name.",
         call. = FALSE)
  }
  data.frame(
    query_id = paste0("plant_", seq_along(species)),
    input_order = seq_along(species),
    query_plant = raw_species,
    query_plant_clean = .plant_clean_name(raw_species),
    species = species,
    genus = genus,
    family = family,
    species_slug = .plant_slug(species),
    taxon_fallback = paste(.uaf_non_empty(taxon_fallback), collapse = "; "),
    stringsAsFactors = FALSE
  )
}

.plant_name_resolution = function(plant_queries) {
  out = data.frame(
    query_id = plant_queries$query_id,
    query_plant = plant_queries$query_plant,
    query_plant_clean = plant_queries$query_plant_clean,
    matched_taxon = plant_queries$species,
    matched_rank = "species",
    species = plant_queries$species,
    genus = plant_queries$genus,
    family = plant_queries$family,
    species_slug = plant_queries$species_slug,
    query_status = ifelse(.plant_is_binomial(plant_queries$species),
                          "parsed_species", "name_needs_review"),
    resolution_source = "input_name_parser",
    resolution_note = ifelse(.plant_is_binomial(plant_queries$species),
                             "Species-like binomial parsed from input.",
                             "Input is not a clear binomial; review taxon name."),
    stringsAsFactors = FALSE
  )
  row.names(out) = NULL
  out
}

.plant_normalize_sources = function(sources) {
  allowed = c("lotus", "knapsack", "npass", "pubchem", "pubmed", "pubtator")
  if (is.null(sources)) return(allowed)
  sources = tolower(.uaf_non_empty(sources))
  if (length(sources) < 1) return(character())
  bad = setdiff(sources, allowed)
  if (length(bad) > 0) {
    warning("Unsupported plant provider(s) skipped: ",
            paste(bad, collapse = ", "), call. = FALSE)
  }
  unique(intersect(sources, allowed))
}

.plant_provider_dispatch = function(plant_queries, plant_aliases, sources,
                                    cache, cache_dir,
                                    throttle, ncbi_email, ncbi_tool,
                                    ncbi_api_key, max_pubmed_records,
                                    max_provider_records, lotus_index,
                                    provider_indexes,
                                    provider_results,
                                    request_fun, pubtator_request_fun,
                                    request_timeout, refresh, progress) {
  occurrence_rows = list()
  literature_rows = list()
  diagnostic_rows = list()
  provenance_rows = list()
  accounting_rows = list()
  source_identity_rows = list()
  pubmed_literature = .uaf_empty_table(.plant_literature_cols())
  dispatch_sources = c(setdiff(sources, "pubtator"),
                       intersect(sources, "pubtator"))
  for (provider in dispatch_sources) {
    result = .plant_query_provider(
      provider = provider,
      plant_queries = plant_queries,
      plant_aliases = plant_aliases,
      cache = cache,
      cache_dir = cache_dir,
      throttle = throttle,
      ncbi_email = ncbi_email,
      ncbi_tool = ncbi_tool,
      ncbi_api_key = ncbi_api_key,
      max_pubmed_records = max_pubmed_records,
      max_provider_records = max_provider_records,
      lotus_index = lotus_index,
      provider_indexes = provider_indexes,
      pubmed_literature = pubmed_literature,
      provider_results = provider_results,
      request_fun = request_fun,
      pubtator_request_fun = pubtator_request_fun,
      request_timeout = request_timeout,
      refresh = refresh,
      progress = progress
    )
    occurrence_rows[[length(occurrence_rows) + 1]] =
      result$PlantCompoundOccurrences
    literature_rows[[length(literature_rows) + 1]] =
      result$LiteratureCandidates
    if (provider == "pubmed") pubmed_literature = result$LiteratureCandidates
    diagnostic_rows[[length(diagnostic_rows) + 1]] =
      result$ProviderDiagnostics
    provenance_rows[[length(provenance_rows) + 1]] = result$Provenance
    accounting_rows[[length(accounting_rows) + 1]] =
      .plant_provider_query_accounting(
        plant_queries, provider, result$PlantCompoundOccurrences,
        result$LiteratureCandidates, result$ProviderDiagnostics
      )
    source_identity_rows[[length(source_identity_rows) + 1]] =
      result$SourceCompoundIdentity
  }
  list(
    PlantCompoundOccurrences = .plant_bind_occurrences(occurrence_rows),
    LiteratureCandidates = .plant_bind_tables(literature_rows,
                                              .plant_literature_cols()),
    ProviderDiagnostics = .plant_bind_tables(diagnostic_rows,
                                             .plant_provider_diagnostic_cols()),
    ProviderQueryAccounting = .plant_bind_tables(
      accounting_rows, .plant_provider_query_accounting_cols()
    ),
    SourceCompoundIdentity = .plant_merge_source_identity_tables(
      source_identity_rows
    ),
    Provenance = .plant_bind_tables(provenance_rows, .plant_provenance_cols())
  )
}

.plant_query_provider = function(provider, plant_queries, plant_aliases,
                                 cache, cache_dir,
                                 throttle, ncbi_email, ncbi_tool,
                                 ncbi_api_key, max_pubmed_records,
                                 max_provider_records, lotus_index,
                                 provider_indexes,
                                 pubmed_literature,
                                 provider_results,
                                 request_fun, pubtator_request_fun,
                                 request_timeout, refresh, progress) {
  started = Sys.time()
  .plant_progress(progress, "uafR plant provider ", provider, ": querying")
  if (!is.null(provider_results) && !is.null(provider_results[[provider]])) {
    result = .plant_provider_from_result(provider, provider_results[[provider]],
                                         plant_queries)
    return(.plant_provider_finish(result, provider, started, progress))
  }
  result = switch(provider,
                  lotus = .plant_query_lotus(
                    plant_queries, cache, cache_dir, throttle, request_fun,
                    max_provider_records, request_timeout, lotus_index,
                    plant_aliases),
                  knapsack = .plant_query_knapsack(
                    plant_queries, cache, cache_dir, throttle, request_fun,
                    max_provider_records, request_timeout, plant_aliases),
                  npass = .plant_query_npass(
                    plant_queries, provider_indexes$npass,
                    max_provider_records, plant_aliases),
                  pubchem = .plant_query_pubchem_occurrences(
                    plant_queries, cache, cache_dir, throttle, ncbi_email,
                    ncbi_tool, ncbi_api_key, max_provider_records,
                    request_fun, request_timeout, plant_aliases),
                  pubmed = .plant_query_pubmed_literature(
                    plant_queries, cache, cache_dir, throttle, ncbi_email,
                    ncbi_tool, ncbi_api_key, max_pubmed_records,
                    request_fun, request_timeout, plant_aliases),
                  pubtator = .plant_query_pubtator_literature(
                    plant_queries, cache, cache_dir, throttle,
                    pubtator_request_fun, request_timeout,
                    pubmed_literature = pubmed_literature,
                    plant_aliases = plant_aliases),
                  .plant_provider_empty_result(provider, plant_queries,
                                               status = "not_implemented",
                                               message = paste0(
                                                 provider,
                                                 " species-first adapter is scaffolded ",
                                                 "but requires a public download/API ",
                                                 "parser or provider_results input."))
  )
  .plant_provider_finish(result, provider, started, progress)
}

.plant_provider_from_result = function(provider, result, plant_queries) {
  source_identity = .plant_empty_source_compound_identity()
  if (is.data.frame(result)) {
    occurrences = .plant_normalize_occurrences(result, source_hint = provider)
    literature = .uaf_empty_table(.plant_literature_cols())
  } else if (is.list(result)) {
    occurrences = .plant_normalize_occurrences(
      result$PlantCompoundOccurrences %||% result$occurrences %||%
        .uaf_empty_table(.plant_occurrence_cols()),
      source_hint = provider
    )
    source_identity = .plant_bind_tables(
      list(result$SourceCompoundIdentity %||% result$source_compound_identity),
      .plant_source_compound_identity_cols()
    )
    literature = .plant_normalize_literature(
      result$LiteratureCandidates %||% result$literature_candidates %||%
        .uaf_empty_table(.plant_literature_cols()),
      source_hint = provider
    )
  } else {
    occurrences = .plant_empty_occurrences()
    literature = .uaf_empty_table(.plant_literature_cols())
  }
  occurrences = .plant_match_occurrences_to_queries(occurrences, plant_queries)
  literature = .plant_match_literature_to_queries(literature, plant_queries)
  list(
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = literature,
    SourceCompoundIdentity = source_identity,
    ProviderDiagnostics = .plant_provider_diagnostics(
      provider, TRUE, TRUE, TRUE, 0, 0,
      nrow(occurrences) + nrow(literature), 0, 0, "ok",
      "Provider rows supplied by provider_results."),
    Provenance = .plant_provenance("provider_dispatch", provider,
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   NA_character_,
                                   nrow(occurrences) + nrow(literature),
                                   "Provider rows supplied locally.")
  )
}

.plant_query_lotus = function(plant_queries, cache, cache_dir, throttle,
                              request_fun, max_records, request_timeout,
                              lotus_index = NULL, plant_aliases = NULL) {
  lotus_index = .plant_lotus_index_or_null(lotus_index)
  if (!is.null(lotus_index)) {
    return(.plant_query_lotus_local(
      plant_queries, lotus_index, max_records, plant_aliases
    ))
  }
  rows = list()
  request_count = 0
  cache_hit_count = 0L
  errors = 0
  error_messages = character()
  consecutive_busy = 0L
  circuit_open = FALSE
  requested_max_records = .plant_max_records(max_records)
  effective_max_records = if (is.null(request_fun)) {
    .plant_lotus_live_record_cap(requested_max_records)
  } else {
    requested_max_records
  }
  for (i in seq_len(nrow(plant_queries))) {
    query_row = plant_queries[i, , drop = FALSE]
    terms = .plant_taxon_query_terms(query_row, plant_aliases)
    terms = terms[terms$rank == "species", , drop = FALSE]
    for (j in seq_len(nrow(terms))) {
      plant = terms$term[[j]]
      url = paste0("https://lotus.naturalproducts.net/api/search/simple?query=",
                   utils::URLencode(plant, reserved = TRUE))
      request_count = request_count + 1
      result = tryCatch({
        .plant_fetch_lotus_json(url, cache, file.path(cache_dir, "lotus"),
                                throttle, request_fun,
                                timeout = request_timeout,
                                max_records = effective_max_records)
      }, error = function(error) {
        errors <<- errors + 1
        error_messages <<- c(error_messages, conditionMessage(error))
        consecutive_busy <<- .plant_provider_busy_update(consecutive_busy,
                                                         error)
        NULL
      })
      if (is.null(result) && .plant_provider_circuit_open(consecutive_busy)) {
        circuit_open = TRUE
        break
      }
      if (!is.null(result)) consecutive_busy = 0L
      if (!is.null(result)) {
        cache_hit_count = cache_hit_count + as.integer(.plant_cache_hit(result))
      }
      parsed = .plant_lotus_rows(
        query_row, result, url, effective_max_records,
        accepted_species = .plant_verified_species_aliases(plant_aliases,
                                                           query_row)
      )
      if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
    }
    if (circuit_open) break
  }
  occurrences = .plant_bind_occurrences(rows)
  status = .plant_provider_status(nrow(occurrences), errors)
  diagnostic_note = paste(
    "LOTUS simple API taxon-oriented candidates; only rows with taxon evidence are retained.",
    .plant_lotus_live_cap_note(requested_max_records, effective_max_records,
                               is.null(request_fun)),
    .plant_provider_circuit_note(circuit_open)
  )
  list(
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    ProviderDiagnostics = .plant_provider_diagnostics(
      "lotus", TRUE, TRUE, TRUE, request_count, cache_hit_count,
      nrow(occurrences), errors, 0, status,
      diagnostic_note,
      error_messages = .plant_error_messages(error_messages)),
    Provenance = .plant_provenance("provider_dispatch", "lotus",
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   "https://lotus.naturalproducts.net/api/search/simple",
                                   nrow(occurrences),
                                   "LOTUS Natural Products Online candidate search.")
  )
}

.plant_query_lotus_local = function(plant_queries, lotus_index, max_records,
                                    plant_aliases = NULL) {
  started = Sys.time()
  errors = 0L
  error_messages = character()
  lookup = .plant_lotus_lookup_info(lotus_index)
  occurrences = if (!is.null(lookup)) {
    tryCatch(
      .plant_query_lotus_lookup_index(
        plant_queries, lookup, max_records, plant_aliases
      ),
      error = function(error) {
        errors <<- errors + 1L
        error_messages <<- conditionMessage(error)
        .plant_empty_occurrences()
      }
    )
  } else {
    index = tryCatch(
      standardizeLotusIndex(lotus_index),
      error = function(error) {
        errors <<- errors + 1L
        error_messages <<- conditionMessage(error)
        .plant_empty_lotus_index()
      }
    )
    if (nrow(index) > 0) {
      .plant_query_lotus_index(
        plant_queries, index, max_records, plant_aliases = plant_aliases
      )
    } else {
      .plant_empty_occurrences()
    }
  }
  source_identity = attr(occurrences, "SourceCompoundIdentity")
  if (!is.data.frame(source_identity)) {
    source_identity = .plant_empty_source_compound_identity()
  }
  status = .plant_provider_status(nrow(occurrences), errors)
  message = if (errors > 0) {
    "Local LOTUS index could not be read; no live LOTUS fallback was attempted because an explicit index was supplied."
  } else if (!is.null(lookup)) {
    paste("LOTUS records were queried from a manifest-backed local lookup index.",
          "Only candidate species/genus/family shards were read.")
  } else {
    paste("LOTUS records were queried from a local standardized index.",
          "This is the recommended path for medium and large plant panels.")
  }
  list(
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    SourceCompoundIdentity = source_identity,
    ProviderDiagnostics = .plant_provider_diagnostics(
      "lotus", TRUE, TRUE, errors < 1, 0, 0,
      nrow(occurrences), errors, 0, status, message,
      elapsed_seconds = round(as.numeric(difftime(Sys.time(), started,
                                                  units = "secs")), 3),
      error_messages = .plant_error_messages(error_messages)),
    Provenance = .plant_provenance(
      "provider_dispatch", "lotus_local_index",
      paste(plant_queries$query_plant, collapse = "; "),
      .plant_lotus_index_source_label(lotus_index), nrow(occurrences),
      "LOTUS species-compound occurrences queried from a local index.")
  )
}

.plant_lotus_live_record_cap = function(max_records) {
  requested = .plant_max_records(max_records)
  cap = suppressWarnings(as.numeric(
    Sys.getenv("UAFR_LOTUS_MAX_RECORDS_PER_SPECIES", "2")
  ))
  if (!is.finite(cap) || cap < 1) cap = 2
  min(requested, as.integer(cap))
}

.plant_lotus_live_cap_note = function(requested_max_records,
                                      effective_max_records,
                                      live_request = TRUE) {
  if (!isTRUE(live_request)) return("")
  requested = .plant_max_records(requested_max_records)
  effective = .plant_max_records(effective_max_records)
  if (!is.finite(requested) || requested > effective) {
    return(paste0(
      "Live LOTUS reads are capped at ", effective,
      " record(s) per species because the unpaged simple API can return",
      " tens of megabytes per plant. Set UAFR_LOTUS_MAX_RECORDS_PER_SPECIES",
      " for a deeper diagnostic run."
    ))
  }
  ""
}

.plant_query_lotus_index = function(plant_queries, lotus_index,
                                    max_records = Inf,
                                    plant_aliases = NULL) {
  lookup = .plant_lotus_lookup_info(lotus_index)
  if (!is.null(lookup)) {
    return(.plant_query_lotus_lookup_index(
      plant_queries, lookup, max_records, plant_aliases
    ))
  }
  index = if (is.data.frame(lotus_index) &&
              all(.plant_lotus_index_cols() %in% names(lotus_index))) {
    lotus_index
  } else {
    standardizeLotusIndex(lotus_index)
  }
  .plant_query_lotus_index_data(
    plant_queries, index, max_records, plant_aliases
  )
}

.plant_query_lotus_index_data = function(plant_queries, index,
                                         max_records = Inf,
                                         plant_aliases = NULL) {
  if (nrow(index) < 1 || nrow(plant_queries) < 1) {
    return(.plant_empty_occurrences())
  }
  rows = list()
  identity_rows = list()
  max_records = .plant_max_records(max_records)
  for (i in seq_len(nrow(plant_queries))) {
    query = plant_queries[i, , drop = FALSE]
    matches = .plant_lotus_index_matches(index, query, plant_aliases)
    if (nrow(matches) < 1) next
    matches$.uaf_rank_order = match(matches$matched_rank,
                                    c("species", "genus", "family",
                                      "unknown"))
    matches = matches[order(matches$.uaf_rank_order,
                            matches$compound_name_clean,
                            matches$source_record_id), , drop = FALSE]
    if (is.finite(max_records)) {
      matches = utils::head(matches, max_records)
    }
    rows[[length(rows) + 1]] =
      .plant_lotus_index_occurrences(query, matches)
    identity_rows[[length(identity_rows) + 1]] =
      .plant_lotus_source_identity_from_matches(matches)
  }
  out = .plant_bind_occurrences(rows)
  attr(out, "SourceCompoundIdentity") = .plant_bind_unique_tables(
    identity_rows, .plant_source_compound_identity_cols(),
    c("source_database", "source_record_id", "source_compound_id",
      "InChIKey", "SMILES")
  )
  out
}

.plant_query_lotus_lookup_index = function(plant_queries, lookup,
                                           max_records = Inf,
                                           plant_aliases = NULL) {
  index = .plant_lotus_lookup_index_rows(
    plant_queries, lookup, plant_aliases
  )
  if (nrow(index) < 1) return(.plant_empty_occurrences())
  .plant_query_lotus_index_data(
    plant_queries, index, max_records, plant_aliases
  )
}

.plant_lotus_lookup_index_rows = function(plant_queries, lookup,
                                          plant_aliases = NULL) {
  lookup = .plant_lotus_lookup_info(lookup)
  if (is.null(lookup) || nrow(plant_queries) < 1) {
    return(.plant_empty_lotus_index())
  }
  keys = .plant_lotus_lookup_query_keys(plant_queries, plant_aliases)
  if (nrow(keys) < 1) return(.plant_empty_lotus_index())
  prefix_length = suppressWarnings(as.integer(
    lookup$manifest$shard_prefix_length[[1]]
  ))
  if (!is.finite(prefix_length) || prefix_length < 1) prefix_length = 2L
  layout = .uaf_first_non_empty_text(lookup$manifest$lookup_layout, "prefix")
  keys$shard_file = mapply(
    .plant_lotus_lookup_shard_file,
    index_key = keys$index_key,
    index_key_type = keys$index_key_type,
    MoreArgs = list(root = lookup$root,
                    prefix_length = prefix_length,
                    layout = layout),
    USE.NAMES = FALSE
  )
  files = unique(keys$shard_file[file.exists(keys$shard_file)])
  if (length(files) < 1) return(.plant_empty_lotus_index())
  key_values = paste(keys$index_key_type, keys$index_key, sep = "\r")
  parts = lapply(files, function(path) {
    x = utils::read.csv(path, stringsAsFactors = FALSE,
                        check.names = FALSE)
    names(x) = .plant_normalize_column_names(names(x))
    if (!all(c("index_key_type", "index_key") %in% names(x))) {
      return(.plant_empty_lotus_lookup_rows())
    }
    hit = paste(x$index_key_type, x$index_key, sep = "\r") %in% key_values
    x[hit, , drop = FALSE]
  })
  rows = .plant_bind_tables(parts, .plant_lotus_lookup_cols())
  if (nrow(rows) < 1) return(.plant_empty_lotus_index())
  .plant_finalize_lotus_index(rows[, .plant_lotus_index_cols(), drop = FALSE])
}

.plant_lotus_lookup_info = function(lotus_index) {
  if (is.list(lotus_index) && !is.null(lotus_index$manifest) &&
      !is.null(lotus_index$root)) {
    return(lotus_index)
  }
  if (!is.character(lotus_index) || length(lotus_index) != 1 ||
      is.na(lotus_index) || !nzchar(lotus_index)) {
    return(NULL)
  }
  manifest_file = NA_character_
  root = NA_character_
  if (dir.exists(lotus_index)) {
    root = lotus_index
    manifest_file = file.path(lotus_index, "manifest.json")
  } else if (file.exists(lotus_index) &&
             tolower(tools::file_ext(lotus_index)) == "json") {
    manifest_file = lotus_index
    root = dirname(lotus_index)
  }
  if (is.na(manifest_file) || !file.exists(manifest_file)) return(NULL)
  manifest = tryCatch(jsonlite::fromJSON(manifest_file, simplifyVector = TRUE),
                      error = function(error) NULL)
  if (is.null(manifest) ||
      !identical(.uaf_first_non_empty_text(manifest$format),
                 "uafR_lotus_lookup_index")) {
    return(NULL)
  }
  list(
    root = normalizePath(root, winslash = "/", mustWork = FALSE),
    manifest_file = normalizePath(manifest_file, winslash = "/",
                                  mustWork = FALSE),
    manifest = manifest
  )
}

.plant_lotus_lookup_cols = function() {
  c("index_key_type", "index_key", .plant_lotus_index_cols())
}

.plant_empty_lotus_lookup_rows = function() {
  .uaf_empty_table(.plant_lotus_lookup_cols())
}

.plant_lotus_lookup_query_keys = function(plant_queries,
                                          plant_aliases = NULL) {
  rows = list()
  for (i in seq_len(nrow(plant_queries))) {
    query = plant_queries[i, , drop = FALSE]
    query_species = .plant_verified_species_aliases(plant_aliases, query)
    query_genus = .uaf_first_non_empty_text(query$genus,
                                            .plant_genus(query_species))
    query_family = .uaf_first_non_empty_text(query$family)
    fallback = .uaf_non_empty(strsplit(
      .uaf_first_non_empty_text(query$taxon_fallback),
      ";",
      fixed = TRUE
    )[[1]])
    fallback = tolower(.uaf_squish_text(fallback))
    add_key = function(type, value) {
      value = .plant_clean_name(value)
      if (is.na(value) || value == "") return(NULL)
      data.frame(index_key_type = type, index_key = value,
                 stringsAsFactors = FALSE)
    }
    for (species_value in query_species) {
      rows[[length(rows) + 1]] = add_key("species", species_value)
    }
    if ("genus" %in% fallback) {
      rows[[length(rows) + 1]] = add_key("genus", query_genus)
    }
    if ("family" %in% fallback) {
      rows[[length(rows) + 1]] = add_key("family", query_family)
    }
  }
  out = .plant_bind_tables(rows, c("index_key_type", "index_key"))
  if (nrow(out) < 1) return(out)
  out = out[!duplicated(out), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_lotus_lookup_prefix = function(index_key, prefix_length = 2L) {
  key = .plant_clean_compound(index_key)
  key[is.na(key) | key == ""] = "_"
  substring(key, 1L, pmax(1L, prefix_length))
}

.plant_lotus_lookup_shard_file = function(root, index_key,
                                          index_key_type = NA_character_,
                                          prefix_length = 2L,
                                          layout = "prefix") {
  prefix = .plant_lotus_lookup_prefix(index_key, prefix_length)
  layout = .uaf_first_non_empty_text(layout, "prefix")
  if (identical(layout, "exact")) {
    key_type = .plant_clean_compound(index_key_type)
    if (is.na(key_type) || key_type == "") key_type = "unknown"
    key_slug = .plant_clean_compound(index_key)
    if (is.na(key_slug) || key_slug == "") key_slug = "_"
    return(file.path(root, "keys", key_type, prefix,
                     paste0(key_slug, ".csv")))
  }
  file.path(root, "shards", paste0(prefix, ".csv"))
}

.plant_lotus_index_matches = function(index, query, plant_aliases = NULL) {
  query_species = .plant_verified_species_aliases(plant_aliases, query)
  query_genus = .uaf_first_non_empty_text(query$genus, .plant_genus(query_species))
  query_family = .uaf_first_non_empty_text(query$family)
  fallback = .uaf_non_empty(strsplit(
    .uaf_first_non_empty_text(query$taxon_fallback),
    ";",
    fixed = TRUE
  )[[1]])
  fallback = tolower(.uaf_squish_text(fallback))

  pieces = list()
  species_hit = .plant_clean_name(index$species) %in%
    .plant_clean_name(query_species)
  species_hit[is.na(species_hit)] = FALSE
  if (any(species_hit, na.rm = TRUE)) {
    hit = index[species_hit, , drop = FALSE]
    hit$matched_rank = "species"
    hit$matched_taxon = hit$species
    pieces[[length(pieces) + 1]] = hit
  }

  if ("genus" %in% fallback && !is.na(query_genus) && query_genus != "") {
    genus_hit = .plant_clean_name(index$genus) == .plant_clean_name(query_genus)
    genus_hit[is.na(genus_hit)] = FALSE
    genus_hit = genus_hit & !species_hit
    if (any(genus_hit, na.rm = TRUE)) {
      hit = index[genus_hit, , drop = FALSE]
      hit$matched_rank = "genus"
      hit$matched_taxon = query_genus
      pieces[[length(pieces) + 1]] = hit
    }
  }

  if ("family" %in% fallback && !is.na(query_family) && query_family != "") {
    family_hit = .plant_clean_name(index$family) ==
      .plant_clean_name(query_family)
    family_hit[is.na(family_hit)] = FALSE
    family_hit = family_hit & !species_hit
    if (any(family_hit, na.rm = TRUE)) {
      hit = index[family_hit, , drop = FALSE]
      hit$matched_rank = "family"
      hit$matched_taxon = query_family
      pieces[[length(pieces) + 1]] = hit
    }
  }

  out = .plant_bind_tables(pieces, c(.plant_lotus_index_cols(),
                                     "matched_rank", "matched_taxon"))
  if (nrow(out) < 1) return(out)
  key_cols = c("matched_rank", "compound_name_clean", "source_record_id",
               "lotus_id", "species")
  out = out[!duplicated(out[, key_cols, drop = FALSE]), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_lotus_index_occurrences = function(query, matches) {
  rows = lapply(seq_len(nrow(matches)), function(i) {
    hit = matches[i, , drop = FALSE]
    rank = .uaf_first_non_empty_text(hit$matched_rank, "species")
    data.frame(
      query_plant = query$query_plant,
      query_plant_clean = query$query_plant_clean,
      matched_taxon = .uaf_first_non_empty_text(hit$matched_taxon,
                                                hit$species),
      matched_rank = rank,
      species = query$species,
      genus = query$genus,
      family = query$family,
      compound_name = hit$compound_name,
      compound_name_clean = hit$compound_name_clean,
      compound_id = .plant_lotus_index_compound_id(hit),
      compound_id_type = .plant_lotus_index_compound_id_type(hit),
      source_database = "LOTUS",
      source_record_id = .plant_lotus_index_source_record_id(hit),
      evidence_text = .uaf_first_non_empty_text(
        hit$evidence_text,
        .plant_lotus_index_evidence_text(hit)
      ),
      evidence_url = .plant_lotus_index_evidence_url(hit),
      reference_id = .uaf_first_non_empty_text(hit$reference_id,
                                               hit$lotus_id,
                                               hit$wikidata_id),
      pmid = hit$pmid,
      doi = hit$doi,
      plant_part = hit$plant_part,
      tissue = hit$tissue,
      method = hit$method,
      occurrence_type = "local_lotus_index_record",
      retrieved_at = .plant_timestamp(),
      confidence = ifelse(rank == "species", "high", "medium"),
      curation_flag = ifelse(rank == "species", "source_database_record",
                             "taxon_fallback_review_required"),
      evidence_tier = ifelse(rank == "species",
                             "direct_species_database",
                             paste0(rank, "_database_fallback")),
      stringsAsFactors = FALSE
    )
  })
  .plant_bind_occurrences(rows)
}

.plant_lotus_source_identity_from_matches = function(matches) {
  if (!is.data.frame(matches) || nrow(matches) < 1) {
    return(.plant_empty_source_compound_identity())
  }
  rows = lapply(seq_len(nrow(matches)), function(i) {
    hit = matches[i, , drop = FALSE]
    data.frame(
      compound_name = hit$compound_name,
      compound_name_clean = hit$compound_name_clean,
      source_database = "LOTUS",
      source_record_id = .plant_lotus_index_source_record_id(hit),
      source_compound_id = .plant_lotus_index_compound_id(hit),
      source_compound_id_type = .plant_lotus_index_compound_id_type(hit),
      CID = suppressWarnings(as.integer(.uaf_first_non_empty_text(hit$cid))),
      InChIKey = .uaf_first_non_empty_text(hit$inchikey),
      SMILES = .uaf_first_non_empty_text(hit$smiles),
      MolecularFormula = .uaf_first_non_empty_text(hit$molecular_formula),
      evidence_url = .plant_lotus_index_evidence_url(hit),
      evidence_text = .uaf_first_non_empty_text(
        hit$evidence_text, .plant_lotus_index_evidence_text(hit)
      ),
      identity_status = NA_character_,
      identity_note = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out = .plant_bind_tables(rows, .plant_source_compound_identity_cols())
  .plant_annotate_source_identity_status(unique(out))
}

.plant_lotus_index_or_null = function(lotus_index) {
  if (is.null(lotus_index)) return(NULL)
  if (is.data.frame(lotus_index)) return(lotus_index)
  if (is.character(lotus_index) && length(lotus_index) == 1 &&
      !is.na(lotus_index) && nzchar(lotus_index)) {
    return(lotus_index)
  }
  NULL
}

.plant_lotus_index_source_label = function(lotus_index) {
  if (is.character(lotus_index) && length(lotus_index) == 1 &&
      !is.na(lotus_index) && nzchar(lotus_index)) {
    return(lotus_index)
  }
  if (is.data.frame(lotus_index) && "source_file" %in% names(lotus_index)) {
    return(.uaf_first_non_empty_text(lotus_index$source_file,
                                     "data_frame"))
  }
  "data_frame"
}

.plant_lotus_index_signature = function(lotus_index) {
  lotus_index = .plant_lotus_index_or_null(lotus_index)
  if (is.null(lotus_index)) return("lotus_index:none")
  if (is.character(lotus_index)) {
    if (!file.exists(lotus_index)) return(paste0("lotus_index:missing:",
                                                 lotus_index))
    lookup = .plant_lotus_lookup_info(lotus_index)
    signature_path = if (is.null(lookup)) lotus_index else
      lookup$manifest_file
    info = file.info(signature_path)
    signature_label = if (is.null(lookup)) "path" else "lookup_manifest"
    return(paste("lotus_index", signature_label,
                 normalizePath(signature_path, winslash = "/",
                               mustWork = FALSE),
                 info$size, info$mtime, sep = ":"))
  }
  paste("lotus_index:data_frame", nrow(lotus_index), ncol(lotus_index),
        paste(names(lotus_index), collapse = ","), sep = ":")
}

.plant_query_knapsack = function(plant_queries, cache, cache_dir, throttle,
                                 request_fun, max_records, request_timeout,
                                 plant_aliases = NULL) {
  rows = list()
  request_count = 0
  cache_hit_count = 0L
  errors = 0
  error_messages = character()
  consecutive_busy = 0L
  circuit_open = FALSE
  effective_throttle = if (is.null(request_fun)) max(throttle, 1) else throttle
  for (i in seq_len(nrow(plant_queries))) {
    if (circuit_open) break
    query_row = plant_queries[i, , drop = FALSE]
    terms = .plant_taxon_query_terms(query_row, plant_aliases)
    terms = terms[terms$rank %in% c("species", "genus"), , drop = FALSE]
    for (j in seq_len(nrow(terms))) {
      url = paste0(
        "https://www.knapsackfamily.com/knapsack_core/info.php?sname=organism&word=",
        utils::URLencode(terms$term[[j]], reserved = TRUE)
      )
      request_count = request_count + 1
      html = tryCatch({
        .plant_fetch_text(url, cache, file.path(cache_dir, "knapsack"),
                          effective_throttle, request_fun,
                          timeout = request_timeout)
      }, error = function(error) {
        errors <<- errors + 1
        error_messages <<- c(error_messages, conditionMessage(error))
        consecutive_busy <<- .plant_provider_busy_update(consecutive_busy,
                                                         error)
        NULL
      })
      if (is.null(html) && .plant_provider_circuit_open(consecutive_busy)) {
        circuit_open = TRUE
        break
      }
      if (!is.null(html)) consecutive_busy = 0L
      if (!is.null(html)) {
        cache_hit_count = cache_hit_count + as.integer(.plant_cache_hit(html))
      }
      result_url = .plant_knapsack_result_url(html, terms$term[[j]])
      if (!is.na(result_url) && !identical(result_url, url)) {
        request_count = request_count + 1L
        html = tryCatch({
          .plant_fetch_text(
            result_url, cache, file.path(cache_dir, "knapsack"),
            effective_throttle, request_fun, timeout = request_timeout
          )
        }, error = function(error) {
          errors <<- errors + 1
          error_messages <<- c(error_messages, conditionMessage(error))
          consecutive_busy <<- .plant_provider_busy_update(consecutive_busy,
                                                           error)
          NULL
        })
        if (is.null(html) && .plant_provider_circuit_open(consecutive_busy)) {
          circuit_open = TRUE
          break
        }
        if (!is.null(html)) consecutive_busy = 0L
        if (!is.null(html)) {
          cache_hit_count = cache_hit_count + as.integer(.plant_cache_hit(html))
        }
      } else {
        result_url = url
      }
      parsed = .plant_knapsack_rows(query_row, html, terms$term[[j]],
                                    terms$rank[[j]], result_url, max_records)
      if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
    }
  }
  occurrences = .plant_bind_occurrences(rows)
  status = .plant_provider_status(nrow(occurrences), errors)
  list(
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    ProviderDiagnostics = .plant_provider_diagnostics(
      "knapsack", TRUE, TRUE, TRUE, request_count, cache_hit_count,
      nrow(occurrences), errors, 0, status,
      paste("KNApSAcK documented organism endpoint records were exact-taxon",
            "filtered and normalized from public HTML output.",
            .plant_provider_circuit_note(circuit_open)),
      error_messages = .plant_error_messages(error_messages)),
    Provenance = .plant_provenance("provider_dispatch", "knapsack",
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   "https://www.knapsackfamily.com/knapsack_core/info.php",
                                   nrow(occurrences),
                                   "KNApSAcK organism-metabolite candidate search.")
  )
}

.plant_query_npass = function(plant_queries, npass_index, max_records,
                              plant_aliases = NULL) {
  if (is.null(.plant_npass_index_or_null(npass_index))) {
    return(.plant_provider_empty_result(
      "npass", plant_queries, status = "unavailable",
      message = paste(
        "NPASS requires a local NPASS 3.0 index. Build one with",
        "`buildNpassIndex()` and supply it as `provider_indexes$npass`."
      )
    ))
  }
  started = Sys.time()
  result = tryCatch(
    .plant_npass_query_result(
      plant_queries, npass_index, max_records, plant_aliases
    ),
    error = function(error) error
  )
  if (inherits(result, "error")) {
    return(list(
      PlantCompoundOccurrences = .plant_empty_occurrences(),
      LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
      SourceCompoundIdentity = .plant_empty_source_compound_identity(),
      ProviderDiagnostics = .plant_provider_diagnostics(
        "npass", TRUE, TRUE, FALSE, 0, 0, 0, 1, 0, "error",
        "The local NPASS index query failed.",
        elapsed_seconds = round(as.numeric(difftime(Sys.time(), started,
                                                    units = "secs")), 3),
        error_messages = conditionMessage(result)
      ),
      Provenance = .plant_provenance(
        "provider_dispatch", "npass_local_index",
        paste(plant_queries$query_plant, collapse = "; "),
        .uaf_first_non_empty_text(npass_index), 0,
        "NPASS local index query failed."
      )
    ))
  }
  occurrences = result$PlantCompoundOccurrences
  list(
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    SourceCompoundIdentity = result$SourceCompoundIdentity,
    ProviderDiagnostics = .plant_provider_diagnostics(
      "npass", TRUE, TRUE, TRUE, 0, 0, nrow(occurrences), 0, 0,
      .plant_provider_status(nrow(occurrences), 0),
      paste("NPASS records were queried from a manifest-backed local",
            "NPASS 3.0 species-source index."),
      elapsed_seconds = round(as.numeric(difftime(Sys.time(), started,
                                                  units = "secs")), 3)
    ),
    Provenance = .plant_provenance(
      "provider_dispatch", "npass_local_index",
      paste(plant_queries$query_plant, collapse = "; "),
      .uaf_first_non_empty_text(npass_index), nrow(occurrences),
      "NPASS species-compound occurrences queried from a local index."
    )
  )
}

.plant_query_pubchem_occurrences = function(plant_queries, cache, cache_dir,
                                            throttle, ncbi_email, ncbi_tool,
                                            ncbi_api_key, max_records,
                                            request_fun, request_timeout,
                                            plant_aliases = NULL) {
  rows = list()
  request_count = 0
  cache_hit_count = 0L
  errors = 0
  error_messages = character()
  consecutive_busy = 0L
  circuit_open = FALSE
  effective_throttle = throttle
  if (is.null(request_fun)) {
    has_ncbi_key = length(.uaf_non_empty(ncbi_api_key)) > 0
    effective_throttle = max(throttle, ifelse(has_ncbi_key, 0.10, 0.34))
  }
  for (i in seq_len(nrow(plant_queries))) {
    query_row = plant_queries[i, , drop = FALSE]
    accepted_species = .plant_verified_species_aliases(
      plant_aliases, query_row
    )
    taxon_match = NULL
    for (query_term in accepted_species) {
      request_count = request_count + 1L
      taxid_search = tryCatch({
        .plant_ncbi_taxonomy_search(query_term, cache, cache_dir,
                                    effective_throttle, ncbi_email, ncbi_tool,
                                    ncbi_api_key, request_fun,
                                    request_timeout)
      }, error = function(error) {
        errors <<- errors + 1L
        error_messages <<- c(error_messages, conditionMessage(error))
        consecutive_busy <<- .plant_provider_busy_update(consecutive_busy,
                                                         error)
        NULL
      })
      if (is.null(taxid_search) &&
          .plant_provider_circuit_open(consecutive_busy)) {
        circuit_open = TRUE
        break
      }
      if (!is.null(taxid_search)) consecutive_busy = 0L
      if (!is.null(taxid_search)) {
        cache_hit_count = cache_hit_count +
          as.integer(.plant_cache_hit(taxid_search))
      }
      suggest_matches = .plant_verified_taxonomy_suggest_matches(
        taxid_search, accepted_species
      )
      if (nrow(suggest_matches) > 0) {
        taxon_match = suggest_matches[1, , drop = FALSE]
        break
      }
      if (.plant_is_taxonomy_suggest_result(taxid_search)) next
      taxids = .plant_pubmed_ids(taxid_search)
      if (length(taxids) < 1) next
      request_count = request_count + 1L
      taxid_summary = tryCatch({
        .plant_ncbi_taxonomy_summary(
          taxids, cache, cache_dir, effective_throttle, ncbi_email,
          ncbi_tool, ncbi_api_key, request_fun, request_timeout
        )
      }, error = function(error) {
        errors <<- errors + 1L
        error_messages <<- c(error_messages, conditionMessage(error))
        consecutive_busy <<- .plant_provider_busy_update(consecutive_busy,
                                                         error)
        NULL
      })
      if (is.null(taxid_summary) &&
          .plant_provider_circuit_open(consecutive_busy)) {
        circuit_open = TRUE
        break
      }
      if (!is.null(taxid_summary)) consecutive_busy = 0L
      if (!is.null(taxid_summary)) {
        cache_hit_count = cache_hit_count +
          as.integer(.plant_cache_hit(taxid_summary))
      }
      matches = .plant_verified_taxonomy_matches(
        taxid_summary, taxids, accepted_species
      )
      if (nrow(matches) > 0) {
        taxon_match = matches[1, , drop = FALSE]
        break
      }
    }
    if (circuit_open) break
    if (is.null(taxon_match) || nrow(taxon_match) < 1) next
    taxid = taxon_match$taxid[[1]]
    matched_scientific_name = taxon_match$scientific_name[[1]]
    url = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/taxonomy/",
                 utils::URLencode(taxid, reserved = TRUE), "/JSON")
    request_count = request_count + 1
    result = tryCatch({
      .plant_fetch_json(url, cache, file.path(cache_dir, "pubchem_taxonomy"),
                        effective_throttle, request_fun,
                        timeout = request_timeout)
    }, error = function(error) {
      errors <<- errors + 1
      error_messages <<- c(error_messages, conditionMessage(error))
      consecutive_busy <<- .plant_provider_busy_update(consecutive_busy,
                                                       error)
      NULL
    })
    if (is.null(result) && .plant_provider_circuit_open(consecutive_busy)) {
      circuit_open = TRUE
      break
    }
    if (!is.null(result)) consecutive_busy = 0L
    if (!is.null(result)) {
      cache_hit_count = cache_hit_count + as.integer(.plant_cache_hit(result))
    }
    parsed = .plant_pubchem_taxonomy_rows(query_row, taxid, result, url,
                                          max_records)
    if (nrow(parsed) > 0) parsed$matched_taxon = matched_scientific_name
    parsed_parts = list(parsed)
    remaining = .plant_remaining_records(parsed_parts, max_records)
    specs = .plant_pubchem_taxonomy_external_specs(result)
    if (remaining > 0 && nrow(specs) > 0) {
      specs = specs[!duplicated(paste(specs$collection, specs$srccmpdkind,
                                      specs$view, sep = "\r")), ,
                    drop = FALSE]
      for (j in seq_len(nrow(specs))) {
        remaining = .plant_remaining_records(parsed_parts, max_records)
        if (remaining <= 0) break
        external_url = .plant_pubchem_taxonomy_external_url(
          taxid = taxid,
          srccmpdkind = specs$srccmpdkind[[j]],
          start = 1,
          limit = remaining
        )
        request_count = request_count + 1
        external_result = tryCatch({
          .plant_fetch_json(external_url, cache,
                            file.path(cache_dir, "pubchem_taxonomy_external"),
                            effective_throttle, request_fun,
                            timeout = request_timeout)
        }, error = function(error) {
          errors <<- errors + 1
          error_messages <<- c(error_messages, conditionMessage(error))
          consecutive_busy <<- .plant_provider_busy_update(consecutive_busy,
                                                           error)
          NULL
        })
        if (is.null(external_result) &&
            .plant_provider_circuit_open(consecutive_busy)) {
          circuit_open = TRUE
          break
        }
        if (!is.null(external_result)) consecutive_busy = 0L
        if (!is.null(external_result)) {
          cache_hit_count = cache_hit_count +
            as.integer(.plant_cache_hit(external_result))
        }
        external_rows = .plant_pubchem_taxonomy_external_rows(
          query_row, taxid, specs[j, , drop = FALSE], external_result,
          external_url, remaining
        )
        if (nrow(external_rows) > 0) {
          external_rows$matched_taxon = matched_scientific_name
        }
        if (nrow(external_rows) > 0) {
          parsed_parts[[length(parsed_parts) + 1]] = external_rows
        }
      }
    }
    parsed = .plant_bind_occurrences(parsed_parts)
    max_rows = .plant_max_records(max_records)
    if (is.finite(max_rows) && nrow(parsed) > max_rows) {
      parsed = parsed[seq_len(max_rows), , drop = FALSE]
    }
    if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
    if (circuit_open) break
  }
  occurrences = .plant_bind_occurrences(rows)
  status = .plant_provider_status(nrow(occurrences), errors)
  list(
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    ProviderDiagnostics = .plant_provider_diagnostics(
      "pubchem", TRUE, TRUE, TRUE, request_count, cache_hit_count,
      nrow(occurrences), errors, 0, status,
      paste("PubChem taxonomy PUG-View and external-table chemical annotations; review evidence context before treating as occurrence.",
            .plant_provider_circuit_note(circuit_open)),
      error_messages = .plant_error_messages(error_messages)),
    Provenance = .plant_provenance("provider_dispatch", "pubchem",
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/taxonomy/",
                                   nrow(occurrences),
                                   "PubChem taxonomy annotation search.")
  )
}

.plant_query_pubmed_literature = function(plant_queries, cache, cache_dir,
                                          throttle, ncbi_email, ncbi_tool,
                                          ncbi_api_key, max_pubmed_records,
                                          request_fun, request_timeout,
                                          plant_aliases = NULL) {
  rows = list()
  request_count = 0
  cache_hit_count = 0L
  errors = 0
  error_messages = character()
  consecutive_busy = 0L
  circuit_open = FALSE
  effective_throttle = throttle
  if (is.null(request_fun)) {
    has_ncbi_key = length(.uaf_non_empty(ncbi_api_key)) > 0
    effective_throttle = max(throttle, ifelse(has_ncbi_key, 0.10, 0.34))
  }
  for (i in seq_len(nrow(plant_queries))) {
    query_row = plant_queries[i, , drop = FALSE]
    accepted_species = .plant_verified_species_aliases(
      plant_aliases, query_row
    )
    term = .plant_pubmed_query(accepted_species)
    esearch_url = .plant_ncbi_url(
      endpoint = "esearch.fcgi",
      params = c(db = "pubmed", term = term, retmode = "json",
                 retmax = as.character(max_pubmed_records),
                 tool = ncbi_tool, email = ncbi_email),
      api_key = ncbi_api_key
    )
    request_count = request_count + 1
    search = tryCatch({
      .plant_fetch_json(esearch_url, cache, file.path(cache_dir, "pubmed"),
                        effective_throttle, request_fun,
                        timeout = request_timeout)
    }, error = function(error) {
      errors <<- errors + 1
      error_messages <<- c(error_messages, conditionMessage(error))
      consecutive_busy <<- .plant_provider_busy_update(consecutive_busy,
                                                       error)
      NULL
    })
    if (is.null(search) && .plant_provider_circuit_open(consecutive_busy)) {
      circuit_open = TRUE
      break
    }
    if (!is.null(search)) consecutive_busy = 0L
    if (!is.null(search)) {
      cache_hit_count = cache_hit_count + as.integer(.plant_cache_hit(search))
    }
    ids = .plant_pubmed_ids(search)
    total_hits = .plant_pubmed_total_hits(search, length(ids))
    truncated = is.finite(total_hits) && total_hits > length(ids)
    if (length(ids) < 1) next
    summary_url = .plant_ncbi_url(
      endpoint = "esummary.fcgi",
      params = c(db = "pubmed", id = paste(ids, collapse = ","),
                 retmode = "json", tool = ncbi_tool, email = ncbi_email),
      api_key = ncbi_api_key
    )
    request_count = request_count + 1
    summary = tryCatch({
      .plant_fetch_json(summary_url, cache, file.path(cache_dir, "pubmed"),
                        effective_throttle, request_fun,
                        timeout = request_timeout)
    }, error = function(error) {
      errors <<- errors + 1
      error_messages <<- c(error_messages, conditionMessage(error))
      consecutive_busy <<- .plant_provider_busy_update(consecutive_busy,
                                                       error)
      NULL
    })
    if (is.null(summary) && .plant_provider_circuit_open(consecutive_busy)) {
      circuit_open = TRUE
      break
    }
    if (!is.null(summary)) consecutive_busy = 0L
    if (!is.null(summary)) {
      cache_hit_count = cache_hit_count + as.integer(.plant_cache_hit(summary))
    }
    fetch_url = .plant_ncbi_url(
      endpoint = "efetch.fcgi",
      params = c(db = "pubmed", id = paste(ids, collapse = ","),
                 retmode = "xml", rettype = "abstract", tool = ncbi_tool,
                 email = ncbi_email),
      api_key = ncbi_api_key
    )
    request_count = request_count + 1L
    abstracts = tryCatch({
      xml = .plant_fetch_text(fetch_url, cache,
                              file.path(cache_dir, "pubmed_abstracts"),
                              effective_throttle, request_fun,
                              timeout = request_timeout)
      cache_hit_count = cache_hit_count + as.integer(.plant_cache_hit(xml))
      .plant_pubmed_abstracts(xml)
    }, error = function(error) {
      errors <<- errors + 1L
      error_messages <<- c(error_messages, conditionMessage(error))
      consecutive_busy <<- .plant_provider_busy_update(consecutive_busy, error)
      character()
    })
    if (length(abstracts) < 1 &&
        .plant_provider_circuit_open(consecutive_busy)) {
      circuit_open = TRUE
      break
    }
    if (length(abstracts) > 0) consecutive_busy = 0L
    rows[[length(rows) + 1]] =
      .plant_pubmed_summary_rows(query_row, ids,
                                 summary, abstracts = abstracts,
                                 total_hits = total_hits,
                                 truncated = truncated,
                                 query = term)
  }
  literature = .plant_bind_tables(rows, .plant_literature_cols())
  status = .plant_provider_status(nrow(literature), errors)
  list(
    PlantCompoundOccurrences = .plant_empty_occurrences(),
    LiteratureCandidates = literature,
    ProviderDiagnostics = .plant_provider_diagnostics(
      "pubmed", TRUE, TRUE, TRUE, request_count, cache_hit_count,
      nrow(literature), errors, 0, status,
      paste("PubMed E-utilities literature candidates; not confirmed occurrence.",
            .plant_provider_circuit_note(circuit_open)),
      error_messages = .plant_error_messages(error_messages)),
    Provenance = .plant_provenance("provider_dispatch", "pubmed",
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
                                   nrow(literature),
                                   "PubMed candidate literature search.")
  )
}

.plant_query_pubtator_literature = function(plant_queries, cache, cache_dir,
                                            throttle, request_fun,
                                            request_timeout,
                                            pubmed_literature = NULL,
                                            batch_size = 100,
                                            plant_aliases = NULL) {
  rows = list()
  request_count = 0
  cache_hit_count = 0L
  errors = 0
  error_messages = character()
  no_record_pmids = character()
  consecutive_busy = 0L
  circuit_open = FALSE
  effective_throttle = if (is.null(request_fun)) max(throttle, 0.34) else
    throttle
  pubmed_literature = .plant_normalize_literature(pubmed_literature, "pubmed")
  pmids = unique(.uaf_non_empty(pubmed_literature$pmid))
  invalid_pmids = pmids[!grepl("^[0-9]+$", pmids)]
  if (length(invalid_pmids) > 0L) {
    errors = errors + 1L
    error_messages = c(
      error_messages,
      paste("PubTator received non-numeric PMID values:",
            paste(invalid_pmids, collapse = ", "))
    )
    pmids = setdiff(pmids, invalid_pmids)
  }
  if (length(pmids) > 0) {
    groups = split(pmids, ceiling(seq_along(pmids) / max(1L, batch_size)))
    documents = list()
    for (group in groups) {
      url = paste0(
        "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson?pmids=",
        paste(group, collapse = ",")
      )
      request_count = request_count + 1L
      result = tryCatch({
        .plant_fetch_json(url, cache, file.path(cache_dir, "pubtator_bioc"),
                          effective_throttle, request_fun,
                          timeout = request_timeout)
      }, error = function(error) {
        if (.plant_pubtator_no_records_error(error, group)) {
          no_record_pmids <<- unique(c(no_record_pmids, group))
          .plant_pubtator_cache_no_records(
            url, group, cache, file.path(cache_dir, "pubtator_bioc")
          )
          return(NULL)
        }
        errors <<- errors + 1L
        error_messages <<- c(error_messages, conditionMessage(error))
        consecutive_busy <<- .plant_provider_busy_update(consecutive_busy,
                                                         error)
        NULL
      })
      if (is.null(result) && .plant_provider_circuit_open(consecutive_busy)) {
        circuit_open = TRUE
        break
      }
      if (!is.null(result)) consecutive_busy = 0L
      if (!is.null(result)) {
        cache_hit_count = cache_hit_count + as.integer(.plant_cache_hit(result))
        cached_no_records = tryCatch(
          unlist(result$uafR_no_records$pmids, use.names = FALSE),
          error = function(error) character()
        )
        no_record_pmids = unique(c(no_record_pmids,
                                   .uaf_non_empty(cached_no_records)))
      }
      documents = c(documents, .plant_pubtator_docs(result))
    }
    for (i in seq_len(nrow(plant_queries))) {
      species_pmids = unique(.uaf_non_empty(pubmed_literature$pmid[
        pubmed_literature$species == plant_queries$species[[i]]
      ]))
      docs = documents[vapply(documents, function(doc) {
        .uaf_first_non_empty_text(doc$pmid, doc$id, doc$sourceid) %in%
          species_pmids
      }, logical(1))]
      query_row = plant_queries[i, , drop = FALSE]
      parsed = .plant_pubtator_rows(
        query_row, docs,
        accepted_species = .plant_verified_species_aliases(
          plant_aliases, query_row
        )
      )
      rows[[length(rows) + 1]] = parsed$LiteratureCandidates
    }
  } else for (i in seq_len(nrow(plant_queries))) {
    query_row = plant_queries[i, , drop = FALSE]
    accepted_species = .plant_verified_species_aliases(
      plant_aliases, query_row
    )
    url = paste0("https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/",
                 "?text=", utils::URLencode(
                   paste(paste(accepted_species, collapse = " OR "),
                         "phytochemical metabolite"),
                   reserved = TRUE))
    request_count = request_count + 1
    result = tryCatch({
      .plant_fetch_json(url, cache, file.path(cache_dir, "pubtator"),
                        effective_throttle, request_fun,
                        timeout = request_timeout)
    }, error = function(error) {
      errors <<- errors + 1
      error_messages <<- c(error_messages, conditionMessage(error))
      consecutive_busy <<- .plant_provider_busy_update(consecutive_busy,
                                                       error)
      NULL
    })
    if (is.null(result) && .plant_provider_circuit_open(consecutive_busy)) {
      circuit_open = TRUE
      break
    }
    if (!is.null(result)) consecutive_busy = 0L
    if (!is.null(result)) {
      cache_hit_count = cache_hit_count + as.integer(.plant_cache_hit(result))
    }
    parsed = .plant_pubtator_rows(
      query_row, result, accepted_species = accepted_species
    )
    rows[[length(rows) + 1]] = parsed$LiteratureCandidates
  }
  literature = .plant_bind_tables(rows, .plant_literature_cols())
  occurrences = .plant_occurrences_from_pubtator(literature)
  status = .plant_provider_status(nrow(literature), errors)
  list(
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = literature,
    ProviderDiagnostics = .plant_provider_diagnostics(
      "pubtator", TRUE, TRUE, TRUE, request_count, cache_hit_count,
      nrow(literature), errors, 0, status,
      paste("PubTator BioC annotations for PubMed PMIDs (free-text search",
            "fallback when no PubMed candidates exist); co-mentions are not",
            "confirmed occurrence.", if (length(no_record_pmids) > 0L) {
              paste(length(no_record_pmids),
                    "PMID(s) had no indexed PubTator record.")
            } else "",
            .plant_provider_circuit_note(circuit_open)),
      error_messages = .plant_error_messages(error_messages)),
    Provenance = .plant_provenance("provider_dispatch", "pubtator",
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/",
                                   nrow(literature),
                                   "PubTator chemical/species candidate search.")
  )
}

.plant_pubtator_no_records_error = function(error, pmids) {
  pmids = .uaf_non_empty(pmids)
  length(pmids) > 0L && all(grepl("^[0-9]+$", pmids)) &&
    grepl("(?:HTTP(?: status)?(?: was)?[^0-9]*)?400(?: Bad Request)?",
          conditionMessage(error), ignore.case = TRUE, perl = TRUE)
}

.plant_pubtator_cache_no_records = function(url, pmids, cache, cache_dir) {
  if (!isTRUE(cache)) return(invisible(NULL))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file = file.path(cache_dir, paste0(.pubchem_url_hash(url), ".json"))
  payload = jsonlite::toJSON(
    list(
      PubTator3 = list(),
      uafR_no_records = list(
        reason = paste(
          "HTTP 400 for a numeric PMID export batch;",
          "no indexed PubTator record returned."
        ),
        pmids = as.list(as.character(pmids)),
        recorded_at = .plant_timestamp()
      )
    ),
    auto_unbox = TRUE, null = "null", na = "null"
  )
  .pubchem_atomic_write_text(payload, cache_file)
  invisible(cache_file)
}

.plant_provider_empty_result = function(provider, plant_queries,
                                        status = "no_records",
                                        message = "No provider records returned.") {
  list(
    PlantCompoundOccurrences = .plant_empty_occurrences(),
    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    SourceCompoundIdentity = .plant_empty_source_compound_identity(),
    ProviderDiagnostics = .plant_provider_diagnostics(
      provider = provider,
      enabled = TRUE,
      queried = FALSE,
      available = FALSE,
      request_count = 0,
      cache_hit_count = 0,
      record_count = 0,
      error_count = 0,
      warning_count = ifelse(status %in% c("not_implemented", "not_queried"),
                             1, 0),
      status = status,
      message = message,
      elapsed_seconds = 0,
      error_messages = NA_character_
    ),
    Provenance = .plant_provenance("provider_dispatch", provider,
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   NA_character_, 0, message)
  )
}

.plant_provider_status = function(record_count, errors) {
  if (errors > 0) return("warning")
  if (record_count > 0) return("ok")
  "no_records"
}

.plant_service_busy_message = function(message) {
  message = paste(.uaf_non_empty(message), collapse = " ")
  nzchar(message) && grepl(
    paste0(
      "(^|[^0-9])(408|425|429|500|502|503|504)([^0-9]|$)|",
      "rate[ -]?limit|service[ -]?busy|service unavailable|",
      "temporar(il)?y unavailable|search backend failed|cannot connect to solr|",
      "timed? out|timeout|connection reset|could not resolve host|",
      "failure when receiving data|empty reply from server"
    ),
    message,
    ignore.case = TRUE
  )
}

.plant_provider_busy_update = function(current, error) {
  if (.plant_service_busy_message(conditionMessage(error))) {
    as.integer(current) + 1L
  } else {
    0L
  }
}

.plant_provider_circuit_open = function(consecutive_busy, limit = 2L) {
  is.finite(consecutive_busy) && consecutive_busy >= limit
}

.plant_provider_circuit_note = function(circuit_open) {
  if (!isTRUE(circuit_open)) return("")
  paste(
    "The provider circuit breaker stopped this chunk after two consecutive",
    "rate-limit or service-busy errors; resume later with the existing cache."
  )
}

.plant_provider_diagnostics = function(provider, enabled, queried, available,
                                       request_count, cache_hit_count,
                                       record_count, error_count,
                                       warning_count, status, message,
                                       elapsed_seconds = NA_real_,
                                       error_messages = NA_character_,
                                       timeout_count = NA_integer_,
                                       rate_limit_count = NA_integer_,
                                       no_hit_reason = NA_character_) {
  timeout_count = .plant_provider_event_count(
    timeout_count, error_messages, "timeout"
  )
  rate_limit_count = .plant_provider_event_count(
    rate_limit_count, error_messages, "rate_limit"
  )
  no_hit_reason = .plant_provider_no_hit_reason(
    status, record_count, error_count, no_hit_reason
  )
  warning_message = if (tolower(.uaf_first_non_empty_text(status, "")) %in%
                        c("warning", "error", "failed", "timeout",
                          "timed_out", "rate_limited", "unavailable",
                          "not_queried", "not_implemented")) {
    .plant_error_messages(c(message, error_messages))
  } else NA_character_
  data.frame(
    provider = provider,
    enabled = .uaf_yes_no(enabled),
    queried = .uaf_yes_no(queried),
    available = .uaf_yes_no(available),
    request_count = suppressWarnings(as.integer(request_count)),
    cache_hit_count = suppressWarnings(as.integer(cache_hit_count)),
    record_count = suppressWarnings(as.integer(record_count)),
    error_count = suppressWarnings(as.integer(error_count)),
    warning_count = suppressWarnings(as.integer(warning_count)),
    timeout_count = suppressWarnings(as.integer(timeout_count)),
    rate_limit_count = suppressWarnings(as.integer(rate_limit_count)),
    status = status,
    no_hit_reason = no_hit_reason,
    warning_message = .plant_redact_secrets(warning_message),
    message = .plant_redact_secrets(message),
    retrieved_at = .plant_timestamp(),
    elapsed_seconds = suppressWarnings(as.numeric(elapsed_seconds)),
    error_messages = .uaf_first_non_empty_text(
      .plant_redact_secrets(error_messages)
    ),
    stringsAsFactors = FALSE
  )
}

.plant_provider_event_count = function(value, messages,
                                        type = c("timeout", "rate_limit")) {
  type = match.arg(type)
  numeric_value = suppressWarnings(as.integer(value[[1]]))
  if (is.finite(numeric_value)) return(numeric_value)
  messages = .uaf_non_empty(.plant_redact_secrets(messages))
  if (length(messages) < 1) return(0L)
  pattern = if (type == "timeout") {
    "time[ -]?out|timed out|operation.*timed"
  } else {
    "(^|[^0-9])(429|503)([^0-9]|$)|rate[ -]?limit|service[ -]?busy|service unavailable"
  }
  as.integer(sum(vapply(messages, grepl, logical(1), pattern = pattern,
                        ignore.case = TRUE, perl = TRUE)))
}

.plant_provider_no_hit_reason = function(status, record_count, error_count,
                                          supplied = NA_character_) {
  supplied = .uaf_first_non_empty_text(supplied)
  if (!is.na(supplied)) return(supplied)
  status = tolower(.uaf_first_non_empty_text(status, "unknown"))
  records = suppressWarnings(as.numeric(record_count[[1]]))
  errors = suppressWarnings(as.numeric(error_count[[1]]))
  if (is.finite(records) && records > 0) return(NA_character_)
  if (status %in% c("unavailable", "not_implemented")) {
    return("provider_unavailable")
  }
  if (status %in% c("not_queried", "disabled")) return("provider_not_queried")
  if ((is.finite(errors) && errors > 0) ||
      status %in% c("error", "failed", "warning", "timeout", "timed_out",
                    "rate_limited", "service_unavailable")) {
    return("provider_error_or_incomplete_query")
  }
  "no_records_returned_for_exact_query"
}

.plant_provider_finish = function(result, provider, started, progress) {
  elapsed = round(as.numeric(difftime(Sys.time(), started,
                                      units = "secs")), 3)
  if (!is.list(result)) {
    result = list(
      PlantCompoundOccurrences = .plant_empty_occurrences(),
      LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
      SourceCompoundIdentity = .plant_empty_source_compound_identity(),
      ProviderDiagnostics = .plant_provider_diagnostics(
        provider, TRUE, TRUE, FALSE, 0, 0, 0, 1, 0, "warning",
        "Provider adapter returned a non-list result.",
        error_messages = "Provider adapter returned a non-list result."
      ),
      Provenance = .plant_provenance("provider_dispatch", provider,
                                     NA_character_, NA_character_, 0,
                                     "Provider adapter returned a non-list result.")
    )
  }
  supplied_identity = if (is.data.frame(result$SourceCompoundIdentity)) {
    result$SourceCompoundIdentity
  } else {
    .plant_empty_source_compound_identity()
  }
  result$SourceCompoundIdentity = .plant_merge_source_identity_tables(list(
    supplied_identity,
    .plant_source_identity_from_occurrences(result$PlantCompoundOccurrences)
  ))
  if (is.data.frame(result$ProviderDiagnostics)) {
    result$ProviderDiagnostics$elapsed_seconds = elapsed
  }
  count = nrow(result$PlantCompoundOccurrences) +
    nrow(result$LiteratureCandidates)
  .plant_progress(progress, "uafR plant provider ", provider, ": ",
                  count, " record(s), ", elapsed, " sec")
  result
}

.plant_source_identity_from_occurrences = function(occurrences) {
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) return(.plant_empty_source_compound_identity())
  keep = !is.na(occurrences$compound_id) & occurrences$compound_id != ""
  occurrences = occurrences[keep, , drop = FALSE]
  if (nrow(occurrences) < 1) return(.plant_empty_source_compound_identity())
  cid_type = grepl("pubchem.*cid|^cid$",
                   tolower(.uaf_squish_text(occurrences$compound_id_type)))
  cid = rep(NA_integer_, nrow(occurrences))
  cid[cid_type] = suppressWarnings(as.integer(occurrences$compound_id[cid_type]))
  molecular_formula = rep(NA_character_, nrow(occurrences))
  knapsack = grepl("knapsack",
                   tolower(.uaf_squish_text(occurrences$source_database)))
  molecular_formula[knapsack] = vapply(
    occurrences$evidence_text[knapsack],
    .plant_knapsack_formula_from_evidence,
    character(1)
  )
  out = data.frame(
    compound_name = occurrences$compound_name,
    compound_name_clean = occurrences$compound_name_clean,
    source_database = occurrences$source_database,
    source_record_id = occurrences$source_record_id,
    source_compound_id = occurrences$compound_id,
    source_compound_id_type = occurrences$compound_id_type,
    CID = cid,
    InChIKey = NA_character_,
    SMILES = NA_character_,
    MolecularFormula = molecular_formula,
    evidence_url = occurrences$evidence_url,
    evidence_text = occurrences$evidence_text,
    identity_status = "source_identifier_only",
    identity_note = paste(
      "The provider supplied this identifier; no structure was inferred from",
      "the occurrence row. Cross-source identity validation is still required."
    ),
    stringsAsFactors = FALSE
  )
  .plant_bind_unique_tables(
    list(out), .plant_source_compound_identity_cols(),
    c("source_database", "source_record_id", "source_compound_id")
  )
}

.plant_knapsack_formula_from_evidence = function(text) {
  text = .uaf_first_non_empty_text(text)
  if (is.na(text)) return(NA_character_)
  match = regexec("(?:Molecular_formula|Formula)=([^;]+)", text,
                  ignore.case = TRUE, perl = TRUE)
  value = regmatches(text, match)[[1]]
  if (length(value) < 2) return(NA_character_)
  formula = .uaf_squish_text(value[[2]])
  if (!.plant_formula_like(formula)) NA_character_ else formula
}

.plant_merge_source_identity_tables = function(tables) {
  out = .plant_bind_tables(tables, .plant_source_compound_identity_cols())
  if (nrow(out) < 1) return(out)
  source_key = paste(
    .plant_clean_name(out$source_database),
    .uaf_squish_text(out$source_record_id),
    .uaf_squish_text(out$source_compound_id),
    sep = "\r"
  )
  cid_numeric = suppressWarnings(as.numeric(out$CID))
  has_specific_identity =
    (!is.na(cid_numeric) & is.finite(cid_numeric)) |
    (!is.na(out$InChIKey) & out$InChIKey != "") |
    (!is.na(out$SMILES) & out$SMILES != "") |
    (!is.na(out$MolecularFormula) & out$MolecularFormula != "")
  specific_keys = unique(source_key[has_specific_identity])
  drop_generic = !has_specific_identity & source_key %in% specific_keys
  out = out[!drop_generic, , drop = FALSE]
  key = paste(source_key[!drop_generic],
              .uaf_squish_text(out$InChIKey),
              .uaf_squish_text(out$SMILES), sep = "\r")
  out = out[!duplicated(key), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_error_messages = function(messages) {
  messages = unique(.uaf_non_empty(.plant_redact_secrets(messages)))
  if (length(messages) < 1) return(NA_character_)
  .plant_truncate(paste(messages, collapse = " | "), 1000)
}

.plant_redact_secrets = function(x) {
  x = as.character(x)
  gsub(
    paste0(
      "(?i)(\\b(?:api[ _-]?key|access[ _-]?token|password|secret)",
      "\\s*[:=]\\s*)[^,;&\\s\"'<>]+"
    ),
    "\\1[REDACTED]",
    x,
    perl = TRUE
  )
}

.plant_progress = function(progress, ...) {
  if (isTRUE(progress)) message(...)
  invisible(NULL)
}

.plant_normalize_occurrences = function(x, source_hint = NA_character_) {
  cols = .plant_occurrence_cols()
  if (!is.data.frame(x) || nrow(x) < 1) return(.uaf_empty_table(cols))
  names(x) = .plant_normalize_column_names(names(x))
  for (col in cols) if (!col %in% names(x)) x[[col]] = NA_character_
  x = x[, cols, drop = FALSE]
  x$query_plant = .uaf_squish_text(x$query_plant)
  x$species = .uaf_squish_text(x$species)
  missing_query = is.na(x$query_plant) | x$query_plant == ""
  x$query_plant[missing_query] = x$species[missing_query]
  x$query_plant_clean = .plant_clean_name(
    ifelse(is.na(x$query_plant_clean) | x$query_plant_clean == "",
           x$query_plant, x$query_plant_clean)
  )
  x$genus = ifelse(is.na(x$genus) | x$genus == "", .plant_genus(x$species),
                   .uaf_squish_text(x$genus))
  x$family = .uaf_squish_text(x$family)
  x$matched_taxon = .uaf_squish_text(ifelse(is.na(x$matched_taxon) |
                                              x$matched_taxon == "",
                                            x$species, x$matched_taxon))
  x$matched_rank = tolower(.uaf_squish_text(ifelse(is.na(x$matched_rank) |
                                                     x$matched_rank == "",
                                                   "species", x$matched_rank)))
  x$matched_rank[!x$matched_rank %in% c("species", "genus", "family")] = "unknown"
  x$compound_name = .uaf_squish_text(x$compound_name)
  x$compound_name_clean = .plant_clean_compound(
    ifelse(is.na(x$compound_name) | x$compound_name == "",
           x$compound_name_clean, x$compound_name)
  )
  x$source_database = .uaf_squish_text(ifelse(is.na(x$source_database) |
                                                x$source_database == "",
                                              source_hint, x$source_database))
  x$source_database[is.na(x$source_database) | x$source_database == ""] =
    "unknown"
  x$retrieved_at = ifelse(is.na(x$retrieved_at) | x$retrieved_at == "",
                          .plant_timestamp(), x$retrieved_at)
  x$plant_part = .plant_clean_context_value(x$plant_part)
  x$tissue = .plant_clean_context_value(x$tissue)
  x$method = .plant_clean_context_value(x$method)
  x$evidence_tier = .plant_normalize_evidence_tier(x$evidence_tier,
                                                   x$source_database,
                                                   x$matched_rank)
  x$confidence = .plant_normalize_confidence(x$confidence, x$evidence_tier)
  x$curation_flag = .uaf_squish_text(ifelse(is.na(x$curation_flag) |
                                              x$curation_flag == "",
                                            "unreviewed", x$curation_flag))
  x$plant_part_group = .plant_context_group(x$plant_part,
                                            x$plant_part_group,
                                            "plant_part")
  x$tissue_group = .plant_context_group(x$tissue, x$tissue_group,
                                        "tissue")
  x$method_group = .plant_method_group(x$method, x$method_group,
                                       x$source_database, x$evidence_tier)
  invalid_part = x$plant_part_group %in%
    c("unknown", "other", "extract_unspecified")
  invalid_tissue = x$tissue_group %in%
    c("unknown", "other", "extract_unspecified")
  method_like_part = .plant_context_raw_method_like(x$plant_part)
  method_like_tissue = .plant_context_raw_method_like(x$tissue)
  x$plant_part[method_like_part & !invalid_part] =
    x$plant_part_group[method_like_part & !invalid_part]
  x$tissue[method_like_tissue & !invalid_tissue] =
    x$tissue_group[method_like_tissue & !invalid_tissue]
  x$plant_part[invalid_part] = NA_character_
  x$tissue[invalid_tissue] = NA_character_
  x$method[x$method_group %in%
             c("unknown", "other", "database_record",
               "literature_curation")] = NA_character_
  x$biological_context_status = .plant_context_status(x$plant_part_group,
                                                      x$tissue_group,
                                                      x$method_group)
  x$occurrence_status = .plant_occurrence_status(x$evidence_tier,
                                                 x$matched_rank,
                                                 x$source_database,
                                                 x$curation_flag)
  x$occurrence_basis = .plant_occurrence_basis(x$evidence_tier,
                                               x$matched_rank,
                                               x$source_database)
  x$evidence_quality_score = .plant_evidence_quality_score(
    x$occurrence_status, x$occurrence_basis, x$confidence,
    x$biological_context_status
  )
  x$analysis_ready = .plant_analysis_ready(x$occurrence_status,
                                           x$confidence,
                                           x$evidence_quality_score)
  x = x[!is.na(x$species) & x$species != "" &
          !is.na(x$compound_name) & x$compound_name != "", ,
        drop = FALSE]
  x = unique(x)
  row.names(x) = NULL
  x
}

.plant_normalize_literature = function(x, source_hint = "pubmed") {
  cols = .plant_literature_cols()
  if (!is.data.frame(x) || nrow(x) < 1) return(.uaf_empty_table(cols))
  names(x) = .plant_normalize_column_names(names(x))
  for (col in cols) if (!col %in% names(x)) x[[col]] = NA_character_
  x = x[, cols, drop = FALSE]
  x$query_plant = .uaf_squish_text(x$query_plant)
  x$query_plant_clean = .plant_clean_name(x$query_plant)
  x$species = .uaf_squish_text(x$species)
  x$genus = ifelse(is.na(x$genus) | x$genus == "", .plant_genus(x$species),
                   .uaf_squish_text(x$genus))
  x$source_database = .uaf_squish_text(ifelse(is.na(x$source_database) |
                                                x$source_database == "",
                                              source_hint, x$source_database))
  x$retrieved_at = ifelse(is.na(x$retrieved_at) | x$retrieved_at == "",
                          .plant_timestamp(), x$retrieved_at)
  x$evidence_tier = .plant_normalize_evidence_tier(x$evidence_tier,
                                                   x$source_database,
                                                   "species")
  x$confidence = .plant_normalize_confidence(x$confidence, x$evidence_tier)
  x$curation_flag = .uaf_squish_text(ifelse(is.na(x$curation_flag) |
                                              x$curation_flag == "",
                                            "candidate", x$curation_flag))
  row.names(x) = NULL
  x
}

.plant_match_occurrences_to_queries = function(occurrences, plant_queries,
                                                normalize = TRUE) {
  if (isTRUE(normalize)) {
    occurrences = .plant_normalize_occurrences(occurrences)
  }
  if (nrow(occurrences) < 1) return(occurrences)
  lookup = plant_queries[, c("query_plant", "query_plant_clean", "species",
                             "genus", "family"), drop = FALSE]
  match_index = .plant_query_lookup_match(
    lookup, .plant_clean_name(occurrences$species)
  )
  matched = !is.na(match_index)
  if (!any(matched)) return(occurrences)
  occurrences$query_plant[matched] = lookup$query_plant[match_index[matched]]
  occurrences$query_plant_clean[matched] =
    lookup$query_plant_clean[match_index[matched]]
  occurrences$family = .plant_fill_matched_values(
    occurrences$family, lookup$family, match_index
  )
  occurrences$genus = .plant_fill_matched_values(
    occurrences$genus, lookup$genus, match_index
  )
  occurrences
}

.plant_match_literature_to_queries = function(literature, plant_queries,
                                               normalize = TRUE) {
  if (isTRUE(normalize)) {
    literature = .plant_normalize_literature(literature)
  }
  if (nrow(literature) < 1) return(literature)
  lookup = plant_queries[, c("query_plant", "query_plant_clean", "species",
                             "genus", "family"), drop = FALSE]
  match_index = .plant_query_lookup_match(
    lookup, .plant_clean_name(literature$species)
  )
  matched = !is.na(match_index)
  if (!any(matched)) return(literature)
  literature$query_plant[matched] = lookup$query_plant[match_index[matched]]
  literature$query_plant_clean[matched] =
    lookup$query_plant_clean[match_index[matched]]
  literature$family = .plant_fill_matched_values(
    literature$family, lookup$family, match_index
  )
  literature$genus = .plant_fill_matched_values(
    literature$genus, lookup$genus, match_index
  )
  literature
}

.plant_query_lookup_match = function(lookup, keys) {
  if (!is.data.frame(lookup) || nrow(lookup) < 1L || length(keys) < 1L) {
    return(rep(NA_integer_, length(keys)))
  }
  lookup_query = .plant_clean_name(lookup$query_plant_clean)
  lookup_species = .plant_clean_name(lookup$species)
  lookup_rows = seq_len(nrow(lookup))
  map_key = c(lookup_query, lookup_species)
  map_row = rep(lookup_rows, 2L)
  valid = !is.na(map_key) & map_key != ""
  map_key = map_key[valid]
  map_row = map_row[valid]
  ordered = order(map_row)
  map_key = map_key[ordered]
  map_row = map_row[ordered]
  keep = !duplicated(map_key)
  map_key = map_key[keep]
  map_row = map_row[keep]
  map_row[match(keys, map_key, nomatch = NA_integer_)]
}

.plant_fill_matched_values = function(current, lookup_values, match_index) {
  current = .uaf_squish_text(current)
  lookup_values = .uaf_squish_text(lookup_values)
  matched = !is.na(match_index)
  missing = is.na(current) | current == ""
  fill = matched & missing
  if (any(fill)) {
    replacement = lookup_values[match_index[fill]]
    available = !is.na(replacement) & replacement != ""
    current[which(fill)[available]] = replacement[available]
  }
  current
}

.plant_bind_occurrences = function(parts) {
  .plant_collapse_occurrence_evidence(
    .plant_normalize_occurrences(.plant_bind_tables(parts, .plant_occurrence_cols()))
  )
}

.plant_collapse_occurrence_evidence = function(occurrences) {
  if (!is.data.frame(occurrences) || nrow(occurrences) < 2) {
    return(occurrences)
  }
  dup_cols = .plant_duplicate_key_cols("PlantCompoundOccurrences")
  if (!all(dup_cols %in% names(occurrences))) return(occurrences)
  key_data = occurrences[dup_cols]
  key_data[] = lapply(key_data, function(value) {
    value = .uaf_squish_text(value)
    value[is.na(value)] = ""
    value
  })
  key = do.call(paste, c(key_data, sep = "||"))
  duplicate = duplicated(key) | duplicated(key, fromLast = TRUE)
  if (!any(duplicate)) {
    out = occurrences[order(key), , drop = FALSE]
    row.names(out) = NULL
    return(out)
  }
  groups = split(which(duplicate), key[duplicate])
  rows = lapply(groups, function(idx) {
    .plant_merge_occurrence_group(occurrences[idx, , drop = FALSE])
  })
  merged = .plant_normalize_occurrences(do.call(rbind, rows))
  unique_idx = which(!duplicate)
  out = rbind(occurrences[unique_idx, , drop = FALSE], merged)
  out_key = c(key[unique_idx], names(groups))
  out = out[order(out_key), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_merge_occurrence_group = function(group) {
  scores = suppressWarnings(as.numeric(group$evidence_quality_score))
  scores[!is.finite(scores)] = 0
  rank_score = .plant_occurrence_rank_score(group$occurrence_status) +
    .plant_confidence_score(group$confidence) + scores
  best = which.max(rank_score)
  row = group[best, , drop = FALSE]
  for (col in names(row)) {
    if (col %in% c("evidence_text", "curation_flag")) next
    if (length(.uaf_non_empty(row[[col]])) > 0) next
    values = .uaf_non_empty(group[[col]])
    if (length(values) > 0) {
      row[[col]] = values[[1]]
      if (identical(col, "plant_part") && "plant_part_group" %in% names(row)) {
        row$plant_part_group = NA_character_
      }
      if (identical(col, "tissue") && "tissue_group" %in% names(row)) {
        row$tissue_group = NA_character_
      }
      if (identical(col, "method") && "method_group" %in% names(row)) {
        row$method_group = NA_character_
      }
    }
  }
  row$evidence_text = .plant_truncate(
    .pubchem_collapse(unique(.uaf_non_empty(group$evidence_text))),
    1200
  )
  row$curation_flag = .pubchem_collapse(unique(.uaf_non_empty(group$curation_flag)))
  row
}

.plant_occurrence_rank_score = function(occurrence_status) {
  status = .plant_matrix_key(occurrence_status)
  out = rep(0, length(status))
  out[status == "curated_reported"] = 5
  out[status == "direct_reported"] = 4
  out[status == "literature_reported"] = 3
  out[status == "taxon_fallback"] = 2
  out[status == "candidate"] = 1
  out
}

.plant_bind_tables = function(parts, cols) {
  parts = parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) < 1) return(.uaf_empty_table(cols))
  normalized = lapply(parts, function(x) {
    if (!is.data.frame(x)) return(.uaf_empty_table(cols))
    for (col in cols) if (!col %in% names(x)) x[[col]] = NA_character_
    x[, cols, drop = FALSE]
  })
  out = do.call(rbind, normalized)
  if (is.null(out) || nrow(out) < 1) return(.uaf_empty_table(cols))
  row.names(out) = NULL
  out
}

.plant_compound_resolution = function(occurrences, categorate_result) {
  occurrences = .plant_normalize_occurrences(occurrences)
  compounds = sort(unique(.uaf_non_empty(occurrences$compound_name)))
  if (length(compounds) < 1) {
    return(.uaf_empty_table(.plant_compound_resolution_cols()))
  }
  props = .plant_categorate_properties(categorate_result)
  resolved_source = .plant_enrichment_source(categorate_result)
  unresolved_note = if (identical(resolved_source, "pubchemProfile_identity")) {
    "No PubChem identity hit."
  } else {
    "No compound enrichment hit."
  }
  compound_keys = .plant_clean_compound(compounds)
  occurrence_keys = occurrences$compound_name_clean
  occurrence_keys[is.na(occurrence_keys)] = ""
  occurrence_counts = table(occurrence_keys)
  query_count = as.integer(occurrence_counts[compound_keys])
  query_count[is.na(query_count)] = 0L

  props_keys = .plant_clean_compound(props$Query)
  props_keys[is.na(props_keys)] = ""
  hit_keys = unique(props_keys)
  has_hit = compound_keys %in% hit_keys
  first_non_empty_by_key = function(values) {
    values = .uaf_squish_text(values)
    keep = !is.na(values) & values != "" & props_keys != ""
    candidates = which(keep)
    first = candidates[!duplicated(props_keys[candidates])]
    if (length(first) < 1L) return(rep(NA_character_, length(compound_keys)))
    values[first][match(compound_keys, props_keys[first])]
  }
  cid = first_non_empty_by_key(props$CID)
  inchikey = first_non_empty_by_key(props$InChIKey)
  molecular_formula = first_non_empty_by_key(props$MolecularFormula)
  isomeric_smiles = first_non_empty_by_key(props$IsomericSMILES)
  canonical_smiles = first_non_empty_by_key(props$CanonicalSMILES)
  general_smiles = first_non_empty_by_key(props$SMILES)
  connectivity_smiles = first_non_empty_by_key(props$ConnectivitySMILES)
  smiles = isomeric_smiles
  fill = is.na(smiles) | smiles == ""
  smiles[fill] = canonical_smiles[fill]
  fill = is.na(smiles) | smiles == ""
  smiles[fill] = general_smiles[fill]
  fill = is.na(smiles) | smiles == ""
  smiles[fill] = connectivity_smiles[fill]
  resolved = has_hit & (!is.na(cid) | !is.na(inchikey) | !is.na(smiles))

  out = data.frame(
    compound_name = compounds,
    compound_name_clean = compound_keys,
    query_count = query_count,
    resolved = resolved,
    CID = cid,
    InChIKey = inchikey,
    SMILES = smiles,
    MolecularFormula = molecular_formula,
    resolution_source = ifelse(has_hit, resolved_source, "unresolved"),
    notes = ifelse(has_hit, "", unresolved_note),
    stringsAsFactors = FALSE
  )
  .plant_collapse_compound_resolution_keys(out)
}

.plant_collapse_compound_resolution_keys = function(resolution) {
  resolution = .plant_bind_tables(list(resolution),
                                  .plant_compound_resolution_cols())
  if (nrow(resolution) < 2) return(resolution)
  key = resolution$compound_name_clean
  key[is.na(key)] = ""
  duplicate_key = duplicated(key) | duplicated(key, fromLast = TRUE)
  if (!any(duplicate_key)) {
    out = resolution[order(key), , drop = FALSE]
    row.names(out) = NULL
    return(out)
  }
  unique_rows = resolution[!duplicate_key, , drop = FALSE]
  duplicate_rows = which(duplicate_key)
  groups = split(duplicate_rows, key[duplicate_rows])
  rows = lapply(groups, function(idx) {
    group = resolution[idx, , drop = FALSE]
    score = as.numeric(group$resolved %in% TRUE) * 100
    score = score + as.numeric(!is.na(group$SMILES) & group$SMILES != "") * 20
    score = score + as.numeric(!is.na(group$InChIKey) &
                                 group$InChIKey != "") * 10
    score = score + as.numeric(!is.na(group$CID) & group$CID != "") * 5
    score = score + as.numeric(!is.na(group$MolecularFormula) &
                                 group$MolecularFormula != "") * 2
    score = score - (nchar(group$compound_name) / 10000)
    best = which.max(score)
    row = group[best, , drop = FALSE]
    counts = suppressWarnings(as.integer(group$query_count))
    row$query_count = ifelse(any(is.finite(counts)), max(counts, na.rm = TRUE),
                             row$query_count)
    cids = unique(.uaf_non_empty(group$CID))
    aliases = setdiff(unique(.uaf_non_empty(group$compound_name)),
                      row$compound_name)
    notes = unique(.uaf_non_empty(c(
      row$notes,
      if (length(aliases) > 0) {
        paste0("Collapsed aliases: ", .pubchem_collapse(aliases))
      },
      if (length(cids) > 1) {
        paste0("Multiple PubChem CIDs observed for this normalized key: ",
               .pubchem_collapse(cids))
      }
    )))
    row$notes = .plant_truncate(.pubchem_collapse(notes), 1200)
    row
  })
  out = rbind(unique_rows, do.call(rbind, rows))
  out = out[order(out$compound_name_clean), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_occurrences_from_input = function(x) {
  if ((inherits(x, "uaf_plant_phytochemistry") ||
       (is.list(x) && !is.data.frame(x))) &&
      is.data.frame(x$PlantCompoundOccurrences)) {
    return(.plant_normalize_occurrences(x$PlantCompoundOccurrences))
  }
  if (is.data.frame(x) && all(.plant_occurrence_cols() %in% names(x))) {
    return(.plant_normalize_occurrences(x))
  }
  standardizePlantCompoundIntake(x)
}

.plant_resolve_compound_identities = function(occurrences, cache, cache_dir,
                                              throttle, batch_size, resume,
                                              progress, pubchem_fun,
                                              request_fun,
                                              lotus_index = NULL,
                                              source_identity = NULL,
                                              source_only = FALSE) {
  occurrences = .plant_normalize_occurrences(occurrences)
  compounds = sort(unique(.uaf_non_empty(occurrences$compound_name)))
  source_identity = .plant_merge_source_identity_tables(list(
    source_identity,
    .plant_source_compound_identity(occurrences, lotus_index)
  ))
  source_identity = .plant_annotate_source_identity_status(source_identity)
  source_resolution = .plant_source_compound_resolution(occurrences,
                                                        source_identity)
  if (length(compounds) < 1) {
    categorate_result = .plant_empty_identity_categorate()
    categorate_result$SourceCompoundIdentity = source_identity
    return(list(
      CategorateResult = categorate_result,
      CompoundResolution = .plant_compound_resolution(occurrences,
                                                      categorate_result),
      Provenance = .plant_provenance("compound_identity", "not_run",
                                     NA_character_, NA_character_, 0,
                                     "No compounds were available for identity resolution.")
    ))
  }
  source_resolved = source_resolution$compound_name_clean[
    source_resolution$resolved %in% TRUE
  ]
  pubchem_compounds = if (isTRUE(source_only)) {
    character()
  } else {
    compounds[!.plant_clean_compound(compounds) %in% source_resolved]
  }
  categorate_result = if (length(pubchem_compounds) > 0) {
    .plant_pubchem_identity_categorate(
      compounds = pubchem_compounds,
      cache = cache,
      cache_dir = cache_dir,
      throttle = throttle,
      batch_size = batch_size,
      resume = resume,
      progress = progress,
      pubchem_fun = pubchem_fun,
      request_fun = request_fun
    )
  } else if (isTRUE(source_only)) {
    .plant_empty_identity_categorate("source_identity_only")
  } else {
    .plant_empty_identity_categorate("source_identity_only")
  }
  categorate_result$SourceCompoundIdentity = source_identity
  pubchem_resolution = .plant_compound_resolution(occurrences,
                                                  categorate_result)
  if (nrow(source_identity) > 0) {
    categorate_result$EnrichmentMode = if (length(pubchem_compounds) > 0) {
      "source_pubchem_identity"
    } else {
      "source_identity_only"
    }
  }
  source_resolved_resolution = source_resolution[
    source_resolution$resolved %in% TRUE, , drop = FALSE
  ]
  source_review_resolution = source_resolution[
    !(source_resolution$resolved %in% TRUE), , drop = FALSE
  ]
  compound_resolution = .plant_merge_compound_resolution(
    pubchem_resolution, source_resolved_resolution
  )
  unresolved_after_merge = compound_resolution$compound_name_clean[
    !(compound_resolution$resolved %in% TRUE)
  ]
  source_review_resolution = source_review_resolution[
    source_review_resolution$compound_name_clean %in% unresolved_after_merge, ,
    drop = FALSE
  ]
  compound_resolution = .plant_merge_compound_resolution(
    compound_resolution, source_review_resolution
  )
  compound_resolution = .plant_add_source_identity_review_notes(
    compound_resolution, source_resolution
  )
  list(
    CategorateResult = categorate_result,
    CompoundResolution = compound_resolution,
    Provenance = .plant_provenance(
      "compound_identity", .plant_enrichment_source(categorate_result),
      paste(compounds, collapse = "; "), NA_character_, length(compounds),
      .plant_enrichment_note(categorate_result)
    )
  )
}

.plant_source_identity_from_input = function(x) {
  supplied = if (is.list(x) && !is.data.frame(x) &&
                 is.data.frame(x$SourceCompoundIdentity)) {
    x$SourceCompoundIdentity
  } else {
    attr(x, "SourceCompoundIdentity", exact = TRUE)
  }
  if (!is.data.frame(supplied)) return(.plant_empty_source_compound_identity())
  .plant_merge_source_identity_tables(list(supplied))
}

.plant_source_compound_identity = function(occurrences, lotus_index = NULL) {
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1 || is.null(.plant_lotus_index_or_null(lotus_index))) {
    return(.plant_empty_source_compound_identity())
  }
  lotus_occurrences = occurrences[
    tolower(.uaf_squish_text(occurrences$source_database)) == "lotus", ,
    drop = FALSE
  ]
  if (nrow(lotus_occurrences) < 1) {
    return(.plant_empty_source_compound_identity())
  }
  index = .plant_lotus_identity_index_for_occurrences(lotus_occurrences,
                                                      lotus_index)
  if (nrow(index) < 1) return(.plant_empty_source_compound_identity())
  occurrence_ids = unique(.uaf_non_empty(c(lotus_occurrences$source_record_id,
                                           lotus_occurrences$compound_id)))
  index_ids = .plant_lotus_identity_record_ids(index)
  keep = if (length(occurrence_ids) > 0) {
    index_ids %in% occurrence_ids
  } else {
    index$compound_name_clean %in% lotus_occurrences$compound_name_clean
  }
  index = index[keep, , drop = FALSE]
  if (nrow(index) < 1) return(.plant_empty_source_compound_identity())
  rows = lapply(seq_len(nrow(index)), function(i) {
    hit = index[i, , drop = FALSE]
    data.frame(
      compound_name = hit$compound_name,
      compound_name_clean = hit$compound_name_clean,
      source_database = "LOTUS",
      source_record_id = .plant_lotus_index_source_record_id(hit),
      source_compound_id = .plant_lotus_index_compound_id(hit),
      source_compound_id_type = .plant_lotus_index_compound_id_type(hit),
      CID = suppressWarnings(as.integer(.uaf_first_non_empty_text(hit$cid))),
      InChIKey = .uaf_first_non_empty_text(hit$inchikey),
      SMILES = .uaf_first_non_empty_text(hit$smiles),
      MolecularFormula = .uaf_first_non_empty_text(hit$molecular_formula),
      evidence_url = .plant_lotus_index_evidence_url(hit),
      evidence_text = .uaf_first_non_empty_text(
        hit$evidence_text,
        .plant_lotus_index_evidence_text(hit)
      ),
      identity_status = NA_character_,
      identity_note = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out = .plant_bind_tables(rows, .plant_source_compound_identity_cols())
  if (nrow(out) < 1) return(out)
  out = unique(out)
  out = .plant_annotate_source_identity_status(out)
  row.names(out) = NULL
  out
}

.plant_lotus_identity_index_for_occurrences = function(occurrences,
                                                       lotus_index) {
  lookup = .plant_lotus_lookup_info(lotus_index)
  if (!is.null(lookup)) {
    plant_queries = .plant_queries_from_occurrences(occurrences)
    return(.plant_lotus_lookup_index_rows(plant_queries, lookup))
  }
  index = tryCatch(
    standardizeLotusIndex(lotus_index),
    error = function(error) .plant_empty_lotus_index()
  )
  if (nrow(index) < 1) return(index)
  occurrence_ids = unique(.uaf_non_empty(c(occurrences$source_record_id,
                                           occurrences$compound_id)))
  index_ids = .plant_lotus_identity_record_ids(index)
  keep = if (length(occurrence_ids) > 0) {
    index_ids %in% occurrence_ids
  } else {
    index$compound_name_clean %in% occurrences$compound_name_clean
  }
  index[keep, , drop = FALSE]
}

.plant_lotus_identity_record_ids = function(index) {
  ids = mapply(
    function(source_record_id, lotus_id, wikidata_id, cid) {
      .uaf_first_non_empty_text(source_record_id, lotus_id, wikidata_id,
                                ifelse(!is.na(cid) & cid != "",
                                       paste0("cid_", cid), NA_character_),
                                cid)
    },
    source_record_id = index$source_record_id,
    lotus_id = index$lotus_id,
    wikidata_id = index$wikidata_id,
    cid = index$cid,
    USE.NAMES = FALSE
  )
  .uaf_squish_text(ids)
}

.plant_queries_from_occurrences = function(occurrences) {
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) return(.plant_queries(character(), "species"))
  x = unique(occurrences[, intersect(c("query_plant", "species", "genus",
                                       "family"), names(occurrences)),
                         drop = FALSE])
  for (col in c("query_plant", "species", "genus", "family")) {
    if (!col %in% names(x)) x[[col]] = NA_character_
  }
  plant_names = ifelse(!is.na(.uaf_squish_text(x$species)) &
                         .uaf_squish_text(x$species) != "",
                       .uaf_squish_text(x$species),
                       .uaf_squish_text(x$query_plant))
  plant_names = unique(.uaf_non_empty(plant_names))
  out = .plant_queries(plant_names, taxon_fallback = c("species", "genus",
                                                       "family"))
  if (nrow(out) < 1) return(out)
  for (i in seq_len(nrow(out))) {
    hit = x[.plant_clean_name(x$species) == out$query_plant_clean |
              .plant_clean_name(x$query_plant) == out$query_plant_clean, ,
            drop = FALSE]
    if (nrow(hit) < 1) next
    out$genus[[i]] = .uaf_first_non_empty_text(out$genus[[i]], hit$genus)
    out$family[[i]] = .uaf_first_non_empty_text(out$family[[i]], hit$family)
  }
  out
}

.plant_annotate_source_identity_status = function(identity) {
  identity = .plant_bind_tables(list(identity),
                                .plant_source_compound_identity_cols())
  if (nrow(identity) < 1) return(identity)
  groups = split(seq_len(nrow(identity)), identity$compound_name_clean)
  for (idx in groups) {
    group = identity[idx, , drop = FALSE]
    structures = unique(.uaf_non_empty(group$SMILES))
    inchikeys = unique(.uaf_non_empty(group$InChIKey))
    cids = unique(.uaf_non_empty(as.character(group$CID)))
    formulas = unique(.uaf_non_empty(group$MolecularFormula))
    sources = unique(.uaf_non_empty(group$source_database))
    source_label = .pubchem_collapse(sources)
    conflicting = length(structures) > 1 || length(inchikeys) > 1 ||
      length(cids) > 1
    if (!conflicting && (length(structures) == 1 || length(inchikeys) == 1)) {
      identity$identity_status[idx] = "source_structure_unique"
      identity$identity_note[idx] =
        paste0("Resolved from one non-conflicting source-backed structure for ",
               "this normalized compound key (source: ", source_label, ").")
    } else if (conflicting) {
      identity$identity_status[idx] = "source_structure_ambiguous"
      identity$identity_note[idx] =
        paste0("Conflicting source structures, InChIKeys, or CIDs map to this ",
               "normalized compound key; review source records before assigning ",
               "one identity (sources: ", source_label, ").")
    } else if (length(formulas) == 1) {
      identity$identity_status[idx] = "source_formula_only"
      identity$identity_note[idx] =
        paste0("A source provides one formula but no unique structure; keep ",
               "unresolved for structure-required workflows (source: ",
               source_label, ").")
    } else {
      identity$identity_status[idx] = "source_identity_incomplete"
      identity$identity_note[idx] =
        paste0("The source record did not provide a usable structure or ",
               "formula (source: ", source_label, ").")
    }
  }
  identity
}

.plant_source_compound_resolution = function(occurrences, source_identity) {
  occurrences = .plant_normalize_occurrences(occurrences)
  source_identity = .plant_bind_tables(list(source_identity),
                                       .plant_source_compound_identity_cols())
  source_identity = .plant_annotate_source_identity_status(source_identity)
  if (nrow(occurrences) < 1 || nrow(source_identity) < 1) {
    return(.uaf_empty_table(.plant_compound_resolution_cols()))
  }
  keys = sort(unique(.uaf_non_empty(source_identity$compound_name_clean)))
  rows = lapply(keys, function(key) {
    hits = source_identity[source_identity$compound_name_clean == key, ,
                           drop = FALSE]
    unique_status = any(hits$identity_status == "source_structure_unique")
    hit = hits[order(!(hits$identity_status == "source_structure_unique"),
                     is.na(hits$SMILES) | hits$SMILES == "",
                     is.na(hits$InChIKey) | hits$InChIKey == ""), ,
               drop = FALSE][1, , drop = FALSE]
    data.frame(
      compound_name = .uaf_first_non_empty_text(hit$compound_name),
      compound_name_clean = key,
      query_count = sum(occurrences$compound_name_clean == key),
      resolved = isTRUE(unique_status),
      CID = if (isTRUE(unique_status)) {
        .uaf_first_non_empty_text(hit$CID)
      } else {
        NA_character_
      },
      InChIKey = if (isTRUE(unique_status)) {
        .uaf_first_non_empty_text(hit$InChIKey)
      } else {
        NA_character_
      },
      SMILES = if (isTRUE(unique_status)) {
        .uaf_first_non_empty_text(hit$SMILES)
      } else {
        NA_character_
      },
      MolecularFormula = if (isTRUE(unique_status)) {
        .uaf_first_non_empty_text(hit$MolecularFormula)
      } else {
        NA_character_
      },
      resolution_source = if (isTRUE(unique_status)) {
        paste0(
          paste(sort(unique(.uaf_non_empty(hits$source_database))),
                collapse = "+"),
          "_source_identity"
        )
      } else {
        "source_identity_review_required"
      },
      notes = .plant_source_identity_resolution_note(hits),
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  .plant_collapse_compound_resolution_keys(out)
}

.plant_source_identity_resolution_note = function(hits) {
  statuses = unique(.uaf_non_empty(hits$identity_status))
  source_ids = unique(.uaf_non_empty(hits$source_record_id))
  notes = unique(.uaf_non_empty(c(
    hits$identity_note,
    if (length(source_ids) > 0) {
      paste0("Source records: ", .pubchem_collapse(utils::head(source_ids, 20)))
    },
    if (length(source_ids) > 20) {
      paste0("Additional source records omitted from note: ",
             length(source_ids) - 20)
    },
    if (length(statuses) > 0) {
      paste0("Source identity statuses: ", .pubchem_collapse(statuses))
    }
  )))
  .plant_truncate(.pubchem_collapse(notes), 1200)
}

.plant_add_source_identity_review_notes = function(resolution,
                                                   source_resolution) {
  resolution = .plant_bind_tables(list(resolution),
                                  .plant_compound_resolution_cols())
  source_resolution = .plant_bind_tables(list(source_resolution),
                                         .plant_compound_resolution_cols())
  if (nrow(resolution) < 1 || nrow(source_resolution) < 1) {
    return(resolution)
  }
  review = source_resolution[!(source_resolution$resolved %in% TRUE), ,
                             drop = FALSE]
  if (nrow(review) < 1) return(resolution)
  for (i in seq_len(nrow(resolution))) {
    hit = review[review$compound_name_clean ==
                   resolution$compound_name_clean[[i]], , drop = FALSE]
    if (nrow(hit) < 1) next
    notes = unique(.uaf_non_empty(c(resolution$notes[[i]], hit$notes)))
    resolution$notes[[i]] = .plant_truncate(.pubchem_collapse(notes), 1200)
  }
  resolution
}

.plant_identity_issue = function(row) {
  name = .uaf_first_non_empty_text(row$compound_name)
  clean = .uaf_first_non_empty_text(row$compound_name_clean,
                                    .plant_clean_compound(name))
  text = tolower(.uaf_first_non_empty_text(name, clean))
  notes = tolower(.uaf_first_non_empty_text(row$notes))
  resolved = row$resolved %in% TRUE
  if (.plant_identity_broad_class_label(text)) {
    return(.plant_identity_issue_row(
      if (resolved) "resolved_broad_class_label" else
        "broad_class_or_family_label",
      if (resolved) "review_label_before_matrix_use" else
        "exclude_from_structure_required_exports",
      NA_character_, "curated class metadata",
      if (resolved) {
        "A source-backed structure is available, but the displayed value is a class or family label; review the source record and replace the display name or keep as class metadata only."
      } else {
        "The value appears to describe a chemical class or family rather than one discrete compound."
      }
    ))
  }
  if (.plant_identity_source_label(text)) {
    return(.plant_identity_issue_row(
      if (resolved) "resolved_source_or_product_label" else
        "source_or_product_label",
      if (resolved) "review_label_before_matrix_use" else
        "exclude_from_compound_identity",
      NA_character_, "source record review",
      if (resolved) {
        "A source-backed structure is available, but the displayed value appears to be a plant, product, material, or source label; review and replace the display name before interpretation."
      } else {
        "The value appears to be a plant, product, material, or common-source label rather than a discrete compound."
      }
    ))
  }
  if (grepl("isomer|isomers|mixture|derivative|derivatives|fraction|extract",
            text, perl = TRUE)) {
    return(.plant_identity_issue_row(
      if (resolved) "resolved_mixture_or_isomer_group" else
        "mixture_or_isomer_group",
      "review_record_level_structure", name, "source record review",
      if (resolved) {
        "A source-backed structure is available, but the label describes a mixture, isomer group, derivative group, or extract fraction; review before using one structure as representative."
      } else {
        "The value describes a mixture, isomer group, derivative group, or extract fraction; it should not be forced to one structure."
      }
    ))
  }
  if (grepl("multiple lotus source structures|source_structure_ambiguous",
            notes, perl = TRUE)) {
    return(.plant_identity_issue_row(
      if (resolved) "resolved_source_structure_ambiguous" else
        "ambiguous_source_structure",
      "review_record_level_structure", name,
      "LOTUS source record",
      if (resolved) {
        "The compound resolved through PubChem or another identity source, but multiple LOTUS source structures map to the same normalized name; review if record-level stereochemistry or derivatives matter."
      } else {
        "A normalized compound key maps to more than one source-backed LOTUS structure; do not assign a single structure without record-level review."
      }
    ))
  }
  if (row$resolved %in% TRUE) {
    return(.plant_identity_issue_row(
      "resolved", "accept_resolved_identity", name, row$resolution_source,
      "Compound already has a usable CID, InChIKey, SMILES, or formula."
    ))
  }
  if (grepl(.plant_greek_letter_regex(), text, perl = TRUE)) {
    return(.plant_identity_issue_row(
      "remaining_greek_alias", "curate_synonym_or_source_structure",
      .pubchem_transliterate_name(name), "PubChem or source record",
      "A Greek-letter alias remains unresolved after deterministic alias retries; review source identifiers or add a vetted synonym only if supported."
    ))
  }
  if (.plant_identity_complex_name(text)) {
    return(.plant_identity_issue_row(
      "complex_structure_name_unresolved", "prefer_source_structure_or_keep_unresolved",
      name, "LOTUS source record or structure database",
      "The value appears to be a long systematic chemical name. Prefer source-backed SMILES/InChIKey; otherwise keep unresolved."
    ))
  }
  if (.plant_identity_short_name(text)) {
    return(.plant_identity_issue_row(
      "short_name_synonym_or_spelling_review", "curate_synonym_if_source_backed",
      name, "PubChem synonym/source record",
      "The value is a short chemical-like name. It may be a spelling variant, synonym, or uncommon trivial name; review before adding a mapping."
    ))
  }
  .plant_identity_issue_row(
    "unclassified_unresolved_identity", "manual_review_required", name,
    "source record review",
    "No conservative automated rule explains the unresolved identity."
  )
}

.plant_greek_letter_regex = function() {
  paste0("[", paste(vapply(c(0x03b1, 0x03b2, 0x03b3, 0x03b4,
                            0x03ba, 0x03bb, 0x03bc),
                          intToUtf8, character(1)),
                    collapse = ""), "]")
}

.plant_identity_review_required_types = function() {
  c("resolved_broad_class_label", "resolved_source_or_product_label",
    "resolved_mixture_or_isomer_group",
    "resolved_source_structure_ambiguous")
}

.plant_identity_issue_row = function(type, decision, suggested_query,
                                     suggested_source, reason) {
  list(type = type,
       decision = decision,
       suggested_query = .uaf_first_non_empty_text(suggested_query),
       suggested_source = .uaf_first_non_empty_text(suggested_source),
       reason = reason)
}

.plant_identity_broad_class_label = function(text) {
  grepl(paste(c(
    "\\balkaloids?\\b", "\\bflavonoids?\\b", "\\bterpenes?\\b",
    "\\bterpenoids?\\b", "\\bmonoterpenes?\\b", "\\bsesquiterpenes?\\b",
    "\\bditerpenes?\\b", "\\btriterpenes?\\b", "\\bphenolics?\\b",
    "\\bpolyphenols?\\b", "\\bphytosterols?\\b", "\\bsteroids?\\b",
    "\\bsaponins?\\b", "\\bglycosides?\\b", "\\blipids?\\b",
    "\\bfatty acids?\\b", "\\bcoumarins?\\b", "\\blignans?\\b",
    "\\bcarotenoids?\\b", "\\banthocyanins?\\b", "\\btannins?\\b"
  ), collapse = "|"), text, perl = TRUE)
}

.plant_identity_source_label = function(text) {
  grepl(paste(c(
    "\\bchamomile\\b", "\\btarragon\\b", "\\bessential oil\\b",
    "\\boil\\b", "\\bresin\\b", "\\bextract\\b", "\\bfraction\\b",
    "\\broot\\b", "\\bleaf\\b", "\\bflower\\b", "\\bseed\\b",
    "\\bwood\\b", "\\bbark\\b"
  ), collapse = "|"), text, perl = TRUE)
}

.plant_identity_complex_name = function(text) {
  nchar(text) > 120 ||
    grepl("\\[|\\]|\\{|\\}|\\([0-9rsze,+-]+\\)|lambda|cyclo|oxan|chromen|benzopyran|pentacyclo|tetracyclo",
          text, perl = TRUE)
}

.plant_identity_short_name = function(text) {
  nchar(text) <= 80 &&
    grepl("^[a-z0-9+/,.'() -]+$", text, perl = TRUE) &&
    grepl("[a-z]", text, perl = TRUE)
}

.plant_identity_review_priority = function(row, occ, issue_type) {
  if (row$resolved %in% TRUE &&
      !issue_type %in% .plant_identity_review_required_types()) {
    return("resolved")
  }
  count = suppressWarnings(as.integer(row$query_count))
  if (!is.finite(count)) count = nrow(occ)
  species_count = length(unique(.uaf_non_empty(occ$species)))
  if (issue_type %in% .plant_identity_review_required_types()) {
    if (count >= 3 || species_count >= 3) return("high")
    return("medium")
  }
  high_types = c("ambiguous_source_structure",
                 "short_name_synonym_or_spelling_review",
                 "remaining_greek_alias")
  low_types = c("broad_class_or_family_label", "source_or_product_label",
                "mixture_or_isomer_group")
  if (issue_type %in% high_types && (count >= 3 || species_count >= 3)) {
    return("high")
  }
  if (issue_type %in% low_types) return("low")
  if (count >= 5 || species_count >= 5) return("high")
  if (count >= 2 || species_count >= 2) return("medium")
  "low"
}

.plant_identity_review_priority_vector = function(resolved, query_count,
                                                   occurrence_count,
                                                   species_count,
                                                   issue_type) {
  n = length(issue_type)
  resolved = rep_len(resolved %in% TRUE, n)
  count = suppressWarnings(as.integer(rep_len(query_count, n)))
  occurrence_count = suppressWarnings(as.integer(rep_len(occurrence_count, n)))
  species_count = suppressWarnings(as.integer(rep_len(species_count, n)))
  count[!is.finite(count)] = occurrence_count[!is.finite(count)]
  count[!is.finite(count)] = 0L
  species_count[!is.finite(species_count)] = 0L
  issue_type = rep_len(as.character(issue_type), n)
  required = issue_type %in% .plant_identity_review_required_types()
  high_types = issue_type %in% c(
    "ambiguous_source_structure",
    "short_name_synonym_or_spelling_review",
    "remaining_greek_alias"
  )
  low_types = issue_type %in% c(
    "broad_class_or_family_label",
    "source_or_product_label",
    "mixture_or_isomer_group"
  )
  threshold_3 = count >= 3L | species_count >= 3L
  out = rep("low", n)
  out[count >= 2L | species_count >= 2L] = "medium"
  out[count >= 5L | species_count >= 5L] = "high"
  out[low_types] = "low"
  out[high_types & threshold_3] = "high"
  out[required] = "medium"
  out[required & threshold_3] = "high"
  out[resolved & !required] = "resolved"
  out
}

.plant_compound_identity_review_decision = function(x) {
  decision = .plant_matrix_key(.uaf_first_non_empty_text(x))
  if (is.na(decision) || decision == "") return("needs_review")
  if (decision %in% c("accept", "accepted", "accept_resolved",
                      "accept_resolved_identity", "accept_identity",
                      "keep_resolved")) {
    return("accept_resolved")
  }
  if (decision %in% c("update", "update_identity", "replace",
                      "replace_identity", "manual_identity",
                      "curate_identity", "promote_identity")) {
    return("update_identity")
  }
  if (decision %in% c("reject", "exclude", "remove", "drop")) {
    return("reject")
  }
  if (decision %in% c("keep", "needs_review", "review_later")) {
    return(decision)
  }
  "needs_review"
}

.plant_compound_identity_review_stamp = function(review, reviewer) {
  reviewer = .uaf_first_non_empty_text(review$reviewed_by, reviewer)
  reviewed_at = .uaf_first_non_empty_text(review$reviewed_at,
                                          .plant_timestamp())
  note = .uaf_first_non_empty_text(review$review_note)
  pieces = c(
    if (!is.na(reviewer)) paste0("reviewer: ", reviewer),
    if (!is.na(reviewed_at)) paste0("reviewed_at: ", reviewed_at),
    if (!is.na(note)) paste0("note: ", note)
  )
  if (length(pieces) < 1) return("")
  paste0("[", paste(pieces, collapse = "; "), "]")
}

.plant_reviewed_identity_values = function(review) {
  list(
    compound_name = .uaf_first_non_empty_text(review$proposed_compound_name),
    CID = .uaf_first_non_empty_text(review$proposed_cid),
    InChIKey = .uaf_first_non_empty_text(review$proposed_inchikey),
    SMILES = .uaf_first_non_empty_text(review$proposed_smiles),
    MolecularFormula = .uaf_first_non_empty_text(
      review$proposed_molecular_formula
    ),
    resolution_source = .uaf_first_non_empty_text(
      review$proposed_resolution_source
    )
  )
}

.plant_append_note = function(existing, note) {
  pieces = unique(.uaf_non_empty(c(existing, note)))
  .plant_truncate(.pubchem_collapse(pieces), 1200)
}

.plant_pubchem_identity_categorate = function(compounds, cache, cache_dir,
                                              throttle, batch_size, resume,
                                              progress, pubchem_fun,
                                              request_fun) {
  batch_size = .plant_enrichment_batch_size(batch_size, length(compounds))
  if (batch_size < length(compounds)) {
    return(.plant_pubchem_identity_categorate_batched(
      compounds = compounds,
      cache = cache,
      cache_dir = cache_dir,
      throttle = throttle,
      batch_size = batch_size,
      resume = resume,
      progress = progress,
      pubchem_fun = pubchem_fun,
      request_fun = request_fun
    ))
  }
  pubchem_fun = pubchem_fun %||% pubchemProfile
  pubchem_cache_dir = if (is.null(cache_dir)) NULL else
    file.path(cache_dir, "pubchem")
  pubchem = tryCatch(
    suppressWarnings(pubchem_fun(compounds = compounds,
                                 profile = "minimal",
                                 cache = cache,
                                 cache_dir = pubchem_cache_dir,
                                 throttle = throttle,
                                 include_annotations = FALSE,
                                 request_fun = request_fun)),
    error = function(error) {
      warning("PubChem identity resolution failed: ",
              conditionMessage(error), call. = FALSE)
      .categorate_empty_pubchem_profile(compounds, "minimal")
    }
  )
  .plant_identity_categorate_from_pubchem(pubchem, "pubchem_identity")
}

.plant_pubchem_identity_categorate_batched = function(compounds, cache,
                                                      cache_dir, throttle,
                                                      batch_size, resume,
                                                      progress, pubchem_fun,
                                                      request_fun) {
  groups = split(compounds, ceiling(seq_along(compounds) / batch_size))
  batch_dir = if (is.null(cache_dir)) {
    file.path(.plant_cache_dir(NULL), "compound_identity",
              "pubchem_identity_batches")
  } else {
    file.path(cache_dir, "pubchem_identity_batches")
  }
  if (isTRUE(cache)) dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
  results = vector("list", length(groups))
  for (i in seq_along(groups)) {
    batch = groups[[i]]
    cache_file = file.path(
      batch_dir,
      paste0("identity_batch_", sprintf("%04d", i), "_",
             .pubchem_url_hash(paste(batch, collapse = "\r")), ".rds")
    )
    if (isTRUE(cache) && isTRUE(resume) && file.exists(cache_file)) {
      cached_result = .plant_read_categorate_checkpoint(cache_file, batch)
      if (!is.null(cached_result)) {
        if (isTRUE(progress)) {
          message("uafR plant identity batch ", i, "/", length(groups),
                  ": using cached result")
        }
        results[[i]] = cached_result
        next
      }
      if (isTRUE(progress)) {
        message("uafR plant identity batch ", i, "/", length(groups),
                ": cached result is incomplete or unreadable; rebuilding")
      }
    }
    if (isTRUE(progress)) {
      message("uafR plant identity batch ", i, "/", length(groups),
              ": resolving ", length(batch), " compound(s)")
    }
    result = .plant_pubchem_identity_categorate(
      compounds = batch,
      cache = cache,
      cache_dir = cache_dir,
      throttle = throttle,
      batch_size = Inf,
      resume = FALSE,
      progress = FALSE,
      pubchem_fun = pubchem_fun,
      request_fun = request_fun
    )
    if (isTRUE(cache) &&
        .plant_categorate_checkpoint_cacheable(result, batch)) {
      .plant_atomic_save_rds(result, cache_file)
    } else if (isTRUE(cache) && isTRUE(progress)) {
      message("uafR plant identity batch ", i, "/", length(groups),
              ": result was not checkpointed because PubChem did not",
              " complete the identity pass")
    }
    results[[i]] = result
  }
  .plant_merge_identity_categorate_results(results, compounds,
                                           "pubchem_identity_batched")
}

.plant_identity_categorate_from_pubchem = function(pubchem,
                                                   enrichment_mode) {
  identity = if (is.list(pubchem) && is.data.frame(pubchem$identity)) {
    pubchem$identity
  } else {
    .uaf_empty_table(c("Query", "CID", "MatchStatus", "QueriedName",
                       "SourceURL"))
  }
  properties = if (is.list(pubchem) && is.data.frame(pubchem$properties)) {
    pubchem$properties
  } else {
    .uaf_empty_table(c("Query", "CID", "MolecularFormula", "InChIKey",
                       "CanonicalSMILES", "IsomericSMILES"))
  }
  synonyms = if (is.list(pubchem) && is.data.frame(pubchem$synonyms)) {
    pubchem$synonyms
  } else {
    .uaf_empty_table(c("Query", "CID", "Synonym"))
  }
  list(
    PubChemIdentity = identity,
    PubChemProperties = properties,
    PubChemSynonyms = synonyms,
    EnrichmentMode = enrichment_mode,
    ValidationSummary = data.frame(Status = "identity_only",
                                   stringsAsFactors = FALSE)
  )
}

.plant_merge_identity_categorate_results = function(results, compounds,
                                                    enrichment_mode) {
  results = results[vapply(results, is.list, logical(1))]
  if (length(results) < 1) return(.plant_empty_identity_categorate())
  out = list(
    PubChemIdentity = .plant_bind_flexible_tables(lapply(results, `[[`,
                                                        "PubChemIdentity")),
    PubChemProperties = .plant_bind_flexible_tables(lapply(results, `[[`,
                                                          "PubChemProperties")),
    PubChemSynonyms = .plant_bind_flexible_tables(lapply(results, `[[`,
                                                        "PubChemSynonyms")),
    SourceCompoundIdentity = .plant_bind_tables(lapply(results, `[[`,
                                                       "SourceCompoundIdentity"),
                                                .plant_source_compound_identity_cols()),
    EnrichmentMode = enrichment_mode,
    BatchCount = length(results),
    BatchCompoundCount = length(compounds),
    ValidationSummary = data.frame(Status = "identity_only",
                                   stringsAsFactors = FALSE)
  )
  out
}

.plant_empty_identity_categorate = function(enrichment_mode = "pubchem_identity") {
  list(
    PubChemIdentity = .uaf_empty_table(c("Query", "CID", "MatchStatus",
                                         "QueriedName", "SourceURL")),
    PubChemProperties = .uaf_empty_table(c("Query", "CID",
                                           "MolecularFormula", "InChIKey",
                                           "CanonicalSMILES",
                                           "IsomericSMILES")),
    PubChemSynonyms = .uaf_empty_table(c("Query", "CID", "Synonym")),
    SourceCompoundIdentity = .plant_empty_source_compound_identity(),
    EnrichmentMode = enrichment_mode,
    ValidationSummary = data.frame(Status = "identity_only",
                                   stringsAsFactors = FALSE)
  )
}

.plant_merge_compound_resolution = function(base, update) {
  cols = .plant_compound_resolution_cols()
  base = .plant_bind_tables(list(base), cols)
  update = .plant_bind_tables(list(update), cols)
  if (nrow(base) < 1) return(update)
  if (nrow(update) < 1) return(base)
  key = update$compound_name_clean
  keep = !base$compound_name_clean %in% key
  out = rbind(base[keep, , drop = FALSE], update)
  out = out[order(out$compound_name_clean), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_mark_unattempted_compounds = function(resolution,
                                            resolution_occurrences) {
  resolution = .plant_bind_tables(list(resolution),
                                  .plant_compound_resolution_cols())
  attempted = unique(.uaf_non_empty(resolution_occurrences$compound_name_clean))
  if (nrow(resolution) < 1) return(resolution)
  not_attempted = !resolution$compound_name_clean %in% attempted
  resolution$resolution_source[not_attempted] = "not_attempted"
  resolution$notes[not_attempted] =
    "Not selected for compound resolution in this run."
  resolution
}

.plant_mark_all_compounds_not_attempted = function(resolution) {
  resolution = .plant_bind_tables(list(resolution),
                                  .plant_compound_resolution_cols())
  if (nrow(resolution) < 1) return(resolution)
  resolution$resolution_source = "not_attempted"
  resolution$notes = "Compound resolution was not requested."
  resolution
}

.plant_pubchem_only_enrichment = function(compounds, detail, cache, cache_dir,
                                          throttle, pubchem_fun = NULL,
                                          request_fun = NULL,
                                          batch_size = Inf,
                                          resume = TRUE,
                                          progress = FALSE,
                                          ...) {
  dots = list(...)
  batch_size = .plant_enrichment_batch_size(batch_size, length(compounds))
  if (batch_size < length(compounds)) {
    return(.plant_pubchem_only_enrichment_batched(
      compounds = compounds,
      detail = detail,
      cache = cache,
      cache_dir = cache_dir,
      throttle = throttle,
      pubchem_fun = pubchem_fun,
      request_fun = request_fun,
      batch_size = batch_size,
      resume = resume,
      progress = progress,
      ...
    ))
  }
  pubchem_fun = pubchem_fun %||% pubchemProfile
  pubchem_profile = if (detail == "full") "full" else "safety"
  pubchem_sections = if (detail == "full") {
    c("Literature", "Patents")
  } else {
    c("Mass Spectrometry", "GC-MS", "MS-MS",
      "Names and Identifiers", "Chemical and Physical Properties",
      "Pharmacology and Biochemistry", "Drug and Medication Information")
  }
  pubchem_sources = c("LOTUS - the natural products occurrence database",
                      "Flavor and Extract Manufacturers Association (FEMA)",
                      "FDA/SPL Indexing Data",
                      "Medical Subject Headings (MeSH)")
  assay_detail_limit = dots$assay_detail_limit %||% 50
  trait_matrix_profile = dots$trait_matrix_profile %||% "core"
  trait_matrix_mode = dots$trait_matrix_mode %||% "binary"
  trait_matrix_min_confidence = dots$trait_matrix_min_confidence %||% 0
  trait_matrix_max_traits = dots$trait_matrix_max_traits %||% Inf

  pubchem_cache_dir = if (is.null(cache_dir)) NULL else
    file.path(cache_dir, "pubchem")
  pubchem = tryCatch(
    pubchem_fun(compounds = compounds,
                profile = pubchem_profile,
                sections = pubchem_sections,
                sources = pubchem_sources,
                cache = cache,
                cache_dir = pubchem_cache_dir,
                throttle = throttle,
                assay_detail_limit = assay_detail_limit,
                request_fun = request_fun),
    error = function(error) {
      warning("PubChem-only plant compound enrichment failed: ",
              conditionMessage(error), call. = FALSE)
      .categorate_empty_pubchem_profile(compounds, pubchem_profile)
    }
  )
  pubchem_profiles = .pubchem_extract_profiles(pubchem)
  kegg = .categorate_empty_kegg_profile()
  data_list = .plant_empty_categorate_data_list()
  normalized = .categorate_normalized_outputs(
    compounds = compounds,
    data_list = data_list,
    pubchem = pubchem,
    kegg = kegg,
    pubchem_profiles = pubchem_profiles,
    trait_matrix_profile = trait_matrix_profile,
    trait_matrix_mode = trait_matrix_mode,
    trait_matrix_min_confidence = trait_matrix_min_confidence,
    trait_matrix_max_traits = trait_matrix_max_traits
  )
  source_coverage = .categorate_source_coverage(
    compounds = compounds,
    data_list = data_list,
    pubchem = pubchem,
    kegg = kegg,
    pubchem_profiles = pubchem_profiles
  )
  derived_groups = .categorate_derived_groups(
    compounds = compounds,
    data_list = data_list,
    pubchem = pubchem,
    kegg = kegg,
    source_coverage = source_coverage,
    pubchem_profiles = pubchem_profiles,
    normalized = normalized
  )
  provenance = .categorate_research_provenance(pubchem = pubchem, kegg = kegg)
  result = list(
    PubChemIdentity = pubchem$identity,
    PubChemProperties = pubchem$properties,
    PubChemSynonyms = pubchem$synonyms,
    PubChemAnnotations = pubchem$annotations,
    PubChemSourceAnnotations = pubchem$source_annotations,
    PubChemSpectra = pubchem$spectra,
    PubChemSafety = pubchem$safety,
    PubChemExperimental = pubchem$experimental,
    PubChemBioactivity = pubchem$bioactivity,
    PubChemBioAssayDetails = pubchem$bioassay_details,
    PubChemIdentifiers = pubchem_profiles$PubChemIdentifierProfile,
    SafetyProfile = pubchem_profiles$SafetyProfile,
    FEMAProfile = pubchem_profiles$FEMAProfile,
    FDA_SPL_Profile = pubchem_profiles$FDA_SPL_Profile,
    LOTUSProfile = pubchem_profiles$LOTUSProfile,
    PubChemClassifications = pubchem_profiles$PubChemClassificationProfile,
    MeSHProfile = pubchem_profiles$MeSHProfile,
    LiteratureProfile = pubchem_profiles$LiteratureProfile,
    ChemicalTerms = normalized$ChemicalTerms,
    ChemicalTraits = normalized$ChemicalTraits,
    ChemicalTraitOntology = normalized$ChemicalTraitOntology,
    ChemicalTraitMatrix = normalized$ChemicalTraitMatrix,
    ChemicalTraitOntologyMatrix = normalized$ChemicalTraitOntologyMatrix,
    ChemicalTraitEvidence = normalized$ChemicalTraitEvidence,
    ChemicalTraitReport = normalized$ChemicalTraitReport,
    ChemicalTraitSummary = normalized$ChemicalTraitSummary,
    ChemicalTraitSimilarity = normalized$ChemicalTraitSimilarity,
    ChemicalClasses = normalized$ChemicalClasses,
    ChemicalMeasurements = normalized$ChemicalMeasurements,
    ChemicalMeasurementSummary = normalized$ChemicalMeasurementSummary,
    ChemicalHazards = normalized$ChemicalHazards,
    ChemicalUses = normalized$ChemicalUses,
    ChemicalBioassays = normalized$ChemicalBioassays,
    ChemicalBioactivities = normalized$ChemicalBioactivities,
    ChemicalTargets = normalized$ChemicalTargets,
    ChemicalPotencies = normalized$ChemicalPotencies,
    ChemicalTaxonomy = normalized$ChemicalTaxonomy,
    ChemicalOccurrences = normalized$ChemicalOccurrences,
    ChemicalPathwayRoles = normalized$ChemicalPathwayRoles,
    KEGGReactionParticipants = normalized$KEGGReactionParticipants,
    KEGGMatches = kegg$matches,
    KEGGSearchCandidates = kegg$search_candidates,
    KEGGRecords = kegg$records,
    KEGGIdentifiers = kegg$identifiers,
    KEGGPathways = kegg$pathways,
    KEGGReactions = kegg$reactions,
    KEGGEnzymes = kegg$enzymes,
    KEGGModules = kegg$modules,
    KEGGLinks = kegg$links,
    KEGGLinkMetadata = kegg$link_metadata,
    KEGGClassifications = kegg$classifications,
    SourceCoverage = source_coverage,
    DerivedGroups = derived_groups,
    Provenance = provenance,
    EnrichmentMode = "pubchem_only"
  )
  validation = validateCategorateResult(result)
  c(result, list(
    DataDictionary = validation$DataDictionary,
    TableQuality = validation$TableQuality,
    SourceDiagnostics = validation$SourceDiagnostics,
    ValidationIssues = validation$Issues,
    ValidationSummary = validation$Summary
  ))
}

.plant_enrichment_batch_size = function(batch_size, n_compounds) {
  batch_size = suppressWarnings(as.numeric(batch_size[[1]]))
  if (!is.finite(batch_size) || batch_size < 1) return(n_compounds)
  max(1, as.integer(batch_size))
}

.plant_pubchem_only_enrichment_batched = function(compounds, detail, cache,
                                                  cache_dir, throttle,
                                                  pubchem_fun, request_fun,
                                                  batch_size, resume,
                                                  progress, ...) {
  groups = split(compounds, ceiling(seq_along(compounds) / batch_size))
  results = vector("list", length(groups))
  batch_dir = if (is.null(cache_dir)) {
    file.path(.plant_cache_dir(NULL), "compound_enrichment",
              "pubchem_only_batches")
  } else {
    file.path(cache_dir, "pubchem_only_batches")
  }
  if (isTRUE(cache)) dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(groups)) {
    batch = groups[[i]]
    cache_file = file.path(
      batch_dir,
      paste0("batch_", sprintf("%04d", i), "_",
             .pubchem_url_hash(paste(c(detail, batch), collapse = "\r")),
             ".rds")
    )
    if (isTRUE(cache) && isTRUE(resume) && file.exists(cache_file)) {
      cached_result = .plant_read_categorate_checkpoint(cache_file, batch)
      if (!is.null(cached_result)) {
        if (isTRUE(progress)) {
          message("uafR plant enrichment batch ", i, "/", length(groups),
                  ": using cached result")
        }
        results[[i]] = cached_result
        next
      }
      if (isTRUE(progress)) {
        message("uafR plant enrichment batch ", i, "/", length(groups),
                ": cached result is incomplete or unreadable; rebuilding")
      }
    }
    if (isTRUE(progress)) {
      message("uafR plant enrichment batch ", i, "/", length(groups),
              ": resolving ", length(batch), " compound(s)")
    }
    result = .plant_pubchem_only_enrichment(
      compounds = batch,
      detail = detail,
      cache = cache,
      cache_dir = cache_dir,
      throttle = throttle,
      pubchem_fun = pubchem_fun,
      request_fun = request_fun,
      batch_size = Inf,
      resume = FALSE,
      progress = FALSE,
      ...
    )
    if (isTRUE(cache) &&
        .plant_categorate_checkpoint_cacheable(result, batch)) {
      .plant_atomic_save_rds(result, cache_file)
    } else if (isTRUE(cache) && isTRUE(progress)) {
      message("uafR plant enrichment batch ", i, "/", length(groups),
              ": result was not checkpointed because PubChem did not",
              " complete the identity pass")
    }
    results[[i]] = result
  }
  .plant_merge_categorate_results(results, compounds,
                                  enrichment_mode = "pubchem_only_batched")
}

.plant_merge_categorate_results = function(results, compounds,
                                           enrichment_mode) {
  results = results[vapply(results, is.list, logical(1))]
  if (length(results) < 1) return(NULL)
  names_all = unique(unlist(lapply(results, names), use.names = FALSE))
  merged = list()
  for (name in names_all) {
    values = lapply(results, `[[`, name)
    values = values[!vapply(values, is.null, logical(1))]
    if (length(values) < 1) next
    if (all(vapply(values, is.data.frame, logical(1)))) {
      merged[[name]] = .plant_bind_flexible_tables(values)
    } else {
      merged[[name]] = values[[1]]
    }
  }
  merged$EnrichmentMode = enrichment_mode
  merged$BatchCount = length(results)
  merged$BatchCompoundCount = length(compounds)
  validation = validateCategorateResult(merged)
  c(merged, list(
    DataDictionary = validation$DataDictionary,
    TableQuality = validation$TableQuality,
    SourceDiagnostics = validation$SourceDiagnostics,
    ValidationIssues = validation$Issues,
    ValidationSummary = validation$Summary
  ))
}

.plant_bind_flexible_tables = function(parts) {
  parts = parts[vapply(parts, is.data.frame, logical(1))]
  if (length(parts) < 1) return(data.frame())
  cols = unique(unlist(lapply(parts, names), use.names = FALSE))
  normalized = lapply(parts, function(part) {
    for (col in cols) if (!col %in% names(part)) part[[col]] = NA_character_
    part[, cols, drop = FALSE]
  })
  out = do.call(rbind, normalized)
  out = unique(out)
  row.names(out) = NULL
  out
}

.plant_batch_checkpoint_version = function() "3.3.0"

.plant_batch_subset_provider_results = function(provider_results,
                                                 species) {
  if (is.null(provider_results) || !is.list(provider_results)) {
    return(provider_results)
  }
  allowed = unique(.plant_clean_name(.uaf_non_empty(species)))
  if (length(allowed) < 1L) return(provider_results)
  out = provider_results
  for (provider in names(out)) {
    value = out[[provider]]
    if (is.data.frame(value)) {
      out[[provider]] = .plant_batch_subset_species_table(value, allowed)
      next
    }
    if (!is.list(value)) next
    occurrences = value$PlantCompoundOccurrences %||% value$occurrences
    literature = value$LiteratureCandidates %||% value$literature_candidates
    source_identity = value$SourceCompoundIdentity %||%
      value$source_compound_identity
    occurrences = .plant_batch_subset_species_table(occurrences, allowed)
    literature = .plant_batch_subset_species_table(literature, allowed)
    source_identity = .plant_batch_subset_source_identity(
      source_identity, occurrences
    )
    if (!is.null(value$PlantCompoundOccurrences)) {
      value$PlantCompoundOccurrences = occurrences
    } else if (!is.null(value$occurrences)) {
      value$occurrences = occurrences
    }
    if (!is.null(value$LiteratureCandidates)) {
      value$LiteratureCandidates = literature
    } else if (!is.null(value$literature_candidates)) {
      value$literature_candidates = literature
    }
    if (!is.null(value$SourceCompoundIdentity)) {
      value$SourceCompoundIdentity = source_identity
    } else if (!is.null(value$source_compound_identity)) {
      value$source_compound_identity = source_identity
    }
    out[[provider]] = value
  }
  out
}

.plant_batch_subset_species_table = function(x, allowed_species) {
  if (!is.data.frame(x) || nrow(x) < 1L) return(x)
  key_columns = intersect(
    c("species", "query_plant", "query_plant_clean"), names(x)
  )
  if (length(key_columns) < 1L) return(x)
  keep = rep(FALSE, nrow(x))
  for (column in key_columns) {
    values = if (identical(column, "query_plant_clean")) {
      .plant_clean_name(x[[column]])
    } else {
      .plant_clean_name(x[[column]])
    }
    keep = keep | (!is.na(values) & values %in% allowed_species)
  }
  x[keep, , drop = FALSE]
}

.plant_batch_subset_source_identity = function(source_identity,
                                                occurrences) {
  if (!is.data.frame(source_identity)) return(source_identity)
  if (nrow(source_identity) < 1L || !is.data.frame(occurrences) ||
      nrow(occurrences) < 1L) {
    return(source_identity[FALSE, , drop = FALSE])
  }
  keep = rep(FALSE, nrow(source_identity))
  if ("source_record_id" %in% names(source_identity) &&
      "source_record_id" %in% names(occurrences)) {
    occurrence_ids = .uaf_non_empty(occurrences$source_record_id)
    keep = keep | source_identity$source_record_id %in% occurrence_ids
  }
  if ("compound_name_clean" %in% names(source_identity)) {
    occurrence_names = if ("compound_name_clean" %in% names(occurrences)) {
      .uaf_non_empty(occurrences$compound_name_clean)
    } else if ("compound_name" %in% names(occurrences)) {
      .plant_clean_compound(occurrences$compound_name)
    } else {
      character()
    }
    keep = keep | source_identity$compound_name_clean %in% occurrence_names
  }
  source_identity[keep, , drop = FALSE]
}

.plant_batch_chunk_inputs = function(plants, plant_queries, species_chunks) {
  if (!is.data.frame(plants)) return(species_chunks)
  input = as.data.frame(plants, stringsAsFactors = FALSE)
  names(input) = .plant_normalize_column_names(names(input))
  species_field = if ("species" %in% names(input)) "species" else names(input)[[1]]
  input_species = .plant_canonical_taxon_name(input[[species_field]])
  lapply(species_chunks, function(chunk_species) {
    rows = input[input_species %in% chunk_species, , drop = FALSE]
    if (nrow(rows) < 1) {
      return(data.frame(species = chunk_species, stringsAsFactors = FALSE))
    }
    order_index = match(.plant_canonical_taxon_name(rows[[species_field]]),
                        chunk_species)
    rows = rows[order(order_index, na.last = TRUE), , drop = FALSE]
    row.names(rows) = NULL
    rows
  })
}

.plant_batch_query_signature = function(plant_queries, plant_aliases) {
  values = function(x) {
    if (!is.data.frame(x) || nrow(x) < 1) return("empty")
    x[] = lapply(x, function(value) {
      value = .uaf_squish_text(value)
      value[is.na(value)] = ""
      value
    })
    paste(c(names(x), unlist(x, use.names = FALSE)), collapse = "\r")
  }
  paste0("plant_query_contract:", .pubchem_url_hash(paste(
    values(plant_queries), values(plant_aliases), sep = "\n--aliases--\n"
  )))
}

.plant_batch_provider_results_signature = function(provider_results) {
  if (is.null(provider_results)) return("provider_results:none")
  path = tempfile("uafR_provider_results_", fileext = ".rds")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  serialized = tryCatch({
    saveRDS(provider_results, path, version = 3, compress = FALSE)
    TRUE
  }, error = function(error) FALSE)
  if (isTRUE(serialized) && file.exists(path)) {
    return(paste0(
      "provider_results:rds_md5=", unname(tools::md5sum(path)[[1L]]),
      ":bytes=", file.info(path)$size
    ))
  }
  paste0("provider_results:class=",
         paste(class(provider_results), collapse = "/"),
         ":length=", length(provider_results))
}

.plant_batch_run_signature = function(species, sources, taxon_fallback,
                                       provider_index_signature,
                                       max_pubmed_records,
                                       max_provider_records,
                                       provider_results,
                                       query_signature = NA_character_) {
  payload = c(
    "uafR_plant_discovery",
    .plant_batch_checkpoint_version(),
    paste(species, collapse = "\n"),
    paste(sort(unique(.uaf_non_empty(sources))), collapse = ";"),
    paste(.uaf_non_empty(taxon_fallback), collapse = ";"),
    provider_index_signature,
    .uaf_first_non_empty_text(query_signature,
                              "plant_query_contract:not_supplied"),
    paste0("max_pubmed_records=", max_pubmed_records),
    paste0("max_provider_records=", max_provider_records),
    .plant_batch_provider_results_signature(provider_results)
  )
  .pubchem_url_hash(paste(payload, collapse = "\r"))
}

.plant_batch_chunk_signature = function(run_signature, chunk_index, species) {
  .pubchem_url_hash(paste(c(run_signature, chunk_index, species),
                          collapse = "\r"))
}

.plant_batch_manifest_plan = function(species_chunks, checkpoint_files,
                                       run_signature, out_dir, cache) {
  starts = cumsum(c(1L, vapply(species_chunks, length, integer(1))))
  starts = starts[seq_along(species_chunks)]
  ends = starts + vapply(species_chunks, length, integer(1)) - 1L
  output_files = if (isTRUE(cache)) {
    if (is.null(out_dir)) {
      normalizePath(checkpoint_files, winslash = "/", mustWork = FALSE)
    } else {
      file.path("checkpoints", basename(checkpoint_files))
    }
  } else {
    rep(NA_character_, length(species_chunks))
  }
  data.frame(
    batch_index = seq_along(species_chunks),
    chunk_id = seq_along(species_chunks),
    query_start = starts,
    query_end = ends,
    query_count = vapply(species_chunks, length, integer(1)),
    query_label = vapply(species_chunks, paste, collapse = "; ",
                         FUN.VALUE = character(1)),
    species_count = vapply(species_chunks, length, integer(1)),
    species = vapply(species_chunks, paste, collapse = "; ",
                     FUN.VALUE = character(1)),
    status = "not_started",
    previous_status = NA_character_,
    started_at = NA_character_,
    finished_at = NA_character_,
    completed_at = NA_character_,
    elapsed_seconds = 0,
    cache_hit_count = 0L,
    request_count = 0L,
    retry_count = 0L,
    error_count = 0L,
    warning_count = 0L,
    error_message = NA_character_,
    occurrence_count = 0L,
    literature_candidate_count = 0L,
    provider_diagnostic_count = 0L,
    output_file = output_files,
    output_path = output_files,
    checkpoint_file = output_files,
    checkpoint_read_status = "not_checked",
    checkpoint_status = "not_checked",
    checkpoint_cache_hit = "No",
    checkpoint_version = .plant_batch_checkpoint_version(),
    run_signature = run_signature,
    chunk_signature = vapply(seq_along(species_chunks), function(i) {
      .plant_batch_chunk_signature(run_signature, i, species_chunks[[i]])
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

.plant_batch_merge_prior_manifest = function(manifest, out_dir,
                                              run_signature, resume,
                                              overwrite) {
  if (is.null(out_dir) || !isTRUE(resume)) return(manifest)
  path = file.path(out_dir, "discovery_chunk_manifest.csv")
  if (!file.exists(path)) return(manifest)
  prior = tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(condition) NULL
  )
  if (!is.data.frame(prior) || nrow(prior) < 1) return(manifest)
  prior_signatures = if ("run_signature" %in% names(prior)) {
    unique(.uaf_non_empty(prior$run_signature))
  } else {
    character()
  }
  if (length(prior_signatures) > 0 &&
      !identical(prior_signatures, run_signature)) {
    if (!isTRUE(overwrite)) {
      stop(
        "Existing discovery manifest was created for different inputs. Use a ",
        "new `out_dir`, or set `overwrite = TRUE` after reviewing the prior run.",
        call. = FALSE
      )
    }
    return(manifest)
  }
  if (!"chunk_signature" %in% names(prior)) return(manifest)
  idx = match(manifest$chunk_signature, prior$chunk_signature)
  matched = which(!is.na(idx))
  if (length(matched) < 1) return(manifest)
  if ("status" %in% names(prior)) {
    manifest$previous_status[matched] = as.character(prior$status[idx[matched]])
  }
  if ("retry_count" %in% names(prior)) {
    retries = suppressWarnings(as.integer(prior$retry_count[idx[matched]]))
    retries[is.na(retries)] = 0L
    manifest$retry_count[matched] = retries
  }
  manifest
}

.plant_valid_discovery_result = function(result) {
  required = c("PlantQueries", "ProviderDiagnostics",
               "PlantCompoundOccurrences", "LiteratureCandidates")
  is.list(result) && all(required %in% names(result)) &&
    all(vapply(result[required], is.data.frame, logical(1)))
}

.plant_batch_checkpoint_envelope = function(result, run_signature,
                                             chunk_signature, species) {
  list(
    format = "uafR_plant_discovery_checkpoint",
    checkpoint_version = .plant_batch_checkpoint_version(),
    run_signature = run_signature,
    chunk_signature = chunk_signature,
    query_species = as.character(species),
    status = "completed",
    created_at = .plant_timestamp(),
    result = result
  )
}

.plant_read_batch_checkpoint = function(path, expected_run_signature,
                                        expected_chunk_signature,
                                        expected_species) {
  if (!file.exists(path)) {
    return(list(ok = FALSE, reason = "missing", result = NULL))
  }
  object = tryCatch(readRDS(path), error = function(condition) condition)
  if (inherits(object, "condition")) {
    return(list(ok = FALSE, reason = "unreadable", result = NULL))
  }
  if (!is.list(object) ||
      !identical(object$format, "uafR_plant_discovery_checkpoint")) {
    return(list(ok = FALSE, reason = "format_mismatch", result = NULL))
  }
  if (!identical(as.character(object$checkpoint_version),
                 .plant_batch_checkpoint_version())) {
    return(list(ok = FALSE, reason = "version_mismatch", result = NULL))
  }
  if (!identical(as.character(object$run_signature),
                 as.character(expected_run_signature)) ||
      !identical(as.character(object$chunk_signature),
                 as.character(expected_chunk_signature))) {
    return(list(ok = FALSE, reason = "signature_mismatch", result = NULL))
  }
  if (!identical(as.character(object$query_species),
                 as.character(expected_species))) {
    return(list(ok = FALSE, reason = "species_mismatch", result = NULL))
  }
  if (!identical(as.character(object$status), "completed") ||
      !.plant_valid_discovery_result(object$result)) {
    return(list(ok = FALSE, reason = "invalid_result", result = NULL))
  }
  list(ok = TRUE, reason = "valid", result = object$result)
}

.plant_atomic_replace = function(temp_file, path) {
  renamed = suppressWarnings(file.rename(temp_file, path))
  if (!isTRUE(renamed)) {
    copied = file.copy(temp_file, path, overwrite = TRUE, copy.mode = TRUE,
                       copy.date = TRUE)
    if (!isTRUE(copied)) {
      stop("Could not atomically replace file: ", path, call. = FALSE)
    }
  }
  invisible(path)
}

.plant_atomic_save_rds = function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_file = tempfile(paste0(".", basename(path), "_"),
                       tmpdir = dirname(path))
  on.exit(if (file.exists(temp_file)) unlink(temp_file), add = TRUE)
  saveRDS(object, temp_file)
  .plant_atomic_replace(temp_file, path)
}

.plant_atomic_write_csv = function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_file = tempfile(paste0(".", basename(path), "_"),
                       tmpdir = dirname(path))
  on.exit(if (file.exists(temp_file)) unlink(temp_file), add = TRUE)
  utils::write.csv(x, temp_file, row.names = FALSE, na = "")
  .plant_atomic_replace(temp_file, path)
}

.plant_atomic_write_json = function(x, path, dataframe = "rows",
                                     pretty = TRUE, auto_unbox = TRUE,
                                     na = "null") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_file = tempfile(paste0(".", basename(path), "_"),
                       tmpdir = dirname(path))
  on.exit(if (file.exists(temp_file)) unlink(temp_file), add = TRUE)
  jsonlite::write_json(
    x, temp_file, dataframe = dataframe, pretty = pretty,
    auto_unbox = auto_unbox, na = na
  )
  .plant_atomic_replace(temp_file, path)
}

.plant_atomic_write_text = function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_file = tempfile(paste0(".", basename(path), "_"),
                       tmpdir = dirname(path))
  on.exit(if (file.exists(temp_file)) unlink(temp_file), add = TRUE)
  writeLines(enc2utf8(as.character(x)), temp_file, useBytes = TRUE)
  .plant_atomic_replace(temp_file, path)
}

.plant_categorate_checkpoint_cacheable = function(result,
                                                   expected_compounds) {
  if (!is.list(result) || !is.data.frame(result$PubChemIdentity) ||
      nrow(result$PubChemIdentity) < 1) {
    return(FALSE)
  }
  identity = result$PubChemIdentity
  if (!"MatchStatus" %in% names(identity)) return(FALSE)
  status = tolower(.uaf_non_empty(identity$MatchStatus))
  if (length(status) < 1 ||
      all(status %in% c("not_run", "error", "failed", "timeout",
                        "rate_limited", "service_unavailable"))) {
    return(FALSE)
  }
  expected = length(unique(.uaf_non_empty(expected_compounds)))
  observed = if ("Query" %in% names(identity)) {
    length(unique(.uaf_non_empty(identity$Query)))
  } else {
    nrow(identity)
  }
  observed >= expected
}

.plant_read_categorate_checkpoint = function(path, expected_compounds) {
  object = tryCatch(readRDS(path), error = function(condition) NULL)
  if (!.plant_categorate_checkpoint_cacheable(object, expected_compounds)) {
    return(NULL)
  }
  object
}

.plant_batch_failed_queries = function(manifest) {
  failed_status = c("failed", "error", "timeout", "timed_out",
                    "rate_limited", "stopped")
  out = manifest[tolower(manifest$status) %in% failed_status, , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_write_batch_operational_files = function(manifest, out_dir) {
  if (is.null(out_dir)) return(invisible(NULL))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  failed = .plant_batch_failed_queries(manifest)
  retry = writePlantChemistryRetryQueue(manifest)
  .plant_atomic_write_csv(
    manifest, file.path(out_dir, "discovery_chunk_manifest.csv")
  )
  .plant_atomic_write_csv(failed, file.path(out_dir, "failed_queries.csv"))
  .plant_atomic_write_csv(retry, file.path(out_dir, "retry_queue.csv"))
  invisible(list(Manifest = manifest, FailedQueries = failed,
                 RetryQueue = retry))
}

.plant_batch_result_counts = function(result) {
  diagnostics = if (is.list(result) &&
                    is.data.frame(result$ProviderDiagnostics)) {
    result$ProviderDiagnostics
  } else {
    data.frame()
  }
  sum_column = function(name) {
    if (!name %in% names(diagnostics)) return(0L)
    value = suppressWarnings(as.numeric(diagnostics[[name]]))
    as.integer(sum(value[is.finite(value)], na.rm = TRUE))
  }
  list(
    occurrence_count = if (is.list(result) &&
                              is.data.frame(result$PlantCompoundOccurrences)) {
      nrow(result$PlantCompoundOccurrences)
    } else 0L,
    literature_candidate_count = if (is.list(result) &&
                                        is.data.frame(result$LiteratureCandidates)) {
      nrow(result$LiteratureCandidates)
    } else 0L,
    provider_diagnostic_count = nrow(diagnostics),
    request_count = sum_column("request_count"),
    cache_hit_count = sum_column("cache_hit_count"),
    error_count = sum_column("error_count"),
    warning_count = sum_column("warning_count")
  )
}

.plant_batch_result_health = function(result) {
  counts = .plant_batch_result_counts(result)
  diagnostics = result$ProviderDiagnostics
  messages = character()
  if (is.data.frame(diagnostics) && nrow(diagnostics) > 0) {
    error_rows = rep(FALSE, nrow(diagnostics))
    if ("error_count" %in% names(diagnostics)) {
      error_rows = suppressWarnings(as.numeric(diagnostics$error_count)) > 0
      error_rows[is.na(error_rows)] = FALSE
    }
    if ("status" %in% names(diagnostics)) {
      error_rows = error_rows | tolower(as.character(diagnostics$status)) %in%
        c("error", "failed", "timeout", "timed_out", "rate_limited")
    }
    for (column in intersect(c("error_messages", "message"),
                             names(diagnostics))) {
      messages = c(messages, as.character(diagnostics[[column]][error_rows]))
    }
  }
  messages = unique(.uaf_non_empty(messages))
  error_message = if (length(messages) > 0) {
    paste(messages, collapse = " | ")
  } else if (counts$error_count > 0) {
    "One or more selected providers reported an error."
  } else {
    NA_character_
  }
  service_busy_text = paste(messages, collapse = " ")
  service_busy = counts$error_count > 0 &&
    .plant_service_busy_message(service_busy_text)
  list(error_count = counts$error_count,
       warning_count = counts$warning_count,
       service_busy = service_busy,
       error_message = error_message)
}

.plant_empty_batch_result = function(species, taxon_fallback,
                                     error_message = NA_character_) {
  plant_queries = .plant_queries(species, taxon_fallback)
  occurrences = .plant_empty_occurrences()
  diagnostics = .plant_provider_diagnostics(
    "batch", TRUE, TRUE, FALSE, 0, 0, 0,
    ifelse(is.na(error_message), 0, 1), 0,
    ifelse(is.na(error_message), "no_records", "error"),
    .uaf_first_non_empty_text(error_message, "No provider records found.")
  )
	  out = list(
	    PlantQueries = plant_queries,
	    PlantQueryAliases = .plant_query_aliases(species, plant_queries),
	    PlantNameResolution = .plant_name_resolution(plant_queries),
	    ProviderDiagnostics = diagnostics,
	    ProviderQueryAccounting =
	      .uaf_empty_table(.plant_provider_query_accounting_cols()),
	    ProviderResourceManifest =
	      .uaf_empty_table(.plant_provider_resource_manifest_cols()),
	    PlantCompoundOccurrences = occurrences,
	    PlantContextEvidence = .uaf_empty_table(.plant_context_evidence_cols()),
	    ProviderContextAudit =
	      .uaf_empty_table(.plant_provider_context_audit_cols()),
	    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    SourceCompoundIdentity = .plant_empty_source_compound_identity(),
    CompoundResolution = .plant_compound_resolution(occurrences, NULL),
    CompoundIdentityReview = .plant_empty_compound_identity_review(),
    CategorateResult = NULL,
    SpeciesChemistrySummary = summarizePlantPhytochemistry(
      plant_compounds = occurrences,
      categorate_result = NULL,
      compound_resolution = .plant_compound_resolution(occurrences, NULL),
      plant_queries = plant_queries,
      provider_diagnostics = diagnostics
    ),
    SpeciesChemistryMatrix = .uaf_empty_table(c("species")),
    ChemistryComparability = .uaf_empty_table(.plant_comparability_cols()),
    ComparableChemistryMatrix = .uaf_empty_table(c("species")),
    TraitEvidence = .uaf_empty_table(.plant_trait_evidence_cols()),
    Validation = list(Summary = data.frame(),
                      TableQuality = data.frame(),
                      ProviderDiagnostics = data.frame(),
                      Issues = data.frame(),
                      DataDictionary = data.frame()),
    DataDictionary = plantPhytochemistrySchema(),
    Provenance = .plant_provenance(
      "batch_discovery", "batch",
      paste(plant_queries$query_plant, collapse = "; "),
      NA_character_, 0,
      .uaf_first_non_empty_text(error_message,
                                "Empty batch result assembled.")
    )
  )
  class(out) = c("uaf_plant_phytochemistry", class(out))
  out$Validation = validatePlantPhytochemistryResult(out)
  out
}

.plant_combine_batch_results = function(results, plant_queries,
                                        taxon_fallback,
                                        link_source_context = FALSE,
                                        defer_derived = FALSE,
                                        validate_result = TRUE) {
  results = results[vapply(results, is.list, logical(1))]
  if (length(results) < 1) {
    return(.plant_empty_batch_result(plant_queries$species, taxon_fallback))
  }
  plant_queries = unique(.plant_bind_tables(list(plant_queries),
                                            .plant_query_cols()))
  occurrences = .plant_bind_occurrences(
    lapply(results, `[[`, "PlantCompoundOccurrences")
  )
  occurrences = .plant_match_occurrences_to_queries(occurrences,
                                                    plant_queries,
                                                    normalize = FALSE)
  literature = .plant_bind_literature(
    lapply(results, `[[`, "LiteratureCandidates")
  )
  literature = .plant_match_literature_to_queries(
    literature, plant_queries, normalize = FALSE
  )
  context_evidence = .plant_bind_tables(
    lapply(results, `[[`, "PlantContextEvidence"),
    .plant_context_evidence_cols()
  )
  if (!isTRUE(defer_derived)) {
    derived_context = plantContextEvidence(occurrences)
    context_evidence = .plant_normalize_context_evidence(
      .plant_bind_tables(
        list(context_evidence, derived_context),
        .plant_context_evidence_cols()
      )
    )
  }
  if (!isTRUE(defer_derived) && isTRUE(link_source_context) &&
      nrow(literature) > 0) {
    source_context = .plant_source_context_evidence(
      occurrences = occurrences,
      context_sources = literature,
      min_confidence = "low"
    )
    context_evidence = .plant_normalize_context_evidence(
      .plant_bind_tables(
        list(context_evidence, source_context),
        .plant_context_evidence_cols()
      )
    )
  }
  if (isTRUE(defer_derived)) {
    provider_context_audit =
      .uaf_empty_table(.plant_provider_context_audit_cols())
  } else {
    occurrences = .plant_apply_context_evidence(occurrences, context_evidence)
    occurrences = .plant_clean_context_conflicts(occurrences)
    provider_context_audit = plantProviderContextAudit(
      list(PlantCompoundOccurrences = occurrences,
           PlantContextEvidence = context_evidence)
    )
  }
  diagnostics = .plant_bind_tables(lapply(results, `[[`, "ProviderDiagnostics"),
                                   .plant_provider_diagnostic_cols())
  provenance = .plant_bind_tables(lapply(results, `[[`, "Provenance"),
                                  .plant_provenance_cols())
  compound_resolution = .plant_compound_resolution(occurrences, NULL)
  if (isTRUE(defer_derived)) {
    summary = .uaf_empty_table(.plant_summary_cols())
    matrix = .uaf_empty_table(c("species"))
    comparability = .uaf_empty_table(.plant_comparability_cols())
    comparable_matrix = .uaf_empty_table(c("species"))
    compound_identity_review = .plant_empty_compound_identity_review()
  } else {
    comparability = plantChemistryComparability(
      list(PlantCompoundOccurrences = occurrences, CategorateResult = NULL),
      min_confidence = "low"
    )
    summary = summarizePlantPhytochemistry(
      plant_compounds = occurrences,
      categorate_result = NULL,
      compound_resolution = compound_resolution,
      plant_queries = plant_queries,
      provider_diagnostics = diagnostics,
      comparability = comparability
    )
    matrix = plantPhytochemistryMatrix(
      list(PlantCompoundOccurrences = occurrences, CategorateResult = NULL),
      level = "species",
      profile = "core",
      mode = "binary",
      min_confidence = "medium",
      max_traits = .plant_automatic_matrix_max_traits()
    )
    comparable_matrix = plantComparableChemistryMatrix(
      list(ChemistryComparability = comparability),
      comparison_scope = "specialized_metabolites",
      level = "species",
      mode = "binary",
      min_comparability_confidence = "medium"
    )
    compound_identity_review = plantCompoundIdentityReviewTable(
      list(PlantCompoundOccurrences = occurrences,
           CompoundResolution = compound_resolution)
    )
  }
  out = list(
    PlantQueries = plant_queries,
    PlantQueryAliases = .plant_merge_alias_tables(results, plant_queries),
    PlantNameResolution = .plant_name_resolution(plant_queries),
    ProviderDiagnostics = diagnostics,
    ProviderQueryAccounting = .plant_merge_accounting_tables(
      results, plant_queries
    ),
    ProviderResourceManifest = .plant_bind_unique_tables(
      lapply(results, `[[`, "ProviderResourceManifest"),
      .plant_provider_resource_manifest_cols(),
      c("provider", "resource_type", "resource_id", "md5")
    ),
    PlantCompoundOccurrences = occurrences,
    PlantContextEvidence = context_evidence,
    ProviderContextAudit = provider_context_audit,
    LiteratureCandidates = literature,
    SourceCompoundIdentity = .plant_merge_source_identity_tables(
      lapply(results, `[[`, "SourceCompoundIdentity")
    ),
    CompoundResolution = compound_resolution,
    CompoundIdentityReview = compound_identity_review,
    CategorateResult = NULL,
    SpeciesChemistrySummary = summary,
    SpeciesChemistryMatrix = matrix,
    ChemistryComparability = comparability,
    ComparableChemistryMatrix = comparable_matrix,
    TraitEvidence = .uaf_empty_table(.plant_trait_evidence_cols()),
    Validation = list(Summary = data.frame(),
                      TableQuality = data.frame(),
                      ProviderDiagnostics = data.frame(),
                      Issues = data.frame(),
                      DataDictionary = data.frame()),
    DataDictionary = plantPhytochemistrySchema(),
    Provenance = .plant_bind_tables(list(
      provenance,
      .plant_provenance("batch_discovery", "uafR",
                        paste(plant_queries$query_plant, collapse = "; "),
                        NA_character_, nrow(occurrences),
                        "Chunked species-first discovery results combined.")
    ), .plant_provenance_cols())
  )
  class(out) = c("uaf_plant_phytochemistry", class(out))
  if (!isTRUE(defer_derived) && isTRUE(validate_result)) {
    out$Validation = validatePlantPhytochemistryResult(out)
  }
  out
}

.plant_bind_literature = function(parts) {
  out = .plant_normalize_literature(
    .plant_bind_tables(parts, .plant_literature_cols())
  )
  if (nrow(out) < 2L) return(out)
  key_cols = setdiff(names(out), "retrieved_at")
  key_data = out[key_cols]
  key_data[] = lapply(key_data, function(value) {
    value = .uaf_squish_text(value)
    value[is.na(value)] = ""
    value
  })
  key = do.call(paste, c(key_data, sep = "||"))
  out = out[!duplicated(key), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_batch_resolution_occurrences = function(occurrences, occurrence_status,
                                               analysis_ready, min_confidence,
                                               max_compounds_per_species,
                                               max_unique_compounds) {
  occurrences = filterPlantPhytochemistryEvidence(
    occurrences,
    occurrence_status = occurrence_status,
    analysis_ready = analysis_ready,
    min_confidence = min_confidence
  )
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) return(occurrences)
  score = .plant_occurrence_rank_score(occurrences$occurrence_status) +
    .plant_confidence_score(occurrences$confidence) +
    suppressWarnings(as.numeric(occurrences$evidence_quality_score))
  score[!is.finite(score)] = 0
  occurrences$.uaf_resolution_score = score
  occurrences = occurrences[order(occurrences$species,
                                  -occurrences$.uaf_resolution_score,
                                  occurrences$compound_name_clean), ,
                            drop = FALSE]
  occurrences = occurrences[!duplicated(paste(occurrences$species,
                                              occurrences$compound_name_clean,
                                              sep = "\r")), ,
                            drop = FALSE]
  max_per_species = suppressWarnings(as.numeric(max_compounds_per_species[[1]]))
  if (is.finite(max_per_species) && max_per_species >= 1) {
    parts = lapply(split(occurrences, occurrences$species), function(part) {
      part[seq_len(min(nrow(part), as.integer(max_per_species))), ,
           drop = FALSE]
    })
    occurrences = do.call(rbind, parts)
  }
  max_unique = suppressWarnings(as.numeric(max_unique_compounds[[1]]))
  if (is.finite(max_unique) && max_unique >= 1) {
    ranked = occurrences[order(-occurrences$.uaf_resolution_score,
                               occurrences$compound_name_clean), ,
                         drop = FALSE]
    keep_compounds = unique(ranked$compound_name_clean)
    keep_compounds = keep_compounds[seq_len(min(length(keep_compounds),
                                            as.integer(max_unique)))]
    occurrences = occurrences[occurrences$compound_name_clean %in%
                                keep_compounds, , drop = FALSE]
  }
  occurrences$.uaf_resolution_score = NULL
  row.names(occurrences) = NULL
  occurrences
}

.plant_batch_run_manifest = function(result, chunk_manifest, started,
                                     completed, sources,
                                     compound_resolution_profile,
                                     resolution_occurrences, out_dir,
                                     cache_dir, discovery_complete,
                                     compound_stage_status, run_signature,
                                     pause_reason = NA_character_) {
  elapsed = round(as.numeric(difftime(completed, started, units = "secs")), 3)
  resolution = if (is.data.frame(result$CompoundResolution)) {
    result$CompoundResolution
  } else {
    .uaf_empty_table(.plant_compound_resolution_cols())
  }
  resolved = sum(resolution$resolved %in% TRUE)
  not_attempted = sum(resolution$resolution_source == "not_attempted",
                      na.rm = TRUE)
  attempted = nrow(resolution) - not_attempted
  attempted_unresolved = sum(!(resolution$resolved %in% TRUE) &
                               resolution$resolution_source != "not_attempted",
                             na.rm = TRUE)
  comparability = if (is.data.frame(result$ChemistryComparability)) {
    result$ChemistryComparability
  } else {
    .uaf_empty_table(.plant_comparability_cols())
  }
  data.frame(
    run_id = paste0("plant_batch_", format(started, "%Y%m%d_%H%M%S")),
    run_signature = run_signature,
    checkpoint_version = .plant_batch_checkpoint_version(),
    run_status = ifelse(discovery_complete, "completed", "incomplete"),
    discovery_complete = .uaf_yes_no(discovery_complete),
    completed_chunk_count = sum(chunk_manifest$status == "completed"),
    failed_chunk_count = sum(chunk_manifest$status %in%
                               c("failed", "error", "timeout", "timed_out",
                                 "rate_limited", "stopped")),
    pending_chunk_count = sum(!chunk_manifest$status %in%
                                c("completed", "failed", "error", "timeout",
                                  "timed_out", "rate_limited", "stopped")),
    pause_reason = .uaf_first_non_empty_text(pause_reason),
    started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    completed_at = format(completed, "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_seconds = elapsed,
    plant_count = nrow(result$PlantQueries),
    chunk_count = nrow(chunk_manifest),
    source_count = length(unique(.uaf_non_empty(sources))),
    sources = paste(unique(.uaf_non_empty(sources)), collapse = "; "),
    occurrence_count = nrow(result$PlantCompoundOccurrences),
    analysis_ready_occurrence_count = sum(result$PlantCompoundOccurrences$analysis_ready == "Yes"),
    context_evidence_count = if (is.data.frame(result$PlantContextEvidence)) {
      nrow(result$PlantContextEvidence)
    } else {
      0
    },
    literature_candidate_count = nrow(result$LiteratureCandidates),
    review_required_count = nrow(plantPhytochemistryReviewTable(result)),
    unique_compound_count = length(unique(.uaf_non_empty(
      result$PlantCompoundOccurrences$compound_name_clean))),
    comparable_occurrence_count = sum(comparability$comparable_for_matrix ==
                                        "Yes", na.rm = TRUE),
    comparison_scope_count = length(unique(.uaf_non_empty(
      comparability$comparison_scope[comparability$comparison_scope !=
                                       "unknown"]))),
    resolution_input_occurrence_count = nrow(resolution_occurrences),
    compound_resolution_profile = compound_resolution_profile,
    compound_stage_status = compound_stage_status,
    attempted_compound_count = attempted,
    resolved_compound_count = resolved,
    attempted_unresolved_compound_count = attempted_unresolved,
    not_attempted_compound_count = not_attempted,
    unresolved_compound_count = sum(!(resolution$resolved %in% TRUE)),
    out_dir = .uaf_first_non_empty_text(out_dir),
    cache_dir = cache_dir,
    validation_status = if (is.list(result$Validation) &&
                            is.data.frame(result$Validation$Summary) &&
                            "Status" %in% names(result$Validation$Summary)) {
      .uaf_first_non_empty_text(result$Validation$Summary$Status)
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )
}

.plant_write_batch_outputs = function(result, out_dir, overwrite,
                                       include_analysis = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(include_analysis)) {
    analysis_ready = filterPlantPhytochemistryEvidence(
      result,
      occurrence_status = c("direct_reported", "curated_reported"),
      analysis_ready = TRUE,
      min_confidence = "medium"
    )$PlantCompoundOccurrences
    review_required = plantPhytochemistryReviewTable(result)
    context_report = .plant_context_coverage_report(
      result$PlantCompoundOccurrences
    )
    tables = list(
      all_occurrences = result$PlantCompoundOccurrences,
      plant_context_evidence = result$PlantContextEvidence,
      provider_context_audit = result$ProviderContextAudit,
      analysis_ready_occurrences = analysis_ready,
      review_required_occurrences = review_required,
      compound_identity_resolution = result$CompoundResolution,
      compound_identity_review = result$CompoundIdentityReview,
      species_chemistry_summary = result$SpeciesChemistrySummary,
      species_chemistry_matrix = result$SpeciesChemistryMatrix,
      chemistry_comparability = result$ChemistryComparability,
      comparable_chemistry_matrix = result$ComparableChemistryMatrix,
      provider_diagnostics = result$ProviderDiagnostics,
      context_coverage_report = context_report,
      batch_run_manifest = result$BatchRunManifest,
      batch_chunk_manifest = result$BatchChunkManifest
    )
  } else {
    tables = list(
      all_occurrences = result$PlantCompoundOccurrences,
      plant_context_evidence = result$PlantContextEvidence,
      literature_candidates = result$LiteratureCandidates,
      source_compound_identity = result$SourceCompoundIdentity,
      provider_diagnostics = result$ProviderDiagnostics,
      provider_query_accounting = result$ProviderQueryAccounting,
      provider_resource_manifest = result$ProviderResourceManifest,
      compound_identity_resolution = result$CompoundResolution,
      batch_run_manifest = result$BatchRunManifest,
      batch_chunk_manifest = result$BatchChunkManifest
    )
  }
  manifest_rows = list()
  for (name in names(tables)) {
    file = file.path(out_dir, paste0(name, ".csv"))
    .plant_write_csv_file(tables[[name]], file, overwrite)
    manifest_rows[[length(manifest_rows) + 1]] = data.frame(
      artifact = name,
      file = file,
      row_count = nrow(tables[[name]]),
      column_count = ncol(tables[[name]]),
      stringsAsFactors = FALSE
    )
  }
  operational = list(
    failed_queries = result$FailedQueries,
    retry_queue = result$RetryQueue,
    discovery_chunk_manifest = result$BatchChunkManifest
  )
  for (name in names(operational)) {
    file = file.path(out_dir, paste0(name, ".csv"))
    if (!file.exists(file)) .plant_atomic_write_csv(operational[[name]], file)
    manifest_rows[[length(manifest_rows) + 1]] = data.frame(
      artifact = name,
      file = file,
      row_count = nrow(operational[[name]]),
      column_count = ncol(operational[[name]]),
      stringsAsFactors = FALSE
    )
  }
  json_file = file.path(out_dir, "run_manifest.json")
  .plant_write_json_file(result$BatchRunManifest, json_file, overwrite)
  manifest_rows[[length(manifest_rows) + 1]] = data.frame(
    artifact = "run_manifest_json",
    file = json_file,
    row_count = nrow(result$BatchRunManifest),
    column_count = ncol(result$BatchRunManifest),
    stringsAsFactors = FALSE
  )
  manifest = do.call(rbind, manifest_rows)
  manifest_file = file.path(out_dir, "batch_export_manifest.csv")
  .plant_write_csv_file(manifest, manifest_file, overwrite)
  row.names(manifest) = NULL
  manifest
}

.plant_pilot_summary = function(result) {
  summary = if (is.data.frame(result$SpeciesChemistrySummary)) {
    result$SpeciesChemistrySummary
  } else {
    .uaf_empty_table(.plant_summary_cols())
  }
  cols = .plant_pilot_summary_cols()
  if (nrow(summary) < 1) return(.uaf_empty_table(cols))
  context = .plant_context_coverage_report(result$PlantCompoundOccurrences)
  review = plantPhytochemistryReviewTable(result)
  review_count = .plant_count_by_species(review, "review_required_count")
  provider_count = .plant_count_by_species(result$PlantCompoundOccurrences,
                                           "provider_occurrence_count")
  out = merge(summary, context[, c("species", "context_known_fraction"),
                               drop = FALSE],
              by = "species", all.x = TRUE)
  out = merge(out, review_count, by = "species", all.x = TRUE)
  out = merge(out, provider_count, by = "species", all.x = TRUE)
  zero_cols = c("compound_count", "analysis_ready_compound_count",
                "direct_species_compound_count",
                "candidate_occurrence_count",
                "taxon_fallback_occurrence_count",
                "resolved_compound_count", "unresolved_compound_count",
                "context_known_occurrence_count",
                "review_required_count", "provider_occurrence_count",
                "comparable_specialized_compound_count",
                "comparable_volatile_compound_count",
                "comparable_primary_compound_count")
  for (col in zero_cols) {
    if (!col %in% names(out)) out[[col]] = 0
    value = suppressWarnings(as.numeric(out[[col]]))
    value[is.na(value)] = 0
    out[[col]] = value
  }
  if (!"context_known_fraction" %in% names(out)) {
    out$context_known_fraction = 0
  }
  out$context_known_fraction = suppressWarnings(
    as.numeric(out$context_known_fraction)
  )
  out$context_known_fraction[is.na(out$context_known_fraction)] = 0
  if (!"source_coverage_score" %in% names(out)) {
    out$source_coverage_score = 0
  }
  out$source_coverage_score = suppressWarnings(
    as.numeric(out$source_coverage_score)
  )
  out$source_coverage_score[is.na(out$source_coverage_score)] = 0
  if (!"uafR_validation_status" %in% names(out)) {
    out$uafR_validation_status = NA_character_
  }
  out$pilot_status = .plant_pilot_status(out)
  out$recommended_next_step = .plant_pilot_next_step(out)
  for (col in cols) if (!col %in% names(out)) out[[col]] = NA
  out = out[, cols, drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_pilot_summary_cols = function() {
  c("species", "query_status", "compound_count",
    "analysis_ready_compound_count", "direct_species_compound_count",
    "candidate_occurrence_count", "taxon_fallback_occurrence_count",
    "resolved_compound_count", "unresolved_compound_count",
    "context_known_occurrence_count", "context_known_fraction",
    "provider_occurrence_count", "review_required_count",
    "source_database_count", "source_databases",
    "comparable_specialized_compound_count",
    "comparable_volatile_compound_count",
    "comparable_primary_compound_count", "comparison_scopes",
    "source_coverage_score", "uafR_validation_status", "pilot_status",
    "recommended_next_step")
}

.plant_count_by_species = function(x, count_col) {
  if (!is.data.frame(x) || nrow(x) < 1 || !"species" %in% names(x)) {
    return(.uaf_empty_table(c("species", count_col)))
  }
  species = .uaf_non_empty(x$species)
  if (length(species) < 1) return(.uaf_empty_table(c("species", count_col)))
  tab = as.data.frame(table(species), stringsAsFactors = FALSE)
  names(tab) = c("species", count_col)
  tab[[count_col]] = as.integer(tab[[count_col]])
  tab
}

.plant_pilot_status = function(summary) {
  out = rep("ready_for_review", nrow(summary))
  out[summary$compound_count < 1] = "no_public_records"
  out[summary$compound_count > 0 &
        summary$analysis_ready_compound_count < 1] = "review_only"
  out[summary$analysis_ready_compound_count > 0 &
        summary$context_known_fraction < 0.25] = "context_sparse"
  out[summary$analysis_ready_compound_count > 0 &
        summary$resolved_compound_count < 1] = "identity_unresolved"
  out[summary$review_required_count > summary$analysis_ready_compound_count &
        summary$analysis_ready_compound_count > 0] = "review_heavy"
  out[summary$analysis_ready_compound_count > 0 &
        summary$context_known_fraction >= 0.25 &
        summary$resolved_compound_count > 0 &
        summary$review_required_count <= summary$analysis_ready_compound_count] =
    "pilot_ready"
  out
}

.plant_pilot_next_step = function(summary) {
  status = summary$pilot_status
  out = rep("Inspect species summary, occurrences, context evidence, and matrices before scaling.",
            length(status))
  out[status == "no_public_records"] =
    "No provider records were found; try alternate names or curated intake."
  out[status == "review_only"] =
    "Review candidate/fallback records before using this species analytically."
  out[status == "context_sparse"] =
    "Prioritize sources or manual review that can confirm plant part, tissue, or method."
  out[status == "identity_unresolved"] =
    "Resolve compound names or provide identifiers before matrix interpretation."
  out[status == "review_heavy"] =
    "Review the flagged occurrence rows before expanding this species set."
  out[status == "pilot_ready"] =
    "Use context-aware comparable matrices and spot-check provenance before scaling."
  out
}

.plant_qa_row = function(check, status, value, threshold, details,
                         recommendation) {
  data.frame(
    check = check,
    status = status,
    value = as.character(value),
    threshold = threshold,
    details = details,
    recommendation = recommendation,
    stringsAsFactors = FALSE
  )
}

.plant_qa_threshold_status = function(value, threshold) {
  value = suppressWarnings(as.numeric(value))
  threshold = suppressWarnings(as.numeric(threshold))
  if (!is.finite(value) || !is.finite(threshold)) return("warning")
  ifelse(value >= threshold, "pass", "warning")
}

.plant_pilot_matrices = function(result, min_confidence, feature, mode,
                                 max_features) {
  definitions = list(
    matrix_specialized_metabolites = list(
      comparison_scope = "specialized_metabolites",
      plant_part_group = NULL,
      require_context = FALSE,
      description = "Comparable specialized-metabolite groups."
    ),
    matrix_volatile_specialized_metabolites = list(
      comparison_scope = "volatile_specialized_metabolites",
      plant_part_group = NULL,
      require_context = FALSE,
      description = "Comparable volatile specialized chemistry."
    ),
    matrix_leaf_associated_chemistry = list(
      comparison_scope = "all_classified",
      plant_part_group = "leaf",
      require_context = TRUE,
      description = "Classified chemistry reported from leaf context."
    ),
    matrix_root_exudate_associated_chemistry = list(
      comparison_scope = "all_classified",
      plant_part_group = c("root_belowground", "exudate_rhizosphere"),
      require_context = TRUE,
      description = "Classified chemistry reported from root or rhizosphere/exudate context."
    ),
    matrix_primary_metabolites = list(
      comparison_scope = "primary_metabolites",
      plant_part_group = NULL,
      require_context = FALSE,
      description = "Comparable primary-metabolism chemistry."
    )
  )
  matrices = lapply(definitions, function(definition) {
    plantComparableChemistryMatrix(
      result,
      comparison_scope = definition$comparison_scope,
      feature = feature,
      mode = mode,
      min_comparability_confidence = min_confidence,
      plant_part_group = definition$plant_part_group,
      require_context = definition$require_context,
      max_features = max_features
    )
  })
  attr(matrices, "definitions") = definitions
  matrices
}

.plant_write_pilot_outputs = function(result, pilot_summary, review_needed,
                                      pilot_matrices, qa_report, out_dir,
                                      overwrite) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rows = list()
  add_artifact = function(artifact, table, file_name, description) {
    file = file.path(out_dir, file_name)
    .plant_write_csv_file(table, file, overwrite)
    rows[[length(rows) + 1]] <<- data.frame(
      artifact = artifact,
      file = file,
      row_count = nrow(table),
      column_count = ncol(table),
      description = description,
      stringsAsFactors = FALSE
    )
  }
  add_artifact("pilot_summary", pilot_summary, "pilot_summary.csv",
               "Compact species-level pilot status and next-step summary.")
  add_artifact("pilot_review_needed", review_needed,
               "pilot_review_needed.csv",
               "Candidate, fallback, unresolved, or review-required occurrence rows.")
	  add_artifact("pilot_qa_report", qa_report, "pilot_qa_report.csv",
	               "Overall pilot readiness checks and recommendations.")
	  add_artifact("pilot_provider_context_audit", result$ProviderContextAudit,
	               "pilot_provider_context_audit.csv",
	               "Provider-level biological-context coverage and review burden.")
	  definitions = attr(pilot_matrices, "definitions")
  for (name in names(pilot_matrices)) {
    description = if (is.list(definitions) &&
                       !is.null(definitions[[name]]$description)) {
      definitions[[name]]$description
    } else {
      "Pilot comparable-chemistry matrix."
    }
    add_artifact(name, pilot_matrices[[name]], paste0(name, ".csv"),
                 description)
  }
  run_manifest = .plant_pilot_run_manifest(result, pilot_summary,
                                           pilot_matrices)
  json_file = file.path(out_dir, "pilot_run_manifest.json")
  .plant_write_json_file(run_manifest, json_file, overwrite)
  rows[[length(rows) + 1]] = data.frame(
    artifact = "pilot_run_manifest_json",
    file = json_file,
    row_count = nrow(run_manifest),
    column_count = ncol(run_manifest),
    description = "One-row pilot run summary for automation and review.",
    stringsAsFactors = FALSE
  )
  manifest = do.call(rbind, rows)
  manifest_file = file.path(out_dir, "pilot_export_manifest.csv")
  .plant_write_csv_file(manifest, manifest_file, overwrite)
  row.names(manifest) = NULL
  manifest
}

.plant_pilot_run_manifest = function(result, pilot_summary, pilot_matrices) {
  validation_status = if (is.list(result$Validation) &&
                          is.data.frame(result$Validation$Summary) &&
                          "Status" %in% names(result$Validation$Summary)) {
    .uaf_first_non_empty_text(result$Validation$Summary$Status)
  } else {
    NA_character_
  }
  data.frame(
    created_at = .plant_timestamp(),
    plant_count = if (is.data.frame(result$PlantQueries)) {
      nrow(result$PlantQueries)
    } else {
      0
    },
    occurrence_count = if (is.data.frame(result$PlantCompoundOccurrences)) {
      nrow(result$PlantCompoundOccurrences)
    } else {
      0
    },
    analysis_ready_occurrence_count =
      if (is.data.frame(result$PlantCompoundOccurrences) &&
          "analysis_ready" %in% names(result$PlantCompoundOccurrences)) {
        sum(result$PlantCompoundOccurrences$analysis_ready == "Yes",
            na.rm = TRUE)
      } else {
        0
      },
    review_required_count = if (is.data.frame(result$PilotReviewNeeded)) {
      nrow(result$PilotReviewNeeded)
    } else {
      sum(pilot_summary$review_required_count, na.rm = TRUE)
    },
    pilot_ready_species_count = sum(pilot_summary$pilot_status ==
                                      "pilot_ready", na.rm = TRUE),
    matrix_count = length(pilot_matrices),
    validation_status = validation_status,
    stringsAsFactors = FALSE
  )
}

.plant_context_coverage_report = function(occurrences) {
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) {
    return(.uaf_empty_table(c("species", "occurrence_count",
                              "context_known_count", "context_known_fraction",
                              "plant_part_known_count", "tissue_known_count",
                              "method_known_count")))
  }
  parts = lapply(split(occurrences, occurrences$species), function(group) {
    part_known = !group$plant_part_group %in%
      c("unknown", "extract_unspecified")
    tissue_known = !group$tissue_group %in%
      c("unknown", "extract_unspecified")
    method_known = !group$method_group %in%
      c("unknown", "database_record", "literature_curation")
    context_known = group$biological_context_status != "context_missing"
    data.frame(
      species = group$species[[1]],
      occurrence_count = nrow(group),
      context_known_count = sum(context_known),
      context_known_fraction = round(sum(context_known) / nrow(group), 3),
      plant_part_known_count = sum(part_known),
      tissue_known_count = sum(tissue_known),
      method_known_count = sum(method_known),
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, parts)
  row.names(out) = NULL
  out
}

.plant_provider_context_audit = function(occurrences, context_evidence) {
  cols = .plant_provider_context_audit_cols()
  occurrences = .plant_normalize_occurrences(occurrences)
  context_evidence = .plant_normalize_context_evidence(context_evidence)
  sources = unique(.uaf_non_empty(c(occurrences$source_database,
                                    context_evidence$source_database)))
  if (length(sources) < 1) return(.uaf_empty_table(cols))
  parts = lapply(sources, function(source) {
    group = occurrences[occurrences$source_database == source, , drop = FALSE]
    context = context_evidence[context_evidence$source_database == source, ,
                               drop = FALSE]
    occurrence_count = nrow(group)
    part_known = if (occurrence_count > 0) {
      !group$plant_part_group %in% c("unknown", "extract_unspecified")
    } else logical()
    tissue_known = if (occurrence_count > 0) {
      !group$tissue_group %in% c("unknown", "extract_unspecified")
    } else logical()
    method_known = if (occurrence_count > 0) {
      !group$method_group %in%
        c("unknown", "database_record", "literature_curation")
    } else logical()
    context_known = if (occurrence_count > 0) {
      group$biological_context_status != "context_missing"
    } else logical()
    evidence_count = nrow(context)
    review_count = if (evidence_count > 0) {
      sum(context$requires_review == "Yes", na.rm = TRUE)
    } else 0
    context_fraction = .plant_safe_fraction(sum(context_known), occurrence_count)
    review_fraction = .plant_safe_fraction(review_count, evidence_count)
    missing_fraction = .plant_safe_fraction(sum(!context_known),
                                            occurrence_count)
    candidate_count = if (occurrence_count > 0) {
      sum(group$occurrence_status == "candidate", na.rm = TRUE)
    } else 0
    audit_status = .plant_provider_context_audit_status(
      occurrence_count = occurrence_count,
      candidate_count = candidate_count,
      context_fraction = context_fraction,
      review_fraction = review_fraction
    )
    data.frame(
      source_database = source,
      occurrence_count = occurrence_count,
      species_count = length(unique(.uaf_non_empty(group$species))),
      compound_count = length(unique(.uaf_non_empty(group$compound_name_clean))),
      direct_reported_count = sum(group$occurrence_status == "direct_reported",
                                  na.rm = TRUE),
      candidate_count = candidate_count,
      taxon_fallback_count = sum(group$occurrence_status == "taxon_fallback",
                                 na.rm = TRUE),
      analysis_ready_count = sum(group$analysis_ready == "Yes", na.rm = TRUE),
      context_known_count = sum(context_known),
      context_known_fraction = round(context_fraction, 3),
      plant_part_known_count = sum(part_known),
      tissue_known_count = sum(tissue_known),
      method_known_count = sum(method_known),
      context_evidence_count = evidence_count,
      high_confidence_context_count =
        sum(context$context_confidence == "high", na.rm = TRUE),
      medium_confidence_context_count =
        sum(context$context_confidence == "medium", na.rm = TRUE),
      low_confidence_context_count =
        sum(context$context_confidence == "low", na.rm = TRUE),
      review_required_context_count = review_count,
      review_required_context_fraction = round(review_fraction, 3),
      context_missing_count = sum(!context_known),
      context_missing_fraction = round(missing_fraction, 3),
      top_context_groups = .plant_top_terms(context$normalized_context),
      top_extraction_rules = .plant_top_terms(context$extraction_rule),
      audit_status = audit_status,
      recommended_action = .plant_provider_context_recommendation(audit_status,
                                                                  source),
      stringsAsFactors = FALSE
    )
  })
  out = .plant_bind_tables(parts, cols)
  out = out[order(out$audit_status, out$source_database), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_provider_context_audit_status = function(occurrence_count,
                                                candidate_count,
                                                context_fraction,
                                                review_fraction) {
  if (occurrence_count < 1) return("no_occurrences")
  if (candidate_count >= occurrence_count) return("candidate_only")
  if (!is.finite(context_fraction) || context_fraction < 0.10) {
    return("context_sparse")
  }
  if (is.finite(review_fraction) && review_fraction > 0.75) {
    return("review_heavy")
  }
  "context_ready"
}

.plant_provider_context_recommendation = function(audit_status, source) {
  provider_key = .plant_provider_key(source)
  if (audit_status == "no_occurrences") {
    return("No occurrence records were returned; inspect provider diagnostics or alternate plant names.")
  }
  if (audit_status == "candidate_only") {
    return("Provider output is candidate evidence only; do not treat context or occurrence as confirmed without manual review.")
  }
  if (audit_status == "context_sparse") {
    if (provider_key == "knapsack") {
      return(paste("KNApSAcK organism rows usually support species-compound",
                   "association and record-level review links, but often do",
                   "not expose plant part, tissue, or method in the organism",
                   "table. Use evidence_url/source records or curated intake",
                   "before context-aware analyses."))
    }
    if (provider_key == "lotus") {
      return(paste("LOTUS simple-search rows usually support taxon-compound",
                   "association and chemistry classification evidence, but",
                   "may not expose plant part, tissue, or analytical method.",
                   "Use source record review or curated intake before",
                   "context-aware analyses."))
    }
    return(paste("Context is sparse for", source,
                 "rows; inspect source records for plant part, tissue, and method before context-aware analyses."))
  }
  if (audit_status == "review_heavy") {
    return(paste("Most extracted context for", source,
                 "requires review; curate or harden provider-specific rules before scaling."))
  }
  "Provider has usable context coverage for pilot interpretation; spot-check source evidence before scaling."
}

.plant_safe_fraction = function(numerator, denominator) {
  numerator = suppressWarnings(as.numeric(numerator))
  denominator = suppressWarnings(as.numeric(denominator))
  if (!is.finite(numerator) || !is.finite(denominator) || denominator <= 0) {
    return(0)
  }
  numerator / denominator
}

.plant_write_csv_file = function(x, file, overwrite) {
  if (file.exists(file) && !isTRUE(overwrite)) {
    stop("Output file exists and `overwrite = FALSE`: ", file,
         call. = FALSE)
  }
  utils::write.csv(x, file, row.names = FALSE, na = "")
}

.plant_write_json_file = function(x, file, overwrite) {
  if (file.exists(file) && !isTRUE(overwrite)) {
    stop("Output file exists and `overwrite = FALSE`: ", file,
         call. = FALSE)
  }
  jsonlite::write_json(x, path = file, dataframe = "rows",
                       pretty = TRUE, na = "null")
}

.plant_empty_categorate_data_list = function() {
  list(
    reactives = .uaf_empty_table(c("Chemical", "reactives")),
    LOTUS = .uaf_empty_table(c("Chemical", "LOTUS")),
    KEGG = .uaf_empty_table(c("Chemical", "KEGG")),
    FEMA = .uaf_empty_table(c("Chemical", "FEMA")),
    FDA_SPL = .uaf_empty_table(c("Chemical", "FDA_SPL"))
  )
}

.plant_enrichment_source = function(categorate_result) {
  mode = if (is.list(categorate_result)) {
    .uaf_first_non_empty_text(categorate_result$EnrichmentMode)
  } else {
    NA_character_
  }
  if (!is.na(mode) && mode %in% c("pubchem_identity",
                                  "pubchem_identity_batched")) {
    return("pubchemProfile_identity")
  }
  if (!is.na(mode) && mode %in% c("source_identity_only",
                                  "source_pubchem_identity")) {
    return("source_compound_identity")
  }
  if (!is.na(mode) && mode %in% c("pubchem_only", "pubchem_only_batched")) {
    return("pubchemProfile")
  }
  "categorate"
}

.plant_enrichment_note = function(categorate_result) {
  mode = if (is.list(categorate_result)) {
    .uaf_first_non_empty_text(categorate_result$EnrichmentMode)
  } else {
    NA_character_
  }
  if (identical(mode, "pubchem_identity_batched")) {
    return("PubChem identity-only compound resolution completed in resumable batches.")
  }
  if (identical(mode, "pubchem_identity")) {
    return("PubChem identity-only compound resolution completed.")
  }
  if (identical(mode, "source_identity_only")) {
    return("Source-backed compound identity resolution completed without PubChem lookups.")
  }
  if (identical(mode, "source_pubchem_identity")) {
    return("Source-backed compound identity resolution completed before PubChem fallback.")
  }
  if (identical(mode, "pubchem_only_batched")) {
    return("PubChem-only compound enrichment completed in resumable batches.")
  }
  if (identical(mode, "pubchem_only")) {
    return("PubChem-only compound enrichment completed.")
  }
  "Compound enrichment completed or supplied."
}

.plant_categorate_properties = function(categorate_result) {
  cols = c("Query", "CID", "MolecularFormula", "InChIKey",
           "CanonicalSMILES", "IsomericSMILES", "SMILES",
           "ConnectivitySMILES")
  identity = if (is.list(categorate_result) &&
                 is.data.frame(categorate_result$PubChemIdentity)) {
    categorate_result$PubChemIdentity
  } else if (is.list(categorate_result) &&
             is.data.frame(categorate_result$identity)) {
    categorate_result$identity
  } else {
    data.frame()
  }
  if (is.list(categorate_result) &&
      is.data.frame(categorate_result$PubChemProperties)) {
    props = categorate_result$PubChemProperties
  } else if (is.list(categorate_result) &&
             is.data.frame(categorate_result$properties)) {
    props = categorate_result$properties
  } else {
    props = .uaf_empty_table(cols)
  }
  for (col in cols) {
    if (!col %in% names(props)) props[[col]] = rep(NA_character_, nrow(props))
  }
  props = props[, cols, drop = FALSE]
  if (is.data.frame(identity) && nrow(identity) > 0 &&
      all(c("Query", "CID") %in% names(identity))) {
    identity = identity[!is.na(identity$CID) & identity$CID != "", ,
                        drop = FALSE]
    if (nrow(identity) > 0) {
      identity_rows = data.frame(
        Query = identity$Query,
        CID = identity$CID,
        MolecularFormula = NA_character_,
        InChIKey = NA_character_,
        CanonicalSMILES = NA_character_,
        IsomericSMILES = NA_character_,
        SMILES = NA_character_,
        ConnectivitySMILES = NA_character_,
        stringsAsFactors = FALSE
      )
      prop_key = paste(.plant_clean_compound(props$Query), props$CID,
                       sep = "\r")
      identity_key = paste(.plant_clean_compound(identity_rows$Query),
                           identity_rows$CID, sep = "\r")
      identity_rows = identity_rows[!identity_key %in% prop_key, ,
                                    drop = FALSE]
      props = rbind(props, identity_rows)
    }
  }
  props
}

.plant_trait_evidence = function(occurrences, categorate_result) {
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1 || !is.list(categorate_result) ||
      !is.data.frame(categorate_result$ChemicalTraitEvidence)) {
    return(.uaf_empty_table(.plant_trait_evidence_cols()))
  }
  evidence = categorate_result$ChemicalTraitEvidence
  if (!"Query" %in% names(evidence)) {
    return(.uaf_empty_table(.plant_trait_evidence_cols()))
  }
  rows = list()
  for (i in seq_len(nrow(occurrences))) {
    occ = occurrences[i, , drop = FALSE]
    hit = evidence[.plant_clean_compound(evidence$Query) ==
                     occ$compound_name_clean, , drop = FALSE]
    if (nrow(hit) < 1) next
    rows[[length(rows) + 1]] = data.frame(
      species = occ$species,
      species_slug = .plant_slug(occ$species),
      compound_name = occ$compound_name,
      compound_name_clean = occ$compound_name_clean,
      EvidenceType = .plant_col_or_default(hit, "EvidenceType", NA_character_),
      AnalysisKey = .plant_col_or_default(hit, "AnalysisKey", NA_character_),
      SourceDatabase = .plant_col_or_default(hit, "SourceDatabase", NA_character_),
      EvidenceText = .plant_col_or_default(hit, "EvidenceText", NA_character_),
      EvidenceURL = .plant_col_or_default(hit, "EvidenceURL", NA_character_),
      Confidence = .plant_col_or_default(hit, "Confidence", NA_character_),
      stringsAsFactors = FALSE
    )
  }
  .plant_bind_tables(rows, .plant_trait_evidence_cols())
}

.plant_summary_trait_info = function(group, categorate_result) {
  compounds = unique(.uaf_non_empty(group$compound_name_clean))
  traits = if (is.list(categorate_result) &&
               is.data.frame(categorate_result$ChemicalTraits)) {
    categorate_result$ChemicalTraits
  } else {
    .uaf_empty_table(c("Query", "TraitType", "TraitGroup", "TraitValue",
                       "Confidence", "SourceDatabase"))
  }
  classes = if (is.list(categorate_result) &&
                is.data.frame(categorate_result$ChemicalClasses)) {
    categorate_result$ChemicalClasses
  } else {
    .uaf_empty_table(c("Query", "ClassName", "ClassGroup",
                       "Superclass", "Class", "Subclass"))
  }
  keep_traits = if ("Query" %in% names(traits)) {
    .plant_clean_compound(traits$Query) %in% compounds
  } else logical(nrow(traits))
  keep_classes = if ("Query" %in% names(classes)) {
    .plant_clean_compound(classes$Query) %in% compounds
  } else logical(nrow(classes))
  traits = traits[keep_traits, , drop = FALSE]
  classes = classes[keep_classes, , drop = FALSE]
  text = tolower(.pubchem_collapse(c(group$compound_name, group$evidence_text,
                                     traits$TraitType, traits$TraitGroup,
                                     traits$TraitValue, classes$ClassName,
                                     classes$ClassGroup)))
  if (nrow(traits) < 1L && nrow(classes) < 1L) {
    compound_index = split(seq_len(nrow(group)), group$compound_name_clean)
    compound_text = vapply(compound_index, function(idx) {
      tolower(.pubchem_collapse(group$compound_name[idx]))
    }, character(1), USE.NAMES = FALSE)
    fraction = function(pattern) {
      if (length(compound_text) < 1L) return(0)
      mean(grepl(pattern, compound_text, perl = TRUE))
    }
    return(list(
      superclasses = NA_character_,
      classes = NA_character_,
      subclasses = NA_character_,
      dominant_classes = NA_character_,
      kingdoms = ifelse(grepl("plantae|plant", text),
                        "Plantae", NA_character_),
      volatile_fraction =
        fraction("volatile|odor|aroma|essential oil|gc-ms|gc ms"),
      lipophilic_fraction =
        fraction("lipophilic|xlogp|logp|terpene|fatty|alkane|lipid"),
      oxygenated_fraction =
        fraction("oxygen|hydroxy|phenol|acid|ester|alcohol|ketone|aldehyde"),
      nitrogenous_fraction =
        fraction("nitrogen|alkaloid|amine|amide|nitrile"),
      sulfur_fraction = fraction("sulfur|thiol|sulfide|sulfoxide"),
      halogen_fraction = fraction("chloro|bromo|fluoro|iodo|halogen"),
      kegg_pathway_groups = NA_character_,
      chemical_trait_count = 0L,
      high_confidence_trait_count = 0L
    ))
  }
  fraction = function(pattern) {
    if (length(compounds) < 1) return(0)
    hits = vapply(compounds, function(compound) {
      compound_text = tolower(.pubchem_collapse(c(
        group$compound_name[group$compound_name_clean == compound],
        traits$TraitValue[.plant_clean_compound(traits$Query) == compound],
        traits$TraitGroup[.plant_clean_compound(traits$Query) == compound]
      )))
      grepl(pattern, compound_text, perl = TRUE)
    }, logical(1))
    mean(hits)
  }
  list(
    superclasses = .pubchem_collapse(.plant_col_or_default(classes, "Superclass", NA_character_)),
    classes = .pubchem_collapse(c(.plant_col_or_default(classes, "Class", NA_character_),
                                  .plant_col_or_default(classes, "ClassName", NA_character_))),
    subclasses = .pubchem_collapse(.plant_col_or_default(classes, "Subclass", NA_character_)),
    dominant_classes = .plant_top_terms(c(.plant_col_or_default(classes, "ClassName", NA_character_),
                                          .plant_col_or_default(traits, "TraitGroup", NA_character_))),
    kingdoms = ifelse(grepl("plantae|plant", text), "Plantae", NA_character_),
    volatile_fraction = fraction("volatile|odor|aroma|essential oil|gc-ms|gc ms"),
    lipophilic_fraction = fraction("lipophilic|xlogp|logp|terpene|fatty|alkane|lipid"),
    oxygenated_fraction = fraction("oxygen|hydroxy|phenol|acid|ester|alcohol|ketone|aldehyde"),
    nitrogenous_fraction = fraction("nitrogen|alkaloid|amine|amide|nitrile"),
    sulfur_fraction = fraction("sulfur|thiol|sulfide|sulfoxide"),
    halogen_fraction = fraction("chloro|bromo|fluoro|iodo|halogen"),
    kegg_pathway_groups = .pubchem_collapse(.plant_col_or_default(traits, "PathwayGroup", NA_character_)),
    chemical_trait_count = nrow(traits),
    high_confidence_trait_count = sum(tolower(.plant_col_or_default(traits, "Confidence", "")) == "high")
  )
}

.plant_wide_matrix = function(long, id_col, mode, max_traits) {
  if (nrow(long) < 1) return(.uaf_empty_table(c(id_col)))
  long = long[!is.na(long$id) & long$id != "" &
                !is.na(long$trait) & long$trait != "", , drop = FALSE]
  if (nrow(long) < 1) return(.uaf_empty_table(c(id_col)))
  traits = names(sort(table(long$trait), decreasing = TRUE))
  max_traits = suppressWarnings(as.numeric(max_traits[[1]]))
  if (is.finite(max_traits)) traits = utils::head(traits, max_traits)
  ids = sort(unique(.uaf_non_empty(long$id)))
  long = long[long$trait %in% traits & long$id %in% ids, , drop = FALSE]
  if (nrow(long) < 1 || length(ids) < 1 || length(traits) < 1) {
    return(.uaf_empty_table(c(id_col)))
  }
  long$id = factor(long$id, levels = ids)
  long$trait = factor(long$trait, levels = traits)
  if (mode == "binary") {
    keyed = unique(long[, c("id", "trait"), drop = FALSE])
    keyed$value = 1
  } else if (mode == "count") {
    keyed = long[, c("id", "trait"), drop = FALSE]
    keyed$value = 1
  } else {
    values = suppressWarnings(as.numeric(long$value))
    values[!is.finite(values)] = 0
    keyed = stats::aggregate(values,
                             by = list(id = long$id, trait = long$trait),
                             FUN = max, na.rm = TRUE)
    names(keyed)[names(keyed) == "x"] = "value"
  }
  mat = stats::xtabs(value ~ id + trait, data = keyed,
                     drop.unused.levels = FALSE)
  out = data.frame(stats::setNames(list(rownames(mat)), id_col),
                   as.data.frame.matrix(mat),
                   stringsAsFactors = FALSE,
                   check.names = FALSE)
  row.names(out) = NULL
  out
}

.plant_chemical_trait_matrix_long = function(occurrences, categorate_result,
                                             level, profile, min_confidence) {
  cols = c("id", "trait", "value")
  if (!is.list(categorate_result) ||
      !is.data.frame(categorate_result$ChemicalTraits) ||
      nrow(categorate_result$ChemicalTraits) < 1 ||
      nrow(occurrences) < 1) {
    return(.uaf_empty_table(cols))
  }
  traits = .normalized_prepare_traits(categorate_result$ChemicalTraits)
  traits = .normalized_filter_trait_matrix_source(
    chemical_traits = traits,
    profile = .plant_chemical_trait_profile(profile),
    min_confidence = .normalized_trait_confidence_threshold(min_confidence)
  )
  if (nrow(traits) < 1) return(.uaf_empty_table(cols))
  traits$query_clean = .plant_clean_compound(traits$Query)
  rows = list()
  for (i in seq_len(nrow(occurrences))) {
    occurrence = occurrences[i, , drop = FALSE]
    id = .uaf_first_non_empty_text(occurrence[[level]],
                                   occurrence$species)
    if (is.na(id) || id == "") next
    hits = traits[traits$query_clean == occurrence$compound_name_clean, ,
                  drop = FALSE]
    if (nrow(hits) < 1) next
    rows[[length(rows) + 1]] = data.frame(
      id = id,
      trait = paste0("chemtrait__", hits$MatrixKey),
      compound_name_clean = occurrence$compound_name_clean,
      value = suppressWarnings(as.numeric(hits$ConfidenceScore)),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = do.call(rbind, rows)
  out$value[is.na(out$value)] = 0
  out = unique(out[, c("id", "trait", "compound_name_clean", "value"),
                   drop = FALSE])
  out = out[, cols, drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_context_matrix_long = function(occurrences, level, group_col, prefix,
                                      exclude = c("unknown")) {
  cols = c("id", "trait", "value")
  if (!is.data.frame(occurrences) || nrow(occurrences) < 1 ||
      !group_col %in% names(occurrences)) {
    return(.uaf_empty_table(cols))
  }
  keep = !is.na(occurrences[[group_col]]) &
    occurrences[[group_col]] != "" &
    !occurrences[[group_col]] %in% exclude
  rows = occurrences[keep, , drop = FALSE]
  if (nrow(rows) < 1) return(.uaf_empty_table(cols))
  id = rows[[level]]
  id[is.na(id) | id == ""] = rows$species[is.na(id) | id == ""]
  value = suppressWarnings(as.numeric(rows$evidence_quality_score))
  missing_value = !is.finite(value)
  value[missing_value] = .plant_confidence_score(rows$confidence[missing_value])
  out = data.frame(
    id = id,
    trait = paste0(prefix, "__", .plant_matrix_key(rows[[group_col]])),
    value = value,
    stringsAsFactors = FALSE
  )
  out = unique(out)
  row.names(out) = NULL
  out
}

.plant_chemical_trait_profile = function(profile) {
  if (identical(profile, "metabolism")) return("kegg")
  profile
}

.plant_comparability_traits = function(categorate_result) {
  if (is.list(categorate_result) &&
      is.data.frame(categorate_result$ChemicalTraits)) {
    traits = categorate_result$ChemicalTraits
  } else {
    traits = .uaf_empty_table(c("Query", "TraitType", "TraitGroup",
                                "TraitValue", "Confidence",
                                "SourceDatabase", "PathwayGroup"))
  }
  traits$query_clean = if ("Query" %in% names(traits)) {
    .plant_clean_compound(traits$Query)
  } else {
    rep(NA_character_, nrow(traits))
  }
  traits
}

.plant_comparability_classes = function(categorate_result) {
  if (is.list(categorate_result) &&
      is.data.frame(categorate_result$ChemicalClasses)) {
    classes = categorate_result$ChemicalClasses
  } else {
    classes = .uaf_empty_table(c("Query", "ClassName", "ClassGroup",
                                 "Superclass", "Class", "Subclass"))
  }
  classes$query_clean = if ("Query" %in% names(classes)) {
    .plant_clean_compound(classes$Query)
  } else {
    rep(NA_character_, nrow(classes))
  }
  classes
}

.plant_comparability_class_sources = function(categorate_result) {
  cols = c("query_clean", "source", "source_table", "source_field",
           "source_value", "source_confidence", "priority")
  if (!is.list(categorate_result)) return(.uaf_empty_table(cols))
  rows = list()
  add_rows = function(table_name, table, query_col, source_field,
                     value_cols, confidence_col = NULL, priority = 50,
                     source = table_name) {
    if (!is.data.frame(table) || nrow(table) < 1 ||
        !query_col %in% names(table)) {
      return(invisible(NULL))
    }
    value_cols = intersect(value_cols, names(table))
    if (length(value_cols) < 1) return(invisible(NULL))
    for (i in seq_len(nrow(table))) {
      row = table[i, , drop = FALSE]
      source_value = .pubchem_collapse(unlist(row[value_cols],
                                              use.names = FALSE))
      if (length(.uaf_non_empty(source_value)) < 1) next
      source_table = .uaf_first_non_empty_text(
        if ("SourceTable" %in% names(row)) row$SourceTable else NA_character_,
        table_name
      )
      source_name = .uaf_first_non_empty_text(
        if (is.function(source)) source(row) else source,
        table_name
      )
      confidence = .uaf_first_non_empty_text(
        if (!is.null(confidence_col) && confidence_col %in% names(row)) {
          row[[confidence_col]]
        } else {
          NA_character_
        },
        NA_character_
      )
      rows[[length(rows) + 1]] <<- data.frame(
        query_clean = .plant_clean_compound(row[[query_col]]),
        source = source_name,
        source_table = source_table,
        source_field = source_field,
        source_value = source_value,
        source_confidence = confidence,
        priority = suppressWarnings(as.numeric(
          if (is.function(priority)) priority(row) else priority
        )),
        stringsAsFactors = FALSE
      )
    }
    invisible(NULL)
  }

  add_rows(
    "ChemicalClasses", categorate_result$ChemicalClasses, "Query",
    "ClassSystem/ClassType/ClassName/ClassGroup/ClassPath",
    c("ClassSystem", "ClassType", "ClassName", "ClassGroup", "ClassPath",
      "EvidenceText", "ExtractionRule"),
    confidence_col = "Confidence",
    priority = .plant_class_source_priority,
    source = function(row) paste(.uaf_non_empty(c(
      "ChemicalClasses",
      .uaf_first_non_empty_text(.plant_col_or_default(row, "ClassSystem",
                                                      NA_character_)),
      .uaf_first_non_empty_text(.plant_col_or_default(row, "ClassType",
                                                      NA_character_))
    )), collapse = ":")
  )
  add_rows(
    "LOTUSProfile", categorate_result$LOTUSProfile, "Query",
    "NaturalProductClass", c("NaturalProductClass", "RawValue"),
    priority = 96,
    source = "LOTUSProfile:NaturalProductClass"
  )
  add_rows(
    "PubChemClassifications", categorate_result$PubChemClassifications,
    "Query", "ClassPath/ClassName",
    c("Source", "TreeName", "TreeType", "ClassName", "ParentClass",
      "ClassPath"),
    priority = .plant_pubchem_class_source_priority,
    source = function(row) paste(.uaf_non_empty(c(
      "PubChemClassifications",
      .uaf_first_non_empty_text(.plant_col_or_default(row, "Source",
                                                      NA_character_)),
      .uaf_first_non_empty_text(.plant_col_or_default(row, "TreeName",
                                                      NA_character_))
    )), collapse = ":")
  )
  add_rows(
    "ChemicalTraitOntology", categorate_result$ChemicalTraitOntology,
    "Query", "OntologyDomain/OntologyGroup/OntologyLabel",
    c("OntologyDomain", "OntologyGroup", "OntologyTerm", "OntologyLabel",
      "SourceTraitType", "SourceTraitGroup", "SourceTraitValue",
      "SourceDatabase"),
    confidence_col = "Confidence",
    priority = .plant_ontology_source_priority,
    source = function(row) paste(.uaf_non_empty(c(
      "ChemicalTraitOntology",
      .uaf_first_non_empty_text(.plant_col_or_default(row, "OntologyDomain",
                                                      NA_character_)),
      .uaf_first_non_empty_text(.plant_col_or_default(row, "OntologyGroup",
                                                      NA_character_))
    )), collapse = ":")
  )
  add_rows(
    "ChemicalTerms", categorate_result$ChemicalTerms, "Query",
    "Domain/TermType/TermClean",
    c("Domain", "TermType", "TermRaw", "TermClean", "TermGroup",
      "SourceTable", "Source", "EvidenceText"),
    confidence_col = "Confidence",
    priority = .plant_terms_source_priority,
    source = function(row) paste(.uaf_non_empty(c(
      "ChemicalTerms",
      .uaf_first_non_empty_text(.plant_col_or_default(row, "Domain",
                                                      NA_character_)),
      .uaf_first_non_empty_text(.plant_col_or_default(row, "TermType",
                                                      NA_character_))
    )), collapse = ":")
  )
  add_rows(
    "KEGGClassifications", categorate_result$KEGGClassifications, "Query",
    "ClassificationType/Classification",
    c("ClassificationType", "Classification", "Evidence"),
    priority = 78,
    source = "KEGGClassifications"
  )
  add_rows(
    "KEGGPathways", categorate_result$KEGGPathways, "Query",
    "PathwayGroup/PathwayName",
    c("PathwayGroup", "PathwayName", "Evidence"),
    priority = 78,
    source = "KEGGPathways"
  )
  add_rows(
    "DerivedGroups", categorate_result$DerivedGroups, "Query",
    "derived_contexts",
    c("natural_product_superclasses", "natural_product_classes",
      "natural_product_subclasses", "metabolic_context",
      "ecological_context", "analytical_context", "volatility_proxy",
      "lipophilicity_bin", "polarity_bin"),
    priority = 62,
    source = "DerivedGroups"
  )

  out = .plant_bind_flexible_tables(rows)
  if (!is.data.frame(out) || nrow(out) < 1) return(.uaf_empty_table(cols))
  for (col in cols) if (!col %in% names(out)) out[[col]] = NA_character_
  out = out[, cols, drop = FALSE]
  out$priority = suppressWarnings(as.numeric(out$priority))
  out$priority[!is.finite(out$priority)] = 50
  out$source_confidence = .plant_comparability_clean_confidence(
    out$source_confidence
  )
  out = out[!is.na(out$query_clean) & out$query_clean != "" &
              !is.na(out$source_value) & out$source_value != "", ,
            drop = FALSE]
  out = unique(out)
  out = out[order(out$query_clean, -out$priority, out$source,
                  out$source_field), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_normalize_comparability = function(x) {
  if (!is.data.frame(x)) return(.uaf_empty_table(.plant_comparability_cols()))
  cols = .plant_comparability_cols()
  for (col in cols) if (!col %in% names(x)) x[[col]] = NA_character_
  x = x[, cols, drop = FALSE]
  char_cols = names(x)[!names(x) %in% "evidence_quality_score"]
  x[char_cols] = lapply(x[char_cols], .uaf_squish_text)
  x$species_slug[is.na(x$species_slug) | x$species_slug == ""] =
    .plant_slug(x$species[is.na(x$species_slug) | x$species_slug == ""])
  x$metabolism_domain[is.na(x$metabolism_domain) |
                        x$metabolism_domain == ""] = "unknown"
  x$biosynthetic_family[is.na(x$biosynthetic_family) |
                          x$biosynthetic_family == ""] = "unknown"
  x$chemical_behavior[is.na(x$chemical_behavior) |
                        x$chemical_behavior == ""] = "unknown"
  x$comparison_scope[is.na(x$comparison_scope) |
                       x$comparison_scope == ""] = "unknown"
  x$comparison_group[is.na(x$comparison_group) |
                       x$comparison_group == ""] = "unknown"
  x$comparison_subgroup[is.na(x$comparison_subgroup) |
                          x$comparison_subgroup == ""] = "unknown"
  x$comparability_confidence[is.na(x$comparability_confidence) |
                               x$comparability_confidence == ""] = "unknown"
  x$classification_source_confidence[
    is.na(x$classification_source_confidence) |
      x$classification_source_confidence == ""
  ] = "unknown"
  x$comparable_for_matrix[is.na(x$comparable_for_matrix) |
                            x$comparable_for_matrix == ""] = "No"
  row.names(x) = NULL
  x
}

.plant_comparability_unenriched_table = function(occurrences) {
  n = nrow(occurrences)
  if (n < 1L) {
    return(.uaf_empty_table(.plant_comparability_classification_cols()))
  }
  compound_value = .plant_row_collapse_values(list(
    occurrences$compound_name, occurrences$compound_name_clean
  ))
  evidence_value = .plant_row_collapse_values(list(
    occurrences$evidence_text, occurrences$occurrence_type,
    occurrences$plant_part, occurrences$tissue, occurrences$method,
    occurrences$method_group, occurrences$source_database
  ))
  compound_text = tolower(compound_value)
  evidence_text = tolower(evidence_value)

  family = rep("unknown", n)
  family_pattern = rep(NA_character_, n)
  family_source = rep("none", n)
  family_source_value = rep(NA_character_, n)
  family_patterns = .plant_comparability_family_patterns()
  for (family_name in names(family_patterns)) {
    patterns = family_patterns[[family_name]]
    combined = paste(patterns, collapse = "|")
    available = family == "unknown"
    compound_hit = available & grepl(combined, compound_text, perl = TRUE)
    compound_hit[is.na(compound_hit)] = FALSE
    if (any(compound_hit)) {
      family[compound_hit] = family_name
      family_source[compound_hit] = "compound_name_rule"
      family_source_value[compound_hit] = compound_value[compound_hit]
      family_pattern[compound_hit] = .plant_first_pattern_vector(
        compound_text[compound_hit], patterns
      )
    }
    available = family == "unknown"
    evidence_hit = available & grepl(combined, evidence_text, perl = TRUE)
    evidence_hit[is.na(evidence_hit)] = FALSE
    if (any(evidence_hit)) {
      family[evidence_hit] = family_name
      family_source[evidence_hit] = "occurrence_evidence_text_rule"
      family_source_value[evidence_hit] = evidence_value[evidence_hit]
      family_pattern[evidence_hit] = .plant_first_pattern_vector(
        evidence_text[evidence_hit], patterns
      )
    }
  }
  broad_patterns = c("\\bnatural product", "\\bmetabolite",
                     "\\bphytochemical", "\\bsecondary metabolite")
  broad_combined = paste(broad_patterns, collapse = "|")
  available = family == "unknown"
  compound_hit = available & grepl(broad_combined, compound_text, perl = TRUE)
  compound_hit[is.na(compound_hit)] = FALSE
  if (any(compound_hit)) {
    family[compound_hit] = "broad_or_uncertain"
    family_source[compound_hit] = "compound_name_rule"
    family_source_value[compound_hit] = compound_value[compound_hit]
    family_pattern[compound_hit] = .plant_first_pattern_vector(
      compound_text[compound_hit], broad_patterns
    )
  }
  available = family == "unknown"
  evidence_hit = available & grepl(broad_combined, evidence_text, perl = TRUE)
  evidence_hit[is.na(evidence_hit)] = FALSE
  if (any(evidence_hit)) {
    family[evidence_hit] = "broad_or_uncertain"
    family_source[evidence_hit] = "occurrence_evidence_text_rule"
    family_source_value[evidence_hit] = evidence_value[evidence_hit]
    family_pattern[evidence_hit] = .plant_first_pattern_vector(
      evidence_text[evidence_hit], broad_patterns
    )
  }

  family_basis = rep(NA_character_, n)
  family_matched = !is.na(family_pattern) & family_pattern != ""
  family_field = ifelse(family_source == "compound_name_rule",
                        "compound_name", "evidence_text/occurrence_context")
  family_basis[family_matched] = paste0(
    family_source[family_matched], " matched `",
    family_pattern[family_matched], "` in ", family_field[family_matched], "."
  )

  volatile_patterns = .plant_comparability_volatile_patterns()
  volatile_combined = paste(volatile_patterns, collapse = "|")
  behavior = rep("nonvolatile_or_unspecified", n)
  behavior_source = rep("default", n)
  behavior_pattern = rep(NA_character_, n)
  compound_volatile = grepl(volatile_combined, compound_text, perl = TRUE)
  compound_volatile[is.na(compound_volatile)] = FALSE
  if (any(compound_volatile)) {
    behavior[compound_volatile] = "volatile_semivolatile"
    behavior_source[compound_volatile] = "compound_name_rule"
    behavior_pattern[compound_volatile] = .plant_first_pattern_vector(
      compound_text[compound_volatile], volatile_patterns
    )
  }
  evidence_volatile = behavior_source == "default" &
    grepl(volatile_combined, evidence_text, perl = TRUE)
  evidence_volatile[is.na(evidence_volatile)] = FALSE
  if (any(evidence_volatile)) {
    behavior[evidence_volatile] = "volatile_semivolatile"
    behavior_source[evidence_volatile] = "occurrence_evidence_text_rule"
    behavior_pattern[evidence_volatile] = .plant_first_pattern_vector(
      evidence_text[evidence_volatile], volatile_patterns
    )
  }
  behavior_basis = rep(
    "No volatile or semivolatile evidence pattern matched.", n
  )
  behavior_matched = !is.na(behavior_pattern) & behavior_pattern != ""
  behavior_field = ifelse(behavior_source == "compound_name_rule",
                          "compound_name", "evidence_text/occurrence_context")
  behavior_basis[behavior_matched] = paste0(
    behavior_source[behavior_matched], " matched `",
    behavior_pattern[behavior_matched], "` in ",
    behavior_field[behavior_matched], "."
  )

  domain = rep("unknown", n)
  domain[family %in% c("carbohydrate", "amino_acid", "organic_acid",
                       "nucleoside_nucleotide")] = "primary_metabolism"
  domain[family == "fatty_acid_lipid"] = "lipid_metabolism"
  domain[family == "plant_hormone_signal"] = "plant_hormone_signaling"
  specialized = c(
    "terpenoid", "phenolic_phenylpropanoid", "flavonoid",
    "alkaloid_nitrogenous", "organosulfur", "glucosinolate",
    "benzoxazinoid", "cyanogenic_glycoside", "saponin",
    "steroid_triterpenoid", "polyketide"
  )
  domain[family %in% specialized] = "specialized_metabolism"
  domain[family == "broad_or_uncertain"] = "broad_or_uncertain"

  scope = rep("unknown", n)
  scope[domain == "broad_or_uncertain"] = "broad_or_uncertain"
  scope[domain == "primary_metabolism"] = "primary_metabolites"
  scope[domain == "lipid_metabolism"] = "lipids_fatty_acids"
  scope[domain == "plant_hormone_signaling"] = "plant_hormone_signaling"
  scope[domain == "specialized_metabolism"] = "specialized_metabolites"
  scope[domain == "specialized_metabolism" &
          behavior == "volatile_semivolatile"] =
    "volatile_specialized_metabolites"
  group = family
  volatile_scope = scope == "volatile_specialized_metabolites" &
    !family %in% c("unknown", "broad_or_uncertain")
  group[volatile_scope] = paste0("volatile_", family[volatile_scope])
  comparable = scope %in% c(
    "specialized_metabolites", "volatile_specialized_metabolites",
    "primary_metabolites", "lipids_fatty_acids",
    "plant_hormone_signaling"
  ) & !family %in% c("unknown", "broad_or_uncertain")

  comparability_confidence = rep("unknown", n)
  comparability_confidence[family_source == "compound_name_rule"] = "medium"
  comparability_confidence[
    family_source == "occurrence_evidence_text_rule"
  ] = "low"
  behavior_only = family_source == "none" &
    behavior_source %in% c("compound_name_rule",
                            "occurrence_evidence_text_rule")
  comparability_confidence[behavior_only] = "low"
  source_table = rep(NA_character_, n)
  source_field = rep(NA_character_, n)
  source_confidence = rep("unknown", n)
  source_backed = family_source != "none"
  source_table[source_backed] = "PlantCompoundOccurrences"
  source_field[family_source == "compound_name_rule"] = "compound_name"
  source_field[family_source == "occurrence_evidence_text_rule"] =
    "evidence_text/occurrence_context"
  source_confidence[family_source == "compound_name_rule"] = "medium"
  source_confidence[family_source == "occurrence_evidence_text_rule"] = "low"

  data.frame(
    metabolism_domain = domain,
    biosynthetic_family = family,
    chemical_behavior = behavior,
    comparison_scope = scope,
    comparison_group = group,
    comparison_subgroup = family,
    comparability_confidence = comparability_confidence,
    comparability_basis = .plant_row_collapse_values(
      list(family_basis, behavior_basis)
    ),
    classification_source = family_source,
    classification_source_table = source_table,
    classification_source_field = source_field,
    classification_source_value = family_source_value,
    classification_source_confidence = source_confidence,
    comparable_for_matrix = .uaf_yes_no(comparable),
    comparison_caveat = .plant_vector_comparison_caveat(scope, behavior),
    stringsAsFactors = FALSE
  )
}

.plant_comparability_classification_cols = function() {
  c("metabolism_domain", "biosynthetic_family", "chemical_behavior",
    "comparison_scope", "comparison_group", "comparison_subgroup",
    "comparability_confidence", "comparability_basis",
    "classification_source", "classification_source_table",
    "classification_source_field", "classification_source_value",
    "classification_source_confidence", "comparable_for_matrix",
    "comparison_caveat")
}

.plant_row_collapse_values = function(values) {
  if (!is.list(values) || length(values) < 1L) return(character())
  n = max(lengths(values))
  values = lapply(values, function(value) {
    rep(.uaf_squish_text(value), length.out = n)
  })
  out = rep(NA_character_, n)
  previous = list()
  for (value in values) {
    available = !is.na(value) & value != ""
    duplicate = rep(FALSE, n)
    if (length(previous) > 0L) {
      for (seen in previous) {
        duplicate = duplicate |
          (available & !is.na(seen) & seen != "" & value == seen)
      }
    }
    add = available & !duplicate
    first = add & (is.na(out) | out == "")
    out[first] = value[first]
    append = add & !first
    out[append] = paste(out[append], value[append], sep = "; ")
    previous[[length(previous) + 1L]] = value
  }
  out
}

.plant_first_pattern_vector = function(text, patterns) {
  out = rep(NA_character_, length(text))
  for (pattern in patterns) {
    hit = is.na(out) & grepl(pattern, text, perl = TRUE)
    hit[is.na(hit)] = FALSE
    out[hit] = pattern
  }
  out
}

.plant_vector_comparison_caveat = function(scope, behavior) {
  key = paste(scope, behavior, sep = "||")
  unique_key = unique(key)
  values = vapply(unique_key, function(item) {
    parts = strsplit(item, "||", fixed = TRUE)[[1]]
    .plant_comparison_caveat(parts[[1]], parts[[2]])
  }, character(1), USE.NAMES = FALSE)
  values[match(key, unique_key)]
}

.plant_comparability_classify = function(occurrence, traits, classes,
                                         class_sources) {
  compound_clean = occurrence$compound_name_clean[[1]]
  trait_hits = traits[traits$query_clean == compound_clean, , drop = FALSE]
  class_hits = classes[classes$query_clean == compound_clean, , drop = FALSE]
  source_hits = class_sources[class_sources$query_clean == compound_clean, ,
                              drop = FALSE]
  text_sources = .plant_comparability_text_sources(occurrence, trait_hits,
                                                   class_hits, source_hits)
  family_hit = .plant_comparability_family_hit(text_sources)
  behavior_hit = .plant_comparability_behavior_hit(text_sources)
  family = family_hit$value
  behavior = behavior_hit$value
  domain = .plant_metabolism_domain_for_family(family)
  scope = .plant_comparison_scope_for_family(family, behavior, domain)
  group = .plant_comparison_group_for_scope(family, scope)
  comparable = scope %in% c("specialized_metabolites",
                            "volatile_specialized_metabolites",
                            "primary_metabolites", "lipids_fatty_acids",
                            "plant_hormone_signaling") &&
    !family %in% c("unknown", "broad_or_uncertain")
  confidence = .plant_comparability_confidence(family_hit, behavior_hit)
  basis = .pubchem_collapse(c(family_hit$basis, behavior_hit$basis))
  if (length(.uaf_non_empty(basis)) < 1) {
    basis = "No source-backed chemistry scope pattern was matched."
  }
  list(
    metabolism_domain = domain,
    biosynthetic_family = family,
    chemical_behavior = behavior,
    comparison_scope = scope,
    comparison_group = group,
    comparison_subgroup = .plant_comparison_subgroup(family, trait_hits,
                                                     class_hits),
    comparability_confidence = confidence,
    comparability_basis = basis,
    classification_source = .uaf_first_non_empty_text(family_hit$source,
                                                      NA_character_),
    classification_source_table =
      .uaf_first_non_empty_text(family_hit$source_table, NA_character_),
    classification_source_field =
      .uaf_first_non_empty_text(family_hit$source_field, NA_character_),
    classification_source_value =
      .uaf_first_non_empty_text(family_hit$source_value, NA_character_),
    classification_source_confidence =
      .plant_comparability_clean_confidence(family_hit$source_confidence),
    comparable_for_matrix = .uaf_yes_no(comparable),
    comparison_caveat = .plant_comparison_caveat(scope, behavior)
  )
}

.plant_comparability_text_sources = function(occurrence, traits, classes,
                                             source_hits) {
  cols = c("source_key", "source", "source_table", "source_field",
           "source_value", "source_confidence", "priority", "text")
  rows = list()
  if (is.data.frame(source_hits) && nrow(source_hits) > 0) {
    for (i in seq_len(nrow(source_hits))) {
      hit = source_hits[i, , drop = FALSE]
      rows[[length(rows) + 1]] = data.frame(
        source_key = "source_class",
        source = hit$source,
        source_table = hit$source_table,
        source_field = hit$source_field,
        source_value = hit$source_value,
        source_confidence = hit$source_confidence,
        priority = suppressWarnings(as.numeric(hit$priority)),
        text = tolower(hit$source_value),
        stringsAsFactors = FALSE
      )
    }
  }
  trait_class_value = .pubchem_collapse(c(
      .plant_col_or_default(traits, "TraitType", NA_character_),
      .plant_col_or_default(traits, "TraitGroup", NA_character_),
      .plant_col_or_default(traits, "TraitValue", NA_character_),
      .plant_col_or_default(traits, "SourceDatabase", NA_character_),
      .plant_col_or_default(traits, "PathwayGroup", NA_character_),
      .plant_col_or_default(classes, "ClassName", NA_character_),
      .plant_col_or_default(classes, "ClassGroup", NA_character_),
      .plant_col_or_default(classes, "Superclass", NA_character_),
      .plant_col_or_default(classes, "Class", NA_character_),
      .plant_col_or_default(classes, "Subclass", NA_character_)
  ))
  if (length(.uaf_non_empty(trait_class_value)) > 0) {
    rows[[length(rows) + 1]] = data.frame(
      source_key = "trait_class",
      source = "ChemicalTraits/ChemicalClasses",
      source_table = "ChemicalTraits; ChemicalClasses",
      source_field = "TraitType/TraitGroup/TraitValue/ClassName/ClassGroup",
      source_value = trait_class_value,
      source_confidence = .plant_comparability_clean_confidence(
        .uaf_first_non_empty_text(
          .plant_col_or_default(traits, "Confidence", NA_character_),
          .plant_col_or_default(classes, "Confidence", NA_character_),
          "medium"
        )
      ),
      priority = 56,
      text = tolower(trait_class_value),
      stringsAsFactors = FALSE
    )
  }
  compound_value = .pubchem_collapse(c(occurrence$compound_name,
                                       occurrence$compound_name_clean))
  rows[[length(rows) + 1]] = data.frame(
    source_key = "compound_name",
    source = "compound_name_rule",
    source_table = "PlantCompoundOccurrences",
    source_field = "compound_name",
    source_value = compound_value,
    source_confidence = "medium",
    priority = 36,
    text = tolower(compound_value),
    stringsAsFactors = FALSE
  )
  evidence_value = .pubchem_collapse(c(
    occurrence$evidence_text, occurrence$occurrence_type,
    occurrence$plant_part, occurrence$tissue, occurrence$method,
    occurrence$method_group, occurrence$source_database
  ))
  if (length(.uaf_non_empty(evidence_value)) > 0) {
    rows[[length(rows) + 1]] = data.frame(
      source_key = "occurrence_evidence",
      source = "occurrence_evidence_text_rule",
      source_table = "PlantCompoundOccurrences",
      source_field = "evidence_text/occurrence_context",
      source_value = evidence_value,
      source_confidence = "low",
      priority = 24,
      text = tolower(evidence_value),
      stringsAsFactors = FALSE
    )
  }
  out = .plant_bind_tables(rows, cols)
  if (nrow(out) < 1) return(.uaf_empty_table(cols))
  out$priority = suppressWarnings(as.numeric(out$priority))
  out$priority[!is.finite(out$priority)] = 0
  out = out[order(-out$priority, out$source_key), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_comparability_family_hit = function(text_sources) {
  patterns = .plant_comparability_family_patterns()
  for (family in names(patterns)) {
    hit = .plant_comparability_pattern_hit(text_sources, patterns[[family]])
    if (!is.null(hit)) {
      hit$value = family
      return(hit)
    }
  }
  broad = .plant_comparability_pattern_hit(
    text_sources,
    c("\\bnatural product", "\\bmetabolite", "\\bphytochemical",
      "\\bsecondary metabolite")
  )
  if (!is.null(broad)) {
    broad$value = "broad_or_uncertain"
    return(broad)
  }
  list(value = "unknown", basis = NA_character_, source = "none",
       source_table = NA_character_, source_field = NA_character_,
       source_value = NA_character_, source_confidence = "unknown",
       priority = 0)
}

.plant_comparability_behavior_hit = function(text_sources) {
  volatile_patterns = .plant_comparability_volatile_patterns()
  hit = .plant_comparability_pattern_hit(text_sources, volatile_patterns)
  if (!is.null(hit)) {
    hit$value = "volatile_semivolatile"
    return(hit)
  }
  list(value = "nonvolatile_or_unspecified",
       basis = "No volatile or semivolatile evidence pattern matched.",
       source = "default")
}

.plant_comparability_volatile_patterns = function() {
  c(
    "\\bvolatile", "\\bsemi[- ]?volatile", "\\bessential oil",
    "\\baroma", "\\bodou?r", "\\bolfactory", "\\bgc[- ]?ms",
    "\\bmonoterp", "\\bsesquiterp", "\\bpinene", "\\blimonene",
    "\\blinalool", "\\bterpin", "\\bcineol", "\\bcineole",
    "\\beucalyptol",
    "\\bcamphor", "\\bcarvone", "\\bmenthol", "\\bgeraniol",
    "\\bcitral", "\\bthymol", "\\bcarvacrol"
  )
}

.plant_comparability_pattern_hit = function(text_sources, patterns) {
  if (!is.data.frame(text_sources) || nrow(text_sources) < 1) return(NULL)
  text_sources = text_sources[order(-suppressWarnings(as.numeric(
    text_sources$priority
  )), text_sources$source_key), , drop = FALSE]
  for (i in seq_len(nrow(text_sources))) {
    row = text_sources[i, , drop = FALSE]
    text = row$text[[1]]
    if (is.na(text) || text == "") next
    matched = patterns[vapply(patterns, grepl, logical(1), x = text,
                              perl = TRUE)]
    if (length(matched) > 0) {
      return(list(
        value = NA_character_,
        basis = paste0(row$source[[1]], " matched `", matched[[1]],
                       "` in ", row$source_field[[1]], "."),
        source = row$source[[1]],
        source_table = row$source_table[[1]],
        source_field = row$source_field[[1]],
        source_value = row$source_value[[1]],
        source_confidence = row$source_confidence[[1]],
        priority = suppressWarnings(as.numeric(row$priority[[1]]))
      ))
    }
  }
  NULL
}

.plant_class_source_priority = function(row) {
  system = tolower(.uaf_first_non_empty_text(
    .plant_col_or_default(row, "ClassSystem", NA_character_), ""
  ))
  type = tolower(.uaf_first_non_empty_text(
    .plant_col_or_default(row, "ClassType", NA_character_), ""
  ))
  source_table = tolower(.uaf_first_non_empty_text(
    .plant_col_or_default(row, "SourceTable", NA_character_), ""
  ))
  if (grepl("lotus", system) && grepl("natural_product", type)) return(100)
  if (grepl("lotus", source_table) && grepl("natural_product", type)) {
    return(98)
  }
  if (grepl("kegg", system) && grepl("pathway|enzyme", type)) return(82)
  if (grepl("mesh", system)) return(54)
  72
}

.plant_pubchem_class_source_priority = function(row) {
  text = tolower(.pubchem_collapse(c(
    .plant_col_or_default(row, "Source", NA_character_),
    .plant_col_or_default(row, "TreeName", NA_character_),
    .plant_col_or_default(row, "TreeType", NA_character_),
    .plant_col_or_default(row, "ClassPath", NA_character_)
  )))
  if (grepl("lotus|natural product|npclassifier|classyfire", text)) {
    return(94)
  }
  if (grepl("chemical", text)) return(84)
  66
}

.plant_ontology_source_priority = function(row) {
  text = tolower(.pubchem_collapse(c(
    .plant_col_or_default(row, "OntologyDomain", NA_character_),
    .plant_col_or_default(row, "OntologyGroup", NA_character_),
    .plant_col_or_default(row, "SourceTraitType", NA_character_),
    .plant_col_or_default(row, "SourceTraitGroup", NA_character_)
  )))
  if (grepl("natural_product|chemical_class", text)) return(88)
  if (grepl("pathway_group|enzyme_class|metabolism", text)) return(80)
  if (grepl("physicochemical|sensory", text)) return(64)
  48
}

.plant_terms_source_priority = function(row) {
  text = tolower(.pubchem_collapse(c(
    .plant_col_or_default(row, "Domain", NA_character_),
    .plant_col_or_default(row, "TermType", NA_character_),
    .plant_col_or_default(row, "TermGroup", NA_character_),
    .plant_col_or_default(row, "SourceTable", NA_character_)
  )))
  if (grepl("classification|natural_product|chemical_class", text)) return(72)
  if (grepl("pathway|metabolism", text)) return(66)
  if (grepl("sensory|analytical|physicochemical", text)) return(58)
  38
}

.plant_comparability_clean_confidence = function(x) {
  x = tolower(.uaf_squish_text(x))
  out = ifelse(x %in% c("high", "medium", "low", "unknown"), x, NA_character_)
  out[is.na(out)] = "unknown"
  out
}

.plant_comparability_family_patterns = function() {
  list(
    plant_hormone_signal = c(
      "\\babscisic acid\\b", "\\bjasmonic acid\\b", "\\bsalicylic acid\\b",
      "\\bgibberellin", "\\bauxin", "\\bindole[- ]?3[- ]?acetic",
      "\\bcytokinin", "\\bbrassinosteroid", "\\bethylene\\b",
      "\\bplant hormone", "\\bsignaling hormone"
    ),
    glucosinolate = c("\\bglucosinolate", "\\bisothiocyanate",
                      "\\bthiocyanate"),
    cyanogenic_glycoside = c("\\bcyanogenic", "\\bamygdalin",
                             "\\blinamarin", "\\bdhurrin"),
    benzoxazinoid = c("\\bbenzoxaz", "\\bdimboa\\b", "\\bmboa\\b"),
    saponin = c("\\bsaponin", "\\bsapogenin"),
    steroid_triterpenoid = c("\\bsteroid", "\\bsterol", "\\btriterpen",
                             "\\bsteroids and steroid derivatives"),
    polyketide = c("\\bpolyketide", "\\bmacrolide", "\\banthraquinone",
                   "\\btetracycline", "\\btype i polyketide",
                   "\\btype ii polyketide"),
    terpenoid = c(
      "\\bterpen", "\\bmonoterp", "\\bsesquiterp", "\\bditerp",
      "\\btriterp", "\\bpinene", "\\blimonene", "\\blinalool",
      "\\bterpin", "\\bcineol", "\\bcineole", "\\beucalyptol", "\\bcamphor",
      "\\bcarvone", "\\bmenthol", "\\bmenthone", "\\bgeraniol",
      "\\bcitral", "\\bthymol", "\\bcarvacrol", "\\bisopren",
      "\\bprenol lipid", "\\bisoprenoid", "\\bmeroterpenoid"
    ),
    flavonoid = c("\\bflavonoid", "\\bflavone", "\\bflavonol",
                  "\\banthocyan", "\\bcatechin", "\\bquercetin",
                  "\\bkaempferol", "\\bnaringenin", "\\brutin\\b",
                  "\\bflavan", "\\bchalcone", "\\bisoflav"),
    phenolic_phenylpropanoid = c(
      "\\bphenolic", "\\bphenol", "\\bphenylpropanoid",
      "\\bshikimate", "\\bshikimates and phenylpropanoids",
      "\\bcaffeic acid\\b", "\\bcoumaric acid\\b", "\\bferulic acid\\b",
      "\\bchlorogenic acid\\b", "\\btannin", "\\blignin",
      "\\bsalicin\\b", "\\bsalicyl", "\\bbenzoic acid\\b",
      "\\bcinnamic acid\\b", "\\bcoumarin", "\\blignan",
      "\\bstilbene", "\\bhydroxycinnamic", "\\bhydroxybenzoic",
      "\\bsimple phenolic acid"
    ),
    alkaloid_nitrogenous = c(
      "\\balkaloid", "\\bcaffeine\\b", "\\bnicotine\\b",
      "\\btheobromine\\b", "\\btheophylline\\b", "\\bindole alkaloid",
      "\\bquinoline", "\\bpyridine", "\\bpiperidine"
    ),
    organosulfur = c("\\borganosulfur", "\\balliin\\b", "\\ballicin\\b",
                     "\\bajoene\\b", "\\bisoalliin\\b", "\\bpropiin\\b",
                     "\\bsulfide", "\\bsulfoxide", "\\bsulfinyl",
                     "\\bsulfanyl", "\\bsulfen", "\\bdisulfen",
                     "\\bthiol", "\\bthio", "\\bdithio",
                     "\\bthiosulfinate", "\\bsulfur"),
    fatty_acid_lipid = c(
      "\\bfatty acid", "\\blipid", "\\blinoleic", "\\blinolenic",
      "\\boleic", "\\bpalmitic", "\\bstearic", "\\bphospholipid",
      "\\btriacyl", "\\btriglyceride", "\\bwax ester",
      "\\bfatty acyl", "\\bglycerolipid", "\\bglycerophospholipid",
      "\\bsphingolipid", "\\bsaccharolipid", "\\blipid metabolism",
      "\\bfatty acid biosynthesis", "\\bfatty acid degradation"
    ),
    carbohydrate = c("\\bcarbohydrate", "\\bsaccharide", "\\bsugar",
                     "\\bglucose\\b", "\\bfructose\\b", "\\bsucrose\\b",
                     "\\bcellulose\\b", "\\bstarch\\b", "\\bmannitol\\b",
                     "\\bcarbohydrate metabolism", "\\bglycolysis",
                     "\\bstarch and sucrose metabolism"),
    amino_acid = c(
      "\\bamino acid", "\\balanine\\b", "\\barginine\\b",
      "\\basparagine\\b", "\\baspartic acid\\b", "\\bcysteine\\b",
      "\\bglutamic acid\\b", "\\bglutamine\\b", "\\bglycine\\b",
      "\\bhistidine\\b", "\\bisoleucine\\b", "\\bleucine\\b",
      "\\blysine\\b", "\\bmethionine\\b", "\\bphenylalanine\\b",
      "\\bproline\\b", "\\bserine\\b", "\\bthreonine\\b",
      "\\btryptophan\\b", "\\btyrosine\\b", "\\bvaline\\b",
      "\\bamino acids peptides and analogues",
      "\\bamino acid metabolism", "\\bbiosynthesis of amino acids"
    ),
    organic_acid = c("\\borganic acid", "\\bcitric acid\\b",
                     "\\bmalic acid\\b", "\\boxalic acid\\b",
                     "\\bsuccinic acid\\b", "\\bfumaric acid\\b",
                     "\\bpyruvic acid\\b", "\\blactic acid\\b",
                     "\\btricarboxylic acid", "\\bcitrate cycle",
                     "\\btca cycle", "\\bdicarboxylic acid",
                     "\\bhydroxy acid", "\\bketo acid"),
    nucleoside_nucleotide = c("\\bnucleoside", "\\bnucleotide",
                              "\\badenosine\\b", "\\bguanosine\\b",
                              "\\buridine\\b", "\\bcytidine\\b")
  )
}

.plant_metabolism_domain_for_family = function(family) {
  if (family %in% c("carbohydrate", "amino_acid", "organic_acid",
                    "nucleoside_nucleotide")) {
    return("primary_metabolism")
  }
  if (family == "fatty_acid_lipid") return("lipid_metabolism")
  if (family == "plant_hormone_signal") return("plant_hormone_signaling")
  if (family %in% c("terpenoid", "phenolic_phenylpropanoid", "flavonoid",
                    "alkaloid_nitrogenous", "organosulfur",
                    "glucosinolate", "benzoxazinoid",
                    "cyanogenic_glycoside", "saponin",
                    "steroid_triterpenoid", "polyketide")) {
    return("specialized_metabolism")
  }
  if (family == "broad_or_uncertain") return("broad_or_uncertain")
  "unknown"
}

.plant_comparison_scope_for_family = function(family, behavior, domain) {
  if (family == "unknown" || domain == "unknown") return("unknown")
  if (family == "broad_or_uncertain" || domain == "broad_or_uncertain") {
    return("broad_or_uncertain")
  }
  if (domain == "primary_metabolism") return("primary_metabolites")
  if (domain == "lipid_metabolism") return("lipids_fatty_acids")
  if (domain == "plant_hormone_signaling") {
    return("plant_hormone_signaling")
  }
  if (domain == "specialized_metabolism" &&
      behavior == "volatile_semivolatile") {
    return("volatile_specialized_metabolites")
  }
  if (domain == "specialized_metabolism") return("specialized_metabolites")
  "unknown"
}

.plant_comparison_group_for_scope = function(family, scope) {
  if (family %in% c("unknown", "broad_or_uncertain")) return(family)
  if (scope == "volatile_specialized_metabolites") {
    return(paste0("volatile_", family))
  }
  family
}

.plant_comparison_subgroup = function(family, traits, classes) {
  terms = .uaf_non_empty(c(
    .plant_col_or_default(classes, "Subclass", NA_character_),
    .plant_col_or_default(classes, "ClassName", NA_character_),
    .plant_col_or_default(traits, "TraitValue", NA_character_)
  ))
  if (length(terms) > 0) return(terms[[1]])
  family
}

.plant_comparability_confidence = function(family_hit, behavior_hit) {
  source = family_hit$source
  priority = suppressWarnings(as.numeric(family_hit$priority))
  source_confidence = .plant_comparability_clean_confidence(
    family_hit$source_confidence
  )
  if (!identical(source, "none") &&
      !source %in% c("compound_name_rule", "occurrence_evidence_text_rule") &&
      is.finite(priority) && priority >= 70) {
    if (source_confidence %in% c("high", "medium", "low")) {
      return(source_confidence)
    }
    return(ifelse(priority >= 80, "high", "medium"))
  }
  if (identical(source, "ChemicalTraits/ChemicalClasses")) {
    if (source_confidence %in% c("high", "medium", "low")) {
      return(source_confidence)
    }
    return("medium")
  }
  if (identical(source, "compound_name_rule")) return("medium")
  if (identical(source, "occurrence_evidence_text_rule")) return("low")
  if (!is.null(behavior_hit$source) &&
      behavior_hit$source %in% c("compound_name_rule",
                                 "occurrence_evidence_text_rule")) {
    return("low")
  }
  "unknown"
}

.plant_comparison_caveat = function(scope, behavior) {
  if (scope == "volatile_specialized_metabolites") {
    return(paste("Volatility describes an analytical or physicochemical",
                 "fraction. Compare these rows only against the same",
                 "volatile-specialized scope or a declared volatile subfamily."))
  }
  if (scope == "specialized_metabolites") {
    return(paste("Specialized-metabolite comparisons should not be mixed with",
                 "primary metabolites, broad lipid pools, or hormone signals",
                 "without an explicit mixed-scope analysis plan."))
  }
  if (scope == "primary_metabolites") {
    return(paste("Primary-metabolite rows support primary-metabolism",
                 "comparisons, not direct comparison against specialized",
                 "defense, aroma, or signaling chemistry."))
  }
  if (scope == "lipids_fatty_acids") {
    return(paste("Lipid and fatty-acid rows are separated because they can",
                 "reflect structural, storage, cuticular, or signaling pools."))
  }
  if (scope == "plant_hormone_signaling") {
    return(paste("Plant hormone and signaling rows should be interpreted as",
                 "signaling chemistry, not general primary or secondary",
                 "metabolite abundance."))
  }
  if (behavior == "volatile_semivolatile") {
    return(paste("A volatile signal was detected but no reliable chemistry",
                 "family was assigned; review before matrix use."))
  }
  "Chemistry scope is broad, uncertain, or unknown; review before comparison."
}

.plant_summary_comparability_info = function(group, categorate_result,
                                             comparability = NULL) {
  if (is.null(comparability) || !is.data.frame(comparability)) {
    comparability = plantChemistryComparability(
      list(PlantCompoundOccurrences = group,
           CategorateResult = categorate_result),
      min_confidence = "low"
    )
  }
  if (nrow(comparability) < 1) {
    return(list(specialized_count = 0, volatile_count = 0,
                primary_count = 0, lipid_count = 0, unknown_count = 0,
                comparison_scopes = NA_character_,
                dominant_groups = NA_character_))
  }
  unique_rows = unique(comparability[, c("compound_name_clean",
                                         "comparison_scope",
                                         "comparison_group",
                                         "comparable_for_matrix"),
                                      drop = FALSE])
  count_scope = function(scope) {
    length(unique(.uaf_non_empty(unique_rows$compound_name_clean[
      unique_rows$comparison_scope %in% scope
    ])))
  }
  list(
    specialized_count = count_scope(c("specialized_metabolites",
                                      "volatile_specialized_metabolites")),
    volatile_count = count_scope("volatile_specialized_metabolites"),
    primary_count = count_scope("primary_metabolites"),
    lipid_count = count_scope("lipids_fatty_acids"),
    unknown_count = count_scope(c("unknown", "broad_or_uncertain")),
    comparison_scopes = .plant_count_summary(unique_rows$comparison_scope),
    dominant_groups = .plant_top_terms(unique_rows$comparison_group[
      unique_rows$comparable_for_matrix == "Yes"
    ])
  )
}

.plant_export_tables = function(x, tables, include_empty, max_cell_chars) {
  default = c("PlantQueries", "PlantQueryAliases", "PlantNameResolution",
              "ProviderDiagnostics", "ProviderQueryAccounting",
              "ProviderResourceManifest",
              "PlantCompoundOccurrences", "PlantContextEvidence",
              "ProviderContextAudit",
              "LiteratureCandidates",
              "SourceCompoundIdentity",
              "CompoundResolution", "CompoundIdentityReview",
              "SpeciesChemistrySummary",
              "SpeciesChemistryMatrix", "ChemistryComparability",
              "ComparableChemistryMatrix", "TraitEvidence",
              "ValidationSummary", "ValidationIssues", "DataDictionary",
              "Provenance",
              "BatchRunManifest", "BatchChunkManifest",
              "FailedQueries", "RetryQueue",
              "BatchExportManifest")
  available = list(
    PlantQueries = x$PlantQueries,
    PlantQueryAliases = x$PlantQueryAliases,
    PlantNameResolution = x$PlantNameResolution,
    ProviderDiagnostics = x$ProviderDiagnostics,
    ProviderQueryAccounting = x$ProviderQueryAccounting,
    ProviderResourceManifest = x$ProviderResourceManifest,
    PlantCompoundOccurrences = .plant_clean_context_conflicts(
      x$PlantCompoundOccurrences
    ),
    PlantContextEvidence = x$PlantContextEvidence,
    ProviderContextAudit = x$ProviderContextAudit,
    LiteratureCandidates = x$LiteratureCandidates,
    SourceCompoundIdentity = x$SourceCompoundIdentity,
    CompoundResolution = x$CompoundResolution,
    CompoundIdentityReview = x$CompoundIdentityReview,
    SpeciesChemistrySummary = x$SpeciesChemistrySummary,
    SpeciesChemistryMatrix = x$SpeciesChemistryMatrix,
    ChemistryComparability = x$ChemistryComparability,
    ComparableChemistryMatrix = x$ComparableChemistryMatrix,
    TraitEvidence = x$TraitEvidence,
    ValidationSummary = if (is.list(x$Validation)) x$Validation$Summary else NULL,
    ValidationIssues = if (is.list(x$Validation)) x$Validation$Issues else NULL,
    DataDictionary = x$DataDictionary,
    Provenance = x$Provenance,
    BatchRunManifest = x$BatchRunManifest,
    BatchChunkManifest = x$BatchChunkManifest,
    FailedQueries = x$FailedQueries,
    RetryQueue = x$RetryQueue,
    BatchExportManifest = x$BatchExportManifest
  )
  requested = .uaf_non_empty(tables)
  table_names = if (length(requested) > 0) requested else default
  out = list()
  for (name in table_names) {
    table = available[[name]]
    if (!is.data.frame(table)) next
    if (!isTRUE(include_empty) && nrow(table) < 1) next
    out[[name]] = .categorate_export_prepare_table(table, max_cell_chars)
  }
  out
}

.plant_apply_optional_filter = function(keep, values, allowed) {
  allowed = .uaf_non_empty(allowed)
  if (length(allowed) < 1 || "all" %in% tolower(allowed)) return(keep)
  keep & values %in% allowed
}

.plant_review_cols = function() {
  c("review_id", "review_decision", "reviewed_by", "reviewed_at",
    "review_note", "proposed_evidence_tier", "proposed_confidence",
    "proposed_source_database", "proposed_citation_or_url",
    "proposed_plant_part", "proposed_tissue", "proposed_method",
    "review_key", .plant_occurrence_cols())
}

.plant_occurrence_key = function(occurrences) {
  dup_cols = .plant_duplicate_key_cols("PlantCompoundOccurrences")
  key_data = occurrences[dup_cols]
  key_data[] = lapply(key_data, function(value) {
    value = .uaf_squish_text(value)
    value[is.na(value)] = ""
    value
  })
  do.call(paste, c(key_data, sep = "||"))
}

.plant_review_match_indices = function(occurrences, review) {
  keys = .plant_occurrence_key(occurrences)
  review_key = if ("review_key" %in% names(review)) {
    as.character(review$review_key[[1]])
  } else {
    NA_character_
  }
  if (is.na(review_key) || review_key == "") review_key = NA_character_
  if (!is.na(review_key)) return(which(keys == review_key))
  dup_cols = .plant_duplicate_key_cols("PlantCompoundOccurrences")
  if (!all(dup_cols %in% names(review))) return(integer())
  review_data = review[dup_cols]
  review_data[] = lapply(review_data, function(value) {
    value = .uaf_squish_text(value)
    value[is.na(value)] = ""
    value
  })
  key = do.call(paste, c(review_data, sep = "||"))
  which(keys == key)
}

.plant_review_decision = function(x) {
  decision = tolower(.uaf_first_non_empty_text(x))
  if (is.na(decision)) return("needs_review")
  decision = gsub("[^a-z0-9]+", "_", decision)
  decision = gsub("^_|_$", "", decision)
  if (decision %in% c("accept", "promote", "promote_reported",
                      "accept_curated", "curated")) {
    return("promote_curated")
  }
  if (decision %in% c("remove", "drop")) return("reject")
  if (decision %in% c("update", "context", "edit_context")) {
    return("update_context")
  }
  if (decision %in% c("candidate", "keep_as_candidate")) {
    return("keep_candidate")
  }
  decision
}

.plant_apply_review_context = function(occurrences, idx, review) {
  update = function(col, value) {
    value = .uaf_first_non_empty_text(value)
    if (!is.na(value) && col %in% names(occurrences)) {
      occurrences[[col]][idx] <<- value
      TRUE
    } else {
      FALSE
    }
  }
  if (update("plant_part", review$proposed_plant_part) &&
      "plant_part_group" %in% names(occurrences)) {
    occurrences$plant_part_group[idx] = NA_character_
  }
  if (update("tissue", review$proposed_tissue) &&
      "tissue_group" %in% names(occurrences)) {
    occurrences$tissue_group[idx] = NA_character_
  }
  if (update("method", review$proposed_method) &&
      "method_group" %in% names(occurrences)) {
    occurrences$method_group[idx] = NA_character_
  }
  update("confidence", review$proposed_confidence)
  update("source_database", review$proposed_source_database)
  update("evidence_url", review$proposed_citation_or_url)
  update("evidence_tier", review$proposed_evidence_tier)
  occurrences
}

.plant_rebuild_filtered_result = function(x, occurrences) {
  occurrences = .plant_normalize_occurrences(occurrences)
  context_evidence = plantContextEvidence(occurrences)
  occurrences = .plant_apply_context_evidence(occurrences, context_evidence)
  out = x
  out$PlantCompoundOccurrences = occurrences
  out$PlantContextEvidence = context_evidence
  out$ProviderContextAudit = plantProviderContextAudit(
    list(PlantCompoundOccurrences = occurrences,
         PlantContextEvidence = context_evidence)
  )
  out$CompoundResolution = .plant_compound_resolution(occurrences,
                                                      x$CategorateResult)
  out$TraitEvidence = .plant_trait_evidence(occurrences, x$CategorateResult)
  out$ChemistryComparability = plantChemistryComparability(
    list(PlantCompoundOccurrences = occurrences,
         CategorateResult = x$CategorateResult),
    min_confidence = "low"
  )
  out$SpeciesChemistrySummary = summarizePlantPhytochemistry(
    plant_compounds = occurrences,
    categorate_result = x$CategorateResult,
    compound_resolution = out$CompoundResolution,
    plant_queries = x$PlantQueries,
    provider_diagnostics = x$ProviderDiagnostics,
    comparability = out$ChemistryComparability
  )
  out$SpeciesChemistryMatrix = plantPhytochemistryMatrix(
    list(PlantCompoundOccurrences = occurrences,
         CategorateResult = x$CategorateResult),
    level = "species",
    profile = "core",
    mode = "binary",
    min_confidence = "medium",
    max_traits = .plant_automatic_matrix_max_traits()
  )
  out$ComparableChemistryMatrix = plantComparableChemistryMatrix(
    list(ChemistryComparability = out$ChemistryComparability),
    comparison_scope = "specialized_metabolites",
    level = "species",
    mode = "binary",
    min_comparability_confidence = "medium"
  )
  out$Validation = validatePlantPhytochemistryResult(out)
  class(out) = unique(c("uaf_plant_phytochemistry", class(x)))
  out
}

.plant_validation_result = function(table_quality, provider_diagnostics,
                                    issues, dictionary) {
  summary = .plant_validation_summary(table_quality, issues)
  out = list(
    Summary = summary,
    TableQuality = table_quality,
    ProviderDiagnostics = provider_diagnostics,
    Issues = issues,
    DataDictionary = dictionary
  )
  class(out) = c("uaf_plant_phytochemistry_validation", class(out))
  out
}

.plant_validation_summary = function(table_quality, issues) {
  errors = sum(issues$severity == "error")
  warnings = sum(issues$severity == "warning")
  data.frame(
    Status = ifelse(errors > 0, "fail", ifelse(warnings > 0, "warning", "pass")),
    TablesExpected = nrow(table_quality),
    TablesPresent = sum(table_quality$Present == "Yes"),
    IssueCount = nrow(issues),
    ErrorCount = errors,
    WarningCount = warnings,
    stringsAsFactors = FALSE
  )
}

.plant_table_quality = function(table_name, table, dictionary, strict) {
  if (table_name == "Validation" && is.list(table)) {
    missing_required = setdiff(dictionary$Column[dictionary$Required],
                               names(table))
    return(data.frame(
      Table = table_name,
      Present = "Yes",
      RowCount = NA_integer_,
      ColumnCount = length(table),
      MissingRequiredColumns = .pubchem_collapse(missing_required),
      MissingOptionalColumns = "",
      DuplicateKeyColumns = "",
      Status = ifelse(length(missing_required) > 0, "fail", "pass"),
      stringsAsFactors = FALSE
    ))
  }
  present = is.data.frame(table)
  required = dictionary$Column[dictionary$Required]
  missing_required = if (present) setdiff(required, names(table)) else required
  optional = dictionary$Column[!dictionary$Required]
  missing_optional = if (present && isTRUE(strict)) {
    setdiff(optional, names(table))
  } else character()
  data.frame(
    Table = table_name,
    Present = .uaf_yes_no(present),
    RowCount = if (present) nrow(table) else NA_integer_,
    ColumnCount = if (present) ncol(table) else NA_integer_,
    MissingRequiredColumns = .pubchem_collapse(missing_required),
    MissingOptionalColumns = .pubchem_collapse(missing_optional),
    DuplicateKeyColumns = .pubchem_collapse(.plant_duplicate_key_cols(table_name)),
    Status = ifelse(!present || length(missing_required) > 0, "fail", "pass"),
    stringsAsFactors = FALSE
  )
}

.plant_validate_table = function(table_name, table, dictionary, strict) {
  rows = list()
  if (table_name == "Validation" && is.list(table)) {
    missing_required = setdiff(dictionary$Column[dictionary$Required],
                               names(table))
    for (col in missing_required) {
      rows[[length(rows) + 1]] =
        .plant_validation_issue("error", table_name, col,
                                "Required validation component is missing",
                                "present", "absent", NA_integer_,
                                NA_character_)
    }
    return(.plant_bind_tables(rows, .plant_validation_issue_cols()))
  }
  if (!is.data.frame(table)) {
    rows[[length(rows) + 1]] =
      .plant_validation_issue("error", table_name, NA_character_,
                              "Expected table is missing", "data.frame",
                              class(table)[[1]], NA_integer_, NA_character_)
    return(.plant_bind_tables(rows, .plant_validation_issue_cols()))
  }
  required = dictionary$Column[dictionary$Required]
  missing_required = setdiff(required, names(table))
  for (col in missing_required) {
    rows[[length(rows) + 1]] =
      .plant_validation_issue("error", table_name, col,
                              "Required column is missing", "present",
                              "absent", nrow(table), NA_character_)
  }
  if (isTRUE(strict)) {
    missing_optional = setdiff(dictionary$Column[!dictionary$Required],
                               names(table))
    for (col in missing_optional) {
      rows[[length(rows) + 1]] =
        .plant_validation_issue("warning", table_name, col,
                                "Optional documented column is missing",
                                "present", "absent", nrow(table),
                                NA_character_)
    }
  }
  present_dictionary = dictionary[dictionary$Column %in% names(table), ,
                                  drop = FALSE]
  for (i in seq_len(nrow(present_dictionary))) {
    allowed = .uaf_non_empty(strsplit(
      .uaf_first_non_empty_text(present_dictionary$AllowedValues[[i]]),
      ";", fixed = TRUE
    )[[1]])
    if (length(allowed) < 1) next
    col = present_dictionary$Column[[i]]
    bad = setdiff(unique(.uaf_non_empty(table[[col]])), allowed)
    if (length(bad) > 0) {
      rows[[length(rows) + 1]] =
        .plant_validation_issue("warning", table_name, col,
                                "Column contains values outside documented set",
                                .pubchem_collapse(allowed),
                                .pubchem_collapse(utils::head(bad, 8)),
                                nrow(table), NA_character_)
    }
  }
  dup_cols = .plant_duplicate_key_cols(table_name)
  if (length(dup_cols) > 0 && all(dup_cols %in% names(table)) &&
      nrow(table) > 0) {
    key = do.call(paste, c(table[dup_cols], sep = "\r"))
    key = key[!is.na(key) & key != ""]
    dup = unique(key[duplicated(key)])
    if (length(dup) > 0) {
      rows[[length(rows) + 1]] =
        .plant_validation_issue("warning", table_name,
                                .pubchem_collapse(dup_cols),
                                "Duplicate evidence keys were found",
                                "unique rows",
                                paste(length(dup), "duplicate key(s)"),
                                nrow(table),
                                .pubchem_collapse(utils::head(dup, 5)))
    }
  }
  .plant_bind_tables(rows, .plant_validation_issue_cols())
}

.plant_extra_validation_issues = function(x, issues) {
  rows = list(issues)
  if (is.data.frame(x$PlantNameResolution) &&
      any(x$PlantNameResolution$query_status == "name_needs_review")) {
    rows[[length(rows) + 1]] =
      .plant_validation_issue("warning", "PlantNameResolution",
                              "query_status",
                              "Some plant names need taxonomic review",
                              "parsed_species",
                              "name_needs_review",
                              nrow(x$PlantNameResolution),
                              .pubchem_collapse(x$PlantNameResolution$query_plant[
                                x$PlantNameResolution$query_status ==
                                  "name_needs_review"
                              ]))
  }
  if (is.data.frame(x$CompoundResolution) &&
      any(!(x$CompoundResolution$resolved %in% TRUE))) {
    rows[[length(rows) + 1]] =
      .plant_validation_issue("warning", "CompoundResolution", "resolved",
                              "Some compounds are unresolved",
                              "all TRUE",
                              paste(sum(!(x$CompoundResolution$resolved %in% TRUE)),
                                    "unresolved"),
                              nrow(x$CompoundResolution),
                              .pubchem_collapse(utils::head(
                                x$CompoundResolution$compound_name[
                                  !(x$CompoundResolution$resolved %in% TRUE)
                                ], 8
                              )))
  }
  if (is.data.frame(x$PlantCompoundOccurrences) &&
      "evidence_quality_score" %in% names(x$PlantCompoundOccurrences)) {
    score = suppressWarnings(as.numeric(
      x$PlantCompoundOccurrences$evidence_quality_score
    ))
    bad_score = is.na(score) | score < 0 | score > 1
    if (any(bad_score)) {
      rows[[length(rows) + 1]] =
        .plant_validation_issue("warning", "PlantCompoundOccurrences",
                                "evidence_quality_score",
                                "Evidence quality scores must be between 0 and 1",
                                "0 <= score <= 1",
                                .pubchem_collapse(utils::head(score[bad_score], 8)),
                                nrow(x$PlantCompoundOccurrences),
                                NA_character_)
    }
    bad_ready = x$PlantCompoundOccurrences$analysis_ready == "Yes" &
      !x$PlantCompoundOccurrences$occurrence_status %in%
      c("direct_reported", "curated_reported")
    if (any(bad_ready, na.rm = TRUE)) {
      rows[[length(rows) + 1]] =
        .plant_validation_issue("warning", "PlantCompoundOccurrences",
                                "analysis_ready",
                                "Analysis-ready rows should be direct or curated reported occurrences",
                                "direct_reported or curated_reported",
                                .pubchem_collapse(unique(
                                  x$PlantCompoundOccurrences$occurrence_status[bad_ready]
                                )),
                                nrow(x$PlantCompoundOccurrences),
                                .pubchem_collapse(utils::head(
                                  x$PlantCompoundOccurrences$compound_name[bad_ready],
                                  8
                                )))
    }
  }
  if (is.data.frame(x$ProviderDiagnostics) &&
      any(x$ProviderDiagnostics$status %in% c("not_implemented",
                                              "not_queried",
                                              "warning"))) {
    rows[[length(rows) + 1]] =
      .plant_validation_issue("warning", "ProviderDiagnostics", "status",
                              "One or more providers did not return usable records",
                              "ok or no_records",
                              .pubchem_collapse(unique(x$ProviderDiagnostics$status)),
                              nrow(x$ProviderDiagnostics),
                              .pubchem_collapse(x$ProviderDiagnostics$provider[
                                x$ProviderDiagnostics$status %in%
                                  c("not_implemented", "not_queried",
                                    "warning")
                              ]))
  }
  if (is.data.frame(x$BatchChunkManifest) &&
      "status" %in% names(x$BatchChunkManifest)) {
    status = tolower(as.character(x$BatchChunkManifest$status))
    incomplete = !status %in% c("completed", "complete", "success",
                                "succeeded", "ok", "pass")
    if (any(incomplete)) {
      failed = status %in% c("failed", "error", "timeout", "timed_out",
                             "rate_limited", "stopped")
      labels = if ("query_label" %in% names(x$BatchChunkManifest)) {
        x$BatchChunkManifest$query_label
      } else if ("species" %in% names(x$BatchChunkManifest)) {
        x$BatchChunkManifest$species
      } else {
        rep(NA_character_, nrow(x$BatchChunkManifest))
      }
      rows[[length(rows) + 1]] =
        .plant_validation_issue(
          ifelse(any(failed), "error", "warning"),
          "BatchChunkManifest", "status",
          paste("Plant discovery is incomplete; unprocessed or failed chunks",
                "must be resumed before compound enrichment or analysis."),
          "all completed",
          .pubchem_collapse(unique(status[incomplete])),
          nrow(x$BatchChunkManifest),
          .pubchem_collapse(utils::head(labels[incomplete], 5))
        )
    }
  }
  .plant_bind_tables(rows, .plant_validation_issue_cols())
}

.plant_validation_issue = function(severity, table, column, issue, expected,
                                   observed, row_count, examples) {
  data.frame(
    severity = severity,
    table = table,
    column = column,
    issue = issue,
    expected = expected,
    observed = observed,
    row_count = row_count,
    examples = examples,
    stringsAsFactors = FALSE
  )
}

.plant_duplicate_key_cols = function(table_name) {
  switch(table_name,
         PlantCompoundOccurrences = c("species", "compound_name_clean",
                                      "source_database", "source_record_id",
                                      "pmid", "evidence_url"),
         PlantContextEvidence = c("species", "compound_name_clean",
                                  "context_type", "normalized_context",
                                  "source_database", "source_record_id",
                                  "pmid", "evidence_url", "source_field"),
         LiteratureCandidates = c("species", "source_database", "pmid",
                                  "chemical_mention"),
         ChemistryComparability = c("species", "compound_name_clean",
                                    "source_database", "source_record_id",
                                    "pmid", "evidence_url"),
         character())
}

.plant_lotus_index_cols = function() {
  c("species", "species_clean", "genus", "family", "compound_name",
    "compound_name_clean", "lotus_id", "wikidata_id", "cid", "smiles",
    "inchikey", "molecular_formula", "source_database", "source_record_id",
    "evidence_text", "evidence_url", "reference_id", "pmid", "doi",
    "plant_part", "tissue", "method", "chemical_class_pathway",
    "chemical_class_superclass", "chemical_class_class",
    "all_chem_classifications", "all_taxa", "retrieved_at", "source_file")
}

.plant_empty_lotus_index = function() {
  .uaf_empty_table(.plant_lotus_index_cols())
}

.plant_lotus_build_manifest_cols = function() {
  c("source_id", "source_label", "source_type", "source_path",
    "raw_row_count", "index_row_count", "status", "error_message",
    "elapsed_seconds", "processed_at")
}

.plant_empty_occurrences = function() .uaf_empty_table(.plant_occurrence_cols())

.plant_provenance = function(step, source, query, source_url, record_count,
                             notes, cache_file = NA_character_) {
  data.frame(
    step = step,
    source = source,
    query = query,
    source_url = source_url,
    cache_file = cache_file,
    retrieved_at = .plant_timestamp(),
    record_count = suppressWarnings(as.integer(record_count)),
    notes = notes,
    stringsAsFactors = FALSE
  )
}

.plant_pubmed_query = function(species) {
  species = unique(.uaf_non_empty(species))
  if (length(species) < 1) return(NA_character_)
  taxon_query = paste0('"', species, '"[Title/Abstract]', collapse = " OR ")
  paste0('(', taxon_query, ') AND (phytochemical* OR metabolite* OR ',
         '"secondary metabolite*" OR "chemical profile" OR "metabolic profile" OR ',
         '"natural product*" OR GC-MS OR "GC MS" OR LC-MS OR "LC MS" OR HPLC OR ',
         '"essential oil" OR volatile* OR flavonoid* OR alkaloid* OR terpene* OR ',
         'phenolic* OR "root exudate*")')
}

.plant_ncbi_url = function(endpoint, params, api_key = "") {
  params = params[!is.na(params) & params != ""]
  if (length(.uaf_non_empty(api_key)) > 0) {
    params = c(params, api_key = api_key)
  }
  paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", endpoint, "?",
         paste(names(params), utils::URLencode(params, reserved = TRUE),
               sep = "=", collapse = "&"))
}

.plant_fetch_json = function(url, cache, cache_dir, throttle, request_fun,
                             timeout = 30) {
  cache_file = file.path(cache_dir, paste0(.pubchem_url_hash(url), ".json"))
  if (isTRUE(cache) && file.exists(cache_file)) {
    txt = paste(readLines(cache_file, warn = FALSE, encoding = "UTF-8"),
                collapse = "\n")
    cached = tryCatch(
      jsonlite::fromJSON(txt, simplifyVector = FALSE),
      error = function(error) NULL
    )
    if (!is.null(cached) && is.na(.plant_json_payload_error(cached))) {
      return(.plant_mark_cache_hit(cached, TRUE))
    }
  }
  result = if (is.null(request_fun)) {
    .plant_live_json_request(url, timeout, throttle)
  } else {
    request_fun(url)
  }
  if (is.character(result)) {
    txt = paste(result, collapse = "\n")
    parsed = jsonlite::fromJSON(txt, simplifyVector = FALSE)
    payload_error = .plant_json_payload_error(parsed)
    if (!is.na(payload_error)) {
      stop("Remote JSON error for ", .plant_redact_secrets(url), ": ",
           payload_error, call. = FALSE)
    }
    if (isTRUE(cache)) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      .pubchem_atomic_write_text(txt, cache_file)
    }
    if (!is.na(throttle) && throttle > 0) Sys.sleep(throttle)
    .plant_mark_cache_hit(parsed, FALSE)
  } else {
    payload_error = .plant_json_payload_error(result)
    if (!is.na(payload_error)) {
      stop("Remote JSON error for ", .plant_redact_secrets(url), ": ",
           payload_error, call. = FALSE)
    }
    .plant_mark_cache_hit(result, FALSE)
  }
}

.plant_fetch_text = function(url, cache, cache_dir, throttle, request_fun,
                             timeout = 30) {
  cache_file = file.path(cache_dir, paste0(.pubchem_url_hash(url), ".txt"))
  if (isTRUE(cache) && file.exists(cache_file)) {
    cached = paste(readLines(cache_file, warn = FALSE, encoding = "UTF-8"),
                   collapse = "\n")
    if (is.na(.plant_text_payload_error(cached))) {
      return(.plant_mark_cache_hit(cached, TRUE))
    }
  }
  result = if (is.null(request_fun)) {
    .plant_live_text_request(url, timeout, throttle)
  } else {
    request_fun(url)
  }
  txt = if (is.list(result)) {
    jsonlite::toJSON(result, auto_unbox = TRUE)
  } else {
    paste(as.character(result), collapse = "\n")
  }
  cacheable = is.null(request_fun) || is.character(result)
  if (isTRUE(cache) && isTRUE(cacheable)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    .pubchem_atomic_write_text(txt, cache_file)
  }
  if (!is.na(throttle) && throttle > 0) Sys.sleep(throttle)
  .plant_mark_cache_hit(txt, FALSE)
}

.plant_live_json_request = function(url, timeout, throttle, attempts = 4L) {
  .plant_live_request_retry(
    url = url, timeout = timeout, throttle = throttle, attempts = attempts,
    validator = function(txt) {
      parsed = jsonlite::fromJSON(txt, simplifyVector = FALSE)
      payload_error = .plant_json_payload_error(parsed)
      if (!is.na(payload_error)) {
        stop("Remote JSON error: ", payload_error, call. = FALSE)
      }
      TRUE
    }
  )
}

.plant_live_text_request = function(url, timeout, throttle, attempts = 4L) {
  .plant_live_request_retry(
    url = url, timeout = timeout, throttle = throttle, attempts = attempts,
    validator = function(txt) {
      payload_error = .plant_text_payload_error(txt)
      if (!is.na(payload_error)) {
        stop("Remote text error: ", payload_error, call. = FALSE)
      }
      TRUE
    }
  )
}

.plant_live_request_retry = function(url, timeout, throttle, attempts,
                                     validator) {
  attempts = suppressWarnings(as.integer(attempts[[1]]))
  if (!is.finite(attempts) || attempts < 1L) attempts = 1L
  last_error = NULL
  attempted = 0L
  for (attempt in seq_len(attempts)) {
    attempted = attempt
    result = tryCatch({
      txt = .plant_read_url_text(url, timeout)
      validator(txt)
      txt
    }, error = function(error) error)
    if (!inherits(result, "error")) return(result)
    last_error = result
    if (attempt >= attempts ||
        !.plant_service_busy_message(conditionMessage(result))) break
    delay = .plant_live_retry_delay(attempt, throttle)
    if (delay > 0) Sys.sleep(delay)
  }
  stop(
    "Live request failed after ", attempted, " attempt(s) for ",
    .plant_redact_secrets(url), ": ", conditionMessage(last_error),
    call. = FALSE
  )
}

.plant_live_retry_delay = function(attempt, throttle) {
  throttle = suppressWarnings(as.numeric(throttle[[1]]))
  if (!is.finite(throttle) || throttle < 0) throttle = 0
  min(8, max(1, throttle) * (2 ^ max(0, as.integer(attempt) - 1L)))
}

.plant_json_payload_error = function(x) {
  if (!is.list(x)) return(NA_character_)
  candidates = list(
    x$ERROR, x$error, x$errors,
    tryCatch(x$esearchresult$ERROR, error = function(error) NULL),
    tryCatch(x$esearchresult$error, error = function(error) NULL),
    tryCatch(x$header$ERROR, error = function(error) NULL),
    tryCatch(x$header$error, error = function(error) NULL)
  )
  values = .uaf_non_empty(unlist(candidates, use.names = FALSE))
  if (length(values) < 1) NA_character_ else .pubchem_collapse(values)
}

.plant_text_payload_error = function(txt) {
  txt = .uaf_first_non_empty_text(txt)
  if (is.na(txt)) return("empty response")
  patterns = c(
    "NCBI/eutils[0-9]+ - WWW Error [0-9]+ Diagnostic",
    "<h1[^>]*>Server Error</h1>",
    "<ERROR>[^<]+</ERROR>",
    "Search Backend failed:.*Status: (500|502|503|504)"
  )
  if (!any(vapply(patterns, grepl, logical(1), x = txt,
                  ignore.case = TRUE, perl = TRUE))) {
    return(NA_character_)
  }
  status = regmatches(txt, regexpr("(408|425|429|500|502|503|504)", txt,
                                   perl = TRUE))
  if (length(status) < 1 || !nzchar(status)) status = "service unavailable"
  paste("provider returned", status)
}

.plant_read_url_text = function(url, timeout = 30) {
  timeout = suppressWarnings(as.numeric(timeout[[1]]))
  if (!is.finite(timeout) || timeout <= 0) timeout = 30
  if (requireNamespace("curl", quietly = TRUE)) {
    handle = curl::new_handle(
      followlocation = TRUE, failonerror = FALSE, timeout = timeout,
      connecttimeout = min(timeout, 20)
    )
    curl::handle_setheaders(handle, `User-Agent` = "uafR public-data client")
    response = curl::curl_fetch_memory(url, handle = handle)
    status = suppressWarnings(as.integer(response$status_code))
    if (!is.na(status) && status >= 400L) {
      stop(
        "HTTP ", status, " returned for ", .plant_redact_secrets(url),
        call. = FALSE
      )
    }
    return(rawToChar(response$content))
  }
  warning_message = NA_character_
  tryCatch(
    withCallingHandlers(
      .plant_with_timeout(timeout, {
        con = base::url(url, open = "rb")
        on.exit(close(con), add = TRUE)
        paste(readLines(con, warn = FALSE, encoding = "UTF-8"),
              collapse = "\n")
      }),
      warning = function(warning) {
        warning_message <<- conditionMessage(warning)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) {
      message = paste(.uaf_non_empty(c(conditionMessage(error),
                                       warning_message)), collapse = "; ")
      stop(.plant_redact_secrets(message), call. = FALSE)
    }
  )
}

.plant_fetch_lotus_json = function(url, cache, cache_dir, throttle,
                                   request_fun, timeout = 30,
                                   max_records = Inf) {
  max_records = .plant_max_records(max_records)
  if (!is.finite(max_records)) {
    return(.plant_fetch_json(url, cache, cache_dir, throttle, request_fun,
                             timeout = timeout))
  }
  cache_file = file.path(
    cache_dir,
    paste0(.pubchem_url_hash(paste(url, max_records, sep = "\n")),
           "_first", max_records, ".json")
  )
  if (isTRUE(cache) && file.exists(cache_file)) {
    txt = paste(readLines(cache_file, warn = FALSE, encoding = "UTF-8"),
                collapse = "\n")
    return(.plant_mark_cache_hit(
      jsonlite::fromJSON(txt, simplifyVector = FALSE), TRUE
    ))
  }
  result = if (is.null(request_fun)) {
    .plant_read_lotus_json_prefix(url, max_records, timeout)
  } else {
    request_fun(url)
  }
  if (is.character(result)) {
    txt = paste(result, collapse = "\n")
    if (isTRUE(cache)) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      writeLines(txt, cache_file, useBytes = TRUE)
    }
    if (!is.na(throttle) && throttle > 0) Sys.sleep(throttle)
    .plant_mark_cache_hit(jsonlite::fromJSON(txt, simplifyVector = FALSE),
                          FALSE)
  } else {
    .plant_mark_cache_hit(result, FALSE)
  }
}

.plant_mark_cache_hit = function(x, cache_hit) {
  if (is.null(x)) return(x)
  attr(x, "uaf_cache_hit") = isTRUE(cache_hit)
  x
}

.plant_cache_hit = function(x) {
  isTRUE(attr(x, "uaf_cache_hit", exact = TRUE))
}

.plant_read_lotus_json_prefix = function(url, max_records, timeout = 30,
                                         chunk_size = 65536) {
  max_records = .plant_max_records(max_records)
  if (!is.finite(max_records)) {
    stop("`max_records` must be finite for LOTUS prefix reads.",
         call. = FALSE)
  }
  .plant_with_timeout(timeout, {
    con = .plant_open_lotus_connection(url, timeout)
    on.exit(close(con), add = TRUE)
    pieces = character()
    repeat {
      chunk = readChar(con, nchars = chunk_size, useBytes = TRUE)
      if (!nzchar(chunk)) break
      pieces[[length(pieces) + 1]] = chunk
      txt = paste(pieces, collapse = "")
      prefix = .plant_lotus_json_prefix(txt, max_records)
      if (!is.null(prefix)) return(prefix)
    }
    txt = paste(pieces, collapse = "")
    if (!nzchar(txt)) {
      stop("LOTUS returned an empty response.", call. = FALSE)
    }
    prefix = .plant_lotus_json_prefix(txt, max_records)
    if (!is.null(prefix)) return(prefix)
    txt
  })
}

.plant_open_lotus_connection = function(url, timeout = 30) {
  timeout = suppressWarnings(as.numeric(timeout[[1]]))
  if (!is.finite(timeout) || timeout <= 0) timeout = 30
  if (requireNamespace("curl", quietly = TRUE)) {
    handle = curl::new_handle(timeout = timeout,
                              connecttimeout = min(timeout, 10))
    return(curl::curl(url, open = "rb", handle = handle))
  }
  base::url(url, open = "rb")
}

.plant_lotus_json_prefix = function(txt, max_records) {
  if (!is.character(txt) || length(txt) < 1 || !nzchar(txt[[1]])) {
    return(NULL)
  }
  max_records = .plant_max_records(max_records)
  key = regexpr('"naturalProducts"\\s*:', txt, perl = TRUE)
  if (key[[1]] < 0) return(NULL)
  after_key = key[[1]] + attr(key, "match.length")
  tail = substring(txt, after_key)
  array_start_rel = regexpr("[", tail, fixed = TRUE)
  if (array_start_rel[[1]] < 0) return(NULL)
  array_start = after_key + array_start_rel[[1]] - 1L
  chars = strsplit(substring(txt, array_start + 1L), "", fixed = TRUE)[[1]]
  if (length(chars) < 1) return(NULL)

  in_string = FALSE
  escaped = FALSE
  depth = 0L
  record_count = 0L
  last_record_end = NA_integer_
  for (i in seq_along(chars)) {
    ch = chars[[i]]
    if (in_string) {
      if (escaped) {
        escaped = FALSE
      } else if (identical(ch, "\\")) {
        escaped = TRUE
      } else if (identical(ch, "\"")) {
        in_string = FALSE
      }
      next
    }
    if (identical(ch, "\"")) {
      in_string = TRUE
      next
    }
    if (identical(ch, "{") || identical(ch, "[")) {
      depth = depth + 1L
      next
    }
    if (identical(ch, "}") || identical(ch, "]")) {
      if (depth == 0L && identical(ch, "]")) {
        return(paste0(substring(txt, 1L, array_start + i), "}"))
      }
      depth = max(0L, depth - 1L)
      if (identical(ch, "}") && depth == 0L) {
        record_count = record_count + 1L
        last_record_end = array_start + i
        if (record_count >= max_records) {
          return(paste0(substring(txt, 1L, last_record_end), "]}"))
        }
      }
    }
  }
  NULL
}

.plant_read_lotus_index_source = function(path) {
  if (!file.exists(path)) {
    stop("LOTUS index file does not exist: ", path, call. = FALSE)
  }
  ext = tolower(tools::file_ext(path))
  if (ext %in% c("zip", "bson")) {
    stop("Raw LOTUS MongoDB downloads are zipped BSON dumps and are not ",
         "parsed by standardizeLotusIndex(). Use ",
         "`tools/flatten_lotus_mongo_dump.py` to create a flat CSV, compact ",
         "CSV, and optional lookup directory, then pass those outputs to ",
         "uafR.", call. = FALSE)
  }
  if (ext == "rds") return(readRDS(path))
  if (ext %in% c("csv")) {
    return(utils::read.csv(path, stringsAsFactors = FALSE,
                           check.names = FALSE))
  }
  if (ext %in% c("tsv", "tab", "txt")) {
    return(utils::read.delim(path, stringsAsFactors = FALSE,
                             check.names = FALSE))
  }
  if (ext %in% c("jsonl", "ndjson")) {
    con = file(path, open = "r", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    return(jsonlite::stream_in(con, verbose = FALSE))
  }
  if (ext == "json") {
    out = jsonlite::fromJSON(path, flatten = TRUE)
    if (is.data.frame(out)) return(out)
    if (is.list(out)) {
      frames = out[vapply(out, is.data.frame, logical(1))]
      if (length(frames) > 0) return(frames[[1]])
    }
  }
  stop("Unsupported LOTUS index file extension: ", ext,
       ". Use CSV, TSV, JSON, JSONL, NDJSON, or RDS.", call. = FALSE)
}

.plant_lotus_build_sources = function(input, recursive = TRUE) {
  if (is.data.frame(input)) {
    return(list(list(type = "data_frame",
                     path = NA_character_,
                     label = "data_frame",
                     data = input)))
  }
  if (!is.character(input) || length(input) < 1) {
    stop("`input` must be a data frame, file path, directory path, or ",
         "character vector of file/directory paths.", call. = FALSE)
  }
  paths = .uaf_non_empty(input)
  out = list()
  for (path in paths) {
    if (dir.exists(path)) {
      files = .plant_lotus_supported_files(path, recursive = recursive)
      for (file in files) {
        out[[length(out) + 1]] = list(
          type = "file",
          path = file,
          label = normalizePath(file, winslash = "/", mustWork = FALSE),
          data = NULL
        )
      }
    } else {
      out[[length(out) + 1]] = list(
        type = "file",
        path = path,
        label = normalizePath(path, winslash = "/", mustWork = FALSE),
        data = NULL
      )
    }
  }
  out
}

.plant_lotus_supported_files = function(path, recursive = TRUE) {
  files = list.files(path, full.names = TRUE, recursive = recursive,
                     no.. = TRUE)
  files[file.exists(files) & !dir.exists(files) &
          tolower(tools::file_ext(files)) %in%
          c(.plant_lotus_supported_flat_ext(), "zip", "bson")]
}

.plant_lotus_supported_flat_ext = function() {
  c("csv", "tsv", "tab", "txt", "json", "jsonl", "ndjson", "rds")
}

.plant_lotus_build_read_source = function(source) {
  if (identical(source$type, "data_frame")) return(source$data)
  .plant_read_lotus_index_source(source$path)
}

.plant_deduplicate_lotus_index = function(index) {
  index = .plant_bind_tables(list(index), .plant_lotus_index_cols())
  if (nrow(index) < 2) return(index)
  key_cols = c("species_clean", "compound_name_clean", "source_record_id",
               "lotus_id", "wikidata_id", "pmid", "doi")
  key = do.call(paste, c(index[key_cols], sep = "\r"))
  out = index[!duplicated(key), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_lotus_index_output_format = function(out_file, format) {
  format = .uaf_first_non_empty_text(format, "auto")
  if (format != "auto") return(format)
  ext = tolower(tools::file_ext(.uaf_first_non_empty_text(out_file, "")))
  if (ext == "rds") return("rds")
  "csv"
}

.plant_write_lotus_index_file = function(index, out_file, format, overwrite) {
  out_file = .uaf_first_non_empty_text(out_file)
  if (is.na(out_file)) return(NA_character_)
  if (file.exists(out_file) && !isTRUE(overwrite)) {
    stop("Output file exists and `overwrite = FALSE`: ", out_file,
         call. = FALSE)
  }
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (format == "rds") {
    saveRDS(index, out_file)
  } else if (format == "csv") {
    utils::write.csv(index, out_file, row.names = FALSE, na = "")
  } else {
    stop("Unsupported LOTUS index output format: ", format, call. = FALSE)
  }
  normalizePath(out_file, winslash = "/", mustWork = FALSE)
}

.plant_lotus_manifest_file = function(manifest_file, out_file) {
  manifest_file = .uaf_first_non_empty_text(manifest_file)
  if (!is.na(manifest_file)) return(manifest_file)
  out_file = .uaf_first_non_empty_text(out_file)
  if (is.na(out_file)) return(NA_character_)
  paste0(tools::file_path_sans_ext(out_file), "_manifest.json")
}

.plant_lotus_index_build_summary = function(index, manifest_rows,
                                            pre_dedup_rows, out_file,
                                            manifest_file) {
  manifest = .plant_bind_tables(manifest_rows,
                                .plant_lotus_build_manifest_cols())
  source_count = nrow(manifest)
  failed_count = sum(manifest$status == "error", na.rm = TRUE)
  no_taxon_count = sum(manifest$status == "no_taxon_records", na.rm = TRUE)
  context_known = !.plant_context_group(index$plant_part, NA_character_,
                                        "plant_part") %in%
    c("unknown", "extract_unspecified") |
    !.plant_context_group(index$tissue, NA_character_, "tissue") %in%
    c("unknown", "extract_unspecified") |
    !.plant_method_group(index$method, NA_character_, "LOTUS",
                         "direct_species_database") %in%
    c("unknown", "database_record")
  reference_known = !is.na(.uaf_squish_text(index$doi)) |
    !is.na(.uaf_squish_text(index$pmid)) |
    !is.na(.uaf_squish_text(index$reference_id))
  data.frame(
    built_at = .plant_timestamp(),
    source_count = source_count,
    failed_source_count = failed_count,
    no_taxon_source_count = no_taxon_count,
    raw_row_count = sum(suppressWarnings(as.integer(manifest$raw_row_count)),
                        na.rm = TRUE),
    pre_dedup_index_row_count = pre_dedup_rows,
    index_row_count = nrow(index),
    duplicate_row_count = max(0L, pre_dedup_rows - nrow(index)),
    species_count = length(unique(.uaf_non_empty(index$species_clean))),
    genus_count = length(unique(.uaf_non_empty(.plant_clean_name(index$genus)))),
    family_count = length(unique(.uaf_non_empty(.plant_clean_name(index$family)))),
    compound_count = length(unique(.uaf_non_empty(index$compound_name_clean))),
    rows_with_lotus_id = sum(!is.na(.uaf_squish_text(index$lotus_id))),
    rows_with_wikidata_id = sum(!is.na(.uaf_squish_text(index$wikidata_id))),
    rows_with_inchikey = sum(!is.na(.uaf_squish_text(index$inchikey))),
    rows_with_smiles = sum(!is.na(.uaf_squish_text(index$smiles))),
    rows_with_reference = sum(reference_known, na.rm = TRUE),
    rows_with_biological_context = sum(context_known, na.rm = TRUE),
    out_file = .uaf_first_non_empty_text(out_file),
    manifest_file = .uaf_first_non_empty_text(manifest_file),
    stringsAsFactors = FALSE
  )
}

.plant_standardize_lotus_index_row = function(row, source_file = NULL) {
  compound_name = .plant_lotus_index_field(
    row,
    c("compound_name", "natural_product_name", "traditional_name",
      "molecule_name", "name", "title", "iupac_name")
  )
  if (is.na(compound_name) || compound_name == "") {
    return(.plant_empty_lotus_index())
  }

  taxon_text = .plant_lotus_index_taxon_text(row)
  explicit_species = .plant_lotus_index_field(
    row,
    c("species", "organism_species", "taxon_species", "matched_species")
  )
  explicit_species = .plant_canonical_taxon_name(explicit_species)
  species_terms = unique(.uaf_non_empty(c(
    if (!is.na(explicit_species) && .plant_is_binomial(explicit_species)) {
      explicit_species
    } else {
      NA_character_
    },
    .lotus_extract_binomial(taxon_text)
  )))

  lineage = .lotus_lineage_ranks(taxon_text,
                                 .uaf_first_non_empty_text(species_terms))
  genus_value = .uaf_first_non_empty_text(
    .plant_lotus_index_field(row, c("genus", "organism_genus",
                                    "taxon_genus")),
    lineage$Genus
  )
  family_value = .uaf_first_non_empty_text(
    .plant_lotus_index_field(row, c("family", "organism_family",
                                    "taxon_family")),
    lineage$Family
  )
  if (length(species_terms) < 1) {
    species_terms = NA_character_
  }

  base = list(
    compound_name = compound_name,
    lotus_id = .plant_lotus_index_field(row, c("lotus_id", "lotus",
                                               "lotus_identifier")),
    wikidata_id = .plant_lotus_index_field(row, c("wikidata_id", "wikidata",
                                                  "wikidata_id_id",
                                                  "all_wikidata_ids")),
    cid = .plant_lotus_index_field(row, c("cid", "pubchem_cid",
                                          "pubchem")),
    smiles = .plant_lotus_index_field(row, c("smiles", "canonical_smiles",
                                             "isomeric_smiles", "smiles2d")),
    inchikey = .plant_lotus_index_field(row, c("inchikey", "inchikey2d",
                                               "inchi_key")),
    molecular_formula = .plant_lotus_index_field(row, c("molecular_formula",
                                                        "formula")),
    evidence_text = .plant_lotus_index_field(row, c("evidence_text",
                                                    "evidence_note",
                                                    "description")),
    evidence_url = .plant_lotus_index_field(row, c("evidence_url",
                                                   "source_url", "url")),
    reference_id = .plant_lotus_index_field(row, c("reference_id",
                                                   "reference",
                                                   "citation")),
    pmid = .plant_lotus_index_field(row, c("pmid", "pubmed_id")),
    doi = .plant_lotus_index_field(row, c("doi")),
    plant_part = .plant_lotus_index_field(row, c("plant_part",
                                                 "part", "organ")),
    tissue = .plant_lotus_index_field(row, c("tissue")),
    method = .plant_lotus_index_field(row, c("method",
                                             "analytical_method",
                                             "instrument")),
    chemical_class_pathway = .plant_lotus_index_field(
      row,
      c("chemical_taxonomy_npclassifier_pathway", "npclassifier_pathway",
        "pathway")
    ),
    chemical_class_superclass = .plant_lotus_index_field(
      row,
      c("chemical_taxonomy_npclassifier_superclass", "npclassifier_superclass",
        "superclass")
    ),
    chemical_class_class = .plant_lotus_index_field(
      row,
      c("chemical_taxonomy_npclassifier_class", "npclassifier_class",
        "chemical_class", "class")
    ),
    all_chem_classifications = .plant_lotus_index_fields_text(
      row,
      "chem.*class|npclassifier|chemical_taxonomy|all_chem"
    ),
    all_taxa = taxon_text,
    source_file = .uaf_first_non_empty_text(source_file)
  )

  rows = lapply(species_terms, function(species) {
    species = .plant_canonical_taxon_name(species)
    genus = .uaf_first_non_empty_text(genus_value, .plant_genus(species))
    data.frame(
      species = species,
      species_clean = .plant_clean_name(species),
      genus = genus,
      family = family_value,
      compound_name = base$compound_name,
      compound_name_clean = .plant_clean_compound(base$compound_name),
      lotus_id = base$lotus_id,
      wikidata_id = base$wikidata_id,
      cid = base$cid,
      smiles = base$smiles,
      inchikey = base$inchikey,
      molecular_formula = base$molecular_formula,
      source_database = "LOTUS",
      source_record_id = .uaf_first_non_empty_text(base$lotus_id,
                                                   base$wikidata_id,
                                                   base$cid),
      evidence_text = .uaf_first_non_empty_text(
        base$evidence_text,
        .plant_lotus_index_evidence_text(base)
      ),
      evidence_url = .uaf_first_non_empty_text(
        base$evidence_url,
        .plant_lotus_index_evidence_url(base)
      ),
      reference_id = .uaf_first_non_empty_text(base$reference_id,
                                               base$lotus_id,
                                               base$wikidata_id),
      pmid = .plant_first_pattern(base$pmid, "\\b\\d{7,9}\\b"),
      doi = .plant_first_pattern(base$doi,
                                 "10\\.\\d{4,9}/[-._;()/:A-Za-z0-9]+"),
      plant_part = base$plant_part,
      tissue = base$tissue,
      method = base$method,
      chemical_class_pathway = base$chemical_class_pathway,
      chemical_class_superclass = base$chemical_class_superclass,
      chemical_class_class = base$chemical_class_class,
      all_chem_classifications = base$all_chem_classifications,
      all_taxa = base$all_taxa,
      retrieved_at = .plant_timestamp(),
      source_file = base$source_file,
      stringsAsFactors = FALSE
    )
  })
  .plant_bind_tables(rows, .plant_lotus_index_cols())
}

.plant_lotus_index_field = function(row, candidates) {
  if (!is.data.frame(row) || nrow(row) < 1) return(NA_character_)
  names_row = names(row)
  candidates = .plant_normalize_column_names(candidates)
  for (candidate in candidates) {
    hit = which(names_row == candidate)
    if (length(hit) > 0) {
      value = .uaf_first_non_empty_text(unlist(row[hit], use.names = FALSE))
      if (!is.na(value)) return(value)
    }
  }
  for (candidate in candidates) {
    hit = grep(candidate, names_row, ignore.case = TRUE)
    if (length(hit) > 0) {
      value = .uaf_first_non_empty_text(unlist(row[hit], use.names = FALSE))
      if (!is.na(value)) return(value)
    }
  }
  NA_character_
}

.plant_lotus_index_fields_text = function(row, pattern) {
  if (!is.data.frame(row) || nrow(row) < 1) return(NA_character_)
  hit = grep(pattern, names(row), ignore.case = TRUE, perl = TRUE)
  if (length(hit) < 1) return(NA_character_)
  .pubchem_collapse(unlist(row[hit], use.names = FALSE))
}

.plant_lotus_index_taxon_text = function(row) {
  if (!is.data.frame(row) || nrow(row) < 1) return(NA_character_)
  hit = grep("all_?taxa|taxonomy|taxon|organism|species|genus|family|lineage",
             names(row), ignore.case = TRUE, perl = TRUE)
  hit = hit[!grepl("chem|npclassifier|class", names(row)[hit],
                   ignore.case = TRUE, perl = TRUE)]
  if (length(hit) < 1) return(NA_character_)
  explicit = .pubchem_collapse(unlist(row[hit], use.names = FALSE))
  .uaf_first_non_empty_text(explicit)
}

.plant_lotus_index_compound_id = function(hit) {
  .uaf_first_non_empty_text(hit$lotus_id, hit$inchikey, hit$wikidata_id,
                            ifelse(!is.na(hit$cid) & hit$cid != "",
                                   paste0("cid_", hit$cid), NA_character_),
                            hit$compound_name_clean)
}

.plant_lotus_index_compound_id_type = function(hit) {
  if (!is.na(.uaf_first_non_empty_text(hit$lotus_id))) return("LOTUS")
  if (!is.na(.uaf_first_non_empty_text(hit$inchikey))) return("InChIKey")
  if (!is.na(.uaf_first_non_empty_text(hit$wikidata_id))) return("Wikidata")
  if (!is.na(.uaf_first_non_empty_text(hit$cid))) return("PubChem_CID")
  "compound_name"
}

.plant_lotus_index_source_record_id = function(hit) {
  .uaf_first_non_empty_text(hit$source_record_id, hit$lotus_id,
                            hit$wikidata_id, hit$cid)
}

.plant_lotus_index_evidence_url = function(hit) {
  url = .uaf_first_non_empty_text(hit$evidence_url)
  if (!is.na(url)) return(url)
  wikidata = .uaf_first_non_empty_text(hit$wikidata_id)
  if (!is.na(wikidata) && grepl("^Q\\d+$", wikidata)) {
    return(paste0("https://www.wikidata.org/wiki/", wikidata))
  }
  "https://lotus.naturalproducts.net/download"
}

.plant_lotus_index_evidence_text = function(hit) {
  pieces = c(
    paste0("LOTUS local index record for ",
           .uaf_first_non_empty_text(hit$compound_name, "unknown compound")),
    if (!is.na(.uaf_first_non_empty_text(hit$lotus_id))) {
      paste0("LOTUS ID: ", hit$lotus_id)
    },
    if (!is.na(.uaf_first_non_empty_text(hit$species))) {
      paste0("taxon: ", hit$species)
    },
    if (!is.na(.uaf_first_non_empty_text(hit$chemical_class_pathway))) {
      paste0("pathway: ", hit$chemical_class_pathway)
    },
    if (!is.na(.uaf_first_non_empty_text(hit$chemical_class_superclass))) {
      paste0("superclass: ", hit$chemical_class_superclass)
    },
    if (!is.na(.uaf_first_non_empty_text(hit$chemical_class_class))) {
      paste0("class: ", hit$chemical_class_class)
    },
    if (!is.na(.uaf_first_non_empty_text(hit$doi))) {
      paste0("DOI: ", hit$doi)
    },
    if (!is.na(.uaf_first_non_empty_text(hit$pmid))) {
      paste0("PMID: ", hit$pmid)
    }
  )
  .plant_truncate(paste(.uaf_non_empty(pieces), collapse = "; "), 800)
}

.plant_with_timeout = function(timeout, expr) {
  timeout = suppressWarnings(as.numeric(timeout[[1]]))
  if (!is.finite(timeout) || timeout <= 0) return(force(expr))
  old_timeout = getOption("timeout")
  options(timeout = max(1, timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  force(expr)
}

.plant_pubmed_ids = function(search) {
  ids = tryCatch(search$esearchresult$idlist, error = function(error) NULL)
  .uaf_non_empty(unlist(ids, use.names = FALSE))
}

.plant_pubmed_total_hits = function(search, fallback = 0L) {
  value = tryCatch(search$esearchresult$count, error = function(error) NULL)
  value = suppressWarnings(as.integer(.uaf_first_non_empty_text(value)))
  if (!is.finite(value)) suppressWarnings(as.integer(fallback)) else value
}

.plant_pubmed_summary_rows = function(query_row, ids, summary,
                                      abstracts = character(),
                                      total_hits = length(ids),
                                      truncated = FALSE,
                                      query = NA_character_) {
  rows = list()
  for (id in ids) {
    item = tryCatch(summary$result[[id]], error = function(error) NULL)
    title = .uaf_first_non_empty_text(item$title)
    doi = .plant_article_id(item, "doi")
    abstract = if (!is.null(names(abstracts)) && id %in% names(abstracts)) {
      .uaf_first_non_empty_text(abstracts[[id]])
    } else NA_character_
    rows[[length(rows) + 1]] = data.frame(
      query_plant = query_row$query_plant,
      query_plant_clean = query_row$query_plant_clean,
      species = query_row$species,
      genus = query_row$genus,
      family = query_row$family,
      source_database = "PubMed",
      source_record_id = id,
      pmid = id,
      doi = doi,
      title = title,
      abstract = abstract,
      chemical_mention = NA_character_,
      species_mention = query_row$species,
      evidence_text = .uaf_first_non_empty_text(abstract, title),
      evidence_url = paste0("https://pubmed.ncbi.nlm.nih.gov/", id, "/"),
      literature_query = query,
      literature_total_hit_count = suppressWarnings(as.integer(total_hits)),
      literature_search_truncated = .uaf_yes_no(truncated),
      abstract_retrieval_status = ifelse(is.na(abstract), "not_available",
                                         "retrieved"),
      retrieved_at = .plant_timestamp(),
      confidence = "low",
      curation_flag = "literature_candidate",
      evidence_tier = "direct_species_literature",
      stringsAsFactors = FALSE
    )
  }
  .plant_bind_tables(rows, .plant_literature_cols())
}

.plant_pubmed_abstracts = function(xml) {
  xml = paste(.uaf_non_empty(xml), collapse = "\n")
  if (!nzchar(xml) || !grepl("<PubmedArticle", xml, fixed = TRUE)) {
    return(character())
  }
  matches = gregexpr(
    "(?s)<PubmedArticle(?:\\s[^>]*)?>.*?</PubmedArticle>", xml,
                     perl = TRUE)
  articles = regmatches(xml, matches)[[1]]
  if (length(articles) < 1 || identical(articles, "")) return(character())
  out = character()
  for (article in articles) {
    pmid = .plant_xml_first_text(article, "PMID")
    abstract_matches = gregexpr(
      "(?s)<AbstractText(?:\\s[^>]*)?>(.*?)</AbstractText>", article,
      perl = TRUE
    )
    abstract_parts = regmatches(article, abstract_matches)[[1]]
    abstract_parts = vapply(abstract_parts, function(value) {
      value = sub("^<AbstractText(?:\\s[^>]*)?>", "", value, perl = TRUE)
      value = sub("</AbstractText>$", "", value, perl = TRUE)
      .plant_xml_text(value)
    }, character(1))
    abstract = .uaf_first_non_empty_text(
      paste(.uaf_non_empty(abstract_parts), collapse = " ")
    )
    if (!is.na(pmid) && !is.na(abstract)) out[[pmid]] = abstract
  }
  out
}

.plant_xml_first_text = function(xml, tag) {
  pattern = paste0("<", tag, "(?:\\s[^>]*)?>(.*?)</", tag, ">")
  hit = regmatches(xml, regexpr(pattern, xml, perl = TRUE))
  if (length(hit) < 1 || identical(hit, "")) return(NA_character_)
  hit = sub(paste0("^<", tag, "(?:\\s[^>]*)?>"), "", hit, perl = TRUE)
  hit = sub(paste0("</", tag, ">$"), "", hit, perl = TRUE)
  .plant_xml_text(hit)
}

.plant_xml_text = function(x) {
  x = gsub("<[^>]+>", " ", x, perl = TRUE)
  replacements = c("&lt;" = "<", "&gt;" = ">", "&amp;" = "&",
                   "&quot;" = "\"", "&apos;" = "'")
  for (from in names(replacements)) x = gsub(from, replacements[[from]], x,
                                               fixed = TRUE)
  .uaf_squish_text(x)
}

.plant_taxon_query_terms = function(query_row, plant_aliases = NULL) {
  fallbacks = .uaf_non_empty(strsplit(
    .uaf_first_non_empty_text(query_row$taxon_fallback),
    ";",
    fixed = TRUE
  )[[1]])
  fallbacks = tolower(.uaf_squish_text(fallbacks))
  ranks = unique(c("species", fallbacks))
  rows = list()
  if ("species" %in% ranks) {
    species_terms = .plant_verified_species_aliases(plant_aliases, query_row)
    for (term in species_terms) {
      rows[[length(rows) + 1]] = data.frame(
        term = term,
        rank = "species",
        stringsAsFactors = FALSE
      )
    }
  }
  if ("genus" %in% ranks && length(.uaf_non_empty(query_row$genus)) > 0) {
    rows[[length(rows) + 1]] = data.frame(
      term = .uaf_first_non_empty_text(query_row$genus),
      rank = "genus",
      stringsAsFactors = FALSE
    )
  }
  if ("family" %in% ranks && length(.uaf_non_empty(query_row$family)) > 0) {
    rows[[length(rows) + 1]] = data.frame(
      term = .uaf_first_non_empty_text(query_row$family),
      rank = "family",
      stringsAsFactors = FALSE
    )
  }
  out = .plant_bind_tables(rows, c("term", "rank"))
  out = out[!is.na(out$term) & out$term != "", , drop = FALSE]
  unique(out)
}

.plant_lotus_rows = function(query_row, result, url, max_records,
                             accepted_species = NULL) {
  items = .plant_json_records(result)
  if (length(items) < 1) return(.plant_empty_occurrences())
  max_records = .plant_max_records(max_records)
  if (is.finite(max_records)) items = utils::head(items, max_records)
  rows = list()
  for (item in items) {
    flat = .plant_flatten_record(item)
    if (length(flat) < 1) next
    evidence_text = .plant_truncate(.pubchem_collapse(flat), 800)
    taxon_values = flat[grepl("organism|taxon|species|genus|family",
                              names(flat), ignore.case = TRUE)]
    taxon_text = .pubchem_collapse(taxon_values)
    matched_rank = .plant_matched_rank_from_text(
      query_row, taxon_text, accepted_species
    )
    if (is.na(matched_rank)) {
      matched_rank = .plant_matched_rank_from_text(
        query_row, evidence_text, accepted_species
      )
    }
    if (is.na(matched_rank)) next
    compound_name = .plant_lotus_compound_name(flat, query_row)
    if (is.na(compound_name)) next
    record_id = .plant_lotus_record_id(flat)
    rows[[length(rows) + 1]] = data.frame(
      query_plant = query_row$query_plant,
      query_plant_clean = query_row$query_plant_clean,
      matched_taxon = ifelse(
        matched_rank == "genus", query_row$genus,
        ifelse(matched_rank == "family", query_row$family,
               .plant_matched_species_from_text(
                 taxon_text, accepted_species, query_row$species
               ))
      ),
      matched_rank = matched_rank,
      species = query_row$species,
      genus = query_row$genus,
      family = query_row$family,
      compound_name = compound_name,
      compound_name_clean = .plant_clean_compound(compound_name),
      compound_id = record_id,
      compound_id_type = ifelse(is.na(record_id), NA_character_,
                                "LOTUS_or_Wikidata"),
      source_database = "LOTUS",
      source_record_id = record_id,
      evidence_text = evidence_text,
      evidence_url = url,
      reference_id = record_id,
      pmid = .plant_first_pattern(flat, "\\b\\d{7,9}\\b"),
      doi = .plant_first_pattern(flat, "10\\.\\d{4,9}/[-._;()/:A-Za-z0-9]+"),
      plant_part = .plant_lotus_plant_part(flat),
      tissue = .plant_lotus_tissue(flat),
      method = .plant_lotus_method(flat),
      occurrence_type = "database_taxon_record",
      retrieved_at = .plant_timestamp(),
      confidence = ifelse(matched_rank == "species", "high", "medium"),
      curation_flag = "database_record_review_recommended",
      evidence_tier = ifelse(matched_rank == "species",
                             "direct_species_database",
                             paste0(matched_rank, "_database_fallback")),
      stringsAsFactors = FALSE
    )
  }
  .plant_bind_occurrences(rows)
}

.plant_knapsack_rows = function(query_row, html, matched_taxon, matched_rank,
                                url, max_records) {
  if (is.null(html) || length(.uaf_non_empty(html)) < 1) {
    return(.plant_empty_occurrences())
  }
  table_rows = .plant_html_row_records(html)
  if (length(table_rows) < 1) return(.plant_empty_occurrences())
  max_records = .plant_max_records(max_records)
  if (is.finite(max_records)) table_rows = utils::head(table_rows, max_records)
  rows = list()
  for (row in table_rows) {
    cells = row$cells
    if (length(cells) < 2) next
    returned_taxon = .plant_knapsack_returned_taxon(cells)
    if (!.plant_knapsack_taxon_matches(returned_taxon, matched_taxon,
                                       matched_rank)) next
    c_id = .uaf_first_non_empty_text(cells[grepl("^C\\d{6,}$", cells)])
    compound_name = .plant_knapsack_compound_name(cells, query_row)
    if (is.na(compound_name)) next
    detail_url = .plant_knapsack_detail_url(row$hrefs, c_id)
    evidence_text = .plant_knapsack_evidence_text(cells, detail_url)
    rows[[length(rows) + 1]] = data.frame(
      query_plant = query_row$query_plant,
      query_plant_clean = query_row$query_plant_clean,
      matched_taxon = returned_taxon,
      matched_rank = matched_rank,
      species = query_row$species,
      genus = query_row$genus,
      family = query_row$family,
      compound_name = compound_name,
      compound_name_clean = .plant_clean_compound(compound_name),
      compound_id = c_id,
      compound_id_type = ifelse(is.na(c_id), NA_character_, "KNApSAcK_C_ID"),
      source_database = "KNApSAcK",
      source_record_id = c_id,
      evidence_text = evidence_text,
      evidence_url = .uaf_first_non_empty_text(detail_url, url),
      reference_id = c_id,
	      pmid = NA_character_,
	      doi = NA_character_,
	      plant_part = .plant_provider_context_cell_value(cells, "plant_part"),
	      tissue = .plant_provider_context_cell_value(cells, "tissue"),
	      method = .plant_provider_context_cell_value(cells, "method"),
	      occurrence_type = "organism_metabolite_record",
      retrieved_at = .plant_timestamp(),
      confidence = ifelse(matched_rank == "species", "high", "medium"),
      curation_flag = ifelse(matched_rank == "species",
                             "source_database_record",
                             "taxon_fallback_review_required"),
      evidence_tier = ifelse(matched_rank == "species",
                             "direct_species_database",
                             paste0(matched_rank, "_database_fallback")),
      stringsAsFactors = FALSE
    )
  }
  .plant_bind_occurrences(rows)
}

.plant_knapsack_result_url = function(html, query_term) {
  html = paste(.uaf_non_empty(html), collapse = "\n")
  if (!nzchar(html) || grepl("<table", html, ignore.case = TRUE)) {
    return(NA_character_)
  }
  pattern = "result[.]php[?][^'\"]+"
  hit = regmatches(html, regexpr(pattern, html, ignore.case = TRUE,
                                perl = TRUE))
  hit = .uaf_first_non_empty_text(hit)
  if (is.na(hit)) {
    return(paste0(
      "https://www.knapsackfamily.com/knapsack_core/result.php?sname=organism&word=",
      utils::URLencode(query_term, reserved = TRUE)
    ))
  }
  hit = gsub("&#32;", "%20", hit, fixed = TRUE)
  .plant_knapsack_absolute_url(hit)
}

.plant_knapsack_returned_taxon = function(cells) {
  cells = .uaf_non_empty(.uaf_squish_text(cells))
  if (length(cells) < 1) return(NA_character_)
  .uaf_first_non_empty_text(utils::tail(cells, 1L))
}

.plant_knapsack_taxon_matches = function(returned_taxon, requested_taxon,
                                         requested_rank) {
  if (requested_rank == "genus") {
    returned = .plant_clean_name(returned_taxon)
    requested = .plant_clean_name(requested_taxon)
    if (is.na(returned) || is.na(requested) || returned == "" ||
        requested == "") return(FALSE)
    return(identical(strsplit(returned, " ", fixed = TRUE)[[1]][[1]],
                     strsplit(requested, " ", fixed = TRUE)[[1]][[1]]))
  }
  if (requested_rank != "species") return(FALSE)
  returned = .plant_clean_name(.plant_taxon_core_name(returned_taxon))
  requested = .plant_clean_name(.plant_taxon_core_name(requested_taxon))
  if (is.na(returned) || is.na(requested) || returned == "" ||
      requested == "") return(FALSE)
  identical(returned, requested)
}

.plant_ncbi_taxonomy_search = function(species, cache, cache_dir, throttle,
                                       ncbi_email, ncbi_tool, ncbi_api_key,
                                       request_fun, request_timeout) {
  if (is.null(request_fun)) {
    url = paste0(
      "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon_suggest/",
      utils::URLencode(species, reserved = TRUE)
    )
    api_key = .uaf_first_non_empty_text(ncbi_api_key)
    if (!is.na(api_key)) {
      url = paste0(url, "?api_key=", utils::URLencode(api_key,
                                                       reserved = TRUE))
    }
    return(.plant_fetch_json(
      url, cache, file.path(cache_dir, "ncbi_taxonomy_suggest"), throttle,
      request_fun = NULL, timeout = request_timeout
    ))
  }
  term = paste0('"', species, '"[Scientific Name]')
  url = .plant_ncbi_url(
    endpoint = "esearch.fcgi",
    params = c(db = "taxonomy", term = term, retmode = "json",
               retmax = "5", tool = ncbi_tool, email = ncbi_email),
    api_key = ncbi_api_key
  )
  .plant_fetch_json(url, cache, file.path(cache_dir, "ncbi_taxonomy"),
                    throttle, request_fun, timeout = request_timeout)
}

.plant_is_taxonomy_suggest_result = function(x) {
  is.list(x) && !is.null(x$sci_name_and_ids)
}

.plant_verified_taxonomy_suggest_matches = function(search,
                                                     accepted_species) {
  columns = c("taxid", "scientific_name", "rank")
  if (!.plant_is_taxonomy_suggest_result(search)) {
    return(.uaf_empty_table(columns))
  }
  accepted = .plant_clean_name(.plant_taxon_core_name(accepted_species))
  accepted = unique(.uaf_non_empty(accepted))
  if (length(accepted) < 1) return(.uaf_empty_table(columns))
  records = search$sci_name_and_ids
  if (is.data.frame(records)) {
    records = lapply(seq_len(nrow(records)), function(i) {
      as.list(records[i, , drop = FALSE])
    })
  }
  if (!is.list(records)) return(.uaf_empty_table(columns))
  rows = lapply(records, function(record) {
    if (!is.list(record)) return(NULL)
    scientific_name = .uaf_first_non_empty_text(
      record$sci_name, record$scientific_name, record$scientificName
    )
    matched_term = .uaf_first_non_empty_text(record$matched_term)
    taxid = .uaf_first_non_empty_text(record$tax_id, record$taxid)
    rank = tolower(.uaf_first_non_empty_text(record$rank, ""))
    names_to_check = .plant_clean_name(.plant_taxon_core_name(
      c(scientific_name, matched_term)
    ))
    names_to_check = unique(.uaf_non_empty(names_to_check))
    rank_ok = !nzchar(rank) || rank %in%
      c("species", "subspecies", "varietas", "forma", "hybrid")
    if (is.na(taxid) || length(intersect(names_to_check, accepted)) < 1 ||
        !rank_ok) {
      return(NULL)
    }
    data.frame(
      taxid = taxid,
      scientific_name = scientific_name,
      rank = ifelse(nzchar(rank), rank, NA_character_),
      stringsAsFactors = FALSE
    )
  })
  out = .plant_bind_tables(rows, columns)
  if (nrow(out) < 1) return(out)
  out[!duplicated(paste(out$taxid, out$scientific_name, sep = "\r")), ,
      drop = FALSE]
}

.plant_ncbi_taxonomy_summary = function(taxids, cache, cache_dir, throttle,
                                         ncbi_email, ncbi_tool, ncbi_api_key,
                                         request_fun, request_timeout) {
  taxids = unique(.uaf_non_empty(as.character(taxids)))
  if (length(taxids) < 1) return(NULL)
  url = .plant_ncbi_url(
    endpoint = "esummary.fcgi",
    params = c(db = "taxonomy", id = paste(taxids, collapse = ","),
               retmode = "json", tool = ncbi_tool, email = ncbi_email),
    api_key = ncbi_api_key
  )
  .plant_fetch_json(url, cache, file.path(cache_dir, "ncbi_taxonomy"),
                    throttle, request_fun, timeout = request_timeout)
}

.plant_verified_taxonomy_matches = function(summary, taxids,
                                             accepted_species) {
  columns = c("taxid", "scientific_name", "rank")
  if (!is.list(summary) || is.null(summary$result)) {
    return(.uaf_empty_table(columns))
  }
  accepted = .plant_clean_name(.plant_taxon_core_name(accepted_species))
  accepted = unique(.uaf_non_empty(accepted))
  if (length(accepted) < 1) return(.uaf_empty_table(columns))
  result = summary$result
  rows = list()
  for (taxid in unique(.uaf_non_empty(as.character(taxids)))) {
    record = result[[taxid]]
    if (is.null(record)) next
    scientific_name = .uaf_first_non_empty_text(
      record$scientificname, record$scientific_name,
      record$scientificName, record$title
    )
    rank = tolower(.uaf_first_non_empty_text(record$rank, ""))
    scientific_key = .plant_clean_name(
      .plant_taxon_core_name(scientific_name)
    )
    rank_ok = !nzchar(rank) || rank %in%
      c("species", "subspecies", "varietas", "forma", "hybrid")
    if (is.na(scientific_key) || !scientific_key %in% accepted || !rank_ok) {
      next
    }
    rows[[length(rows) + 1L]] = data.frame(
      taxid = taxid,
      scientific_name = scientific_name,
      rank = ifelse(nzchar(rank), rank, NA_character_),
      stringsAsFactors = FALSE
    )
  }
  .plant_bind_tables(rows, columns)
}

.plant_pubchem_taxonomy_rows = function(query_row, taxid, result, url,
                                        max_records) {
  annotations = .pubchem_parse_pugview(
    json = result,
    cid = NA_integer_,
    query = query_row$species,
    heading = paste0("Taxonomy: ", query_row$species),
    pubchem_url = url
  )
  if (!is.data.frame(annotations) || nrow(annotations) < 1) {
    return(.plant_empty_occurrences())
  }
  rows = list()
  for (i in seq_len(nrow(annotations))) {
    annotation = annotations[i, , drop = FALSE]
    ids = .plant_pubchem_compound_links(annotation)
    if (nrow(ids) < 1) next
    for (j in seq_len(nrow(ids))) {
      compound_name = .uaf_first_non_empty_text(ids$label[[j]],
                                                annotation$CleanValue,
                                                annotation$Name)
      if (is.na(compound_name) ||
          .plant_text_has_taxon(compound_name, query_row)) {
        next
      }
      rows[[length(rows) + 1]] = data.frame(
        query_plant = query_row$query_plant,
        query_plant_clean = query_row$query_plant_clean,
        matched_taxon = query_row$species,
        matched_rank = "species",
        species = query_row$species,
        genus = query_row$genus,
        family = query_row$family,
        compound_name = compound_name,
        compound_name_clean = .plant_clean_compound(compound_name),
        compound_id = ids$cid[[j]],
        compound_id_type = "PubChem_CID",
        source_database = "PubChem Taxonomy",
        source_record_id = paste("taxonomy", taxid, "cid", ids$cid[[j]],
                                 sep = ":"),
        evidence_text = .plant_truncate(annotation$CleanValue, 800),
        evidence_url = .uaf_first_non_empty_text(annotation$SourceURL,
                                                annotation$PubChemURL,
                                                url),
        reference_id = .uaf_first_non_empty_text(annotation$Source,
                                                annotation$SourceURL),
	        pmid = .plant_explicit_pmid(c(annotation$CleanValue,
	                                      annotation$SourceURL)),
	        doi = .plant_explicit_doi(c(annotation$CleanValue,
	                                   annotation$SourceURL)),
	        plant_part = .plant_provider_text_context_value(
	          annotation$CleanValue, "plant_part"
	        ),
	        tissue = .plant_provider_text_context_value(annotation$CleanValue,
	                                                    "tissue"),
	        method = .plant_provider_text_context_value(annotation$CleanValue,
	                                                    "method"),
        occurrence_type = "pubchem_taxonomy_annotation",
        retrieved_at = .plant_timestamp(),
        confidence = "medium",
        curation_flag = "source_annotation_review_recommended",
        evidence_tier = "direct_species_database",
        stringsAsFactors = FALSE
      )
    }
  }
  occurrences = .plant_bind_occurrences(rows)
  max_records = .plant_max_records(max_records)
  if (is.finite(max_records) && nrow(occurrences) > max_records) {
    occurrences = occurrences[seq_len(max_records), , drop = FALSE]
  }
  occurrences
}

.plant_pubchem_taxonomy_external_specs = function(result) {
  cols = c("collection", "view", "srccmpdkind", "section_heading",
           "section_path", "description", "reference_number",
           "external_table_name", "occurrence_type", "confidence",
           "curation_flag")
  rows = list()
  walk_section = function(section, path = character()) {
    if (is.null(section) || !is.list(section)) return(invisible(NULL))
    heading = .uaf_first_non_empty_text(section$TOCHeading, section$Name)
    path_next = .uaf_non_empty(c(path, heading))
    infos = section$Information
    if (!is.null(infos) && is.list(infos)) {
      if (!is.null(infos$Value)) infos = list(infos)
      for (info in infos) {
        external_name = .uaf_first_non_empty_text(
          tryCatch(info$Value$ExternalTableName, error = function(error) NA_character_),
          tryCatch(info$ExternalTableName, error = function(error) NA_character_)
        )
        if (is.na(external_name)) next
        params = .plant_parse_query_string(external_name)
        collection = .uaf_first_non_empty_text(params$collection)
        kind = .uaf_first_non_empty_text(params$srccmpdkind)
        if (!identical(collection, "consolidatedcompoundtaxonomy") ||
            is.na(kind)) {
          next
        }
        rows[[length(rows) + 1]] <<- data.frame(
          collection = collection,
          view = .uaf_first_non_empty_text(params$view),
          srccmpdkind = kind,
          section_heading = .uaf_first_non_empty_text(heading),
          section_path = .pubchem_collapse(path_next),
          description = .uaf_first_non_empty_text(info$Description,
                                                  section$Description),
          reference_number = .uaf_first_non_empty_text(info$ReferenceNumber),
          external_table_name = external_name,
          occurrence_type = .plant_pubchem_external_occurrence_type(kind),
          confidence = "medium",
          curation_flag = "source_table_review_recommended",
          stringsAsFactors = FALSE
        )
      }
    }
    subsections = section$Section
    if (!is.null(subsections) && is.list(subsections)) {
      if (!is.null(subsections$TOCHeading) || !is.null(subsections$Information)) {
        subsections = list(subsections)
      }
      for (subsection in subsections) walk_section(subsection, path_next)
    }
    invisible(NULL)
  }
  sections = tryCatch(result$Record$Section, error = function(error) NULL)
  if (is.null(sections)) return(.uaf_empty_table(cols))
  if (!is.null(sections$TOCHeading) || !is.null(sections$Information)) {
    sections = list(sections)
  }
  for (section in sections) walk_section(section)
  .plant_bind_tables(rows, cols)
}

.plant_pubchem_external_occurrence_type = function(kind) {
  kind = tolower(.uaf_first_non_empty_text(kind))
  if (is.na(kind)) return("pubchem_taxonomy_external_table")
  if (kind == "metabolite") return("pubchem_taxonomy_metabolite_table")
  if (kind == "natural product") {
    return("pubchem_taxonomy_natural_product_table")
  }
  if (kind == "food compound") {
    return("pubchem_taxonomy_food_compound_table")
  }
  "pubchem_taxonomy_external_table"
}

.plant_pubchem_taxonomy_external_url = function(taxid, srccmpdkind,
                                                start = 1, limit = 100) {
  limit = suppressWarnings(as.integer(limit[[1]]))
  if (!is.finite(limit) || limit < 1) limit = 100L
  limit = min(limit, 500L)
  fields = c("cid", "synonym", "evids", "dsn", "cmpdname", "taxname",
             "pmids", "dois", "evurls", "citations", "srccmpdkind",
             "srccmpd", "srccmpdurl", "srcpart")
  query_json = paste0(
    "{\"select\":", jsonlite::toJSON(fields, auto_unbox = TRUE),
    ",\"collection\":\"consolidatedcompoundtaxonomy\"",
    ",\"order\":[\"cid,asc\"]",
    ",\"start\":", as.integer(start),
    ",\"limit\":", limit,
    ",\"where\":{\"ands\":[{\"taxid\":",
    jsonlite::toJSON(as.character(taxid), auto_unbox = TRUE),
    "},{\"srccmpdkind\":",
    jsonlite::toJSON(.uaf_first_non_empty_text(srccmpdkind),
                     auto_unbox = TRUE),
    "}]},\"width\":1000000}"
  )
  paste0("https://pubchem.ncbi.nlm.nih.gov/sdq/sphinxql.cgi?",
         "infmt=json&outfmt=json&query=",
         utils::URLencode(query_json, reserved = TRUE))
}

.plant_pubchem_taxonomy_external_rows = function(query_row, taxid, spec,
                                                 result, url, max_records) {
  rows_raw = .plant_pubchem_sphinx_rows(result)
  if (length(rows_raw) < 1) return(.plant_empty_occurrences())
  max_records = .plant_max_records(max_records)
  if (is.finite(max_records)) rows_raw = utils::head(rows_raw, max_records)
  rows = list()
  for (row in rows_raw) {
    if (is.data.frame(row)) row = as.list(row[1, , drop = FALSE])
    compound_name = .uaf_first_non_empty_text(row$cmpdname, row$synonym,
                                              row$srccmpd)
    if (is.na(compound_name) ||
        .plant_text_has_taxon(compound_name, query_row)) {
      next
    }
    cid = .uaf_first_non_empty_text(row$cid)
    source_reference = .plant_pubchem_external_reference(
      row = row,
      query_row = query_row,
      compound_name = compound_name
    )
    evidence_parts = c(
      paste0("PubChem taxonomy external table: ",
             .uaf_first_non_empty_text(spec$section_heading,
                                       spec$srccmpdkind)),
      paste0("source_compound_kind=", .uaf_first_non_empty_text(row$srccmpdkind,
                                                               spec$srccmpdkind)),
      paste0("source=", .uaf_first_non_empty_text(row$dsn)),
      paste0("taxon=", .uaf_first_non_empty_text(row$taxname)),
      paste0("evidence_ids=", .uaf_first_non_empty_text(row$evids)),
      paste0("source_pmid_count=", source_reference$pmid_count),
      paste0("source_doi_count=", source_reference$doi_count),
      paste0("source_literature_link_status=", source_reference$status)
    )
    # Only provider-supplied source-part fields support biological context.
    # Collection citations and database names are not plant-part evidence.
    text_for_context = .uaf_first_non_empty_text(row$srcpart)
    rows[[length(rows) + 1]] = data.frame(
      query_plant = query_row$query_plant,
      query_plant_clean = query_row$query_plant_clean,
      matched_taxon = query_row$species,
      matched_rank = "species",
      species = query_row$species,
      genus = query_row$genus,
      family = query_row$family,
      compound_name = compound_name,
      compound_name_clean = .plant_clean_compound(compound_name),
      compound_id = cid,
      compound_id_type = ifelse(is.na(cid), NA_character_, "PubChem_CID"),
      source_database = "PubChem Taxonomy",
      source_record_id = paste("taxonomy", taxid, "external",
                               .plant_matrix_key(spec$srccmpdkind), "cid",
                               .uaf_first_non_empty_text(cid, "unknown"),
                               sep = ":"),
      evidence_text = .plant_truncate(.pubchem_collapse(evidence_parts), 1200),
      evidence_url = .uaf_first_non_empty_text(
        .plant_first_pipe(row$evurls),
        row$srccmpdurl,
        ifelse(is.na(cid), NA_character_,
               paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", cid)),
        url
      ),
      reference_id = .uaf_first_non_empty_text(row$evids, row$dsn,
                                              spec$external_table_name),
      pmid = source_reference$pmid,
      doi = source_reference$doi,
      plant_part = .plant_provider_text_context_value(text_for_context,
                                                      "plant_part"),
      tissue = .plant_provider_text_context_value(text_for_context, "tissue"),
      method = .plant_provider_text_context_value(text_for_context, "method"),
      occurrence_type = spec$occurrence_type,
      retrieved_at = .plant_timestamp(),
      confidence = spec$confidence,
      curation_flag = .pubchem_collapse(c(
        spec$curation_flag,
        if (source_reference$status != "row_text_verified") {
          "source_literature_links_not_independently_assigned"
        }
      )),
      evidence_tier = "direct_species_database",
      stringsAsFactors = FALSE
    )
  }
  .plant_bind_occurrences(rows)
}

.plant_pubchem_external_reference = function(row, query_row, compound_name) {
  pmids = .plant_pipe_values(row$pmids)
  dois = .plant_pipe_values(row$dois)
  citation = .uaf_first_non_empty_text(row$citations)
  species = .uaf_first_non_empty_text(query_row$species,
                                      query_row$query_plant)
  citation_links_row = !is.na(citation) &&
    .plant_text_has_term(citation, species) &&
    .plant_text_has_compound(citation, compound_name)
  unambiguous = isTRUE(citation_links_row) &&
    length(pmids) <= 1L && length(dois) <= 1L &&
    (length(pmids) + length(dois)) > 0L
  list(
    pmid = if (unambiguous && length(pmids) == 1L) pmids[[1]] else
      NA_character_,
    doi = if (unambiguous && length(dois) == 1L) dois[[1]] else
      NA_character_,
    pmid_count = length(pmids),
    doi_count = length(dois),
    status = if (unambiguous) {
      "row_text_verified"
    } else if (length(pmids) + length(dois) > 0L || !is.na(citation)) {
      "not_independently_assigned"
    } else {
      "not_supplied"
    }
  )
}

.plant_pubchem_sphinx_rows = function(result) {
  sets = tryCatch(result$SDQOutputSet, error = function(error) NULL)
  if (is.null(sets)) return(list())
  if (is.data.frame(sets)) return(lapply(seq_len(nrow(sets)), function(i) {
    as.list(sets[i, , drop = FALSE])
  }))
  if (!is.null(sets$rows) || !is.null(sets$status)) sets = list(sets)
  out = list()
  for (set in sets) {
    status_code = suppressWarnings(as.integer(.uaf_first_non_empty_text(
      tryCatch(set$status$code, error = function(error) NA_character_)
    )))
    if (is.finite(status_code) && status_code != 0) next
    rows = tryCatch(set$rows, error = function(error) NULL)
    if (is.null(rows)) next
    if (is.data.frame(rows)) {
      out = c(out, lapply(seq_len(nrow(rows)), function(i) {
        as.list(rows[i, , drop = FALSE])
      }))
    } else if (is.list(rows)) {
      out = c(out, rows)
    }
  }
  out
}

.plant_parse_query_string = function(query) {
  query = .uaf_first_non_empty_text(query)
  if (is.na(query)) return(list())
  pieces = strsplit(query, "&", fixed = TRUE)[[1]]
  out = list()
  for (piece in pieces) {
    pair = strsplit(piece, "=", fixed = TRUE)[[1]]
    if (length(pair) < 2) next
    key = utils::URLdecode(pair[[1]])
    value = utils::URLdecode(paste(pair[-1], collapse = "="))
    value = gsub("^['\"]|['\"]$", "", trimws(value))
    out[[key]] = value
  }
  out
}

.plant_first_pipe = function(x) {
  .uaf_first_non_empty_text(.plant_pipe_values(x))
}

.plant_pipe_values = function(x) {
  text = .pubchem_collapse(x)
  if (is.na(text) || text == "") return(character())
  unique(.uaf_non_empty(unlist(strsplit(text, "\\|"), use.names = FALSE)))
}

.plant_explicit_pmid = function(x) {
  text = .pubchem_collapse(x)
  if (is.na(text) || text == "") return(NA_character_)
  patterns = c(
    "(?i)\\bPMID\\s*[:#]?\\s*(\\d{7,9})\\b",
    "(?i)pubmed\\.ncbi\\.nlm\\.nih\\.gov/(\\d{7,9})(?:/|\\b)"
  )
  for (pattern in patterns) {
    hit = regexec(pattern, text, perl = TRUE)
    values = regmatches(text, hit)[[1]]
    if (length(values) >= 2L) return(values[[2]])
  }
  NA_character_
}

.plant_explicit_doi = function(x) {
  text = .pubchem_collapse(x)
  if (is.na(text) || text == "") return(NA_character_)
  hit = regexpr("10\\.\\d{4,9}/[-._;()/:A-Za-z0-9]+", text,
                ignore.case = TRUE, perl = TRUE)
  if (hit[[1]] < 0) return(NA_character_)
  sub("[.,;:)]+$", "", regmatches(text, hit))
}

.plant_article_id = function(item, id_type) {
  articleids = tryCatch(item$articleids, error = function(error) NULL)
  if (is.null(articleids)) return(NA_character_)
  for (entry in articleids) {
    if (identical(tolower(.uaf_first_non_empty_text(entry$idtype)),
                  tolower(id_type))) {
      return(.uaf_first_non_empty_text(entry$value))
    }
  }
  NA_character_
}

.plant_pubtator_rows = function(query_row, result,
                                accepted_species = NULL) {
  docs = .plant_pubtator_docs(result)
  accepted_species = unique(.uaf_non_empty(c(query_row$species,
                                              accepted_species)))
  rows = list()
  for (doc in docs) {
    pmid = .uaf_first_non_empty_text(doc$pmid, doc$id, doc$sourceid)
    title = .uaf_first_non_empty_text(doc$title,
                                      .plant_pubtator_passage_text(doc,
                                                                  "title"))
    abstract = .uaf_first_non_empty_text(
      doc$abstract, .plant_pubtator_passage_text(doc, "abstract")
    )
    text = .pubchem_collapse(c(title, abstract, doc$text_hl,
                               .plant_pubtator_all_passage_text(doc)))
    chemicals = unique(.uaf_non_empty(c(
      .plant_pubtator_mentions(doc, c("Chemical", "chemical")),
      .plant_pubtator_highlight_mentions(doc, "CHEMICAL")
    )))
    species_mentions = unique(.uaf_non_empty(c(
      .plant_pubtator_mentions(doc, c("Species", "species", "Organism",
                                      "Taxon")),
      .plant_pubtator_highlight_mentions(doc, "SPECIES")
    )))
    if (length(species_mentions) < 1 &&
        .plant_text_contains(text, accepted_species)) {
      species_mentions = accepted_species[vapply(
        accepted_species, function(value) .plant_text_contains(text, value),
        logical(1)
      )]
    }
    if (length(species_mentions) < 1 ||
        !.plant_text_contains(.pubchem_collapse(species_mentions),
                              accepted_species)) {
      next
    }
    matched_species = accepted_species[vapply(
      accepted_species,
      function(value) .plant_text_contains(
        .pubchem_collapse(species_mentions), value
      ), logical(1)
    )]
    if (length(chemicals) < 1) chemicals = NA_character_
    for (chemical in chemicals) {
      rows[[length(rows) + 1]] = data.frame(
        query_plant = query_row$query_plant,
        query_plant_clean = query_row$query_plant_clean,
        species = query_row$species,
        genus = query_row$genus,
        family = query_row$family,
        source_database = "PubTator",
        source_record_id = pmid,
        pmid = pmid,
        doi = NA_character_,
        title = title,
        abstract = abstract,
        chemical_mention = chemical,
        species_mention = .pubchem_collapse(species_mentions),
        evidence_text = .plant_pubtator_evidence_snippet(
          text, .uaf_first_non_empty_text(matched_species,
                                          query_row$species), chemical
        ),
        evidence_url = ifelse(is.na(pmid), NA_character_,
                              paste0("https://pubmed.ncbi.nlm.nih.gov/",
                                     pmid, "/")),
        retrieved_at = .plant_timestamp(),
        confidence = "low",
        curation_flag = "candidate_co_mention",
        evidence_tier = "direct_species_pubtator_candidate",
        stringsAsFactors = FALSE
      )
    }
  }
  list(LiteratureCandidates = .plant_bind_tables(rows, .plant_literature_cols()))
}

.plant_pubtator_docs = function(result) {
  if (is.null(result)) return(list())
  if (!is.null(result$PubTator3) && is.list(result$PubTator3)) {
    return(result$PubTator3)
  }
  if (!is.null(result$documents) && is.list(result$documents)) {
    return(result$documents)
  }
  if (!is.null(result$results) && is.list(result$results)) {
    return(result$results)
  }
  if (!is.null(result$passages) || !is.null(result$annotations)) {
    return(list(result))
  }
  if (is.list(result) && length(result) > 0 && is.null(names(result))) {
    is_doc = vapply(result, function(item) {
      is.list(item) && (!is.null(item$id) || !is.null(item$pmid) ||
                          !is.null(item$passages) || !is.null(item$annotations))
    }, logical(1))
    if (all(is_doc)) return(result)
  }
  list()
}

.plant_pubtator_all_passage_text = function(doc) {
  passages = tryCatch(doc$passages, error = function(error) NULL)
  if (!is.list(passages)) return(NA_character_)
  .pubchem_collapse(vapply(passages, function(passage) {
    .uaf_first_non_empty_text(passage$text)
  }, character(1)))
}

.plant_pubtator_passage_text = function(doc, type) {
  passages = tryCatch(doc$passages, error = function(error) NULL)
  if (!is.list(passages)) return(NA_character_)
  values = vapply(passages, function(passage) {
    passage_type = tolower(.uaf_first_non_empty_text(
      passage$infons$type, passage$infons$section_type, passage$type, ""
    ))
    if (!identical(passage_type, tolower(type))) return(NA_character_)
    .uaf_first_non_empty_text(passage$text)
  }, character(1))
  .uaf_first_non_empty_text(values)
}

.plant_pubtator_evidence_snippet = function(text, species, chemical) {
  text = .uaf_first_non_empty_text(text)
  if (is.na(text)) return(NA_character_)
  sentences = unlist(strsplit(text, "(?<=[.!?])\\s+", perl = TRUE),
                     use.names = FALSE)
  has_species = vapply(sentences, .plant_text_contains, logical(1),
                       value = species)
  if (!is.na(chemical)) {
    has_chemical = vapply(sentences, .plant_text_contains, logical(1),
                          value = chemical)
    both = sentences[has_species & has_chemical]
    if (length(.uaf_non_empty(both)) > 0) {
      return(.plant_truncate(.uaf_first_non_empty_text(both), 1200))
    }
  }
  species_sentences = sentences[has_species]
  .plant_truncate(.uaf_first_non_empty_text(species_sentences, text), 1200)
}

.plant_pubtator_mentions = function(doc, types) {
  annotations = doc$annotations
  if (is.null(annotations) && !is.null(doc$passages)) {
    annotations = unlist(lapply(doc$passages, function(passage) {
      passage$annotations
    }), recursive = FALSE)
  }
  if (is.null(annotations)) return(character())
  values = character()
  for (ann in annotations) {
    ann_type = .uaf_first_non_empty_text(ann$infons$type, ann$type,
                                         ann$infons$identifier)
    if (!ann_type %in% types) next
    values = c(values, .uaf_first_non_empty_text(ann$text))
  }
  unique(.uaf_non_empty(values))
}

.plant_pubtator_highlight_mentions = function(doc, type) {
  text = .uaf_first_non_empty_text(doc$text_hl, doc$highlight,
                                   doc$snippet)
  if (is.na(text) || text == "") return(character())
  pattern = paste0("@", type, "_[^[:space:]]+\\s+@@@(.*?)@@@")
  matches = gregexpr(pattern, text, perl = TRUE, ignore.case = TRUE)
  raw = regmatches(text, matches)[[1]]
  if (length(raw) < 1 || identical(raw, "")) return(character())
  values = vapply(raw, function(item) {
    out = sub(pattern, "\\1", item, perl = TRUE, ignore.case = TRUE)
    .plant_html_text(out)
  }, character(1))
  unique(.uaf_non_empty(values))
}

.plant_occurrences_from_pubtator = function(literature) {
  literature = .plant_normalize_literature(literature, "pubtator")
  if (nrow(literature) < 1) return(.plant_empty_occurrences())
  keep = !is.na(literature$chemical_mention) & literature$chemical_mention != ""
  literature = literature[keep, , drop = FALSE]
  if (nrow(literature) < 1) return(.plant_empty_occurrences())
  out = data.frame(
    query_plant = literature$query_plant,
    query_plant_clean = literature$query_plant_clean,
    matched_taxon = literature$species,
    matched_rank = "species",
    species = literature$species,
    genus = literature$genus,
    family = literature$family,
    compound_name = literature$chemical_mention,
    compound_name_clean = .plant_clean_compound(literature$chemical_mention),
    compound_id = NA_character_,
    compound_id_type = NA_character_,
    source_database = "PubTator",
    source_record_id = literature$source_record_id,
    evidence_text = literature$evidence_text,
    evidence_url = literature$evidence_url,
    reference_id = literature$source_record_id,
    pmid = literature$pmid,
    doi = literature$doi,
    plant_part = NA_character_,
    tissue = NA_character_,
    method = NA_character_,
    occurrence_type = "candidate_co_mention",
    retrieved_at = literature$retrieved_at,
    confidence = "low",
    curation_flag = "candidate_co_mention",
    evidence_tier = "direct_species_pubtator_candidate",
    stringsAsFactors = FALSE
  )
  .plant_normalize_occurrences(out, "pubtator")
}

.plant_json_records = function(x, depth = 0) {
  if (depth > 5 || is.null(x)) return(list())
  if (is.data.frame(x)) {
    return(lapply(seq_len(nrow(x)), function(i) {
      as.list(x[i, , drop = FALSE])
    }))
  }
  if (!is.list(x)) return(list())
  names_x = .plant_normalize_column_names(names(x))
  record_keys = c("name", "title", "compound", "molecule", "smiles", "inchi",
                  "inchikey", "lotus", "wikidata", "organism", "taxon",
                  "species", "genus", "family")
  if (any(grepl(paste(record_keys, collapse = "|"), names_x,
                ignore.case = TRUE))) {
    return(list(x))
  }
  out = list()
  for (item in x) {
    out = c(out, .plant_json_records(item, depth + 1))
  }
  out
}

.plant_flatten_record = function(x, prefix = character(), depth = 0) {
  if (depth > 5 || is.null(x)) return(character())
  if (is.data.frame(x)) {
    if (nrow(x) < 1) return(character())
    x = as.list(x[1, , drop = FALSE])
  }
  if (is.list(x)) {
    out = character()
    names_x = names(x)
    if (is.null(names_x)) names_x = paste0("item", seq_along(x))
    for (i in seq_along(x)) {
      name = .plant_normalize_column_names(names_x[[i]])
      if (is.na(name) || name == "") name = paste0("item", i)
      child = .plant_flatten_record(x[[i]], c(prefix, name), depth + 1)
      out = c(out, child)
    }
    return(out)
  }
  value = .uaf_squish_text(as.character(unlist(x, use.names = FALSE)))
  value = .uaf_non_empty(value)
  if (length(value) < 1) return(character())
  key = paste(prefix, collapse = "_")
  if (key == "") key = "value"
  stats::setNames(value, rep(key, length(value)))
}

.plant_lotus_compound_name = function(flat, query_row) {
  priority = c("compound_name", "molecule_name", "natural_product_name",
               "traditional_name", "iupac_name", "name", "title")
  keys = names(flat)
  for (pattern in priority) {
    values = flat[grepl(pattern, keys, ignore.case = TRUE)]
    values = .uaf_non_empty(values)
    values = values[!vapply(values, .plant_text_has_taxon, logical(1),
                            query_row = query_row)]
    values = values[!grepl("^LTS\\d+|^Q\\d+$", values, ignore.case = TRUE)]
    if (length(values) > 0) return(values[[1]])
  }
  NA_character_
}

.plant_lotus_record_id = function(flat) {
  values = flat[grepl("lotus|wikidata|record|identifier|^id$",
                      names(flat), ignore.case = TRUE)]
  values = .uaf_non_empty(values)
  hit = values[grepl("^LTS\\d+|^Q\\d+$|^https?://", values,
                     ignore.case = TRUE)]
  .uaf_first_non_empty_text(hit, values)
}

.plant_knapsack_compound_name = function(cells, query_row) {
  cells = .uaf_non_empty(.uaf_squish_text(cells))
  if (length(cells) < 1) return(NA_character_)
  reject = grepl(
    paste(c("^C\\d{6,}$", "^CAS$", "^C_ID$", "metabolite", "organism",
            "formula", "molecular", "^mw$", "exact", "inchi", "smiles",
            "search", "image", "link", "name", "taxonomy", "all rights"),
          collapse = "|"),
    cells,
    ignore.case = TRUE
  )
  taxon_cell = vapply(cells, .plant_text_has_taxon, logical(1),
                      query_row = query_row)
  reject = reject | .plant_formula_like(cells) | taxon_cell
  candidates = cells[!reject]
  candidates = candidates[grepl("[A-Za-z]", candidates)]
  candidates = candidates[nchar(candidates) >= 3]
  .uaf_first_non_empty_text(candidates)
}

.plant_knapsack_detail_url = function(hrefs, c_id) {
  hrefs = .uaf_non_empty(.uaf_squish_text(hrefs))
  c_id = .uaf_first_non_empty_text(c_id)
  if (length(hrefs) > 0) {
    hits = hrefs[grepl("information\\.php", hrefs, ignore.case = TRUE)]
    if (!is.na(c_id) && length(hits) > 0) {
      c_hits = hits[grepl(c_id, hits, fixed = TRUE)]
      if (length(c_hits) > 0) hits = c_hits
    }
    hit = .uaf_first_non_empty_text(hits)
    if (!is.na(hit)) return(.plant_knapsack_absolute_url(hit))
  }
  if (is.na(c_id)) return(NA_character_)
  paste0("https://www.knapsackfamily.com/knapsack_core/information.php?word=",
         utils::URLencode(c_id, reserved = TRUE))
}

.plant_knapsack_absolute_url = function(url) {
  url = .uaf_first_non_empty_text(url)
  if (is.na(url)) return(NA_character_)
  url = gsub("&amp;", "&", url, ignore.case = TRUE)
  if (grepl("^https?://", url, ignore.case = TRUE)) return(url)
  if (startsWith(url, "/")) {
    return(paste0("https://www.knapsackfamily.com", url))
  }
  paste0("https://www.knapsackfamily.com/knapsack_core/", url)
}

.plant_knapsack_evidence_text = function(cells, detail_url = NA_character_) {
  cells = .uaf_non_empty(.uaf_squish_text(cells))
  if (length(cells) < 1) return(NA_character_)
  labels = switch(
    as.character(length(cells)),
    "4" = c("C_ID", "Metabolite", "Molecular_formula", "Organism"),
    "5" = c("C_ID", "CAS_ID", "Metabolite", "Molecular_formula",
            "Organism"),
    "6" = c("C_ID", "CAS_ID", "Metabolite", "Molecular_formula",
            "Molecular_weight", "Organism"),
    paste0("cell_", seq_along(cells))
  )
  if (length(labels) != length(cells)) {
    labels = paste0("cell_", seq_along(cells))
  }
  fields = paste0(labels, "=", cells)
  detail_url = .uaf_first_non_empty_text(detail_url)
  if (!is.na(detail_url)) fields = c(fields, paste0("detail_url=", detail_url))
  .plant_truncate(paste(fields, collapse = "; "), 800)
}

.plant_provider_context_cell_value = function(cells,
                                              kind = c("plant_part", "tissue",
                                                       "method")) {
  kind = match.arg(kind)
  cells = .uaf_non_empty(.uaf_squish_text(cells))
  if (length(cells) < 1) return(NA_character_)
  for (cell in cells) {
    if (kind == "method") {
      group = .plant_method_group(cell)
      reject = c("unknown", "other", "database_record", "literature_curation")
    } else {
      group = .plant_context_group(cell, kind = kind)
      reject = c("unknown", "other", "extract_unspecified")
    }
    if (!group %in% reject) return(cell)
  }
  NA_character_
}

.plant_provider_text_context_value = function(text,
                                              kind = c("plant_part", "tissue",
                                                       "method")) {
  kind = match.arg(kind)
  text = .uaf_first_non_empty_text(text)
  if (is.na(text) || text == "") return(NA_character_)
  patterns = .plant_context_patterns()
  patterns = patterns[patterns$context_type == kind, , drop = FALSE]
  for (i in seq_len(nrow(patterns))) {
    hit = .plant_regex_match_text(text, patterns$pattern[[i]])
    if (is.na(hit) || hit == "") next
    if (kind == "method") {
      group = .plant_method_group(hit)
      reject = c("unknown", "other", "database_record", "literature_curation")
    } else {
      group = .plant_context_group(hit, kind = kind)
      method_group = .plant_method_group(hit)
      method_like = !method_group %in% c("unknown", "other",
                                         "database_record",
                                         "literature_curation")
      reject = c("unknown", "other", "extract_unspecified")
      if (method_like) next
    }
    if (!group %in% reject) return(hit)
  }
  NA_character_
}

.plant_context_raw_method_like = function(x) {
  x = .uaf_squish_text(x)
  method_group = .plant_method_group(x)
  !is.na(x) & x != "" &
    !method_group %in% c("unknown", "other", "database_record",
                         "literature_curation")
}

.plant_clean_context_conflicts = function(occurrences) {
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) return(occurrences)
  invalid_part = occurrences$plant_part_group %in%
    c("unknown", "other", "extract_unspecified")
  invalid_tissue = occurrences$tissue_group %in%
    c("unknown", "other", "extract_unspecified")
  method_like_part = .plant_context_raw_method_like(occurrences$plant_part)
  method_like_tissue = .plant_context_raw_method_like(occurrences$tissue)
  replace_part = method_like_part & !invalid_part
  replace_tissue = method_like_tissue & !invalid_tissue
  occurrences$plant_part[replace_part] =
    occurrences$plant_part_group[replace_part]
  occurrences$tissue[replace_tissue] =
    occurrences$tissue_group[replace_tissue]
  occurrences$plant_part[method_like_part & invalid_part] = NA_character_
  occurrences$tissue[method_like_tissue & invalid_tissue] = NA_character_
  occurrences$biological_context_status = .plant_context_status(
    occurrences$plant_part_group,
    occurrences$tissue_group,
    occurrences$method_group
  )
  occurrences
}

.plant_html_rows = function(html) {
  records = .plant_html_row_records(html)
  lapply(records, function(record) record$cells)
}

.plant_html_row_records = function(html) {
  html = paste(as.character(html), collapse = "\n")
  rows = regmatches(html, gregexpr("<tr[^>]*>.*?</tr>", html,
                                   perl = TRUE, ignore.case = TRUE))[[1]]
  out = list()
  if (length(rows) > 0 && !identical(rows, "")) {
    for (row in rows) {
      cells = regmatches(row, gregexpr("<t[dh][^>]*>.*?</t[dh]>", row,
                                       perl = TRUE, ignore.case = TRUE))[[1]]
      cells = .uaf_non_empty(.plant_html_text(cells))
      if (length(cells) > 0) {
        out[[length(out) + 1]] = list(
          cells = cells,
          hrefs = .plant_html_hrefs(row)
        )
      }
    }
  }
  if (length(out) > 0) return(out)
  text_rows = strsplit(.plant_html_text(html), "\n", fixed = TRUE)[[1]]
  text_rows = .uaf_non_empty(text_rows)
  for (row in text_rows) {
    cells = .uaf_non_empty(strsplit(row, "\\s{2,}|\\t", perl = TRUE)[[1]])
    if (length(cells) > 1) {
      out[[length(out) + 1]] = list(cells = cells, hrefs = character())
    }
  }
  out
}

.plant_html_hrefs = function(html) {
  matches = regmatches(
    html,
    gregexpr("href\\s*=\\s*['\"][^'\"]+['\"]|href\\s*=\\s*[^[:space:]>]+",
             html, perl = TRUE, ignore.case = TRUE)
  )[[1]]
  if (length(matches) < 1 || identical(matches, "")) return(character())
  matches = sub("^href\\s*=\\s*", "", matches, ignore.case = TRUE)
  matches = gsub("^[\"']|[\"']$", "", matches)
  matches = gsub("&amp;", "&", matches, ignore.case = TRUE)
  .uaf_non_empty(.uaf_squish_text(matches))
}

.plant_html_text = function(x) {
  x = gsub("<script[^>]*>.*?</script>", " ", x, ignore.case = TRUE,
           perl = TRUE)
  x = gsub("<style[^>]*>.*?</style>", " ", x, ignore.case = TRUE,
           perl = TRUE)
  x = gsub("<[^>]+>", " ", x, perl = TRUE)
  x = gsub("&nbsp;", " ", x, ignore.case = TRUE)
  x = gsub("&amp;", "&", x, ignore.case = TRUE)
  x = gsub("&quot;", "\"", x, ignore.case = TRUE)
  x = gsub("&#39;|&apos;", "'", x, ignore.case = TRUE)
  .uaf_squish_text(x)
}

.plant_pubchem_compound_links = function(annotation) {
  urls = .pubchem_split_collapsed(annotation$MarkupURL)
  labels = .pubchem_split_collapsed(annotation$MarkupText)
  rows = list()
  for (i in seq_along(urls)) {
    cid = .plant_pubchem_cid_from_url(urls[[i]])
    if (is.na(cid)) next
    label = if (length(labels) >= i) labels[[i]] else NA_character_
    rows[[length(rows) + 1]] = data.frame(
      cid = cid,
      label = label,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) < 1) return(.uaf_empty_table(c("cid", "label")))
  unique(do.call(rbind, rows))
}

.plant_pubchem_cid_from_url = function(url) {
  url = .uaf_first_non_empty_text(url)
  if (is.na(url)) return(NA_character_)
  patterns = c("/compound/(\\d+)", "[?&]cid=(\\d+)", "/cid/(\\d+)",
               "CID(\\d+)")
  for (pattern in patterns) {
    hit = regexec(pattern, url, ignore.case = TRUE)
    match = regmatches(url, hit)[[1]]
    if (length(match) > 1) return(match[[2]])
  }
  NA_character_
}

.plant_matched_rank_from_text = function(query_row, text,
                                         accepted_species = NULL) {
  species = unique(.uaf_non_empty(c(query_row$species, accepted_species)))
  if (.plant_text_contains(text, species)) return("species")
  if (.plant_text_contains(text, query_row$genus)) return("genus")
  if (.plant_text_contains(text, query_row$family)) return("family")
  NA_character_
}

.plant_matched_species_from_text = function(text, accepted_species,
                                            fallback) {
  candidates = unique(.uaf_non_empty(c(accepted_species, fallback)))
  if (length(candidates) < 1) return(NA_character_)
  hit = vapply(candidates, function(value) {
    .plant_text_contains(text, value)
  }, logical(1))
  if (!any(hit)) return(.uaf_first_non_empty_text(fallback))
  candidates[hit][order(nchar(candidates[hit]), decreasing = TRUE)][[1]]
}

.plant_text_has_taxon = function(text, query_row) {
  .plant_text_contains(text, query_row$species) ||
    .plant_text_contains(text, query_row$genus) ||
    .plant_text_contains(text, query_row$family)
}

.plant_text_contains = function(text, value) {
  text = .plant_clean_name(.pubchem_collapse(text))
  value = .plant_clean_name(value)
  value = .uaf_non_empty(value)
  if (length(value) < 1 || is.na(text) || text == "") return(FALSE)
  any(vapply(value, function(v) {
    !is.na(v) && v != "" && grepl(v, text, fixed = TRUE)
  }, logical(1)))
}

.plant_formula_like = function(x) {
  x = .uaf_squish_text(x)
  grepl("^[A-Z][a-z]?\\d*([A-Z][a-z]?\\d*)+$", x) |
    grepl("^[0-9.]+$", x)
}

.plant_first_named = function(flat, pattern) {
  values = flat[grepl(pattern, names(flat), ignore.case = TRUE)]
  .uaf_first_non_empty_text(values)
}

.plant_first_named_filtered = function(flat, include, exclude = NA_character_) {
  nm = names(flat)
  keep = grepl(include, nm, ignore.case = TRUE, perl = TRUE)
  if (!is.na(exclude)) {
    keep = keep & !grepl(exclude, nm, ignore.case = TRUE, perl = TRUE)
  }
  values = flat[keep]
  if (length(values) < 1) return(NA_character_)
  for (value in unname(values)) {
    value = .uaf_first_non_empty_text(value)
    if (is.na(value) || value == "") next
    if (grepl("^[0-9]+$", value)) next
    if (grepl("^https?://", value, ignore.case = TRUE)) next
    return(value)
  }
  NA_character_
}

.plant_lotus_context_exclude = function() {
  paste(c("organism", "taxonomy", "taxon", "kingdom", "phylum", "classx",
          "family", "genus", "species", "cleaned_organism_id", "wikidata",
          "reference", "pubchemfingerprint", "fragment",
          "chemicaltaxonomy"), collapse = "|")
}

.plant_lotus_plant_part = function(flat) {
  .plant_first_named_filtered(
    flat,
    paste(c("(^|[._])plant[_ ]?part([._]|$)",
            "(^|[._])plantpart([._]|$)",
            "(^|[._])part[_ ]?used([._]|$)",
            "(^|[._])partused([._]|$)",
            "(^|[._])sample[_ ]?organ([._]|$)",
            "(^|[._])sample[_ ]?part([._]|$)",
            "(^|[._])material([._]|$)"), collapse = "|"),
    .plant_lotus_context_exclude()
  )
}

.plant_lotus_tissue = function(flat) {
  .plant_first_named_filtered(
    flat,
    paste(c("(^|[._])tissue([._]|$)",
            "(^|[._])sample[_ ]?tissue([._]|$)"), collapse = "|"),
    .plant_lotus_context_exclude()
  )
}

.plant_lotus_method = function(flat) {
  .plant_first_named_filtered(
    flat,
    paste(c("(^|[._])method([._]|$)",
            "(^|[._])analysis([._]|$)",
            "(^|[._])technique([._]|$)",
            "(^|[._])instrument([._]|$)"), collapse = "|"),
    "chemicaltaxonomy|pubchemfingerprint|fragment"
  )
}

.plant_first_pattern = function(values, pattern) {
  text = .pubchem_collapse(unlist(values, use.names = FALSE))
  if (is.na(text) || text == "") return(NA_character_)
  hit = regexpr(pattern, text, perl = TRUE, ignore.case = TRUE)
  if (is.na(hit[[1]]) || hit[[1]] < 0) return(NA_character_)
  regmatches(text, hit)
}

.plant_max_records = function(max_records) {
  max_records = suppressWarnings(as.numeric(max_records[[1]]))
  if (!is.finite(max_records) || max_records < 1) return(Inf)
  max_records
}

.plant_remaining_records = function(parts, max_records) {
  max_records = .plant_max_records(max_records)
  if (!is.finite(max_records)) return(Inf)
  used = sum(vapply(parts, function(x) {
    if (is.data.frame(x)) nrow(x) else 0L
  }, integer(1)), na.rm = TRUE)
  max(0L, as.integer(max_records - used))
}

.plant_truncate = function(x, n) {
  x = .uaf_first_non_empty_text(x)
  if (is.na(x) || nchar(x) <= n) return(x)
  paste0(substr(x, 1, n - 3), "...")
}

.plant_normalize_column_names = function(x) {
  x = gsub("[^A-Za-z0-9]+", "_", x)
  x = gsub("_+", "_", x)
  x = gsub("^_|_$", "", x)
  tolower(x)
}

.plant_col_or_default = function(table, column, default) {
  n = if (is.data.frame(table)) nrow(table) else 0
  if (!is.data.frame(table) || !column %in% names(table)) {
    return(rep(default, n))
  }
  value = table[[column]]
  value = .uaf_squish_text(value)
  value[is.na(value) | value == ""] = default
  value
}

.plant_clean_name = function(x) {
  x = tolower(.uaf_squish_text(x))
  x = gsub("\\s+", " ", x)
  x
}

.plant_clean_compound = function(x) {
  x = .uaf_squish_text(x)
  plus_minus = intToUtf8(0x00b1)
  x = gsub(plus_minus, " plus_minus ", x, fixed = TRUE)
  x = gsub("\\(\\s*(\\+-|\\+\\s*/\\s*-|-\\s*/\\s*\\+)\\s*\\)",
           " plus_minus ", x, perl = TRUE)
  x = gsub("\\(\\s*\\+\\s*\\)", " plus ", x, perl = TRUE)
  x = gsub("\\(\\s*-\\s*\\)", " minus ", x, perl = TRUE)
  x = tolower(x)
  greek = stats::setNames(
    c(" alpha ", " beta ", " gamma ", " delta ", " epsilon ", " zeta ",
      " eta ", " theta ", " lambda ", " mu ", " pi ", " sigma ",
      " omega "),
    vapply(c(0x03b1, 0x03b2, 0x03b3, 0x03b4, 0x03b5, 0x03b6,
             0x03b7, 0x03b8, 0x03bb, 0x03bc, 0x03c0, 0x03c3,
             0x03c9), intToUtf8, character(1))
  )
  for (pattern in names(greek)) {
    x = gsub(pattern, greek[[pattern]], x, fixed = TRUE)
  }
  x = gsub("[^a-z0-9]+", "_", x)
  x = gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

.plant_slug = function(x) {
  x = .plant_clean_compound(x)
  x[is.na(x) | x == ""] = NA_character_
  x
}

.plant_genus = function(species) {
  species = .uaf_squish_text(species)
  vapply(strsplit(species, "\\s+"), function(parts) {
    if (length(parts) < 1) return(NA_character_)
    parts[[1]]
  }, character(1))
}

.plant_canonical_taxon_name = function(x) {
  x = .uaf_squish_text(x)
  vapply(x, function(item) {
    if (is.na(item) || !nzchar(item)) return(NA_character_)
    parts = strsplit(item, "\\s+")[[1]]
    if (length(parts) < 1) return(NA_character_)
    parts = parts[nzchar(parts)]
    if (length(parts) < 1) return(NA_character_)
    parts[[1]] = .plant_title_word(parts[[1]])
    if (length(parts) > 1) {
      parts[-1] = tolower(parts[-1])
    }
    paste(parts, collapse = " ")
  }, character(1), USE.NAMES = FALSE)
}

.plant_title_word = function(x) {
  x = .uaf_squish_text(x)
  out = x
  ok = !is.na(x) & nzchar(x)
  out[ok] = paste0(toupper(substr(x[ok], 1, 1)),
                   tolower(substr(x[ok], 2, nchar(x[ok]))))
  out
}

.plant_is_binomial = function(x) {
  x = .uaf_squish_text(x)
  vapply(strsplit(x, "\\s+"), function(parts) {
    parts = parts[nzchar(parts)]
    if (length(parts) < 2) return(FALSE)
    genus_ok = grepl("^[A-Z][A-Za-z-]+$", parts[[1]])
    if (!genus_ok) return(FALSE)
    marker = tolower(gsub("[.]$", "", parts[[2]]))
    if (marker %in% c("sp", "spp", "cf", "aff")) return(FALSE)
    epithet = parts[[2]]
    if (identical(marker, "x")) {
      if (length(parts) < 3) return(FALSE)
      epithet = parts[[3]]
    }
    epithet = tolower(gsub("[.]$", "", epithet))
    grepl("^[a-z][a-z-]+$", epithet) &&
      !(epithet %in% c("sp", "spp", "cf", "aff"))
  }, logical(1), USE.NAMES = FALSE)
}

.plant_cache_dir = function(cache_dir) {
  if (!is.null(cache_dir) && length(.uaf_non_empty(cache_dir)) > 0) {
    return(.uaf_non_empty(cache_dir)[[1]])
  }
  out = tryCatch(tools::R_user_dir("uafR", "cache"),
                 error = function(error) NULL)
  if (is.null(out)) out = file.path(tempdir(), "uafR-plant-cache")
  file.path(out, "plant_phytochemistry")
}

.plant_timestamp = function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

.plant_normalize_evidence_tier = function(evidence_tier, source_database,
                                          matched_rank) {
  evidence_tier = .uaf_squish_text(evidence_tier)
  source_database = tolower(.uaf_squish_text(source_database))
  matched_rank = tolower(.uaf_squish_text(matched_rank))
  out = evidence_tier
  missing = is.na(out) | out == ""
  out[missing & source_database %in% c("manual", "curated")] = "manual_curated"
  out[missing & source_database %in% c("pubtator")] =
    "direct_species_pubtator_candidate"
  out[missing & source_database %in% c("pubmed")] = "direct_species_literature"
  out[missing & matched_rank == "genus"] = "genus_database_fallback"
  out[missing & matched_rank == "family"] = "family_database_fallback"
  out[missing] = "direct_species_database"
  out[!out %in% .plant_evidence_tiers()] = "unresolved"
  out
}

.plant_evidence_tier_confidence = function(evidence_tier) {
  tier = .plant_normalize_evidence_tier(evidence_tier, "", "species")
  ifelse(tier %in% c("direct_species_database", "manual_curated"), "high",
         ifelse(tier %in% c("direct_species_literature",
                            "genus_database_fallback",
                            "genus_literature_fallback"), "medium", "low"))
}

.plant_normalize_confidence = function(confidence, evidence_tier) {
  confidence = tolower(.uaf_squish_text(confidence))
  fallback = .plant_evidence_tier_confidence(evidence_tier)
  confidence[is.na(confidence) | confidence == ""] =
    fallback[is.na(confidence) | confidence == ""]
  confidence[!confidence %in% .plant_confidence_values()] = "unknown"
  confidence
}

.plant_confidence_score = function(confidence) {
  confidence = tolower(.uaf_squish_text(confidence))
  if (length(confidence) == 1 && confidence %in% c("low", "medium", "high",
                                                   "unknown")) {
    return(switch(confidence, low = 0.35, medium = 0.65, high = 1,
                  unknown = 0))
  }
  out = suppressWarnings(as.numeric(confidence))
  out[confidence == "low"] = 0.35
  out[confidence == "medium"] = 0.65
  out[confidence == "high"] = 1
  out[confidence == "unknown" | is.na(confidence)] = 0
  out[is.na(out)] = 0
  out
}

.plant_clean_context_value = function(x) {
  out = .uaf_squish_text(x)
  lower = tolower(out)
  out[is.na(out) | out == "" |
        lower %in% c("na", "n/a", "none", "unknown", "unspecified",
                     "not available", "not reported")] = NA_character_
  out[grepl("^[0-9]+$", out)] = NA_character_
  out[grepl("^https?://", out, ignore.case = TRUE)] = NA_character_
  out
}

.plant_normalize_context_evidence = function(x) {
  cols = .plant_context_evidence_cols()
  if (!is.data.frame(x) || nrow(x) < 1) return(.uaf_empty_table(cols))
  names(x) = .plant_normalize_column_names(names(x))
  for (col in cols) if (!col %in% names(x)) x[[col]] = NA_character_
  x = x[, cols, drop = FALSE]
  x[cols] = lapply(x[cols], .uaf_squish_text)
  missing_slug = is.na(x$species_slug) | x$species_slug == ""
  x$species_slug[missing_slug] = .plant_slug(x$species[missing_slug])
  x$compound_name_clean = .plant_clean_compound(
    ifelse(is.na(x$compound_name) | x$compound_name == "",
           x$compound_name_clean, x$compound_name)
  )
  x$context_type = .plant_matrix_key(x$context_type)
  x$context_type[!x$context_type %in% .plant_context_type_values()] =
    "plant_part"
  x$normalized_context = .plant_matrix_key(x$normalized_context)
  x$context_confidence =
    .plant_comparability_clean_confidence(x$context_confidence)
  review = x$context_confidence == "low" |
    x$normalized_context %in% c("other", "extract_unspecified", "unknown")
  x$requires_review = .uaf_yes_no(review)
  x$retrieved_at[is.na(x$retrieved_at) | x$retrieved_at == ""] =
    .plant_timestamp()
  x = x[!is.na(x$species) & x$species != "" &
          !is.na(x$compound_name_clean) & x$compound_name_clean != "" &
          !is.na(x$normalized_context) & x$normalized_context != "" &
          x$normalized_context != "unknown", , drop = FALSE]
  x = unique(x)
  x = .plant_collapse_context_evidence_keys(x)
  row.names(x) = NULL
  x
}

.plant_collapse_context_evidence_keys = function(x) {
  if (!is.data.frame(x) || nrow(x) < 2) return(x)
  key_cols = .plant_duplicate_key_cols("PlantContextEvidence")
  if (!all(key_cols %in% names(x))) return(x)
  key_data = x[key_cols]
  key_data[] = lapply(key_data, function(value) {
    value = .uaf_squish_text(value)
    value[is.na(value)] = ""
    value
  })
  key = do.call(paste, c(key_data, sep = "||"))
  duplicate = duplicated(key) | duplicated(key, fromLast = TRUE)
  if (!any(duplicate)) {
    out = x[order(key), , drop = FALSE]
    row.names(out) = NULL
    return(out)
  }
  groups = split(which(duplicate), key[duplicate])
  rows = lapply(groups, function(idx) {
    .plant_merge_context_evidence_group(x[idx, , drop = FALSE])
  })
  unique_idx = which(!duplicate)
  merged = do.call(rbind, rows)
  out = rbind(x[unique_idx, , drop = FALSE], merged)
  out_key = c(key[unique_idx], names(groups))
  out = out[order(out_key), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_merge_context_evidence_group = function(group) {
  basis_rank = .plant_context_basis_rank(group$evidence_basis)
  confidence_score = .plant_confidence_score(group$context_confidence)
  best = which.max(confidence_score + basis_rank)
  row = group[best, , drop = FALSE]
  row$raw_context_text = .plant_truncate(
    .pubchem_collapse(unique(.uaf_non_empty(group$raw_context_text))),
    800
  )
  row$extraction_rule = .plant_truncate(
    .pubchem_collapse(unique(.uaf_non_empty(group$extraction_rule))),
    800
  )
  row$evidence_basis = .plant_truncate(
    .pubchem_collapse(unique(.uaf_non_empty(group$evidence_basis))),
    800
  )
  row$context_confidence = .plant_best_context_confidence(
    group$context_confidence
  )
  review = row$context_confidence == "low" |
    row$normalized_context %in% c("other", "extract_unspecified", "unknown")
  row$requires_review = .uaf_yes_no(review)
  row
}

.plant_context_basis_rank = function(evidence_basis) {
  basis = .plant_matrix_key(evidence_basis)
  out = rep(0, length(basis))
  out[basis == "reported_context_field"] = 5
  out[basis == "source_backed_species_compound_sentence"] = 4
  out[basis == "source_backed_species_chemical_context_sentence"] = 3
  out[basis == "provider_row_text_regex"] = 2.5
  out[basis == "provider_record_text_regex"] = 2.5
  out[basis == "source_backed_document_context_sentence"] = 2
  out[basis == "source_backed_compound_context_sentence"] = 1
  out
}

.plant_best_context_confidence = function(confidence) {
  confidence = .plant_comparability_clean_confidence(confidence)
  scores = .plant_confidence_score(confidence)
  if (length(scores) < 1 || all(!is.finite(scores))) return("unknown")
  confidence[[which.max(scores)]]
}

.plant_context_evidence_table = function(occurrences, context_type,
                                          raw_context_text,
                                          normalized_context, source_field,
                                          extraction_rule,
                                          context_confidence,
                                          evidence_basis) {
  if (!is.data.frame(occurrences) || nrow(occurrences) < 1L) {
    return(.uaf_empty_table(.plant_context_evidence_cols()))
  }
  n = nrow(occurrences)
  recycle = function(x) rep(x, length.out = n)
  column = function(name, default = NA_character_) {
    .plant_col_or_default(occurrences, name, default)
  }
  species = column("species")
  species_slug = column("species_slug")
  missing_slug = is.na(species_slug) | species_slug == ""
  species_slug[missing_slug] = .plant_slug(species[missing_slug])
  normalized_context = recycle(normalized_context)
  context_confidence = recycle(context_confidence)
  data.frame(
    species = species,
    species_slug = species_slug,
    genus = column("genus"),
    family = column("family"),
    compound_name = column("compound_name"),
    compound_name_clean = column("compound_name_clean"),
    source_database = column("source_database"),
    source_record_id = column("source_record_id"),
    pmid = column("pmid"),
    evidence_url = column("evidence_url"),
    context_type = recycle(context_type),
    raw_context_text = recycle(raw_context_text),
    normalized_context = normalized_context,
    source_field = recycle(source_field),
    extraction_rule = recycle(extraction_rule),
    context_confidence = context_confidence,
    evidence_basis = recycle(evidence_basis),
    requires_review = .uaf_yes_no(
      context_confidence == "low" |
        normalized_context %in% c("other", "extract_unspecified")
    ),
    retrieved_at = rep(.plant_timestamp(), n),
    stringsAsFactors = FALSE
  )
}

.plant_direct_context_confidence_vector = function(occurrences) {
  confidence = .plant_comparability_clean_confidence(
    .plant_col_or_default(occurrences, "confidence", "unknown")
  )
  tier_confidence = .plant_evidence_tier_confidence(
    .plant_col_or_default(occurrences, "evidence_tier", "unresolved")
  )
  out = tier_confidence
  supplied = confidence %in% c("high", "medium")
  out[supplied] = confidence[supplied]
  out[!out %in% c("high", "medium")] = "low"
  out
}

.plant_direct_context_table = function(occurrences) {
  rows = list()
  fields = c("plant_part", "tissue", "method")
  direct_confidence = .plant_direct_context_confidence_vector(occurrences)
  for (field in fields) {
    value = .plant_col_or_default(occurrences, field, NA_character_)
    idx = which(!is.na(value) & value != "")
    if (length(idx) < 1L) next
    if (field == "method") {
      normalized = .plant_method_group(
        value[idx],
        .plant_col_or_default(occurrences, "method_group", "unknown")[idx],
        .plant_col_or_default(occurrences, "source_database", NA_character_)[idx],
        .plant_col_or_default(occurrences, "evidence_tier", "unresolved")[idx]
      )
      context_type = "method"
    } else {
      normalized = .plant_context_group(
        value[idx],
        .plant_col_or_default(
          occurrences, paste0(field, "_group"), "unknown"
        )[idx],
        field
      )
      context_type = field
    }
    keep = !normalized %in%
      c("unknown", "other", "database_record", "literature_curation")
    if (!any(keep)) next
    selected = idx[keep]
    rows[[length(rows) + 1L]] = .plant_context_evidence_table(
      occurrences[selected, , drop = FALSE],
      context_type = context_type,
      raw_context_text = value[selected],
      normalized_context = normalized[keep],
      source_field = field,
      extraction_rule = paste0("direct_column:", field),
      context_confidence = direct_confidence[selected],
      evidence_basis = "reported_context_field"
    )
  }
  .plant_bind_tables(rows, .plant_context_evidence_cols())
}

.plant_context_text_source_table = function(occurrences) {
  n = nrow(occurrences)
  if (n < 1L) {
    return(.uaf_empty_table(c("occurrence_index", "source_field", "text",
                              "rule_prefix", "evidence_basis")))
  }
  source = vapply(
    .plant_col_or_default(occurrences, "source_database", NA_character_),
    .plant_provider_key,
                  character(1))
  evidence_text = .plant_col_or_default(
    occurrences, "evidence_text", NA_character_
  )
  literature = source %in% c("pubtator", "pubmed") &
    !is.na(evidence_text) & evidence_text != ""
  if (any(literature)) {
    literature_idx = which(literature)
    evidence_text[literature_idx] = vapply(literature_idx, function(i) {
      .plant_literature_context_sentences(
        occurrences[i, , drop = FALSE], evidence_text[[i]]
      )
    }, character(1))
  }
  rule_prefix = rep("evidence_text_regex", n)
  evidence_basis = rep("source_text_regex", n)
  rule_prefix[source == "knapsack"] = "knapsack_row_text_regex"
  evidence_basis[source == "knapsack"] = "provider_row_text_regex"
  rule_prefix[source == "lotus"] = "lotus_record_text_regex"
  evidence_basis[source == "lotus"] = "provider_record_text_regex"
  rule_prefix[source == "pubchem_taxonomy"] =
    "pubchem_taxonomy_annotation_regex"
  evidence_basis[source == "pubchem_taxonomy"] =
    "provider_annotation_text_regex"
  literature_source = source %in% c("pubtator", "pubmed")
  rule_prefix[literature_source] = paste0(
    source[literature_source], "_co_mention_sentence_regex"
  )
  evidence_basis[literature_source] = "provider_sentence_regex"
  evidence_idx = which(!is.na(evidence_text) & evidence_text != "")
  evidence = data.frame(
    occurrence_index = evidence_idx,
    source_field = rep("evidence_text", length(evidence_idx)),
    text = evidence_text[evidence_idx],
    rule_prefix = rule_prefix[evidence_idx],
    evidence_basis = evidence_basis[evidence_idx],
    stringsAsFactors = FALSE
  )
  occurrence_type = .plant_col_or_default(
    occurrences, "occurrence_type", NA_character_
  )
  occurrence_idx = which(!is.na(occurrence_type) & occurrence_type != "")
  occurrence = data.frame(
    occurrence_index = occurrence_idx,
    source_field = rep("occurrence_type", length(occurrence_idx)),
    text = occurrence_type[occurrence_idx],
    rule_prefix = rep("occurrence_type_regex", length(occurrence_idx)),
    evidence_basis = rep("occurrence_type_regex", length(occurrence_idx)),
    stringsAsFactors = FALSE
  )
  out = rbind(evidence, occurrence)
  row.names(out) = NULL
  out
}

.plant_text_context_confidence_vector = function(occurrences,
                                                  source_field) {
  base = .plant_direct_context_confidence_vector(occurrences)
  tier = .plant_normalize_evidence_tier(
    .plant_col_or_default(occurrences, "evidence_tier", "unresolved"),
    .plant_col_or_default(occurrences, "source_database", NA_character_),
    .plant_col_or_default(occurrences, "matched_rank", "unresolved")
  )
  out = rep("low", nrow(occurrences))
  occurrence_type = source_field == "occurrence_type"
  out[occurrence_type & base == "high"] = "medium"
  evidence_text = !occurrence_type
  out[evidence_text & tier %in%
        c("manual_curated", "direct_species_database") & base == "high"] =
    "medium"
  out[evidence_text & tier %in%
        c("direct_species_literature", "genus_database_fallback",
          "genus_literature_fallback") & base != "unknown"] = "medium"
  out
}

.plant_prune_context_table = function(x, group) {
  if (!is.data.frame(x) || nrow(x) < 2L) return(x)
  drop = rep(FALSE, nrow(x))
  groups_with = function(type, values) {
    unique(group[x$context_type == type &
                   x$normalized_context %in% values])
  }
  exudate = groups_with("plant_part", "exudate_rhizosphere")
  if (length(exudate) > 0L) {
    drop = drop | (group %in% exudate & x$context_type == "plant_part" &
                     x$normalized_context %in%
                       c("root_belowground", "extract_unspecified"))
  }
  specific_parts = setdiff(
    .plant_context_group_values(),
    c("unknown", "other", "extract_unspecified", "whole_plant")
  )
  specific_part_groups = groups_with("plant_part", specific_parts)
  if (length(specific_part_groups) > 0L) {
    drop = drop | (group %in% specific_part_groups &
                     x$context_type == "plant_part" &
                     x$normalized_context == "extract_unspecified")
  }
  specific_methods = c("gc_ms", "lc_ms", "hplc", "nmr")
  specific_method_groups = groups_with("method", specific_methods)
  if (length(specific_method_groups) > 0L) {
    drop = drop | (group %in% specific_method_groups &
                     x$context_type == "method" &
                     x$normalized_context %in%
                       c("mass_spectrometry", "chromatography",
                         "spectroscopy"))
  }
  mass_only_groups = setdiff(
    groups_with("method", "mass_spectrometry"), specific_method_groups
  )
  if (length(mass_only_groups) > 0L) {
    drop = drop | (group %in% mass_only_groups &
                     x$context_type == "method" &
                     x$normalized_context %in%
                       c("chromatography", "spectroscopy"))
  }
  x[!drop, , drop = FALSE]
}

.plant_text_context_table = function(occurrences) {
  text_sources = .plant_context_text_source_table(occurrences)
  if (nrow(text_sources) < 1L) {
    return(.uaf_empty_table(.plant_context_evidence_cols()))
  }
  patterns = .plant_context_patterns()
  rows = list()
  for (i in seq_len(nrow(patterns))) {
    hit = regexpr(patterns$pattern[[i]], text_sources$text,
                  ignore.case = TRUE, perl = TRUE)
    idx = which(hit > 0L)
    if (length(idx) < 1L) next
    match_length = attr(hit, "match.length")[idx]
    raw = substring(text_sources$text[idx], hit[idx],
                    hit[idx] + match_length - 1L)
    occurrence_idx = text_sources$occurrence_index[idx]
    selected = occurrences[occurrence_idx, , drop = FALSE]
    part = .plant_context_evidence_table(
      selected,
      context_type = patterns$context_type[[i]],
      raw_context_text = .uaf_squish_text(raw),
      normalized_context = patterns$normalized_context[[i]],
      source_field = text_sources$source_field[idx],
      extraction_rule = paste0(text_sources$rule_prefix[idx], ":",
                               patterns$rule_name[[i]]),
      context_confidence = .plant_text_context_confidence_vector(
        selected, text_sources$source_field[idx]
      ),
      evidence_basis = text_sources$evidence_basis[idx]
    )
    part$.text_source_id = idx
    rows[[length(rows) + 1L]] = part
  }
  if (length(rows) < 1L) {
    return(.uaf_empty_table(.plant_context_evidence_cols()))
  }
  out = do.call(rbind, rows)
  out = .plant_prune_context_table(out, out$.text_source_id)
  out$.text_source_id = NULL
  row.names(out) = NULL
  out
}

.plant_direct_context_rows = function(occurrence) {
  rows = list()
  fields = list(
    plant_part = "plant_part",
    tissue = "tissue",
    method = "method"
  )
  for (field in names(fields)) {
    value = .uaf_first_non_empty_text(occurrence[[field]])
    if (is.na(value) || value == "") next
    if (field == "method") {
      normalized = .plant_method_group(
        value, occurrence$method_group, occurrence$source_database,
        occurrence$evidence_tier
      )
      context_type = "method"
    } else {
      normalized = .plant_context_group(
        value, occurrence[[paste0(field, "_group")]], fields[[field]]
      )
      context_type = fields[[field]]
    }
    if (length(normalized) < 1 || normalized[[1]] %in%
        c("unknown", "other", "database_record", "literature_curation")) {
      next
    }
    rows[[length(rows) + 1]] = .plant_context_evidence_row(
      occurrence = occurrence,
      context_type = context_type,
      raw_context_text = value,
      normalized_context = normalized[[1]],
      source_field = field,
      extraction_rule = paste0("direct_column:", field),
      context_confidence = .plant_direct_context_confidence(occurrence),
      evidence_basis = "reported_context_field"
    )
  }
  rows
}

.plant_text_context_rows = function(occurrence) {
  rows = list()
  patterns = .plant_context_patterns()
  text_sources = .plant_context_text_sources(occurrence)
  if (nrow(text_sources) < 1) return(rows)
  for (j in seq_len(nrow(text_sources))) {
    field = text_sources$source_field[[j]]
    text = .uaf_first_non_empty_text(text_sources$text[[j]])
    if (is.na(text) || text == "") next
    field_rows = list()
    for (i in seq_len(nrow(patterns))) {
      hit = .plant_regex_match_text(text, patterns$pattern[[i]])
      if (is.na(hit) || hit == "") next
      field_rows[[length(field_rows) + 1]] = .plant_context_evidence_row(
        occurrence = occurrence,
        context_type = patterns$context_type[[i]],
        raw_context_text = hit,
        normalized_context = patterns$normalized_context[[i]],
        source_field = field,
        extraction_rule = paste0(text_sources$rule_prefix[[j]], ":",
                                 patterns$rule_name[[i]]),
        context_confidence = .plant_text_context_confidence(
          occurrence, patterns$context_type[[i]], field
        ),
        evidence_basis = text_sources$evidence_basis[[j]]
      )
    }
    rows = c(rows, .plant_prune_context_rows(field_rows))
  }
  rows
}

.plant_context_text_sources = function(occurrence) {
  source = .plant_provider_key(occurrence$source_database)
  evidence_text = .uaf_first_non_empty_text(occurrence$evidence_text)
  occurrence_type = .uaf_first_non_empty_text(occurrence$occurrence_type)
  rows = list()
  add = function(field, text, rule_prefix, evidence_basis) {
    text = .uaf_first_non_empty_text(text)
    if (is.na(text) || text == "") return(invisible(FALSE))
    rows[[length(rows) + 1]] <<- data.frame(
      source_field = field,
      text = text,
      rule_prefix = rule_prefix,
      evidence_basis = evidence_basis,
      stringsAsFactors = FALSE
    )
    invisible(TRUE)
  }
  if (source %in% c("pubtator", "pubmed")) {
    relevant = .plant_literature_context_sentences(occurrence, evidence_text)
    add("evidence_text", relevant, paste0(source, "_co_mention_sentence_regex"),
        "provider_sentence_regex")
  } else if (source == "knapsack") {
    add("evidence_text", evidence_text, "knapsack_row_text_regex",
        "provider_row_text_regex")
  } else if (source == "lotus") {
    add("evidence_text", evidence_text, "lotus_record_text_regex",
        "provider_record_text_regex")
  } else if (source == "pubchem_taxonomy") {
    add("evidence_text", evidence_text, "pubchem_taxonomy_annotation_regex",
        "provider_annotation_text_regex")
  } else {
    add("evidence_text", evidence_text, "evidence_text_regex",
        "source_text_regex")
  }
  add("occurrence_type", occurrence_type, "occurrence_type_regex",
      "occurrence_type_regex")
  if (length(rows) < 1) {
    return(.uaf_empty_table(c("source_field", "text", "rule_prefix",
                              "evidence_basis")))
  }
  do.call(rbind, rows)
}

.plant_literature_context_sentences = function(occurrence, text) {
  text = .uaf_first_non_empty_text(text)
  if (is.na(text) || text == "") return(NA_character_)
  sentences = .plant_split_sentences(text)
  if (length(sentences) < 1) return(NA_character_)
  species = .uaf_first_non_empty_text(occurrence$species,
                                      occurrence$query_plant)
  compound = .uaf_first_non_empty_text(occurrence$compound_name)
  keep_both = vapply(sentences, function(sentence) {
    .plant_text_contains(sentence, species) &&
      .plant_text_contains(sentence, compound)
  }, logical(1))
  if (any(keep_both)) {
    return(.plant_truncate(.pubchem_collapse(sentences[keep_both]), 1200))
  }
  keep_compound = vapply(sentences, .plant_text_contains, logical(1),
                         value = compound)
  if (any(keep_compound)) {
    return(.plant_truncate(.pubchem_collapse(sentences[keep_compound]), 1200))
  }
  NA_character_
}

.plant_split_sentences = function(text) {
  text = .uaf_first_non_empty_text(text)
  if (is.na(text) || text == "") return(character())
  sentences = strsplit(text, "(?<=[.!?])\\s+", perl = TRUE)[[1]]
  .uaf_non_empty(.uaf_squish_text(sentences))
}

.plant_provider_key = function(source_database) {
  source = .plant_matrix_key(.uaf_first_non_empty_text(source_database))
  if (source %in% c("pubchem_taxonomy", "pubchem_taxonomy_annotations")) {
    return("pubchem_taxonomy")
  }
  if (grepl("pubchem.*taxonomy|taxonomy.*pubchem", source)) {
    return("pubchem_taxonomy")
  }
  source
}

.plant_prune_context_rows = function(rows) {
  rows = rows[vapply(rows, is.data.frame, logical(1))]
  if (length(rows) < 2) return(rows)
  table = do.call(rbind, rows)
  drop = rep(FALSE, nrow(table))
  has_context = function(type, values) {
    any(table$context_type == type & table$normalized_context %in% values)
  }
  specific_parts = setdiff(.plant_context_group_values(),
                           c("unknown", "other", "extract_unspecified",
                             "whole_plant"))
  if (has_context("plant_part", "exudate_rhizosphere")) {
    drop = drop | (table$context_type == "plant_part" &
                     table$normalized_context %in%
                       c("root_belowground", "extract_unspecified"))
  }
  if (has_context("plant_part", specific_parts)) {
    drop = drop | (table$context_type == "plant_part" &
                     table$normalized_context == "extract_unspecified")
  }
  specific_methods = c("gc_ms", "lc_ms", "hplc", "nmr")
  if (has_context("method", specific_methods)) {
    drop = drop | (table$context_type == "method" &
                     table$normalized_context %in%
                       c("mass_spectrometry", "chromatography",
                         "spectroscopy"))
  } else if (has_context("method", "mass_spectrometry")) {
    drop = drop | (table$context_type == "method" &
                     table$normalized_context %in%
                       c("chromatography", "spectroscopy"))
  }
  table = table[!drop, , drop = FALSE]
  lapply(seq_len(nrow(table)), function(i) table[i, , drop = FALSE])
}

.plant_context_evidence_row = function(occurrence, context_type,
                                       raw_context_text, normalized_context,
                                       source_field, extraction_rule,
                                       context_confidence, evidence_basis) {
  data.frame(
    species = .uaf_first_non_empty_text(occurrence$species),
    species_slug = .uaf_first_non_empty_text(occurrence$species_slug,
                                            .plant_slug(occurrence$species)),
    genus = .uaf_first_non_empty_text(occurrence$genus),
    family = .uaf_first_non_empty_text(occurrence$family),
    compound_name = .uaf_first_non_empty_text(occurrence$compound_name),
    compound_name_clean =
      .uaf_first_non_empty_text(occurrence$compound_name_clean),
    source_database = .uaf_first_non_empty_text(occurrence$source_database),
    source_record_id = .uaf_first_non_empty_text(occurrence$source_record_id),
    pmid = .uaf_first_non_empty_text(occurrence$pmid),
    evidence_url = .uaf_first_non_empty_text(occurrence$evidence_url),
    context_type = context_type,
    raw_context_text = raw_context_text,
    normalized_context = normalized_context,
    source_field = source_field,
    extraction_rule = extraction_rule,
    context_confidence = context_confidence,
    evidence_basis = evidence_basis,
    requires_review = .uaf_yes_no(context_confidence == "low" ||
                                    normalized_context %in%
                                      c("other", "extract_unspecified")),
    retrieved_at = .plant_timestamp(),
    stringsAsFactors = FALSE
  )
}

.plant_context_patterns = function() {
  data.frame(
    context_type = c(
      rep("plant_part", 15),
      rep("tissue", 6),
      rep("method", 10)
    ),
    normalized_context = c(
      "exudate_rhizosphere", "root_belowground", "leaf", "stem_shoot",
      "bark_wood", "flower", "fruit_seed", "aerial", "whole_plant",
      "extract_unspecified", "root_belowground", "fruit_seed",
      "fruit_seed", "leaf", "stem_shoot",
      "vascular", "secretory", "epidermal", "microbial", "whole_plant",
      "root_belowground",
      "gc_ms", "lc_ms", "hplc", "nmr", "mass_spectrometry",
      "chromatography", "spectroscopy", "lc_ms", "gc_ms",
      "chromatography"
    ),
    rule_name = c(
      "exudate_or_rhizosphere", "root_belowground", "leaf_or_foliar",
      "stem_or_shoot", "bark_or_wood", "flower_or_floral",
      "fruit_or_seed", "aerial_tissue", "whole_plant",
      "extract_or_essential_oil", "hairy_root_or_root_culture",
      "cotyledon_or_endosperm", "husk_hull_or_pericarp",
      "leaf_blade_or_foliage", "culm_or_stalk",
      "vascular_tissue", "secretory_tissue", "epidermal_tissue",
      "microbial_context", "callus_or_cell_culture",
      "hairy_root_tissue",
      "gc_ms", "lc_ms", "hplc_uplc", "nmr", "mass_spectrometry",
      "chromatography", "spectroscopy", "tandem_lc_ms",
      "tandem_gc_ms", "headspace_spme"
    ),
    pattern = c(
      "\\b(root exudates?|exudates?|exudation|rhizosphere)\\b",
      "\\b(roots?|rhizomes?|tubers?|bulbs?|below[- ]?ground)\\b",
      "\\b(leaves|leaf|foliar|needles?)\\b",
      "\\b(stems?|shoots?|twigs?|branches?|culms?)\\b",
      "\\b(bark|wood|xylem|phloem)\\b",
      "\\b(flowers?|floral|inflorescences?|petals?|pollen|stamens?)\\b",
      "\\b(fruits?|seeds?|grains?|kernels?|berries|pods?|achenes?|nuts?)\\b",
      "\\b(aerial|above[- ]?ground)\\b",
      "\\b(whole plant|whole-plant|entire plant|whole organism)\\b",
      "\\b(extracts?|essential oils?|volatile oils?)\\b",
      "\\b(hairy roots?|root cultures?)\\b",
      "\\b(cotyledons?|endosperm|embryos?)\\b",
      "\\b(husks?|hulls?|pericarp|peels?|rinds?)\\b",
      "\\b(foliage|leaf blades?|lamina)\\b",
      "\\b(culms?|stalks?)\\b",
      "\\b(vascular|veins?|vasculature)\\b",
      "\\b(glands?|secretory|trichomes?|resin ducts?)\\b",
      "\\b(epidermis|epidermal|cuticle|surface)\\b",
      "\\b(microbial|microbiome|endophytes?|rhizobia|rhizobacterial)\\b",
      "\\b(callus|calli|cell suspension|cell cultures?|in vitro cultures?)\\b",
      "\\b(hairy roots?|root cultures?)\\b",
      "\\b(gc[- ]?ms|gas chromatograph[^.;,]*mass)\\b",
      "\\b(lc[- ]?ms|liquid chromatograph[^.;,]*mass)\\b",
      "\\b(hplc|uplc)\\b",
      "\\b(nmr|nuclear magnetic resonance)\\b",
      "\\b(mass spectrometry|mass spectrometric|\\bms\\b)\\b",
      "\\b(chromatography|chromatographic)\\b",
      "\\b(spectroscopy|spectroscopic|uv[- ]?vis|infrared|ftir)\\b",
      "\\b(lc[- ]?ms/ms|uplc[- ]?ms/ms|uhplc[- ]?ms|uhplc[- ]?qtof|uplc[- ]?qtof)\\b",
      "\\b(gc[- ]?ms/ms|gc[- ]?tof[- ]?ms|gc[- ]?fid|gas chromatography[- ]flame ionization)\\b",
      "\\b(hs[- ]?spme|headspace|solid[- ]phase microextraction)\\b"
    ),
    stringsAsFactors = FALSE
  )
}

.plant_regex_match_text = function(text, pattern) {
  text = .uaf_first_non_empty_text(text)
  if (is.na(text) || text == "") return(NA_character_)
  hit = regexpr(pattern, text, ignore.case = TRUE, perl = TRUE)
  if (hit[[1]] < 0) return(NA_character_)
  .uaf_squish_text(regmatches(text, hit))
}

.plant_direct_context_confidence = function(occurrence) {
  confidence = .plant_comparability_clean_confidence(occurrence$confidence)
  if (confidence %in% c("high", "medium")) return(confidence)
  tier_conf = .plant_evidence_tier_confidence(occurrence$evidence_tier)
  if (tier_conf %in% c("high", "medium")) return(tier_conf)
  "low"
}

.plant_text_context_confidence = function(occurrence, context_type,
                                          source_field) {
  base = .plant_direct_context_confidence(occurrence)
  if (source_field == "occurrence_type") {
    return(ifelse(base == "high", "medium", "low"))
  }
  tier = .plant_normalize_evidence_tier(occurrence$evidence_tier,
                                        occurrence$source_database,
                                        occurrence$matched_rank)
  if (tier %in% c("manual_curated", "direct_species_database") &&
      base == "high") {
    return("medium")
  }
  if (tier %in% c("direct_species_literature",
                  "genus_database_fallback",
                  "genus_literature_fallback") && base != "unknown") {
    return("medium")
  }
  "low"
}

.plant_apply_context_evidence = function(occurrences, context_evidence) {
  occurrences = .plant_normalize_occurrences(occurrences)
  context_evidence = .plant_normalize_context_evidence(context_evidence)
  if (nrow(occurrences) < 1 || nrow(context_evidence) < 1) {
    return(occurrences)
  }
  best = .plant_best_context_evidence(context_evidence)
  if (nrow(best) < 1) return(occurrences)
  occurrence_keys = .plant_context_occurrence_key(occurrences)
  evidence_keys = .plant_context_occurrence_key(best)
  occurrence_relaxed_keys = .plant_context_occurrence_key(
    occurrences, include_source_detail = FALSE
  )
  evidence_relaxed_keys = .plant_context_occurrence_key(
    best, include_source_detail = FALSE
  )
  for (context_type in c("plant_part", "tissue", "method")) {
    evidence_idx = which(best$context_type == context_type)
    if (length(evidence_idx) < 1L) next
    selected = rep(NA_integer_, nrow(occurrences))
    exact_match = match(occurrence_keys, evidence_keys[evidence_idx])
    exact = !is.na(exact_match)
    selected[exact] = evidence_idx[exact_match[exact]]

    relaxed_score = .plant_confidence_score(
      best$context_confidence[evidence_idx]
    ) + .plant_context_specificity_score(
      best$context_type[evidence_idx],
      best$normalized_context[evidence_idx]
    )
    relaxed_order = order(
      evidence_relaxed_keys[evidence_idx], -relaxed_score,
      best$requires_review[evidence_idx],
      best$normalized_context[evidence_idx], evidence_idx
    )
    relaxed_idx = evidence_idx[relaxed_order]
    relaxed_idx = relaxed_idx[
      !duplicated(evidence_relaxed_keys[relaxed_idx])
    ]
    unmatched = which(is.na(selected))
    relaxed_match = match(
      occurrence_relaxed_keys[unmatched],
      evidence_relaxed_keys[relaxed_idx]
    )
    relaxed = !is.na(relaxed_match)
    selected[unmatched[relaxed]] = relaxed_idx[relaxed_match[relaxed]]

    has_evidence = !is.na(selected)
    if (!any(has_evidence)) next
    normalized = rep(NA_character_, nrow(occurrences))
    raw = rep(NA_character_, nrow(occurrences))
    normalized[has_evidence] = best$normalized_context[selected[has_evidence]]
    raw[has_evidence] = .uaf_squish_text(
      best$raw_context_text[selected[has_evidence]]
    )
    missing_raw = has_evidence & (is.na(raw) | raw == "")
    raw[missing_raw] = normalized[missing_raw]
    if (context_type == "plant_part") {
      missing = has_evidence & (occurrences$plant_part_group %in%
        c("unknown", "other", "extract_unspecified") |
        .plant_context_raw_method_like(occurrences$plant_part))
      occurrences$plant_part[missing] = raw[missing]
      occurrences$plant_part_group[missing] = normalized[missing]
    } else if (context_type == "tissue") {
      missing = has_evidence & (occurrences$tissue_group %in%
        c("unknown", "other", "extract_unspecified") |
        .plant_context_raw_method_like(occurrences$tissue))
      occurrences$tissue[missing] = raw[missing]
      occurrences$tissue_group[missing] = normalized[missing]
    } else if (context_type == "method") {
      missing = has_evidence & occurrences$method_group %in%
        c("unknown", "other", "database_record", "literature_curation")
      occurrences$method[missing] = raw[missing]
      occurrences$method_group[missing] = normalized[missing]
    }
  }
  .plant_normalize_occurrences(occurrences)
}

.plant_best_context_evidence = function(context_evidence) {
  context_evidence = .plant_normalize_context_evidence(context_evidence)
  if (nrow(context_evidence) < 1) return(context_evidence)
  key = paste(.plant_context_occurrence_key(context_evidence),
              context_evidence$context_type, sep = "||")
  score = .plant_confidence_score(context_evidence$context_confidence) +
    .plant_context_specificity_score(
      context_evidence$context_type,
      context_evidence$normalized_context
    )
  ordered = order(
    key, -score, context_evidence$requires_review,
    context_evidence$normalized_context, seq_len(nrow(context_evidence))
  )
  keep = ordered[!duplicated(key[ordered])]
  out = context_evidence[keep, , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_context_specificity_score = function(context_type, normalized_context) {
  context_type = .plant_matrix_key(context_type)
  normalized_context = .plant_matrix_key(normalized_context)
  out = rep(0.05, length(normalized_context))
  out[normalized_context %in% c("exudate_rhizosphere", "secretory",
                                "epidermal", "vascular", "microbial",
                                "gc_ms", "lc_ms", "hplc", "nmr")] = 0.20
  out[normalized_context %in% c("root_belowground", "leaf", "stem_shoot",
                                "bark_wood", "flower", "fruit_seed",
                                "aerial")] = 0.15
  out[normalized_context %in% c("whole_plant", "mass_spectrometry",
                                "chromatography", "spectroscopy")] = 0.10
  out[normalized_context %in% c("extract_unspecified",
                                "literature_curation",
                                "database_record", "other", "unknown")] = 0
  out[context_type == "method" & normalized_context %in%
        c("gc_ms", "lc_ms", "hplc", "nmr")] =
    out[context_type == "method" & normalized_context %in%
          c("gc_ms", "lc_ms", "hplc", "nmr")] + 0.05
  out
}

.plant_context_occurrence_key = function(x, include_source_detail = TRUE) {
  cols = c("species", "compound_name_clean", "source_database",
           "source_record_id")
  if (isTRUE(include_source_detail)) {
    cols = c(cols, "pmid", "evidence_url")
  }
  for (col in cols) if (!col %in% names(x)) x[[col]] = NA_character_
  key_data = x[, cols, drop = FALSE]
  key_data[] = lapply(key_data, function(value) {
    value = .uaf_squish_text(value)
    value[is.na(value)] = ""
    value
  })
  do.call(paste, c(key_data, sep = "||"))
}

.plant_context_group = function(value, supplied = NA_character_,
                                kind = c("plant_part", "tissue")) {
  kind = match.arg(kind)
  supplied = rep(.plant_matrix_key(supplied), length.out = length(value))
  valid = .plant_context_group_values()
  supplied[supplied == "unknown"] = NA_character_
  supplied[!supplied %in% valid] = NA_character_
  computed = .plant_first_regex_class(
    value,
    patterns = c(
      "rhizosphere|exudate|exudation",
      "root|rhizome|tuber|bulb|belowground|hairy root",
      "leaf|leaves|foliar|needle|foliage|leaf blade|lamina",
      "stem|shoot|twig|branch|culm|stalk",
      "bark|wood|xylem|phloem",
      "flower|floral|inflorescence|petal|pollen|stamen",
      "seedling|seedlings",
      paste0(
        "fruit|\\bseed\\b|\\bseeds\\b|grain|kernel|berry|pod|achene|",
        "nut|cotyledon|endosperm|embryo|husk|hull|pericarp|peel|rind"
      ),
      "aerial|aboveground|above-ground",
      "whole plant|whole-plant|whole organism|entire plant",
      "vascular|vein|vasculature",
      "gland|secretory|trichome|resin duct",
      "epiderm|cuticle|surface",
      "microb|endophy|rhizob",
      "callus|calli|cell suspension|cell culture|in vitro culture",
      "extract|essential oil|oil"
    ),
    classes = c(
      "exudate_rhizosphere", "root_belowground", "leaf", "stem_shoot",
      "bark_wood", "flower", "whole_plant", "fruit_seed", "aerial",
      "whole_plant", "vascular", "secretory", "epidermal", "microbial",
      "whole_plant", "extract_unspecified"
    )
  )
  out = ifelse(!is.na(supplied) & supplied != "", supplied, computed)
  out[!out %in% valid] = "unknown"
  out
}

.plant_first_regex_class = function(value, patterns, classes,
                                    default = "other") {
  if (length(patterns) != length(classes)) {
    stop("`patterns` and `classes` must have the same length.",
         call. = FALSE)
  }
  text = tolower(.plant_clean_context_value(value))
  computed = rep(default, length(text))
  available = !is.na(text) & text != ""
  computed[!available] = "unknown"
  unmatched = available
  for (i in seq_along(patterns)) {
    hit = unmatched & grepl(patterns[[i]], text)
    hit[is.na(hit)] = FALSE
    computed[hit] = classes[[i]]
    unmatched[hit] = FALSE
  }
  computed
}

.plant_method_group = function(value, supplied = NA_character_,
                               source_database = NA_character_,
                               evidence_tier = NA_character_) {
  n = length(value)
  supplied = rep(.plant_matrix_key(supplied), length.out = n)
  source_database = rep(source_database, length.out = n)
  evidence_tier = rep(evidence_tier, length.out = n)
  valid = .plant_method_group_values()
  supplied[supplied == "unknown"] = NA_character_
  supplied[!supplied %in% valid] = NA_character_
  computed = .plant_first_regex_class(
    value,
    patterns = c(
      "gc[- ]?ms|gc[- ]?tof|gc[- ]?fid|gas chromat",
      "lc[- ]?ms|uplc[- ]?ms|uhplc[- ]?ms|liquid chromat.*mass",
      "hplc|uplc",
      "\\bnmr\\b|nuclear magnetic resonance",
      "mass spectrom|\\bms\\b",
      "chromat",
      "spectroscop|uv-vis|infrared|ftir",
      "literature|paper|publication|curat",
      "database record|provider record|source table"
    ),
    classes = c(
      "gc_ms", "lc_ms", "hplc", "nmr", "mass_spectrometry",
      "chromatography", "spectroscopy", "literature_curation",
      "database_record"
    )
  )
  missing_method = computed == "unknown"
  source = tolower(.uaf_squish_text(source_database))
  tier = .plant_normalize_evidence_tier(evidence_tier, source_database,
                                        "species")
  computed[missing_method & source %in% c("lotus", "knapsack", "npass",
                                          "pubchem taxonomy")] =
    "database_record"
  computed[missing_method & grepl("literature|pubtator", tier)] =
    "literature_curation"
  out = ifelse(!is.na(supplied) & supplied != "", supplied, computed)
  out[!out %in% valid] = "unknown"
  out
}

.plant_context_status = function(plant_part_group, tissue_group, method_group) {
  part_known = !.plant_matrix_key(plant_part_group) %in%
    c("unknown", "extract_unspecified")
  tissue_known = !.plant_matrix_key(tissue_group) %in%
    c("unknown", "extract_unspecified")
  method_known = !.plant_matrix_key(method_group) %in%
    c("unknown", "database_record", "literature_curation")
  ifelse((part_known | tissue_known) & method_known,
         "plant_part_and_method_known",
         ifelse(part_known | tissue_known, "plant_part_known",
                ifelse(method_known, "method_known", "context_missing")))
}

.plant_occurrence_status = function(evidence_tier, matched_rank,
                                    source_database, curation_flag) {
  tier = .plant_normalize_evidence_tier(evidence_tier, source_database,
                                        matched_rank)
  rank = tolower(.uaf_squish_text(matched_rank))
  status = rep("unknown", length(tier))
  status[tier == "manual_curated"] = "curated_reported"
  status[tier == "direct_species_database" & rank == "species"] =
    "direct_reported"
  status[tier == "direct_species_literature"] = "literature_reported"
  status[tier == "direct_species_pubtator_candidate"] = "candidate"
  status[grepl("fallback", tier)] = "taxon_fallback"
  status[tier == "unresolved"] = "unresolved"
  status
}

.plant_occurrence_basis = function(evidence_tier, matched_rank,
                                   source_database) {
  tier = .plant_normalize_evidence_tier(evidence_tier, source_database,
                                        matched_rank)
  rank = tolower(.uaf_squish_text(matched_rank))
  basis = rep("unknown", length(tier))
  basis[tier == "manual_curated"] = "species_curated_record"
  basis[tier == "direct_species_database" & rank == "species"] =
    "species_database_record"
  basis[tier == "direct_species_literature"] = "species_literature_record"
  basis[tier == "direct_species_pubtator_candidate"] = "candidate_co_mention"
  basis[tier == "genus_database_fallback" | rank == "genus"] =
    "genus_fallback_record"
  basis[tier == "family_database_fallback" | rank == "family"] =
    "family_fallback_record"
  basis[tier == "unresolved"] = "unresolved"
  basis
}

.plant_evidence_quality_score = function(occurrence_status, occurrence_basis,
                                         confidence,
                                         biological_context_status) {
  status = .plant_matrix_key(occurrence_status)
  base = rep(0, length(status))
  base[status == "curated_reported"] = 0.95
  base[status == "direct_reported"] = 0.85
  base[status == "literature_reported"] = 0.55
  base[status == "taxon_fallback"] = 0.45
  base[status == "candidate"] = 0.20
  base[status == "unresolved"] = 0
  conf = .plant_confidence_score(confidence)
  context = .plant_matrix_key(biological_context_status)
  bonus = rep(0, length(status))
  bonus[context == "plant_part_and_method_known"] = 0.05
  bonus[context == "plant_part_known"] = 0.03
  bonus[context == "method_known"] = 0.02
  round(pmin(1, base * conf + bonus), 3)
}

.plant_analysis_ready = function(occurrence_status, confidence,
                                 evidence_quality_score) {
  status = .plant_matrix_key(occurrence_status)
  conf_ok = .plant_confidence_score(confidence) >= .plant_confidence_score("medium")
  score = suppressWarnings(as.numeric(evidence_quality_score))
  ifelse(status %in% c("direct_reported", "curated_reported") &
           conf_ok & !is.na(score) & score >= 0.50, "Yes", "No")
}

.plant_matrix_key = function(x) {
  key = .plant_clean_compound(x)
  key[is.na(key) | key == ""] = "unknown"
  key
}

.plant_count_summary = function(x) {
  x = .uaf_non_empty(x)
  if (length(x) < 1) return(NA_character_)
  tab = sort(table(x), decreasing = TRUE)
  paste(paste(names(tab), as.integer(tab), sep = "="), collapse = "; ")
}

.plant_top_terms = function(x, n = 5) {
  x = .uaf_non_empty(x)
  if (length(x) < 1) return(NA_character_)
  tab = sort(table(x), decreasing = TRUE)
  .pubchem_collapse(utils::head(names(tab), n))
}

.plant_source_coverage_score = function(source_count, provider_diagnostics) {
  denom = if (is.data.frame(provider_diagnostics) &&
              nrow(provider_diagnostics) > 0) {
    max(1, sum(provider_diagnostics$enabled == "Yes"))
  } else {
    max(1, source_count)
  }
  round(min(1, source_count / denom), 3)
}

.plant_categorate_validation_status = function(categorate_result) {
  if (is.list(categorate_result) &&
      is.data.frame(categorate_result$ValidationSummary) &&
      "Status" %in% names(categorate_result$ValidationSummary)) {
    return(.uaf_first_non_empty_text(categorate_result$ValidationSummary$Status))
  }
  if (is.null(categorate_result)) return("not_run")
  "not_available"
}

"%||%" = function(x, y) if (is.null(x)) y else x
