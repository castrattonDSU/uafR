# Production provider contracts and local species-metabolite indexes.

.plant_npass_index_version = function() "3"

.plant_npass_index_cols = function() {
  c(
    "index_key_type", "index_key", "species", "species_clean", "genus",
    "family", "org_id", "org_tax_level", "compound_name",
    "compound_name_clean", "np_id", "pubchem_cid", "inchikey", "smiles",
    "inchi", "molecular_formula", "source_record_id", "source_pair_id",
    "new_compound_found",
    "plant_part", "collection_location", "collection_time", "reference_type",
    "reference_id", "reference_id_type", "reference_url", "source_file"
  )
}

.plant_empty_npass_index = function() {
  .uaf_empty_table(.plant_npass_index_cols())
}

.plant_query_alias_cols = function() {
  c("query_id", "query_plant", "species", "alias", "alias_clean",
    "alias_type", "authority", "verified", "source_field")
}

.plant_provider_query_accounting_cols = function() {
  c("query_id", "query_plant", "species", "provider", "query_status",
    "occurrence_count", "literature_candidate_count", "direct_record_count",
    "fallback_record_count", "candidate_record_count", "provider_status",
    "provider_queried", "provider_available", "total_hit_count",
    "query_truncated", "retry_required", "message", "retrieved_at")
}

.plant_provider_resource_manifest_cols = function() {
  c("provider", "resource_type", "resource_id", "resource_path",
    "source_url", "version", "md5", "sha256", "file_size_bytes",
    "availability_status", "redistribution_status", "checked_at", "notes")
}

#' Build a local NPASS species-metabolite index
#'
#' @description
#' Builds a manifest-backed, sharded lookup index from the official NPASS 3.0
#' natural-product general-information, structure, species-source, and species
#' taxonomy TSV files. The large species-source table is streamed in chunks.
#' No records or chemical identifiers are inferred: rows without a source
#' organism or named compound remain excluded from the lookup and are counted
#' in the build summary.
#'
#' @param general_info Path to `NPASS3.0_naturalproducts_generalinfo.txt`.
#' @param structures Path to `NPASS3.0_naturalproducts_structure.txt`.
#' @param species_pairs Path to `NPASS3.0_naturalproducts_species_pair.txt`.
#' @param species_info Path to `NPASS3.0_species_info.txt`.
#' @param out_dir Output directory for sharded RDS files and `manifest.json`.
#' @param overwrite Logical. Replace an existing NPASS index directory.
#' @param chunk_size Number of species-pair lines parsed per streaming chunk.
#' @param progress Logical. Print chunk progress.
#'
#' @return A list with `BuildSummary`, `SourceManifest`, `ShardManifest`, and
#' `IndexPath`.
#'
#' @export
buildNpassIndex = function(general_info,
                           structures,
                           species_pairs,
                           species_info,
                           out_dir,
                           overwrite = FALSE,
                           chunk_size = 100000,
                           progress = interactive()) {
  source_paths = c(
    general_info = general_info,
    structures = structures,
    species_pairs = species_pairs,
    species_info = species_info
  )
  source_paths = vapply(source_paths, .plant_npass_source_path, character(1))
  missing = names(source_paths)[!file.exists(source_paths)]
  if (length(missing) > 0) {
    stop("Missing NPASS source file(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  if ((dir.exists(out_dir) || file.exists(out_dir)) && !isTRUE(overwrite)) {
    stop("NPASS index output exists. Use `overwrite = TRUE`: ", out_dir,
         call. = FALSE)
  }
  chunk_size = suppressWarnings(as.integer(chunk_size[[1]]))
  if (is.na(chunk_size) || chunk_size < 1) {
    stop("`chunk_size` must be a positive integer.", call. = FALSE)
  }

  stage_dir = tempfile("uafR_npass_index_", tmpdir = dirname(out_dir))
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(if (dir.exists(stage_dir)) unlink(stage_dir, recursive = TRUE,
                                            force = TRUE), add = TRUE)
  temp_shards = file.path(stage_dir, "temp_shards")
  final_shards = file.path(stage_dir, "shards")
  dir.create(temp_shards, recursive = TRUE, showWarnings = FALSE)
  dir.create(final_shards, recursive = TRUE, showWarnings = FALSE)

  .plant_progress(progress, "uafR NPASS index: reading compound maps")
  general = .plant_npass_read_tsv(source_paths[["general_info"]])
  structure = .plant_npass_read_tsv(source_paths[["structures"]])
  species = .plant_npass_read_tsv(source_paths[["species_info"]])
  .plant_npass_require_cols(
    general, c("np_id", "pref_name", "iupac_name", "pubchem_id"),
    "general information"
  )
  .plant_npass_require_cols(
    structure, c("np_id", "inchi", "inchikey", "smiles"), "structures"
  )
  .plant_npass_require_cols(
    species, c("org_id", "org_name", "org_tax_level", "genus_name",
               "family_name"), "species information"
  )
  general = .plant_npass_unique_map(general, "np_id")
  structure = .plant_npass_unique_map(structure, "np_id")
  species = .plant_npass_unique_map(species, "org_id")

  source_manifest = .plant_npass_source_manifest(source_paths)
  pair_connection = file(source_paths[["species_pairs"]], open = "r",
                         encoding = "UTF-8")
  pair_connection_open = TRUE
  on.exit(if (isTRUE(pair_connection_open)) close(pair_connection), add = TRUE)
  header = readLines(pair_connection, n = 1L, warn = FALSE)
  if (length(header) != 1L) {
    stop("The NPASS species-pair file is empty.", call. = FALSE)
  }
  pair_columns = .plant_npass_column_names(
    strsplit(header, "\t", fixed = TRUE)[[1]]
  )
  required_pairs = c("src_org_record_id", "src_org_pair_id", "org_id",
                     "np_id", "org_isolation_part", "ref_type", "ref_id",
                     "ref_id_type", "ref_url")
  missing_pairs = setdiff(required_pairs, pair_columns)
  if (length(missing_pairs) > 0) {
    stop("NPASS species-pair file is missing columns: ",
         paste(missing_pairs, collapse = ", "), call. = FALSE)
  }

  chunks = 0L
  pair_rows = 0L
  usable_rows = 0L
  omitted_rows = 0L
  repeat {
    lines = readLines(pair_connection, n = chunk_size, warn = FALSE)
    if (length(lines) < 1) break
    chunks = chunks + 1L
    pair_rows = pair_rows + length(lines)
    .plant_progress(progress, "uafR NPASS index chunk ", chunks, ": ",
                    format(length(lines), big.mark = ","), " pair row(s)")
    pairs = .plant_npass_parse_lines(header, lines)
    normalized = .plant_npass_join_chunk(
      pairs, general, structure, species,
      source_file = basename(source_paths[["species_pairs"]])
    )
    usable_rows = usable_rows + nrow(normalized)
    omitted_rows = omitted_rows + nrow(pairs) - nrow(normalized)
    expanded = .plant_npass_expand_keys(normalized)
    .plant_npass_append_shards(expanded, temp_shards)
  }
  close(pair_connection)
  pair_connection_open = FALSE

  .plant_progress(progress, "uafR NPASS index: finalizing shards")
  shard_manifest = .plant_npass_finalize_shards(temp_shards, final_shards)
  if (nrow(shard_manifest) < 1 || usable_rows < 1) {
    stop("No usable NPASS species-compound rows were indexed.", call. = FALSE)
  }
  build_summary = data.frame(
    provider = "NPASS",
    index_version = .plant_npass_index_version(),
    source_version = "NPASS 3.0 / NPASS-2026",
    pair_row_count = pair_rows,
    usable_occurrence_row_count = usable_rows,
    omitted_pair_row_count = omitted_rows,
    indexed_key_row_count = sum(shard_manifest$row_count),
    shard_count = nrow(shard_manifest),
    chunk_count = chunks,
    built_at = .plant_timestamp(),
    stringsAsFactors = FALSE
  )
  manifest = list(
    Format = "uafR_npass_lookup_index",
    IndexVersion = .plant_npass_index_version(),
    ProviderVersion = "NPASS 3.0 / NPASS-2026",
    TemporaryShardFormat = "chunked_rds",
    BuildSummary = build_summary,
    SourceManifest = source_manifest,
    ShardManifest = shard_manifest
  )
  jsonlite::write_json(
    manifest, file.path(stage_dir, "manifest.json"), dataframe = "rows",
    pretty = TRUE, na = "null", auto_unbox = TRUE
  )
  if (dir.exists(out_dir) || file.exists(out_dir)) {
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }
  if (!file.rename(stage_dir, out_dir)) {
    stop("Could not publish the completed NPASS index: ", out_dir,
         call. = FALSE)
  }
  stage_dir = NA_character_
  out = list(
    BuildSummary = build_summary,
    SourceManifest = source_manifest,
    ShardManifest = shard_manifest,
    IndexPath = out_dir
  )
  class(out) = c("uaf_npass_index_build", "list")
  out
}

#' Query a local NPASS species-metabolite index
#'
#' @param plants Character vector or data frame of plant names.
#' @param npass_index Manifest-backed directory created by
#' `buildNpassIndex()`, or a standardized NPASS index data frame.
#' @param taxon_fallback Taxonomic ranks to query. Species is always queried;
#' genus and family rows remain explicit fallback evidence.
#' @param max_records Maximum records retained per input plant.
#'
#' @return A normalized `PlantCompoundOccurrences` table. Source structure
#' evidence is attached as the `SourceCompoundIdentity` attribute.
#'
#' @export
queryNpassIndex = function(plants, npass_index,
                           taxon_fallback = c("species", "genus"),
                           max_records = Inf) {
  queries = .plant_queries(plants, taxon_fallback)
  aliases = .plant_query_aliases(plants, queries)
  result = .plant_npass_query_result(
    queries, npass_index, max_records, plant_aliases = aliases
  )
  occurrences = result$PlantCompoundOccurrences
  attr(occurrences, "SourceCompoundIdentity") = result$SourceCompoundIdentity
  occurrences
}

#' Inspect plant-provider readiness
#'
#' @description
#' Reports whether each requested species-first provider has the local index or
#' live-service configuration required for a run. This is a configuration and
#' resource check by default; it does not make live requests.
#'
#' @param sources Provider names.
#' @param provider_indexes Named list of local provider indexes. Supported names
#' currently include `lotus` and `npass`.
#' @param lotus_index Backward-compatible LOTUS index argument.
#' @param probe_live Logical. If `TRUE`, perform a lightweight request to live
#' provider landing/API endpoints.
#' @param request_fun Optional request function used when `probe_live = TRUE`.
#' @param request_timeout Maximum live probe time in seconds.
#'
#' @return A provider resource manifest data frame.
#'
#' @export
plantProviderAvailability = function(
    sources = c("lotus", "knapsack", "npass", "pubchem", "pubmed",
                "pubtator"),
    provider_indexes = NULL,
    lotus_index = Sys.getenv("UAFR_LOTUS_INDEX", ""),
    probe_live = FALSE,
    request_fun = NULL,
    request_timeout = 15) {
  sources = .plant_normalize_sources(sources)
  indexes = .plant_provider_indexes(provider_indexes, lotus_index)
  rows = lapply(sources, function(provider) {
    .plant_provider_resource_row(
      provider, indexes[[provider]], probe_live, request_fun, request_timeout
    )
  })
  .plant_bind_tables(rows, .plant_provider_resource_manifest_cols())
}

#' Merge plant phytochemistry discovery results
#'
#' @description
#' Deterministically combines separately completed provider or species batches.
#' Occurrence evidence is preserved, exact duplicate evidence keys are removed,
#' and per-species provider accounting prefers completed record/no-record states
#' over earlier retry-required states.
#'
#' @param ... Plant phytochemistry result objects, or one list containing them.
#' @param strict Logical passed to result validation.
#'
#' @return A merged `uaf_plant_phytochemistry` result.
#'
#' @export
mergePlantPhytochemistryResults = function(..., strict = FALSE) {
  results = list(...)
  if (length(results) == 1L && is.list(results[[1]]) &&
      !inherits(results[[1]], "uaf_plant_phytochemistry") &&
      is.null(results[[1]]$PlantQueries)) {
    results = results[[1]]
  }
  results = results[vapply(results, .plant_valid_discovery_result, logical(1))]
  if (length(results) < 1) {
    stop("Supply at least one valid plant phytochemistry result.", call. = FALSE)
  }
  queries = .plant_merge_query_tables(results)
  merged = .plant_combine_batch_results(results, queries,
                                        unique(queries$taxon_fallback))
  merged$PlantQueryAliases = .plant_merge_alias_tables(results, queries)
  merged$ProviderQueryAccounting = .plant_merge_accounting_tables(
    results, queries
  )
  merged$ProviderResourceManifest = .plant_bind_unique_tables(
    lapply(results, `[[`, "ProviderResourceManifest"),
    .plant_provider_resource_manifest_cols(),
    c("provider", "resource_type", "resource_id", "md5")
  )
  merged$SourceCompoundIdentity = .plant_merge_source_identity_tables(
    lapply(results, `[[`, "SourceCompoundIdentity")
  )
  merged$DataDictionary = plantPhytochemistrySchema()
  merged$Validation = validatePlantPhytochemistryResult(merged,
                                                        strict = strict)
  merged
}

.plant_npass_source_path = function(path) {
  path = .uaf_first_non_empty_text(path)
  if (is.na(path)) return(NA_character_)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.plant_npass_column_names = function(x) {
  x = tolower(.uaf_squish_text(x))
  x = gsub("[^a-z0-9]+", "_", x)
  x = gsub("^_+|_+$", "", x)
  make.unique(x, sep = "_")
}

.plant_npass_read_tsv = function(path) {
  out = utils::read.delim(
    path, stringsAsFactors = FALSE, check.names = FALSE, quote = "",
    comment.char = "", fill = TRUE, na.strings = character(),
    fileEncoding = "UTF-8"
  )
  names(out) = .plant_npass_column_names(names(out))
  out[] = lapply(out, .plant_npass_clean_value)
  out
}

.plant_npass_clean_value = function(x) {
  x = .uaf_squish_text(x)
  missing = is.na(x) | x == "" | tolower(x) %in%
    c("n.a.", "n.a", "na", "n/a", "null", "none", "-")
  x[missing] = NA_character_
  x
}

.plant_npass_require_cols = function(x, required, label) {
  missing = setdiff(required, names(x))
  if (length(missing) > 0) {
    stop("NPASS ", label, " file is missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.plant_npass_unique_map = function(x, key) {
  x = x[!is.na(x[[key]]) & x[[key]] != "", , drop = FALSE]
  x[!duplicated(x[[key]]), , drop = FALSE]
}

.plant_npass_parse_lines = function(header, lines) {
  text = paste(c(header, lines), collapse = "\n")
  out = utils::read.delim(
    text = text, stringsAsFactors = FALSE, check.names = FALSE, quote = "",
    comment.char = "", fill = TRUE, na.strings = character()
  )
  names(out) = .plant_npass_column_names(names(out))
  out[] = lapply(out, .plant_npass_clean_value)
  out
}

.plant_npass_join_chunk = function(pairs, general, structure, species,
                                   source_file) {
  general_idx = match(pairs$np_id, general$np_id)
  structure_idx = match(pairs$np_id, structure$np_id)
  species_idx = match(pairs$org_id, species$org_id)
  pick = function(table, idx, column) {
    if (!column %in% names(table)) return(rep(NA_character_, length(idx)))
    out = rep(NA_character_, length(idx))
    hit = !is.na(idx)
    out[hit] = as.character(table[[column]][idx[hit]])
    .plant_npass_clean_value(out)
  }
  first = function(...) {
    values = list(...)
    out = rep(NA_character_, nrow(pairs))
    for (value in values) {
      value = .plant_npass_clean_value(value)
      fill = (is.na(out) | out == "") & !is.na(value) & value != ""
      out[fill] = value[fill]
    }
    out
  }
  species_name = first(pick(species, species_idx, "species_name"),
                       pick(species, species_idx, "org_name"))
  genus = first(pick(species, species_idx, "genus_name"),
                .plant_genus(species_name))
  family = pick(species, species_idx, "family_name")
  compound_name = first(pick(general, general_idx, "pref_name"),
                        pick(general, general_idx, "iupac_name"))
  inchikey = first(pick(structure, structure_idx, "inchikey"),
                   pick(general, general_idx, "inchikey"))
  out = data.frame(
    species = species_name,
    species_clean = .plant_clean_name(species_name),
    genus = genus,
    family = family,
    org_id = pairs$org_id,
    org_tax_level = pick(species, species_idx, "org_tax_level"),
    compound_name = compound_name,
    compound_name_clean = .plant_clean_compound(compound_name),
    np_id = pairs$np_id,
    pubchem_cid = pick(general, general_idx, "pubchem_id"),
    inchikey = inchikey,
    smiles = pick(structure, structure_idx, "smiles"),
    inchi = pick(structure, structure_idx, "inchi"),
    molecular_formula = first(
      pick(general, general_idx, "molecular_formula"),
      pick(general, general_idx, "formula")
    ),
    source_record_id = pairs$src_org_record_id,
    source_pair_id = pairs$src_org_pair_id,
    new_compound_found = .plant_npass_pair_col(pairs, "new_cp_found"),
    plant_part = .plant_npass_pair_col(pairs, "org_isolation_part"),
    collection_location = .plant_npass_pair_col(pairs,
                                                "org_collect_location"),
    collection_time = .plant_npass_pair_col(pairs, "org_collect_time"),
    reference_type = .plant_npass_pair_col(pairs, "ref_type"),
    reference_id = .plant_npass_pair_col(pairs, "ref_id"),
    reference_id_type = .plant_npass_pair_col(pairs, "ref_id_type"),
    reference_url = .plant_npass_pair_col(pairs, "ref_url"),
    source_file = source_file,
    stringsAsFactors = FALSE
  )
  usable = !is.na(out$species_clean) & out$species_clean != "" &
    !is.na(out$compound_name_clean) & out$compound_name_clean != ""
  out = out[usable, , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_npass_pair_col = function(x, column) {
  if (!column %in% names(x)) return(rep(NA_character_, nrow(x)))
  .plant_npass_clean_value(x[[column]])
}

.plant_npass_expand_keys = function(x) {
  if (!is.data.frame(x) || nrow(x) < 1) return(.plant_empty_npass_index())
  pieces = list()
  add = function(type, key) {
    keep = !is.na(key) & key != ""
    if (!any(keep)) return(NULL)
    rows = x[keep, , drop = FALSE]
    rows$index_key_type = type
    rows$index_key = .plant_clean_name(key[keep])
    rows[, .plant_npass_index_cols(), drop = FALSE]
  }
  pieces[[1]] = add("species", x$species_clean)
  pieces[[2]] = add("genus", x$genus)
  pieces[[3]] = add("family", x$family)
  .plant_bind_tables(pieces, .plant_npass_index_cols())
}

.plant_npass_shard_prefix = function(x) {
  key = .plant_clean_compound(x)
  key[is.na(key) | key == ""] = "_"
  prefix = substring(key, 1L, 2L)
  gsub("[^a-z0-9_]", "_", prefix)
}

.plant_npass_append_shards = function(x, root) {
  if (!is.data.frame(x) || nrow(x) < 1) return(invisible(NULL))
  group = paste(x$index_key_type, .plant_npass_shard_prefix(x$index_key),
                sep = "/")
  groups = split(seq_len(nrow(x)), group)
  for (name in names(groups)) {
    directory = file.path(root, name)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    file = tempfile("chunk_", tmpdir = directory, fileext = ".rds")
    saveRDS(
      x[groups[[name]], .plant_npass_index_cols(), drop = FALSE],
      file, compress = FALSE, version = 3
    )
  }
  invisible(NULL)
}

.plant_npass_finalize_shards = function(temp_root, final_root) {
  files = list.files(temp_root, pattern = "[.]rds$", recursive = TRUE,
                     full.names = TRUE)
  relative = substring(files, nchar(temp_root) + 2L)
  relative = gsub("\\\\", "/", relative)
  shard_group = sub("/[^/]+$", "", relative)
  grouped_files = split(files, shard_group)
  rows = list()
  for (i in seq_along(grouped_files)) {
    table = .plant_bind_tables(
      lapply(grouped_files[[i]], readRDS), .plant_npass_index_cols()
    )
    key = do.call(paste, c(table[c("index_key_type", "index_key",
                                  "source_record_id", "np_id")], sep = "\r"))
    table = table[!duplicated(key), , drop = FALSE]
    table = table[order(table$index_key, table$species,
                        table$compound_name_clean, table$source_record_id),
                  , drop = FALSE]
    relative_out = paste0(names(grouped_files)[[i]], ".rds")
    out_file = file.path(final_root, relative_out)
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(table, out_file, compress = "gzip")
    rows[[i]] = data.frame(
      index_key_type = table$index_key_type[[1]],
      prefix = sub("[.]rds$", "", basename(out_file)),
      file = file.path("shards", relative_out),
      row_count = nrow(table),
      key_count = length(unique(table$index_key)),
      md5 = unname(tools::md5sum(out_file)[[1]]),
      file_size_bytes = file.info(out_file)$size,
      stringsAsFactors = FALSE
    )
  }
  .plant_bind_tables(rows, c("index_key_type", "prefix", "file",
                             "row_count", "key_count", "md5",
                             "file_size_bytes"))
}

.plant_npass_source_manifest = function(paths) {
  info = file.info(paths)
  data.frame(
    source_role = names(paths),
    source_file = basename(paths),
    source_path = unname(paths),
    source_url = paste0(
      "https://bidd.group/NPASS/downloadFiles/", basename(paths)
    ),
    file_size_bytes = as.numeric(info$size),
    md5 = unname(tools::md5sum(paths)),
    sha256 = vapply(paths, .plant_sha256_file, character(1)),
    checked_at = .plant_timestamp(),
    stringsAsFactors = FALSE
  )
}

.plant_sha256_file = function(path) {
  executable = Sys.which(c("sha256sum", "shasum"))
  if (nzchar(executable[["sha256sum"]])) {
    out = suppressWarnings(system2(executable[["sha256sum"]], path,
                                   stdout = TRUE, stderr = TRUE))
  } else if (nzchar(executable[["shasum"]])) {
    out = suppressWarnings(system2(executable[["shasum"]], c("-a", "256", path),
                                   stdout = TRUE, stderr = TRUE))
  } else {
    return(NA_character_)
  }
  value = strsplit(.uaf_first_non_empty_text(out), "\\s+")[[1]][1]
  if (is.na(value) || !grepl("^[a-fA-F0-9]{64}$", value)) NA_character_ else
    tolower(value)
}

.plant_npass_lookup_info = function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !dir.exists(path)) return(NULL)
  manifest_file = file.path(path, "manifest.json")
  if (!file.exists(manifest_file)) return(NULL)
  manifest = tryCatch(jsonlite::read_json(manifest_file, simplifyVector = TRUE),
                      error = function(error) NULL)
  if (!is.list(manifest) ||
      !identical(manifest$Format, "uafR_npass_lookup_index")) return(NULL)
  list(root = normalizePath(path, winslash = "/", mustWork = TRUE),
       manifest_file = manifest_file, manifest = manifest)
}

.plant_npass_index_or_null = function(x) {
  if (is.data.frame(x)) return(x)
  if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) return(x)
  NULL
}

.plant_npass_query_result = function(plant_queries, npass_index, max_records,
                                     plant_aliases = NULL) {
  npass_index = .plant_npass_index_or_null(npass_index)
  if (is.null(npass_index)) {
    stop("A local NPASS index is required. Build one with `buildNpassIndex()`.",
         call. = FALSE)
  }
  if (is.data.frame(npass_index)) {
    index = .plant_bind_tables(list(npass_index), .plant_npass_index_cols())
    lookup = NULL
  } else {
    lookup = .plant_npass_lookup_info(npass_index)
    if (is.null(lookup)) {
      stop("`npass_index` is not a valid manifest-backed NPASS index: ",
           npass_index, call. = FALSE)
    }
    index = NULL
  }
  occurrence_rows = list()
  identity_rows = list()
  for (i in seq_len(nrow(plant_queries))) {
    query = plant_queries[i, , drop = FALSE]
    keys = .plant_lotus_lookup_query_keys(query, plant_aliases)
    matches = if (is.null(lookup)) {
      .plant_npass_index_matches(index, query, plant_aliases)
    } else {
      .plant_npass_lookup_matches(lookup, keys, query, plant_aliases)
    }
    if (nrow(matches) < 1) next
    rank_order = match(matches$index_key_type,
                       c("species", "genus", "family"))
    matches = matches[order(rank_order, matches$compound_name_clean,
                            matches$source_record_id), , drop = FALSE]
    evidence_key = paste(matches$source_record_id, matches$np_id,
                         matches$species_clean, sep = "\r")
    matches = matches[!duplicated(evidence_key), , drop = FALSE]
    limit = .plant_max_records(max_records)
    if (is.finite(limit) && nrow(matches) > limit) {
      matches = matches[seq_len(limit), , drop = FALSE]
    }
    occurrence_rows[[length(occurrence_rows) + 1L]] =
      .plant_npass_occurrences(query, matches)
    identity_rows[[length(identity_rows) + 1L]] =
      .plant_npass_source_identity(matches)
  }
  list(
    PlantCompoundOccurrences = .plant_bind_occurrences(occurrence_rows),
    SourceCompoundIdentity = .plant_bind_unique_tables(
      identity_rows, .plant_source_compound_identity_cols(),
      c("source_database", "source_record_id", "source_compound_id",
        "InChIKey", "SMILES")
    )
  )
}

.plant_npass_lookup_matches = function(lookup, keys, query,
                                       plant_aliases = NULL) {
  species_keys = .plant_clean_name(
    .plant_verified_species_aliases(plant_aliases, query)
  )
  rows = list()
  for (i in seq_len(nrow(keys))) {
    type = keys$index_key_type[[i]]
    key = keys$index_key[[i]]
    file = file.path(lookup$root, "shards", type,
                     paste0(.plant_npass_shard_prefix(key), ".rds"))
    if (!file.exists(file)) next
    table = tryCatch(readRDS(file), error = function(error) NULL)
    if (!is.data.frame(table)) next
    table = .plant_bind_tables(list(table), .plant_npass_index_cols())
    hit = table$index_key_type == type & table$index_key == key
    if (type != "species") {
      hit = hit & !table$species_clean %in% species_keys
    }
    if (any(hit, na.rm = TRUE)) rows[[length(rows) + 1L]] = table[hit, , drop = FALSE]
  }
  .plant_bind_tables(rows, .plant_npass_index_cols())
}

.plant_npass_index_matches = function(index, query, plant_aliases = NULL) {
  keys = .plant_lotus_lookup_query_keys(query, plant_aliases)
  species_keys = .plant_clean_name(
    .plant_verified_species_aliases(plant_aliases, query)
  )
  rows = list()
  for (i in seq_len(nrow(keys))) {
    type = keys$index_key_type[[i]]
    key = keys$index_key[[i]]
    hit = index$index_key_type == type & index$index_key == key
    if (type != "species") {
      hit = hit & !index$species_clean %in% species_keys
    }
    if (any(hit, na.rm = TRUE)) rows[[length(rows) + 1L]] = index[hit, , drop = FALSE]
  }
  .plant_bind_tables(rows, .plant_npass_index_cols())
}

.plant_npass_occurrences = function(query, matches) {
  rows = lapply(seq_len(nrow(matches)), function(i) {
    hit = matches[i, , drop = FALSE]
    rank = hit$index_key_type[[1]]
    pmid = if (.plant_npass_id_type(hit$reference_id_type) == "pmid") {
      hit$reference_id[[1]]
    } else NA_character_
    doi = if (.plant_npass_id_type(hit$reference_id_type) == "doi") {
      hit$reference_id[[1]]
    } else NA_character_
    data.frame(
      query_plant = query$query_plant,
      query_plant_clean = query$query_plant_clean,
      matched_taxon = hit$index_key,
      matched_rank = rank,
      species = query$species,
      genus = query$genus,
      family = .uaf_first_non_empty_text(query$family, hit$family),
      compound_name = hit$compound_name,
      compound_name_clean = hit$compound_name_clean,
      compound_id = hit$np_id,
      compound_id_type = "NPASS_NP_ID",
      source_database = "NPASS",
      source_record_id = hit$source_record_id,
      evidence_text = .plant_npass_evidence_text(hit),
      evidence_url = .plant_npass_evidence_url(hit),
      reference_id = hit$reference_id,
      pmid = pmid,
      doi = doi,
      plant_part = hit$plant_part,
      tissue = NA_character_,
      method = NA_character_,
      occurrence_type = "npass_species_source_record",
      retrieved_at = .plant_timestamp(),
      confidence = ifelse(rank == "species", "high", "medium"),
      curation_flag = ifelse(rank == "species", "source_database_record",
                             "taxon_fallback_review_required"),
      evidence_tier = ifelse(rank == "species", "direct_species_database",
                             paste0(rank, "_database_fallback")),
      stringsAsFactors = FALSE
    )
  })
  .plant_bind_occurrences(rows)
}

.plant_npass_source_identity = function(matches) {
  rows = lapply(seq_len(nrow(matches)), function(i) {
    hit = matches[i, , drop = FALSE]
    structure_available = length(.uaf_non_empty(c(hit$inchikey, hit$smiles))) > 0
    data.frame(
      compound_name = hit$compound_name,
      compound_name_clean = hit$compound_name_clean,
      source_database = "NPASS",
      source_record_id = hit$source_record_id,
      source_compound_id = hit$np_id,
      source_compound_id_type = "NPASS_NP_ID",
      CID = suppressWarnings(as.integer(hit$pubchem_cid)),
      InChIKey = hit$inchikey,
      SMILES = hit$smiles,
      MolecularFormula = hit$molecular_formula,
      evidence_url = .plant_npass_evidence_url(hit),
      evidence_text = .plant_npass_evidence_text(hit),
      identity_status = ifelse(structure_available,
                               "source_structure_available",
                               "source_identifier_only"),
      identity_note = paste(
        "NPASS source identifiers and structures are retained as published;",
        "cross-source identity agreement still requires validation."
      ),
      stringsAsFactors = FALSE
    )
  })
  .plant_bind_tables(rows, .plant_source_compound_identity_cols())
}

.plant_npass_id_type = function(x) {
  value = tolower(.uaf_first_non_empty_text(x, ""))
  if (grepl("pmid|pubmed", value)) return("pmid")
  if (grepl("doi", value)) return("doi")
  value
}

.plant_npass_evidence_url = function(hit) {
  .uaf_first_non_empty_text(
    hit$reference_url,
    paste0("https://bidd.group/NPASS/compound.php?compoundID=", hit$np_id)
  )
}

.plant_npass_evidence_text = function(hit) {
  value = function(column) .uaf_first_non_empty_text(hit[[column]])
  values = c(
    paste0("NPASS species source: ", value("species")),
    if (!is.na(value("plant_part"))) paste0("isolation part: ",
                                             value("plant_part")),
    if (!is.na(value("collection_location"))) paste0(
      "collection location: ", value("collection_location")
    ),
    if (!is.na(value("collection_time"))) paste0("collection time: ",
                                                 value("collection_time")),
    if (!is.na(value("reference_type"))) paste0("reference type: ",
                                               value("reference_type")),
    if (!is.na(value("reference_id"))) paste0("reference: ",
                                             value("reference_id"))
  )
  .plant_truncate(paste(.uaf_non_empty(values), collapse = "; "), 800)
}

.plant_provider_indexes = function(provider_indexes, lotus_index) {
  indexes = if (is.list(provider_indexes)) provider_indexes else list()
  if (length(indexes) > 0) {
    index_names = names(indexes)
    if (is.null(index_names) || any(is.na(index_names) | index_names == "")) {
      stop("`provider_indexes` must be a named list.", call. = FALSE)
    }
    names(indexes) = tolower(index_names)
  }
  lotus = .plant_lotus_index_or_null(lotus_index)
  if (is.null(indexes$lotus) && !is.null(lotus)) indexes$lotus = lotus
  indexes
}

.plant_provider_indexes_signature = function(indexes, sources) {
  sources = sort(unique(.plant_normalize_sources(sources)))
  pieces = vapply(sources, function(provider) {
    index = indexes[[provider]]
    if (provider == "lotus") return(.plant_lotus_index_signature(index))
    if (provider == "npass") {
      info = .plant_npass_lookup_info(index)
      if (!is.null(info)) {
        stat = file.info(info$manifest_file)
        return(paste("npass_index", normalizePath(info$manifest_file,
                                                   winslash = "/"),
                     stat$size, stat$mtime,
                     unname(tools::md5sum(info$manifest_file)[[1]]), sep = ":"))
      }
      if (is.data.frame(index)) {
        return(paste("npass_index:data_frame", nrow(index), ncol(index),
                     paste(names(index), collapse = ","), sep = ":"))
      }
      return("npass_index:none")
    }
    paste0(provider, ":live_service")
  }, character(1))
  paste0("provider_indexes:", .pubchem_url_hash(paste(pieces,
                                                       collapse = "\r")))
}

.plant_provider_resource_row = function(provider, index, probe_live,
                                         request_fun, request_timeout) {
  urls = c(
    lotus = "https://lotus.naturalproducts.net/",
    knapsack = "https://www.knapsackfamily.com/knapsack_core/top.php",
    npass = "https://bidd.group/NPASS/downloadnpass.html",
    pubchem = "https://pubchem.ncbi.nlm.nih.gov/",
    pubmed = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
    pubtator = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/"
  )
  probe_urls = urls
  probe_urls[["pubmed"]] = paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi",
    "?db=pubmed&retmode=json"
  )
  local_required = provider %in% c("lotus", "npass")
  resource_path = if (is.character(index) && length(index) == 1L) index else
    NA_character_
  available = if (provider == "lotus") {
    !is.null(.plant_lotus_index_or_null(index))
  } else if (provider == "npass") {
    is.data.frame(index) || !is.null(.plant_npass_lookup_info(index))
  } else TRUE
  version = NA_character_
  md5 = NA_character_
  sha256 = NA_character_
  checksum_path = NA_character_
  if (provider == "npass" && !is.null(.plant_npass_lookup_info(index))) {
    info = .plant_npass_lookup_info(index)
    version = .uaf_first_non_empty_text(info$manifest$ProviderVersion)
    checksum_path = info$manifest_file
  } else if (provider == "lotus" && is.character(index) && file.exists(index)) {
    lookup = .plant_lotus_lookup_info(index)
    manifest = if (is.null(lookup)) index else lookup$manifest_file
    checksum_path = manifest
    version = if (is.null(lookup)) {
      "local LOTUS index"
    } else {
      .uaf_first_non_empty_text(
        lookup$manifest$source_version,
        lookup$manifest$provider_version,
        lookup$manifest$release,
        "manifest-backed local LOTUS index"
      )
    }
  }
  if (!is.na(checksum_path) && file.exists(checksum_path)) {
    md5 = unname(tools::md5sum(checksum_path)[[1]])
    sha256 = .plant_sha256_file(checksum_path)
  }
  probe_note = NA_character_
  if (isTRUE(probe_live) && !local_required) {
    probe = tryCatch({
      if (is.function(request_fun)) request_fun(probe_urls[[provider]]) else
        .plant_fetch_text(
          probe_urls[[provider]], cache = FALSE, cache_dir = tempdir(),
          throttle = 0, request_fun = NULL, timeout = request_timeout
        )
      TRUE
    }, error = function(error) {
      probe_note <<- .plant_redact_secrets(conditionMessage(error))
      FALSE
    })
    available = isTRUE(probe)
  }
  data.frame(
    provider = provider,
    resource_type = ifelse(local_required, "local_index", "live_service"),
    resource_id = ifelse(local_required, paste0(provider, "_index"),
                         paste0(provider, "_service")),
    resource_path = resource_path,
    source_url = urls[[provider]],
    version = version,
    md5 = md5,
    sha256 = sha256,
    file_size_bytes = if (!is.na(checksum_path) && file.exists(checksum_path)) {
      as.numeric(file.info(checksum_path)$size)
    } else if (!is.na(resource_path) && file.exists(resource_path) &&
               !dir.exists(resource_path)) {
      as.numeric(file.info(resource_path)$size)
    } else NA_real_,
    availability_status = ifelse(available, "available", "unavailable"),
    redistribution_status = ifelse(provider == "knapsack",
                                   "raw_response_not_for_redistribution",
                                   "review_provider_terms_before_redistribution"),
    checked_at = .plant_timestamp(),
    notes = .uaf_first_non_empty_text(
      probe_note,
      ifelse(local_required && !available,
             paste0("Configure a valid local ", provider, " index."),
             "Resource configuration is ready; no biological hit is implied.")
    ),
    stringsAsFactors = FALSE
  )
}

.plant_query_aliases = function(plants, plant_queries) {
  rows = list()
  add_alias = function(query, alias, type, source_field,
                       authority = NA_character_, verified = FALSE) {
    alias = .uaf_first_non_empty_text(alias)
    if (is.na(alias)) return(NULL)
    data.frame(
      query_id = query$query_id,
      query_plant = query$query_plant,
      species = query$species,
      alias = alias,
      alias_clean = .plant_clean_name(alias),
      alias_type = type,
      authority = authority,
      verified = .uaf_yes_no(verified),
      source_field = source_field,
      stringsAsFactors = FALSE
    )
  }
  input = if (is.data.frame(plants)) plants else NULL
  if (!is.null(input)) names(input) = .plant_normalize_column_names(names(input))
  input_species = if (!is.null(input)) {
    field = if ("species" %in% names(input)) "species" else names(input)[[1]]
    .plant_canonical_taxon_name(input[[field]])
  } else character()
  input_raw_clean = if (!is.null(input)) {
    field = if ("species" %in% names(input)) "species" else names(input)[[1]]
    .plant_clean_name(input[[field]])
  } else character()
  for (i in seq_len(nrow(plant_queries))) {
    query = plant_queries[i, , drop = FALSE]
    canonical_key = .plant_clean_name(query$species)
    rows[[length(rows) + 1L]] = add_alias(query, query$query_plant,
                                         "submitted_name", "query_plant",
                                         verified = identical(
                                           .plant_clean_name(query$query_plant),
                                           canonical_key
                                         ))
    rows[[length(rows) + 1L]] = add_alias(query, query$species,
                                         "normalized_species", "species",
                                         verified = TRUE)
    input_i = if (!is.null(input)) {
      which(input_species == query$species |
              input_raw_clean == query$query_plant_clean)[1]
    } else NA_integer_
    if (!is.na(input_i)) {
      for (field in intersect(c("original_species", "canonical_input_species",
                                "accepted_species_name", "wfo_full_name",
                                "verified_synonym"), names(input))) {
        type = ifelse(field == "accepted_species_name", "accepted_name",
                      ifelse(field == "wfo_full_name",
                             "accepted_name_with_authority",
                             ifelse(field == "verified_synonym",
                                    "verified_synonym", "submitted_alias")))
        authority = if (field == "wfo_full_name" && "wfo_id" %in% names(input)) {
          input$wfo_id[[input_i]]
        } else NA_character_
        field_value = input[[field]][[input_i]]
        wfo_verified = field == "wfo_full_name" &&
          length(.uaf_non_empty(authority)) > 0 &&
          (!"wfo_match_status" %in% names(input) ||
             grepl("unambiguous|accepted|exact",
                   .uaf_first_non_empty_text(input$wfo_match_status[[input_i]],
                                             ""),
                   ignore.case = TRUE))
        is_verified = field %in% c("accepted_species_name",
                                   "verified_synonym") || wfo_verified
        if (field %in% c("original_species", "canonical_input_species")) {
          is_verified = identical(.plant_clean_name(field_value),
                                  canonical_key)
        }
        values = .plant_split_taxon_aliases(field_value)
        for (value in values) {
          rows[[length(rows) + 1L]] = add_alias(
            query, value, type, field, authority, is_verified
          )
          if (field == "wfo_full_name" && isTRUE(is_verified)) {
            core = .plant_taxon_core_name(value)
            if (length(.uaf_non_empty(core)) > 0 &&
                .plant_clean_name(core) != .plant_clean_name(value)) {
              rows[[length(rows) + 1L]] = add_alias(
                query, core, "accepted_name", paste0(field, "_taxon_core"),
                authority, TRUE
              )
            }
          }
        }
      }
    }
  }
  out = .plant_bind_tables(rows, .plant_query_alias_cols())
  if (nrow(out) < 1) return(out)
  key = paste(out$query_id, out$alias_clean, sep = "\r")
  out = out[!duplicated(key), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_split_taxon_aliases = function(x) {
  x = .uaf_first_non_empty_text(x)
  if (is.na(x)) return(character())
  parts = unlist(strsplit(x, "[;|]", perl = TRUE), use.names = FALSE)
  unique(.uaf_non_empty(.uaf_squish_text(parts)))
}

.plant_taxon_core_name = function(x) {
  x = .uaf_squish_text(x)
  vapply(x, function(value) {
    if (is.na(value) || !nzchar(value)) return(NA_character_)
    value = gsub(intToUtf8(0x00d7), " x ", value, fixed = TRUE)
    parts = .uaf_non_empty(strsplit(.uaf_squish_text(value), "\\s+")[[1]])
    if (length(parts) < 2 ||
        !grepl("^[A-Za-z][A-Za-z-]+$", parts[[1]])) return(NA_character_)
    genus = .plant_title_word(parts[[1]])
    position = 2L
    marker = tolower(gsub("[.]$", "", parts[[position]]))
    if (marker == "x") {
      if (length(parts) < 3L) return(NA_character_)
      epithet = tolower(gsub("[^A-Za-z-]", "", parts[[3]]))
      if (!grepl("^[a-z][a-z-]+$", epithet)) return(NA_character_)
      core = c(genus, "x", epithet)
      position = 4L
    } else {
      epithet = tolower(gsub("[^A-Za-z-]", "", parts[[position]]))
      if (!grepl("^[a-z][a-z-]+$", epithet) ||
          epithet %in% c("sp", "spp", "cf", "aff")) return(NA_character_)
      core = c(genus, epithet)
      position = 3L
    }
    if (length(parts) >= position + 1L) {
      rank = tolower(gsub("[.]$", "", parts[[position]]))
      ranks = c("subsp", "ssp", "subspecies", "var", "variety",
                "f", "forma")
      if (rank %in% ranks) {
        infra = tolower(gsub("[^A-Za-z-]", "", parts[[position + 1L]]))
        if (grepl("^[a-z][a-z-]+$", infra)) {
          rank = if (rank %in% c("ssp", "subspecies")) "subsp" else
            if (rank == "variety") "var" else if (rank == "forma") "f" else
              rank
          core = c(core, rank, infra)
        }
      }
    }
    paste(core, collapse = " ")
  }, character(1), USE.NAMES = FALSE)
}

.plant_verified_species_aliases = function(plant_aliases, query_row) {
  canonical = .uaf_first_non_empty_text(query_row$species,
                                        query_row$query_plant)
  values = canonical
  if (is.data.frame(plant_aliases) && nrow(plant_aliases) > 0) {
    id = .uaf_first_non_empty_text(query_row$query_id)
    rows = if (!is.na(id) && "query_id" %in% names(plant_aliases)) {
      plant_aliases[plant_aliases$query_id == id, , drop = FALSE]
    } else {
      plant_aliases[
        .plant_clean_name(plant_aliases$species) ==
          .plant_clean_name(query_row$species), , drop = FALSE
      ]
    }
    allowed = c("normalized_species", "accepted_name", "verified_synonym")
    keep = tolower(.uaf_squish_text(rows$verified)) %in% c("yes", "true", "1") &
      rows$alias_type %in% allowed
    values = c(values, rows$alias[keep])
  }
  values = .plant_taxon_core_name(values)
  unique(.uaf_non_empty(values))
}

.plant_provider_query_accounting = function(plant_queries, provider,
                                             occurrences, literature,
                                             diagnostics) {
  occurrences = .plant_normalize_occurrences(occurrences)
  literature = .plant_normalize_literature(literature, source_hint = provider)
  diagnostic = diagnostics[diagnostics$provider == provider, , drop = FALSE]
  if (nrow(diagnostic) < 1) {
    diagnostic = .plant_provider_diagnostics(
      provider, TRUE, FALSE, FALSE, 0, 0, 0, 1, 0, "missing_diagnostic",
      "Provider returned no diagnostic row."
    )
  }
  rows = lapply(seq_len(nrow(plant_queries)), function(i) {
    query = plant_queries[i, , drop = FALSE]
    occ = occurrences$species == query$species |
      occurrences$query_plant_clean == query$query_plant_clean
    lit = literature$species == query$species |
      literature$query_plant_clean == query$query_plant_clean
    occ[is.na(occ)] = FALSE
    lit[is.na(lit)] = FALSE
    occurrence_count = sum(occ)
    literature_count = sum(lit)
    record_count = max(occurrence_count, literature_count)
    provider_status = tolower(.uaf_first_non_empty_text(diagnostic$status,
                                                        "unknown"))
    retry = record_count < 1 && provider_status %in%
      c("warning", "error", "failed", "timeout", "timed_out",
        "rate_limited", "service_unavailable", "not_queried",
        "not_implemented", "unavailable", "missing_diagnostic")
    status = if (record_count > 0) "records" else if (retry) {
      "retry_required"
    } else "no_records"
    data.frame(
      query_id = query$query_id,
      query_plant = query$query_plant,
      species = query$species,
      provider = provider,
      query_status = status,
      occurrence_count = occurrence_count,
      literature_candidate_count = literature_count,
      direct_record_count = sum(
        occ & occurrences$matched_rank == "species" &
          occurrences$occurrence_status %in% c("direct_reported",
                                                "curated_reported"),
        na.rm = TRUE
      ),
      fallback_record_count = sum(occ & occurrences$matched_rank %in%
                                    c("genus", "family"), na.rm = TRUE),
      candidate_record_count = max(
        sum(occ & occurrences$occurrence_status == "candidate", na.rm = TRUE),
        literature_count
      ),
      provider_status = provider_status,
      provider_queried = .uaf_first_non_empty_text(diagnostic$queried, "No"),
      provider_available = .uaf_first_non_empty_text(diagnostic$available,
                                                     "No"),
      total_hit_count = .plant_accounting_total_hits(literature[lit, ,
                                                                  drop = FALSE],
                                                     record_count),
      query_truncated = .plant_accounting_truncated(literature[lit, ,
                                                                 drop = FALSE]),
      retry_required = .uaf_yes_no(retry),
      message = .uaf_first_non_empty_text(diagnostic$message),
      retrieved_at = .uaf_first_non_empty_text(diagnostic$retrieved_at,
                                               .plant_timestamp()),
      stringsAsFactors = FALSE
    )
  })
  .plant_bind_tables(rows, .plant_provider_query_accounting_cols())
}

.plant_accounting_total_hits = function(literature, fallback) {
  if (is.data.frame(literature) &&
      "literature_total_hit_count" %in% names(literature)) {
    values = suppressWarnings(as.integer(literature$literature_total_hit_count))
    values = values[is.finite(values)]
    if (length(values) > 0) return(max(values))
  }
  suppressWarnings(as.integer(fallback))
}

.plant_accounting_truncated = function(literature) {
  if (is.data.frame(literature) &&
      "literature_search_truncated" %in% names(literature)) {
    values = tolower(.uaf_non_empty(literature$literature_search_truncated))
    if (any(values %in% c("yes", "true", "1"))) return("Yes")
  }
  "No"
}

.plant_bind_unique_tables = function(tables, cols, key_cols) {
  out = .plant_bind_tables(tables, cols)
  if (nrow(out) < 1) return(out)
  key_cols = intersect(key_cols, names(out))
  if (length(key_cols) > 0) {
    key = do.call(paste, c(lapply(out[key_cols], function(x) {
      x = as.character(x)
      x[is.na(x)] = ""
      x
    }), sep = "\r"))
    out = out[!duplicated(key), , drop = FALSE]
  }
  row.names(out) = NULL
  out
}

.plant_merge_query_tables = function(results) {
  out = .plant_bind_tables(lapply(results, `[[`, "PlantQueries"),
                           .plant_query_cols())
  if (nrow(out) < 1) return(out)
  key = .plant_clean_name(out$query_plant)
  missing_key = is.na(key) | !nzchar(key)
  key[missing_key] = .plant_clean_name(out$species[missing_key])
  missing_key = is.na(key) | !nzchar(key)
  key[missing_key] = paste0("query_row_", which(missing_key))
  keys = unique(key)
  rows = lapply(seq_along(keys), function(i) {
    candidates = out[key == keys[[i]], , drop = FALSE]
    row = candidates[1L, , drop = FALSE]
    for (name in setdiff(names(row),
                         c("query_id", "input_order", "taxon_fallback"))) {
      values = candidates[[name]]
      if (is.character(values)) {
        values = values[!is.na(values) & nzchar(.uaf_squish_text(values))]
      } else {
        values = values[!is.na(values)]
      }
      if (length(values) > 0) row[[name]] = values[[1L]]
    }
    fallbacks = unlist(strsplit(
      .uaf_non_empty(candidates$taxon_fallback), ";", fixed = TRUE
    ), use.names = FALSE)
    fallbacks = unique(tolower(.uaf_non_empty(fallbacks)))
    rank_order = c("species", "genus", "family")
    fallbacks = c(intersect(rank_order, fallbacks),
                  setdiff(fallbacks, rank_order))
    row$query_id = paste0("plant_", i)
    row$input_order = as.integer(i)
    row$taxon_fallback = if (length(fallbacks) > 0) {
      paste(fallbacks, collapse = "; ")
    } else {
      NA_character_
    }
    row
  })
  out = .plant_bind_tables(rows, .plant_query_cols())
  row.names(out) = NULL
  out
}

.plant_merge_alias_tables = function(results, queries) {
  supplied = lapply(results, `[[`, "PlantQueryAliases")
  out = .plant_bind_tables(supplied, .plant_query_alias_cols())
  if (nrow(out) > 0) {
    out = .plant_remap_query_ids(out, queries)
    out = .plant_bind_unique_tables(
      list(out), .plant_query_alias_cols(), c("query_id", "alias_clean")
    )
    return(out)
  }
  .plant_query_aliases(queries$species, queries)
}

.plant_merge_accounting_tables = function(results, queries = NULL) {
  out = .plant_bind_tables(lapply(results, `[[`, "ProviderQueryAccounting"),
                           .plant_provider_query_accounting_cols())
  if (nrow(out) < 1) return(out)
  if (is.data.frame(queries) && nrow(queries) > 0) {
    out = .plant_remap_query_ids(out, queries)
  }
  rank = match(out$query_status, c("retry_required", "no_records", "records"))
  rank[is.na(rank)] = 0L
  time = suppressWarnings(as.numeric(as.POSIXct(out$retrieved_at)))
  time[!is.finite(time)] = 0
  order_idx = order(out$query_id, out$provider, -rank, -time)
  out = out[order_idx, , drop = FALSE]
  key = paste(out$query_id, out$provider, sep = "\r")
  out = out[!duplicated(key), , drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_remap_query_ids = function(x, queries) {
  if (!is.data.frame(x) || nrow(x) < 1 || !is.data.frame(queries) ||
      nrow(queries) < 1) return(x)
  query_species = .plant_clean_name(queries$species)
  row_species = if ("species" %in% names(x)) {
    .plant_clean_name(x$species)
  } else {
    rep(NA_character_, nrow(x))
  }
  index = match(row_species, query_species)
  if ("query_plant" %in% names(x)) {
    missing = is.na(index)
    index[missing] = match(
      .plant_clean_name(x$query_plant[missing]),
      .plant_clean_name(queries$query_plant)
    )
  }
  matched = !is.na(index)
  x$query_id[matched] = queries$query_id[index[matched]]
  if ("query_plant" %in% names(x)) {
    x$query_plant[matched] = queries$query_plant[index[matched]]
  }
  if ("species" %in% names(x)) {
    x$species[matched] = queries$species[index[matched]]
  }
  x
}
