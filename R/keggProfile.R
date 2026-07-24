#' Build a detailed KEGG profile for query chemicals
#'
#' @description
#' `keggProfile()` resolves query chemicals and/or supplied KEGG IDs against
#' KEGG, retrieves flat-file entries, follows KEGG cross-reference links, and
#' returns tidy tables that are useful for grouping compounds by pathway,
#' reaction, enzyme, module, identifier, and metabolic context.
#'
#' @param compounds Optional character vector of chemical names to search in
#' KEGG COMPOUND and KEGG DRUG.
#' @param kegg_ids Optional character vector of known KEGG identifiers. Values
#' can be bare IDs such as `"C01405"` or source strings such as
#' `"KEGG: C01405"`.
#' @param pubchem_profile Optional object returned by `pubchemProfile()`.
#' KEGG IDs found in PubChem annotations are reused.
#' @param cache Logical. If `TRUE`, raw KEGG text responses are cached under
#' `cache_dir`.
#' @param cache_dir Directory for cached KEGG responses. Defaults to a
#' user-cache location from `tools::R_user_dir()` when available, otherwise a
#' temporary directory.
#' @param throttle Seconds to wait between uncached KEGG requests. KEGG asks
#' users to keep API calls at or below 3 requests per second.
#' @param max_matches_per_query Maximum number of KEGG name-search matches to
#' expand per query chemical. Name searches can be broad for terms such as
#' `"glucose"`, so the default keeps the most exact matches. Use `Inf` for an
#' exhaustive expansion.
#' @param link_targets KEGG link databases to follow for each resolved KEGG ID.
#' Defaults to pathway, module, reaction, enzyme, BRITE, drug, disease, and
#' PubMed links. Use a smaller vector such as `c("pathway", "reaction",
#' "enzyme")` for faster smoke tests.
#' @param resolve_link_metadata Logical. If `TRUE`, linked KEGG IDs are fetched
#' with `get` so output tables include names, definitions, and equations when
#' KEGG provides them.
#' @param max_link_metadata Maximum number of linked target records to resolve.
#' This prevents broad queries from turning into very large metadata requests.
#' @param request_fun Optional function used to retrieve a URL. This is intended
#' for tests and advanced users. It should return text from the requested URL.
#'
#' @returns A list with tidy data frames: `matches`, `records`, `identifiers`,
#' `pathways`, `reactions`, `enzymes`, `modules`, `links`, `link_metadata`,
#' `classifications`, and `provenance`. The returned object has class
#' `"uaf_kegg_profile"`.
#'
#' @examples
#' \dontrun{
#' kegg = keggProfile(c("aspirin", "glucose"))
#' kegg$pathways
#' kegg$classifications
#' }
#'
#' @export
keggProfile = function(compounds = NULL,
                       kegg_ids = NULL,
                       pubchem_profile = NULL,
                       cache = TRUE,
                       cache_dir = NULL,
                       throttle = 0.35,
                       max_matches_per_query = 3,
                       link_targets = c("pathway", "module", "reaction",
                                        "enzyme", "brite", "drug",
                                        "disease", "pubmed"),
                       resolve_link_metadata = TRUE,
                       max_link_metadata = 100,
                       request_fun = NULL) {
  compound_queries = character()
  if (!is.null(compounds)) {
    compound_queries = unique(.uaf_non_empty(compounds))
  }

  explicit_ids = .kegg_normalize_ids(kegg_ids)
  pubchem_ids = .kegg_ids_from_pubchem_profile(pubchem_profile)
  if (length(compound_queries) < 1 &&
      length(explicit_ids) < 1 &&
      length(pubchem_ids) < 1) {
    stop("Provide at least one compound, KEGG ID, or PubChem profile.",
         call. = FALSE)
  }

  retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  max_matches_per_query = .kegg_validate_max_matches(max_matches_per_query)
  link_targets = .kegg_validate_link_targets(link_targets)
  max_link_metadata = .kegg_validate_max_matches(max_link_metadata)
  fetch = .kegg_fetcher(cache = cache,
                        cache_dir = cache_dir,
                        throttle = throttle,
                        request_fun = request_fun)

  matches = .kegg_search_queries(compound_queries = compound_queries,
                                 fetch = fetch,
                                 retrieved_at = retrieved_at,
                                 max_matches_per_query = max_matches_per_query)
  all_ids = unique(c(matches$KEGG_ID, explicit_ids, pubchem_ids))
  all_ids = all_ids[!is.na(all_ids) & all_ids != ""]
  query_map = .kegg_query_map(ids = all_ids,
                              matches = matches,
                              explicit_ids = explicit_ids,
                              pubchem_profile = pubchem_profile,
                              fallback_queries = compound_queries)

  records = .kegg_fetch_records(ids = all_ids,
                                fetch = fetch,
                                query_map = query_map,
                                retrieved_at = retrieved_at)
  links = .kegg_fetch_links(ids = all_ids,
                            fetch = fetch,
                            query_map = query_map,
                            retrieved_at = retrieved_at,
                            link_targets = link_targets)
  link_metadata = if (isTRUE(resolve_link_metadata)) {
    .kegg_fetch_link_metadata(links = links,
                              fetch = fetch,
                              retrieved_at = retrieved_at,
                              max_link_metadata = max_link_metadata)
  } else {
    .kegg_empty_link_metadata()
  }
  pathways = .kegg_pathways(records = records,
                            links = links,
                            link_metadata = link_metadata)
  reactions = .kegg_reactions(records = records,
                              links = links,
                              link_metadata = link_metadata)
  enzymes = .kegg_enzymes(records = records,
                          links = links,
                          link_metadata = link_metadata)
  modules = .kegg_modules(records = records,
                          links = links,
                          link_metadata = link_metadata)
  identifiers = .kegg_identifiers(records = records)
  classifications = .kegg_classifications(pathways = pathways,
                                          enzymes = enzymes,
                                          records = records)

  out = list(
    matches = matches,
    records = records,
    identifiers = identifiers,
    pathways = pathways,
    reactions = reactions,
    enzymes = enzymes,
    modules = modules,
    links = links,
    link_metadata = link_metadata,
    classifications = classifications,
    provenance = .kegg_profile_provenance(matches = matches,
                                          records = records,
                                          links = links,
                                          link_metadata = link_metadata),
    retrieved_at = retrieved_at
  )
  class(out) = c("uaf_kegg_profile", class(out))
  out
}

#' @export
print.uaf_kegg_profile = function(x, ...) {
  cat("uaf KEGG profile\n")
  cat("  matches: ", nrow(x$matches), " rows\n", sep = "")
  cat("  records: ", length(unique(x$records$KEGG_ID)), " KEGG IDs\n", sep = "")
  cat("  pathways: ", nrow(x$pathways), " rows\n", sep = "")
  cat("  reactions: ", nrow(x$reactions), " rows\n", sep = "")
  cat("  enzymes: ", nrow(x$enzymes), " rows\n", sep = "")
  invisible(x)
}

.kegg_base_url = function() {
  "https://rest.kegg.jp"
}

.kegg_fetcher = function(cache, cache_dir, throttle, request_fun) {
  if (is.null(cache_dir)) {
    cache_dir = .kegg_default_cache_dir()
  }
  force(cache)
  force(cache_dir)
  force(throttle)
  force(request_fun)

  function(url) {
    cache_file = file.path(cache_dir, paste0(.kegg_url_hash(url), ".txt"))
    if (isTRUE(cache) && file.exists(cache_file)) {
      return(paste(readLines(cache_file, warn = FALSE, encoding = "UTF-8"),
                   collapse = "\n"))
    }

    result = tryCatch({
      if (is.null(request_fun)) {
        con = base::url(url, open = "rb")
        on.exit(close(con), add = TRUE)
        txt = paste(readLines(con, warn = FALSE, encoding = "UTF-8"),
                    collapse = "\n")
        if (isTRUE(cache)) {
          dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
          writeLines(txt, cache_file, useBytes = TRUE)
        }
        if (!is.na(throttle) && throttle > 0) Sys.sleep(throttle)
        txt
      } else {
        request_result = request_fun(url)
        txt = paste(as.character(request_result), collapse = "\n")
        if (isTRUE(cache)) {
          dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
          writeLines(txt, cache_file, useBytes = TRUE)
        }
        txt
      }
    }, error = function(error) {
      warning("KEGG request failed: ", conditionMessage(error),
              "\nURL: ", url, call. = FALSE)
      NULL
    })

    result
  }
}

.kegg_default_cache_dir = function() {
  cache_dir = tryCatch(tools::R_user_dir("uafR", "cache"),
                       error = function(error) NULL)
  if (is.null(cache_dir)) {
    cache_dir = file.path(tempdir(), "uafR-cache")
  }
  file.path(cache_dir, "kegg")
}

.kegg_url_hash = function(url) {
  ints = utf8ToInt(url)
  val = sum((ints * seq_along(ints)) %% .Machine$integer.max)
  paste0(nchar(url), "_", sprintf("%08x", as.integer(val %% .Machine$integer.max)))
}

.kegg_search_queries = function(compound_queries, fetch, retrieved_at,
                                max_matches_per_query) {
  cols = c("Query", "KEGG_ID", "Database", "MatchName", "MatchStatus",
           "MatchScore", "MatchRank", "SourceURL", "RetrievedAt")
  if (length(compound_queries) < 1) return(.uaf_empty_table(cols))

  rows = list()
  for (query in compound_queries) {
    for (database in c("compound", "drug")) {
      url = paste0(.kegg_base_url(), "/find/", database, "/",
                   utils::URLencode(query, reserved = TRUE))
      txt = fetch(url)
      parsed = .kegg_parse_find(txt = txt,
                                query = query,
                                database = database,
                                url = url,
                                retrieved_at = retrieved_at)
      if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
    }
  }

  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out = unique(out)
  .kegg_limit_matches(out, max_matches_per_query)
}

.kegg_parse_find = function(txt, query, database, url, retrieved_at) {
  cols = c("Query", "KEGG_ID", "Database", "MatchName", "MatchStatus",
           "MatchScore", "MatchRank", "SourceURL", "RetrievedAt")
  if (is.null(txt) || length(txt) < 1 || is.na(txt[[1]]) || txt[[1]] == "") {
    return(.uaf_empty_table(cols))
  }
  txt = paste(as.character(txt), collapse = "\n")
  lines = unlist(strsplit(txt, "\n", fixed = TRUE), use.names = FALSE)
  rows = list()
  for (line in lines) {
    pieces = strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(pieces) < 2) next
    kegg_id = .kegg_normalize_ids(pieces[[1]])
    if (length(kegg_id) < 1) next
    rows[[length(rows) + 1]] = data.frame(
      Query = query,
      KEGG_ID = kegg_id[[1]],
      Database = .kegg_database(kegg_id[[1]]),
      MatchName = .uaf_squish_text(pieces[[2]]),
      MatchStatus = ifelse(database == .kegg_database(kegg_id[[1]]),
                           "name_match", paste0(database, "_name_match")),
      MatchScore = .kegg_match_score(query, pieces[[2]]),
      MatchRank = NA_integer_,
      SourceURL = url,
      RetrievedAt = retrieved_at,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.kegg_validate_max_matches = function(max_matches_per_query) {
  if (length(max_matches_per_query) != 1 ||
      is.na(max_matches_per_query) ||
      max_matches_per_query <= 0) {
    stop("`max_matches_per_query` must be a positive number or Inf.",
         call. = FALSE)
  }
  if (is.infinite(max_matches_per_query)) return(Inf)
  as.integer(max_matches_per_query)
}

.kegg_validate_link_targets = function(link_targets) {
  allowed = c("pathway", "module", "reaction", "enzyme", "brite",
              "drug", "disease", "pubmed")
  if (is.null(link_targets) || length(link_targets) < 1) return(character())
  link_targets = unique(.uaf_non_empty(link_targets))
  unknown = setdiff(link_targets, allowed)
  if (length(unknown) > 0) {
    stop("Unsupported KEGG link target(s): ",
         paste(unknown, collapse = ", "),
         ". Supported targets are: ",
         paste(allowed, collapse = ", "),
         call. = FALSE)
  }
  link_targets
}

.kegg_limit_matches = function(matches, max_matches_per_query) {
  if (!is.data.frame(matches) || nrow(matches) < 1) return(matches)
  matches$.OriginalOrder = seq_len(nrow(matches))
  database_rank = ifelse(matches$Database == "compound", 0,
                         ifelse(matches$Database == "drug", 1, 2))
  matches = matches[order(matches$Query,
                          matches$MatchScore,
                          database_rank,
                          matches$.OriginalOrder), , drop = FALSE]

  kept = list()
  for (query in unique(matches$Query)) {
    query_rows = matches[matches$Query == query, , drop = FALSE]
    if (!is.infinite(max_matches_per_query)) {
      query_rows = utils::head(query_rows, max_matches_per_query)
    }
    query_rows$MatchRank = seq_len(nrow(query_rows))
    kept[[length(kept) + 1]] = query_rows
  }
  out = do.call(rbind, kept)
  out$.OriginalOrder = NULL
  row.names(out) = NULL
  out
}

.kegg_match_score = function(query, match_name) {
  query_clean = tolower(.uaf_squish_text(query))
  query_key = .kegg_match_key(query_clean)
  names = unlist(strsplit(.uaf_squish_text(match_name), ";", fixed = TRUE),
                 use.names = FALSE)
  names = tolower(.uaf_non_empty(names))
  name_keys = .kegg_match_key(names)

  if (length(names) < 1) return(4)
  if (any(names == query_clean) || any(name_keys == query_key)) return(0)
  if (any(startsWith(names, query_clean)) ||
      any(startsWith(name_keys, query_key))) return(1)
  if (grepl(.kegg_regex_escape(query_clean),
            paste(names, collapse = " "),
            ignore.case = TRUE)) return(2)
  3
}

.kegg_match_key = function(x) {
  gsub("[^a-z0-9]", "", tolower(x))
}

.kegg_regex_escape = function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

.kegg_normalize_ids = function(x) {
  x = .uaf_non_empty(x)
  if (length(x) < 1) return(character())
  ids = unlist(lapply(x, .uaf_extract_pattern,
                      pattern = "\\b[CDGMR]\\d{5}\\b"),
               use.names = FALSE)
  unique(ids[!is.na(ids) & ids != ""])
}

.kegg_ids_from_pubchem_profile = function(pubchem_profile) {
  if (is.null(pubchem_profile) ||
      !is.list(pubchem_profile) ||
      is.null(pubchem_profile$annotations) ||
      !is.data.frame(pubchem_profile$annotations) ||
      nrow(pubchem_profile$annotations) < 1) {
    return(character())
  }
  text = paste(pubchem_profile$annotations$Name,
               pubchem_profile$annotations$Value,
               pubchem_profile$annotations$Source,
               collapse = " ")
  .kegg_normalize_ids(text)
}

.kegg_query_map = function(ids, matches, explicit_ids, pubchem_profile,
                           fallback_queries) {
  out = stats::setNames(rep(NA_character_, length(ids)), ids)

  if (is.data.frame(matches) && nrow(matches) > 0) {
    for (id in ids) {
      match_rows = matches[matches$KEGG_ID == id, , drop = FALSE]
      if (nrow(match_rows) > 0) out[[id]] = match_rows$Query[[1]]
    }
  }

  if (!is.null(pubchem_profile) &&
      is.data.frame(pubchem_profile$identity) &&
      length(.uaf_non_empty(pubchem_profile$identity$Query)) > 0) {
    fallback = pubchem_profile$identity$Query[[1]]
  } else if (length(fallback_queries) > 0) {
    fallback = fallback_queries[[1]]
  } else {
    fallback = NA_character_
  }

  out[is.na(out) | out == ""] = fallback
  out
}

.kegg_fetch_records = function(ids, fetch, query_map, retrieved_at) {
  cols = c("Query", "KEGG_ID", "Database", "Field", "Value", "CleanValue",
           "ValueNumeric", "UnitClean", "EvidenceURL", "RetrievedAt")
  ids = unique(ids)
  ids = ids[!is.na(ids) & ids != ""]
  if (length(ids) < 1) return(.uaf_empty_table(cols))

  rows = list()
  chunks = split(ids, ceiling(seq_along(ids) / 10))
  for (chunk in chunks) {
    url = paste0(.kegg_base_url(), "/get/", paste(chunk, collapse = "+"))
    txt = fetch(url)
    parsed = .kegg_parse_flat_records(txt = txt,
                                      url = url,
                                      query_map = query_map,
                                      retrieved_at = retrieved_at)
    if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
  }

  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.kegg_parse_flat_records = function(txt, url, query_map, retrieved_at) {
  cols = c("Query", "KEGG_ID", "Database", "Field", "Value", "CleanValue",
           "ValueNumeric", "UnitClean", "EvidenceURL", "RetrievedAt")
  if (is.null(txt) || is.na(txt) || txt == "") return(.uaf_empty_table(cols))
  raw_records = unlist(strsplit(txt, "\n///", fixed = TRUE), use.names = FALSE)
  rows = list()

  for (record in raw_records) {
    lines = unlist(strsplit(record, "\n", fixed = TRUE), use.names = FALSE)
    lines = lines[lines != ""]
    if (length(lines) < 1) next

    entry_line = lines[grepl("^ENTRY", lines)][1]
    kegg_id = .kegg_normalize_ids(entry_line)
    if (length(kegg_id) < 1) next
    kegg_id = kegg_id[[1]]
    current_field = NA_character_

    for (line in lines) {
      field = trimws(substr(line, 1, 12))
      value = trimws(substr(line, 13, nchar(line)))
      if (field == "") {
        field = current_field
      } else {
        current_field = field
      }
      if (is.na(field) || field == "" || value == "") next
      measurement = .uaf_parse_measurement(value)
      rows[[length(rows) + 1]] = data.frame(
        Query = unname(query_map[[kegg_id]]),
        KEGG_ID = kegg_id,
        Database = .kegg_database(kegg_id),
        Field = field,
        Value = value,
        CleanValue = .uaf_squish_text(gsub(";+\\s*$", "", value)),
        ValueNumeric = measurement$ValueNumeric[[1]],
        UnitClean = measurement$UnitClean[[1]],
        EvidenceURL = url,
        RetrievedAt = retrieved_at,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.kegg_fetch_links = function(ids, fetch, query_map, retrieved_at,
                             link_targets) {
  cols = c("Query", "KEGG_ID", "Database", "TargetDatabase", "TargetID",
           "TargetPrefix", "SourceEntry", "EvidenceURL", "RetrievedAt")
  ids = unique(ids)
  ids = ids[!is.na(ids) & ids != ""]
  if (length(ids) < 1) return(.uaf_empty_table(cols))

  targets = link_targets
  if (length(targets) < 1) return(.uaf_empty_table(cols))
  rows = list()
  for (id in ids) {
    source_entry = .kegg_prefixed_id(id)
    for (target in targets) {
      url = paste0(.kegg_base_url(), "/link/", target, "/", source_entry)
      txt = fetch(url)
      parsed = .kegg_parse_links(txt = txt,
                                 id = id,
                                 target = target,
                                 source_entry = source_entry,
                                 url = url,
                                 query_map = query_map,
                                 retrieved_at = retrieved_at)
      if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
    }
  }

  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  unique(out)
}

.kegg_parse_links = function(txt, id, target, source_entry, url, query_map,
                             retrieved_at) {
  cols = c("Query", "KEGG_ID", "Database", "TargetDatabase", "TargetID",
           "TargetPrefix", "SourceEntry", "EvidenceURL", "RetrievedAt")
  if (is.null(txt) || is.na(txt) || txt == "") return(.uaf_empty_table(cols))
  lines = unlist(strsplit(txt, "\n", fixed = TRUE), use.names = FALSE)
  rows = list()
  for (line in lines) {
    pieces = strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(pieces) < 2) next
    target_entry = pieces[[2]]
    target_pieces = strsplit(target_entry, ":", fixed = TRUE)[[1]]
    prefix = ifelse(length(target_pieces) > 1, target_pieces[[1]], NA_character_)
    target_id = ifelse(length(target_pieces) > 1, target_pieces[[2]], target_entry)
    rows[[length(rows) + 1]] = data.frame(
      Query = unname(query_map[[id]]),
      KEGG_ID = id,
      Database = .kegg_database(id),
      TargetDatabase = target,
      TargetID = target_id,
      TargetPrefix = prefix,
      SourceEntry = source_entry,
      EvidenceURL = url,
      RetrievedAt = retrieved_at,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.kegg_empty_link_metadata = function() {
  .uaf_empty_table(c("TargetDatabase", "TargetID", "TargetEntry", "Name",
                     "Definition", "Equation", "Class", "PathwayGroup",
                     "EvidenceURL", "RetrievedAt"))
}

.kegg_fetch_link_metadata = function(links, fetch, retrieved_at,
                                    max_link_metadata) {
  if (!is.data.frame(links) || nrow(links) < 1) {
    return(.kegg_empty_link_metadata())
  }

  supported = c("pathway", "module", "reaction", "enzyme", "drug",
                "disease", "brite")
  targets = links[links$TargetDatabase %in% supported, , drop = FALSE]
  if (nrow(targets) < 1) return(.kegg_empty_link_metadata())

  targets$TargetEntry = mapply(.kegg_link_target_entry,
                               targets$TargetDatabase,
                               targets$TargetID,
                               targets$TargetPrefix,
                               USE.NAMES = FALSE)
  targets = targets[!is.na(targets$TargetEntry) &
                      targets$TargetEntry != "", , drop = FALSE]
  targets = unique(targets[, c("TargetDatabase", "TargetID", "TargetEntry"),
                           drop = FALSE])
  if (nrow(targets) < 1) return(.kegg_empty_link_metadata())

  if (!is.infinite(max_link_metadata) && nrow(targets) > max_link_metadata) {
    warning("Resolving only the first ", max_link_metadata,
            " KEGG linked target records. Increase `max_link_metadata` for exhaustive metadata.",
            call. = FALSE)
    targets = utils::head(targets, max_link_metadata)
  }

  rows = list()
  for (i in seq_len(nrow(targets))) {
    target = targets[i, , drop = FALSE]
    url = paste0(.kegg_base_url(), "/get/", target$TargetEntry)
    txt = fetch(url)
    parsed = .kegg_parse_link_metadata(txt = txt,
                                       target = target,
                                       url = url,
                                       retrieved_at = retrieved_at)
    if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
  }

  if (length(rows) < 1) return(.kegg_empty_link_metadata())
  out = unique(do.call(rbind, rows))
  row.names(out) = NULL
  out
}

.kegg_parse_link_metadata = function(txt, target, url, retrieved_at) {
  if (is.null(txt) || length(txt) < 1 || is.na(txt[[1]]) || txt[[1]] == "") {
    return(.kegg_empty_link_metadata())
  }
  fields = .kegg_flat_record_fields(txt)
  if (length(fields) < 1) return(.kegg_empty_link_metadata())

  name = .kegg_clean_names(fields[["NAME"]])
  definition = .uaf_first_non_empty_text(fields[["DEFINITION"]],
                                         fields[["DESCRIPTION"]])
  equation = .uaf_first_non_empty_text(fields[["EQUATION"]])
  class = .kegg_clean_names(c(fields[["CLASS"]], fields[["BRITE"]]))

  data.frame(
    TargetDatabase = target$TargetDatabase,
    TargetID = target$TargetID,
    TargetEntry = target$TargetEntry,
    Name = name,
    Definition = definition,
    Equation = equation,
    Class = class,
    PathwayGroup = .kegg_pathway_group(.uaf_first_non_empty_text(name,
                                                                 definition,
                                                                 class)),
    EvidenceURL = url,
    RetrievedAt = retrieved_at,
    stringsAsFactors = FALSE
  )
}

.kegg_flat_record_fields = function(txt) {
  txt = paste(as.character(txt), collapse = "\n")
  records = unlist(strsplit(txt, "\n///", fixed = TRUE), use.names = FALSE)
  record = records[[1]]
  lines = unlist(strsplit(record, "\n", fixed = TRUE), use.names = FALSE)
  lines = lines[lines != ""]
  fields = list()
  current_field = NA_character_

  for (line in lines) {
    field = trimws(substr(line, 1, 12))
    value = trimws(substr(line, 13, nchar(line)))
    if (field == "") {
      field = current_field
    } else {
      current_field = field
    }
    if (is.na(field) || field == "" || value == "") next
    fields[[field]] = c(fields[[field]], value)
  }
  fields
}

.kegg_clean_names = function(values) {
  values = .uaf_non_empty(values)
  if (length(values) < 1) return(NA_character_)
  values = gsub(";+$", "", values)
  values = .uaf_squish_text(values)
  values = values[!is.na(values) & values != ""]
  if (length(values) < 1) return(NA_character_)
  paste(unique(values), collapse = "; ")
}

.kegg_link_target_entry = function(target_database, target_id, target_prefix) {
  target_id = .uaf_squish_text(target_id)
  target_prefix = .uaf_squish_text(target_prefix)
  if (is.na(target_id) || target_id == "") return(NA_character_)

  if (!is.na(target_prefix) && target_prefix != "") {
    return(paste0(target_prefix, ":", target_id))
  }

  prefix = switch(target_database,
                  pathway = "path",
                  module = "md",
                  reaction = "rn",
                  enzyme = "ec",
                  drug = "dr",
                  disease = "ds",
                  brite = "br",
                  NA_character_)
  if (is.na(prefix)) return(NA_character_)
  paste0(prefix, ":", target_id)
}

.kegg_metadata_lookup = function(link_metadata, target_database, target_id) {
  empty = list(Name = NA_character_,
               Definition = NA_character_,
               Equation = NA_character_,
               Class = NA_character_,
               PathwayGroup = NA_character_)
  if (!is.data.frame(link_metadata) ||
      nrow(link_metadata) < 1 ||
      !"TargetDatabase" %in% colnames(link_metadata) ||
      !"TargetID" %in% colnames(link_metadata)) {
    return(empty)
  }
  hit = link_metadata[link_metadata$TargetDatabase == target_database &
                        link_metadata$TargetID == target_id, , drop = FALSE]
  if (nrow(hit) < 1) return(empty)
  hit = hit[1, , drop = FALSE]
  list(Name = .uaf_first_non_empty_text(hit$Name),
       Definition = .uaf_first_non_empty_text(hit$Definition),
       Equation = .uaf_first_non_empty_text(hit$Equation),
       Class = .uaf_first_non_empty_text(hit$Class),
       PathwayGroup = .uaf_first_non_empty_text(hit$PathwayGroup))
}

.kegg_database = function(id) {
  first = substr(id, 1, 1)
  switch(first,
         C = "compound",
         D = "drug",
         G = "glycan",
         M = "module",
         R = "reaction",
         "kegg")
}

.kegg_prefixed_id = function(id) {
  first = substr(id, 1, 1)
  prefix = switch(first,
                  C = "cpd",
                  D = "dr",
                  G = "gl",
                  M = "md",
                  R = "rn",
                  "kegg")
  paste0(prefix, ":", id)
}

.kegg_identifiers = function(records) {
  cols = c("Query", "KEGG_ID", "Database", "IdentifierType", "Identifier",
           "SourceField", "EvidenceURL", "RetrievedAt")
  if (!is.data.frame(records) || nrow(records) < 1) return(.uaf_empty_table(cols))

  rows = list()
  for (i in seq_len(nrow(records))) {
    record = records[i, , drop = FALSE]
    if (record$Field == "DBLINKS") {
      pieces = strsplit(record$CleanValue, ":", fixed = TRUE)[[1]]
      if (length(pieces) > 1) {
        source = .uaf_squish_text(pieces[[1]])
        identifiers = unlist(strsplit(.uaf_squish_text(paste(pieces[-1],
                                                             collapse = ":")),
                                       "[,;[:space:]]+"),
                             use.names = FALSE)
        identifiers = identifiers[identifiers != ""]
        for (identifier in identifiers) {
          rows[[length(rows) + 1]] = data.frame(
            Query = record$Query,
            KEGG_ID = record$KEGG_ID,
            Database = record$Database,
            IdentifierType = source,
            Identifier = identifier,
            SourceField = record$Field,
            EvidenceURL = record$EvidenceURL,
            RetrievedAt = record$RetrievedAt,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    extracted = .uaf_extract_identifiers(record$CleanValue)
    if (nrow(extracted) < 1) next
    for (j in seq_len(nrow(extracted))) {
      rows[[length(rows) + 1]] = data.frame(
        Query = record$Query,
        KEGG_ID = record$KEGG_ID,
        Database = record$Database,
        IdentifierType = extracted$IdentifierType[[j]],
        Identifier = extracted$Identifier[[j]],
        SourceField = record$Field,
        EvidenceURL = record$EvidenceURL,
        RetrievedAt = record$RetrievedAt,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = unique(do.call(rbind, rows))
  row.names(out) = NULL
  out
}

.kegg_pathways = function(records, links, link_metadata) {
  cols = c("Query", "KEGG_ID", "Database", "PathwayID", "PathwayName",
           "PathwayGroup", "Evidence", "EvidenceURL", "RetrievedAt")
  rows = list()

  if (is.data.frame(records) && nrow(records) > 0) {
    pathway_rows = records[records$Field == "PATHWAY", , drop = FALSE]
    for (i in seq_len(nrow(pathway_rows))) {
      record = pathway_rows[i, , drop = FALSE]
      pathway_id = .uaf_extract_pattern(record$CleanValue,
                                        "\\b(?:map|[a-z]{2,4})\\d{5}\\b")
      pathway_name = record$CleanValue
      if (length(pathway_id) > 0) {
        pathway_name = .uaf_squish_text(sub(pathway_id[[1]], "",
                                            pathway_name, fixed = TRUE))
      } else {
        pathway_id = NA_character_
      }
      metadata = .kegg_metadata_lookup(link_metadata, "pathway",
                                       pathway_id[[1]])
      pathway_name = .uaf_first_non_empty_text(metadata$Name, pathway_name)
      pathway_group = .uaf_first_non_empty_text(metadata$PathwayGroup,
                                                .kegg_pathway_group(pathway_name))
      rows[[length(rows) + 1]] = data.frame(
        Query = record$Query,
        KEGG_ID = record$KEGG_ID,
        Database = record$Database,
        PathwayID = pathway_id[[1]],
        PathwayName = pathway_name,
        PathwayGroup = pathway_group,
        Evidence = "record",
        EvidenceURL = record$EvidenceURL,
        RetrievedAt = record$RetrievedAt,
        stringsAsFactors = FALSE
      )
    }
  }

  if (is.data.frame(links) && nrow(links) > 0) {
    link_rows = links[links$TargetDatabase == "pathway", , drop = FALSE]
    for (i in seq_len(nrow(link_rows))) {
      link = link_rows[i, , drop = FALSE]
      metadata = .kegg_metadata_lookup(link_metadata, "pathway",
                                       link$TargetID)
      pathway_name = metadata$Name
      pathway_group = .uaf_first_non_empty_text(metadata$PathwayGroup,
                                                .kegg_pathway_group(pathway_name))
      rows[[length(rows) + 1]] = data.frame(
        Query = link$Query,
        KEGG_ID = link$KEGG_ID,
        Database = link$Database,
        PathwayID = link$TargetID,
        PathwayName = pathway_name,
        PathwayGroup = pathway_group,
        Evidence = "link",
        EvidenceURL = link$EvidenceURL,
        RetrievedAt = link$RetrievedAt,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = unique(do.call(rbind, rows))
  row.names(out) = NULL
  out
}

.kegg_reactions = function(records, links, link_metadata) {
  cols = c("Query", "KEGG_ID", "Database", "ReactionID", "ReactionName",
           "ReactionDefinition", "Equation", "Evidence", "EvidenceURL",
           "RetrievedAt")
  rows = list()
  if (is.data.frame(records) && nrow(records) > 0) {
    reaction_rows = records[records$Field == "REACTION", , drop = FALSE]
    for (i in seq_len(nrow(reaction_rows))) {
      record = reaction_rows[i, , drop = FALSE]
      reaction_ids = .uaf_extract_pattern(record$CleanValue, "\\bR\\d{5}\\b")
      if (length(reaction_ids) < 1) reaction_ids = NA_character_
      for (reaction_id in reaction_ids) {
        metadata = .kegg_metadata_lookup(link_metadata, "reaction",
                                         reaction_id)
        rows[[length(rows) + 1]] = data.frame(
          Query = record$Query,
          KEGG_ID = record$KEGG_ID,
          Database = record$Database,
          ReactionID = reaction_id,
          ReactionName = metadata$Name,
          ReactionDefinition = metadata$Definition,
          Equation = metadata$Equation,
          Evidence = "record",
          EvidenceURL = record$EvidenceURL,
          RetrievedAt = record$RetrievedAt,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (is.data.frame(links) && nrow(links) > 0) {
    link_rows = links[links$TargetDatabase == "reaction", , drop = FALSE]
    for (i in seq_len(nrow(link_rows))) {
      link = link_rows[i, , drop = FALSE]
      metadata = .kegg_metadata_lookup(link_metadata, "reaction",
                                       link$TargetID)
      rows[[length(rows) + 1]] = data.frame(
        Query = link$Query,
        KEGG_ID = link$KEGG_ID,
        Database = link$Database,
        ReactionID = link$TargetID,
        ReactionName = metadata$Name,
        ReactionDefinition = metadata$Definition,
        Equation = metadata$Equation,
        Evidence = "link",
        EvidenceURL = link$EvidenceURL,
        RetrievedAt = link$RetrievedAt,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = unique(do.call(rbind, rows))
  row.names(out) = NULL
  out
}

.kegg_enzymes = function(records, links, link_metadata) {
  cols = c("Query", "KEGG_ID", "Database", "ECNumber", "EnzymeName",
           "EnzymeClass", "Evidence", "EvidenceURL", "RetrievedAt")
  rows = list()
  if (is.data.frame(records) && nrow(records) > 0) {
    enzyme_rows = records[records$Field == "ENZYME", , drop = FALSE]
    for (i in seq_len(nrow(enzyme_rows))) {
      record = enzyme_rows[i, , drop = FALSE]
      ecs = .uaf_extract_pattern(record$CleanValue,
                                 "\\b[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9A-Za-z-]+\\b")
      if (length(ecs) < 1) ecs = NA_character_
      for (ec in ecs) {
        metadata = .kegg_metadata_lookup(link_metadata, "enzyme", ec)
        rows[[length(rows) + 1]] = data.frame(
          Query = record$Query,
          KEGG_ID = record$KEGG_ID,
          Database = record$Database,
          ECNumber = ec,
          EnzymeName = metadata$Name,
          EnzymeClass = .kegg_enzyme_class(ec),
          Evidence = "record",
          EvidenceURL = record$EvidenceURL,
          RetrievedAt = record$RetrievedAt,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (is.data.frame(links) && nrow(links) > 0) {
    link_rows = links[links$TargetDatabase == "enzyme", , drop = FALSE]
    for (i in seq_len(nrow(link_rows))) {
      link = link_rows[i, , drop = FALSE]
      metadata = .kegg_metadata_lookup(link_metadata, "enzyme",
                                       link$TargetID)
      rows[[length(rows) + 1]] = data.frame(
        Query = link$Query,
        KEGG_ID = link$KEGG_ID,
        Database = link$Database,
        ECNumber = link$TargetID,
        EnzymeName = metadata$Name,
        EnzymeClass = .kegg_enzyme_class(link$TargetID),
        Evidence = "link",
        EvidenceURL = link$EvidenceURL,
        RetrievedAt = link$RetrievedAt,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = unique(do.call(rbind, rows))
  row.names(out) = NULL
  out
}

.kegg_modules = function(records, links, link_metadata) {
  cols = c("Query", "KEGG_ID", "Database", "ModuleID", "ModuleName",
           "ModuleDefinition", "Evidence", "EvidenceURL", "RetrievedAt")
  rows = list()
  if (is.data.frame(records) && nrow(records) > 0) {
    module_rows = records[records$Field == "MODULE", , drop = FALSE]
    for (i in seq_len(nrow(module_rows))) {
      record = module_rows[i, , drop = FALSE]
      module_ids = .uaf_extract_pattern(record$CleanValue, "\\bM\\d{5}\\b")
      module_name = record$CleanValue
      if (length(module_ids) < 1) module_ids = NA_character_
      for (module_id in module_ids) {
        metadata = .kegg_metadata_lookup(link_metadata, "module",
                                         module_id)
        rows[[length(rows) + 1]] = data.frame(
          Query = record$Query,
          KEGG_ID = record$KEGG_ID,
          Database = record$Database,
          ModuleID = module_id,
          ModuleName = .uaf_first_non_empty_text(
            metadata$Name,
            .uaf_squish_text(sub(module_id, "", module_name, fixed = TRUE))
          ),
          ModuleDefinition = metadata$Definition,
          Evidence = "record",
          EvidenceURL = record$EvidenceURL,
          RetrievedAt = record$RetrievedAt,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (is.data.frame(links) && nrow(links) > 0) {
    link_rows = links[links$TargetDatabase == "module", , drop = FALSE]
    for (i in seq_len(nrow(link_rows))) {
      link = link_rows[i, , drop = FALSE]
      metadata = .kegg_metadata_lookup(link_metadata, "module",
                                       link$TargetID)
      rows[[length(rows) + 1]] = data.frame(
        Query = link$Query,
        KEGG_ID = link$KEGG_ID,
        Database = link$Database,
        ModuleID = link$TargetID,
        ModuleName = metadata$Name,
        ModuleDefinition = metadata$Definition,
        Evidence = "link",
        EvidenceURL = link$EvidenceURL,
        RetrievedAt = link$RetrievedAt,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = unique(do.call(rbind, rows))
  row.names(out) = NULL
  out
}

.kegg_link_rows = function(links, target, id_col, name_col, cols) {
  if (!is.data.frame(links) || nrow(links) < 1) return(list())
  link_rows = links[links$TargetDatabase == target, , drop = FALSE]
  rows = list()
  for (i in seq_len(nrow(link_rows))) {
    link = link_rows[i, , drop = FALSE]
    row = as.list(stats::setNames(rep(NA_character_, length(cols)), cols))
    row$Query = link$Query
    row$KEGG_ID = link$KEGG_ID
    row$Database = link$Database
    row[[id_col]] = link$TargetID
    row[[name_col]] = NA_character_
    row$Evidence = "link"
    row$EvidenceURL = link$EvidenceURL
    row$RetrievedAt = link$RetrievedAt
    rows[[length(rows) + 1]] = as.data.frame(row, stringsAsFactors = FALSE)
  }
  rows
}

.kegg_classifications = function(pathways, enzymes, records) {
  cols = c("Query", "KEGG_ID", "Database", "ClassificationType",
           "Classification", "Evidence", "EvidenceURL", "RetrievedAt")
  rows = list()

  if (is.data.frame(pathways) && nrow(pathways) > 0) {
    for (i in seq_len(nrow(pathways))) {
      pathway = pathways[i, , drop = FALSE]
      if (is.na(pathway$PathwayGroup) || pathway$PathwayGroup == "") next
      rows[[length(rows) + 1]] = data.frame(
        Query = pathway$Query,
        KEGG_ID = pathway$KEGG_ID,
        Database = pathway$Database,
        ClassificationType = "KEGG pathway group",
        Classification = pathway$PathwayGroup,
        Evidence = pathway$PathwayName,
        EvidenceURL = pathway$EvidenceURL,
        RetrievedAt = pathway$RetrievedAt,
        stringsAsFactors = FALSE
      )
    }
  }

  if (is.data.frame(enzymes) && nrow(enzymes) > 0) {
    for (i in seq_len(nrow(enzymes))) {
      enzyme = enzymes[i, , drop = FALSE]
      if (is.na(enzyme$EnzymeClass) || enzyme$EnzymeClass == "") next
      rows[[length(rows) + 1]] = data.frame(
        Query = enzyme$Query,
        KEGG_ID = enzyme$KEGG_ID,
        Database = enzyme$Database,
        ClassificationType = "EC enzyme class",
        Classification = enzyme$EnzymeClass,
        Evidence = enzyme$ECNumber,
        EvidenceURL = enzyme$EvidenceURL,
        RetrievedAt = enzyme$RetrievedAt,
        stringsAsFactors = FALSE
      )
    }
  }

  if (is.data.frame(records) && nrow(records) > 0) {
    brite_rows = records[records$Field == "BRITE", , drop = FALSE]
    for (i in seq_len(nrow(brite_rows))) {
      brite = brite_rows[i, , drop = FALSE]
      rows[[length(rows) + 1]] = data.frame(
        Query = brite$Query,
        KEGG_ID = brite$KEGG_ID,
        Database = brite$Database,
        ClassificationType = "KEGG BRITE",
        Classification = brite$CleanValue,
        Evidence = brite$Field,
        EvidenceURL = brite$EvidenceURL,
        RetrievedAt = brite$RetrievedAt,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = unique(do.call(rbind, rows))
  row.names(out) = NULL
  out
}

.kegg_pathway_group = function(pathway_name) {
  pathway_name = .uaf_squish_text(pathway_name)
  if (is.na(pathway_name) || pathway_name == "") return(NA_character_)
  patterns = c(
    "Carbohydrate metabolism" = "glycolysis|gluconeogenesis|pentose|starch|sucrose|galactose|fructose|mannose|carbohydrate",
    "Lipid metabolism" = "lipid|fatty acid|arachidonic|glycerolipid|glycerophospholipid|sphingolipid|steroid|bile",
    "Amino acid metabolism" = "amino acid|alanine|aspartate|glutamate|lysine|tryptophan|tyrosine|phenylalanine|histidine|arginine|proline|glycine|serine|threonine|cysteine|methionine|valine|leucine|isoleucine",
    "Energy metabolism" = "carbon fixation|methane|nitrogen|sulfur|oxidative phosphorylation",
    "Nucleotide metabolism" = "purine|pyrimidine|nucleotide",
    "Cofactor and vitamin metabolism" = "cofactor|vitamin|porphyrin|folate|nicotinate|riboflavin|biotin|lipoic",
    "Secondary metabolism" = "secondary metabolite|terpenoid|polyketide|phenylpropanoid|flavonoid|alkaloid",
    "Xenobiotic metabolism" = "xenobiotic|drug|cytochrome p450|dioxin|chloro|toluene|benzoate|polycyclic|ddt",
    "Signaling and regulation" = "signaling|hormone|neuroactive|synapse|inflammation"
  )
  for (label in names(patterns)) {
    if (grepl(patterns[[label]], pathway_name, ignore.case = TRUE)) {
      return(label)
    }
  }
  "Other KEGG pathway"
}

.kegg_enzyme_class = function(ec) {
  if (is.na(ec) || ec == "") return(NA_character_)
  first = strsplit(ec, "\\.")[[1]][[1]]
  switch(first,
         "1" = "Oxidoreductases",
         "2" = "Transferases",
         "3" = "Hydrolases",
         "4" = "Lyases",
         "5" = "Isomerases",
         "6" = "Ligases",
         "7" = "Translocases",
         NA_character_)
}

.kegg_profile_provenance = function(...) {
  tables = list(...)
  rows = list()
  for (name in names(tables)) {
    table = tables[[name]]
    if (!is.data.frame(table)) next
    url_col = intersect(c("SourceURL", "EvidenceURL"), colnames(table))
    if (length(url_col) < 1) next
    urls = unique(unlist(table[url_col], use.names = FALSE))
    urls = urls[!is.na(urls) & urls != ""]
    if (length(urls) < 1) next
    rows[[length(rows) + 1]] = data.frame(
      Table = name,
      SourceURL = urls,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) < 1) return(.uaf_empty_table(c("Table", "SourceURL")))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}
