#' Build a detailed PubChem profile for query chemicals
#'
#' @description
#' `pubchemProfile()` retrieves structured PubChem information for query
#' chemicals and returns tidy tables that can be reused by downstream uafR
#' workflows. It separates PubChem enrichment from mass spectrometry processing
#' so results can be cached, inspected, and tested independently.
#'
#' @param compounds A character vector of chemical names to resolve on PubChem.
#' @param profile One of `"minimal"`, `"ms"`, `"safety"`, `"bioactivity"`, or
#' `"full"`. Profiles control which PubChem properties and PUG-View headings
#' are requested.
#' @param sections Optional PUG-View headings to request in addition to the
#' selected profile.
#' @param sources Optional PUG-View source names to request in addition to the
#' selected profile. Use this for source-specific records such as LOTUS, FEMA,
#' FDA/SPL, or MeSH.
#' @param cache Logical. If `TRUE`, raw PubChem JSON responses are cached under
#' `cache_dir`.
#' @param cache_dir Directory for cached PubChem responses. Defaults to a
#' user-cache location from `tools::R_user_dir()` when available, otherwise a
#' temporary directory.
#' @param throttle Seconds to wait between uncached PubChem requests.
#' @param assay_detail_limit Maximum number of unique PubChem BioAssay AIDs per
#' query for which assay-description metadata is fetched when `profile` is
#' `"bioactivity"` or `"full"`. Active, numeric, and target-bearing assays are
#' prioritized. Set to `0` to skip assay-description requests.
#' @param include_annotations Logical. If `TRUE`, fetch PubChem PUG-View
#' annotation sections selected by `profile`, `sections`, and `sources`. Set to
#' `FALSE` when only identity and property fields are needed.
#' @param query_overrides Optional named character vector or two-column data
#' frame mapping each displayed compound name to a source-backed PubChem query,
#' such as `"cid:2519"` or a full InChIKey. The displayed `Query` remains the
#' compound name. Overrides do not fabricate identity and should come from a
#' reviewed source record.
#' @param request_fun Optional function used to retrieve a URL. This is intended
#' for tests and advanced users. It should return either JSON text or a parsed
#' list.
#'
#' @returns A list with tidy data frames: `identity`, `properties`, `synonyms`,
#' `annotations`, `source_annotations`, `taxonomy`, `classifications`,
#' `spectra`, `safety`, `experimental`, `bioactivity`,
#' `bioassay_details`, and `provenance`. Annotation-derived tables include raw
#' `Value`, normalized `CleanValue`, parsed `ValueNumeric`, and normalized
#' `UnitClean` columns when those fields can be extracted. The `identity` table
#' includes the exact `QueriedName` sent to PubChem and marks conservative
#' deterministic alias matches, such as Greek-letter transliterations and
#' trailing-punctuation cleanup, with `MatchStatus = "resolved_alias"`. The
#' `bioactivity`
#' table stores one PubChem BioAssay summary row per assay with activity
#' outcome, assay name/type, target identifiers, and numeric activity values
#' when PubChem reports them. The `bioassay_details` table stores bounded
#' assay-description metadata for selected AIDs. The returned object has class
#' `"uaf_pubchem_profile"`.
#'
#' @examples
#' \dontrun{
#' profile = pubchemProfile(c("methyl salicylate", "octanal"), profile = "ms")
#' profile$properties
#' profile$spectra
#' }
#'
#' @importFrom jsonlite fromJSON
#' @export
pubchemProfile = function(compounds,
                          profile = c("minimal", "ms", "safety",
                                      "bioactivity", "full"),
                          sections = NULL,
                          sources = NULL,
                          cache = TRUE,
                          cache_dir = NULL,
                          throttle = 0.2,
                          assay_detail_limit = 50,
                          include_annotations = TRUE,
                          query_overrides = NULL,
                          request_fun = NULL) {
  profile = match.arg(profile)
  compounds = .uaf_clean_compounds(compounds)
  fetch = .pubchem_fetcher(cache = cache,
                           cache_dir = cache_dir,
                           throttle = throttle,
                           request_fun = request_fun)

  query_overrides = .pubchem_query_override_map(query_overrides, compounds)
  identity = .pubchem_resolve_cids(compounds, fetch, query_overrides)
  cids = identity$CID[!is.na(identity$CID)]
  cid_query = stats::setNames(identity$Query[!is.na(identity$CID)],
                              paste0(identity$CID[!is.na(identity$CID)]))

  properties = .pubchem_fetch_properties(
    cids = cids,
    properties = .pubchem_profile_properties(profile),
    fetch = fetch,
    cid_query = cid_query
  )

  synonyms = .pubchem_fetch_synonyms(cids = cids,
                                     fetch = fetch,
                                     cid_query = cid_query)

  if (isTRUE(include_annotations)) {
    headings = unique(c(.pubchem_profile_headings(profile), sections))
    heading_annotations = .pubchem_fetch_annotations(cids = cids,
                                                     headings = headings,
                                                     fetch = fetch,
                                                     cid_query = cid_query)
    source_names = unique(c(.pubchem_profile_sources(profile), sources))
    source_annotations = .pubchem_fetch_source_annotations(
      cids = cids,
      sources = source_names,
      fetch = fetch,
      cid_query = cid_query
    )
  } else {
    heading_annotations = .pubchem_empty_table(.pubchem_annotation_cols())
    source_annotations = .pubchem_empty_table(.pubchem_annotation_cols())
  }
  taxonomy = .pubchem_fetch_taxonomy_records(source_annotations, fetch)
  classifications = .pubchem_fetch_classification_records(source_annotations,
                                                          fetch)
  annotations = .pubchem_bind_tables(heading_annotations, source_annotations)

  bioactivity = if (profile %in% c("bioactivity", "full")) {
    .pubchem_fetch_bioactivity(cids = cids,
                               fetch = fetch,
                               cid_query = cid_query)
  } else {
    .pubchem_empty_table(.pubchem_bioactivity_cols())
  }
  bioassay_details = if (profile %in% c("bioactivity", "full")) {
    .pubchem_fetch_bioassay_details(bioactivity = bioactivity,
                                    fetch = fetch,
                                    limit = assay_detail_limit)
  } else {
    .pubchem_empty_table(.pubchem_bioassay_detail_cols())
  }

  out = list(
    identity = identity,
    properties = properties,
    synonyms = synonyms,
    annotations = annotations,
    source_annotations = source_annotations,
    taxonomy = taxonomy,
    classifications = classifications,
    spectra = .pubchem_filter_annotations(annotations, "spectra"),
    safety = .pubchem_filter_annotations(annotations, "safety"),
    experimental = .pubchem_filter_annotations(annotations, "experimental"),
    bioactivity = bioactivity,
    bioassay_details = bioassay_details,
    provenance = .pubchem_profile_provenance(identity = identity,
                                             properties = properties,
                                             synonyms = synonyms,
                                             annotations = annotations,
                                             source_annotations = source_annotations,
                                             taxonomy = taxonomy,
                                             classifications = classifications,
                                             bioactivity = bioactivity,
                                             bioassay_details = bioassay_details),
    profile = profile,
    retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  class(out) = c("uaf_pubchem_profile", class(out))
  out
}

#' @export
print.uaf_pubchem_profile = function(x, ...) {
  cat("uaf PubChem profile\n")
  cat("  profile: ", x$profile, "\n", sep = "")
  cat("  queries: ", nrow(x$identity), "\n", sep = "")
  cat("  resolved CIDs: ", sum(!is.na(x$identity$CID)), "\n", sep = "")
  cat("  properties: ", nrow(x$properties), " rows\n", sep = "")
  cat("  synonyms: ", nrow(x$synonyms), " rows\n", sep = "")
  cat("  annotations: ", nrow(x$annotations), " rows\n", sep = "")
  if (is.data.frame(x$taxonomy)) {
    cat("  taxonomy: ", nrow(x$taxonomy), " rows\n", sep = "")
  }
  if (is.data.frame(x$classifications)) {
    cat("  classifications: ", nrow(x$classifications), " rows\n", sep = "")
  }
  if (is.data.frame(x$bioactivity)) {
    cat("  bioactivity: ", nrow(x$bioactivity), " rows\n", sep = "")
  }
  if (is.data.frame(x$bioassay_details)) {
    cat("  bioassay details: ", nrow(x$bioassay_details), " rows\n", sep = "")
  }
  invisible(x)
}

.uaf_clean_compounds = function(compounds) {
  if (missing(compounds) || is.null(compounds)) {
    stop("`compounds` must be a character vector of chemical names.",
         call. = FALSE)
  }
  compounds = unique(trimws(as.character(compounds)))
  compounds = compounds[!is.na(compounds) & compounds != ""]
  if (length(compounds) < 1) {
    stop("`compounds` must contain at least one non-empty chemical name.",
         call. = FALSE)
  }
  compounds
}

.pubchem_base_url = function() {
  "https://pubchem.ncbi.nlm.nih.gov/rest"
}

.pubchem_fetcher = function(cache, cache_dir, throttle, request_fun,
                            service_busy_limit = Inf,
                            event_fun = NULL,
                            max_attempts = NULL) {
  if (is.null(cache_dir)) {
    cache_dir = .pubchem_default_cache_dir()
  }
  force(cache)
  force(cache_dir)
  force(throttle)
  force(request_fun)
  force(service_busy_limit)
  force(event_fun)
  force(max_attempts)

  service_busy_limit = suppressWarnings(as.numeric(service_busy_limit))
  if (length(service_busy_limit) != 1L || is.na(service_busy_limit) ||
      service_busy_limit < 1) {
    stop("`service_busy_limit` must be a positive number or Inf.",
         call. = FALSE)
  }
  state = new.env(parent = emptyenv())
  state$consecutive_service_busy = 0L

  notify = function(event, url, status_code = NA_integer_, attempt = NA_integer_,
                    attempts = NA_integer_, cache_hit = FALSE,
                    detail = NA_character_) {
    if (!is.function(event_fun)) return(invisible(FALSE))
    payload = list(
      phase = "pubchem_request",
      event = event,
      event_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      url = url,
      status_code = status_code,
      attempt = attempt,
      attempts = attempts,
      cache_hit = isTRUE(cache_hit),
      consecutive_service_busy = state$consecutive_service_busy,
      detail = detail
    )
    tryCatch(event_fun(payload), error = function(error) {
      warning("PubChem progress callback failed: ", conditionMessage(error),
              call. = FALSE)
      invisible(FALSE)
    })
    invisible(TRUE)
  }

  function(url) {
    cache_stem = file.path(cache_dir, .pubchem_url_hash(url))
    cache_file = paste0(cache_stem, ".json")
    cache_url_file = paste0(cache_stem, ".url")
    cache_bound_to_request = isTRUE(cache) && file.exists(cache_file) &&
      file.exists(cache_url_file) && identical(
        paste(readLines(cache_url_file, warn = FALSE, encoding = "UTF-8"),
              collapse = "\n"),
        enc2utf8(url)
      )
    if (cache_bound_to_request) {
      txt = paste(readLines(cache_file, warn = FALSE, encoding = "UTF-8"),
                  collapse = "\n")
      state$consecutive_service_busy = 0L
      notify("cache_hit", url, status_code = 200L, cache_hit = TRUE)
      return(jsonlite::fromJSON(txt, simplifyVector = FALSE))
    }

    attempts = if (!is.null(max_attempts)) {
      suppressWarnings(as.integer(max_attempts))
    } else if (is.null(request_fun)) {
      .pubchem_env_integer("UAFR_PUBCHEM_MAX_ATTEMPTS", 8L)
    } else {
      1L
    }
    if (length(attempts) != 1L || is.na(attempts) || attempts < 1L) {
      stop("`max_attempts` must be a positive integer when supplied.",
           call. = FALSE)
    }
    effective_throttle = if (is.null(request_fun)) {
      max(throttle, .pubchem_env_number("UAFR_PUBCHEM_MIN_DELAY", 1))
    } else {
      throttle
    }
    last_error = NA_character_
    last_status = NA_integer_
    for (attempt in seq_len(attempts)) {
      if (is.null(request_fun)) {
        .pubchem_rate_wait(effective_throttle)
      }
      notify("request_attempt", url, attempt = attempt, attempts = attempts)
      fetched = tryCatch({
        if (is.null(request_fun)) {
          .pubchem_http_get(url)
        } else {
          request_result = request_fun(url)
          if (is.character(request_result)) {
            list(ok = TRUE,
                 status_code = 200L,
                 text = paste(request_result, collapse = "\n"),
                 headers = list())
          } else {
            list(ok = TRUE,
                 status_code = 200L,
                 parsed_direct = request_result,
                 headers = list())
          }
        }
      }, error = function(error) {
        last_error <<- conditionMessage(error)
        last_status <<- .pubchem_status_from_error(last_error)
        list(ok = FALSE,
             status_code = last_status,
             text = NA_character_,
             headers = list(),
             error = last_error)
      })

      if (isTRUE(fetched$ok)) {
        if (!is.null(fetched$parsed_direct)) {
          state$consecutive_service_busy = 0L
          notify("request_success", url, status_code = fetched$status_code,
                 attempt = attempt, attempts = attempts)
          return(fetched$parsed_direct)
        }
        parsed = tryCatch(
          jsonlite::fromJSON(fetched$text, simplifyVector = FALSE),
          error = function(error) {
            last_error <<- conditionMessage(error)
            NULL
          }
        )
        if (!is.null(parsed)) {
          if (isTRUE(cache)) {
            dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
            .pubchem_atomic_write_text(fetched$text, cache_file)
            .pubchem_atomic_write_text(enc2utf8(url), cache_url_file)
          }
          if (is.null(request_fun)) {
            .pubchem_rate_record(fetched$headers, effective_throttle)
          }
          state$consecutive_service_busy = 0L
          notify("request_success", url, status_code = fetched$status_code,
                 attempt = attempt, attempts = attempts)
          return(parsed)
        }
        state$consecutive_service_busy = 0L
        last_status = fetched$status_code
      } else {
        if (!is.null(fetched$error) && length(fetched$error) > 0 &&
            !all(is.na(fetched$error))) {
          last_error = fetched$error
        }
        last_status = fetched$status_code
      }

      service_busy = !is.na(last_status) &&
        as.integer(last_status) %in% c(429L, 503L)
      if (service_busy) {
        state$consecutive_service_busy =
          state$consecutive_service_busy + 1L
        notify("service_busy_response", url, status_code = last_status,
               attempt = attempt, attempts = attempts,
               detail = "PubChem returned HTTP 429 or 503.")
        if (is.finite(service_busy_limit) &&
            state$consecutive_service_busy >= service_busy_limit) {
          delay = .pubchem_retry_delay(
            attempt = attempt,
            status_code = last_status,
            headers = fetched$headers,
            throttle = effective_throttle
          )
          .pubchem_set_next_request_time(delay)
          notify("service_busy_limit_reached", url,
                 status_code = last_status, attempt = attempt,
                 attempts = attempts,
                 detail = paste("Pause before retrying for at least",
                                round(delay, 1), "seconds."))
          condition = structure(
            list(
              message = paste0(
                "PubChem service-busy circuit breaker opened after ",
                state$consecutive_service_busy,
                " consecutive HTTP 429/503 responses. Cached successes are ",
                "preserved; pause and resume later with the same cache."
              ),
              call = NULL,
              status_code = as.integer(last_status),
              url = url,
              consecutive_service_busy = state$consecutive_service_busy,
              retry_after_seconds = delay
            ),
            class = c("uaf_pubchem_service_busy", "error", "condition")
          )
          stop(condition)
        }
      } else {
        state$consecutive_service_busy = 0L
      }

      if (!.pubchem_should_retry(last_status) || attempt >= attempts) {
        break
      }
      delay = .pubchem_retry_delay(
        attempt = attempt,
        status_code = last_status,
        headers = fetched$headers,
        throttle = effective_throttle
      )
      .pubchem_log_retry(url = url,
                         attempt = attempt,
                         attempts = attempts,
                         status_code = last_status,
                         delay = delay)
      .pubchem_set_next_request_time(delay)
      if (is.null(request_fun)) Sys.sleep(delay)
    }

    if (!is.na(last_status) && identical(as.integer(last_status), 404L)) {
      state$consecutive_service_busy = 0L
      notify("request_no_hit", url, status_code = 404L,
             attempt = attempts, attempts = attempts)
      return(NULL)
    }
    notify("request_exhausted", url, status_code = last_status,
           attempt = attempts, attempts = attempts, detail = last_error)
    if (!is.na(last_status) && .pubchem_should_retry(last_status)) {
      warning("PubChem request failed after ", attempts, " attempt(s)",
              " with HTTP status ", last_status, ". This usually means ",
              "PubChem is rate-limiting or temporarily busy; rerun later and ",
              "keep cache enabled.\nURL: ", url, call. = FALSE)
    } else {
      warning("PubChem request failed after ", attempts, " attempt(s): ",
              last_error, "\nURL: ", url, call. = FALSE)
    }
    NULL
  }
}

.pubchem_rate_state = new.env(parent = emptyenv())
.pubchem_rate_state$next_request_time = as.POSIXct(0, origin = "1970-01-01")

.pubchem_http_get = function(url) {
  if (requireNamespace("curl", quietly = TRUE)) {
    handle = curl::new_handle(
      useragent = .pubchem_user_agent(),
      timeout = .pubchem_env_number("UAFR_PUBCHEM_TIMEOUT", 60),
      connecttimeout = .pubchem_env_number("UAFR_PUBCHEM_CONNECT_TIMEOUT", 20)
    )
    curl::handle_setheaders(handle,
                            Accept = "application/json",
                            `Accept-Encoding` = "gzip, deflate")
    response = curl::curl_fetch_memory(url, handle = handle)
    headers = curl::parse_headers_list(response$headers)
    text = rawToChar(response$content)
    return(list(ok = response$status_code >= 200 &&
                  response$status_code < 300,
                status_code = response$status_code,
                text = text,
                headers = headers,
                error = if (response$status_code >= 200 &&
                            response$status_code < 300) {
                  NA_character_
                } else {
                  paste("HTTP status", response$status_code)
                }))
  }

  con = base::url(url, open = "rb")
  on.exit(close(con), add = TRUE)
  text = paste(readLines(con, warn = FALSE, encoding = "UTF-8"),
               collapse = "\n")
  list(ok = TRUE,
       status_code = 200L,
       text = text,
       headers = list(),
       error = NA_character_)
}

.pubchem_user_agent = function() {
  version = tryCatch(as.character(utils::packageVersion("uafR")),
                     error = function(error) "development")
  Sys.getenv(
    "UAFR_PUBCHEM_USER_AGENT",
    paste0("uafR/", version,
           " (https://github.com/castrattonDSU/uafR)")
  )
}

.pubchem_rate_wait = function(throttle) {
  now = Sys.time()
  next_time = .pubchem_rate_state$next_request_time
  if (!inherits(next_time, "POSIXct")) {
    next_time = as.POSIXct(0, origin = "1970-01-01")
  }
  wait = as.numeric(difftime(next_time, now, units = "secs"))
  if (is.finite(wait) && wait > 0) Sys.sleep(wait)
  throttle = max(throttle, 0, na.rm = TRUE)
  .pubchem_rate_state$next_request_time = Sys.time() + throttle
}

.pubchem_rate_record = function(headers, throttle) {
  control = .pubchem_header_value(headers, "x-throttling-control")
  delay = throttle
  if (!is.na(control)) {
    control_lower = tolower(control)
    if (grepl("black", control_lower, fixed = TRUE)) {
      delay = max(delay, .pubchem_env_number("UAFR_PUBCHEM_BLACK_DELAY", 300))
    } else if (grepl("red", control_lower, fixed = TRUE)) {
      delay = max(delay, .pubchem_env_number("UAFR_PUBCHEM_RED_DELAY", 60))
    } else if (grepl("yellow", control_lower, fixed = TRUE)) {
      delay = max(delay, .pubchem_env_number("UAFR_PUBCHEM_YELLOW_DELAY", 10))
    }
  }
  .pubchem_set_next_request_time(delay)
}

.pubchem_set_next_request_time = function(delay) {
  delay = max(delay, 0, na.rm = TRUE)
  target = Sys.time() + delay
  current = .pubchem_rate_state$next_request_time
  if (!inherits(current, "POSIXct") || target > current) {
    .pubchem_rate_state$next_request_time = target
  }
}

.pubchem_should_retry = function(status_code) {
  if (length(status_code) != 1 || is.na(status_code)) return(TRUE)
  as.integer(status_code) %in% c(408L, 425L, 429L, 500L, 502L, 503L, 504L)
}

.pubchem_retry_delay = function(attempt, status_code, headers, throttle) {
  retry_after = suppressWarnings(as.numeric(
    .pubchem_header_value(headers, "retry-after")
  ))
  if (length(retry_after) == 1 && is.finite(retry_after) &&
      retry_after > 0) {
    return(retry_after)
  }
  base = if (!is.na(status_code) &&
             as.integer(status_code) %in% c(429L, 503L)) {
    .pubchem_env_number("UAFR_PUBCHEM_BUSY_BACKOFF", 60)
  } else {
    .pubchem_env_number("UAFR_PUBCHEM_BASE_BACKOFF", 10)
  }
  max_delay = .pubchem_env_number("UAFR_PUBCHEM_MAX_BACKOFF", 900)
  delay = min(max_delay, max(base, throttle) * 2^(attempt - 1))
  jitter = stats::runif(1, min = 0, max = min(5, delay * 0.1))
  delay + jitter
}

.pubchem_log_retry = function(url, attempt, attempts, status_code, delay) {
  verbose = Sys.getenv("UAFR_PUBCHEM_VERBOSE", "TRUE")
  if (!tolower(verbose) %in% c("true", "t", "1", "yes", "y")) {
    return(invisible(FALSE))
  }
  status_label = if (length(status_code) == 1 && !is.na(status_code)) {
    paste0("HTTP ", status_code)
  } else {
    "request error"
  }
  message("PubChem retry backoff: ", status_label,
          "; attempt ", attempt, "/", attempts,
          "; waiting ", round(delay, 1), " seconds before retrying ",
          .pubchem_shorten_url(url))
  invisible(TRUE)
}

.pubchem_shorten_url = function(url, width = 120) {
  url = as.character(url)
  if (nchar(url) <= width) return(url)
  paste0(substr(url, 1, width - 3), "...")
}

.pubchem_header_value = function(headers, name) {
  if (!is.list(headers) || length(headers) < 1) return(NA_character_)
  names_lower = tolower(names(headers))
  hit = which(names_lower == tolower(name))
  if (length(hit) < 1) return(NA_character_)
  as.character(headers[[hit[[1]]]])
}

.pubchem_status_from_error = function(message) {
  if (is.na(message)) return(NA_integer_)
  status = sub(".*HTTP status was ['\"]?([0-9]{3}).*", "\\1", message)
  if (identical(status, message)) {
    status = sub(".*HTTP status ([0-9]{3}).*", "\\1", message)
  }
  if (identical(status, message)) return(NA_integer_)
  suppressWarnings(as.integer(status))
}

.pubchem_env_number = function(name, default) {
  value = Sys.getenv(name, unset = NA_character_)
  value = suppressWarnings(as.numeric(value))
  if (length(value) != 1 || !is.finite(value)) return(default)
  value
}

.pubchem_env_integer = function(name, default) {
  value = .pubchem_env_number(name, default)
  if (!is.finite(value)) return(default)
  as.integer(value)
}

.pubchem_default_cache_dir = function() {
  cache_dir = tryCatch(tools::R_user_dir("uafR", "cache"),
                       error = function(error) NULL)
  if (is.null(cache_dir)) {
    cache_dir = file.path(tempdir(), "uafR-pubchem-cache")
  }
  cache_dir
}

.pubchem_url_hash = function(url) {
  if (length(url) != 1L || is.na(url) || !nzchar(url)) {
    stop("PubChem cache URLs must be one non-empty string.", call. = FALSE)
  }
  temp_file = tempfile("uafr_pubchem_url_")
  on.exit(unlink(temp_file, force = TRUE), add = TRUE)
  writeBin(charToRaw(enc2utf8(as.character(url))), temp_file)
  paste0("v2_", unname(tools::md5sum(temp_file)[[1]]))
}

.pubchem_atomic_write_text = function(text, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_file = tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temp_file, force = TRUE), add = TRUE)
  writeLines(enc2utf8(as.character(text)), temp_file, useBytes = TRUE)
  if (!file.rename(temp_file, path)) {
    stop("Could not atomically write PubChem cache file: ", path,
         call. = FALSE)
  }
  invisible(path)
}

.pubchem_encode_path = function(x) {
  utils::URLencode(x, reserved = TRUE)
}

.pubchem_resolve_cids = function(compounds, fetch, query_overrides = NULL) {
  rows = lapply(compounds, function(compound) {
    override = .uaf_first_non_empty_text(query_overrides[[compound]])
    aliases = unique(.uaf_non_empty(c(
      override, .pubchem_query_aliases(compound)
    )))
    cid = NA_character_
    url = NA_character_
    queried_name = NA_character_
    query_endpoint = NA_character_
    for (alias in aliases) {
      cid_alias = .pubchem_cid_alias(alias)
      endpoint = if (!is.na(cid_alias)) {
        "cid"
      } else if (.uaf_is_inchikey(alias)) {
        "inchikey"
      } else {
        "name"
      }
      query_value = if (identical(endpoint, "cid")) cid_alias else alias
      url = paste0(.pubchem_base_url(), "/pug/compound/", endpoint, "/",
                   .pubchem_encode_path(query_value), "/cids/JSON")
      if (identical(endpoint, "cid")) {
        cid = cid_alias
      } else {
        json = fetch(url)
        cid = tryCatch(json$IdentifierList$CID[[1]],
                       error = function(error) NA)
        cid = .uaf_first_non_empty_text(cid)
      }
      queried_name = alias
      query_endpoint = endpoint
      if (!is.na(cid)) break
    }
    used_override = !is.na(override) && identical(queried_name, override)
    match_status = if (is.na(cid)) {
      "not_found"
    } else if (used_override && identical(query_endpoint, "inchikey")) {
      "resolved_source_inchikey"
    } else if (used_override && identical(query_endpoint, "cid")) {
      "resolved_source_cid"
    } else if (identical(query_endpoint, "inchikey")) {
      ifelse(identical(queried_name, compound), "resolved_inchikey",
             "resolved_inchikey_alias")
    } else if (identical(query_endpoint, "cid")) {
      ifelse(identical(queried_name, compound), "resolved_cid",
             "resolved_cid_alias")
    } else if (identical(queried_name, compound)) {
      "resolved"
    } else {
      "resolved_alias"
    }
    data.frame(Query = compound,
               CID = suppressWarnings(as.integer(cid)),
               MatchStatus = match_status,
               QueriedName = queried_name,
               SourceURL = url,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

.pubchem_query_override_map = function(x, compounds) {
  out = stats::setNames(rep(NA_character_, length(compounds)), compounds)
  if (is.null(x)) return(out)
  if (is.data.frame(x)) {
    if (ncol(x) < 2) {
      stop("`query_overrides` data frames require query and override columns.",
           call. = FALSE)
    }
    query_col = intersect(c("Query", "query", "compound_name"), names(x))
    override_col = intersect(c("PubChemQuery", "pubchem_query", "override"),
                             names(x))
    query_col = if (length(query_col) > 0) query_col[[1]] else names(x)[[1]]
    override_col = if (length(override_col) > 0) override_col[[1]] else
      names(x)[[2]]
    values = as.character(x[[override_col]])
    names(values) = as.character(x[[query_col]])
    x = values
  }
  if (!is.character(x) || is.null(names(x))) {
    stop("`query_overrides` must be a named character vector or data frame.",
         call. = FALSE)
  }
  x_names = names(x)
  x = .uaf_squish_text(x)
  names(x) = .uaf_squish_text(x_names)
  x = x[!is.na(names(x)) & names(x) != "" & !duplicated(names(x))]
  hit = intersect(names(out), names(x))
  out[hit] = x[hit]
  out
}

.pubchem_cid_alias = function(x) {
  x = .uaf_squish_text(x)
  if (is.na(x) || x == "") return(NA_character_)
  x = sub("^cid\\s*[:=]\\s*", "", x, ignore.case = TRUE, perl = TRUE)
  if (!grepl("^[0-9]+$", x)) return(NA_character_)
  x
}

.pubchem_query_aliases = function(compound) {
  compound = .uaf_squish_text(compound)
  if (is.na(compound) || compound == "") return(character())
  aliases = c(compound)
  no_terminal_punctuation = gsub("[,;:.]+$", "", compound, perl = TRUE)
  aliases = c(aliases, no_terminal_punctuation)
  aliases = c(aliases, .pubchem_transliterate_name(no_terminal_punctuation))
  aliases = c(aliases, .pubchem_transliterate_name(compound))
  aliases = c(aliases, gsub("[<>]", "", aliases, perl = TRUE))
  unique(.uaf_non_empty(.uaf_squish_text(aliases)))
}

.pubchem_transliterate_name = function(x) {
  x = .uaf_squish_text(x)
  plus_minus = intToUtf8(0x00b1)
  x = gsub(plus_minus, "+/-", x, fixed = TRUE)
  greek = stats::setNames(
    c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta",
      "theta", "lambda", "mu", "pi", "sigma", "omega"),
    vapply(c(0x03b1, 0x03b2, 0x03b3, 0x03b4, 0x03b5, 0x03b6,
             0x03b7, 0x03b8, 0x03bb, 0x03bc, 0x03c0, 0x03c3,
             0x03c9), intToUtf8, character(1))
  )
  for (pattern in names(greek)) {
    x = gsub(pattern, greek[[pattern]], x, fixed = TRUE)
  }
  x
}

.pubchem_profile_properties = function(profile) {
  minimal = c("Title", "MolecularFormula", "MolecularWeight",
              "IUPACName", "InChI", "InChIKey", "CanonicalSMILES",
              "IsomericSMILES", "SMILES", "ConnectivitySMILES",
              "ExactMass", "MonoisotopicMass")
  descriptors = c("XLogP", "TPSA", "Complexity", "Charge",
                  "HBondDonorCount", "HBondAcceptorCount",
                  "RotatableBondCount", "HeavyAtomCount",
                  "IsotopeAtomCount", "AtomStereoCount",
                  "DefinedAtomStereoCount", "UndefinedAtomStereoCount",
                  "BondStereoCount", "DefinedBondStereoCount",
                  "UndefinedBondStereoCount", "CovalentUnitCount")
  three_d = c("Volume3D", "XStericQuadrupole3D", "YStericQuadrupole3D",
              "ZStericQuadrupole3D", "FeatureCount3D",
              "FeatureAcceptorCount3D", "FeatureDonorCount3D",
              "FeatureAnionCount3D", "FeatureCationCount3D",
              "FeatureRingCount3D", "FeatureHydrophobeCount3D",
              "ConformerModelRMSD3D", "EffectiveRotorCount3D",
              "ConformerCount3D")

  switch(profile,
         minimal = minimal,
         ms = unique(c(minimal, descriptors)),
         safety = unique(c(minimal, descriptors)),
         bioactivity = unique(c(minimal, descriptors)),
         full = unique(c(minimal, descriptors, three_d)))
}

.pubchem_profile_headings = function(profile) {
  switch(profile,
         minimal = character(),
         ms = c("Mass Spectrometry", "GC-MS", "MS-MS"),
         safety = c("Safety and Hazards", "Experimental Properties"),
         bioactivity = c("Pharmacology and Biochemistry", "Drug and Medication Information"),
         full = c("Mass Spectrometry", "GC-MS", "MS-MS",
                  "Safety and Hazards", "Experimental Properties",
                  "Pharmacology and Biochemistry",
                  "Drug and Medication Information",
                  "Names and Identifiers", "Chemical and Physical Properties",
                  "Literature", "Patents"))
}

.pubchem_profile_sources = function(profile) {
  switch(profile,
         minimal = character(),
         ms = character(),
         safety = character(),
         bioactivity = character(),
         full = c("LOTUS - the natural products occurrence database",
                  "Flavor and Extract Manufacturers Association (FEMA)",
                  "FDA/SPL Indexing Data",
                  "Medical Subject Headings (MeSH)"))
}

.pubchem_fetch_properties = function(cids, properties, fetch, cid_query) {
  cols = c("Query", "CID", properties, "SourceURL")
  if (length(cids) < 1) return(.pubchem_empty_table(cols))

  chunks = split(cids, ceiling(seq_along(cids) / 100))
  rows = list()
  for (chunk in chunks) {
    fetched = .pubchem_fetch_property_chunk(chunk = chunk,
                                            properties = properties,
                                            fetch = fetch)
    if (length(fetched) < 1) next
    for (item in fetched) {
      record = item$record
      row = as.list(stats::setNames(rep(NA_character_, length(cols)), cols))
      row$CID = suppressWarnings(as.integer(record$CID))
      row$Query = unname(cid_query[paste0(row$CID)])
      row$SourceURL = item$url
      for (property in properties) {
        row[[property]] = .pubchem_scalar(record[[property]])
      }
      rows[[length(rows) + 1]] = as.data.frame(row, stringsAsFactors = FALSE)
    }
  }
  if (length(rows) < 1) return(.pubchem_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.pubchem_fetch_property_chunk = function(chunk, properties, fetch) {
  prop_string = paste(properties, collapse = ",")
  fetch_once = function(cids) {
    cid_string = paste(cids, collapse = ",")
    url = paste0(.pubchem_base_url(), "/pug/compound/cid/",
                 cid_string, "/property/", prop_string, "/JSON")
    json = fetch(url)
    records = tryCatch(json$PropertyTable$Properties,
                       error = function(error) NULL)
    if (is.null(records)) return(list())
    lapply(records, function(record) list(record = record, url = url))
  }
  out = fetch_once(chunk)
  if (length(out) > 0 || length(chunk) <= 1) return(out)
  out = unlist(lapply(chunk, fetch_once), recursive = FALSE)
  if (length(out) < 1) list() else out
}

.pubchem_fetch_synonyms = function(cids, fetch, cid_query) {
  cols = c("Query", "CID", "Synonym", "SourceURL")
  if (length(cids) < 1) return(.pubchem_empty_table(cols))

  rows = list()
  chunks = split(cids, ceiling(seq_along(cids) / 100))
  for (chunk in chunks) {
    url = paste0(.pubchem_base_url(), "/pug/compound/cid/",
                 paste(chunk, collapse = ","), "/synonyms/JSON")
    json = fetch(url)
    records = tryCatch(json$InformationList$Information, error = function(error) NULL)
    if (is.null(records)) next
    for (record in records) {
      cid = suppressWarnings(as.integer(record$CID))
      synonyms = unlist(record$Synonym, use.names = FALSE)
      if (length(synonyms) < 1) next
      for (synonym in unique(as.character(synonyms))) {
        rows[[length(rows) + 1]] = data.frame(
          Query = unname(cid_query[paste0(cid)]),
          CID = cid,
          Synonym = synonym,
          SourceURL = url,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) < 1) return(.pubchem_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out = unique(out)
  row.names(out) = NULL
  out
}

.pubchem_fetch_annotations = function(cids, headings, fetch, cid_query) {
  cols = .pubchem_annotation_cols()
  if (length(cids) < 1 || length(headings) < 1) {
    return(.pubchem_empty_table(cols))
  }

  rows = list()
  for (cid in cids) {
    for (heading in headings) {
      url = paste0(.pubchem_base_url(), "/pug_view/data/compound/",
                   cid, "/JSON?heading=",
                   utils::URLencode(heading, reserved = TRUE))
      json = fetch(url)
      parsed = .pubchem_parse_pugview(json = json,
                                      cid = cid,
                                      query = unname(cid_query[paste0(cid)]),
                                      heading = heading,
                                      pubchem_url = url)
      if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
    }
  }
  if (length(rows) < 1) return(.pubchem_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  .pubchem_dedupe_annotations(out)
}

.pubchem_fetch_source_annotations = function(cids, sources, fetch, cid_query) {
  cols = .pubchem_annotation_cols()
  if (length(cids) < 1 || length(sources) < 1) {
    return(.pubchem_empty_table(cols))
  }

  rows = list()
  for (cid in cids) {
    for (source_name in sources) {
      url = paste0(.pubchem_base_url(), "/pug_view/data/compound/",
                   cid, "/JSON?source=",
                   utils::URLencode(source_name, reserved = TRUE))
      json = fetch(url)
      parsed = .pubchem_parse_pugview(
        json = json,
        cid = cid,
        query = unname(cid_query[paste0(cid)]),
        heading = paste0("Source: ", source_name),
        pubchem_url = url
      )
      if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
    }
  }
  if (length(rows) < 1) return(.pubchem_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  .pubchem_dedupe_annotations(out)
}

.pubchem_fetch_taxonomy_records = function(annotations, fetch) {
  cols = c("Query", "CID", "TaxonomyID", "Organism", "CommonName",
           "Rank", "Domain", "Kingdom", "Phylum", "Class", "Order",
           "Family", "Genus", "Species", "Lineage", "SourceAnnotation",
           "Source", "SourceURL", "TaxonomyURL", "PubChemURL")
  links = .pubchem_taxonomy_links(annotations)
  if (nrow(links) < 1) return(.pubchem_empty_table(cols))

  detail_cache = list()
  rows = list()
  for (i in seq_len(nrow(links))) {
    link = links[i, , drop = FALSE]
    taxid_key = as.character(link$TaxonomyID[[1]])
    if (is.null(detail_cache[[taxid_key]])) {
      url = paste0(.pubchem_base_url(), "/pug_view/data/taxonomy/",
                   taxid_key, "/JSON")
      json = fetch(url)
      detail_cache[[taxid_key]] = .pubchem_parse_taxonomy_record(json, url)
    }
    detail = detail_cache[[taxid_key]]
    organism = .uaf_first_non_empty_text(detail$Organism, link$Organism)
    lineage = detail$Lineage
    ranks = .lotus_lineage_ranks(lineage, organism)
    rows[[length(rows) + 1]] = data.frame(
      Query = link$Query,
      CID = link$CID,
      TaxonomyID = taxid_key,
      Organism = organism,
      CommonName = detail$CommonName,
      Rank = detail$Rank,
      Domain = .uaf_first_non_empty_text(detail$Domain, ranks$Domain),
      Kingdom = ranks$Kingdom,
      Phylum = ranks$Phylum,
      Class = ranks$Class,
      Order = ranks$Order,
      Family = ranks$Family,
      Genus = ranks$Genus,
      Species = ranks$Species,
      Lineage = lineage,
      SourceAnnotation = link$SourceAnnotation,
      Source = link$Source,
      SourceURL = link$SourceURL,
      TaxonomyURL = paste0(.pubchem_base_url(), "/pug_view/data/taxonomy/",
                           taxid_key, "/JSON"),
      PubChemURL = link$PubChemURL,
      stringsAsFactors = FALSE
    )
  }

  out = do.call(rbind, rows)
  row.names(out) = NULL
  unique(out)
}

.pubchem_classification_cols = function() {
  c("Query", "CID", "Source", "TreeID", "TreeName", "TreeType",
    "RootHNID", "HNID", "NodeID", "ParentNodeID", "ClassName",
    "ParentClass", "ClassPath", "ClassDepth", "CompoundCount",
    "TaxonomyCount", "DOICount", "PubMedCount", "SourceURL",
    "ClassificationURL", "PubChemURL")
}

.pubchem_fetch_classification_records = function(annotations, fetch) {
  cols = .pubchem_classification_cols()
  links = .pubchem_classification_links(annotations)
  if (nrow(links) < 1) return(.pubchem_empty_table(cols))

  detail_cache = list()
  rows = list()
  for (i in seq_len(nrow(links))) {
    link = links[i, , drop = FALSE]
    hid = as.character(link$TreeID[[1]])
    cid = as.character(link$CID[[1]])
    if (is.na(hid) || hid == "" || is.na(cid) || cid == "") next
    cache_key = paste(hid, cid, sep = "\r")
    if (is.null(detail_cache[[cache_key]])) {
      url = .pubchem_classification_url(hid = hid, cid = cid)
      json = fetch(url)
      detail_cache[[cache_key]] = .pubchem_parse_classification_response(
        json = json,
        link = link,
        url = url
      )
    }
    detail = detail_cache[[cache_key]]
    if (is.data.frame(detail) && nrow(detail) > 0) {
      rows[[length(rows) + 1]] = detail
    }
  }

  if (length(rows) < 1) return(.pubchem_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  unique(out)
}

.pubchem_classification_links = function(annotations) {
  cols = c("Query", "CID", "TreeID", "TreeName", "TreeType", "Source",
           "SourceURL", "PubChemURL")
  if (!is.data.frame(annotations) || nrow(annotations) < 1) {
    return(.pubchem_empty_table(cols))
  }

  rows = list()
  for (i in seq_len(nrow(annotations))) {
    annotation = annotations[i, , drop = FALSE]
    context = paste(annotation$Heading, annotation$HeadingPath,
                    annotation$Name, annotation$CleanValue,
                    annotation$Source, annotation$PubChemURL)
    if (!grepl("classification", context, ignore.case = TRUE)) next
    if (!grepl("LOTUS", context, ignore.case = TRUE)) next
    name = .uaf_first_non_empty_text(annotation$Name)
    if (is.na(name) || !grepl("^HID$", name, ignore.case = TRUE)) next

    hid = .uaf_first_non_empty_text(annotation$ValueNumeric,
                                    annotation$CleanValue,
                                    annotation$Value)
    hid = sub("\\.0$", "", hid)
    if (is.na(hid) || !grepl("^[0-9]+$", hid)) next

    tree_name = .uaf_first_non_empty_text(
      sub("^.*>\\s*", "", annotation$HeadingPath),
      annotation$HeadingPath
    )
    tree_type = if (grepl("biological", tree_name, ignore.case = TRUE)) {
      "biological"
    } else if (grepl("chemical", tree_name, ignore.case = TRUE)) {
      "chemical"
    } else {
      "other"
    }

    rows[[length(rows) + 1]] = data.frame(
      Query = annotation$Query,
      CID = annotation$CID,
      TreeID = hid,
      TreeName = tree_name,
      TreeType = tree_type,
      Source = annotation$Source,
      SourceURL = annotation$SourceURL,
      PubChemURL = annotation$PubChemURL,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) < 1) return(.pubchem_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  unique(out)
}

.pubchem_classification_url = function(hid, cid) {
  paste0("https://pubchem.ncbi.nlm.nih.gov/classification/cgi/",
         "classifications.fcgi?format=json&hid=",
         utils::URLencode(as.character(hid), reserved = TRUE),
         "&search_uid=",
         utils::URLencode(as.character(cid), reserved = TRUE),
         "&search_uid_type=cid&search_type=list")
}

.pubchem_parse_classification_response = function(json, link, url) {
  cols = .pubchem_classification_cols()
  if (is.null(json) || is.null(json$Hierarchies)) {
    return(.pubchem_empty_table(cols))
  }

  hierarchies = .pubchem_classification_as_list(
    json$Hierarchies$Hierarchy,
    object_fields = c("Node", "Information", "SourceName", "SourceID", "HID")
  )
  if (length(hierarchies) < 1) return(.pubchem_empty_table(cols))

  rows = list()
  for (hierarchy in hierarchies) {
    root_info = hierarchy$Information
    nodes = .pubchem_classification_node_rows(hierarchy$Node)
    if (nrow(nodes) < 1) next

    match_index = which(nodes$Match %in% TRUE)
    if (length(match_index) < 1) match_index = 1L
    match_index = match_index[[1]]
    leaf_to_root = .pubchem_classification_leaf_order(nodes, match_index)
    match_index = which(leaf_to_root$Match %in% TRUE)
    if (length(match_index) < 1) match_index = 1L
    match_index = match_index[[1]]

    matched = leaf_to_root[match_index, , drop = FALSE]
    parent = if (match_index < nrow(leaf_to_root)) {
      leaf_to_root[match_index + 1L, , drop = FALSE]
    } else {
      NULL
    }
    path_terms = rev(.uaf_non_empty(leaf_to_root$ClassName))
    class_path = .pubchem_collapse(path_terms)

    rows[[length(rows) + 1]] = data.frame(
      Query = link$Query,
      CID = link$CID,
      Source = .uaf_first_non_empty_text(link$Source,
                                         hierarchy$SourceName,
                                         hierarchy$SourceID,
                                         "PubChem Classification"),
      TreeID = .uaf_first_non_empty_text(link$TreeID, hierarchy$HID,
                                         root_info$HID),
      TreeName = .uaf_first_non_empty_text(link$TreeName,
                                           .pubchem_classification_info_name(root_info)),
      TreeType = link$TreeType,
      RootHNID = .uaf_first_non_empty_text(root_info$HNID),
      HNID = matched$HNID,
      NodeID = matched$NodeID,
      ParentNodeID = matched$ParentNodeID,
      ClassName = matched$ClassName,
      ParentClass = if (is.null(parent)) NA_character_ else parent$ClassName,
      ClassPath = class_path,
      ClassDepth = length(path_terms),
      CompoundCount = matched$CompoundCount,
      TaxonomyCount = matched$TaxonomyCount,
      DOICount = matched$DOICount,
      PubMedCount = matched$PubMedCount,
      SourceURL = link$SourceURL,
      ClassificationURL = url,
      PubChemURL = link$PubChemURL,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) < 1) return(.pubchem_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  unique(out)
}

.pubchem_classification_node_rows = function(nodes) {
  cols = c("NodeID", "ParentNodeID", "HNID", "ClassName", "Match",
           "CompoundCount", "TaxonomyCount", "DOICount", "PubMedCount")
  node_list = .pubchem_classification_as_list(
    nodes,
    object_fields = c("NodeID", "ParentID", "ParentNodeID", "Information")
  )
  if (length(node_list) < 1) return(.pubchem_empty_table(cols))

  rows = lapply(node_list, function(node) {
    info = node$Information
    data.frame(
      NodeID = .uaf_first_non_empty_text(node$NodeID, info$NodeID),
      ParentNodeID = .uaf_first_non_empty_text(node$ParentID,
                                               node$ParentNodeID,
                                               info$ParentID,
                                               info$ParentNodeID),
      HNID = .uaf_first_non_empty_text(info$HNID),
      ClassName = .pubchem_classification_info_name(info),
      Match = .pubchem_truthy(info$Match),
      CompoundCount = .pubchem_classification_count(info, "Compound"),
      TaxonomyCount = .pubchem_classification_count(info, "Taxonomy"),
      DOICount = .pubchem_classification_count(info, "DOI"),
      PubMedCount = .pubchem_classification_count(info, "PubMed"),
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.pubchem_classification_leaf_order = function(nodes, match_index) {
  if (nrow(nodes) < 2) return(nodes)

  matched_first = match_index == 1L
  matched_last = match_index == nrow(nodes)
  if (matched_last && !matched_first) {
    return(nodes[rev(seq_len(nrow(nodes))), , drop = FALSE])
  }

  first_parent = nodes$ParentNodeID[[1]]
  second_node = nodes$NodeID[[2]]
  second_parent = nodes$ParentNodeID[[2]]
  first_node = nodes$NodeID[[1]]
  if (!is.na(second_parent) && !is.na(first_node) &&
      second_parent == first_node) {
    return(nodes[rev(seq_len(nrow(nodes))), , drop = FALSE])
  }
  if (!is.na(first_parent) && !is.na(second_node) &&
      first_parent == second_node) {
    return(nodes)
  }
  nodes
}

.pubchem_classification_as_list = function(x, object_fields = character()) {
  if (is.null(x)) return(list())
  if (is.data.frame(x)) {
    return(lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE])))
  }
  if (!is.list(x)) return(list(x))
  if (length(x) < 1) return(list())
  if (length(intersect(names(x), object_fields)) > 0) return(list(x))
  x
}

.pubchem_classification_info_name = function(info) {
  name = .uaf_first_non_empty_text(
    tryCatch(info$Name$StringWithMarkup[[1]]$String,
             error = function(error) NA_character_),
    tryCatch(info$Name$StringWithMarkup$String,
             error = function(error) NA_character_),
    tryCatch(info$Name$String,
             error = function(error) NA_character_),
    tryCatch(info$Name,
             error = function(error) NA_character_)
  )
  name
}

.pubchem_classification_count = function(info, type) {
  counts = info$Counts
  if (is.null(counts)) return(NA_integer_)
  if (!is.null(counts[[type]])) {
    return(suppressWarnings(as.integer(.pubchem_scalar(counts[[type]]))))
  }
  items = .pubchem_classification_as_list(
    counts,
    object_fields = c("Type", "Count")
  )
  for (item in items) {
    item_type = .pubchem_scalar(item$Type)
    if (!is.na(item_type) && tolower(item_type) == tolower(type)) {
      return(suppressWarnings(as.integer(.pubchem_scalar(item$Count))))
    }
  }
  NA_integer_
}

.pubchem_truthy = function(x) {
  value = tolower(.pubchem_scalar(x))
  !is.na(value) && value %in% c("true", "t", "1", "yes", "y")
}

.pubchem_taxonomy_links = function(annotations) {
  cols = c("Query", "CID", "TaxonomyID", "Organism", "SourceAnnotation",
           "Source", "SourceURL", "PubChemURL")
  if (!is.data.frame(annotations) || nrow(annotations) < 1 ||
      !"MarkupURL" %in% colnames(annotations)) {
    return(.pubchem_empty_table(cols))
  }

  rows = list()
  for (i in seq_len(nrow(annotations))) {
    annotation = annotations[i, , drop = FALSE]
    context = paste(annotation$Heading, annotation$HeadingPath,
                    annotation$Source, annotation$PubChemURL)
    if (!grepl("LOTUS", context, ignore.case = TRUE)) next

    urls = .pubchem_split_collapsed(annotation$MarkupURL)
    labels = .pubchem_split_collapsed(annotation$MarkupText)
    if (length(urls) < 1) next
    for (j in seq_along(urls)) {
      taxid = .pubchem_taxonomy_id_from_url(urls[[j]])
      if (is.na(taxid)) next
      label = if (length(labels) >= j) labels[[j]] else NA_character_
      if (.lotus_is_junk_taxon(label)) label = NA_character_
      rows[[length(rows) + 1]] = data.frame(
        Query = annotation$Query,
        CID = annotation$CID,
        TaxonomyID = taxid,
        Organism = label,
        SourceAnnotation = annotation$CleanValue,
        Source = annotation$Source,
        SourceURL = annotation$SourceURL,
        PubChemURL = annotation$PubChemURL,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) < 1) return(.pubchem_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  unique(out)
}

.pubchem_parse_taxonomy_record = function(json, source_url) {
  empty = list(Organism = NA_character_, CommonName = NA_character_,
               Rank = NA_character_, Domain = NA_character_,
               Lineage = NA_character_)
  record_number = tryCatch(json$Record$RecordNumber, error = function(error) NA)
  title = tryCatch(json$Record$RecordTitle, error = function(error) NA)
  parsed = .pubchem_parse_pugview(json = json,
                                  cid = record_number,
                                  query = title,
                                  heading = "Taxonomy",
                                  pubchem_url = source_url)
  if (!is.data.frame(parsed) || nrow(parsed) < 1) return(empty)

  list(
    Organism = .pubchem_taxonomy_field(parsed, "Scientific Name"),
    CommonName = .pubchem_taxonomy_field(parsed, "Common Name"),
    Rank = .pubchem_taxonomy_field(parsed, "Rank"),
    Domain = .pubchem_taxonomy_field(parsed, "Domain"),
    Lineage = .pubchem_taxonomy_field(parsed, "Lineage")
  )
}

.pubchem_taxonomy_field = function(parsed, field) {
  keep = grepl(paste0("(^| > )", .lotus_regex_escape(field), "$"),
               parsed$HeadingPath, ignore.case = TRUE, perl = TRUE)
  keep[is.na(keep)] = FALSE
  .uaf_first_non_empty_text(parsed$CleanValue[keep])
}

.pubchem_taxonomy_id_from_url = function(url) {
  url = .uaf_squish_text(url)
  if (is.na(url) || url == "") return(NA_character_)
  match = regmatches(url, regexec("/taxonomy/([0-9]+)", url,
                                  perl = TRUE))[[1]]
  if (length(match) < 2) return(NA_character_)
  match[[2]]
}

.pubchem_split_collapsed = function(x) {
  x = .uaf_non_empty(x)
  if (length(x) < 1) return(character())
  .uaf_non_empty(unlist(strsplit(paste(x, collapse = "; "),
                                 "\\s*;\\s*", perl = TRUE),
                        use.names = FALSE))
}

.pubchem_parse_pugview = function(json, cid, query, heading, pubchem_url) {
  cols = .pubchem_annotation_cols()
  if (is.null(json) || is.null(json$Record)) return(.pubchem_empty_table(cols))

  references = .pubchem_reference_map(json$Record$Reference)
  sections = json$Record$Section
  rows = .pubchem_parse_sections(sections = sections,
                                 path = character(),
                                 cid = cid,
                                 query = query,
                                 heading = heading,
                                 references = references,
                                 pubchem_url = pubchem_url)
  if (length(rows) < 1) return(.pubchem_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.pubchem_parse_sections = function(sections, path, cid, query, heading,
                                   references, pubchem_url) {
  if (is.null(sections)) return(list())
  if (!is.list(sections) || (!is.null(sections$TOCHeading) || !is.null(sections$Information))) {
    sections = list(sections)
  }

  rows = list()
  for (section in sections) {
    section_heading = .pubchem_scalar(section$TOCHeading)
    section_path = c(path, section_heading)
    section_path = section_path[!is.na(section_path) & section_path != ""]

    info_rows = .pubchem_parse_information(section$Information,
                                           cid = cid,
                                           query = query,
                                           heading = heading,
                                           heading_path = paste(section_path,
                                                                collapse = " > "),
                                           references = references,
                                           pubchem_url = pubchem_url)
    rows = c(rows, info_rows)

    child_rows = .pubchem_parse_sections(section$Section,
                                         path = section_path,
                                         cid = cid,
                                         query = query,
                                         heading = heading,
                                         references = references,
                                         pubchem_url = pubchem_url)
    rows = c(rows, child_rows)
  }
  rows
}

.pubchem_parse_information = function(information, cid, query, heading,
                                      heading_path, references, pubchem_url) {
  if (is.null(information)) return(list())
  if (!is.list(information) || !is.null(information$Value)) {
    information = list(information)
  }

  rows = list()
  for (info in information) {
    name = .pubchem_scalar(info$Name)
    ref_number = .pubchem_scalar(info$ReferenceNumber)
    reference = references[[paste0(ref_number)]]
    values = .pubchem_value_rows(info$Value)
    if (length(values) < 1) {
      values = list(list(
        value = NA_character_,
        unit = NA_character_,
        markup_text = NA_character_,
        markup_url = NA_character_,
        markup_extra = NA_character_
      ))
    }

    for (value in values) {
      value_text = .pubchem_scalar(value$value)
      value_unit = .pubchem_scalar(value$unit)
      measurement = .uaf_parse_measurement(value_text)
      unit_clean = .uaf_first_non_empty_text(value_unit,
                                             measurement$UnitClean[[1]])
      rows[[length(rows) + 1]] = data.frame(
        Query = query,
        CID = suppressWarnings(as.integer(cid)),
        Heading = heading,
        HeadingPath = heading_path,
        Name = name,
        Value = value_text,
        CleanValue = .uaf_squish_text(value_text),
        ValueNumeric = measurement$ValueNumeric[[1]],
        Unit = value_unit,
        UnitClean = unit_clean,
        MarkupText = .pubchem_scalar(value$markup_text),
        MarkupURL = .pubchem_scalar(value$markup_url),
        MarkupExtra = .pubchem_scalar(value$markup_extra),
        Source = .pubchem_reference_field(reference, "Source"),
        SourceURL = .pubchem_reference_field(reference, "URL"),
        PubChemURL = pubchem_url,
        stringsAsFactors = FALSE
      )
    }
  }
  rows
}

.pubchem_reference_map = function(references) {
  out = list()
  if (is.null(references)) return(out)
  if (!is.list(references) || !is.null(references$ReferenceNumber)) {
    references = list(references)
  }
  for (reference in references) {
    number = .pubchem_scalar(reference$ReferenceNumber)
    if (is.na(number)) next
    out[[paste0(number)]] = list(
      Source = .pubchem_first_non_empty(reference$SourceName,
                                        reference$Name,
                                        reference$SourceID),
      URL = .pubchem_first_non_empty(reference$URL,
                                     reference$SourceURL)
    )
  }
  out
}

.pubchem_value_rows = function(value) {
  if (is.null(value)) return(list())
  unit = .pubchem_scalar(value$Unit)

  if (!is.null(value$StringWithMarkup)) {
    return(.pubchem_string_with_markup_rows(value$StringWithMarkup, unit))
  }

  if (!is.null(value$String)) {
    strings = unlist(value$String, use.names = FALSE)
    return(lapply(as.character(strings), function(x) {
      list(value = x, unit = unit, markup_text = NA_character_,
           markup_url = NA_character_, markup_extra = NA_character_)
    }))
  }

  if (!is.null(value$Number)) {
    numbers = unlist(value$Number, use.names = FALSE)
    return(lapply(as.character(numbers), function(x) {
      list(value = x, unit = unit, markup_text = NA_character_,
           markup_url = NA_character_, markup_extra = NA_character_)
    }))
  }

  flattened = unlist(value, use.names = FALSE)
  flattened = flattened[!is.na(flattened)]
  if (length(flattened) < 1) return(list())
  list(list(value = paste(unique(as.character(flattened)), collapse = " | "),
            unit = unit, markup_text = NA_character_,
            markup_url = NA_character_, markup_extra = NA_character_))
}

.pubchem_string_with_markup = function(string_with_markup) {
  rows = .pubchem_string_with_markup_rows(string_with_markup,
                                          unit = NA_character_)
  vapply(rows, function(row) row$value, character(1))
}

.pubchem_string_with_markup_rows = function(string_with_markup, unit) {
  if (is.null(string_with_markup)) return(character())
  if (!is.list(string_with_markup) || !is.null(string_with_markup$String)) {
    string_with_markup = list(string_with_markup)
  }
  rows = lapply(string_with_markup, function(x) {
    value = .pubchem_scalar(x$String)
    markup = .pubchem_markup_summary(value, x$Markup)
    list(value = value,
         unit = unit,
         markup_text = markup$Text,
         markup_url = markup$URL,
         markup_extra = markup$Extra)
  })
  rows[vapply(rows, function(row) {
    !is.na(row$value) && row$value != ""
  }, logical(1))]
}

.pubchem_markup_summary = function(text, markup) {
  empty = list(Text = NA_character_, URL = NA_character_, Extra = NA_character_)
  text = .uaf_squish_text(text)
  if (is.null(markup) || is.na(text) || text == "") return(empty)
  if (!is.list(markup) || !is.null(markup$URL) || !is.null(markup$Start)) {
    markup = list(markup)
  }
  rows = lapply(markup, function(item) {
    start = suppressWarnings(as.integer(.pubchem_scalar(item$Start)))
    length = suppressWarnings(as.integer(.pubchem_scalar(item$Length)))
    label = NA_character_
    if (!is.na(start) && !is.na(length) && length > 0) {
      label = substr(text, start + 1, start + length)
    }
    c(Text = label,
      URL = .pubchem_scalar(item$URL),
      Extra = .pubchem_scalar(item$Extra))
  })
  if (length(rows) < 1) return(empty)
  values = as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  list(Text = .pubchem_collapse(values$Text),
       URL = .pubchem_collapse(values$URL),
       Extra = .pubchem_collapse(values$Extra))
}

.pubchem_filter_annotations = function(annotations, group) {
  cols = colnames(annotations)
  if (nrow(annotations) < 1) return(.pubchem_empty_table(cols))
  pattern = switch(group,
                   spectra = "Mass Spectrometry|GC-MS|MS-MS|Spectral|Spectrum|m/z",
                   safety = "Safety|Hazard|GHS|Toxicity|Fire|Exposure|Handling",
                   experimental = "Experimental|Boiling|Melting|Vapor|Solubility|Density|Flash|Henry|LogP|pKa",
                   group)
  keep = grepl(pattern,
               paste(annotations$Heading, annotations$HeadingPath,
                     annotations$Name, annotations$Value),
               ignore.case = TRUE)
  annotations[keep, , drop = FALSE]
}

.pubchem_fetch_bioactivity = function(cids, fetch, cid_query) {
  cols = .pubchem_bioactivity_cols()
  if (length(cids) < 1) return(.pubchem_empty_table(cols))

  rows = list()
  for (cid in cids) {
    url = paste0(.pubchem_base_url(), "/pug/compound/cid/",
                 cid, "/assaysummary/JSON")
    json = fetch(url)
    parsed = .pubchem_parse_assay_summary(
      json = json,
      cid = cid,
      query = unname(cid_query[paste0(cid)]),
      source_url = url
    )
    if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
  }
  if (length(rows) < 1) return(.pubchem_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  unique(out)
}

.pubchem_bioactivity_cols = function() {
  c("Query", "CID", "Dataset", "AID", "PanelMemberID", "SID",
    "ActivityOutcome", "ActivityClass", "ActivityValue", "ActivityUnit",
    "ActivityName", "AssayName", "AssayType", "TargetAccession",
    "TargetGeneID", "PubMedID", "RNAi", "Key", "Value", "SourceURL")
}

.pubchem_bioassay_detail_cols = function() {
  c("Query", "CID", "AID", "AssayName", "AssaySourceName",
    "AssaySourceID", "AssayDescription", "AssayProtocol", "AssayComment",
    "ActivityOutcomeMethod", "ProjectCategory", "TargetName",
    "TargetDescription", "TargetType", "TargetAccession", "TargetOtherID",
    "TargetGeneID", "TargetTaxonomyID", "TargetOrganism",
    "TargetCommonName", "EndpointNames", "EndpointCount", "ResultNames",
    "PubMedID", "ExternalURL", "AssayDetailURL", "SourceURL")
}

.pubchem_parse_assay_summary = function(json, cid, query, source_url) {
  cols = .pubchem_bioactivity_cols()
  if (is.null(json) || is.null(json$Table) || is.null(json$Table$Row)) {
    return(.pubchem_empty_table(cols))
  }

  column_names = .pubchem_assay_columns(json$Table$Columns)
  row_list = .pubchem_assay_rows(json$Table$Row)
  if (length(row_list) < 1) return(.pubchem_empty_table(cols))

  rows = list()
  for (row in row_list) {
    values = .pubchem_assay_row_values(row, column_names)
    activity_value_col = .pubchem_assay_activity_value_col(names(values))
    activity_value = if (is.na(activity_value_col)) {
      NA_character_
    } else {
      values[[activity_value_col]]
    }
    activity_unit = if (is.na(activity_value_col)) {
      NA_character_
    } else {
      .pubchem_assay_activity_unit(activity_value_col)
    }
    outcome = .pubchem_assay_field(values, c("Activity Outcome",
                                             "ActivityOutcome",
                                             "Outcome"))

    rows[[length(rows) + 1]] = data.frame(
      Query = query,
      CID = suppressWarnings(as.integer(.uaf_first_non_empty_text(
        .pubchem_assay_field(values, "CID"),
        cid
      ))),
      Dataset = "assaysummary",
      AID = suppressWarnings(as.integer(.pubchem_assay_field(values, "AID"))),
      PanelMemberID = .pubchem_assay_field(values, c("Panel Member ID",
                                                     "PanelMemberID")),
      SID = suppressWarnings(as.integer(.pubchem_assay_field(values, "SID"))),
      ActivityOutcome = outcome,
      ActivityClass = .pubchem_activity_class(outcome),
      ActivityValue = suppressWarnings(as.numeric(activity_value)),
      ActivityUnit = activity_unit,
      ActivityName = .pubchem_assay_field(values, c("Activity Name",
                                                    "ActivityName")),
      AssayName = .pubchem_assay_field(values, c("Assay Name", "AssayName")),
      AssayType = .pubchem_assay_field(values, c("Assay Type", "AssayType")),
      TargetAccession = .pubchem_assay_field(values, c("Target Accession",
                                                       "TargetAccession")),
      TargetGeneID = .pubchem_assay_field(values, c("Target GeneID",
                                                    "Target Gene ID",
                                                    "TargetGeneID")),
      PubMedID = .pubchem_assay_field(values, c("PubMed ID", "PubMedID",
                                                "PMID")),
      RNAi = .pubchem_assay_field(values, "RNAi"),
      Key = "Activity Outcome",
      Value = outcome,
      SourceURL = source_url,
      stringsAsFactors = FALSE
    )
  }

  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.pubchem_assay_columns = function(columns) {
  values = tryCatch(columns$Column, error = function(error) NULL)
  values = unlist(values, use.names = FALSE)
  .uaf_non_empty(values)
}

.pubchem_assay_rows = function(rows) {
  if (is.null(rows)) return(list())
  if (!is.list(rows) || !is.null(rows$Cell) || !is.null(rows$AID)) {
    return(list(rows))
  }
  rows
}

.pubchem_assay_row_values = function(row, column_names) {
  if (!is.null(row$Cell)) {
    cells = as.character(unlist(row$Cell, use.names = FALSE))
    if (length(column_names) < length(cells)) {
      column_names = paste0("Column", seq_along(cells))
    }
    length(cells) = length(column_names)
    return(stats::setNames(cells, column_names))
  }

  values = unlist(row, recursive = TRUE, use.names = FALSE)
  values = as.character(values)
  names(values) = names(row)
  values
}

.pubchem_assay_field = function(values, names) {
  if (length(values) < 1) return(NA_character_)
  hit = match(tolower(names), tolower(names(values)))
  hit = hit[!is.na(hit)]
  if (length(hit) < 1) return(NA_character_)
  .uaf_first_non_empty_text(values[[hit[[1]]]])
}

.pubchem_assay_activity_value_col = function(names) {
  hit = grep("^Activity Value", names, ignore.case = TRUE, value = TRUE)
  if (length(hit) < 1) return(NA_character_)
  hit[[1]]
}

.pubchem_assay_activity_unit = function(name) {
  name = .uaf_squish_text(name)
  if (is.na(name) || name == "") return(NA_character_)
  match = regmatches(name, regexec("\\[([^]]+)\\]", name, perl = TRUE))[[1]]
  if (length(match) < 2) return(NA_character_)
  match[[2]]
}

.pubchem_activity_class = function(outcome) {
  outcome = tolower(.uaf_squish_text(outcome))
  if (is.na(outcome) || outcome == "") return("unspecified")
  if (grepl("active", outcome) && !grepl("inactive", outcome)) return("active")
  if (grepl("inactive", outcome)) return("inactive")
  if (grepl("inconclusive|unspecified|probe|not tested", outcome)) {
    return("inconclusive")
  }
  .normalized_clean_label(outcome)
}

.pubchem_fetch_bioassay_details = function(bioactivity, fetch, limit = 50) {
  cols = .pubchem_bioassay_detail_cols()
  if (!is.data.frame(bioactivity) || nrow(bioactivity) < 1 ||
      !"AID" %in% colnames(bioactivity)) {
    return(.pubchem_empty_table(cols))
  }
  limit = suppressWarnings(as.integer(limit))
  if (is.na(limit) || limit < 1) return(.pubchem_empty_table(cols))

  selected = .pubchem_select_bioassay_detail_rows(bioactivity, limit)
  if (nrow(selected) < 1) return(.pubchem_empty_table(cols))

  aids = unique(.uaf_non_empty(selected$AID))
  detail_rows = list()
  for (chunk in split(aids, ceiling(seq_along(aids) / 25))) {
    url = paste0(.pubchem_base_url(), "/pug/assay/aid/",
                 paste(chunk, collapse = ","),
                 "/description/JSON")
    json = fetch(url)
    parsed = .pubchem_parse_bioassay_detail_response(json, url)
    if (nrow(parsed) > 0) detail_rows[[length(detail_rows) + 1]] = parsed
  }
  if (length(detail_rows) < 1) return(.pubchem_empty_table(cols))

  details = do.call(rbind, detail_rows)
  rows = list()
  for (i in seq_len(nrow(selected))) {
    selected_row = selected[i, , drop = FALSE]
    aid_rows = details[details$AID == selected_row$AID[[1]], , drop = FALSE]
    if (nrow(aid_rows) < 1) next
    for (j in seq_len(nrow(aid_rows))) {
      detail = aid_rows[j, , drop = FALSE]
      detail$Query = selected_row$Query
      detail$CID = selected_row$CID
      rows[[length(rows) + 1]] = detail
    }
  }

  if (length(rows) < 1) return(.pubchem_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  unique(out[, cols, drop = FALSE])
}

.pubchem_select_bioassay_detail_rows = function(bioactivity, limit) {
  cols = c("Query", "CID", "AID")
  for (col in setdiff(c(cols, "ActivityClass", "ActivityValue",
                        "TargetAccession", "TargetGeneID"),
                      colnames(bioactivity))) {
    bioactivity[[col]] = NA_character_
  }
  source = bioactivity[!is.na(bioactivity$AID) & bioactivity$AID != "",
                       ,
                       drop = FALSE]
  if (nrow(source) < 1) return(.pubchem_empty_table(cols))

  source$Priority = 5L
  source$Priority[source$ActivityClass == "active"] = 1L
  source$Priority[!is.na(source$ActivityValue) &
                    source$Priority > 2L] = 2L
  has_target = (!is.na(source$TargetAccession) &
                  source$TargetAccession != "") |
    (!is.na(source$TargetGeneID) & source$TargetGeneID != "")
  source$Priority[has_target & source$Priority > 3L] = 3L
  source = source[order(source$Query, source$Priority, source$AID),
                  ,
                  drop = FALSE]

  rows = list()
  query_keys = unique(.pubchem_key_value(source$Query))
  for (query_key in query_keys) {
    query_rows = source[.pubchem_key_value(source$Query) == query_key,
                        ,
                        drop = FALSE]
    key = paste(.pubchem_key_value(query_rows$Query),
                .pubchem_key_value(query_rows$CID),
                .pubchem_key_value(query_rows$AID),
                sep = "\r")
    query_rows = query_rows[!duplicated(key), , drop = FALSE]
    query_rows = utils::head(query_rows, limit)
    rows[[length(rows) + 1]] = query_rows[, cols, drop = FALSE]
  }

  out = do.call(rbind, rows)
  row.names(out) = NULL
  unique(out)
}

.pubchem_parse_bioassay_detail_response = function(json, source_url) {
  cols = .pubchem_bioassay_detail_cols()
  containers = tryCatch(json$PC_AssayContainer, error = function(error) NULL)
  containers = .pubchem_assay_container_list(containers)
  if (length(containers) < 1) return(.pubchem_empty_table(cols))

  rows = list()
  for (container in containers) {
    descr = tryCatch(container$assay$descr, error = function(error) NULL)
    if (is.null(descr)) next
    aid = suppressWarnings(as.integer(.pubchem_scalar(descr$aid$id)))
    source = descr$aid_source$db
    xrefs = .pubchem_assay_xrefs(descr$xref)
    targets = .pubchem_assay_target_list(descr$target)
    if (length(targets) < 1) targets = list(NULL)

    result_names = .pubchem_assay_result_names(descr$results)
    endpoint_names = .pubchem_assay_endpoint_names(descr$results)
    assay_detail_url = if (!is.na(aid)) {
      paste0(.pubchem_base_url(), "/pug/assay/aid/", aid,
             "/description/JSON")
    } else {
      source_url
    }

    for (target in targets) {
      target_info = .pubchem_assay_target_info(target, xrefs)
      rows[[length(rows) + 1]] = data.frame(
        Query = NA_character_,
        CID = NA_integer_,
        AID = aid,
        AssayName = .pubchem_scalar(descr$name),
        AssaySourceName = .pubchem_scalar(source$name),
        AssaySourceID = .uaf_first_non_empty_text(
          .pubchem_scalar(source$source_id$str),
          .pubchem_scalar(source$source_id$id)
        ),
        AssayDescription = .pubchem_collapse(descr$description),
        AssayProtocol = .pubchem_collapse(descr$protocol),
        AssayComment = .pubchem_collapse(descr$comment),
        ActivityOutcomeMethod = .pubchem_scalar(descr$activity_outcome_method),
        ProjectCategory = .pubchem_scalar(descr$project_category),
        TargetName = target_info$TargetName,
        TargetDescription = target_info$TargetDescription,
        TargetType = target_info$TargetType,
        TargetAccession = target_info$TargetAccession,
        TargetOtherID = target_info$TargetOtherID,
        TargetGeneID = target_info$TargetGeneID,
        TargetTaxonomyID = target_info$TargetTaxonomyID,
        TargetOrganism = target_info$TargetOrganism,
        TargetCommonName = target_info$TargetCommonName,
        EndpointNames = endpoint_names,
        EndpointCount = length(.uaf_non_empty(result_names)),
        ResultNames = result_names,
        PubMedID = .pubchem_collapse(xrefs$PubMedID),
        ExternalURL = .pubchem_collapse(xrefs$ExternalURL),
        AssayDetailURL = assay_detail_url,
        SourceURL = source_url,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) < 1) return(.pubchem_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  unique(out[, cols, drop = FALSE])
}

.pubchem_assay_container_list = function(containers) {
  if (is.null(containers)) return(list())
  if (!is.list(containers) || !is.null(containers$assay)) {
    return(list(containers))
  }
  containers
}

.pubchem_assay_target_list = function(targets) {
  if (is.null(targets)) return(list())
  if (!is.list(targets) || !is.null(targets$name) || !is.null(targets$mol_id) ||
      !is.null(targets$tid)) {
    return(list(targets))
  }
  targets
}

.pubchem_assay_xrefs = function(xref) {
  items = .pubchem_assay_xref_list(xref)
  rows = list()
  for (item in items) {
    value = item$xref
    if (is.null(value)) next
    rows[[length(rows) + 1]] = data.frame(
      PubMedID = .pubchem_scalar(value$pmid),
      GeneID = .pubchem_scalar(value$gene),
      TaxonomyID = .pubchem_scalar(value$taxonomy),
      TaxonomyName = if (!is.na(.pubchem_scalar(value$taxonomy))) {
        .pubchem_scalar(item$comment)
      } else {
        NA_character_
      },
      ExternalURL = .uaf_first_non_empty_text(.pubchem_scalar(value$dburl),
                                              .pubchem_scalar(value$asurl)),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) < 1) {
    return(.uaf_empty_table(c("PubMedID", "GeneID", "TaxonomyID",
                              "TaxonomyName", "ExternalURL")))
  }
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.pubchem_assay_xref_list = function(xref) {
  if (is.null(xref)) return(list())
  if (!is.list(xref) || !is.null(xref$xref)) return(list(xref))
  xref
}

.pubchem_assay_target_info = function(target, xrefs) {
  empty = list(TargetName = NA_character_,
               TargetDescription = NA_character_,
               TargetType = NA_character_,
               TargetAccession = NA_character_,
               TargetOtherID = NA_character_,
               TargetGeneID = .uaf_first_non_empty_text(xrefs$GeneID),
               TargetTaxonomyID = .uaf_first_non_empty_text(xrefs$TaxonomyID),
               TargetOrganism = NA_character_,
               TargetCommonName = NA_character_)
  if (is.null(target)) return(empty)

  mol_id = target$mol_id
  org = target$organism$org
  taxonomy = .uaf_first_non_empty_text(
    .pubchem_assay_taxonomy_from_db(org$db),
    empty$TargetTaxonomyID
  )
  accession = .uaf_first_non_empty_text(
    .pubchem_scalar(mol_id$protein_accession),
    .pubchem_scalar(mol_id$nucleotide_accession)
  )
  other_id = .pubchem_scalar(mol_id$other)
  target_type = if (!is.na(accession)) {
    "protein_or_gene"
  } else if (!is.na(other_id) && grepl("cell-line", other_id,
                                       ignore.case = TRUE)) {
    "cell_line"
  } else if (!is.na(.pubchem_scalar(target$name))) {
    "assay_target"
  } else {
    NA_character_
  }

  list(TargetName = .pubchem_scalar(target$name),
       TargetDescription = .pubchem_scalar(target$descr),
       TargetType = target_type,
       TargetAccession = accession,
       TargetOtherID = other_id,
       TargetGeneID = empty$TargetGeneID,
       TargetTaxonomyID = taxonomy,
       TargetOrganism = .uaf_first_non_empty_text(.pubchem_scalar(org$taxname),
                                                  xrefs$TaxonomyName),
       TargetCommonName = .pubchem_scalar(org$common))
}

.pubchem_assay_taxonomy_from_db = function(db) {
  items = .pubchem_assay_target_list(db)
  for (item in items) {
    name = .pubchem_scalar(item$db)
    if (!is.na(name) && name == "taxon") {
      return(.pubchem_scalar(item$tag$id))
    }
  }
  NA_character_
}

.pubchem_assay_result_names = function(results) {
  items = .pubchem_assay_target_list(results)
  .pubchem_collapse(vapply(items, function(item) {
    .pubchem_scalar(item$name)
  }, character(1)))
}

.pubchem_assay_endpoint_names = function(results) {
  items = .pubchem_assay_target_list(results)
  names = vapply(items, function(item) .pubchem_scalar(item$name),
                 character(1))
  keep = grepl("potency|ic50|ec50|ac50|gi50|lc50|ki|kd|efficacy|phenotype|activity outcome|response",
               names, ignore.case = TRUE, perl = TRUE)
  .pubchem_collapse(names[keep])
}

.pubchem_flatten_json = function(json) {
  cols = c("Key", "Value")
  if (is.null(json)) return(.pubchem_empty_table(cols))
  flat = unlist(json, recursive = TRUE, use.names = TRUE)
  if (length(flat) < 1) return(.pubchem_empty_table(cols))
  out = data.frame(Key = names(flat),
                   Value = as.character(flat),
                   stringsAsFactors = FALSE)
  row.names(out) = NULL
  out
}

.pubchem_profile_provenance = function(...) {
  tables = list(...)
  rows = list()
  for (name in names(tables)) {
    table = tables[[name]]
    if (!is.data.frame(table) || !"SourceURL" %in% colnames(table)) next
    urls = unique(table$SourceURL[!is.na(table$SourceURL) & table$SourceURL != ""])
    if (length(urls) < 1) next
    rows[[length(rows) + 1]] = data.frame(
      Table = name,
      SourceURL = urls,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) < 1) return(.pubchem_empty_table(c("Table", "SourceURL")))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.pubchem_bind_tables = function(...) {
  tables = list(...)
  tables = tables[vapply(tables, is.data.frame, logical(1))]
  if (length(tables) < 1) return(data.frame())
  cols = unique(unlist(lapply(tables, colnames), use.names = FALSE))
  aligned = lapply(tables, function(table) {
    missing_cols = setdiff(cols, colnames(table))
    for (col in missing_cols) table[[col]] = NA_character_
    table[, cols, drop = FALSE]
  })
  out = do.call(rbind, aligned)
  row.names(out) = NULL
  .pubchem_dedupe_annotations(out)
}

.pubchem_scalar = function(x) {
  if (is.null(x) || length(x) < 1) return(NA_character_)
  as.character(unlist(x, use.names = FALSE)[1])
}

.pubchem_first_non_empty = function(...) {
  vals = list(...)
  for (val in vals) {
    val = .pubchem_scalar(val)
    if (!is.na(val) && val != "") return(val)
  }
  NA_character_
}

.pubchem_reference_field = function(reference, field) {
  if (is.null(reference)) return(NA_character_)
  .pubchem_scalar(reference[[field]])
}

.pubchem_empty_table = function(cols) {
  out = as.data.frame(stats::setNames(rep(list(character()), length(cols)), cols),
                      stringsAsFactors = FALSE)
  out[0, , drop = FALSE]
}
