.uaf_squish_text = function(x) {
  x = as.character(x)
  x[is.na(x)] = NA_character_
  x = gsub("[\r\n\t]+", " ", x)
  x = gsub("[[:space:]]+", " ", x)
  trimws(x)
}

.uaf_non_empty = function(x) {
  x = .uaf_squish_text(x)
  x[!is.na(x) & x != ""]
}

.uaf_first_non_empty_text = function(...) {
  vals = .uaf_non_empty(unlist(list(...), use.names = FALSE))
  if (length(vals) < 1) return(NA_character_)
  vals[[1]]
}

.uaf_empty_table = function(cols) {
  out = as.data.frame(stats::setNames(rep(list(character()), length(cols)), cols),
                      stringsAsFactors = FALSE)
  out[0, , drop = FALSE]
}

.uaf_extract_pattern = function(text, pattern) {
  text = paste(.uaf_non_empty(text), collapse = " ")
  if (is.na(text) || text == "") return(character())
  matches = gregexpr(pattern, text, perl = TRUE, ignore.case = FALSE)
  out = regmatches(text, matches)[[1]]
  unique(out[out != ""])
}

.uaf_extract_identifiers = function(text) {
  patterns = list(
    CAS = "\\b[0-9]{2,7}-[0-9]{2}-[0-9]\\b",
    KEGG = "\\b[CDGMR]\\d{5}\\b",
    EC = "\\b[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9A-Za-z-]+\\b",
    ChEBI = "\\b(?:CHEBI|ChEBI):\\s*\\d+\\b",
    PubMed = "\\bPMID:?\\s*\\d+\\b|\\bPubMed:?\\s*\\d+\\b",
    DOI = "\\b10\\.\\d{4,9}/[-._;()/:A-Za-z0-9]+\\b",
    UNII = "\\b[A-Z0-9]{10}\\b",
    GHS = "\\bH[0-9]{3}[A-Za-z]?\\b|\\bP[0-9]{3}[A-Za-z]?\\b"
  )

  rows = list()
  for (identifier_type in names(patterns)) {
    values = .uaf_extract_pattern(text, patterns[[identifier_type]])
    if (length(values) < 1) next
    values = gsub("^(CHEBI|ChEBI):\\s*", "CHEBI:", values)
    values = gsub("^(PMID|PubMed):?\\s*", "PMID:", values)
    values = unique(values)
    rows[[length(rows) + 1]] = data.frame(
      IdentifierType = identifier_type,
      Identifier = values,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) < 1) {
    return(.uaf_empty_table(c("IdentifierType", "Identifier")))
  }
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.uaf_parse_measurement = function(x) {
  x = .uaf_squish_text(x)
  pattern = "([-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?)\\s*([^,;|()]*)"
  matches = regexec(pattern, x, perl = TRUE)
  pieces = regmatches(x, matches)

  values = vapply(pieces, function(piece) {
    if (length(piece) < 2) return(NA_real_)
    suppressWarnings(as.numeric(piece[[2]]))
  }, numeric(1))
  units = vapply(pieces, function(piece) {
    if (length(piece) < 3) return(NA_character_)
    unit = .uaf_squish_text(piece[[3]])
    if (is.na(unit) || unit == "") return(NA_character_)
    unit
  }, character(1))

  data.frame(ValueNumeric = values,
             UnitClean = units,
             stringsAsFactors = FALSE)
}

.uaf_parse_numeric = function(x) {
  .uaf_parse_measurement(x)$ValueNumeric
}

.uaf_yes_no = function(x) {
  ifelse(isTRUE(x), "Yes", "No")
}

.uaf_extract_cid = function(cid_result) {
  if (is.null(cid_result) || length(cid_result) < 1) return(NA_character_)

  value = tryCatch({
    if (is.data.frame(cid_result) || is.matrix(cid_result)) {
      if (nrow(cid_result) < 1 || ncol(cid_result) < 2) return(NA_character_)
      cid_result[[1, 2]]
    } else if (is.list(cid_result) && length(cid_result) >= 2) {
      cid_result[[2]][[1]]
    } else {
      cid_result[[1]]
    }
  }, error = function(error) NA_character_)

  value = .uaf_squish_text(value[[1]])
  if (length(value) < 1 ||
      is.na(value) ||
      value %in% c("", "NA", "NaN", "NULL", "Limit Met")) {
    return(NA_character_)
  }
  value
}

.uaf_get_cid = function(query, from = NULL) {
  cid_result = tryCatch({
    if (is.null(from)) {
      webchem::get_cid(query)
    } else {
      webchem::get_cid(query, from = from)
    }
  }, error = function(error) NA_character_)
  .uaf_extract_cid(cid_result)
}
