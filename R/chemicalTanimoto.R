#' Pairwise PubChem fingerprint Tanimoto similarity for chemicals
#'
#' @description
#' `chemicalTanimotoSimilarity()` resolves chemicals to PubChem CIDs, fetches
#' PubChem Fingerprint2D bit strings, decodes them with `ChemmineR`, and
#' computes pairwise Tanimoto similarity. It accepts a plain character vector or
#' a data frame with compound identifiers. When a grouping column is supplied,
#' it also returns group-level summaries that are suitable for joining to other
#' sample, species, treatment, or site similarity tables.
#'
#' @param compounds Character vector or data frame. Character vectors are
#' treated as compound names. Data frames can contain compound names, PubChem
#' CIDs, InChIKeys, SMILES, and optional grouping columns.
#' @param compound_col Column containing compound names when `compounds` is a
#' data frame. If `NULL`, common names such as `compound_name`, `Query`, or
#' `Chemical` are detected.
#' @param cid_col Column containing PubChem CIDs. If `NULL`, common names such
#' as `CID`, `pubchem_cid`, or `cid` are detected.
#' @param inchikey_col Column containing InChIKeys. If `NULL`, common names such
#' as `InChIKey` or `inchikey` are detected.
#' @param smiles_col Optional SMILES column retained as metadata. PubChem
#' Fingerprint2D values are still fetched from PubChem so the Tanimoto
#' definition is consistent across workflows.
#' @param compound_id_col Optional stable compound identifier column. If absent,
#' identifiers are derived from InChIKey, PubChem CID, or compound name.
#' @param group_cols Optional grouping columns. For plants, use `"species"` or
#' call `plantChemicalTanimotoSimilarity()`.
#' @param thresholds Numeric Tanimoto thresholds to count in group summaries.
#' @param top_n_pairs Number of top compound-pair explanations to include in
#' group summaries.
#' @param return_compound_pairs Logical. If `TRUE`, include or write the
#' compound-compound pair table.
#' @param return_group_compound_pairs Logical. If `TRUE` and `group_cols` are
#' supplied, include or write all cross-group compound-pair rows.
#' @param out_dir Optional output directory. When supplied, large pairwise
#' tables are written as compressed CSV files and the result records their
#' paths.
#' @param max_in_memory_pairs Maximum pair rows allowed in memory when
#' `out_dir` is `NULL`.
#' @param pair_block_size Number of focal rows per pairwise computation block.
#' @param cache Logical. If `TRUE`, PubChem JSON responses are cached.
#' @param cache_dir Cache directory. Defaults to the uafR user cache.
#' @param throttle Seconds to wait after uncached PubChem requests.
#' @param refresh Logical. If `TRUE`, ignore cached PubChem responses.
#' @param request_fun Optional request function for tests. It receives a URL and
#' returns parsed JSON or JSON text.
#' @param name_fallback Logical. If `TRUE`, unresolved compounds may be queried
#' by compound name after CID and InChIKey attempts fail.
#'
#' @return A list with class `"uaf_tanimoto_similarity"` containing
#' `CompoundInput`, `CompoundResolution`, `PubChemFingerprints`,
#' `CompoundTanimoto`, `GroupCompoundMembership`,
#' `GroupPairTanimotoSummary`, `GroupCompoundTanimoto`, `ExportManifest`, and
#' `Provenance`.
#'
#' @examples
#' \dontrun{
#' compounds = data.frame(
#'   species = c("Plant A", "Plant A", "Plant B"),
#'   compound_name = c("caffeine", "theobromine", "aspirin")
#' )
#' sim = chemicalTanimotoSimilarity(compounds, group_cols = "species")
#' sim$GroupPairTanimotoSummary
#' }
#'
#' @export
chemicalTanimotoSimilarity = function(compounds,
                                      compound_col = NULL,
                                      cid_col = NULL,
                                      inchikey_col = NULL,
                                      smiles_col = NULL,
                                      compound_id_col = NULL,
                                      group_cols = NULL,
                                      thresholds = c(0.5, 0.7, 0.85, 0.95),
                                      top_n_pairs = 5,
                                      return_compound_pairs = TRUE,
                                      return_group_compound_pairs = FALSE,
                                      out_dir = NULL,
                                      max_in_memory_pairs = 1000000,
                                      pair_block_size = 250,
                                      cache = TRUE,
                                      cache_dir = NULL,
                                      throttle = 0.2,
                                      refresh = FALSE,
                                      request_fun = NULL,
                                      name_fallback = TRUE) {
  thresholds = suppressWarnings(as.numeric(thresholds))
  thresholds = thresholds[!is.na(thresholds)]
  if (length(thresholds) < 1) thresholds = c(0.5, 0.7, 0.85, 0.95)
  input = .tanimoto_normalize_input(
    compounds = compounds,
    compound_col = compound_col,
    cid_col = cid_col,
    inchikey_col = inchikey_col,
    smiles_col = smiles_col,
    compound_id_col = compound_id_col,
    group_cols = group_cols
  )
  compound_input = input$compounds
  membership = input$membership

  if (nrow(compound_input) < 2) {
    stop("At least two unique compounds are required for Tanimoto similarity.",
         call. = FALSE)
  }

  if (is.null(cache_dir)) cache_dir = .pubchem_default_cache_dir()
  if (isTRUE(refresh) && isTRUE(cache)) {
    unlink(file.path(cache_dir, "tanimoto"), recursive = TRUE, force = TRUE)
  }
  fetch = .pubchem_fetcher(cache = cache,
                           cache_dir = file.path(cache_dir, "tanimoto"),
                           throttle = throttle,
                           request_fun = request_fun)

  resolution = .tanimoto_resolve_cids(compound_input, fetch,
                                      name_fallback = name_fallback)
  fingerprints = .tanimoto_fetch_fingerprints(resolution, fetch)
  fingerprinted = fingerprints[
    !is.na(fingerprints$Fingerprint2D) & fingerprints$Fingerprint2D != "",
    , drop = FALSE
  ]
  if (nrow(fingerprinted) < 2) {
    stop("Fewer than two compounds resolved to PubChem Fingerprint2D values.",
         call. = FALSE)
  }

  bits = .tanimoto_decode_fingerprints(fingerprinted)
  fingerprinted = fingerprinted[match(rownames(bits), fingerprinted$compound_id),
                                , drop = FALSE]

  manifest = .tanimoto_empty_manifest()
  compound_pairs = .uaf_empty_table(.tanimoto_compound_pair_cols())
  if (isTRUE(return_compound_pairs)) {
    pair_count = choose(nrow(bits), 2)
    if (!is.null(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      path = file.path(out_dir, "compound_pair_tanimoto.csv.gz")
      .tanimoto_write_compound_pairs(bits, fingerprinted, path,
                                     block_size = pair_block_size)
      manifest = rbind(manifest, .tanimoto_manifest_row(
        "CompoundTanimoto", path, pair_count,
        "Pairwise PubChem Fingerprint2D Tanimoto between unique compounds."
      ))
      compound_pairs = .uaf_empty_table(.tanimoto_compound_pair_cols())
    } else {
      if (pair_count > max_in_memory_pairs) {
        stop("Compound pair table would contain ", pair_count,
             " rows. Supply `out_dir` to stream the compressed CSV output, ",
             "or increase `max_in_memory_pairs` deliberately.",
             call. = FALSE)
      }
      compound_pairs = .tanimoto_compound_pair_table(bits, fingerprinted,
                                                     block_size = pair_block_size)
    }
  }

  group_summary = .uaf_empty_table(.tanimoto_group_summary_cols(group_cols))
  group_pairs = .uaf_empty_table(.tanimoto_group_pair_cols(group_cols))
  if (length(group_cols) > 0 && nrow(membership) > 0) {
    membership = membership[membership$compound_id %in% rownames(bits), ,
                            drop = FALSE]
    group_summary = .tanimoto_group_pair_summary(
      bits = bits,
      membership = membership,
      compounds = fingerprinted,
      group_cols = group_cols,
      thresholds = thresholds,
      top_n_pairs = top_n_pairs
    )
    if (isTRUE(return_group_compound_pairs)) {
      group_pair_count = .tanimoto_group_pair_count(membership)
      if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        path = file.path(out_dir, "group_compound_pair_tanimoto.csv.gz")
        .tanimoto_write_group_compound_pairs(bits, membership, fingerprinted,
                                             path,
                                             block_size = pair_block_size)
        manifest = rbind(manifest, .tanimoto_manifest_row(
          "GroupCompoundTanimoto", path, group_pair_count,
          "Cross-group compound-pair PubChem Fingerprint2D Tanimoto table."
        ))
      } else {
        if (group_pair_count > max_in_memory_pairs) {
          stop("Group compound pair table would contain ", group_pair_count,
               " rows. Supply `out_dir` to stream the compressed CSV output, ",
               "or increase `max_in_memory_pairs` deliberately.",
               call. = FALSE)
        }
        group_pairs = .tanimoto_group_compound_pair_table(
          bits = bits,
          membership = membership,
          compounds = fingerprinted,
          block_size = pair_block_size
        )
      }
    }
  }

  result = list(
    CompoundInput = compound_input,
    CompoundResolution = resolution,
    PubChemFingerprints = fingerprints,
    CompoundTanimoto = compound_pairs,
    GroupCompoundMembership = membership,
    GroupPairTanimotoSummary = group_summary,
    GroupCompoundTanimoto = group_pairs,
    ExportManifest = manifest,
    Provenance = data.frame(
      Function = "chemicalTanimotoSimilarity",
      FingerprintSource = "PubChem Fingerprint2D",
      RetrievedAt = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      CacheEnabled = isTRUE(cache),
      CacheDir = cache_dir,
      stringsAsFactors = FALSE
    )
  )
  class(result) = c("uaf_tanimoto_similarity", "list")
  result
}

#' Plant-labeled PubChem fingerprint Tanimoto similarity
#'
#' @description
#' `plantChemicalTanimotoSimilarity()` converts a plant phytochemistry result
#' from `resolvePlantPhytochemistry()` or `runPlantPhytochemistryBatch()` into a
#' species-compound membership table, then calls
#' `chemicalTanimotoSimilarity()`. The returned tables preserve plant labels,
#' evidence tiers, plant-part context, and source identifiers so pairwise
#' chemistry outputs can be joined directly to downstream ecological,
#' phylogenetic, remediation, or trait analyses.
#'
#' @param plant_chemistry A plant phytochemistry result list containing
#' `PlantCompoundOccurrences` and `CompoundResolution`, or a data frame already
#' containing plant-compound rows.
#' @param level Plant grouping level. Defaults to `"species"`.
#' @param include_review_required Logical. If `TRUE`, carry
#' `CompoundIdentityReview` flags into the membership table when present.
#' @param ... Additional arguments passed to `chemicalTanimotoSimilarity()`,
#' including `out_dir`, `cache_dir`, `return_group_compound_pairs`, and
#' `request_fun`.
#'
#' @return A `"uaf_plant_tanimoto_similarity"` list containing the generic
#' Tanimoto tables plus plant-labeled aliases: `PlantCompoundMembership`,
#' `PlantPairTanimotoSummary`, and `PlantCompoundTanimoto`.
#'
#' @examples
#' \dontrun{
#' phyto = resolvePlantPhytochemistry(c("Salix nigra", "Zea mays"))
#' sim = plantChemicalTanimotoSimilarity(
#'   phyto,
#'   out_dir = "plant_tanimoto",
#'   return_group_compound_pairs = TRUE
#' )
#' sim$PlantPairTanimotoSummary
#' }
#'
#' @export
plantChemicalTanimotoSimilarity = function(plant_chemistry,
                                           level = "species",
                                           include_review_required = TRUE,
                                           ...) {
  table = .tanimoto_plant_membership_table(
    plant_chemistry,
    level = level,
    include_review_required = include_review_required
  )
  if (!level %in% names(table)) {
    stop("Plant chemistry input does not contain grouping level `", level,
         "`.", call. = FALSE)
  }
  result = chemicalTanimotoSimilarity(
    table,
    compound_col = "compound_name",
    cid_col = "CID",
    inchikey_col = "InChIKey",
    smiles_col = "SMILES",
    compound_id_col = "compound_id",
    group_cols = level,
    ...
  )

  result$PlantCompoundMembership = result$GroupCompoundMembership
  result$PlantPairTanimotoSummary = .tanimoto_label_plant_summary(
    result$GroupPairTanimotoSummary,
    level = level
  )
  result$PlantCompoundTanimoto = .tanimoto_label_plant_pairs(
    result$GroupCompoundTanimoto,
    level = level
  )
  result$Provenance$Function = "plantChemicalTanimotoSimilarity"
  class(result) = c("uaf_plant_tanimoto_similarity",
                    "uaf_tanimoto_similarity", "list")
  result
}

.tanimoto_normalize_input = function(compounds, compound_col, cid_col,
                                     inchikey_col, smiles_col,
                                     compound_id_col, group_cols) {
  if (is.character(compounds)) {
    raw = data.frame(compound_name = compounds, stringsAsFactors = FALSE)
  } else if (is.data.frame(compounds)) {
    raw = compounds
  } else {
    stop("`compounds` must be a character vector or data frame.",
         call. = FALSE)
  }
  if (nrow(raw) < 1) stop("`compounds` contains no rows.", call. = FALSE)

  compound_col = .tanimoto_detect_col(raw, compound_col,
                                      c("compound_name", "Query", "Chemical",
                                        "name", "compound"))
  cid_col = .tanimoto_detect_col(raw, cid_col,
                                 c("CID", "pubchem_cid", "cid",
                                   "PubChemCID"))
  inchikey_col = .tanimoto_detect_col(raw, inchikey_col,
                                      c("InChIKey", "inchikey",
                                        "pubchem_inchikey"))
  smiles_col = .tanimoto_detect_col(raw, smiles_col,
                                    c("SMILES", "smiles", "CanonicalSMILES",
                                      "IsomericSMILES"))
  compound_id_col = .tanimoto_detect_col(raw, compound_id_col,
                                         c("compound_id", "CompoundID",
                                           "compound_key"))
  if (is.null(compound_col) && is.null(cid_col) && is.null(inchikey_col)) {
    stop("Could not detect a compound name, CID, or InChIKey column.",
         call. = FALSE)
  }

  group_cols = .tanimoto_validate_group_cols(raw, group_cols)
  compound_name = if (!is.null(compound_col)) {
    .uaf_squish_text(raw[[compound_col]])
  } else {
    rep(NA_character_, nrow(raw))
  }
  cid = if (!is.null(cid_col)) .uaf_squish_text(raw[[cid_col]]) else {
    rep(NA_character_, nrow(raw))
  }
  inchikey = if (!is.null(inchikey_col)) {
    .uaf_squish_text(raw[[inchikey_col]])
  } else {
    rep(NA_character_, nrow(raw))
  }
  smiles = if (!is.null(smiles_col)) {
    .uaf_squish_text(raw[[smiles_col]])
  } else {
    rep(NA_character_, nrow(raw))
  }
  compound_id = if (!is.null(compound_id_col)) {
    .uaf_squish_text(raw[[compound_id_col]])
  } else {
    .tanimoto_make_compound_id(inchikey, cid, compound_name)
  }
  missing_id = is.na(compound_id) | compound_id == ""
  compound_id[missing_id] = .tanimoto_make_compound_id(
    inchikey[missing_id], cid[missing_id], compound_name[missing_id]
  )

  normalized = data.frame(
    compound_id = compound_id,
    compound_name = compound_name,
    CID = cid,
    InChIKey = inchikey,
    SMILES = smiles,
    stringsAsFactors = FALSE
  )
  normalized = normalized[!is.na(normalized$compound_id) &
                            normalized$compound_id != "", , drop = FALSE]
  normalized = normalized[!duplicated(normalized$compound_id), , drop = FALSE]
  row.names(normalized) = NULL

  membership = .uaf_empty_table(c(group_cols, "group_id", "compound_id",
                                  "compound_name"))
  if (length(group_cols) > 0) {
    membership = raw[, group_cols, drop = FALSE]
    for (col in group_cols) membership[[col]] = .uaf_squish_text(membership[[col]])
    membership$group_id = .tanimoto_group_id(membership, group_cols)
    membership$compound_id = compound_id
    membership$compound_name = compound_name
    extra_cols = setdiff(names(raw), names(membership))
    for (col in extra_cols) membership[[col]] = raw[[col]]
    membership = membership[!is.na(membership$group_id) &
                              !is.na(membership$compound_id), , drop = FALSE]
    membership = membership[!duplicated(membership[, c("group_id",
                                                       "compound_id")]),
                            , drop = FALSE]
    row.names(membership) = NULL
  }
  list(compounds = normalized, membership = membership)
}

.tanimoto_detect_col = function(data, explicit, candidates) {
  if (!is.null(explicit)) {
    if (!explicit %in% names(data)) {
      stop("Column `", explicit, "` was not found.", call. = FALSE)
    }
    return(explicit)
  }
  hit = candidates[candidates %in% names(data)]
  if (length(hit) > 0) return(hit[[1]])
  lower = tolower(names(data))
  idx = match(tolower(candidates), lower, nomatch = 0)
  idx = idx[idx > 0]
  if (length(idx) > 0) names(data)[[idx[[1]]]] else NULL
}

.tanimoto_validate_group_cols = function(data, group_cols) {
  if (is.null(group_cols) || length(group_cols) < 1) return(character())
  missing = setdiff(group_cols, names(data))
  if (length(missing) > 0) {
    stop("Grouping columns were not found: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  group_cols
}

.tanimoto_make_compound_id = function(inchikey, cid, name) {
  id = .tanimoto_sanitize_id(inchikey)
  missing = is.na(id) | id == ""
  cid_id = .tanimoto_sanitize_id(cid)
  id[missing & !is.na(cid_id)] = paste0("cid_", cid_id[missing & !is.na(cid_id)])
  missing = is.na(id) | id == ""
  name_id = .tanimoto_sanitize_id(name)
  id[missing & !is.na(name_id)] = paste0("name_", name_id[missing & !is.na(name_id)])
  id
}

.tanimoto_sanitize_id = function(x) {
  x = .uaf_squish_text(x)
  x = gsub("[^A-Za-z0-9]+", "_", x, perl = TRUE)
  x = gsub("^_+|_+$", "", x, perl = TRUE)
  x[x == ""] = NA_character_
  x
}

.tanimoto_group_id = function(data, group_cols) {
  if (length(group_cols) == 1) return(.uaf_squish_text(data[[group_cols]]))
  apply(data[, group_cols, drop = FALSE], 1, function(row) {
    paste(paste(group_cols, .uaf_squish_text(row), sep = "="),
          collapse = "|")
  })
}

.tanimoto_resolve_cids = function(compounds, fetch, name_fallback) {
  rows = lapply(seq_len(nrow(compounds)), function(i) {
    row = compounds[i, , drop = FALSE]
    cid = suppressWarnings(as.integer(row$CID))
    status = if (!is.na(cid)) "input_cid" else "not_queried"
    source = if (!is.na(cid)) "input_cid" else NA_character_
    url = NA_character_
    error = NA_character_
    if (is.na(cid) && !is.na(row$InChIKey) && row$InChIKey != "") {
      url = paste0(.pubchem_base_url(), "/pug/compound/inchikey/",
                   .pubchem_encode_path(row$InChIKey), "/cids/JSON")
      json = fetch(url)
      if (is.null(json)) {
        status = "inchikey_error"
        error = "PubChem InChIKey request failed."
      } else {
        cid = .tanimoto_json_cid(json)
        status = ifelse(is.na(cid), "inchikey_no_hit", "inchikey_resolved")
        source = ifelse(is.na(cid), NA_character_, "pubchem_inchikey")
      }
    }
    if (is.na(cid) && isTRUE(name_fallback) &&
        !is.na(row$compound_name) && row$compound_name != "") {
      url = paste0(.pubchem_base_url(), "/pug/compound/name/",
                   .pubchem_encode_path(row$compound_name), "/cids/JSON")
      json = fetch(url)
      if (is.null(json)) {
        status = "name_error"
        error = "PubChem name request failed."
      } else {
        cid = .tanimoto_json_cid(json)
        status = ifelse(is.na(cid), "name_no_hit", "name_resolved")
        source = ifelse(is.na(cid), NA_character_, "pubchem_name")
      }
    }
    data.frame(
      compound_id = row$compound_id,
      compound_name = row$compound_name,
      input_CID = row$CID,
      input_InChIKey = row$InChIKey,
      input_SMILES = row$SMILES,
      pubchem_cid = cid,
      cid_resolution_status = status,
      cid_resolution_source = source,
      cid_resolution_url = url,
      cid_resolution_error = error,
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.tanimoto_json_cid = function(json) {
  cid = tryCatch(json$IdentifierList$CID[[1]], error = function(error) NA)
  cid = suppressWarnings(as.integer(cid))
  if (length(cid) != 1 || is.na(cid)) NA_integer_ else cid
}

.tanimoto_fetch_fingerprints = function(resolution, fetch, chunk_size = 100) {
  cids = unique(resolution$pubchem_cid[!is.na(resolution$pubchem_cid)])
  cols = c("compound_id", "compound_name", "pubchem_cid", "Title",
           "MolecularFormula", "InChIKey", "CanonicalSMILES",
           "IsomericSMILES", "Fingerprint2D", "fingerprint_status",
           "fingerprint_url")
  if (length(cids) < 1) return(.uaf_empty_table(cols))
  chunks = split(cids, ceiling(seq_along(cids) / chunk_size))
  props = paste(c("Fingerprint2D", "CanonicalSMILES", "IsomericSMILES",
                  "InChIKey", "MolecularFormula", "Title"), collapse = ",")
  fetched = list()
  for (chunk in chunks) {
    url = paste0(.pubchem_base_url(), "/pug/compound/cid/",
                 paste(chunk, collapse = ","), "/property/", props, "/JSON")
    json = fetch(url)
    rows = tryCatch(json$PropertyTable$Properties,
                    error = function(error) list())
    if (length(rows) < 1) next
    fetched[[length(fetched) + 1L]] = do.call(rbind, lapply(rows, function(x) {
      data.frame(
        pubchem_cid = suppressWarnings(as.integer(.uaf_first_non_empty_text(x$CID))),
        Title = .uaf_first_non_empty_text(x$Title),
        MolecularFormula = .uaf_first_non_empty_text(x$MolecularFormula),
        InChIKey = .uaf_first_non_empty_text(x$InChIKey),
        CanonicalSMILES = .uaf_first_non_empty_text(x$CanonicalSMILES),
        IsomericSMILES = .uaf_first_non_empty_text(x$IsomericSMILES),
        Fingerprint2D = .uaf_first_non_empty_text(x$Fingerprint2D),
        fingerprint_status = ifelse(
          is.na(.uaf_first_non_empty_text(x$Fingerprint2D)),
          "fingerprint_missing", "fingerprint_resolved"
        ),
        fingerprint_url = url,
        stringsAsFactors = FALSE
      )
    }))
  }
  if (length(fetched) < 1) return(.uaf_empty_table(cols))
  props = do.call(rbind, fetched)
  out = merge(resolution, props, by = "pubchem_cid", all.x = TRUE)
  out$fingerprint_status[is.na(out$fingerprint_status)] =
    "fingerprint_not_available"
  out = out[, intersect(cols, names(out)), drop = FALSE]
  row.names(out) = NULL
  out
}

.tanimoto_decode_fingerprints = function(fingerprints) {
  if (!requireNamespace("ChemmineR", quietly = TRUE)) {
    stop("ChemmineR is required to decode PubChem Fingerprint2D values.",
         call. = FALSE)
  }
  fp = fingerprints$Fingerprint2D
  names(fp) = fingerprints$compound_id
  keep = !is.na(fp) & fp != "" & !is.na(names(fp)) & names(fp) != ""
  fp = fp[keep]
  if (length(fp) < 2) {
    stop("Fewer than two fingerprints were available.", call. = FALSE)
  }
  bits = ChemmineR::fp2bit(fp, type = 2)
  storage.mode(bits) = "numeric"
  bits
}

.tanimoto_compound_pair_cols = function() {
  c("compound_pair_id", "compound_id_a", "compound_name_a",
    "pubchem_cid_a", "inchikey_a", "compound_id_b", "compound_name_b",
    "pubchem_cid_b", "inchikey_b", "tanimoto", "intersection_bits",
    "union_bits", "fingerprint_source")
}

.tanimoto_compound_pair_table = function(bits, compounds, block_size) {
  n = nrow(bits)
  counts = rowSums(bits)
  rows = list()
  starts = seq(1L, n - 1L, by = block_size)
  for (start in starts) {
    end = min(start + block_size - 1L, n - 1L)
    block_idx = start:end
    inter = bits[block_idx, , drop = FALSE] %*% t(bits)
    for (k in seq_along(block_idx)) {
      i = block_idx[[k]]
      j = seq.int(i + 1L, n)
      rows[[length(rows) + 1L]] =
        .tanimoto_compound_pair_rows(compounds, counts, inter[k, j],
                                     i = i, j = j)
    }
  }
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.tanimoto_write_compound_pairs = function(bits, compounds, path, block_size) {
  n = nrow(bits)
  counts = rowSums(bits)
  con = gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste(.tanimoto_compound_pair_cols(), collapse = ","), con)
  starts = seq(1L, n - 1L, by = block_size)
  for (start in starts) {
    end = min(start + block_size - 1L, n - 1L)
    block_idx = start:end
    inter = bits[block_idx, , drop = FALSE] %*% t(bits)
    rows = lapply(seq_along(block_idx), function(k) {
      i = block_idx[[k]]
      j = seq.int(i + 1L, n)
      .tanimoto_compound_pair_rows(compounds, counts, inter[k, j],
                                   i = i, j = j)
    })
    utils::write.table(do.call(rbind, rows), con, sep = ",",
                       row.names = FALSE, col.names = FALSE, quote = TRUE,
                       na = "")
  }
  invisible(path)
}

.tanimoto_compound_pair_rows = function(compounds, counts, intersections,
                                        i, j) {
  intersection = as.numeric(intersections)
  union = counts[[i]] + counts[j] - intersection
  tanimoto = ifelse(union > 0, intersection / union, NA_real_)
  data.frame(
    compound_pair_id = paste(compounds$compound_id[[i]],
                             compounds$compound_id[j], sep = "__"),
    compound_id_a = compounds$compound_id[[i]],
    compound_name_a = compounds$compound_name[[i]],
    pubchem_cid_a = compounds$pubchem_cid[[i]],
    inchikey_a = compounds$InChIKey[[i]],
    compound_id_b = compounds$compound_id[j],
    compound_name_b = compounds$compound_name[j],
    pubchem_cid_b = compounds$pubchem_cid[j],
    inchikey_b = compounds$InChIKey[j],
    tanimoto = round(tanimoto, 6),
    intersection_bits = intersection,
    union_bits = union,
    fingerprint_source = "PubChem_Fingerprint2D",
    stringsAsFactors = FALSE
  )
}

.tanimoto_group_summary_cols = function(group_cols) {
  c("group_pair_id", "group_a", "group_b", "group_level",
    "compounds_a", "compounds_b", "compound_pair_count",
    "shared_compound_count", "mean_tanimoto", "median_tanimoto",
    "p95_tanimoto", "max_tanimoto", "top_compound_pair_ids",
    "top_compound_pairs", "fingerprint_source")
}

.tanimoto_group_pair_summary = function(bits, membership, compounds,
                                        group_cols, thresholds,
                                        top_n_pairs) {
  if (nrow(membership) < 2) {
    return(.uaf_empty_table(c(.tanimoto_group_summary_cols(group_cols),
                              .tanimoto_threshold_cols(thresholds))))
  }
  membership$fp_index = match(membership$compound_id, rownames(bits))
  membership = membership[!is.na(membership$fp_index) &
                            !is.na(membership$group_id), , drop = FALSE]
  groups = sort(unique(membership$group_id))
  if (length(groups) < 2) {
    return(.uaf_empty_table(c(.tanimoto_group_summary_cols(group_cols),
                              .tanimoto_threshold_cols(thresholds))))
  }
  counts = rowSums(bits)
  compound_names = stats::setNames(compounds$compound_name,
                                   compounds$compound_id)
  rows = list()
  for (i in seq_len(length(groups) - 1L)) {
    a = membership[membership$group_id == groups[[i]], , drop = FALSE]
    for (j in seq.int(i + 1L, length(groups))) {
      b = membership[membership$group_id == groups[[j]], , drop = FALSE]
      inter = bits[a$fp_index, , drop = FALSE] %*%
        t(bits[b$fp_index, , drop = FALSE])
      denom = outer(counts[a$fp_index], counts[b$fp_index], "+") - inter
      sim = ifelse(denom > 0, inter / denom, NA_real_)
      sim_vec = as.numeric(sim)
      sim_vec = sim_vec[!is.na(sim_vec)]
      if (length(sim_vec) < 1) next
      top = .tanimoto_top_pair_labels(sim, a, b, compound_names, top_n_pairs)
      thresholds_df = as.data.frame(stats::setNames(
        lapply(thresholds, function(th) sum(sim_vec >= th, na.rm = TRUE)),
        .tanimoto_threshold_cols(thresholds)
      ), stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] = cbind(data.frame(
        group_pair_id = paste(groups[[i]], groups[[j]], sep = "__"),
        group_a = groups[[i]],
        group_b = groups[[j]],
        group_level = paste(group_cols, collapse = "+"),
        compounds_a = nrow(a),
        compounds_b = nrow(b),
        compound_pair_count = length(sim_vec),
        shared_compound_count = length(intersect(a$compound_id,
                                                 b$compound_id)),
        mean_tanimoto = round(mean(sim_vec), 6),
        median_tanimoto = round(stats::median(sim_vec), 6),
        p95_tanimoto = round(unname(stats::quantile(sim_vec, 0.95,
                                                    na.rm = TRUE)), 6),
        max_tanimoto = round(max(sim_vec), 6),
        top_compound_pair_ids = paste(top$ids, collapse = "; "),
        top_compound_pairs = paste(top$labels, collapse = "; "),
        fingerprint_source = "PubChem_Fingerprint2D",
        stringsAsFactors = FALSE
      ), thresholds_df)
    }
  }
  if (length(rows) < 1) {
    return(.uaf_empty_table(c(.tanimoto_group_summary_cols(group_cols),
                              .tanimoto_threshold_cols(thresholds))))
  }
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out[order(-out$max_tanimoto, -out$mean_tanimoto, out$group_a,
            out$group_b), , drop = FALSE]
}

.tanimoto_threshold_cols = function(thresholds) {
  paste0("compound_pair_count_ge_", gsub("\\.", "_", thresholds))
}

.tanimoto_top_pair_labels = function(sim, a, b, compound_names, top_n_pairs) {
  ord = order(as.numeric(sim), decreasing = TRUE, na.last = NA)
  if (length(ord) < 1) return(list(ids = character(), labels = character()))
  arr = arrayInd(utils::head(ord, top_n_pairs), .dim = dim(sim))
  ids = labels = character(nrow(arr))
  for (z in seq_len(nrow(arr))) {
    ca = a$compound_id[[arr[z, 1]]]
    cb = b$compound_id[[arr[z, 2]]]
    score = round(sim[arr[z, 1], arr[z, 2]], 3)
    ids[[z]] = paste0(ca, " | ", cb, "=", score)
    labels[[z]] = paste0(compound_names[[ca]], " | ",
                         compound_names[[cb]], "=", score)
  }
  list(ids = ids, labels = labels)
}

.tanimoto_group_pair_cols = function(group_cols) {
  c("group_compound_pair_id", "group_a", "compound_id_a",
    "compound_name_a", "group_b", "compound_id_b", "compound_name_b",
    "tanimoto", "intersection_bits", "union_bits", "fingerprint_source")
}

.tanimoto_group_pair_count = function(membership) {
  if (nrow(membership) < 2) return(0)
  total = choose(nrow(membership), 2)
  within = sum(vapply(split(membership$compound_id, membership$group_id),
                      function(x) choose(length(x), 2), numeric(1)))
  total - within
}

.tanimoto_group_compound_pair_table = function(bits, membership, compounds,
                                               block_size) {
  tf = tempfile(fileext = ".csv.gz")
  .tanimoto_write_group_compound_pairs(bits, membership, compounds, tf,
                                       block_size = block_size)
  utils::read.csv(gzfile(tf), stringsAsFactors = FALSE, check.names = FALSE)
}

.tanimoto_write_group_compound_pairs = function(bits, membership, compounds,
                                                path, block_size) {
  membership$fp_index = match(membership$compound_id, rownames(bits))
  membership = membership[!is.na(membership$fp_index) &
                            !is.na(membership$group_id), , drop = FALSE]
  membership = membership[order(membership$group_id, membership$compound_id),
                          , drop = FALSE]
  counts = rowSums(bits)
  compound_names = stats::setNames(compounds$compound_name,
                                   compounds$compound_id)
  con = gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste(.tanimoto_group_pair_cols(character()), collapse = ","),
             con)
  n = nrow(membership)
  starts = seq(1L, n - 1L, by = block_size)
  for (start in starts) {
    end = min(start + block_size - 1L, n - 1L)
    block_idx = start:end
    inter = bits[membership$fp_index[block_idx], , drop = FALSE] %*%
      t(bits[membership$fp_index, , drop = FALSE])
    rows = list()
    for (k in seq_along(block_idx)) {
      i = block_idx[[k]]
      j = seq.int(i + 1L, n)
      j = j[membership$group_id[j] != membership$group_id[[i]]]
      if (length(j) < 1) next
      intersection = as.numeric(inter[k, j])
      union = counts[membership$fp_index[[i]]] +
        counts[membership$fp_index[j]] - intersection
      tanimoto = ifelse(union > 0, intersection / union, NA_real_)
      ca = membership$compound_id[[i]]
      cb = membership$compound_id[j]
      rows[[length(rows) + 1L]] = data.frame(
        group_compound_pair_id = paste(
          paste(membership$group_id[[i]], ca, sep = "::"),
          paste(membership$group_id[j], cb, sep = "::"),
          sep = "__"
        ),
        group_a = membership$group_id[[i]],
        compound_id_a = ca,
        compound_name_a = compound_names[[ca]],
        group_b = membership$group_id[j],
        compound_id_b = cb,
        compound_name_b = unname(compound_names[cb]),
        tanimoto = round(tanimoto, 6),
        intersection_bits = intersection,
        union_bits = union,
        fingerprint_source = "PubChem_Fingerprint2D",
        stringsAsFactors = FALSE
      )
    }
    if (length(rows) > 0) {
      utils::write.table(do.call(rbind, rows), con, sep = ",",
                         row.names = FALSE, col.names = FALSE, quote = TRUE,
                         na = "")
    }
  }
  invisible(path)
}

.tanimoto_plant_membership_table = function(plant_chemistry, level,
                                            include_review_required) {
  if (is.data.frame(plant_chemistry)) return(plant_chemistry)
  if (!is.list(plant_chemistry) ||
      !is.data.frame(plant_chemistry$PlantCompoundOccurrences) ||
      !is.data.frame(plant_chemistry$CompoundResolution)) {
    stop("`plant_chemistry` must be a plant phytochemistry result or a ",
         "plant-compound data frame.", call. = FALSE)
  }
  occ = plant_chemistry$PlantCompoundOccurrences
  cr = plant_chemistry$CompoundResolution
  for (col in c("compound_name_clean", "compound_name", "CID", "InChIKey",
                "SMILES", "MolecularFormula", "resolved",
                "resolution_source")) {
    if (!col %in% names(cr)) cr[[col]] = NA
  }
  cr$compound_id = .tanimoto_make_compound_id(cr$InChIKey, cr$CID,
                                              cr$compound_name_clean)
  cr = cr[cr$resolved %in% TRUE & !is.na(cr$compound_id), , drop = FALSE]
  if (isTRUE(include_review_required) &&
      is.data.frame(plant_chemistry$CompoundIdentityReview) &&
      nrow(plant_chemistry$CompoundIdentityReview) > 0) {
    review = plant_chemistry$CompoundIdentityReview
    if (!"compound_name_clean" %in% names(review)) {
      review$compound_name_clean = NA_character_
    }
    if (!"review_required" %in% names(review)) {
      review$review_required = TRUE
    }
    flags = stats::aggregate(
      list(identity_review_required = review$review_required %in% TRUE),
      by = list(compound_name_clean = .uaf_squish_text(review$compound_name_clean)),
      FUN = function(x) any(x %in% TRUE, na.rm = TRUE)
    )
    cr = merge(cr, flags, by = "compound_name_clean", all.x = TRUE)
    cr$identity_review_required[is.na(cr$identity_review_required)] = FALSE
  } else {
    cr$identity_review_required = FALSE
  }
  keep = c("compound_name_clean", "compound_id", "CID", "InChIKey", "SMILES",
           "MolecularFormula", "resolution_source",
           "identity_review_required")
  resolved = cr[, keep, drop = FALSE]
  names(resolved)[names(resolved) == "compound_id"] = "resolved_compound_id"
  joined = merge(occ, resolved, by = "compound_name_clean")
  if ("compound_id" %in% names(joined)) {
    joined$source_compound_id = joined$compound_id
  }
  joined$compound_id = joined$resolved_compound_id
  if (!"compound_name" %in% names(joined) && "compound_name.x" %in% names(joined)) {
    joined$compound_name = joined$compound_name.x
  }
  if (!level %in% names(joined)) {
    stop("Plant occurrence table does not contain `", level, "`.",
         call. = FALSE)
  }
  joined = joined[!is.na(joined[[level]]) & !is.na(joined$compound_id), ,
                  drop = FALSE]
  joined = joined[!duplicated(joined[, c(level, "compound_id")]), ,
                  drop = FALSE]
  row.names(joined) = NULL
  joined
}

.tanimoto_label_plant_summary = function(summary, level) {
  if (!is.data.frame(summary) || nrow(summary) < 1) return(summary)
  names(summary)[names(summary) == "group_a"] = paste0(level, "_a")
  names(summary)[names(summary) == "group_b"] = paste0(level, "_b")
  names(summary)[names(summary) == "group_pair_id"] =
    paste0(level, "_pair_id")
  summary
}

.tanimoto_label_plant_pairs = function(pairs, level) {
  if (!is.data.frame(pairs) || nrow(pairs) < 1) return(pairs)
  names(pairs)[names(pairs) == "group_a"] = paste0(level, "_a")
  names(pairs)[names(pairs) == "group_b"] = paste0(level, "_b")
  names(pairs)[names(pairs) == "group_compound_pair_id"] =
    paste0(level, "_compound_pair_id")
  pairs
}

.tanimoto_empty_manifest = function() {
  .uaf_empty_table(c("Table", "Path", "Rows", "Description"))
}

.tanimoto_manifest_row = function(table, path, rows, description) {
  data.frame(Table = table,
             Path = normalizePath(path, winslash = "/", mustWork = FALSE),
             Rows = rows,
             Description = description,
             stringsAsFactors = FALSE)
}
