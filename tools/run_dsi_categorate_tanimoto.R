#!/usr/bin/env Rscript

file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path = if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[[1]])
} else {
  file.path(getwd(), "tools", "run_dsi_categorate_tanimoto.R")
}
repo_root = normalizePath(file.path(dirname(script_path), ".."),
                          winslash = "/", mustWork = FALSE)

load_uafr = function(repo_root) {
  if (file.exists(file.path(repo_root, "DESCRIPTION")) &&
      requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(repo_root, quiet = TRUE)
    return(invisible(TRUE))
  }
  if (requireNamespace("uafR", quietly = TRUE)) {
    suppressPackageStartupMessages(library(uafR))
    return(invisible(TRUE))
  }
  stop("Could not load uafR. Install the current package or run this script ",
       "from the uafR repository with devtools installed.", call. = FALSE)
}

parse_args = function(args) {
  out = list()
  i = 1L
  while (i <= length(args)) {
    arg = args[[i]]
    if (!grepl("^--", arg)) {
      stop("Unexpected argument: ", arg, call. = FALSE)
    }
    key = sub("^--", "", arg)
    value = TRUE
    if (grepl("=", key, fixed = TRUE)) {
      parts = strsplit(key, "=", fixed = TRUE)[[1]]
      key = parts[[1]]
      value = paste(parts[-1], collapse = "=")
    } else if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
      value = args[[i + 1L]]
      i = i + 1L
    }
    key = gsub("-", "_", key)
    out[[key]] = value
    i = i + 1L
  }
  out
}

usage = function() {
  cat(
    "Usage:\n",
    "Rscript tools/run_dsi_categorate_tanimoto.R \\\n",
    "  --plant-result-rds lotus_cache/exports/dsi_lotus_identity_source_resolved_20260611/dsi_lotus_phyto_identity_source_resolved_result.rds \\\n",
    "  --out-dir lotus_cache/exports/dsi_categorate_tanimoto_20260611 \\\n",
    "  --cache-dir lotus_cache/exports/dsi_categorate_tanimoto_20260611/cache \\\n",
    "  --write-all-plant-compound-pairs\n\n",
    "Optional categorate pass:\n",
    "  --run-categorate --skip-tanimoto --categorate-query-mode inchikey \\\n",
    "    --categorate-detail research --categorate-batch-size 50\n\n",
    "Notes:\n",
    "- The Tanimoto path uses PubChem Fingerprint2D when a CID can be resolved.\n",
    "- The all plant-compound pair file can be large; use it for final DSI runs.\n",
    "- The species-pair summary is compact and is the primary file to join ",
    "with phylogenetic similarity.\n",
    sep = ""
  )
}

required_arg = function(args, name) {
  value = args[[name]]
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    stop("Missing required argument `--", gsub("_", "-", name), "`.",
         call. = FALSE)
  }
  value
}

optional_arg = function(args, name, default = NULL) {
  value = args[[name]]
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    return(default)
  }
  value
}

flag_arg = function(args, name, default = FALSE) {
  value = args[[name]]
  if (is.null(value)) return(default)
  if (identical(value, TRUE)) return(TRUE)
  tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

numeric_arg = function(args, name, default) {
  value = optional_arg(args, name, NULL)
  if (is.null(value)) return(default)
  if (tolower(value) %in% c("inf", "infinite")) return(Inf)
  out = suppressWarnings(as.numeric(value))
  if (length(out) != 1 || is.na(out) ||
      (!is.finite(out) && !is.infinite(out))) {
    stop("Argument `--", gsub("_", "-", name), "` must be numeric.",
         call. = FALSE)
  }
  out
}

integer_arg = function(args, name, default) {
  as.integer(numeric_arg(args, name, default))
}

split_arg = function(value, default = character()) {
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    return(default)
  }
  pieces = unlist(strsplit(value, "\\s*[,;]\\s*", perl = TRUE),
                  use.names = FALSE)
  pieces = trimws(pieces)
  pieces[nzchar(pieces)]
}

clean_text = function(x) {
  x = as.character(x)
  x = gsub("[[:space:]]+", " ", x, perl = TRUE)
  x = trimws(x)
  x[x == ""] = NA_character_
  x
}

first_non_empty = function(...) {
  values = clean_text(unlist(list(...), use.names = FALSE))
  values = values[!is.na(values)]
  if (length(values) < 1) return(NA_character_)
  values[[1]]
}

safe_col = function(x, col, default = NA_character_) {
  if (is.data.frame(x) && col %in% names(x)) return(x[[col]])
  rep(default, if (is.data.frame(x)) nrow(x) else 0)
}

hash_text = function(x) {
  x = paste0(x, collapse = "\n")
  ints = utf8ToInt(x)
  val = sum((ints * seq_along(ints)) %% .Machine$integer.max)
  paste0(nchar(x), "_", sprintf("%08x", as.integer(val %% .Machine$integer.max)))
}

sanitize_id = function(x) {
  x = clean_text(x)
  x = gsub("[^A-Za-z0-9]+", "_", x, perl = TRUE)
  x = gsub("^_+|_+$", "", x, perl = TRUE)
  x[nchar(x) < 1] = NA_character_
  x
}

make_compound_id = function(inchikey, cid, name_clean) {
  out = sanitize_id(inchikey)
  missing = is.na(out)
  cid_clean = sanitize_id(cid)
  out[missing & !is.na(cid_clean)] = paste0("cid_", cid_clean[missing & !is.na(cid_clean)])
  missing = is.na(out)
  name_clean = sanitize_id(name_clean)
  out[missing & !is.na(name_clean)] = paste0("name_", name_clean[missing & !is.na(name_clean)])
  out
}

make_pair_id = function(a, b) {
  paste(a, b, sep = "__")
}

read_plant_result = function(path) {
  if (!file.exists(path)) {
    stop("Plant result RDS does not exist: ", path, call. = FALSE)
  }
  result = readRDS(path)
  required = c("PlantCompoundOccurrences", "CompoundResolution")
  missing = setdiff(required, names(result))
  if (length(missing) > 0) {
    stop("Plant result is missing required tables: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  result
}

prepare_compounds = function(result, max_compounds = Inf,
                              include_review_required = TRUE) {
  cr = result$CompoundResolution
  for (col in c("compound_name", "compound_name_clean", "resolved", "CID",
                "InChIKey", "SMILES", "MolecularFormula",
                "resolution_source", "notes")) {
    if (!col %in% names(cr)) cr[[col]] = NA
  }
  cr$compound_name = clean_text(cr$compound_name)
  cr$compound_name_clean = clean_text(cr$compound_name_clean)
  cr$CID = clean_text(cr$CID)
  cr$InChIKey = clean_text(cr$InChIKey)
  cr$SMILES = clean_text(cr$SMILES)
  cr$MolecularFormula = clean_text(cr$MolecularFormula)
  cr$resolution_source = clean_text(cr$resolution_source)
  cr$notes = clean_text(cr$notes)
  cr$resolved = cr$resolved %in% TRUE
  cr$compound_id = make_compound_id(cr$InChIKey, cr$CID,
                                    cr$compound_name_clean)

  if (include_review_required &&
      is.data.frame(result$CompoundIdentityReview) &&
      nrow(result$CompoundIdentityReview) > 0) {
    review = result$CompoundIdentityReview
    if (!"compound_name_clean" %in% names(review)) {
      review$compound_name_clean = NA_character_
    }
    if (!"review_required" %in% names(review)) {
      review$review_required = TRUE
    }
    review_flags = aggregate(
      list(identity_review_required = review$review_required %in% TRUE),
      by = list(compound_name_clean = clean_text(review$compound_name_clean)),
      FUN = function(x) any(x %in% TRUE, na.rm = TRUE)
    )
    cr = merge(cr, review_flags, by = "compound_name_clean", all.x = TRUE)
    cr$identity_review_required[is.na(cr$identity_review_required)] = FALSE
  } else {
    cr$identity_review_required = FALSE
  }

  keep = cr$resolved %in% TRUE & !is.na(cr$compound_id) &
    (!is.na(cr$SMILES) | !is.na(cr$CID) | !is.na(cr$InChIKey))
  cr = cr[keep, , drop = FALSE]
  cr = cr[!duplicated(cr$compound_id), , drop = FALSE]
  if (is.finite(max_compounds) && nrow(cr) > max_compounds) {
    cr = cr[seq_len(max_compounds), , drop = FALSE]
  }
  row.names(cr) = NULL
  cr
}

prepare_membership = function(result, compounds) {
  occ = result$PlantCompoundOccurrences
  for (col in c("species", "genus", "family", "matched_rank",
                "compound_name", "compound_name_clean", "source_database",
                "source_record_id", "evidence_tier", "confidence",
                "occurrence_status", "occurrence_basis", "plant_part_group",
                "tissue_group", "method_group", "biological_context_status",
                "evidence_quality_score", "evidence_url", "doi", "pmid")) {
    if (!col %in% names(occ)) occ[[col]] = NA
  }
  occ$species = clean_text(occ$species)
  occ$compound_name_clean = clean_text(occ$compound_name_clean)
  compound_keys = compounds[, c("compound_name_clean", "compound_id", "CID",
                                "InChIKey", "SMILES", "MolecularFormula",
                                "resolution_source",
                                "identity_review_required")]
  names(compound_keys)[names(compound_keys) == "compound_id"] =
    "resolved_compound_id"
  joined = merge(
    occ,
    compound_keys,
    by = "compound_name_clean",
    all.x = FALSE,
    all.y = FALSE
  )
  if ("compound_id" %in% names(joined)) {
    joined$source_compound_id = joined$compound_id
  }
  joined$compound_id = joined$resolved_compound_id
  keep_cols = c("species", "genus", "family", "matched_rank",
                "compound_id", "compound_name", "compound_name_clean",
                "CID", "InChIKey", "SMILES", "MolecularFormula",
                "source_database", "source_record_id", "evidence_tier",
                "confidence", "occurrence_status", "occurrence_basis",
                "plant_part_group", "tissue_group", "method_group",
                "biological_context_status", "evidence_quality_score",
                "evidence_url", "doi", "pmid", "resolution_source",
                "identity_review_required")
  joined = joined[, intersect(keep_cols, names(joined)), drop = FALSE]
  joined = joined[!is.na(joined$species) & !is.na(joined$compound_id), ,
                  drop = FALSE]
  joined = joined[!duplicated(joined[, c("species", "compound_id")]), ,
                  drop = FALSE]
  row.names(joined) = NULL
  joined[order(joined$species, joined$compound_id), , drop = FALSE]
}

fetch_json_cached = function(url, cache_dir, throttle = 0.2,
                             refresh = FALSE) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file = file.path(cache_dir, paste0(hash_text(url), ".json"))
  if (file.exists(cache_file) && !isTRUE(refresh)) {
    txt = paste(readLines(cache_file, warn = FALSE), collapse = "\n")
    return(jsonlite::fromJSON(txt, simplifyVector = FALSE))
  }
  if (is.finite(throttle) && throttle > 0) Sys.sleep(throttle)
  txt = tryCatch(
    paste(readLines(url, warn = FALSE), collapse = "\n"),
    error = function(error) {
      structure(conditionMessage(error), class = c("uaf_fetch_error",
                                                  "character"))
    }
  )
  if (inherits(txt, "uaf_fetch_error")) {
    return(list(.error = as.character(txt), .url = url))
  }
  writeLines(txt, cache_file, useBytes = TRUE)
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

pubchem_cid_from_json = function(json) {
  cid = tryCatch(json$IdentifierList$CID[[1]], error = function(error) NA)
  cid = suppressWarnings(as.integer(cid))
  if (length(cid) != 1 || is.na(cid)) return(NA_integer_)
  cid
}

resolve_pubchem_cids = function(compounds, cache_dir, throttle = 0.2,
                                refresh = FALSE, name_fallback = FALSE) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  rows = vector("list", nrow(compounds))
  base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
  for (i in seq_len(nrow(compounds))) {
    cid = suppressWarnings(as.integer(compounds$CID[[i]]))
    status = if (!is.na(cid)) "input_cid" else "not_queried"
    source = if (!is.na(cid)) "existing_compound_resolution" else NA_character_
    url = NA_character_
    error = NA_character_

    if (is.na(cid) && !is.na(compounds$InChIKey[[i]])) {
      url = paste0(base, "/inchikey/",
                   utils::URLencode(compounds$InChIKey[[i]], reserved = TRUE),
                   "/cids/JSON")
      json = fetch_json_cached(url, cache_dir = cache_dir,
                               throttle = throttle, refresh = refresh)
      if (!is.null(json$.error)) {
        error = json$.error
        status = "inchikey_error"
      } else {
        cid = pubchem_cid_from_json(json)
        status = ifelse(is.na(cid), "inchikey_no_hit",
                        "inchikey_resolved")
        source = ifelse(is.na(cid), NA_character_, "pubchem_inchikey")
      }
    }

    if (is.na(cid) && isTRUE(name_fallback) &&
        !is.na(compounds$compound_name[[i]])) {
      url = paste0(base, "/name/",
                   utils::URLencode(compounds$compound_name[[i]],
                                    reserved = TRUE),
                   "/cids/JSON")
      json = fetch_json_cached(url, cache_dir = cache_dir,
                               throttle = throttle, refresh = refresh)
      if (!is.null(json$.error)) {
        error = json$.error
        status = "name_error"
      } else {
        cid = pubchem_cid_from_json(json)
        status = ifelse(is.na(cid), "name_no_hit", "name_resolved")
        source = ifelse(is.na(cid), NA_character_, "pubchem_name")
      }
    }

    rows[[i]] = data.frame(
      compound_id = compounds$compound_id[[i]],
      compound_name = compounds$compound_name[[i]],
      compound_name_clean = compounds$compound_name_clean[[i]],
      source_CID = compounds$CID[[i]],
      pubchem_cid = cid,
      cid_resolution_status = status,
      cid_resolution_source = source,
      cid_resolution_url = url,
      cid_resolution_error = error,
      stringsAsFactors = FALSE
    )
    if (i %% 100 == 0) message("CID resolution: ", i, "/", nrow(compounds))
  }
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

fetch_pubchem_fingerprints = function(cid_table, cache_dir, throttle = 0.2,
                                      chunk_size = 100, refresh = FALSE) {
  cids = unique(suppressWarnings(as.integer(cid_table$pubchem_cid)))
  cids = cids[!is.na(cids)]
  if (length(cids) < 1) return(data.frame())
  chunks = split(cids, ceiling(seq_along(cids) / chunk_size))
  rows = list()
  base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
  props = paste(c("Fingerprint2D", "CanonicalSMILES", "IsomericSMILES",
                  "InChIKey", "MolecularFormula", "Title"), collapse = ",")
  for (i in seq_along(chunks)) {
    cid_string = paste(chunks[[i]], collapse = ",")
    url = paste0(base, cid_string, "/property/", props, "/JSON")
    json = fetch_json_cached(url, cache_dir = cache_dir,
                             throttle = throttle, refresh = refresh)
    if (!is.null(json$.error)) {
      rows[[length(rows) + 1L]] = data.frame(
        CID = chunks[[i]],
        Fingerprint2D = NA_character_,
        CanonicalSMILES = NA_character_,
        IsomericSMILES = NA_character_,
        InChIKey = NA_character_,
        MolecularFormula = NA_character_,
        Title = NA_character_,
        fingerprint_status = "fetch_error",
        fingerprint_url = url,
        fingerprint_error = json$.error,
        stringsAsFactors = FALSE
      )
      next
    }
    props_rows = tryCatch(json$PropertyTable$Properties,
                          error = function(error) list())
    if (length(props_rows) < 1) next
    parsed = do.call(rbind, lapply(props_rows, function(row) {
      data.frame(
        CID = suppressWarnings(as.integer(first_non_empty(row$CID))),
        Fingerprint2D = first_non_empty(row$Fingerprint2D),
        CanonicalSMILES = first_non_empty(row$CanonicalSMILES),
        IsomericSMILES = first_non_empty(row$IsomericSMILES),
        InChIKey = first_non_empty(row$InChIKey),
        MolecularFormula = first_non_empty(row$MolecularFormula),
        Title = first_non_empty(row$Title),
        fingerprint_status =
          ifelse(is.na(first_non_empty(row$Fingerprint2D)),
                 "fingerprint_missing", "fingerprint_resolved"),
        fingerprint_url = url,
        fingerprint_error = NA_character_,
        stringsAsFactors = FALSE
      )
    }))
    rows[[length(rows) + 1L]] = parsed
    message("Fingerprint chunk: ", i, "/", length(chunks))
  }
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

decode_pubchem_fingerprints = function(fingerprints) {
  if (!requireNamespace("ChemmineR", quietly = TRUE)) {
    stop("ChemmineR is required to decode PubChem Fingerprint2D values.",
         call. = FALSE)
  }
  fp = clean_text(fingerprints$Fingerprint2D)
  names(fp) = fingerprints$compound_id
  keep = !is.na(fp) & !is.na(names(fp)) & names(fp) != ""
  fp = fp[keep]
  if (length(fp) < 2) {
    stop("Fewer than two PubChem fingerprints were available.", call. = FALSE)
  }
  bits = ChemmineR::fp2bit(fp, type = 2)
  storage.mode(bits) = "numeric"
  bits
}

write_compound_pair_tanimoto = function(bits, compounds, out_file,
                                        block_size = 100) {
  n = nrow(bits)
  if (n < 2) stop("Need at least two fingerprinted compounds.", call. = FALSE)
  compounds = compounds[match(rownames(bits), compounds$compound_id), ,
                        drop = FALSE]
  counts = rowSums(bits)
  con = gzfile(out_file, "wt")
  on.exit(close(con), add = TRUE)
  header = c("compound_pair_id", "compound_id_a", "compound_name_a", "cid_a",
             "inchikey_a", "compound_id_b", "compound_name_b", "cid_b",
             "inchikey_b", "tanimoto", "intersection_bits", "union_bits",
             "fingerprint_source")
  writeLines(paste(header, collapse = ","), con)
  starts = seq(1L, n - 1L, by = block_size)
  for (start in starts) {
    end = min(start + block_size - 1L, n - 1L)
    block_idx = start:end
    inter = bits[block_idx, , drop = FALSE] %*% t(bits)
    rows = vector("list", length(block_idx))
    for (k in seq_along(block_idx)) {
      i = block_idx[[k]]
      j = seq.int(i + 1L, n)
      intersection = as.numeric(inter[k, j])
      union = counts[[i]] + counts[j] - intersection
      tanimoto = ifelse(union > 0, intersection / union, NA_real_)
      rows[[k]] = data.frame(
        compound_pair_id = make_pair_id(compounds$compound_id[[i]],
                                        compounds$compound_id[j]),
        compound_id_a = compounds$compound_id[[i]],
        compound_name_a = compounds$compound_name[[i]],
        cid_a = compounds$pubchem_cid[[i]],
        inchikey_a = compounds$InChIKey[[i]],
        compound_id_b = compounds$compound_id[j],
        compound_name_b = compounds$compound_name[j],
        cid_b = compounds$pubchem_cid[j],
        inchikey_b = compounds$InChIKey[j],
        tanimoto = round(tanimoto, 6),
        intersection_bits = intersection,
        union_bits = union,
        fingerprint_source = "PubChem_Fingerprint2D",
        stringsAsFactors = FALSE
      )
    }
    utils::write.table(do.call(rbind, rows), con, sep = ",", row.names = FALSE,
                       col.names = FALSE, quote = TRUE, na = "")
    message("Compound Tanimoto block: ", end, "/", n)
  }
  invisible(out_file)
}

species_pair_summary = function(bits, membership, compounds,
                                top_n_pairs = 5,
                                thresholds = c(0.5, 0.7, 0.85, 0.95)) {
  membership = membership[membership$compound_id %in% rownames(bits), ,
                          drop = FALSE]
  membership$fp_index = match(membership$compound_id, rownames(bits))
  membership = membership[!is.na(membership$species) &
                            !is.na(membership$fp_index), , drop = FALSE]
  species = sort(unique(membership$species))
  counts = rowSums(bits)
  compound_names = stats::setNames(compounds$compound_name,
                                   compounds$compound_id)
  rows = list()
  for (i in seq_len(length(species) - 1L)) {
    sp_a = species[[i]]
    a = membership[membership$species == sp_a, , drop = FALSE]
    for (j in seq.int(i + 1L, length(species))) {
      sp_b = species[[j]]
      b = membership[membership$species == sp_b, , drop = FALSE]
      inter = bits[a$fp_index, , drop = FALSE] %*%
        t(bits[b$fp_index, , drop = FALSE])
      denom = outer(counts[a$fp_index], counts[b$fp_index], "+") - inter
      sim = ifelse(denom > 0, inter / denom, NA_real_)
      sim_vec = as.numeric(sim)
      sim_vec = sim_vec[!is.na(sim_vec)]
      if (length(sim_vec) < 1) next
      ord = order(as.numeric(sim), decreasing = TRUE, na.last = NA)
      top = character()
      top_ids = character()
      if (length(ord) > 0) {
        top_idx = utils::head(ord, top_n_pairs)
        arr = arrayInd(top_idx, .dim = dim(sim))
        top = vapply(seq_len(nrow(arr)), function(z) {
          ca = a$compound_id[[arr[z, 1]]]
          cb = b$compound_id[[arr[z, 2]]]
          paste0(compound_names[[ca]], " | ", compound_names[[cb]],
                 "=", round(sim[arr[z, 1], arr[z, 2]], 3))
        }, character(1))
        top_ids = vapply(seq_len(nrow(arr)), function(z) {
          ca = a$compound_id[[arr[z, 1]]]
          cb = b$compound_id[[arr[z, 2]]]
          paste0(ca, " | ", cb, "=",
                 round(sim[arr[z, 1], arr[z, 2]], 3))
        }, character(1))
      }
      threshold_counts = stats::setNames(
        lapply(thresholds, function(th) sum(sim_vec >= th, na.rm = TRUE)),
        paste0("compound_pair_count_ge_", gsub("\\.", "_", thresholds))
      )
      rows[[length(rows) + 1L]] = data.frame(
        species_a = sp_a,
        species_b = sp_b,
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
        top_compound_pair_ids = paste(top_ids, collapse = "; "),
        top_compound_pairs = paste(top, collapse = "; "),
        stringsAsFactors = FALSE
      )
      rows[[length(rows)]] = cbind(rows[[length(rows)]],
                                   as.data.frame(threshold_counts,
                                                 stringsAsFactors = FALSE))
    }
    message("Species-pair summary: ", i, "/", length(species))
  }
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out[order(-out$max_tanimoto, -out$mean_tanimoto, out$species_a,
            out$species_b), , drop = FALSE]
}

write_plant_compound_pair_tanimoto = function(bits, membership, compounds,
                                              out_file, block_size = 100,
                                              cross_species_only = TRUE) {
  membership = membership[membership$compound_id %in% rownames(bits), ,
                          drop = FALSE]
  membership$fp_index = match(membership$compound_id, rownames(bits))
  membership = membership[!is.na(membership$species) &
                            !is.na(membership$fp_index), , drop = FALSE]
  membership = membership[order(membership$species, membership$compound_id),
                          , drop = FALSE]
  n = nrow(membership)
  if (n < 2) stop("Need at least two species-compound memberships.",
                  call. = FALSE)
  counts = rowSums(bits)
  compound_names = stats::setNames(compounds$compound_name,
                                   compounds$compound_id)
  con = gzfile(out_file, "wt")
  on.exit(close(con), add = TRUE)
  header = c("plant_compound_pair_id", "species_a", "compound_id_a",
             "compound_name_a", "evidence_tier_a", "plant_part_group_a",
             "species_b", "compound_id_b", "compound_name_b",
             "evidence_tier_b", "plant_part_group_b", "tanimoto",
             "intersection_bits", "union_bits", "fingerprint_source")
  writeLines(paste(header, collapse = ","), con)
  starts = seq(1L, n - 1L, by = block_size)
  for (start in starts) {
    end = min(start + block_size - 1L, n - 1L)
    block_idx = start:end
    block_bits = bits[membership$fp_index[block_idx], , drop = FALSE]
    all_bits = bits[membership$fp_index, , drop = FALSE]
    inter = block_bits %*% t(all_bits)
    rows = list()
    for (k in seq_along(block_idx)) {
      i = block_idx[[k]]
      j = seq.int(i + 1L, n)
      if (isTRUE(cross_species_only)) {
        j = j[membership$species[j] != membership$species[[i]]]
      }
      if (length(j) < 1) next
      intersection = as.numeric(inter[k, j])
      union = counts[membership$fp_index[[i]]] +
        counts[membership$fp_index[j]] - intersection
      tanimoto = ifelse(union > 0, intersection / union, NA_real_)
      ca = membership$compound_id[[i]]
      cb = membership$compound_id[j]
      rows[[length(rows) + 1L]] = data.frame(
        plant_compound_pair_id =
          make_pair_id(paste(membership$species[[i]], ca, sep = "::"),
                       paste(membership$species[j], cb, sep = "::")),
        species_a = membership$species[[i]],
        compound_id_a = ca,
        compound_name_a = compound_names[[ca]],
        evidence_tier_a = membership$evidence_tier[[i]],
        plant_part_group_a = membership$plant_part_group[[i]],
        species_b = membership$species[j],
        compound_id_b = cb,
        compound_name_b = unname(compound_names[cb]),
        evidence_tier_b = membership$evidence_tier[j],
        plant_part_group_b = membership$plant_part_group[j],
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
    message("Plant-compound Tanimoto block: ", end, "/", n)
  }
  invisible(out_file)
}

run_categorate_batches = function(compounds, out_dir, cache_dir,
                                  detail = "research", batch_size = 50,
                                  throttle = 0.2, max_batches = Inf) {
  if (!exists("categorate")) {
    stop("categorate() was not loaded.", call. = FALSE)
  }
  data("library_data", package = "uafR", envir = environment())
  query_col = if ("categorate_query" %in% names(compounds)) {
    "categorate_query"
  } else {
    "compound_name"
  }
  queries = compounds[[query_col]]
  queries = queries[!is.na(queries) & queries != ""]
  queries = unique(queries)
  chunks = split(queries, ceiling(seq_along(queries) / batch_size))
  if (is.finite(max_batches)) chunks = chunks[seq_len(min(length(chunks),
                                                        max_batches))]
  batch_dir = file.path(out_dir, "categorate_batches")
  dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
  rows = list()
  for (i in seq_along(chunks)) {
    batch_file = file.path(batch_dir, sprintf("categorate_batch_%04d.rds", i))
    if (file.exists(batch_file)) {
      result = readRDS(batch_file)
      status = if (inherits(result, "error")) "cached_error" else "cached"
    } else {
      result = NULL
      status = "new"
    }
    if (identical(status, "new") || identical(status, "cached_error")) {
      result = tryCatch(
        categorate(
          compounds = chunks[[i]],
          chemical_library = library_data,
          input_format = "wide",
          detail = detail,
          cache = TRUE,
          cache_dir = file.path(cache_dir, "categorate"),
          throttle = throttle,
          trait_matrix_profile = "core",
          trait_matrix_min_confidence = "medium"
        ),
        error = function(error) error
      )
      status = if (inherits(result, "error")) "error" else "completed"
      saveRDS(result, batch_file)
    }
    validation = if (!inherits(result, "error") &&
                     exists("validateCategorateResult")) {
      tryCatch(validateCategorateResult(result)$ValidationSummary,
               error = function(error) data.frame())
    } else {
      data.frame()
    }
    rows[[i]] = data.frame(
      batch_index = i,
      batch_file = batch_file,
      query_count = length(chunks[[i]]),
      status = status,
      error = if (inherits(result, "error")) conditionMessage(result) else NA_character_,
      validation_status = if (is.data.frame(validation) &&
                              "Status" %in% names(validation) &&
                              nrow(validation) > 0) {
        paste(unique(validation$Status), collapse = "; ")
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
    message("categorate batch: ", i, "/", length(chunks), " (", status, ")")
  }
  out = do.call(rbind, rows)
  utils::write.csv(out, file.path(out_dir, "dsi_categorate_batch_summary.csv"),
                   row.names = FALSE, na = "")
  out
}

prepare_categorate_query_map = function(compounds, fp_ready,
                                        mode = c("inchikey",
                                                 "pubchem_title",
                                                 "compound_name")) {
  mode = match.arg(mode)
  fp_cols = c("compound_id", "pubchem_cid", "Title", "InChIKey.y",
              "CanonicalSMILES", "IsomericSMILES")
  fp_cols = intersect(fp_cols, names(fp_ready))
  map = merge(compounds, fp_ready[, fp_cols, drop = FALSE],
              by = "compound_id", all.x = TRUE, suffixes = c("", "_pubchem"))
  pubchem_inchikey_col = if ("InChIKey.y" %in% names(map)) {
    "InChIKey.y"
  } else if ("InChIKey_pubchem" %in% names(map)) {
    "InChIKey_pubchem"
  } else {
    NA_character_
  }
  pubchem_inchikey = if (!is.na(pubchem_inchikey_col)) {
    clean_text(map[[pubchem_inchikey_col]])
  } else {
    rep(NA_character_, nrow(map))
  }
  title = if ("Title" %in% names(map)) {
    clean_text(map$Title)
  } else {
    rep(NA_character_, nrow(map))
  }
  if (identical(mode, "inchikey")) {
    map$categorate_query = clean_text(map$InChIKey)
    map$categorate_query[is.na(map$categorate_query)] =
      pubchem_inchikey[is.na(map$categorate_query)]
    map$categorate_query[is.na(map$categorate_query)] =
      title[is.na(map$categorate_query)]
    map$categorate_query[is.na(map$categorate_query)] =
      map$compound_name[is.na(map$categorate_query)]
  } else if (identical(mode, "pubchem_title")) {
    map$categorate_query = title
    map$categorate_query[is.na(map$categorate_query)] =
      clean_text(map$InChIKey)[is.na(map$categorate_query)]
    map$categorate_query[is.na(map$categorate_query)] =
      map$compound_name[is.na(map$categorate_query)]
  } else {
    map$categorate_query = map$compound_name
  }
  map$categorate_query_mode = mode
  map$categorate_query_source = ifelse(
    !is.na(map$categorate_query) &
      map$categorate_query == clean_text(map$InChIKey),
    "inchikey",
    ifelse(!is.na(map$categorate_query) &
             map$categorate_query == title,
           "pubchem_title", "compound_name")
  )
  keep = c("compound_id", "compound_name", "compound_name_clean",
           "categorate_query", "categorate_query_mode",
           "categorate_query_source", "CID", "pubchem_cid", "InChIKey",
           "Title", "SMILES", "MolecularFormula", "resolution_source",
           "identity_review_required")
  map[, intersect(keep, names(map)), drop = FALSE]
}

write_summary_json = function(path, summary) {
  jsonlite::write_json(summary, path, auto_unbox = TRUE, pretty = TRUE,
                       null = "null")
}

main = function() {
  args = commandArgs(trailingOnly = TRUE)
  if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
    usage()
    quit(status = ifelse(length(args) == 0, 1L, 0L))
  }
  parsed = parse_args(args)
  load_uafr(repo_root)

  plant_result_rds = required_arg(parsed, "plant_result_rds")
  out_dir = required_arg(parsed, "out_dir")
  cache_dir = optional_arg(parsed, "cache_dir", file.path(out_dir, "cache"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  result = read_plant_result(plant_result_rds)
  max_compounds = numeric_arg(parsed, "max_compounds", Inf)
  throttle = numeric_arg(parsed, "throttle", 0.2)
  refresh = flag_arg(parsed, "refresh", FALSE)

  compounds = prepare_compounds(result, max_compounds = max_compounds)
  membership = prepare_membership(result, compounds)
  utils::write.csv(compounds, file.path(out_dir,
                                        "dsi_resolved_compounds_for_categorate.csv"),
                   row.names = FALSE, na = "")
  utils::write.csv(membership, file.path(out_dir,
                                         "dsi_species_compound_membership.csv"),
                   row.names = FALSE, na = "")

  cid_table = resolve_pubchem_cids(
    compounds = compounds,
    cache_dir = file.path(cache_dir, "pubchem_cids"),
    throttle = throttle,
    refresh = refresh,
    name_fallback = flag_arg(parsed, "name_cid_fallback", FALSE)
  )
  utils::write.csv(cid_table, file.path(out_dir, "dsi_pubchem_cid_resolution.csv"),
                   row.names = FALSE, na = "")

  fp_raw = fetch_pubchem_fingerprints(
    cid_table = cid_table,
    cache_dir = file.path(cache_dir, "pubchem_fingerprints"),
    throttle = throttle,
    chunk_size = integer_arg(parsed, "fingerprint_chunk_size", 100),
    refresh = refresh
  )
  if (nrow(fp_raw) < 1) {
    stop("No PubChem fingerprints were fetched. Check CID resolution and ",
         "network/cache status.", call. = FALSE)
  }
  fp_joined = merge(cid_table, fp_raw, by.x = "pubchem_cid", by.y = "CID",
                    all.x = TRUE)
  fp_joined = merge(fp_joined, compounds, by = c("compound_id",
                                                 "compound_name_clean",
                                                 "compound_name"),
                    all.x = TRUE, suffixes = c("", "_source"))
  utils::write.csv(fp_joined, file.path(out_dir, "dsi_pubchem_fingerprints.csv"),
                   row.names = FALSE, na = "")

  fp_ready = fp_joined[!is.na(fp_joined$Fingerprint2D) &
                         fp_joined$Fingerprint2D != "", , drop = FALSE]
  compound_pair_file = file.path(out_dir, "dsi_compound_pair_tanimoto.csv.gz")
  plant_pair_file = NA_character_
  species_summary_file = file.path(out_dir,
                                   "dsi_species_pair_tanimoto_summary.csv")
  if (!flag_arg(parsed, "skip_tanimoto", FALSE)) {
    bits = decode_pubchem_fingerprints(fp_ready)
    write_compound_pair_tanimoto(
      bits = bits,
      compounds = fp_ready,
      out_file = compound_pair_file,
      block_size = integer_arg(parsed, "compound_pair_block_size", 100)
    )

    species_summary = species_pair_summary(
      bits = bits,
      membership = membership,
      compounds = fp_ready,
      top_n_pairs = integer_arg(parsed, "top_n_pairs", 5),
      thresholds = split_arg(optional_arg(parsed, "thresholds", NULL),
                             c("0.5", "0.7", "0.85", "0.95"))
    )
    utils::write.csv(species_summary, species_summary_file,
                     row.names = FALSE, na = "")

    if (flag_arg(parsed, "write_all_plant_compound_pairs", FALSE)) {
      plant_pair_file = file.path(out_dir,
                                  "dsi_plant_compound_pair_tanimoto.csv.gz")
      write_plant_compound_pair_tanimoto(
        bits = bits,
        membership = membership,
        compounds = fp_ready,
        out_file = plant_pair_file,
        block_size = integer_arg(parsed, "plant_pair_block_size", 100),
        cross_species_only = !flag_arg(parsed, "include_within_species", FALSE)
      )
    }
  } else {
    species_summary = if (file.exists(species_summary_file)) {
      utils::read.csv(species_summary_file, stringsAsFactors = FALSE,
                      check.names = FALSE)
    } else {
      data.frame()
    }
    if (file.exists(file.path(out_dir,
                              "dsi_plant_compound_pair_tanimoto.csv.gz"))) {
      plant_pair_file = file.path(out_dir,
                                  "dsi_plant_compound_pair_tanimoto.csv.gz")
    }
  }

  categorate_summary = data.frame()
  if (flag_arg(parsed, "run_categorate", FALSE)) {
    categorate_map = prepare_categorate_query_map(
      compounds = compounds,
      fp_ready = fp_ready,
      mode = optional_arg(parsed, "categorate_query_mode", "inchikey")
    )
    utils::write.csv(categorate_map,
                     file.path(out_dir, "dsi_categorate_query_map.csv"),
                     row.names = FALSE, na = "")
    categorate_summary = run_categorate_batches(
      compounds = categorate_map,
      out_dir = out_dir,
      cache_dir = cache_dir,
      detail = optional_arg(parsed, "categorate_detail", "research"),
      batch_size = integer_arg(parsed, "categorate_batch_size", 50),
      throttle = throttle,
      max_batches = integer_arg(parsed, "categorate_max_batches", Inf)
    )
  }

  summary = list(
    plant_result_rds = normalizePath(plant_result_rds, winslash = "/",
                                     mustWork = FALSE),
    out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    resolved_compounds_input = nrow(compounds),
    species_compound_memberships = nrow(membership),
    pubchem_cids_resolved = sum(!is.na(cid_table$pubchem_cid)),
    pubchem_fingerprints_resolved = nrow(fp_ready),
    species_pair_rows = nrow(species_summary),
    compound_pair_tanimoto_file = compound_pair_file,
    plant_compound_pair_tanimoto_file = plant_pair_file,
    categorate_batches = if (is.data.frame(categorate_summary)) {
      nrow(categorate_summary)
    } else {
      0
    },
    notes = c(
      "Tanimoto values use PubChem Fingerprint2D bit vectors.",
      "Species-pair summaries compare all fingerprinted compounds reported for each species pair.",
      "Plant-compound pair output is cross-species by default unless --include-within-species is used.",
      "categorate() is optional and runs in resumable batches because full DSI enrichment is network-heavy."
    )
  )
  write_summary_json(file.path(out_dir, "dsi_categorate_tanimoto_summary.json"),
                     summary)

  cat("DSI categorate/Tanimoto preparation complete.\n")
  cat("Output directory: ", normalizePath(out_dir, winslash = "/",
                                         mustWork = FALSE), "\n", sep = "")
  cat("Resolved compounds: ", nrow(compounds), "\n", sep = "")
  cat("Species-compound memberships: ", nrow(membership), "\n", sep = "")
  cat("PubChem fingerprints: ", nrow(fp_ready), "\n", sep = "")
  cat("Species-pair summary: ",
      file.path(out_dir, "dsi_species_pair_tanimoto_summary.csv"), "\n",
      sep = "")
}

main()
