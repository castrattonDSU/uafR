#' Export aiNsect molecular-olfaction inputs from wide EO GC-MS tables
#'
#' @description
#' `exportAiNsectMolOlfInputs()` converts paired wide essential-oil GC-MS
#' identity and abundance tables into molecular-olfaction input files for
#' aiNsect. Compound structures are resolved through `pubchemProfile()`; the
#' exporter does not fabricate SMILES, InChIKeys, PubChem CIDs, formulas, or
#' abundance values.
#'
#' @param chem_id_csv Path to the wide compound-identity CSV. Columns are
#' treatments and rows are ranked compound names.
#' @param chem_quant_csv Path to the matching wide abundance CSV. Columns are
#' treatments and rows are abundance values for the same treatment/rank cells.
#' @param out_dir Directory where aiNsect input files will be written.
#' @param cache_dir Directory used by `pubchemProfile()` for PubChem response
#' caching. If `NULL`, `pubchemProfile()` uses its default cache directory.
#' @param profile PubChem enrichment profile passed to `pubchemProfile()`.
#' Defaults to `"ms"`.
#' @param throttle Seconds to wait between uncached PubChem requests.
#' @param min_relative_abundance Minimum relative abundance retained per
#' treatment after unresolved compounds are removed and duplicate compounds are
#' aggregated. Defaults to `0`.
#' @param source_label Source label written to `uafR_compounds.csv`.
#' @param uafR_run_id Optional run identifier. If `NULL`, a timestamped ID is
#' generated.
#' @param treatment_aliases Optional named character vector mapping treatment
#' aliases to canonical treatment names. Built-in aliases handle the known EO
#' table cleanup issues such as `"Clary Sagae"` to `"Clary Sage"`.
#' @param write_pubchem_audit Logical. If `TRUE`, write raw PubChem identity
#' and property audit CSVs when those tables are available.
#' @param profile_fun Advanced/testing hook. Defaults to `pubchemProfile`.
#'
#' @return A list containing output paths, exported data frames, unresolved
#' compounds, and summary metadata. Files written to `out_dir` include
#' `uafR_compounds.csv`, `treatment_compound_abundance.csv`,
#' `uafR_compounds_unresolved.csv`, and
#' `uafR_mololf_export_summary.json`.
#'
#' @examples
#' \dontrun{
#' exportAiNsectMolOlfInputs(
#'   chem_id_csv = "20240612-EO-gcms-data_all.csv",
#'   chem_quant_csv = "20240612-EO-gcms-quant_all.csv",
#'   out_dir = "mololf_export",
#'   cache_dir = "pubchem_cache",
#'   profile = "ms",
#'   throttle = 0.2
#' )
#' }
#'
#' @export
exportAiNsectMolOlfInputs = function(chem_id_csv,
                                     chem_quant_csv,
                                     out_dir,
                                     cache_dir = NULL,
                                     profile = "ms",
                                     throttle = 0.2,
                                     min_relative_abundance = 0,
                                     source_label = "uafR_pubchem_eo_gcms",
                                     uafR_run_id = NULL,
                                     treatment_aliases = NULL,
                                     write_pubchem_audit = TRUE,
                                     profile_fun = pubchemProfile) {
  chem_id_csv = .ainsect_existing_file(chem_id_csv, "chem_id_csv")
  chem_quant_csv = .ainsect_existing_file(chem_quant_csv, "chem_quant_csv")
  out_dir = .ainsect_required_path(out_dir, "out_dir")
  profile = .uaf_first_non_empty_text(profile)
  source_label = .uaf_first_non_empty_text(source_label,
                                           "uafR_pubchem_eo_gcms")
  throttle = .ainsect_number(throttle, "throttle", min_value = 0,
                             allow_na = TRUE)
  min_relative_abundance = .ainsect_number(min_relative_abundance,
                                           "min_relative_abundance",
                                           min_value = 0,
                                           max_value = 1)
  if (is.null(uafR_run_id) || length(.uaf_non_empty(uafR_run_id)) < 1) {
    uafR_run_id = paste0("uafR_aiNsect_",
                         format(Sys.time(), "%Y%m%d_%H%M%S"))
  } else {
    uafR_run_id = .uaf_non_empty(uafR_run_id)[[1]]
  }
  if (!is.function(profile_fun)) {
    stop("`profile_fun` must be a function.", call. = FALSE)
  }

  identity_wide = .ainsect_read_wide_csv(chem_id_csv, "chem_id_csv")
  quant_wide = .ainsect_read_wide_csv(chem_quant_csv, "chem_quant_csv")
  if (nrow(identity_wide) != nrow(quant_wide)) {
    stop("Identity and abundance tables must have the same number of rows. ",
         "Observed ", nrow(identity_wide), " and ", nrow(quant_wide), ".",
         call. = FALSE)
  }

  aliases = .ainsect_treatment_aliases(treatment_aliases)
  alignment = .ainsect_align_treatments(identity_wide, quant_wide, aliases)
  parsed = .ainsect_long_profile(identity_wide = identity_wide,
                                 quant_wide = quant_wide,
                                 alignment = alignment,
                                 chem_id_csv = chem_id_csv,
                                 chem_quant_csv = chem_quant_csv,
                                 uafR_run_id = uafR_run_id)
  if (nrow(parsed$valid_rows) < 1) {
    stop("No valid compound/abundance pairs were found in the input tables.",
         call. = FALSE)
  }

  compound_names = unique(parsed$valid_rows$compound_name)
  compound_names = compound_names[order(tolower(compound_names))]
  pubchem_profile = .ainsect_call_profile_fun(
    profile_fun = profile_fun,
    compounds = compound_names,
    profile = profile,
    cache_dir = cache_dir,
    throttle = throttle
  )
  compound_resolution = .ainsect_compound_resolution(
    compound_names = compound_names,
    pubchem_profile = pubchem_profile,
    source_label = source_label
  )
  compounds_all = compound_resolution$records
  compounds = .ainsect_compounds_for_export(compounds_all)
  unresolved = .ainsect_unresolved_table(compounds_all, parsed$valid_rows)
  if (nrow(compounds) < 1) {
    stop("No compounds resolved to usable PubChem SMILES. Check PubChem ",
         "connectivity, query names, cache contents, or the selected profile.",
         call. = FALSE)
  }

  abundance = .ainsect_abundance_for_export(
    valid_rows = parsed$valid_rows,
    compound_records = compounds_all,
    min_relative_abundance = min_relative_abundance,
    uafR_run_id = uafR_run_id
  )
  if (nrow(abundance) < 1) {
    stop("No treatment/compound abundance rows remained after unresolved ",
         "compounds and abundance filters were applied.", call. = FALSE)
  }
  missing_ids = setdiff(unique(abundance$compound_id), compounds$compound_id)
  if (length(missing_ids) > 0) {
    stop("Internal export error: abundance rows reference compound IDs that ",
         "are absent from `uafR_compounds.csv`: ",
         paste(missing_ids, collapse = ", "), call. = FALSE)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) {
    stop("Could not create output directory: ", out_dir, call. = FALSE)
  }

  paths = list(
    compounds = file.path(out_dir, "uafR_compounds.csv"),
    abundance = file.path(out_dir, "treatment_compound_abundance.csv"),
    unresolved = file.path(out_dir, "uafR_compounds_unresolved.csv"),
    summary = file.path(out_dir, "uafR_mololf_export_summary.json"),
    pubchem_identity_audit = file.path(out_dir,
                                       "uafR_pubchem_identity_audit.csv"),
    pubchem_properties_audit = file.path(out_dir,
                                         "uafR_pubchem_properties_audit.csv")
  )

  utils::write.csv(compounds, paths$compounds, row.names = FALSE, na = "")
  utils::write.csv(abundance, paths$abundance, row.names = FALSE, na = "")
  utils::write.csv(unresolved, paths$unresolved, row.names = FALSE, na = "")
  if (isTRUE(write_pubchem_audit)) {
    .ainsect_write_audit_table(pubchem_profile$identity,
                               paths$pubchem_identity_audit)
    .ainsect_write_audit_table(pubchem_profile$properties,
                               paths$pubchem_properties_audit)
  } else {
    paths$pubchem_identity_audit = NA_character_
    paths$pubchem_properties_audit = NA_character_
  }

  summary = .ainsect_export_summary(
    uafR_run_id = uafR_run_id,
    source_label = source_label,
    profile = profile,
    throttle = throttle,
    min_relative_abundance = min_relative_abundance,
    chem_id_csv = chem_id_csv,
    chem_quant_csv = chem_quant_csv,
    out_dir = out_dir,
    cache_dir = cache_dir,
    alignment = alignment,
    parsed = parsed,
    compounds_all = compounds_all,
    compounds = compounds,
    unresolved = unresolved,
    abundance = abundance,
    paths = paths
  )
  jsonlite::write_json(summary, paths$summary, pretty = TRUE,
                       auto_unbox = TRUE, null = "null")

  out = list(
    paths = paths,
    compounds = compounds,
    treatment_compound_abundance = abundance,
    unresolved = unresolved,
    summary = summary
  )
  class(out) = c("uaf_ainsect_mololf_export", class(out))
  out
}

.ainsect_existing_file = function(path, arg) {
  path = .ainsect_required_path(path, arg)
  if (!file.exists(path)) {
    stop("`", arg, "` does not exist: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.ainsect_required_path = function(path, arg) {
  values = .uaf_non_empty(path)
  if (length(values) < 1) {
    stop("`", arg, "` is required.", call. = FALSE)
  }
  values[[1]]
}

.ainsect_number = function(x, arg, min_value = -Inf, max_value = Inf,
                           allow_na = FALSE) {
  value = suppressWarnings(as.numeric(x[[1]]))
  if (is.na(value)) {
    if (isTRUE(allow_na)) return(value)
    stop("`", arg, "` must be numeric.", call. = FALSE)
  }
  if (!is.infinite(min_value) && value < min_value) {
    stop("`", arg, "` must be >= ", min_value, ".", call. = FALSE)
  }
  if (!is.infinite(max_value) && value > max_value) {
    stop("`", arg, "` must be <= ", max_value, ".", call. = FALSE)
  }
  value
}

.ainsect_read_wide_csv = function(path, arg) {
  out = tryCatch(
    utils::read.csv(path,
                    check.names = FALSE,
                    stringsAsFactors = FALSE,
                    na.strings = c("", "NA", "N/A", "NULL")),
    error = function(error) {
      stop("Could not read `", arg, "`: ", conditionMessage(error),
           call. = FALSE)
    }
  )
  if (!is.data.frame(out) || nrow(out) < 1 || ncol(out) < 1) {
    stop("`", arg, "` must contain at least one row and one treatment column.",
         call. = FALSE)
  }
  names(out) = .uaf_squish_text(names(out))
  if (any(is.na(names(out)) | names(out) == "")) {
    stop("`", arg, "` contains an empty treatment column name.",
         call. = FALSE)
  }
  out
}

.ainsect_treatment_aliases = function(treatment_aliases = NULL) {
  defaults = c(
    "Melissa " = "Melissa",
    "Melissa" = "Melissa",
    "Clary Sage " = "Clary Sage",
    "Clary Sage" = "Clary Sage",
    "Clary Sagae" = "Clary Sage",
    "Labdanum " = "Labdanum",
    "Labdanum" = "Labdanum",
    "Nootka " = "Nootka",
    "Nootka" = "Nootka",
    "Vanilla " = "Vanilla",
    "Vanilla" = "Vanilla"
  )
  if (is.null(treatment_aliases)) return(defaults)
  if (is.null(names(treatment_aliases)) ||
      any(names(treatment_aliases) == "")) {
    stop("`treatment_aliases` must be a named character vector where names ",
         "are aliases and values are canonical treatment names.", call. = FALSE)
  }
  c(defaults, as.character(treatment_aliases))
}

.ainsect_treatment_key = function(x) {
  x = .uaf_squish_text(x)
  x = tolower(x)
  x = gsub("[^a-z0-9]+", " ", x)
  .uaf_squish_text(x)
}

.ainsect_canonical_treatment = function(x, aliases) {
  cleaned = .uaf_squish_text(x)
  alias_keys = .ainsect_treatment_key(names(aliases))
  alias_values = .uaf_squish_text(unname(aliases))
  lookup = stats::setNames(alias_values, alias_keys)
  key = .ainsect_treatment_key(cleaned)
  out = cleaned
  matched = key %in% names(lookup)
  out[matched] = unname(lookup[key[matched]])
  .uaf_squish_text(out)
}

.ainsect_treatment_map = function(table, table_name, aliases) {
  raw_names = names(table)
  treatment_name = .ainsect_canonical_treatment(raw_names, aliases)
  key = .ainsect_treatment_key(treatment_name)
  out = data.frame(
    table = table_name,
    column_name = raw_names,
    treatment_name = treatment_name,
    treatment_key = key,
    stringsAsFactors = FALSE
  )
  duplicated_keys = unique(out$treatment_key[duplicated(out$treatment_key)])
  if (length(duplicated_keys) > 0) {
    examples = vapply(duplicated_keys, function(key_value) {
      paste(out$column_name[out$treatment_key == key_value], collapse = ", ")
    }, character(1))
    stop("Duplicate treatment columns after cleanup in ", table_name, ": ",
         paste(examples, collapse = "; "), call. = FALSE)
  }
  out
}

.ainsect_align_treatments = function(identity_wide, quant_wide, aliases) {
  id_map = .ainsect_treatment_map(identity_wide, "identity table", aliases)
  quant_map = .ainsect_treatment_map(quant_wide, "abundance table", aliases)
  missing_quant = setdiff(id_map$treatment_key, quant_map$treatment_key)
  missing_identity = setdiff(quant_map$treatment_key, id_map$treatment_key)
  if (length(missing_quant) > 0 || length(missing_identity) > 0) {
    id_missing = id_map$treatment_name[id_map$treatment_key %in% missing_quant]
    quant_missing = quant_map$treatment_name[
      quant_map$treatment_key %in% missing_identity
    ]
    stop("Treatment columns could not be aligned. Missing from abundance ",
         "table: ", paste(id_missing, collapse = ", "),
         ". Missing from identity table: ",
         paste(quant_missing, collapse = ", "), ".", call. = FALSE)
  }

  quant_index = match(id_map$treatment_key, quant_map$treatment_key)
  treatment_ids = .ainsect_make_unique_ids(
    .ainsect_sanitize_id(id_map$treatment_name),
    "treatment"
  )
  data.frame(
    treatment_id = treatment_ids,
    treatment_name = id_map$treatment_name,
    identity_column = id_map$column_name,
    abundance_column = quant_map$column_name[quant_index],
    treatment_key = id_map$treatment_key,
    stringsAsFactors = FALSE
  )
}

.ainsect_long_profile = function(identity_wide, quant_wide, alignment,
                                 chem_id_csv, chem_quant_csv, uafR_run_id) {
  rows = list()
  for (i in seq_len(nrow(alignment))) {
    item = alignment[i, , drop = FALSE]
    compound_name = .uaf_squish_text(identity_wide[[item$identity_column]])
    absolute_abundance = .ainsect_parse_abundance(
      quant_wide[[item$abundance_column]]
    )
    rows[[length(rows) + 1]] = data.frame(
      treatment_id = item$treatment_id,
      treatment_name = item$treatment_name,
      rank = seq_len(nrow(identity_wide)),
      compound_name = compound_name,
      absolute_abundance = absolute_abundance,
      abundance_units = "gcms_relative_area_or_input_units",
      uafR_run_id = uafR_run_id,
      source_file = paste(basename(chem_id_csv), basename(chem_quant_csv),
                          sep = ";"),
      stringsAsFactors = FALSE
    )
  }
  raw_rows = do.call(rbind, rows)
  invalid_abundance = is.na(raw_rows$absolute_abundance)
  missing_compound = is.na(raw_rows$compound_name) | raw_rows$compound_name == ""
  negative_abundance = !is.na(raw_rows$absolute_abundance) &
    raw_rows$absolute_abundance < 0
  if (any(negative_abundance)) {
    stop("Negative abundance values were found. Relative abundance cannot ",
         "be computed from negative GC-MS abundance input.", call. = FALSE)
  }
  valid_rows = raw_rows[!invalid_abundance & !missing_compound, , drop = FALSE]
  row.names(valid_rows) = NULL
  list(
    raw_rows = raw_rows,
    valid_rows = valid_rows,
    invalid_abundance_rows = sum(invalid_abundance),
    missing_compound_rows = sum(missing_compound)
  )
}

.ainsect_parse_abundance = function(x) {
  x = .uaf_squish_text(x)
  x = gsub(",", "", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

.ainsect_compound_resolution = function(compound_names, pubchem_profile,
                                        source_label) {
  identity = .ainsect_profile_table(pubchem_profile, "identity")
  properties = .ainsect_profile_table(pubchem_profile, "properties")
  records = lapply(compound_names, function(compound_name) {
    id_row = .ainsect_match_profile_row(identity, "Query", compound_name)
    cid = .ainsect_first_existing(id_row, "CID")
    if (is.na(cid) || cid == "") {
      prop_row = .ainsect_match_profile_row(properties, "Query", compound_name)
    } else {
      prop_row = .ainsect_match_profile_row(properties, "CID", cid)
      if (nrow(prop_row) < 1) {
        prop_row = .ainsect_match_profile_row(properties, "Query",
                                              compound_name)
      }
    }
    cid = .uaf_first_non_empty_text(
      .ainsect_first_existing(prop_row, "CID"),
      cid
    )
    inchikey = .ainsect_first_existing(prop_row, "InChIKey")
    molecular_formula = .ainsect_first_existing(prop_row, "MolecularFormula")
    smiles = .uaf_first_non_empty_text(
      .ainsect_first_existing(prop_row, "IsomericSMILES"),
      .ainsect_first_existing(prop_row, "CanonicalSMILES"),
      .ainsect_first_existing(prop_row, "SMILES"),
      .ainsect_first_existing(prop_row, "ConnectivitySMILES")
    )
    compound_id = .ainsect_compound_id(compound_name, inchikey, cid)
    resolved_cid = !is.na(cid) && cid != ""
    resolved_smiles = !is.na(smiles) && smiles != ""
    note = .ainsect_resolution_note(compound_id, inchikey, cid,
                                    resolved_cid, resolved_smiles)
    data.frame(
      compound_id = compound_id,
      compound_name = compound_name,
      smiles = smiles,
      inchikey = inchikey,
      pubchem_cid = cid,
      molecular_formula = molecular_formula,
      source = source_label,
      notes = note,
      has_usable_smiles = resolved_smiles,
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, records)
  row.names(out) = NULL
  list(records = out, identity = identity, properties = properties)
}

.ainsect_call_profile_fun = function(profile_fun, compounds, profile, cache_dir,
                                     throttle) {
  args = list(
    compounds = compounds,
    profile = profile,
    cache = TRUE,
    cache_dir = cache_dir,
    throttle = throttle,
    assay_detail_limit = 0,
    include_annotations = FALSE
  )
  formals_names = names(formals(profile_fun))
  if (!is.null(formals_names) && !"..." %in% formals_names) {
    args = args[names(args) %in% formals_names]
  }
  do.call(profile_fun, args)
}

.ainsect_profile_table = function(profile, table_name) {
  table = tryCatch(profile[[table_name]], error = function(error) NULL)
  if (is.data.frame(table)) return(table)
  .uaf_empty_table(character())
}

.ainsect_match_profile_row = function(table, column, value) {
  if (!is.data.frame(table) || nrow(table) < 1 || !column %in% names(table)) {
    return(table[0, , drop = FALSE])
  }
  candidate = as.character(table[[column]])
  index = which(candidate == as.character(value))
  if (length(index) < 1 && column == "Query") {
    index = which(tolower(candidate) == tolower(as.character(value)))
  }
  if (length(index) < 1) return(table[0, , drop = FALSE])
  table[index[[1]], , drop = FALSE]
}

.ainsect_first_existing = function(table, column) {
  if (!is.data.frame(table) || nrow(table) < 1 || !column %in% names(table)) {
    return(NA_character_)
  }
  .uaf_first_non_empty_text(table[[column]])
}

.ainsect_resolution_note = function(compound_id, inchikey, cid, resolved_cid,
                                    resolved_smiles) {
  id_source = if (!is.na(inchikey) && inchikey != "") {
    "compound_id_from_inchikey"
  } else if (!is.na(cid) && cid != "") {
    "compound_id_from_pubchem_cid"
  } else {
    "compound_id_from_query_name"
  }
  status = if (isTRUE(resolved_smiles)) {
    "resolved_with_pubchem_smiles"
  } else if (isTRUE(resolved_cid)) {
    "pubchem_cid_resolved_but_no_usable_smiles"
  } else {
    "pubchem_query_unresolved"
  }
  paste(status, id_source, paste0("compound_id=", compound_id), sep = "; ")
}

.ainsect_compound_id = function(compound_name, inchikey, cid) {
  inchikey = .uaf_first_non_empty_text(inchikey)
  cid = .uaf_first_non_empty_text(cid)
  if (!is.na(inchikey) && inchikey != "") {
    return(.ainsect_sanitize_id(inchikey))
  }
  if (!is.na(cid) && cid != "") {
    cid = sub("\\.0$", "", cid)
    return(paste0("cid_", .ainsect_sanitize_id(cid)))
  }
  paste0("unresolved_", .ainsect_sanitize_id(compound_name))
}

.ainsect_sanitize_id = function(x) {
  x = .uaf_squish_text(x)
  x = gsub("[^A-Za-z0-9]+", "_", x)
  x = gsub("^_+|_+$", "", x)
  x = ifelse(is.na(x) | x == "", "unknown", x)
  x
}

.ainsect_make_unique_ids = function(ids, fallback_prefix) {
  ids = .ainsect_sanitize_id(ids)
  ids[ids == "" | is.na(ids)] = fallback_prefix
  make.unique(ids, sep = "_")
}

.ainsect_compounds_for_export = function(compounds_all) {
  usable = compounds_all[compounds_all$has_usable_smiles %in% TRUE,
                         , drop = FALSE]
  if (nrow(usable) < 1) {
    return(.uaf_empty_table(.ainsect_compound_cols()))
  }
  grouped = split(usable, usable$compound_id)
  rows = lapply(grouped, function(group) {
    notes = unique(.uaf_non_empty(group$notes))
    input_names = unique(.uaf_non_empty(group$compound_name))
    if (length(input_names) > 1) {
      notes = c(notes, paste0("input_names=",
                              paste(input_names, collapse = " | ")))
    }
    data.frame(
      compound_id = group$compound_id[[1]],
      compound_name = .uaf_first_non_empty_text(group$compound_name),
      smiles = .uaf_first_non_empty_text(group$smiles),
      inchikey = .uaf_first_non_empty_text(group$inchikey),
      pubchem_cid = .uaf_first_non_empty_text(group$pubchem_cid),
      molecular_formula = .uaf_first_non_empty_text(group$molecular_formula),
      source = .uaf_first_non_empty_text(group$source),
      notes = paste(notes, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  out = out[order(tolower(out$compound_name)), .ainsect_compound_cols(),
            drop = FALSE]
  row.names(out) = NULL
  out
}

.ainsect_compound_cols = function() {
  c("compound_id", "compound_name", "smiles", "inchikey", "pubchem_cid",
    "molecular_formula", "source", "notes")
}

.ainsect_unresolved_table = function(compounds_all, valid_rows) {
  unresolved = compounds_all[!compounds_all$has_usable_smiles %in% TRUE,
                             , drop = FALSE]
  cols = c("compound_name", "proposed_compound_id", "pubchem_cid",
           "molecular_formula", "inchikey", "smiles", "unresolved_reason",
           "notes", "treatment_count", "absolute_abundance_total")
  if (nrow(unresolved) < 1) return(.uaf_empty_table(cols))
  abundance = stats::aggregate(
    absolute_abundance ~ compound_name,
    data = valid_rows,
    FUN = sum
  )
  treatment_count = stats::aggregate(
    treatment_id ~ compound_name,
    data = unique(valid_rows[c("compound_name", "treatment_id")]),
    FUN = length
  )
  out = merge(unresolved, abundance, by = "compound_name", all.x = TRUE)
  out = merge(out, treatment_count, by = "compound_name", all.x = TRUE)
  out$unresolved_reason = ifelse(
    is.na(out$pubchem_cid) | out$pubchem_cid == "",
    "pubchem_query_unresolved",
    "pubchem_cid_resolved_but_no_usable_smiles"
  )
  out$proposed_compound_id = out$compound_id
  out$treatment_count = out$treatment_id
  out$absolute_abundance_total = out$absolute_abundance
  out = out[order(tolower(out$compound_name)), cols, drop = FALSE]
  row.names(out) = NULL
  out
}

.ainsect_abundance_for_export = function(valid_rows, compound_records,
                                         min_relative_abundance,
                                         uafR_run_id) {
  key = compound_records[, c("compound_name", "compound_id",
                             "has_usable_smiles"), drop = FALSE]
  joined = merge(valid_rows, key, by = "compound_name", all.x = TRUE,
                 sort = FALSE)
  joined = joined[joined$has_usable_smiles %in% TRUE, , drop = FALSE]
  if (nrow(joined) < 1) {
    return(.uaf_empty_table(.ainsect_abundance_cols()))
  }

  aggregate_cols = c("treatment_id", "treatment_name", "compound_id",
                     "abundance_units", "source_file")
  abundance = stats::aggregate(
    absolute_abundance ~ treatment_id + treatment_name + compound_id +
      abundance_units + source_file,
    data = joined[, c(aggregate_cols, "absolute_abundance"), drop = FALSE],
    FUN = sum
  )
  abundance = .ainsect_add_relative_abundance(abundance)
  abundance = abundance[
    abundance$relative_abundance >= min_relative_abundance,
    , drop = FALSE
  ]
  if (nrow(abundance) < 1) {
    return(.uaf_empty_table(.ainsect_abundance_cols()))
  }
  abundance = .ainsect_add_relative_abundance(abundance)
  abundance$uafR_run_id = uafR_run_id
  abundance = abundance[order(abundance$treatment_name,
                              -abundance$absolute_abundance,
                              abundance$compound_id),
                        .ainsect_abundance_cols(), drop = FALSE]
  row.names(abundance) = NULL
  abundance
}

.ainsect_add_relative_abundance = function(abundance) {
  totals = stats::aggregate(
    absolute_abundance ~ treatment_id,
    data = abundance,
    FUN = sum
  )
  names(totals)[names(totals) == "absolute_abundance"] = "treatment_total"
  out = merge(abundance, totals, by = "treatment_id", all.x = TRUE,
              sort = FALSE)
  out = out[out$treatment_total > 0, , drop = FALSE]
  out$relative_abundance = out$absolute_abundance / out$treatment_total
  out$treatment_total = NULL
  out
}

.ainsect_abundance_cols = function() {
  c("treatment_id", "treatment_name", "compound_id", "relative_abundance",
    "absolute_abundance", "abundance_units", "uafR_run_id", "source_file")
}

.ainsect_write_audit_table = function(table, path) {
  if (!is.data.frame(table)) table = .uaf_empty_table(character())
  utils::write.csv(as.data.frame(table, stringsAsFactors = FALSE),
                   path, row.names = FALSE, na = "")
}

.ainsect_export_summary = function(uafR_run_id, source_label, profile,
                                   throttle, min_relative_abundance,
                                   chem_id_csv, chem_quant_csv, out_dir,
                                   cache_dir, alignment, parsed,
                                   compounds_all, compounds, unresolved,
                                   abundance, paths) {
  relative_sums = if (nrow(abundance) > 0) {
    stats::aggregate(relative_abundance ~ treatment_id + treatment_name,
                     data = abundance, FUN = sum)
  } else {
    .uaf_empty_table(c("treatment_id", "treatment_name",
                       "relative_abundance"))
  }
  unresolved_abundance = if (nrow(unresolved) > 0) {
    sum(suppressWarnings(as.numeric(unresolved$absolute_abundance_total)),
        na.rm = TRUE)
  } else {
    0
  }
  list(
    uafR_run_id = uafR_run_id,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    source_label = source_label,
    profile = profile,
    throttle = throttle,
    min_relative_abundance = min_relative_abundance,
    inputs = list(
      chem_id_csv = chem_id_csv,
      chem_quant_csv = chem_quant_csv,
      cache_dir = cache_dir
    ),
    outputs = paths,
    treatment_count = nrow(alignment),
    treatment_names = alignment$treatment_name,
    row_counts = list(
      raw_pair_rows = nrow(parsed$raw_rows),
      valid_pair_rows = nrow(parsed$valid_rows),
      missing_compound_rows = parsed$missing_compound_rows,
      invalid_abundance_rows = parsed$invalid_abundance_rows,
      unique_compounds_queried = nrow(compounds_all),
      compound_queries_with_pubchem_cid = sum(!is.na(compounds_all$pubchem_cid) &
                                                compounds_all$pubchem_cid != ""),
      compound_queries_with_usable_smiles = sum(compounds_all$has_usable_smiles %in%
                                                  TRUE),
      exported_unique_compounds = nrow(compounds),
      unresolved_compounds = nrow(unresolved),
      abundance_rows_written = nrow(abundance)
    ),
    quality = list(
      relative_abundance_sum_min = if (nrow(relative_sums) > 0) {
        min(relative_sums$relative_abundance)
      } else {
        NA_real_
      },
      relative_abundance_sum_max = if (nrow(relative_sums) > 0) {
        max(relative_sums$relative_abundance)
      } else {
        NA_real_
      },
      unresolved_absolute_abundance_total = unresolved_abundance,
      abundance_compound_ids_all_in_compounds = all(
        unique(abundance$compound_id) %in% compounds$compound_id
      )
    )
  )
}
