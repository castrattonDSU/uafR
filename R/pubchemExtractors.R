.pubchem_extract_profiles = function(pubchem) {
  list(
    PubChemIdentifierProfile = .pubchem_extract_identifier_profile(pubchem),
    SafetyProfile = .pubchem_extract_safety_profile(pubchem),
    FEMAProfile = .pubchem_extract_fema_profile(pubchem),
    FDA_SPL_Profile = .pubchem_extract_fda_spl_profile(pubchem),
    LOTUSProfile = .pubchem_extract_lotus_profile(pubchem),
    PubChemClassificationProfile = .pubchem_extract_classification_profile(pubchem),
    MeSHProfile = .pubchem_extract_mesh_profile(pubchem),
    LiteratureProfile = .pubchem_extract_literature_profile(pubchem)
  )
}

.pubchem_profile_annotations = function(pubchem) {
  if (is.null(pubchem) ||
      !is.list(pubchem) ||
      is.null(pubchem$annotations) ||
      !is.data.frame(pubchem$annotations)) {
    return(.pubchem_empty_annotation_table())
  }
  pubchem$annotations
}

.pubchem_dedupe_profile = function(table, key_cols,
                                   collapse_cols = c("RawValue", "Source",
                                                     "SourceURL", "PubChemURL"),
                                   bool_cols = character(),
                                   fallback_cols = c("RawValue")) {
  if (!is.data.frame(table) || nrow(table) < 1) return(table)

  key_cols = intersect(key_cols, colnames(table))
  if (length(key_cols) < 1) {
    out = unique(table)
    row.names(out) = NULL
    return(out)
  }

  key = .pubchem_row_key(table, key_cols, fallback_cols)
  group_keys = unique(key)
  rows = lapply(group_keys, function(group_key) {
    group = table[key == group_key, , drop = FALSE]
    row = group[1, , drop = FALSE]

    for (col in intersect(collapse_cols, colnames(group))) {
      row[[col]] = .pubchem_collapse(group[[col]])
    }

    for (col in intersect(bool_cols, colnames(group))) {
      row[[col]] = any(group[[col]] %in% TRUE, na.rm = TRUE)
    }

    for (col in setdiff(colnames(group),
                        c(key_cols, collapse_cols, bool_cols))) {
      if (is.numeric(group[[col]])) {
        vals = group[[col]][!is.na(group[[col]])]
        if (length(vals) > 0 && is.na(row[[col]][[1]])) row[[col]] = vals[[1]]
      } else {
        vals = .uaf_non_empty(group[[col]])
        if (length(vals) > 0 &&
            (is.na(row[[col]][[1]]) || row[[col]][[1]] == "")) {
          row[[col]] = vals[[1]]
        }
      }
    }

    row
  })

  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.pubchem_dedupe_annotations = function(table) {
  if (!is.data.frame(table) || nrow(table) < 1) return(table)
  .pubchem_dedupe_profile(
    table = table,
    key_cols = c("Query", "CID", "Name", "CleanValue", "ValueNumeric",
                 "UnitClean", "Source"),
    collapse_cols = c("Heading", "HeadingPath", "Value", "Unit", "SourceURL",
                      "PubChemURL", "MarkupText", "MarkupURL", "MarkupExtra"),
    fallback_cols = c("HeadingPath", "Name", "CleanValue", "Value")
  )
}

.pubchem_row_key = function(table, key_cols, fallback_cols = character()) {
  key_frame = as.data.frame(lapply(table[key_cols], .pubchem_key_value),
                            stringsAsFactors = FALSE)
  key = apply(key_frame, 1, paste, collapse = "\r")

  semantic_cols = setdiff(key_cols, c("Query", "CID"))
  if (length(semantic_cols) > 0) {
    semantic_frame = key_frame[semantic_cols]
    has_semantic_value = apply(semantic_frame, 1, function(row) {
      any(!is.na(row) & row != "" & row != "<NA>")
    })
  } else {
    has_semantic_value = rep(FALSE, nrow(table))
  }

  fallback_cols = intersect(fallback_cols, colnames(table))
  if (length(fallback_cols) > 0 && any(!has_semantic_value)) {
    fallback_frame = as.data.frame(lapply(table[fallback_cols],
                                          .pubchem_key_value),
                                  stringsAsFactors = FALSE)
    fallback_key = apply(fallback_frame, 1, paste, collapse = "\r")
    key[!has_semantic_value] = paste(key[!has_semantic_value],
                                     fallback_key[!has_semantic_value],
                                     sep = "\r")
  }

  key
}

.pubchem_key_value = function(x) {
  x = .uaf_squish_text(x)
  x[is.na(x)] = "<NA>"
  x
}

.pubchem_empty_annotation_table = function() {
  .uaf_empty_table(.pubchem_annotation_cols())
}

.pubchem_annotation_cols = function() {
  c("Query", "CID", "Heading", "HeadingPath", "Name",
    "Value", "CleanValue", "ValueNumeric", "Unit",
    "UnitClean", "MarkupText", "MarkupURL", "MarkupExtra",
    "Source", "SourceURL", "PubChemURL")
}

.pubchem_annotation_text = function(annotation) {
  paste(annotation$Heading,
        annotation$HeadingPath,
        annotation$Name,
        annotation$CleanValue,
        annotation$Source,
        sep = " | ")
}

.pubchem_filter_profile_annotations = function(pubchem, pattern) {
  annotations = .pubchem_profile_annotations(pubchem)
  if (nrow(annotations) < 1) return(annotations)
  text = paste(annotations$Heading,
               annotations$HeadingPath,
               annotations$Name,
               annotations$CleanValue,
               annotations$Source)
  annotations[grepl(pattern, text, ignore.case = TRUE), , drop = FALSE]
}

.pubchem_extract_identifier_profile = function(pubchem) {
  cols = c("Query", "CID", "IdentifierType", "Identifier", "SourceField",
           "RawValue", "Source", "SourceURL", "PubChemURL")
  rows = list()
  annotations = .pubchem_profile_annotations(pubchem)

  if (nrow(annotations) > 0) {
    for (i in seq_len(nrow(annotations))) {
      annotation = annotations[i, , drop = FALSE]
      identifiers = .uaf_extract_identifiers(
        paste(annotation$Name, annotation$CleanValue, annotation$Source)
      )
      if (nrow(identifiers) < 1) next
      for (j in seq_len(nrow(identifiers))) {
        rows[[length(rows) + 1]] = data.frame(
          Query = annotation$Query,
          CID = annotation$CID,
          IdentifierType = identifiers$IdentifierType[[j]],
          Identifier = identifiers$Identifier[[j]],
          SourceField = .uaf_first_non_empty_text(annotation$HeadingPath,
                                                  annotation$Name),
          RawValue = annotation$CleanValue,
          Source = annotation$Source,
          SourceURL = annotation$SourceURL,
          PubChemURL = annotation$PubChemURL,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (is.data.frame(pubchem$synonyms) && nrow(pubchem$synonyms) > 0) {
    for (i in seq_len(nrow(pubchem$synonyms))) {
      synonym = pubchem$synonyms[i, , drop = FALSE]
      identifiers = .uaf_extract_identifiers(synonym$Synonym)
      if (nrow(identifiers) < 1) next
      for (j in seq_len(nrow(identifiers))) {
        rows[[length(rows) + 1]] = data.frame(
          Query = synonym$Query,
          CID = synonym$CID,
          IdentifierType = identifiers$IdentifierType[[j]],
          Identifier = identifiers$Identifier[[j]],
          SourceField = "synonym",
          RawValue = synonym$Synonym,
          Source = "PubChem synonyms",
          SourceURL = synonym$SourceURL,
          PubChemURL = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) < 1) return(.uaf_empty_table(cols))
  .pubchem_dedupe_profile(
    table = do.call(rbind, rows),
    key_cols = c("Query", "CID", "IdentifierType", "Identifier"),
    collapse_cols = c("SourceField", "RawValue", "Source", "SourceURL",
                      "PubChemURL"),
    fallback_cols = c("RawValue", "SourceField")
  )
}

.pubchem_extract_safety_profile = function(pubchem) {
  cols = c("Query", "CID", "SignalWord", "HazardCode", "PrecautionCode",
           "HazardClass", "HazardStatement", "RawValue", "Source",
           "SourceURL", "PubChemURL")
  annotations = .pubchem_filter_profile_annotations(
    pubchem,
    "Safety|Hazard|GHS|Toxicity|Toxicological|Fire|Exposure|Handling|Warning|Danger|H[0-9]{3}|P[0-9]{3}"
  )
  if (nrow(annotations) < 1) return(.uaf_empty_table(cols))

  rows = list()
  for (i in seq_len(nrow(annotations))) {
    annotation = annotations[i, , drop = FALSE]
    text = .pubchem_annotation_text(annotation)
    hazard_codes = .uaf_extract_pattern(text, "\\bH[0-9]{3}[A-Za-z]?\\b")
    precaution_codes = .uaf_extract_pattern(text, "\\bP[0-9]{3}[A-Za-z]?\\b")
    signal_word = .uaf_first_non_empty_text(.pubchem_signal_word(text))
    hazard_class = if (grepl("GHS|Hazard|Toxicity|Irritation|Flammable|Corrosion|Carcinogen|Aquatic",
                             text, ignore.case = TRUE)) {
      .uaf_first_non_empty_text(annotation$Name, .pubchem_last_path(annotation$HeadingPath))
    } else {
      NA_character_
    }
    hazard_statement = if (grepl("H[0-9]{3}|hazard|toxic|irritation|flammable|warning|danger|causes|harmful",
                                 text, ignore.case = TRUE)) {
      annotation$CleanValue
    } else {
      NA_character_
    }
    rows[[length(rows) + 1]] = data.frame(
      Query = annotation$Query,
      CID = annotation$CID,
      SignalWord = signal_word,
      HazardCode = .pubchem_collapse(hazard_codes),
      PrecautionCode = .pubchem_collapse(precaution_codes),
      HazardClass = hazard_class,
      HazardStatement = hazard_statement,
      RawValue = annotation$CleanValue,
      Source = annotation$Source,
      SourceURL = annotation$SourceURL,
      PubChemURL = annotation$PubChemURL,
      stringsAsFactors = FALSE
    )
  }
  .pubchem_dedupe_profile(
    table = do.call(rbind, rows),
    key_cols = c("Query", "CID", "SignalWord", "HazardCode",
                 "PrecautionCode", "HazardStatement"),
    collapse_cols = c("HazardClass", "RawValue", "Source", "SourceURL",
                      "PubChemURL"),
    fallback_cols = c("RawValue", "HazardClass")
  )
}

.pubchem_extract_fema_profile = function(pubchem) {
  cols = c("Query", "CID", "FEMANumber", "GRASStatus", "JECFANumber",
           "DescriptorTerms", "FlavorTerms", "OdorTerms", "RawValue",
           "Source", "SourceURL", "PubChemURL")
  annotations = .pubchem_filter_profile_annotations(
    pubchem,
    "FEMA|Flavor and Extract|Flavor|Odor|Taste|GRAS|JECFA"
  )
  if (nrow(annotations) < 1) return(.uaf_empty_table(cols))

  rows = lapply(seq_len(nrow(annotations)), function(i) {
    annotation = annotations[i, , drop = FALSE]
    text = .pubchem_annotation_text(annotation)
    data.frame(
      Query = annotation$Query,
      CID = annotation$CID,
      FEMANumber = .pubchem_collapse(.pubchem_extract_fema_numbers(text)),
      GRASStatus = .pubchem_gras_status(text),
      JECFANumber = .pubchem_collapse(.pubchem_extract_jecfa_numbers(text)),
      DescriptorTerms = .pubchem_collapse(.pubchem_descriptor_terms(text)),
      FlavorTerms = .pubchem_collapse(.pubchem_descriptor_terms(text, "flavor")),
      OdorTerms = .pubchem_collapse(.pubchem_descriptor_terms(text, "odor")),
      RawValue = annotation$CleanValue,
      Source = annotation$Source,
      SourceURL = annotation$SourceURL,
      PubChemURL = annotation$PubChemURL,
      stringsAsFactors = FALSE
    )
  })
  .pubchem_dedupe_profile(
    table = do.call(rbind, rows),
    key_cols = c("Query", "CID", "FEMANumber", "GRASStatus", "JECFANumber",
                 "DescriptorTerms", "FlavorTerms", "OdorTerms"),
    collapse_cols = c("RawValue", "Source", "SourceURL", "PubChemURL"),
    fallback_cols = c("RawValue")
  )
}

.pubchem_extract_fda_spl_profile = function(pubchem) {
  cols = c("Query", "CID", "LabelSection", "ActiveIngredient",
           "PharmacologicClass", "Route", "DosageForm", "HasBoxedWarning",
           "WarningTerms", "RawValue", "Source", "SourceURL", "PubChemURL")
  annotations = .pubchem_filter_profile_annotations(
    pubchem,
    "FDA/SPL|SPL Indexing|Drug and Medication|Therapeutic Uses|Drug|Medication|Clinical|Pharmacologic|Indication|Contraindication|Adverse|Route|Dosage|Ingredient|Boxed Warning"
  )
  if (nrow(annotations) < 1) return(.uaf_empty_table(cols))

  rows = lapply(seq_len(nrow(annotations)), function(i) {
    annotation = annotations[i, , drop = FALSE]
    text = .pubchem_annotation_text(annotation)
    data.frame(
      Query = annotation$Query,
      CID = annotation$CID,
      LabelSection = .uaf_first_non_empty_text(.pubchem_last_path(annotation$HeadingPath),
                                               annotation$Name),
      ActiveIngredient = if (grepl("active ingredient|ingredient",
                                   text, ignore.case = TRUE)) annotation$CleanValue else NA_character_,
      PharmacologicClass = if (grepl("pharmacologic|therapeutic|class",
                                     text, ignore.case = TRUE)) annotation$CleanValue else NA_character_,
      Route = .pubchem_collapse(.pubchem_extract_routes(text)),
      DosageForm = .pubchem_collapse(.pubchem_extract_dosage_forms(text)),
      HasBoxedWarning = grepl("boxed warning|black box", text,
                              ignore.case = TRUE),
      WarningTerms = .pubchem_collapse(.pubchem_warning_terms(text)),
      RawValue = annotation$CleanValue,
      Source = annotation$Source,
      SourceURL = annotation$SourceURL,
      PubChemURL = annotation$PubChemURL,
      stringsAsFactors = FALSE
    )
  })
  .pubchem_dedupe_profile(
    table = do.call(rbind, rows),
    key_cols = c("Query", "CID", "ActiveIngredient", "PharmacologicClass",
                 "Route", "DosageForm", "HasBoxedWarning", "WarningTerms"),
    collapse_cols = c("LabelSection", "RawValue", "Source", "SourceURL",
                      "PubChemURL"),
    bool_cols = "HasBoxedWarning",
    fallback_cols = c("RawValue", "LabelSection")
  )
}

.pubchem_extract_lotus_profile = function(pubchem) {
  cols = c("Query", "CID", "LOTUS_ID", "TaxonomyID", "Organism",
           "CommonName", "TaxonomyRank", "Taxonomy", "TaxonomyLineage",
           "NaturalProductClass", "Occurrence", "ReferenceID", "RawValue",
           "Source", "SourceURL", "PubChemURL")
  annotations = .pubchem_filter_profile_annotations(
    pubchem,
    "LOTUS|natural product|organism|taxonomy|taxon|occurrence|species|genus|family"
  )
  taxonomy = if (is.data.frame(pubchem$taxonomy)) pubchem$taxonomy else NULL
  if (nrow(annotations) < 1 &&
      (is.null(taxonomy) || nrow(taxonomy) < 1)) {
    return(.uaf_empty_table(cols))
  }

  rows = list()
  if (nrow(annotations) > 0) rows = lapply(seq_len(nrow(annotations)), function(i) {
    annotation = annotations[i, , drop = FALSE]
    field_context = paste(annotation$HeadingPath,
                          annotation$Name,
                          sep = " | ")
    text = .pubchem_annotation_text(annotation)
    value = annotation$CleanValue
    markup_text = if ("MarkupText" %in% colnames(annotation)) {
      annotation$MarkupText
    } else {
      NA_character_
    }
    markup_url = if ("MarkupURL" %in% colnames(annotation)) {
      annotation$MarkupURL
    } else {
      NA_character_
    }
    organisms = .lotus_extract_organisms(value, field_context, markup_text)
    data.frame(
      Query = annotation$Query,
      CID = annotation$CID,
      LOTUS_ID = .pubchem_collapse(.pubchem_extract_lotus_ids(text)),
      TaxonomyID = .pubchem_collapse(.pubchem_extract_taxonomy_ids(markup_url)),
      Organism = .pubchem_collapse(organisms),
      CommonName = NA_character_,
      TaxonomyRank = NA_character_,
      Taxonomy = .lotus_extract_taxonomy(value, field_context),
      TaxonomyLineage = NA_character_,
      NaturalProductClass = .pubchem_collapse(.pubchem_np_class_terms(text)),
      Occurrence = if (.lotus_context_is(field_context,
                                         c("occurrence", "biological source",
                                           "source organism", "plant part",
                                           "record description"))) {
        value
      } else {
        NA_character_
      },
      ReferenceID = .pubchem_collapse(.uaf_extract_identifiers(text)$Identifier),
      RawValue = value,
      Source = annotation$Source,
      SourceURL = annotation$SourceURL,
      PubChemURL = annotation$PubChemURL,
      stringsAsFactors = FALSE
    )
  })

  if (!is.null(taxonomy) && nrow(taxonomy) > 0) {
    for (i in seq_len(nrow(taxonomy))) {
      taxon = taxonomy[i, , drop = FALSE]
      if (!grepl("LOTUS", paste(taxon$Source, taxon$PubChemURL),
                 ignore.case = TRUE)) next
      rows[[length(rows) + 1]] = data.frame(
        Query = taxon$Query,
        CID = taxon$CID,
        LOTUS_ID = NA_character_,
        TaxonomyID = taxon$TaxonomyID,
        Organism = taxon$Organism,
        CommonName = taxon$CommonName,
        TaxonomyRank = taxon$Rank,
        Taxonomy = taxon$Lineage,
        TaxonomyLineage = taxon$Lineage,
        NaturalProductClass = NA_character_,
        Occurrence = taxon$SourceAnnotation,
        ReferenceID = taxon$TaxonomyID,
        RawValue = .pubchem_collapse(c(taxon$Organism, taxon$Lineage)),
        Source = taxon$Source,
        SourceURL = taxon$TaxonomyURL,
        PubChemURL = taxon$PubChemURL,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) < 1) return(.uaf_empty_table(cols))
  table = do.call(rbind, rows)
  semantic_cols = c("LOTUS_ID", "TaxonomyID", "Organism", "Taxonomy",
                    "TaxonomyLineage", "NaturalProductClass", "ReferenceID")
  keep = apply(table[semantic_cols], 1, function(row) {
    length(.uaf_non_empty(row)) > 0
  })
  table = table[keep, , drop = FALSE]
  compound_key = paste(.pubchem_key_value(table$Query),
                       .pubchem_key_value(table$CID),
                       sep = "\r")
  taxon_specific = !is.na(table$TaxonomyID) & table$TaxonomyID != "" &
    !grepl(";", table$TaxonomyID, fixed = TRUE) &
    !is.na(table$TaxonomyLineage) & table$TaxonomyLineage != ""
  compounds_with_taxa = unique(compound_key[taxon_specific])
  aggregate_taxon_row = compound_key %in% compounds_with_taxa &
    !is.na(table$TaxonomyID) & grepl(";", table$TaxonomyID, fixed = TRUE) &
    (is.na(table$TaxonomyLineage) | table$TaxonomyLineage == "")
  table = table[!aggregate_taxon_row, , drop = FALSE]
  if (nrow(table) < 1) return(.uaf_empty_table(cols))
  .pubchem_dedupe_profile(
    table = table,
    key_cols = c("Query", "CID", "LOTUS_ID", "TaxonomyID", "Organism",
                 "Taxonomy", "NaturalProductClass", "Occurrence",
                 "ReferenceID"),
    collapse_cols = c("RawValue", "Source", "SourceURL", "PubChemURL"),
    fallback_cols = c("RawValue")
  )
}

.pubchem_extract_classification_profile = function(pubchem) {
  cols = .pubchem_classification_cols()
  if (is.null(pubchem) ||
      !is.list(pubchem) ||
      is.null(pubchem$classifications) ||
      !is.data.frame(pubchem$classifications) ||
      nrow(pubchem$classifications) < 1) {
    return(.uaf_empty_table(cols))
  }

  table = pubchem$classifications
  for (col in setdiff(cols, colnames(table))) table[[col]] = NA_character_
  table = table[, cols, drop = FALSE]
  .pubchem_dedupe_profile(
    table = table,
    key_cols = c("Query", "CID", "Source", "TreeID", "TreeType", "HNID",
                 "ClassName", "ClassPath"),
    collapse_cols = c("TreeName", "RootHNID", "NodeID", "ParentNodeID",
                      "ParentClass", "SourceURL", "ClassificationURL",
                      "PubChemURL"),
    fallback_cols = c("ClassPath", "ClassName")
  )
}

.pubchem_extract_mesh_profile = function(pubchem) {
  cols = c("Query", "CID", "Descriptor", "PharmacologicAction",
           "TreeCategory", "RawValue", "Source", "SourceURL", "PubChemURL")
  annotations = .pubchem_filter_profile_annotations(
    pubchem,
    "MeSH|Medical Subject Headings|Pharmacologic Action|Descriptor|Tree Number"
  )
  if (nrow(annotations) < 1) return(.uaf_empty_table(cols))

  rows = lapply(seq_len(nrow(annotations)), function(i) {
    annotation = annotations[i, , drop = FALSE]
    text = .pubchem_annotation_text(annotation)
    data.frame(
      Query = annotation$Query,
      CID = annotation$CID,
      Descriptor = annotation$CleanValue,
      PharmacologicAction = if (grepl("pharmacologic action|therapeutic|drug",
                                      text, ignore.case = TRUE)) annotation$CleanValue else NA_character_,
      TreeCategory = .pubchem_collapse(.uaf_extract_pattern(text, "\\b[A-Z][0-9]{2}(?:\\.[0-9]+)+\\b")),
      RawValue = annotation$CleanValue,
      Source = annotation$Source,
      SourceURL = annotation$SourceURL,
      PubChemURL = annotation$PubChemURL,
      stringsAsFactors = FALSE
    )
  })
  .pubchem_dedupe_profile(
    table = do.call(rbind, rows),
    key_cols = c("Query", "CID", "Descriptor", "PharmacologicAction",
                 "TreeCategory"),
    collapse_cols = c("RawValue", "Source", "SourceURL", "PubChemURL"),
    fallback_cols = c("RawValue")
  )
}

.pubchem_extract_literature_profile = function(pubchem) {
  cols = c("Query", "CID", "PubMedID", "DOI", "PatentID", "ReferenceType",
           "TitleOrText", "RawValue", "Source", "SourceURL", "PubChemURL")
  annotations = .pubchem_filter_profile_annotations(
    pubchem,
    "Literature|PubMed|PMID|DOI|Patent|Reference|Citation"
  )
  if (nrow(annotations) < 1) return(.uaf_empty_table(cols))

  rows = lapply(seq_len(nrow(annotations)), function(i) {
    annotation = annotations[i, , drop = FALSE]
    text = .pubchem_annotation_text(annotation)
    identifiers = .uaf_extract_identifiers(text)
    pubmed = identifiers$Identifier[identifiers$IdentifierType == "PubMed"]
    doi = identifiers$Identifier[identifiers$IdentifierType == "DOI"]
    data.frame(
      Query = annotation$Query,
      CID = annotation$CID,
      PubMedID = .pubchem_collapse(pubmed),
      DOI = .pubchem_collapse(doi),
      PatentID = .pubchem_collapse(.uaf_extract_pattern(text, "\\b(?:US|EP|WO|JP|CN)[0-9A-Z-]{5,}\\b")),
      ReferenceType = .pubchem_reference_type(text),
      TitleOrText = annotation$CleanValue,
      RawValue = annotation$CleanValue,
      Source = annotation$Source,
      SourceURL = annotation$SourceURL,
      PubChemURL = annotation$PubChemURL,
      stringsAsFactors = FALSE
    )
  })
  .pubchem_dedupe_profile(
    table = do.call(rbind, rows),
    key_cols = c("Query", "CID", "PubMedID", "DOI", "PatentID",
                 "ReferenceType", "TitleOrText"),
    collapse_cols = c("RawValue", "Source", "SourceURL", "PubChemURL"),
    fallback_cols = c("RawValue", "TitleOrText")
  )
}

.pubchem_signal_word = function(text) {
  if (grepl("\\bDanger\\b", text, ignore.case = TRUE)) return("Danger")
  if (grepl("\\bWarning\\b", text, ignore.case = TRUE)) return("Warning")
  NA_character_
}

.pubchem_extract_fema_numbers = function(text) {
  matches = .uaf_extract_pattern(text, "\\bFEMA\\s*(?:No\\.?|Number|#)?\\s*[:#]?\\s*[0-9]{3,5}\\b")
  numbers = .uaf_extract_pattern(matches, "\\b[0-9]{3,5}\\b")
  unique(numbers)
}

.pubchem_extract_jecfa_numbers = function(text) {
  matches = .uaf_extract_pattern(text, "\\bJECFA\\s*(?:No\\.?|Number|#)?\\s*[:#]?\\s*[0-9]{1,5}\\b")
  unique(.uaf_extract_pattern(matches, "\\b[0-9]{1,5}\\b"))
}

.pubchem_extract_lotus_ids = function(text) {
  .uaf_extract_pattern(text, "\\b(?:LTS|LOTUS)[A-Z0-9:-]*[0-9][A-Z0-9:-]*\\b")
}

.pubchem_extract_taxonomy_ids = function(text) {
  text = paste(.uaf_non_empty(text), collapse = "; ")
  if (is.na(text) || text == "") return(character())
  matches = gregexpr("/taxonomy/([0-9]+)", text, perl = TRUE)
  values = regmatches(text, matches)[[1]]
  if (length(values) < 1 || identical(values, character(0))) {
    return(character())
  }
  unique(sub("^/taxonomy/([0-9]+)$", "\\1", values, perl = TRUE))
}

.pubchem_descriptor_terms = function(text, mode = c("all", "flavor", "odor")) {
  mode = match.arg(mode)
  terms = c("citrus", "lemon", "orange", "lime", "floral", "rose",
            "green", "herbal", "mint", "pine", "woody", "earthy",
            "fruity", "berry", "apple", "banana", "sweet", "bitter",
            "sour", "fatty", "waxy", "spicy", "pepper", "smoky",
            "nutty", "almond", "vanilla", "phenolic", "sulfurous",
            "alliaceous", "camphor", "terpene", "fresh")
  text_l = tolower(text)
  found = terms[vapply(terms, function(term) grepl(paste0("\\b", term, "\\b"),
                                                   text_l), logical(1))]
  if (mode == "flavor" && !grepl("flavo[u]?r|taste|aroma", text_l)) return(character())
  if (mode == "odor" && !grepl("odor|odour|aroma|scent|smell", text_l)) return(character())
  unique(found)
}

.pubchem_np_class_terms = function(text) {
  terms = c("alkaloid", "terpenoid", "terpene", "monoterpene",
            "sesquiterpene", "diterpene", "triterpene", "flavonoid",
            "phenylpropanoid", "polyketide", "peptide", "lipid",
            "steroid", "saponin", "coumarin", "lignan", "tannin")
  text_l = tolower(text)
  unique(terms[vapply(terms, function(term) grepl(paste0("\\b", term, "\\b"),
                                                  text_l), logical(1))])
}

.pubchem_extract_routes = function(text) {
  terms = c("oral", "intravenous", "topical", "inhalation",
            "subcutaneous", "intramuscular", "ophthalmic", "nasal",
            "transdermal", "rectal", "sublingual", "buccal",
            "intradermal", "epidural")
  text_l = tolower(text)
  unique(terms[vapply(terms, function(term) grepl(paste0("\\b", term, "\\b"),
                                                  text_l), logical(1))])
}

.pubchem_extract_dosage_forms = function(text) {
  terms = c("tablet", "capsule", "solution", "suspension", "injection",
            "cream", "ointment", "gel", "patch", "spray", "powder",
            "suppository", "syrup")
  text_l = tolower(text)
  unique(terms[vapply(terms, function(term) grepl(paste0("\\b", term, "\\b"),
                                                  text_l), logical(1))])
}

.pubchem_warning_terms = function(text) {
  terms = c("contraindication", "adverse reaction", "boxed warning",
            "warning", "precaution", "interaction", "overdosage",
            "pregnancy", "carcinogenesis", "hypersensitivity")
  text_l = tolower(text)
  unique(terms[vapply(terms, function(term) grepl(term, text_l, fixed = TRUE),
                      logical(1))])
}

.pubchem_gras_status = function(text) {
  if (grepl("\\bGRAS\\b|Generally Recognized as Safe", text,
            ignore.case = TRUE)) {
    return("GRAS")
  }
  NA_character_
}

.pubchem_reference_type = function(text) {
  if (grepl("patent", text, ignore.case = TRUE)) return("Patent")
  if (grepl("PubMed|PMID|DOI|journal|article", text, ignore.case = TRUE)) {
    return("Literature")
  }
  if (grepl("reference|citation", text, ignore.case = TRUE)) return("Reference")
  NA_character_
}

.pubchem_last_path = function(path) {
  path = .uaf_squish_text(path)
  if (is.na(path) || path == "") return(NA_character_)
  pieces = .uaf_non_empty(strsplit(path, ">", fixed = TRUE)[[1]])
  if (length(pieces) < 1) return(NA_character_)
  pieces[[length(pieces)]]
}

.pubchem_collapse = function(x) {
  x = unique(.uaf_non_empty(x))
  if (length(x) < 1) return(NA_character_)
  paste(x, collapse = "; ")
}
