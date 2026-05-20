.categorate_normalized_outputs = function(compounds, data_list, pubchem, kegg,
                                          pubchem_profiles,
                                          trait_matrix_profile = "core",
                                          trait_matrix_mode = "binary",
                                          trait_matrix_min_confidence = 0,
                                          trait_matrix_max_traits = Inf) {
  cid_lookup = .normalized_cid_lookup(pubchem)
  chemical_classes = .normalized_chemical_classes(compounds = compounds,
                                                  data_list = data_list,
                                                  pubchem_profiles = pubchem_profiles,
                                                  kegg = kegg,
                                                  cid_lookup = cid_lookup)
  chemical_measurements = .normalized_chemical_measurements(pubchem = pubchem,
                                                            pubchem_profiles = pubchem_profiles)
  chemical_measurement_summary = chemicalMeasurementSummary(chemical_measurements)
  chemical_hazards = .normalized_chemical_hazards(pubchem_profiles)
  chemical_uses = .normalized_chemical_uses(pubchem_profiles)
  chemical_bioassays = .normalized_chemical_bioassays(pubchem)
  chemical_bioactivities = .normalized_chemical_bioactivities(chemical_bioassays)
  chemical_targets = .normalized_chemical_targets(chemical_bioassays)
  chemical_potencies = .normalized_chemical_potencies(chemical_bioassays)
  chemical_taxonomy = .normalized_chemical_taxonomy(pubchem_profiles)
  chemical_occurrences = .normalized_chemical_occurrences(chemical_taxonomy)
  chemical_pathway_roles = .normalized_chemical_pathway_roles(kegg = kegg,
                                                              cid_lookup = cid_lookup)
  kegg_reaction_participants = .normalized_kegg_reaction_participants(kegg)
  chemical_terms = .normalized_chemical_terms(chemical_classes = chemical_classes,
                                              chemical_hazards = chemical_hazards,
                                              chemical_uses = chemical_uses,
                                              chemical_bioactivities = chemical_bioactivities,
                                              chemical_targets = chemical_targets,
                                              chemical_taxonomy = chemical_taxonomy,
                                              chemical_pathway_roles = chemical_pathway_roles)
  chemical_traits = .normalized_chemical_traits(
    chemical_terms = chemical_terms,
    chemical_classes = chemical_classes,
    chemical_measurements = chemical_measurements,
    chemical_hazards = chemical_hazards,
    chemical_uses = chemical_uses,
    chemical_bioactivities = chemical_bioactivities,
    chemical_targets = chemical_targets,
    chemical_potencies = chemical_potencies,
    chemical_taxonomy = chemical_taxonomy,
    chemical_occurrences = chemical_occurrences,
    chemical_pathway_roles = chemical_pathway_roles,
    kegg_reaction_participants = kegg_reaction_participants
  )
  chemical_trait_ontology = chemicalTraitOntology(
    chemical_traits,
    min_confidence = trait_matrix_min_confidence
  )
  chemical_trait_matrix = chemicalTraitMatrix(
    chemical_traits,
    profile = trait_matrix_profile,
    mode = trait_matrix_mode,
    min_confidence = trait_matrix_min_confidence,
    max_traits = trait_matrix_max_traits
  )
  chemical_trait_ontology_matrix = chemicalTraitOntologyMatrix(
    chemical_trait_ontology,
    mode = trait_matrix_mode,
    min_confidence = trait_matrix_min_confidence,
    max_terms = trait_matrix_max_traits
  )
  chemical_trait_evidence = chemicalTraitEvidence(
    chemical_trait_ontology,
    type = "ontology",
    min_confidence = trait_matrix_min_confidence
  )
  chemical_trait_report = chemicalTraitReport(
    list(ChemicalTraits = chemical_traits,
         ChemicalTraitOntology = chemical_trait_ontology),
    min_confidence = trait_matrix_min_confidence
  )
  chemical_trait_summary = chemicalTraitSummary(
    chemical_traits,
    by = "type",
    min_confidence = trait_matrix_min_confidence
  )
  chemical_trait_similarity = chemicalTraitSimilarity(
    chemical_traits,
    profile = trait_matrix_profile,
    min_confidence = trait_matrix_min_confidence,
    max_traits = trait_matrix_max_traits
  )

  list(
    ChemicalTerms = chemical_terms,
    ChemicalTraits = chemical_traits,
    ChemicalTraitOntology = chemical_trait_ontology,
    ChemicalTraitMatrix = chemical_trait_matrix,
    ChemicalTraitOntologyMatrix = chemical_trait_ontology_matrix,
    ChemicalTraitEvidence = chemical_trait_evidence,
    ChemicalTraitReport = chemical_trait_report,
    ChemicalTraitSummary = chemical_trait_summary,
    ChemicalTraitSimilarity = chemical_trait_similarity,
    ChemicalClasses = chemical_classes,
    ChemicalMeasurements = chemical_measurements,
    ChemicalMeasurementSummary = chemical_measurement_summary,
    ChemicalHazards = chemical_hazards,
    ChemicalUses = chemical_uses,
    ChemicalBioassays = chemical_bioassays,
    ChemicalBioactivities = chemical_bioactivities,
    ChemicalTargets = chemical_targets,
    ChemicalPotencies = chemical_potencies,
    ChemicalTaxonomy = chemical_taxonomy,
    ChemicalOccurrences = chemical_occurrences,
    ChemicalPathwayRoles = chemical_pathway_roles,
    KEGGReactionParticipants = kegg_reaction_participants
  )
}

.normalized_cid_lookup = function(pubchem) {
  if (is.null(pubchem) ||
      !is.list(pubchem) ||
      !is.data.frame(pubchem$identity) ||
      nrow(pubchem$identity) < 1 ||
      !"Query" %in% colnames(pubchem$identity) ||
      !"CID" %in% colnames(pubchem$identity)) {
    return(stats::setNames(character(), character()))
  }
  out = pubchem$identity$CID
  names(out) = pubchem$identity$Query
  out
}

.normalized_lookup_cid = function(query, cid_lookup) {
  if (length(cid_lookup) < 1 || is.na(query) || query == "") return(NA_integer_)
  value = unname(cid_lookup[[query]])
  if (is.null(value) || length(value) < 1 || is.na(value)) return(NA_integer_)
  suppressWarnings(as.integer(value))
}

.normalized_bind_rows = function(rows, cols) {
  rows = rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) < 1) return(.uaf_empty_table(cols))
  aligned = lapply(rows, function(row) {
    for (col in setdiff(cols, colnames(row))) row[[col]] = NA
    row[, cols, drop = FALSE]
  })
  out = unique(do.call(rbind, aligned))
  row.names(out) = NULL
  out
}

.normalized_dedupe_table = function(table, key_cols,
                                    collapse_cols = c("SourceTable", "Source",
                                                      "EvidenceText",
                                                      "EvidenceURL",
                                                      "PubChemURL",
                                                      "ExtractionRule")) {
  if (!is.data.frame(table) || nrow(table) < 2) return(table)
  key_cols = intersect(key_cols, colnames(table))
  if (length(key_cols) < 1) return(unique(table))
  key = .pubchem_row_key(table, key_cols, fallback_cols = key_cols)
  rows = lapply(unique(key), function(group_key) {
    group = table[key == group_key, , drop = FALSE]
    row = group[1, , drop = FALSE]
    for (col in intersect(collapse_cols, colnames(group))) {
      row[[col]] = .pubchem_collapse(group[[col]])
    }
    if ("Confidence" %in% colnames(group)) {
      row$Confidence = .normalized_best_confidence(group$Confidence)
    }
    row
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.normalized_best_confidence = function(confidence) {
  confidence = .uaf_non_empty(confidence)
  if ("high" %in% confidence) return("high")
  if ("medium" %in% confidence) return("medium")
  if ("low" %in% confidence) return("low")
  NA_character_
}

.normalized_term_dictionary = function() {
  data.frame(
    Pattern = c(
      "citrus|lemon|orange|lime",
      "sweet|honey|sugary",
      "floral|rose|jasmine",
      "green|fresh|leafy|grassy",
      "herbal|mint|menthol",
      "pine|woody|wood",
      "earthy|soil|mushroom",
      "fruity|berry|apple|banana",
      "fatty|waxy|oily",
      "spicy|pepper|clove",
      "smoky|smoke|burnt",
      "nutty|almond",
      "vanilla|vanillin",
      "phenolic|phenol",
      "sulfur|sulphur|sulfurous|alliaceous|garlic|onion",
      "camphor",
      "alkaloid",
      "monoterpene",
      "sesquiterpene",
      "diterpene",
      "triterpene",
      "terpenoid|terpene",
      "flavonoid",
      "phenylpropanoid",
      "polyketide",
      "peptide",
      "lipid|fatty acid",
      "steroid",
      "saponin",
      "coumarin",
      "lignan",
      "tannin",
      "anti[- ]?inflammatory|non[- ]?steroidal",
      "cyclooxygenase|cox",
      "analgesic|pain",
      "antipyretic|fever",
      "antimicrobial|antibacterial|antifungal",
      "antioxidant",
      "inhibitor",
      "agonist",
      "antagonist",
      "oral|by mouth|per os",
      "intravenous|\\biv\\b|i\\.v\\.",
      "topical|dermal|cutaneous",
      "inhalation|inhaled",
      "subcutaneous",
      "intramuscular",
      "tablet",
      "capsule",
      "solution",
      "suspension",
      "injection",
      "cream|ointment|gel",
      "patch|transdermal",
      "flammable|combustible",
      "irritat|eye damage|eye irritation",
      "corrosive|corrosion",
      "toxic|toxicity|fatal|harmful",
      "carcinogen|carcinogenic",
      "aquatic|environmental"
    ),
    TermClean = c(
      "citrus", "sweet", "floral", "green/fresh", "herbal/mint",
      "woody/pine", "earthy", "fruity", "fatty/waxy", "spicy",
      "smoky", "nutty", "vanilla", "phenolic", "sulfurous",
      "camphor", "alkaloid", "monoterpene", "sesquiterpene",
      "diterpene", "triterpene", "terpenoid", "flavonoid",
      "phenylpropanoid", "polyketide", "peptide", "lipid",
      "steroid", "saponin", "coumarin", "lignan", "tannin",
      "anti-inflammatory", "cyclooxygenase inhibitor", "analgesic",
      "antipyretic", "antimicrobial", "antioxidant", "inhibitor",
      "agonist", "antagonist", "oral", "intravenous", "topical",
      "inhalation", "subcutaneous", "intramuscular", "tablet",
      "capsule", "solution", "suspension", "injection",
      "semisolid topical", "transdermal patch", "flammable",
      "irritant", "corrosive", "toxic", "carcinogen",
      "aquatic toxicity"
    ),
    TermGroup = c(
      "citrus", "sweet", "floral", "green/fresh", "herbal/mint",
      "woody/pine", "earthy", "fruity", "fatty/waxy", "spicy",
      "smoky", "nutty", "vanilla", "phenolic", "sulfurous",
      "camphoraceous", "alkaloid", "terpenoid", "terpenoid",
      "terpenoid", "terpenoid", "terpenoid", "phenolic",
      "phenylpropanoid", "polyketide", "peptide", "lipid",
      "steroid", "saponin", "coumarin", "lignan", "tannin",
      "anti-inflammatory", "enzyme inhibitor", "analgesic",
      "antipyretic", "antimicrobial", "antioxidant", "bioactivity",
      "bioactivity", "bioactivity", "enteral route",
      "parenteral route", "topical route", "respiratory route",
      "parenteral route", "parenteral route", "solid dosage form",
      "solid dosage form", "liquid dosage form", "liquid dosage form",
      "parenteral dosage form", "topical dosage form",
      "transdermal dosage form", "physical hazard", "health hazard",
      "health hazard", "health hazard", "health hazard",
      "environmental hazard"
    ),
    stringsAsFactors = FALSE
  )
}

.normalized_term_lookup = function(term) {
  term = .uaf_squish_text(term)
  if (is.na(term) || term == "") {
    return(list(TermClean = NA_character_,
                TermGroup = NA_character_,
                Confidence = "low"))
  }
  if (grepl("^H[0-9]{3}[A-Za-z]?$", term, ignore.case = TRUE)) {
    code = toupper(term)
    group = .normalized_hazard_group(code, code)
    return(list(TermClean = tolower(code),
                TermGroup = .uaf_first_non_empty_text(group, "hazard_code"),
                Confidence = "high"))
  }
  if (grepl("anti[- ]?inflammatory|non[- ]?steroidal", term,
            ignore.case = TRUE, perl = TRUE)) {
    return(list(TermClean = "anti-inflammatory",
                TermGroup = "anti-inflammatory",
                Confidence = "high"))
  }
  if (grepl("cyclooxygenase|\\bcox\\b", term,
            ignore.case = TRUE, perl = TRUE)) {
    return(list(TermClean = "cyclooxygenase inhibitor",
                TermGroup = "enzyme inhibitor",
                Confidence = "high"))
  }
  dict = .normalized_term_dictionary()
  text = tolower(term)
  hit = which(vapply(dict$Pattern, function(pattern) {
    grepl(pattern, text, perl = TRUE, ignore.case = TRUE)
  }, logical(1)))
  if (length(hit) > 0) {
    row = dict[hit[[1]], , drop = FALSE]
    return(list(TermClean = row$TermClean,
                TermGroup = row$TermGroup,
                Confidence = "high"))
  }
  clean = .normalized_clean_label(term)
  list(TermClean = clean, TermGroup = clean, Confidence = "medium")
}

.normalized_clean_label = function(x) {
  x = .uaf_squish_text(x)
  x = tolower(x)
  x = gsub("[_]+", " ", x)
  x = gsub("[^a-z0-9+./ -]+", " ", x)
  x = gsub("\\s+", " ", x)
  trimws(x)
}

.normalized_split_terms = function(x) {
  x = .uaf_non_empty(x)
  if (length(x) < 1) return(character())
  terms = unlist(strsplit(paste(x, collapse = "; "), "\\s*(?:;|\\|)\\s*",
                          perl = TRUE), use.names = FALSE)
  terms = .uaf_non_empty(terms)
  terms = terms[!grepl("^[A-Z][0-9]{2}(?:\\.[0-9]+)+$", terms)]
  terms = terms[!grepl("^(?:NA|None|Unknown)$", terms, ignore.case = TRUE)]
  unique(terms)
}

.normalized_split_codes = function(x) {
  x = .uaf_non_empty(x)
  if (length(x) < 1) return(character())
  values = unlist(strsplit(paste(x, collapse = "; "), "\\s*(?:;|\\|)\\s*",
                           perl = TRUE), use.names = FALSE)
  unique(.uaf_non_empty(values))
}

.normalized_classification_path_terms = function(path, fallback = NA_character_) {
  values = .uaf_non_empty(c(path, fallback))
  if (length(values) < 1) return(character())
  terms = unlist(strsplit(values[[1]], "\\s*(?:;|\\|>|>|\\|)\\s*",
                          perl = TRUE), use.names = FALSE)
  terms = .uaf_non_empty(terms)
  terms = terms[!grepl("^(?:biological|chemical) tree$", terms,
                       ignore.case = TRUE)]
  unique(terms)
}

.normalized_source_url = function(row) {
  .uaf_first_non_empty_text(row$SourceURL, row$EvidenceURL, row$PubChemURL)
}

.normalized_add_term = function(rows, query, cid, domain, term_type, term_raw,
                                source_table, source, evidence_text,
                                evidence_url, pubchem_url = NA_character_,
                                extraction_rule = "dictionary") {
  lookup = .normalized_term_lookup(term_raw)
  rows[[length(rows) + 1]] = data.frame(
    Query = query,
    CID = suppressWarnings(as.integer(cid)),
    Domain = domain,
    TermType = term_type,
    TermRaw = .uaf_squish_text(term_raw),
    TermClean = lookup$TermClean,
    TermGroup = lookup$TermGroup,
    SourceTable = source_table,
    Source = source,
    EvidenceText = evidence_text,
    EvidenceURL = evidence_url,
    PubChemURL = pubchem_url,
    ExtractionRule = extraction_rule,
    Confidence = lookup$Confidence,
    stringsAsFactors = FALSE
  )
  rows
}

.normalized_chemical_terms = function(chemical_classes, chemical_hazards,
                                      chemical_uses, chemical_bioactivities,
                                      chemical_targets, chemical_taxonomy,
                                      chemical_pathway_roles) {
  cols = c("Query", "CID", "Domain", "TermType", "TermRaw", "TermClean",
           "TermGroup", "SourceTable", "Source", "EvidenceText",
           "EvidenceURL", "PubChemURL", "ExtractionRule", "Confidence")
  rows = list()

  if (is.data.frame(chemical_classes) && nrow(chemical_classes) > 0) {
    for (i in seq_len(nrow(chemical_classes))) {
      row = chemical_classes[i, , drop = FALSE]
      raw = .uaf_first_non_empty_text(row$ClassName, row$ClassGroup)
      if (is.na(raw)) next
      rows = .normalized_add_term(rows, row$Query, row$CID,
                                  domain = "classification",
                                  term_type = row$ClassType,
                                  term_raw = raw,
                                  source_table = "ChemicalClasses",
                                  source = row$ClassSystem,
                                  evidence_text = row$EvidenceText,
                                  evidence_url = row$EvidenceURL,
                                  extraction_rule = row$ExtractionRule)
    }
  }

  if (is.data.frame(chemical_hazards) && nrow(chemical_hazards) > 0) {
    for (i in seq_len(nrow(chemical_hazards))) {
      row = chemical_hazards[i, , drop = FALSE]
      terms = .normalized_split_terms(c(row$HazardCode, row$HazardCategory,
                                        row$HazardGroup, row$TargetOrgan,
                                        row$SignalWord))
      for (term in terms) {
        rows = .normalized_add_term(rows, row$Query, row$CID,
                                    domain = "safety",
                                    term_type = "hazard",
                                    term_raw = term,
                                    source_table = "ChemicalHazards",
                                    source = row$Source,
                                    evidence_text = row$EvidenceText,
                                    evidence_url = row$EvidenceURL,
                                    pubchem_url = row$PubChemURL,
                                    extraction_rule = row$ExtractionRule)
      }
    }
  }

  if (is.data.frame(chemical_uses) && nrow(chemical_uses) > 0) {
    for (i in seq_len(nrow(chemical_uses))) {
      row = chemical_uses[i, , drop = FALSE]
      rows = .normalized_add_term(rows, row$Query, row$CID,
                                  domain = row$UseDomain,
                                  term_type = row$UseType,
                                  term_raw = row$UseTerm,
                                  source_table = "ChemicalUses",
                                  source = row$Source,
                                  evidence_text = row$EvidenceText,
                                  evidence_url = row$EvidenceURL,
                                  pubchem_url = row$PubChemURL,
                                  extraction_rule = row$ExtractionRule)
    }
  }

  if (is.data.frame(chemical_bioactivities) &&
      nrow(chemical_bioactivities) > 0) {
    for (i in seq_len(nrow(chemical_bioactivities))) {
      row = chemical_bioactivities[i, , drop = FALSE]
      terms = .normalized_split_terms(c(row$ActivityClass,
                                        row$ActivityDirection,
                                        row$BioactivityDomain,
                                        row$AssayType))
      for (term in terms) {
        rows = .normalized_add_term(rows, row$Query, row$CID,
                                    domain = "bioactivity",
                                    term_type = "bioactivity",
                                    term_raw = term,
                                    source_table = "ChemicalBioactivities",
                                    source = "PubChem BioAssay",
                                    evidence_text = row$EvidenceText,
                                    evidence_url = row$EvidenceURL,
                                    extraction_rule = row$ExtractionRule)
      }
    }
  }

  if (is.data.frame(chemical_targets) && nrow(chemical_targets) > 0) {
    for (i in seq_len(nrow(chemical_targets))) {
      row = chemical_targets[i, , drop = FALSE]
      terms = .normalized_split_terms(c(row$TargetName, row$TargetType,
                                        row$TargetGroup,
                                        row$TargetOrganism,
                                        row$TargetCommonName,
                                        row$BioactivityDomains,
                                        row$EndpointNames))
      for (term in terms) {
        rows = .normalized_add_term(rows, row$Query, row$CID,
                                    domain = "bioactivity",
                                    term_type = "target",
                                    term_raw = term,
                                    source_table = "ChemicalTargets",
                                    source = "PubChem BioAssay",
                                    evidence_text = row$EvidenceText,
                                    evidence_url = row$EvidenceURL,
                                    extraction_rule = row$ExtractionRule)
      }
    }
  }

  if (is.data.frame(chemical_taxonomy) && nrow(chemical_taxonomy) > 0) {
    for (i in seq_len(nrow(chemical_taxonomy))) {
      row = chemical_taxonomy[i, , drop = FALSE]
      terms = .normalized_split_terms(c(row$Domain, row$Kingdom, row$Phylum,
                                        row$Class,
                                        row$Order, row$Family, row$Genus,
                                        row$Species, row$TaxonomyGroup,
                                        row$TaxonomyTerms,
                                        row$NaturalProductClass))
      for (term in terms) {
        rows = .normalized_add_term(rows, row$Query, row$CID,
                                    domain = "ecology",
                                    term_type = "taxonomy",
                                    term_raw = term,
                                    source_table = "ChemicalTaxonomy",
                                    source = row$Source,
                                    evidence_text = row$EvidenceText,
                                    evidence_url = row$EvidenceURL,
                                    pubchem_url = row$PubChemURL,
                                    extraction_rule = row$ExtractionRule)
      }
    }
  }

  if (is.data.frame(chemical_pathway_roles) && nrow(chemical_pathway_roles) > 0) {
    for (i in seq_len(nrow(chemical_pathway_roles))) {
      row = chemical_pathway_roles[i, , drop = FALSE]
      terms = .normalized_split_terms(c(row$PathwayGroup, row$EnzymeClass,
                                        row$RoleType))
      for (term in terms) {
        rows = .normalized_add_term(rows, row$Query, row$CID,
                                    domain = "metabolism",
                                    term_type = "pathway_role",
                                    term_raw = term,
                                    source_table = "ChemicalPathwayRoles",
                                    source = "KEGG",
                                    evidence_text = .pubchem_collapse(c(row$PathwayName,
                                                                         row$ReactionName,
                                                                         row$EnzymeName)),
                                    evidence_url = row$EvidenceURL,
                                    extraction_rule = row$ExtractionRule)
      }
    }
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "CID", "Domain", "TermType", "TermClean",
                 "TermGroup"),
    collapse_cols = c("TermRaw", "SourceTable", "Source", "EvidenceText",
                      "EvidenceURL", "PubChemURL", "ExtractionRule")
  )
}

.normalized_chemical_traits = function(chemical_terms, chemical_classes,
                                       chemical_measurements,
                                       chemical_hazards, chemical_uses,
                                       chemical_bioactivities,
                                       chemical_targets,
                                       chemical_potencies,
                                       chemical_taxonomy,
                                       chemical_occurrences,
                                       chemical_pathway_roles,
                                       kegg_reaction_participants) {
  cols = .normalized_trait_cols()
  rows = list()
  cid_lookup = .normalized_trait_cid_lookup(list(
    chemical_terms,
    chemical_classes,
    chemical_measurements,
    chemical_hazards,
    chemical_uses,
    chemical_bioactivities,
    chemical_targets,
    chemical_potencies,
    chemical_taxonomy,
    chemical_occurrences,
    chemical_pathway_roles
  ))

  if (is.data.frame(chemical_classes) && nrow(chemical_classes) > 0) {
    for (i in seq_len(nrow(chemical_classes))) {
      row = chemical_classes[i, , drop = FALSE]
      source_database = .normalized_trait_source_database(row$SourceTable,
                                                          row$ClassSystem)
      rows = .normalized_add_trait(
        rows, row$Query, row$CID,
        trait_type = "chemical_class",
        trait_group = row$ClassType,
        trait_value = .uaf_first_non_empty_text(row$ClassName, row$ClassGroup),
        source_database = source_database,
        source_table = row$SourceTable,
        evidence_id = row$ClassID,
        evidence_text = row$EvidenceText,
        evidence_url = row$EvidenceURL,
        extraction_rule = row$ExtractionRule,
        confidence = row$Confidence
      )
    }
  }

  if (is.data.frame(chemical_measurements) &&
      nrow(chemical_measurements) > 0) {
    for (i in seq_len(nrow(chemical_measurements))) {
      row = chemical_measurements[i, , drop = FALSE]
      measurement_traits = .normalized_measurement_traits(row$Property,
                                                          row$Value)
      if (nrow(measurement_traits) < 1) next
      for (j in seq_len(nrow(measurement_traits))) {
        rows = .normalized_add_trait(
          rows, row$Query, row$CID,
          trait_type = "property",
          trait_group = measurement_traits$TraitGroup[[j]],
          trait_value = measurement_traits$TraitValue[[j]],
          source_database = "PubChem",
          source_table = row$SourceTable,
          evidence_id = row$Property,
          evidence_text = row$EvidenceText,
          evidence_url = row$EvidenceURL,
          extraction_rule = .pubchem_collapse(c(row$ExtractionRule,
                                                "property_bin")),
          confidence = measurement_traits$Confidence[[j]]
        )
      }
    }
  }

  if (is.data.frame(chemical_hazards) && nrow(chemical_hazards) > 0) {
    for (i in seq_len(nrow(chemical_hazards))) {
      row = chemical_hazards[i, , drop = FALSE]
      hazard_values = list(
        hazard_code = row$HazardCode,
        hazard_category = row$HazardCategory,
        hazard_group = row$HazardGroup,
        signal_word = row$SignalWord,
        exposure_route = row$ExposureRoute,
        target_organ = row$TargetOrgan,
        toxicity_metric = row$ToxicityMetric,
        toxicity_species = row$Species
      )
      for (trait_group in names(hazard_values)) {
        rows = .normalized_add_trait(
          rows, row$Query, row$CID,
          trait_type = "hazard",
          trait_group = trait_group,
          trait_value = hazard_values[[trait_group]],
          source_database = "PubChem",
          source_table = "ChemicalHazards",
          evidence_id = row$HazardCode,
          evidence_text = row$EvidenceText,
          evidence_url = row$EvidenceURL,
          extraction_rule = row$ExtractionRule,
          confidence = row$Confidence,
          split = TRUE
        )
      }
    }
  }

  if (is.data.frame(chemical_uses) && nrow(chemical_uses) > 0) {
    for (i in seq_len(nrow(chemical_uses))) {
      row = chemical_uses[i, , drop = FALSE]
      source_database = .normalized_trait_source_database(row$SourceTable,
                                                          row$Source)
      rows = .normalized_add_trait(
        rows, row$Query, row$CID,
        trait_type = "use",
        trait_group = row$UseType,
        trait_value = row$UseTerm,
        source_database = source_database,
        source_table = "ChemicalUses",
        evidence_text = row$EvidenceText,
        evidence_url = row$EvidenceURL,
        extraction_rule = row$ExtractionRule,
        confidence = row$Confidence
      )
      rows = .normalized_add_trait(
        rows, row$Query, row$CID,
        trait_type = "use",
        trait_group = "use_group",
        trait_value = row$UseGroup,
        source_database = source_database,
        source_table = "ChemicalUses",
        evidence_text = row$EvidenceText,
        evidence_url = row$EvidenceURL,
        extraction_rule = row$ExtractionRule,
        confidence = row$Confidence
      )
    }
  }

  if (is.data.frame(chemical_bioactivities) &&
      nrow(chemical_bioactivities) > 0) {
    for (i in seq_len(nrow(chemical_bioactivities))) {
      row = chemical_bioactivities[i, , drop = FALSE]
      bioactivity_values = list(
        activity_class = row$ActivityClass,
        activity_direction = row$ActivityDirection,
        bioactivity_domain = row$BioactivityDomain,
        assay_type = row$AssayType,
        endpoint_metric = row$EndpointNames
      )
      for (trait_group in names(bioactivity_values)) {
        rows = .normalized_add_trait(
          rows, row$Query, row$CID,
          trait_type = "bioactivity",
          trait_group = trait_group,
          trait_value = bioactivity_values[[trait_group]],
          source_database = "PubChem BioAssay",
          source_table = "ChemicalBioactivities",
          evidence_id = row$AID,
          evidence_text = row$EvidenceText,
          evidence_url = row$EvidenceURL,
          extraction_rule = row$ExtractionRule,
          confidence = row$Confidence,
          split = TRUE
        )
      }
    }
  }

  if (is.data.frame(chemical_targets) && nrow(chemical_targets) > 0) {
    for (i in seq_len(nrow(chemical_targets))) {
      row = chemical_targets[i, , drop = FALSE]
      target_values = list(
        target_name = row$TargetName,
        target_type = row$TargetType,
        target_accession = row$TargetAccession,
        target_gene = row$TargetGeneID,
        target_taxonomy = row$TargetTaxonomyID,
        target_organism = row$TargetOrganism,
        target_common_name = row$TargetCommonName,
        target_group = row$TargetGroup,
        endpoint_metric = row$EndpointNames
      )
      for (trait_group in names(target_values)) {
        rows = .normalized_add_trait(
          rows, row$Query, row$CID,
          trait_type = "target",
          trait_group = trait_group,
          trait_value = target_values[[trait_group]],
          source_database = "PubChem BioAssay",
          source_table = "ChemicalTargets",
          evidence_id = .uaf_first_non_empty_text(row$TargetGeneID,
                                                  row$TargetAccession,
                                                  row$TargetTaxonomyID),
          evidence_text = row$EvidenceText,
          evidence_url = row$EvidenceURL,
          extraction_rule = row$ExtractionRule,
          confidence = row$Confidence,
          split = TRUE
        )
      }
    }
  }

  if (is.data.frame(chemical_potencies) && nrow(chemical_potencies) > 0) {
    for (i in seq_len(nrow(chemical_potencies))) {
      row = chemical_potencies[i, , drop = FALSE]
      potency_values = list(
        potency_metric = row$PotencyMetric,
        potency_bucket = row$PotencyBucket,
        target_name = row$TargetName
      )
      for (trait_group in names(potency_values)) {
        rows = .normalized_add_trait(
          rows, row$Query, row$CID,
          trait_type = "potency",
          trait_group = trait_group,
          trait_value = potency_values[[trait_group]],
          source_database = "PubChem BioAssay",
          source_table = "ChemicalPotencies",
          evidence_id = row$AID,
          evidence_text = row$EvidenceText,
          evidence_url = row$EvidenceURL,
          extraction_rule = row$ExtractionRule,
          confidence = row$Confidence
        )
      }
    }
  }

  if (is.data.frame(chemical_taxonomy) && nrow(chemical_taxonomy) > 0) {
    for (i in seq_len(nrow(chemical_taxonomy))) {
      row = chemical_taxonomy[i, , drop = FALSE]
      taxonomy_values = list(
        domain = row$Domain,
        kingdom = row$Kingdom,
        phylum = row$Phylum,
        class = row$Class,
        order = row$Order,
        family = row$Family,
        genus = row$Genus,
        species = row$Species,
        organism = row$Organism,
        common_name = row$CommonName,
        taxonomy_group = row$TaxonomyGroup,
        natural_product_class = row$NaturalProductClass
      )
      for (trait_group in names(taxonomy_values)) {
        rows = .normalized_add_trait(
          rows, row$Query, row$CID,
          trait_type = "taxonomy",
          trait_group = trait_group,
          trait_value = taxonomy_values[[trait_group]],
          source_database = .normalized_trait_source_database("ChemicalTaxonomy",
                                                              row$Source),
          source_table = "ChemicalTaxonomy",
          evidence_id = row$TaxonomyID,
          evidence_text = row$EvidenceText,
          evidence_url = row$EvidenceURL,
          extraction_rule = row$ExtractionRule,
          confidence = row$Confidence,
          split = trait_group == "natural_product_class"
        )
      }
    }
  }

  if (is.data.frame(chemical_occurrences) &&
      nrow(chemical_occurrences) > 0) {
    for (i in seq_len(nrow(chemical_occurrences))) {
      row = chemical_occurrences[i, , drop = FALSE]
      occurrence_values = list(
        occurrence_type = row$OccurrenceType,
        kingdom = row$Kingdom,
        family = row$Family,
        genus = row$Genus,
        species = row$Species,
        organism = row$Organism,
        natural_product_class = row$NaturalProductClass
      )
      for (trait_group in names(occurrence_values)) {
        rows = .normalized_add_trait(
          rows, row$Query, row$CID,
          trait_type = "occurrence",
          trait_group = trait_group,
          trait_value = occurrence_values[[trait_group]],
          source_database = row$SourceDatabase,
          source_table = "ChemicalOccurrences",
          evidence_id = row$TaxonomyID,
          evidence_text = row$EvidenceText,
          evidence_url = row$EvidenceURL,
          extraction_rule = row$ExtractionRule,
          confidence = row$Confidence,
          split = trait_group == "natural_product_class"
        )
      }
    }
  }

  if (is.data.frame(chemical_pathway_roles) &&
      nrow(chemical_pathway_roles) > 0) {
    for (i in seq_len(nrow(chemical_pathway_roles))) {
      row = chemical_pathway_roles[i, , drop = FALSE]
      pathway_values = list(
        pathway_group = row$PathwayGroup,
        pathway_name = row$PathwayName,
        reaction_name = row$ReactionName,
        enzyme_class = row$EnzymeClass,
        ec_number = row$ECNumber,
        role_type = row$RoleType
      )
      for (trait_group in names(pathway_values)) {
        rows = .normalized_add_trait(
          rows, row$Query, row$CID,
          trait_type = "pathway",
          trait_group = trait_group,
          trait_value = pathway_values[[trait_group]],
          source_database = "KEGG",
          source_table = "ChemicalPathwayRoles",
          evidence_id = .uaf_first_non_empty_text(row$PathwayID,
                                                  row$ReactionID,
                                                  row$ECNumber,
                                                  row$KEGG_ID),
          evidence_text = .pubchem_collapse(c(row$PathwayName,
                                               row$ReactionName,
                                               row$EnzymeName)),
          evidence_url = row$EvidenceURL,
          extraction_rule = row$ExtractionRule,
          confidence = row$Confidence
        )
      }
    }
  }

  if (is.data.frame(kegg_reaction_participants) &&
      nrow(kegg_reaction_participants) > 0) {
    for (i in seq_len(nrow(kegg_reaction_participants))) {
      row = kegg_reaction_participants[i, , drop = FALSE]
      reaction_values = list(
        participant_role = row$ParticipantRole,
        participant_side = row$Side,
        participant_id = row$ParticipantID,
        participant_name = row$ParticipantName
      )
      for (trait_group in names(reaction_values)) {
        rows = .normalized_add_trait(
          rows, row$Query, .normalized_lookup_cid(row$Query, cid_lookup),
          trait_type = "reaction",
          trait_group = trait_group,
          trait_value = reaction_values[[trait_group]],
          source_database = "KEGG",
          source_table = "KEGGReactionParticipants",
          evidence_id = .uaf_first_non_empty_text(row$ReactionID,
                                                  row$ParticipantID,
                                                  row$KEGG_ID),
          evidence_text = .pubchem_collapse(c(row$Definition, row$Equation)),
          evidence_url = row$EvidenceURL,
          extraction_rule = row$ExtractionRule,
          confidence = row$Confidence,
          matrix_eligible = trait_group != "participant_name"
        )
      }
    }
  }

  if (is.data.frame(chemical_terms) && nrow(chemical_terms) > 0) {
    for (i in seq_len(nrow(chemical_terms))) {
      row = chemical_terms[i, , drop = FALSE]
      rows = .normalized_add_trait(
        rows, row$Query, row$CID,
        trait_type = "vocabulary",
        trait_group = row$TermType,
        trait_value = row$TermClean,
        source_database = .normalized_trait_source_database(row$SourceTable,
                                                            row$Source),
        source_table = "ChemicalTerms",
        evidence_text = row$EvidenceText,
        evidence_url = row$EvidenceURL,
        extraction_rule = row$ExtractionRule,
        confidence = row$Confidence,
        matrix_eligible = FALSE
      )
    }
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_traits(out)
}

.normalized_trait_cols = function() {
  c("Query", "CID", "TraitType", "TraitGroup", "TraitValue",
    "TraitValueClean", "TraitLabel", "SourceDatabase", "SourceTable",
    "EvidenceID", "EvidenceText", "EvidenceURL", "ExtractionRule",
    "Confidence", "ConfidenceScore", "MatrixKey", "MatrixEligible")
}

.normalized_add_trait = function(rows, query, cid, trait_type, trait_group,
                                 trait_value, source_database, source_table,
                                 evidence_id = NA_character_,
                                 evidence_text = NA_character_,
                                 evidence_url = NA_character_,
                                 extraction_rule = NA_character_,
                                 confidence = "medium",
                                 matrix_eligible = TRUE,
                                 split = FALSE) {
  values = if (split) {
    .normalized_split_terms(trait_value)
  } else {
    .uaf_non_empty(trait_value)
  }
  if (length(values) < 1) return(rows)

  trait_type = .normalized_trait_key_piece(trait_type)
  trait_group = .normalized_trait_key_piece(trait_group)
  source_database = .uaf_first_non_empty_text(source_database, "unknown")
  source_table = .uaf_first_non_empty_text(source_table, NA_character_)
  evidence_id = .uaf_first_non_empty_text(evidence_id)
  confidence = .uaf_first_non_empty_text(confidence, "medium")

  for (value in values) {
    value = .uaf_squish_text(value)
    if (is.na(value) || value == "") next
    if (.normalized_trait_is_descriptive_value(value)) next
    value_clean = .normalized_trait_value_clean(value)
    if (is.na(value_clean) || value_clean == "") next
    is_matrix_trait = isTRUE(matrix_eligible)
    matrix_key = if (is_matrix_trait) {
      .normalized_trait_matrix_key(trait_type, trait_group, value_clean)
    } else {
      NA_character_
    }
    rows[[length(rows) + 1]] = data.frame(
      Query = query,
      CID = suppressWarnings(as.integer(cid)),
      TraitType = trait_type,
      TraitGroup = trait_group,
      TraitValue = value,
      TraitValueClean = value_clean,
      TraitLabel = .normalized_trait_label(trait_type, trait_group, value),
      SourceDatabase = source_database,
      SourceTable = source_table,
      EvidenceID = evidence_id,
      EvidenceText = evidence_text,
      EvidenceURL = evidence_url,
      ExtractionRule = extraction_rule,
      Confidence = confidence,
      ConfidenceScore = .normalized_trait_confidence_score(
        confidence = confidence,
        source_database = source_database,
        evidence_id = evidence_id,
        extraction_rule = extraction_rule
      ),
      MatrixKey = matrix_key,
      MatrixEligible = is_matrix_trait,
      stringsAsFactors = FALSE
    )
  }
  rows
}

.normalized_trait_is_descriptive_value = function(value) {
  value = .uaf_squish_text(value)
  if (is.na(value) || value == "") return(TRUE)
  if (nchar(value) > 140) return(TRUE)
  word_count = length(strsplit(value, "\\s+", perl = TRUE)[[1]])
  if (word_count > 6 && grepl("\\.$", value)) return(TRUE)
  if (word_count > 7 &&
      grepl("^drugs?\\s+that|^compounds?\\s+or\\s+agents|agents?\\s+that|used\\s+to|\\bconvert\\b",
            value,
            ignore.case = TRUE,
            perl = TRUE)) {
    return(TRUE)
  }
  if (word_count > 14 &&
      grepl("\\b(that|which|their|they|during|following|addition|accounts|used to|acts? by|leading to)\\b",
            value,
            ignore.case = TRUE,
            perl = TRUE)) {
    return(TRUE)
  }
  if (nchar(value) > 90 && grepl("[.;]", value)) return(TRUE)
  FALSE
}

.normalized_trait_key_piece = function(x) {
  x = .uaf_first_non_empty_text(x)
  if (is.na(x)) return(NA_character_)
  x = .normalized_clean_label(x)
  x = gsub("[^a-z0-9]+", "_", x, perl = TRUE)
  x = gsub("^_+|_+$", "", x, perl = TRUE)
  if (is.na(x) || x == "") return("unknown")
  x
}

.normalized_trait_value_clean = function(x) {
  x = .uaf_squish_text(x)
  if (is.na(x) || x == "") return(NA_character_)
  .normalized_clean_label(x)
}

.normalized_trait_label = function(trait_type, trait_group, trait_value) {
  paste(.uaf_non_empty(c(trait_type, trait_group, trait_value)),
        collapse = ":")
}

.normalized_trait_matrix_key = function(trait_type, trait_group,
                                        trait_value_clean) {
  pieces = vapply(c(trait_type, trait_group, trait_value_clean),
                  .normalized_trait_key_piece,
                  character(1))
  paste(pieces, collapse = "__")
}

.normalized_trait_confidence_score = function(confidence, source_database,
                                              evidence_id, extraction_rule) {
  confidence = tolower(.uaf_first_non_empty_text(confidence, "medium"))
  score = if (confidence == "high") {
    0.90
  } else if (confidence == "medium") {
    0.65
  } else if (confidence == "low") {
    0.35
  } else {
    0.50
  }

  source_database = .uaf_first_non_empty_text(source_database)
  if (!is.na(source_database) &&
      source_database %in% c("PubChem", "PubChem BioAssay", "KEGG", "LOTUS",
                             "FEMA", "FDA/SPL", "MeSH")) {
    score = score + 0.03
  }
  if (length(.uaf_non_empty(evidence_id)) > 0) score = score + 0.05
  rule = tolower(.pubchem_collapse(extraction_rule))
  if (!is.na(rule) &&
      grepl("id|kegg|pubchem|lotus|mesh|ghs|assay|equation|property",
            rule, perl = TRUE)) {
    score = score + 0.02
  }
  round(min(score, 1), 2)
}

.normalized_trait_source_database = function(source_table,
                                             source = NA_character_) {
  text = paste(.uaf_non_empty(c(source_table, source)), collapse = " ")
  if (grepl("KEGG", text, ignore.case = TRUE)) return("KEGG")
  if (grepl("LOTUS|natural products", text, ignore.case = TRUE)) return("LOTUS")
  if (grepl("FEMA|Flavor and Extract", text, ignore.case = TRUE)) return("FEMA")
  if (grepl("FDA|SPL", text, ignore.case = TRUE)) return("FDA/SPL")
  if (grepl("MeSH|Medical Subject", text, ignore.case = TRUE)) return("MeSH")
  if (grepl("BioAssay", text, ignore.case = TRUE)) return("PubChem BioAssay")
  if (grepl("PubChem|ChemicalHazards|ChemicalMeasurements|reactives",
            text, ignore.case = TRUE)) {
    return("PubChem")
  }
  .uaf_first_non_empty_text(source, source_table, "unknown")
}

.normalized_trait_cid_lookup = function(tables) {
  rows = list()
  for (table in tables) {
    if (!is.data.frame(table) ||
        nrow(table) < 1 ||
        !"Query" %in% colnames(table) ||
        !"CID" %in% colnames(table)) {
      next
    }
    source = table[!is.na(table$Query) & table$Query != "" &
                     !is.na(table$CID) & table$CID != "",
                   c("Query", "CID"),
                   drop = FALSE]
    if (nrow(source) > 0) rows[[length(rows) + 1]] = source
  }
  if (length(rows) < 1) return(stats::setNames(character(), character()))
  out = unique(do.call(rbind, rows))
  out = out[!duplicated(out$Query), , drop = FALSE]
  stats::setNames(out$CID, out$Query)
}

.normalized_measurement_traits = function(property, value) {
  cols = c("TraitGroup", "TraitValue", "Confidence")
  property_clean = .normalized_trait_key_piece(property)
  value = suppressWarnings(as.numeric(value))
  if (is.na(value)) return(.uaf_empty_table(cols))

  trait = NA_character_
  group = NA_character_
  if (property_clean %in% c("xlogp", "logp")) {
    group = "xlogp_bin"
    trait = if (value < 0) {
      "hydrophilic_XLogP_lt_0"
    } else if (value <= 3) {
      "balanced_XLogP_0_3"
    } else if (value <= 5) {
      "lipophilic_XLogP_3_5"
    } else {
      "high_XLogP_gt_5"
    }
  } else if (property_clean %in% c("tpsa", "topological_polar_surface_area")) {
    group = "tpsa_bin"
    trait = if (value < 40) {
      "low_TPSA_lt_40"
    } else if (value <= 90) {
      "moderate_TPSA_40_90"
    } else {
      "high_TPSA_gt_90"
    }
  } else if (property_clean %in% c("molecularweight", "molecular_weight",
                                   "exactmass", "exact_mass",
                                   "monoisotopicmass", "monoisotopic_mass")) {
    group = "molecular_size_bin"
    trait = if (value < 200) {
      "small_molecule_lt_200_Da"
    } else if (value <= 500) {
      "drug_like_size_200_500_Da"
    } else {
      "large_molecule_gt_500_Da"
    }
  } else if (property_clean == "rotatablebondcount") {
    group = "rotatable_bond_bin"
    trait = if (value <= 3) {
      "low_rotatable_bonds_0_3"
    } else if (value <= 10) {
      "moderate_rotatable_bonds_4_10"
    } else {
      "high_rotatable_bonds_gt_10"
    }
  } else if (property_clean == "hbonddonorcount") {
    group = "hbond_donor_bin"
    trait = if (value == 0) {
      "no_hbond_donors"
    } else if (value <= 2) {
      "low_hbond_donors_1_2"
    } else {
      "high_hbond_donors_gt_2"
    }
  } else if (property_clean == "hbondacceptorcount") {
    group = "hbond_acceptor_bin"
    trait = if (value <= 4) {
      "low_hbond_acceptors_0_4"
    } else if (value <= 10) {
      "moderate_hbond_acceptors_5_10"
    } else {
      "high_hbond_acceptors_gt_10"
    }
  }

  if (is.na(group) || is.na(trait)) return(.uaf_empty_table(cols))
  data.frame(TraitGroup = group,
             TraitValue = trait,
             Confidence = "high",
             stringsAsFactors = FALSE)
}

.normalized_dedupe_traits = function(table) {
  cols = .normalized_trait_cols()
  if (!is.data.frame(table) || nrow(table) < 1) return(.uaf_empty_table(cols))
  key_cols = c("Query", "CID", "TraitType", "TraitGroup", "TraitValueClean")
  key = .pubchem_row_key(table, key_cols, fallback_cols = key_cols)
  rows = lapply(unique(key), function(group_key) {
    group = table[key == group_key, , drop = FALSE]
    row = group[1, , drop = FALSE]
    for (col in c("TraitValue", "TraitLabel", "SourceDatabase", "SourceTable",
                  "EvidenceID", "EvidenceText", "EvidenceURL",
                  "ExtractionRule")) {
      row[[col]] = .pubchem_collapse(group[[col]])
    }
    row$Confidence = .normalized_best_confidence(group$Confidence)
    row$ConfidenceScore = max(suppressWarnings(as.numeric(group$ConfidenceScore)),
                              na.rm = TRUE)
    row$MatrixEligible = any(group$MatrixEligible %in% TRUE, na.rm = TRUE)
    row$MatrixKey = if (isTRUE(row$MatrixEligible)) {
      .normalized_trait_matrix_key(row$TraitType, row$TraitGroup,
                                   row$TraitValueClean)
    } else {
      NA_character_
    }
    row
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out[, cols, drop = FALSE]
}

#' Map normalized chemical traits to a controlled ontology
#'
#' @description
#' `chemicalTraitOntology()` converts extracted `ChemicalTraits` into a
#' conservative, source-backed ontology table. The function only maps
#' whitelisted `TraitType`/`TraitGroup` combinations whose values are already
#' discrete in `ChemicalTraits`; it does not infer new biology or parse long
#' descriptive text.
#'
#' @param traits A `ChemicalTraits` data frame or a categorate result list that
#' contains `ChemicalTraits`.
#' @param min_confidence Minimum confidence score for included traits. Accepts a
#' numeric score or `"low"`, `"medium"`, or `"high"`.
#' @param include_unmapped Logical. If `TRUE`, unmapped discrete traits are kept
#' under `OntologyDomain = "unmapped"` for auditing. Defaults to `FALSE`.
#'
#' @return A data frame with controlled ontology domains, groups, terms,
#' external identifiers when available, source trait keys, source evidence, and
#' confidence metadata.
#'
#' @examples
#' \dontrun{
#' result = categorate(compounds, library_data, detail = "research")
#' ontology = chemicalTraitOntology(result, min_confidence = "medium")
#' }
#'
#' @export
chemicalTraitOntology = function(traits,
                                 min_confidence = 0,
                                 include_unmapped = FALSE) {
  traits = .normalized_prepare_traits(.normalized_extract_traits(traits))
  cols = .normalized_ontology_cols()
  if (!is.data.frame(traits) || nrow(traits) < 1) {
    return(.uaf_empty_table(cols))
  }

  threshold = .normalized_trait_confidence_threshold(min_confidence)
  scores = suppressWarnings(as.numeric(traits$ConfidenceScore))
  scores[is.na(scores)] = 0
  traits = traits[scores >= threshold, , drop = FALSE]
  if (nrow(traits) < 1) return(.uaf_empty_table(cols))

  rows = list()
  for (i in seq_len(nrow(traits))) {
    trait = traits[i, , drop = FALSE]
    rule = .normalized_trait_ontology_rule(trait$TraitType[[1]],
                                           trait$TraitGroup[[1]],
                                           trait$TraitValueClean[[1]])
    if (is.null(rule)) {
      if (!isTRUE(include_unmapped)) next
      rule = list(
        domain = "unmapped",
        group = paste(trait$TraitType[[1]], trait$TraitGroup[[1]], sep = "_"),
        rule = "unmapped_trait_passthrough",
        confidence = "low"
      )
    }
    ontology_term = .normalized_trait_value_clean(trait$TraitValue[[1]])
    if (is.na(ontology_term) || ontology_term == "") next
    source_trait_key = .uaf_first_non_empty_text(
      trait$MatrixKey[[1]],
      .normalized_trait_matrix_key(trait$TraitType[[1]],
                                   trait$TraitGroup[[1]],
                                   trait$TraitValueClean[[1]])
    )
    ontology_id = .normalized_ontology_identifier(trait = trait,
                                                  domain = rule$domain,
                                                  group = rule$group,
                                                  term = ontology_term)
    rows[[length(rows) + 1]] = data.frame(
      Query = trait$Query[[1]],
      CID = suppressWarnings(as.integer(trait$CID[[1]])),
      OntologyDomain = rule$domain,
      OntologyGroup = rule$group,
      OntologyTerm = ontology_term,
      OntologyLabel = trait$TraitValue[[1]],
      OntologyKey = .normalized_trait_ontology_key(rule$domain,
                                                   rule$group,
                                                   ontology_term),
      OntologyID = ontology_id$ID,
      OntologyIDSource = ontology_id$Source,
      OntologyIDURL = ontology_id$URL,
      SourceTraitType = trait$TraitType[[1]],
      SourceTraitGroup = trait$TraitGroup[[1]],
      SourceTraitValue = trait$TraitValue[[1]],
      SourceTraitKey = source_trait_key,
      SourceDatabase = trait$SourceDatabase[[1]],
      SourceTable = trait$SourceTable[[1]],
      EvidenceID = trait$EvidenceID[[1]],
      EvidenceText = trait$EvidenceText[[1]],
      EvidenceURL = trait$EvidenceURL[[1]],
      ExtractionRule = trait$ExtractionRule[[1]],
      Confidence = trait$Confidence[[1]],
      ConfidenceScore = suppressWarnings(as.numeric(trait$ConfidenceScore[[1]])),
      MappingRule = rule$rule,
      MappingConfidence = rule$confidence,
      stringsAsFactors = FALSE
    )
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_ontology(out)
}

#' Build a matrix from controlled chemical trait ontology terms
#'
#' @description
#' `chemicalTraitOntologyMatrix()` converts `ChemicalTraitOntology`, a
#' categorate result, or a `ChemicalTraits` table into a wide matrix whose
#' columns are controlled ontology keys. This gives researchers a more stable
#' analysis surface than source-specific trait names while preserving source
#' identifiers, raw evidence, and confidence metadata in `ChemicalTraits` and
#' `ChemicalTraitOntology`.
#'
#' @param ontology A `ChemicalTraitOntology` data frame, a categorate result, or
#' a `ChemicalTraits` data frame.
#' @param mode Matrix value mode. `"binary"` stores 0/1 presence, `"count"`
#' stores ontology-row counts, and `"confidence"` stores the maximum confidence
#' score for each compound-ontology pair.
#' @param min_confidence Minimum confidence score for included ontology rows.
#' Accepts a numeric score or `"low"`, `"medium"`, or `"high"`.
#' @param max_terms Maximum number of ontology columns to include, ranked by
#' prevalence and confidence.
#'
#' @return A data frame with `Query`, `CID`, and one column per selected
#' ontology key.
#'
#' @examples
#' \dontrun{
#' result = categorate(compounds, library_data, detail = "research")
#' ontology_matrix = chemicalTraitOntologyMatrix(result,
#'                                               min_confidence = "medium")
#' }
#'
#' @export
chemicalTraitOntologyMatrix = function(ontology,
                                       mode = c("binary", "count",
                                                "confidence"),
                                       min_confidence = 0,
                                       max_terms = Inf) {
  mode = match.arg(mode)
  ontology = .normalized_extract_ontology(ontology,
                                          min_confidence = min_confidence)
  ontology = .normalized_prepare_ontology(ontology)
  base_cols = c("Query", "CID")
  if (!is.data.frame(ontology) || nrow(ontology) < 1) {
    return(.uaf_empty_table(base_cols))
  }
  base = unique(ontology[, base_cols, drop = FALSE])
  row.names(base) = NULL

  threshold = .normalized_trait_confidence_threshold(min_confidence)
  scores = suppressWarnings(as.numeric(ontology$ConfidenceScore))
  scores[is.na(scores)] = 0
  source = ontology[scores >= threshold &
                      !is.na(ontology$OntologyKey) &
                      ontology$OntologyKey != "",
                    ,
                    drop = FALSE]
  if (nrow(source) < 1) return(base)

  keys = .normalized_select_ontology_keys(source, max_terms)
  source = source[source$OntologyKey %in% keys, , drop = FALSE]
  zero_value = if (mode == "confidence") 0 else 0L
  for (key in keys) base[[key]] = zero_value
  row_key = paste(.pubchem_key_value(base$Query),
                  .pubchem_key_value(base$CID),
                  sep = "\r")
  source_key = paste(.pubchem_key_value(source$Query),
                     .pubchem_key_value(source$CID),
                     sep = "\r")
  for (key in keys) {
    key_source = source[source$OntologyKey == key, , drop = FALSE]
    key_source_key = source_key[source$OntologyKey == key]
    for (compound_key in unique(key_source_key)) {
      hit_row = match(compound_key, row_key)
      if (is.na(hit_row)) next
      compound_rows = key_source[key_source_key == compound_key, ,
                                 drop = FALSE]
      base[[key]][[hit_row]] = .normalized_trait_matrix_value(compound_rows,
                                                              mode)
    }
  }
  base
}

#' Trace trait and ontology keys back to source evidence
#'
#' @description
#' `chemicalTraitEvidence()` returns an audit table for `ChemicalTraits`,
#' `ChemicalTraitOntology`, or a categorate result. It is designed to explain
#' where matrix or ontology keys came from: source trait keys, source databases,
#' evidence IDs, evidence text/URLs, extraction rules, confidence scores, and
#' external identifiers when available.
#'
#' @param x A categorate result, a `ChemicalTraitOntology` data frame, or a
#' `ChemicalTraits` data frame.
#' @param keys Optional vector of ontology keys, trait matrix keys, source trait
#' keys, external IDs, or terms to keep. Use values such as
#' `"safety__ghs_hazard_code__h319"`, `"hazard__hazard_code__h319"`,
#' `"H319"`, or `"map00590"`.
#' @param type Evidence source. `"auto"` uses `ChemicalTraitOntology` when
#' available and otherwise falls back to `ChemicalTraits`. `"ontology"` returns
#' ontology-key evidence. `"trait"` returns source-trait evidence.
#' @param min_confidence Minimum confidence score for included evidence rows.
#' Accepts a numeric score or `"low"`, `"medium"`, or `"high"`.
#'
#' @return A data frame with analysis keys, source traits, source evidence,
#' external identifiers, extraction rules, and confidence metadata.
#'
#' @examples
#' \dontrun{
#' result = categorate(compounds, library_data, detail = "research")
#' chemicalTraitEvidence(result, keys = "safety__ghs_hazard_code__h319")
#' chemicalTraitEvidence(result,
#'                       keys = "hazard__hazard_code__h319",
#'                       type = "trait")
#' }
#'
#' @export
chemicalTraitEvidence = function(x,
                                 keys = NULL,
                                 type = c("auto", "ontology", "trait"),
                                 min_confidence = 0) {
  type = match.arg(type)
  if (type == "auto") {
    type = .normalized_evidence_auto_type(x)
  }

  evidence = if (type == "trait") {
    .normalized_trait_evidence(
      .normalized_prepare_traits(.normalized_extract_traits(x))
    )
  } else {
    .normalized_ontology_evidence(
      .normalized_prepare_ontology(
        .normalized_extract_ontology(x, min_confidence = min_confidence)
      )
    )
  }

  cols = .normalized_evidence_cols()
  if (!is.data.frame(evidence) || nrow(evidence) < 1) {
    return(.uaf_empty_table(cols))
  }

  threshold = .normalized_trait_confidence_threshold(min_confidence)
  scores = suppressWarnings(as.numeric(evidence$ConfidenceScore))
  scores[is.na(scores)] = 0
  evidence = evidence[scores >= threshold, , drop = FALSE]
  evidence = .normalized_filter_evidence_keys(evidence, keys)
  if (nrow(evidence) < 1) return(.uaf_empty_table(cols))
  evidence = evidence[order(evidence$Query, evidence$EvidenceType,
                            evidence$AnalysisDomain,
                            evidence$AnalysisGroup,
                            evidence$AnalysisTerm),
                      cols,
                      drop = FALSE]
  row.names(evidence) = NULL
  evidence
}

.normalized_evidence_cols = function() {
  c("Query", "CID", "EvidenceType", "AnalysisKey", "AnalysisDomain",
    "AnalysisGroup", "AnalysisTerm", "AnalysisLabel", "ExternalID",
    "ExternalIDSource", "ExternalIDURL", "SourceTraitKey",
    "SourceTraitType", "SourceTraitGroup", "SourceTraitValue",
    "SourceDatabase", "SourceTable", "EvidenceID", "EvidenceText",
    "EvidenceURL", "ExtractionRule", "Confidence", "ConfidenceScore",
    "MappingRule", "MappingConfidence")
}

.normalized_evidence_auto_type = function(x) {
  if (is.list(x) && !is.data.frame(x) &&
      is.data.frame(x$ChemicalTraitOntology)) {
    return("ontology")
  }
  if (is.data.frame(x) && "OntologyDomain" %in% colnames(x)) {
    return("ontology")
  }
  "trait"
}

.normalized_ontology_evidence = function(ontology) {
  cols = .normalized_evidence_cols()
  if (!is.data.frame(ontology) || nrow(ontology) < 1) {
    return(.uaf_empty_table(cols))
  }
  for (col in setdiff(.normalized_ontology_cols(), colnames(ontology))) {
    ontology[[col]] = NA
  }
  out = data.frame(
    Query = ontology$Query,
    CID = suppressWarnings(as.integer(ontology$CID)),
    EvidenceType = "ontology",
    AnalysisKey = ontology$OntologyKey,
    AnalysisDomain = ontology$OntologyDomain,
    AnalysisGroup = ontology$OntologyGroup,
    AnalysisTerm = ontology$OntologyTerm,
    AnalysisLabel = ontology$OntologyLabel,
    ExternalID = ontology$OntologyID,
    ExternalIDSource = ontology$OntologyIDSource,
    ExternalIDURL = ontology$OntologyIDURL,
    SourceTraitKey = ontology$SourceTraitKey,
    SourceTraitType = ontology$SourceTraitType,
    SourceTraitGroup = ontology$SourceTraitGroup,
    SourceTraitValue = ontology$SourceTraitValue,
    SourceDatabase = ontology$SourceDatabase,
    SourceTable = ontology$SourceTable,
    EvidenceID = ontology$EvidenceID,
    EvidenceText = ontology$EvidenceText,
    EvidenceURL = ontology$EvidenceURL,
    ExtractionRule = ontology$ExtractionRule,
    Confidence = ontology$Confidence,
    ConfidenceScore = suppressWarnings(as.numeric(ontology$ConfidenceScore)),
    MappingRule = ontology$MappingRule,
    MappingConfidence = ontology$MappingConfidence,
    stringsAsFactors = FALSE
  )
  out[, cols, drop = FALSE]
}

.normalized_trait_evidence = function(traits) {
  cols = .normalized_evidence_cols()
  if (!is.data.frame(traits) || nrow(traits) < 1) {
    return(.uaf_empty_table(cols))
  }
  for (col in setdiff(.normalized_trait_cols(), colnames(traits))) {
    traits[[col]] = NA
  }
  out = data.frame(
    Query = traits$Query,
    CID = suppressWarnings(as.integer(traits$CID)),
    EvidenceType = "trait",
    AnalysisKey = traits$MatrixKey,
    AnalysisDomain = traits$TraitType,
    AnalysisGroup = traits$TraitGroup,
    AnalysisTerm = traits$TraitValueClean,
    AnalysisLabel = traits$TraitValue,
    ExternalID = NA_character_,
    ExternalIDSource = NA_character_,
    ExternalIDURL = NA_character_,
    SourceTraitKey = traits$MatrixKey,
    SourceTraitType = traits$TraitType,
    SourceTraitGroup = traits$TraitGroup,
    SourceTraitValue = traits$TraitValue,
    SourceDatabase = traits$SourceDatabase,
    SourceTable = traits$SourceTable,
    EvidenceID = traits$EvidenceID,
    EvidenceText = traits$EvidenceText,
    EvidenceURL = traits$EvidenceURL,
    ExtractionRule = traits$ExtractionRule,
    Confidence = traits$Confidence,
    ConfidenceScore = suppressWarnings(as.numeric(traits$ConfidenceScore)),
    MappingRule = NA_character_,
    MappingConfidence = NA_character_,
    stringsAsFactors = FALSE
  )
  out[, cols, drop = FALSE]
}

.normalized_filter_evidence_keys = function(evidence, keys = NULL) {
  keys = .uaf_non_empty(keys)
  if (length(keys) < 1) return(evidence)
  normalized_keys = unique(c(keys, vapply(keys,
                                          .normalized_trait_key_piece,
                                          character(1))))
  keep = vapply(seq_len(nrow(evidence)), function(i) {
    key_cols = c("AnalysisKey", "AnalysisTerm", "AnalysisLabel",
                 "ExternalID", "SourceTraitKey", "SourceTraitValue")
    if (identical(evidence$EvidenceType[[i]], "trait")) {
      key_cols = c(key_cols, "EvidenceID")
    }
    row_values = unlist(evidence[i, key_cols, drop = TRUE],
                        use.names = FALSE)
    row_values = unlist(strsplit(.pubchem_collapse(row_values),
                                 "\\s*;\\s*",
                                 perl = TRUE),
                        use.names = FALSE)
    row_values = .uaf_non_empty(row_values)
    if (length(row_values) < 1) return(FALSE)
    row_values_norm = vapply(row_values,
                             .normalized_trait_key_piece,
                             character(1))
    any(unique(c(row_values, row_values_norm)) %in% normalized_keys)
  }, logical(1))
  evidence[keep, , drop = FALSE]
}

#' Build a researcher-facing chemical trait report
#'
#' @description
#' `chemicalTraitReport()` summarizes the normalized trait system into one
#' compact row per compound. It combines trait counts, ontology domains,
#' evidence coverage, source databases, high-confidence terms, external IDs, and
#' nearest ontology neighbors so results are easier to inspect before filtering,
#' clustering, modeling, or exporting.
#'
#' @param x A categorate result, a `ChemicalTraitOntology` data frame, or a
#' `ChemicalTraits` data frame.
#' @param min_confidence Minimum confidence score for included rows. Accepts a
#' numeric score or `"low"`, `"medium"`, or `"high"`.
#' @param top_n Maximum number of terms, IDs, URLs, and shared keys to collapse
#' into each report cell.
#' @param neighbor_count Maximum number of nearest ontology neighbors to include
#' per compound. Set to `0` to skip neighbor summaries.
#'
#' @return A data frame with one row per compound and compact analysis-ready
#' summaries.
#'
#' @examples
#' \dontrun{
#' result = categorate(compounds, library_data, detail = "research")
#' report = chemicalTraitReport(result, min_confidence = "medium")
#' }
#'
#' @export
chemicalTraitReport = function(x,
                               min_confidence = 0,
                               top_n = 8,
                               neighbor_count = 3) {
  cols = .normalized_report_cols()
  traits = .normalized_prepare_traits(.normalized_extract_traits(x))
  ontology = .normalized_prepare_ontology(
    .normalized_extract_ontology(x, min_confidence = min_confidence)
  )
  evidence = chemicalTraitEvidence(ontology,
                                   type = "ontology",
                                   min_confidence = min_confidence)
  threshold = .normalized_trait_confidence_threshold(min_confidence)
  traits = .normalized_filter_report_confidence(traits, threshold)
  ontology = .normalized_filter_report_confidence(ontology, threshold)
  evidence = .normalized_filter_report_confidence(evidence, threshold)

  base = .normalized_report_base(traits, ontology, evidence)
  if (!is.data.frame(base) || nrow(base) < 1) return(.uaf_empty_table(cols))

  neighbor_info = .normalized_report_neighbors(ontology = ontology,
                                               min_confidence = min_confidence,
                                               top_n = top_n,
                                               neighbor_count = neighbor_count)
  rows = lapply(seq_len(nrow(base)), function(i) {
    query = base$Query[[i]]
    cid = base$CID[[i]]
    trait_rows = .normalized_report_compound_rows(traits, query, cid)
    ontology_rows = .normalized_report_compound_rows(ontology, query, cid)
    evidence_rows = .normalized_report_compound_rows(evidence, query, cid)
    neighbor = neighbor_info[[.normalized_report_compound_key(query, cid)]]
    if (is.null(neighbor)) {
      neighbor = list(NearestNeighbors = NA_character_,
                      NearestNeighborSimilarity = NA_real_,
                      SharedOntologyKeys = NA_character_)
    }
    scores = suppressWarnings(as.numeric(evidence_rows$ConfidenceScore))
    scores = scores[!is.na(scores)]
    trait_ids = if (nrow(trait_rows) > 0) {
      paste(.pubchem_key_value(trait_rows$TraitType),
            .pubchem_key_value(trait_rows$TraitGroup),
            .pubchem_key_value(trait_rows$TraitValueClean),
            sep = "\r")
    } else {
      .uaf_non_empty(ontology_rows$SourceTraitKey)
    }
    data.frame(
      Query = query,
      CID = suppressWarnings(as.integer(cid)),
      TraitCount = length(unique(.uaf_non_empty(trait_ids))),
      OntologyTermCount = length(unique(.uaf_non_empty(
        ontology_rows$OntologyKey
      ))),
      OntologyDomainCount = length(unique(.uaf_non_empty(
        ontology_rows$OntologyDomain
      ))),
      EvidenceCount = nrow(evidence_rows),
      HighConfidenceEvidenceCount = sum(scores >= 0.90),
      MeanConfidence = if (length(scores) > 0) round(mean(scores), 3) else NA_real_,
      SourceDatabases = .pubchem_collapse(ontology_rows$SourceDatabase),
      OntologyDomains = .pubchem_collapse(ontology_rows$OntologyDomain),
      SafetySignals = .normalized_report_domain_terms(ontology_rows,
                                                      "safety", top_n),
      PhysicochemicalSignals = .normalized_report_domain_terms(
        ontology_rows, "physicochemical", top_n
      ),
      ExposureSignals = .normalized_report_domain_terms(ontology_rows,
                                                        "exposure", top_n),
      MetabolismSignals = .normalized_report_domain_terms(ontology_rows,
                                                          "metabolism", top_n),
      EcologySignals = .normalized_report_domain_terms(ontology_rows,
                                                       "ecology", top_n),
      BioactivitySignals = .normalized_report_domain_terms(ontology_rows,
                                                           "bioactivity",
                                                           top_n),
      SensorySignals = .normalized_report_domain_terms(ontology_rows,
                                                       "sensory", top_n),
      BiomedicalSignals = .normalized_report_domain_terms(ontology_rows,
                                                          "biomedical",
                                                          top_n),
      RegulatorySignals = .normalized_report_domain_terms(ontology_rows,
                                                          "regulatory",
                                                          top_n),
      ReactivitySignals = .normalized_report_domain_terms(ontology_rows,
                                                          "reactivity",
                                                          top_n),
      ExternalIDs = .normalized_report_external_ids(ontology_rows, top_n),
      EvidenceURLs = .normalized_report_urls(evidence_rows, top_n),
      NearestNeighbors = neighbor$NearestNeighbors,
      NearestNeighborSimilarity = neighbor$NearestNeighborSimilarity,
      SharedOntologyKeys = neighbor$SharedOntologyKeys,
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out[order(out$Query), cols, drop = FALSE]
}

.normalized_report_cols = function() {
  c("Query", "CID", "TraitCount", "OntologyTermCount",
    "OntologyDomainCount", "EvidenceCount", "HighConfidenceEvidenceCount",
    "MeanConfidence", "SourceDatabases", "OntologyDomains", "SafetySignals",
    "PhysicochemicalSignals", "ExposureSignals", "MetabolismSignals",
    "EcologySignals", "BioactivitySignals", "SensorySignals",
    "BiomedicalSignals", "RegulatorySignals", "ReactivitySignals",
    "ExternalIDs", "EvidenceURLs", "NearestNeighbors",
    "NearestNeighborSimilarity", "SharedOntologyKeys")
}

.normalized_filter_report_confidence = function(table, threshold) {
  if (!is.data.frame(table) || nrow(table) < 1 ||
      !"ConfidenceScore" %in% colnames(table)) {
    return(table)
  }
  scores = suppressWarnings(as.numeric(table$ConfidenceScore))
  scores[is.na(scores)] = 0
  table[scores >= threshold, , drop = FALSE]
}

.normalized_report_base = function(traits, ontology, evidence) {
  rows = list()
  for (table in list(traits, ontology, evidence)) {
    if (!is.data.frame(table) || nrow(table) < 1 ||
        !"Query" %in% colnames(table) ||
        !"CID" %in% colnames(table)) {
      next
    }
    rows[[length(rows) + 1]] = table[, c("Query", "CID"), drop = FALSE]
  }
  if (length(rows) < 1) return(.uaf_empty_table(c("Query", "CID")))
  out = unique(do.call(rbind, rows))
  row.names(out) = NULL
  out
}

.normalized_report_compound_key = function(query, cid) {
  paste(.pubchem_key_value(query), .pubchem_key_value(cid), sep = "\r")
}

.normalized_report_compound_rows = function(table, query, cid) {
  if (!is.data.frame(table) || nrow(table) < 1 ||
      !"Query" %in% colnames(table) ||
      !"CID" %in% colnames(table)) {
    return(.uaf_empty_table(colnames(table)))
  }
  key = .normalized_report_compound_key(table$Query, table$CID)
  table[key == .normalized_report_compound_key(query, cid), , drop = FALSE]
}

.normalized_report_domain_terms = function(ontology_rows, domain, top_n) {
  if (!is.data.frame(ontology_rows) || nrow(ontology_rows) < 1) {
    return(NA_character_)
  }
  rows = ontology_rows[ontology_rows$OntologyDomain == domain, , drop = FALSE]
  if (nrow(rows) < 1) return(NA_character_)
  labels = .uaf_first_non_empty_vector(rows$OntologyLabel,
                                       rows$OntologyTerm)
  groups = .uaf_first_non_empty_vector(rows$OntologyGroup,
                                       rep(domain, nrow(rows)))
  values = paste(groups, labels, sep = "=")
  .normalized_top_values(values, top_n)
}

.normalized_report_external_ids = function(ontology_rows, top_n) {
  if (!is.data.frame(ontology_rows) || nrow(ontology_rows) < 1) {
    return(NA_character_)
  }
  ids = .uaf_non_empty(ontology_rows$OntologyID)
  if (length(ids) < 1) return(NA_character_)
  sources = ontology_rows$OntologyIDSource[
    !is.na(ontology_rows$OntologyID) & ontology_rows$OntologyID != ""
  ]
  values = ifelse(!is.na(sources) & sources != "",
                  paste(sources, ids, sep = ":"),
                  ids)
  .normalized_top_values(values, top_n)
}

.normalized_report_urls = function(evidence_rows, top_n) {
  if (!is.data.frame(evidence_rows) || nrow(evidence_rows) < 1) {
    return(NA_character_)
  }
  values = unlist(strsplit(.pubchem_collapse(c(evidence_rows$ExternalIDURL,
                                               evidence_rows$EvidenceURL)),
                           "\\s*;\\s*",
                           perl = TRUE),
                  use.names = FALSE)
  .normalized_top_values(values, top_n)
}

.normalized_report_neighbors = function(ontology, min_confidence, top_n,
                                        neighbor_count) {
  neighbor_count = suppressWarnings(as.integer(neighbor_count[[1]]))
  if (is.na(neighbor_count) || neighbor_count < 1 ||
      !is.data.frame(ontology) || nrow(ontology) < 1) {
    return(list())
  }
  matrix = chemicalTraitOntologyMatrix(ontology,
                                       mode = "binary",
                                       min_confidence = min_confidence)
  similarity = chemicalTraitSimilarity(matrix)
  if (!is.data.frame(similarity) || nrow(similarity) < 1) return(list())
  out = list()
  base = unique(matrix[, c("Query", "CID"), drop = FALSE])
  for (i in seq_len(nrow(base))) {
    query = base$Query[[i]]
    cid = base$CID[[i]]
    compound_key = .normalized_report_compound_key(query, cid)
    key_a = .normalized_report_compound_key(similarity$QueryA,
                                            similarity$CIDA)
    key_b = .normalized_report_compound_key(similarity$QueryB,
                                            similarity$CIDB)
    rows = similarity[
      key_a == compound_key | key_b == compound_key,
      ,
      drop = FALSE
    ]
    if (nrow(rows) < 1) next
    rows = rows[order(-rows$JaccardSimilarity, -rows$SharedTraitCount),
                ,
                drop = FALSE]
    rows = utils::head(rows, neighbor_count)
    neighbor_names = vapply(seq_len(nrow(rows)), function(j) {
      row_key_a = .normalized_report_compound_key(rows$QueryA[[j]],
                                                  rows$CIDA[[j]])
      other = if (row_key_a == compound_key) {
        rows$QueryB[[j]]
      } else {
        rows$QueryA[[j]]
      }
      paste0(other, " (", rows$JaccardSimilarity[[j]], ")")
    }, character(1))
    shared_keys = unique(unlist(strsplit(.pubchem_collapse(
      utils::head(rows$SharedTraits, top_n)
    ), "\\s*;\\s*", perl = TRUE), use.names = FALSE))
    out[[.normalized_report_compound_key(query, cid)]] = list(
      NearestNeighbors = .pubchem_collapse(neighbor_names),
      NearestNeighborSimilarity = rows$JaccardSimilarity[[1]],
      SharedOntologyKeys = .pubchem_collapse(utils::head(.uaf_non_empty(shared_keys),
                                                        top_n))
    )
  }
  out
}

.uaf_first_non_empty_vector = function(primary, fallback) {
  primary = as.character(primary)
  fallback = as.character(fallback)
  missing = is.na(primary) | primary == ""
  primary[missing] = fallback[missing]
  primary
}

.normalized_ontology_cols = function() {
  c("Query", "CID", "OntologyDomain", "OntologyGroup", "OntologyTerm",
    "OntologyLabel", "OntologyKey", "OntologyID", "OntologyIDSource",
    "OntologyIDURL", "SourceTraitType", "SourceTraitGroup",
    "SourceTraitValue", "SourceTraitKey", "SourceDatabase", "SourceTable",
    "EvidenceID", "EvidenceText", "EvidenceURL", "ExtractionRule",
    "Confidence", "ConfidenceScore", "MappingRule", "MappingConfidence")
}

.normalized_trait_ontology_rule = function(trait_type, trait_group,
                                           trait_value_clean) {
  trait_type = .normalized_trait_key_piece(trait_type)
  trait_group = .normalized_trait_key_piece(trait_group)
  trait_value_clean = .normalized_trait_key_piece(trait_value_clean)
  key = paste(trait_type, trait_group, sep = ":")
  direct = list(
    "property:molecular_size_bin" = c("physicochemical", "molecular_size",
                                     "property_bin", "high"),
    "property:xlogp_bin" = c("physicochemical", "lipophilicity",
                             "property_bin", "high"),
    "property:tpsa_bin" = c("physicochemical", "polarity",
                            "property_bin", "high"),
    "property:rotatable_bond_bin" = c("physicochemical", "flexibility",
                                      "property_bin", "high"),
    "property:hbond_donor_bin" = c("physicochemical",
                                   "hbond_donor_capacity",
                                   "property_bin", "high"),
    "property:hbond_acceptor_bin" = c("physicochemical",
                                      "hbond_acceptor_capacity",
                                      "property_bin", "high"),
    "hazard:hazard_code" = c("safety", "ghs_hazard_code",
                             "ghs_code", "high"),
    "hazard:hazard_category" = c("safety", "hazard_category",
                                 "ghs_category", "high"),
    "hazard:hazard_group" = c("safety", "hazard_endpoint",
                              "hazard_group", "high"),
    "hazard:signal_word" = c("safety", "signal_word",
                             "ghs_signal_word", "high"),
    "hazard:exposure_route" = c("exposure", "exposure_route",
                                "route", "high"),
    "hazard:target_organ" = c("safety", "affected_organ",
                              "target_organ", "high"),
    "hazard:toxicity_metric" = c("safety", "toxicity_metric",
                                 "toxicity_measure", "high"),
    "hazard:toxicity_species" = c("safety", "toxicity_species",
                                  "toxicity_species", "high"),
    "use:descriptor" = c("sensory", "sensory_descriptor",
                         "fema_descriptor", "high"),
    "use:flavor" = c("sensory", "flavor",
                     "fema_flavor", "high"),
    "use:odor" = c("sensory", "odor",
                   "fema_odor", "high"),
    "use:regulatory_status" = c("regulatory",
                                "flavor_regulatory_status",
                                "fema_regulatory_status", "high"),
    "use:active_ingredient" = c("biomedical", "active_ingredient",
                                "fda_spl_active_ingredient", "high"),
    "use:pharmacologic_class" = c("biomedical", "pharmacologic_class",
                                  "pharmacologic_class", "high"),
    "use:pharmacologic_action" = c("biomedical", "pharmacologic_action",
                                   "pharmacologic_action", "high"),
    "use:administration_route" = c("biomedical", "administration_route",
                                   "fda_spl_route", "high"),
    "use:dosage_form" = c("biomedical", "dosage_form",
                          "fda_spl_dosage_form", "high"),
    "use:warning" = c("safety", "label_warning",
                      "fda_spl_warning", "high"),
    "use:use_group" = c("use", "use_group",
                        "normalized_use_group", "medium"),
    "sensory:descriptor" = c("sensory", "sensory_descriptor",
                             "sensory_descriptor", "high"),
    "sensory:flavor" = c("sensory", "flavor",
                         "sensory_flavor", "high"),
    "sensory:odor" = c("sensory", "odor",
                       "sensory_odor", "high"),
    "chemical_class:reactive_group" = c("reactivity", "reactive_group",
                                        "pubchem_reactive_group", "high"),
    "chemical_class:natural_product_superclass" = c("ecology",
                                                    "natural_product_superclass",
                                                    "lotus_chemical_tree",
                                                    "high"),
    "chemical_class:natural_product_class" = c("ecology",
                                               "natural_product_class",
                                               "lotus_chemical_tree",
                                               "high"),
    "chemical_class:natural_product_subclass" = c("ecology",
                                                  "natural_product_subclass",
                                                  "lotus_chemical_tree",
                                                  "high"),
    "chemical_class:pathway_group" = c("metabolism", "pathway_group",
                                       "kegg_pathway_group", "high"),
    "chemical_class:enzyme_class" = c("metabolism", "enzyme_class",
                                      "kegg_enzyme_class", "high"),
    "chemical_class:mesh_tree_category" = c("biomedical",
                                            "mesh_tree_category",
                                            "mesh_tree_number", "high"),
    "bioactivity:activity_class" = c("bioactivity", "assay_outcome",
                                     "pubchem_assay_outcome", "high"),
    "bioactivity:activity_direction" = c("bioactivity",
                                         "activity_direction",
                                         "pubchem_assay_direction", "high"),
    "bioactivity:bioactivity_domain" = c("bioactivity", "assay_domain",
                                         "pubchem_assay_domain", "high"),
    "bioactivity:assay_type" = c("bioactivity", "assay_type",
                                 "pubchem_assay_type", "high"),
    "bioactivity:endpoint_metric" = c("bioactivity", "endpoint_metric",
                                      "pubchem_endpoint_metric", "high"),
    "target:target_name" = c("bioactivity", "target_name",
                             "pubchem_target_name", "high"),
    "target:target_type" = c("bioactivity", "target_type",
                             "pubchem_target_type", "high"),
    "target:target_accession" = c("bioactivity", "target_accession",
                                  "pubchem_target_accession", "high"),
    "target:target_gene" = c("bioactivity", "target_gene",
                             "pubchem_target_gene", "high"),
    "target:target_taxonomy" = c("bioactivity", "target_taxonomy",
                                 "pubchem_target_taxonomy", "high"),
    "target:target_organism" = c("bioactivity", "target_organism",
                                 "pubchem_target_organism", "high"),
    "target:target_common_name" = c("bioactivity", "target_common_name",
                                    "pubchem_target_common_name", "high"),
    "target:target_group" = c("bioactivity", "target_group",
                              "pubchem_target_group", "high"),
    "target:endpoint_metric" = c("bioactivity", "endpoint_metric",
                                 "pubchem_endpoint_metric", "high"),
    "potency:potency_metric" = c("bioactivity", "potency_metric",
                                 "pubchem_potency_metric", "high"),
    "potency:potency_bucket" = c("bioactivity", "potency_bucket",
                                 "pubchem_potency_bucket", "high"),
    "potency:target_name" = c("bioactivity", "target_name",
                              "pubchem_potency_target", "high"),
    "pathway:pathway_group" = c("metabolism", "pathway_group",
                                "kegg_pathway_group", "high"),
    "pathway:pathway_name" = c("metabolism", "pathway_name",
                               "kegg_pathway_name", "high"),
    "pathway:reaction_name" = c("metabolism", "reaction_name",
                                "kegg_reaction_name", "high"),
    "pathway:enzyme_class" = c("metabolism", "enzyme_class",
                               "kegg_enzyme_class", "high"),
    "pathway:ec_number" = c("metabolism", "ec_number",
                            "kegg_ec_number", "high"),
    "pathway:role_type" = c("metabolism", "pathway_role",
                            "kegg_pathway_role", "high"),
    "reaction:participant_role" = c("metabolism",
                                    "reaction_participant_role",
                                    "kegg_reaction_participant", "high"),
    "reaction:participant_side" = c("metabolism",
                                    "reaction_participant_side",
                                    "kegg_reaction_side", "high"),
    "reaction:participant_id" = c("metabolism",
                                  "reaction_participant_id",
                                  "kegg_reaction_participant_id", "high")
  )
  hit = direct[[key]]
  if (!is.null(hit)) {
    return(list(domain = hit[[1]], group = hit[[2]], rule = hit[[3]],
                confidence = hit[[4]]))
  }

  taxonomy_ranks = c("domain", "kingdom", "phylum", "class", "order",
                     "family", "genus", "species", "organism",
                     "common_name", "taxonomy_group",
                     "natural_product_class")
  if (trait_type == "taxonomy" && trait_group %in% taxonomy_ranks) {
    return(list(domain = "ecology",
                group = paste0("taxonomic_", trait_group),
                rule = "lotus_taxonomy_rank",
                confidence = "high"))
  }
  occurrence_ranks = c("kingdom", "family", "genus", "species", "organism",
                       "natural_product_class")
  if (trait_type == "occurrence" && trait_group %in% occurrence_ranks) {
    return(list(domain = "ecology",
                group = paste0("taxonomic_", trait_group),
                rule = "lotus_occurrence_rank",
                confidence = "high"))
  }
  if (trait_type == "occurrence" && trait_group == "occurrence_type") {
    return(list(domain = "ecology",
                group = "occurrence_type",
                rule = "lotus_occurrence_type",
                confidence = "high"))
  }
  NULL
}

.normalized_trait_ontology_key = function(domain, group, term) {
  paste(vapply(c(domain, group, term),
               .normalized_trait_key_piece,
               character(1)),
        collapse = "__")
}

.normalized_ontology_identifier = function(trait, domain, group, term) {
  out = list(ID = NA_character_, Source = NA_character_,
             URL = NA_character_)
  source_type = .uaf_first_non_empty_text(trait$TraitType[[1]])
  source_group = .uaf_first_non_empty_text(trait$TraitGroup[[1]])
  source_value = .uaf_first_non_empty_text(trait$TraitValue[[1]])
  evidence_id = .uaf_first_non_empty_text(trait$EvidenceID[[1]])
  evidence_text = .uaf_first_non_empty_text(trait$EvidenceText[[1]])
  text = .pubchem_collapse(c(evidence_id, source_value, evidence_text))

  if (group == "ghs_hazard_code") {
    id = .normalized_first_identifier(text, "\\bH[0-9]{3}[A-Za-z]?\\b")
    if (!is.na(id)) {
      out$ID = toupper(id)
      out$Source = "GHS"
      return(out)
    }
  }

  if (source_group == "mesh_tree_category") {
    id = .normalized_first_identifier(text,
                                      "\\b[A-Z][0-9]{2}(?:\\.[0-9]+)+\\b")
    if (!is.na(id)) {
      out$ID = id
      out$Source = "MeSH tree number"
      return(out)
    }
  }

  if (source_type %in% c("taxonomy", "occurrence")) {
    id = .normalized_first_identifier(evidence_id, "\\b[0-9]{1,12}\\b")
    if (!is.na(id)) {
      out$ID = id
      out$Source = "NCBI Taxonomy"
      out$URL = paste0("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/",
                       "wwwtax.cgi?id=", id)
      return(out)
    }
  }

  if (source_type %in% c("pathway", "reaction", "chemical_class")) {
    ec = .normalized_first_identifier(text, "\\b[0-9]+(?:\\.[0-9-]+){3}\\b")
    if (source_group %in% c("ec_number", "enzyme_class") && !is.na(ec)) {
      out$ID = paste0("ec:", ec)
      out$Source = "EC"
      out$URL = .normalized_kegg_entry_url(out$ID)
      return(out)
    }

    pathway = .normalized_first_identifier(text, "\\b(?:map|ko)[0-9]{5}\\b")
    if (source_group %in% c("pathway_group", "pathway_name") &&
        !is.na(pathway)) {
      out$ID = pathway
      out$Source = "KEGG PATHWAY"
      out$URL = .normalized_kegg_entry_url(pathway)
      return(out)
    }

    reaction = .normalized_first_identifier(text, "\\bR[0-9]{5}\\b")
    if (source_group %in% c("reaction_name", "participant_role",
                            "participant_side") &&
        !is.na(reaction)) {
      out$ID = reaction
      out$Source = "KEGG REACTION"
      out$URL = .normalized_kegg_entry_url(reaction)
      return(out)
    }

    compound = .normalized_first_identifier(source_value, "\\bC[0-9]{5}\\b")
    if (source_group == "participant_id" && !is.na(compound)) {
      out$ID = compound
      out$Source = "KEGG COMPOUND"
      out$URL = .normalized_kegg_entry_url(compound)
      return(out)
    }
  }

  if (source_type %in% c("bioactivity", "target", "potency")) {
    if (source_group == "target_gene") {
      id = .normalized_first_identifier(source_value, "\\b[0-9]{1,12}\\b")
      if (!is.na(id)) {
        out$ID = id
        out$Source = "NCBI Gene"
        out$URL = paste0("https://www.ncbi.nlm.nih.gov/gene/", id)
        return(out)
      }
    }
    if (source_group == "target_taxonomy") {
      id = .normalized_first_identifier(source_value, "\\b[0-9]{1,12}\\b")
      if (!is.na(id)) {
        out$ID = id
        out$Source = "NCBI Taxonomy"
        out$URL = paste0("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/",
                         "wwwtax.cgi?id=", id)
        return(out)
      }
    }
    if (source_group == "target_accession") {
      id = .uaf_first_non_empty_text(source_value)
      if (!is.na(id)) {
        out$ID = id
        out$Source = "PubChem BioAssay target accession"
        return(out)
      }
    }
    aid = .normalized_first_identifier(evidence_id, "\\b[0-9]{1,12}\\b")
    if (source_type %in% c("bioactivity", "potency") &&
        !is.na(aid) &&
        source_group %in% c("activity_class", "activity_direction",
                            "bioactivity_domain", "assay_type",
                            "endpoint_metric", "potency_metric",
                            "potency_bucket")) {
      out$ID = paste0("AID:", aid)
      out$Source = "PubChem BioAssay"
      out$URL = paste0("https://pubchem.ncbi.nlm.nih.gov/bioassay/", aid)
      return(out)
    }
  }

  out
}

.normalized_first_identifier = function(text, pattern) {
  text = .pubchem_collapse(text)
  if (is.na(text) || text == "") return(NA_character_)
  hits = regmatches(text, gregexpr(pattern, text, perl = TRUE,
                                   ignore.case = FALSE))[[1]]
  if (length(hits) < 1 ||
      identical(hits, character(0)) ||
      (length(hits) == 1 && hits[[1]] == "")) {
    return(NA_character_)
  }
  hits[[1]]
}

.normalized_kegg_entry_url = function(id) {
  id = .uaf_first_non_empty_text(id)
  if (is.na(id)) return(NA_character_)
  paste0("https://www.kegg.jp/entry/", utils::URLencode(id, reserved = FALSE))
}

.normalized_dedupe_ontology = function(table) {
  cols = .normalized_ontology_cols()
  if (!is.data.frame(table) || nrow(table) < 1) {
    return(.uaf_empty_table(cols))
  }
  key_cols = c("Query", "CID", "OntologyDomain", "OntologyGroup",
               "OntologyTerm")
  key = .pubchem_row_key(table, key_cols, fallback_cols = key_cols)
  rows = lapply(unique(key), function(group_key) {
    group = table[key == group_key, , drop = FALSE]
    row = group[1, , drop = FALSE]
    for (col in c("OntologyLabel", "SourceTraitType", "SourceTraitGroup",
                  "SourceTraitValue", "SourceTraitKey", "OntologyID",
                  "OntologyIDSource", "OntologyIDURL", "SourceDatabase",
                  "SourceTable", "EvidenceID", "EvidenceText",
                  "EvidenceURL", "ExtractionRule", "MappingRule")) {
      row[[col]] = .pubchem_collapse(group[[col]])
    }
    row$Confidence = .normalized_best_confidence(group$Confidence)
    row$ConfidenceScore = max(suppressWarnings(as.numeric(group$ConfidenceScore)),
                              na.rm = TRUE)
    row$MappingConfidence =
      .normalized_best_confidence(group$MappingConfidence)
    row$OntologyKey = .normalized_trait_ontology_key(row$OntologyDomain,
                                                     row$OntologyGroup,
                                                     row$OntologyTerm)
    row
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out[order(out$Query, out$OntologyDomain, out$OntologyGroup,
            out$OntologyTerm),
      cols,
      drop = FALSE]
}

.normalized_extract_ontology = function(ontology, min_confidence = 0) {
  if (is.list(ontology) && !is.data.frame(ontology)) {
    if (is.data.frame(ontology$ChemicalTraitOntology)) {
      return(ontology$ChemicalTraitOntology)
    }
    if (is.data.frame(ontology$ChemicalTraits)) {
      return(chemicalTraitOntology(ontology$ChemicalTraits,
                                   min_confidence = min_confidence))
    }
  }
  if (is.data.frame(ontology) &&
      "OntologyDomain" %in% colnames(ontology)) {
    return(ontology)
  }
  if (is.data.frame(ontology)) {
    return(chemicalTraitOntology(ontology,
                                 min_confidence = min_confidence))
  }
  .uaf_empty_table(.normalized_ontology_cols())
}

.normalized_prepare_ontology = function(ontology) {
  cols = .normalized_ontology_cols()
  if (!is.data.frame(ontology) || nrow(ontology) < 1) {
    return(.uaf_empty_table(cols))
  }
  for (col in setdiff(cols, colnames(ontology))) ontology[[col]] = NA
  ontology$OntologyDomain = vapply(ontology$OntologyDomain,
                                   .normalized_trait_key_piece,
                                   character(1))
  ontology$OntologyGroup = vapply(ontology$OntologyGroup,
                                  .normalized_trait_key_piece,
                                  character(1))
  ontology$OntologyTerm = vapply(ontology$OntologyTerm,
                                 .normalized_trait_value_clean,
                                 character(1))
  ontology$OntologyLabel = .uaf_squish_text(ontology$OntologyLabel)
  missing_label = is.na(ontology$OntologyLabel) |
    ontology$OntologyLabel == ""
  ontology$OntologyLabel[missing_label] = ontology$OntologyTerm[missing_label]
  missing_key = is.na(ontology$OntologyKey) | ontology$OntologyKey == ""
  if (any(missing_key)) {
    ontology$OntologyKey[missing_key] = mapply(
      .normalized_trait_ontology_key,
      ontology$OntologyDomain[missing_key],
      ontology$OntologyGroup[missing_key],
      ontology$OntologyTerm[missing_key],
      USE.NAMES = FALSE
    )
  }
  ontology$Confidence = vapply(ontology$Confidence, function(x) {
    .uaf_first_non_empty_text(x, "medium")
  }, character(1))
  scores = suppressWarnings(as.numeric(ontology$ConfidenceScore))
  missing_score = is.na(scores)
  if (any(missing_score)) {
    scores[missing_score] = mapply(
      .normalized_trait_confidence_score,
      confidence = ontology$Confidence[missing_score],
      source_database = ontology$SourceDatabase[missing_score],
      evidence_id = ontology$EvidenceID[missing_score],
      extraction_rule = ontology$ExtractionRule[missing_score],
      USE.NAMES = FALSE
    )
  }
  ontology$ConfidenceScore = scores
  ontology[, cols, drop = FALSE]
}

.normalized_select_ontology_keys = function(source, max_terms) {
  keys = unique(source$OntologyKey)
  max_value = suppressWarnings(as.numeric(max_terms[[1]]))
  if (is.na(max_value) || max_value < 1) return(character())
  max_terms = if (is.infinite(max_value)) {
    length(keys)
  } else {
    as.integer(max_value)
  }
  source_key = paste(.pubchem_key_value(source$Query),
                     .pubchem_key_value(source$CID),
                     sep = "\r")
  prevalence = vapply(keys, function(key) {
    length(unique(source_key[source$OntologyKey == key]))
  }, integer(1))
  confidence = vapply(keys, function(key) {
    values = suppressWarnings(as.numeric(source$ConfidenceScore[
      source$OntologyKey == key
    ]))
    values = values[!is.na(values)]
    if (length(values) < 1) return(0)
    mean(values)
  }, numeric(1))
  ranked = data.frame(OntologyKey = keys,
                      Prevalence = prevalence,
                      Confidence = confidence,
                      stringsAsFactors = FALSE)
  ranked = ranked[order(-ranked$Prevalence, -ranked$Confidence,
                        ranked$OntologyKey),
                  ,
                  drop = FALSE]
  utils::head(ranked$OntologyKey, max_terms)
}

#' Build a matrix from normalized chemical traits
#'
#' @description
#' `chemicalTraitMatrix()` converts a `ChemicalTraits` table, or a
#' `categorate(detail = "research")`/`categorate(detail = "full")` result, into
#' a wide matrix for filtering, clustering, ordination, heatmaps, and model
#' inputs. The long `ChemicalTraits` table remains the complete evidence table;
#' this helper lets users choose compact or source-specific matrix views without
#' rerunning web queries.
#'
#' @param traits A `ChemicalTraits` data frame or a categorate result list that
#' contains `ChemicalTraits`.
#' @param profile Matrix profile. `"core"` keeps compact cross-domain grouping
#' traits. `"full"` keeps every eligible trait. Other profiles keep traits for
#' one analysis domain: `"bioactivity"`, `"safety"`, `"ecology"`, `"kegg"`,
#' `"sensory"`, or `"biomedical"`.
#' @param mode Matrix value mode. `"binary"` stores 0/1 presence, `"count"`
#' stores trait-row counts, and `"confidence"` stores the maximum confidence
#' score for each compound-trait pair.
#' @param min_confidence Minimum confidence score for included traits. Accepts a
#' numeric score or `"low"`, `"medium"`, or `"high"`.
#' @param max_traits Maximum number of trait columns to include, ranked by
#' prevalence and confidence.
#'
#' @return A data frame with `Query`, `CID`, and one column per selected trait.
#'
#' @examples
#' \dontrun{
#' result = categorate(compounds, library_data, detail = "full")
#' core_matrix = result$ChemicalTraitMatrix
#' bioactivity_matrix = chemicalTraitMatrix(result, profile = "bioactivity")
#' confidence_matrix = chemicalTraitMatrix(result$ChemicalTraits,
#'                                        profile = "full",
#'                                        mode = "confidence",
#'                                        min_confidence = "high")
#' }
#'
#' @export
chemicalTraitMatrix = function(traits,
                               profile = c("core", "full", "bioactivity",
                                           "safety", "ecology", "kegg",
                                           "sensory", "biomedical"),
                               mode = c("binary", "count", "confidence"),
                               min_confidence = 0,
                               max_traits = Inf) {
  profile = match.arg(profile)
  mode = match.arg(mode)
  traits = .normalized_extract_traits(traits)
  .normalized_chemical_trait_matrix(chemical_traits = traits,
                                    profile = profile,
                                    mode = mode,
                                    min_confidence = min_confidence,
                                    max_traits = max_traits)
}

#' Summarize normalized chemical traits
#'
#' @description
#' `chemicalTraitSummary()` produces compact per-compound summaries from a
#' `ChemicalTraits` table or a categorate result. Use it to see how many traits
#' were extracted by domain or source, which sources contributed evidence, and
#' which trait groups/values dominate each compound.
#'
#' @param traits A `ChemicalTraits` data frame or a categorate result list that
#' contains `ChemicalTraits`.
#' @param by Summary level: `"type"` summarizes by `TraitType`, `"source"`
#' summarizes by `SourceDatabase`, and `"compound"` summarizes all traits per
#' compound.
#' @param min_confidence Minimum confidence score for included traits. Accepts a
#' numeric score or `"low"`, `"medium"`, or `"high"`.
#' @param top_n Number of top trait groups and values to include in collapsed
#' summary columns.
#'
#' @return A data frame with counts, confidence summaries, source coverage, and
#' top trait groups/values.
#'
#' @examples
#' \dontrun{
#' result = categorate(compounds, library_data, detail = "full")
#' chemicalTraitSummary(result)
#' chemicalTraitSummary(result, by = "source", min_confidence = "high")
#' }
#'
#' @export
chemicalTraitSummary = function(traits,
                                by = c("type", "source", "compound"),
                                min_confidence = 0,
                                top_n = 8) {
  by = match.arg(by)
  traits = .normalized_prepare_traits(.normalized_extract_traits(traits))
  cols = c("Query", "CID", "SummaryLevel", "TraitType", "SourceDatabase",
           "TraitCount", "MatrixEligibleCount", "HighConfidenceCount",
           "MeanConfidence", "SourceDatabases", "TopTraitGroups",
           "TopTraitValues")
  if (!is.data.frame(traits) || nrow(traits) < 1) return(.uaf_empty_table(cols))

  threshold = .normalized_trait_confidence_threshold(min_confidence)
  scores = suppressWarnings(as.numeric(traits$ConfidenceScore))
  scores[is.na(scores)] = 0
  traits = traits[scores >= threshold, , drop = FALSE]
  if (nrow(traits) < 1) return(.uaf_empty_table(cols))

  summary_values = if (by == "type") {
    traits$TraitType
  } else if (by == "source") {
    traits$SourceDatabase
  } else {
    rep("all", nrow(traits))
  }
  summary_values = .pubchem_key_value(summary_values)
  group_keys = paste(.pubchem_key_value(traits$Query),
                     .pubchem_key_value(traits$CID),
                     summary_values,
                     sep = "\r")

  rows = lapply(unique(group_keys), function(group_key) {
    group = traits[group_keys == group_key, , drop = FALSE]
    trait_id = paste(.pubchem_key_value(group$TraitType),
                     .pubchem_key_value(group$TraitGroup),
                     .pubchem_key_value(group$TraitValueClean),
                     sep = "\r")
    matrix_keys = unique(.uaf_non_empty(group$MatrixKey[
      group$MatrixEligible %in% TRUE
    ]))
    group_scores = suppressWarnings(as.numeric(group$ConfidenceScore))
    group_scores = group_scores[!is.na(group_scores)]
    data.frame(
      Query = group$Query[[1]],
      CID = suppressWarnings(as.integer(group$CID[[1]])),
      SummaryLevel = by,
      TraitType = if (by == "type") group$TraitType[[1]] else NA_character_,
      SourceDatabase = if (by == "source") {
        group$SourceDatabase[[1]]
      } else {
        NA_character_
      },
      TraitCount = length(unique(trait_id)),
      MatrixEligibleCount = length(matrix_keys),
      HighConfidenceCount = sum(group_scores >= 0.90),
      MeanConfidence = if (length(group_scores) > 0) {
        round(mean(group_scores), 3)
      } else {
        NA_real_
      },
      SourceDatabases = .pubchem_collapse(group$SourceDatabase),
      TopTraitGroups = .normalized_top_values(group$TraitGroup, top_n),
      TopTraitValues = .normalized_top_values(group$TraitValue, top_n),
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out[order(out$Query, out$SummaryLevel, out$TraitType, out$SourceDatabase),
      cols,
      drop = FALSE]
}

#' Compare chemicals by normalized traits
#'
#' @description
#' `chemicalTraitSimilarity()` computes pairwise similarity from
#' `ChemicalTraits`, a categorate result, or an already-built trait matrix. It
#' reports shared and distinct trait counts plus collapsed shared/distinct
#' matrix keys so researchers can quickly see why compounds group together.
#'
#' @param traits A `ChemicalTraits` data frame, a categorate result containing
#' `ChemicalTraits`, or a trait matrix with `Query`/`CID` columns.
#' @param profile Trait matrix profile used when `traits` is not already a
#' matrix. See `chemicalTraitMatrix()`.
#' @param min_confidence Minimum confidence score for included traits. Accepts a
#' numeric score or `"low"`, `"medium"`, or `"high"`.
#' @param max_traits Maximum number of trait columns to include before computing
#' similarity.
#' @param top_n Number of shared/distinct matrix keys to include in collapsed
#' explanatory columns.
#'
#' @return A pairwise data frame with Jaccard and overlap similarities.
#'
#' @examples
#' \dontrun{
#' result = categorate(compounds, library_data, detail = "full")
#' chemicalTraitSimilarity(result)
#' chemicalTraitSimilarity(result, profile = "bioactivity")
#' }
#'
#' @export
chemicalTraitSimilarity = function(traits,
                                   profile = c("core", "full", "bioactivity",
                                               "safety", "ecology", "kegg",
                                               "sensory", "biomedical"),
                                   min_confidence = 0,
                                   max_traits = Inf,
                                   top_n = 12) {
  profile = match.arg(profile)
  matrix = .normalized_similarity_matrix(traits = traits,
                                         profile = profile,
                                         min_confidence = min_confidence,
                                         max_traits = max_traits)
  cols = c("QueryA", "CIDA", "QueryB", "CIDB", "SharedTraitCount",
           "UnionTraitCount", "OnlyA_TraitCount", "OnlyB_TraitCount",
           "JaccardSimilarity", "OverlapCoefficient", "SharedTraits",
           "OnlyA_Traits", "OnlyB_Traits")
  if (!is.data.frame(matrix) || nrow(matrix) < 2 || ncol(matrix) <= 2) {
    return(.uaf_empty_table(cols))
  }

  trait_cols = setdiff(colnames(matrix), c("Query", "CID"))
  rows = list()
  for (i in seq_len(nrow(matrix) - 1L)) {
    a = suppressWarnings(as.numeric(matrix[i, trait_cols, drop = TRUE]))
    a[is.na(a)] = 0
    for (j in seq.int(i + 1L, nrow(matrix))) {
      b = suppressWarnings(as.numeric(matrix[j, trait_cols, drop = TRUE]))
      b[is.na(b)] = 0
      a_present = a > 0
      b_present = b > 0
      shared = a_present & b_present
      union = a_present | b_present
      only_a = a_present & !b_present
      only_b = b_present & !a_present
      shared_count = sum(shared)
      union_count = sum(union)
      min_count = min(sum(a_present), sum(b_present))
      rows[[length(rows) + 1]] = data.frame(
        QueryA = matrix$Query[[i]],
        CIDA = suppressWarnings(as.integer(matrix$CID[[i]])),
        QueryB = matrix$Query[[j]],
        CIDB = suppressWarnings(as.integer(matrix$CID[[j]])),
        SharedTraitCount = shared_count,
        UnionTraitCount = union_count,
        OnlyA_TraitCount = sum(only_a),
        OnlyB_TraitCount = sum(only_b),
        JaccardSimilarity = ifelse(union_count > 0,
                                   round(shared_count / union_count, 4),
                                   NA_real_),
        OverlapCoefficient = ifelse(min_count > 0,
                                    round(shared_count / min_count, 4),
                                    NA_real_),
        SharedTraits = .pubchem_collapse(utils::head(trait_cols[shared],
                                                     top_n)),
        OnlyA_Traits = .pubchem_collapse(utils::head(trait_cols[only_a],
                                                    top_n)),
        OnlyB_Traits = .pubchem_collapse(utils::head(trait_cols[only_b],
                                                    top_n)),
        stringsAsFactors = FALSE
      )
    }
  }
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out[order(-out$JaccardSimilarity, -out$SharedTraitCount, out$QueryA,
            out$QueryB),
      cols,
      drop = FALSE]
}

.normalized_extract_traits = function(traits) {
  if (is.list(traits) &&
      !is.data.frame(traits) &&
      is.data.frame(traits$ChemicalTraits)) {
    return(traits$ChemicalTraits)
  }
  if (is.data.frame(traits)) return(traits)
  .uaf_empty_table(.normalized_trait_cols())
}

.normalized_prepare_traits = function(traits) {
  cols = .normalized_trait_cols()
  if (!is.data.frame(traits) || nrow(traits) < 1) return(.uaf_empty_table(cols))
  original_cols = colnames(traits)
  for (col in setdiff(cols, colnames(traits))) traits[[col]] = NA

  traits$TraitType = vapply(traits$TraitType,
                            .normalized_trait_key_piece,
                            character(1))
  traits$TraitGroup = vapply(traits$TraitGroup,
                             .normalized_trait_key_piece,
                             character(1))
  traits$TraitValue = .uaf_squish_text(traits$TraitValue)
  missing_value = is.na(traits$TraitValue) | traits$TraitValue == ""
  traits$TraitValue[missing_value] =
    traits$TraitValueClean[missing_value]
  traits$TraitValueClean = vapply(traits$TraitValue,
                                  .normalized_trait_value_clean,
                                  character(1))
  traits$TraitLabel = mapply(.normalized_trait_label,
                             traits$TraitType,
                             traits$TraitGroup,
                             traits$TraitValue,
                             USE.NAMES = FALSE)
  traits$SourceDatabase = vapply(seq_len(nrow(traits)), function(i) {
    .uaf_first_non_empty_text(
      traits$SourceDatabase[[i]],
      .normalized_trait_source_database(traits$SourceTable[[i]])
    )
  }, character(1))
  traits$Confidence = vapply(traits$Confidence, function(x) {
    .uaf_first_non_empty_text(x, "medium")
  }, character(1))

  scores = suppressWarnings(as.numeric(traits$ConfidenceScore))
  missing_score = is.na(scores)
  if (any(missing_score)) {
    scores[missing_score] = mapply(
      .normalized_trait_confidence_score,
      confidence = traits$Confidence[missing_score],
      source_database = traits$SourceDatabase[missing_score],
      evidence_id = traits$EvidenceID[missing_score],
      extraction_rule = traits$ExtractionRule[missing_score],
      USE.NAMES = FALSE
    )
  }
  traits$ConfidenceScore = scores

  if ("MatrixEligible" %in% original_cols) {
    eligible = traits$MatrixEligible %in% TRUE |
      tolower(as.character(traits$MatrixEligible)) == "true"
    eligible[is.na(eligible)] = FALSE
  } else {
    eligible = rep(TRUE, nrow(traits))
  }
  descriptive = vapply(traits$TraitValue,
                       .normalized_trait_is_descriptive_value,
                       logical(1))
  eligible[descriptive] = FALSE
  traits$MatrixEligible = eligible

  missing_key = is.na(traits$MatrixKey) | traits$MatrixKey == ""
  if (any(missing_key & traits$MatrixEligible)) {
    key_rows = missing_key & traits$MatrixEligible
    traits$MatrixKey[key_rows] = mapply(
      .normalized_trait_matrix_key,
      traits$TraitType[key_rows],
      traits$TraitGroup[key_rows],
      traits$TraitValueClean[key_rows],
      USE.NAMES = FALSE
    )
  }
  traits$MatrixKey[!traits$MatrixEligible] = NA_character_
  traits[, cols, drop = FALSE]
}

.normalized_similarity_matrix = function(traits, profile, min_confidence,
                                         max_traits) {
  if (is.data.frame(traits) && !"TraitType" %in% colnames(traits)) {
    return(traits)
  }
  if (is.list(traits) &&
      !is.data.frame(traits) &&
      is.data.frame(traits$ChemicalTraitMatrix) &&
      !is.data.frame(traits$ChemicalTraits)) {
    return(traits$ChemicalTraitMatrix)
  }
  chemicalTraitMatrix(traits,
                      profile = profile,
                      mode = "binary",
                      min_confidence = min_confidence,
                      max_traits = max_traits)
}

.normalized_top_values = function(values, top_n = 8) {
  values = .uaf_non_empty(values)
  top_n = suppressWarnings(as.integer(top_n[[1]]))
  if (length(values) < 1 || is.na(top_n) || top_n < 1) return(NA_character_)
  tab = sort(table(values), decreasing = TRUE)
  .pubchem_collapse(utils::head(names(tab), top_n))
}

.normalized_chemical_trait_matrix = function(chemical_traits,
                                             profile = "core",
                                             mode = "binary",
                                             min_confidence = 0,
                                             max_traits = Inf) {
  base_cols = c("Query", "CID")
  if (!is.data.frame(chemical_traits) || nrow(chemical_traits) < 1) {
    return(.uaf_empty_table(base_cols))
  }
  chemical_traits = .normalized_prepare_traits(chemical_traits)
  base = unique(chemical_traits[, base_cols, drop = FALSE])
  row.names(base) = NULL
  threshold = .normalized_trait_confidence_threshold(min_confidence)
  source = .normalized_filter_trait_matrix_source(
    chemical_traits = chemical_traits,
    profile = profile,
    min_confidence = threshold
  )
  if (nrow(source) < 1) return(base)

  keys = .normalized_select_trait_matrix_keys(source, max_traits)
  source = source[source$MatrixKey %in% keys, , drop = FALSE]
  zero_value = if (mode == "confidence") 0 else 0L
  for (key in keys) base[[key]] = zero_value
  row_key = paste(.pubchem_key_value(base$Query),
                  .pubchem_key_value(base$CID),
                  sep = "\r")
  source_key = paste(.pubchem_key_value(source$Query),
                     .pubchem_key_value(source$CID),
                     sep = "\r")
  for (key in keys) {
    key_source = source[source$MatrixKey == key, , drop = FALSE]
    key_source_key = source_key[source$MatrixKey == key]
    for (compound_key in unique(key_source_key)) {
      hit_row = match(compound_key, row_key)
      if (is.na(hit_row)) next
      compound_rows = key_source[key_source_key == compound_key, ,
                                 drop = FALSE]
      base[[key]][[hit_row]] = .normalized_trait_matrix_value(compound_rows,
                                                              mode)
    }
  }
  base
}

.normalized_filter_trait_matrix_source = function(chemical_traits, profile,
                                                 min_confidence) {
  source = chemical_traits[chemical_traits$MatrixEligible %in% TRUE &
                             !is.na(chemical_traits$MatrixKey) &
                             chemical_traits$MatrixKey != "",
                           ,
                           drop = FALSE]
  if (nrow(source) < 1) return(source)
  scores = suppressWarnings(as.numeric(source$ConfidenceScore))
  scores[is.na(scores)] = 0
  source = source[scores >= min_confidence, , drop = FALSE]
  if (nrow(source) < 1) return(source)
  keep = .normalized_trait_profile_keep(source, profile)
  source[keep, , drop = FALSE]
}

.normalized_trait_confidence_threshold = function(min_confidence) {
  if (is.character(min_confidence) && length(min_confidence) > 0) {
    value = tolower(.uaf_first_non_empty_text(min_confidence, "0"))
    if (value == "high") return(0.90)
    if (value == "medium") return(0.65)
    if (value == "low") return(0.35)
  }
  value = suppressWarnings(as.numeric(min_confidence[[1]]))
  if (is.na(value)) return(0)
  max(0, min(value, 1))
}

.normalized_trait_profile_keep = function(traits, profile) {
  if (profile == "full") return(rep(TRUE, nrow(traits)))
  key = paste(traits$TraitType, traits$TraitGroup, sep = ":")
  if (profile == "core") {
    return(key %in% c(
      "property:molecular_size_bin",
      "property:xlogp_bin",
      "property:tpsa_bin",
      "property:rotatable_bond_bin",
      "property:hbond_donor_bin",
      "property:hbond_acceptor_bin",
      "hazard:hazard_code",
      "hazard:hazard_category",
      "hazard:hazard_group",
      "hazard:signal_word",
      "hazard:exposure_route",
      "hazard:target_organ",
      "hazard:toxicity_metric",
      "hazard:toxicity_species",
      "use:descriptor",
      "use:flavor",
      "use:odor",
      "use:regulatory_status",
      "use:pharmacologic_class",
      "use:pharmacologic_action",
      "use:administration_route",
      "use:dosage_form",
      "use:warning",
      "use:use_group",
      "chemical_class:reactive_group",
      "chemical_class:natural_product_superclass",
      "chemical_class:natural_product_class",
      "chemical_class:natural_product_subclass",
      "chemical_class:pathway_group",
      "chemical_class:enzyme_class",
      "chemical_class:mesh_tree_category",
      "taxonomy:kingdom",
      "taxonomy:family",
      "taxonomy:taxonomy_group",
      "taxonomy:natural_product_class",
      "occurrence:occurrence_type",
      "occurrence:kingdom",
      "occurrence:family",
      "occurrence:natural_product_class",
      "pathway:pathway_group",
      "pathway:enzyme_class",
      "pathway:role_type",
      "reaction:participant_role",
      "reaction:participant_side",
      "bioactivity:activity_class",
      "bioactivity:activity_direction",
      "bioactivity:bioactivity_domain",
      "bioactivity:endpoint_metric",
      "target:target_type",
      "target:target_organism",
      "target:target_group",
      "target:endpoint_metric",
      "potency:potency_metric",
      "potency:potency_bucket"
    ))
  }
  if (profile == "bioactivity") {
    return(traits$TraitType %in% c("bioactivity", "target", "potency"))
  }
  if (profile == "safety") {
    return(traits$TraitType == "hazard" |
             key %in% c("use:warning", "use:administration_route",
                        "property:xlogp_bin", "property:tpsa_bin"))
  }
  if (profile == "ecology") {
    return(traits$TraitType %in% c("taxonomy", "occurrence") |
             key %in% c("chemical_class:natural_product_superclass",
                        "chemical_class:natural_product_class",
                        "chemical_class:natural_product_subclass"))
  }
  if (profile == "kegg") {
    return(traits$TraitType %in% c("pathway", "reaction") |
             key %in% c("chemical_class:pathway_group",
                        "chemical_class:enzyme_class"))
  }
  if (profile == "sensory") {
    return(key %in% c("use:descriptor", "use:flavor", "use:odor",
                      "use:regulatory_status", "use:use_group"))
  }
  if (profile == "biomedical") {
    return(traits$TraitType %in% c("bioactivity", "target", "potency") |
             key %in% c("use:active_ingredient",
                        "use:pharmacologic_class",
                        "use:pharmacologic_action",
                        "use:administration_route",
                        "use:dosage_form",
                        "use:use_group",
                        "chemical_class:mesh_tree_category",
                        "property:molecular_size_bin",
                        "property:xlogp_bin",
                        "property:tpsa_bin"))
  }
  rep(FALSE, nrow(traits))
}

.normalized_select_trait_matrix_keys = function(source, max_traits) {
  keys = unique(source$MatrixKey)
  max_value = suppressWarnings(as.numeric(max_traits[[1]]))
  if (is.na(max_value) || max_value < 1) return(character())
  max_traits = if (is.infinite(max_value)) {
    length(keys)
  } else {
    as.integer(max_value)
  }
  source_key = paste(.pubchem_key_value(source$Query),
                     .pubchem_key_value(source$CID),
                     sep = "\r")
  prevalence = vapply(keys, function(key) {
    length(unique(source_key[source$MatrixKey == key]))
  }, integer(1))
  confidence = vapply(keys, function(key) {
    values = suppressWarnings(as.numeric(source$ConfidenceScore[
      source$MatrixKey == key
    ]))
    values = values[!is.na(values)]
    if (length(values) < 1) return(0)
    mean(values)
  }, numeric(1))
  ranked = data.frame(MatrixKey = keys,
                      Prevalence = prevalence,
                      Confidence = confidence,
                      stringsAsFactors = FALSE)
  ranked = ranked[order(-ranked$Prevalence, -ranked$Confidence,
                        ranked$MatrixKey),
                  ,
                  drop = FALSE]
  utils::head(ranked$MatrixKey, max_traits)
}

.normalized_trait_matrix_value = function(rows, mode) {
  if (mode == "count") return(nrow(rows))
  if (mode == "confidence") {
    values = suppressWarnings(as.numeric(rows$ConfidenceScore))
    values = values[!is.na(values)]
    if (length(values) < 1) return(0)
    return(max(values))
  }
  1L
}

.normalized_chemical_classes = function(compounds, data_list, pubchem_profiles,
                                        kegg, cid_lookup) {
  cols = c("Query", "CID", "ClassSystem", "ClassType", "ClassID",
           "ClassName", "ClassGroup", "SourceTable", "EvidenceText",
           "EvidenceURL", "ExtractionRule", "Confidence")
  rows = list()

  reactives = data_list$reactives
  if (is.data.frame(reactives) &&
      "Chemical" %in% colnames(reactives) &&
      "reactives" %in% colnames(reactives)) {
    for (i in seq_len(nrow(reactives))) {
      value = .uaf_first_non_empty_text(reactives$reactives[[i]])
      query = .uaf_first_non_empty_text(reactives$Chemical[[i]])
      if (is.na(value) || value %in% c("None", "NA")) next
      rows[[length(rows) + 1]] = data.frame(
        Query = query,
        CID = .normalized_lookup_cid(query, cid_lookup),
        ClassSystem = "PubChem/reactives",
        ClassType = "reactive_group",
        ClassID = NA_character_,
        ClassName = value,
        ClassGroup = .normalized_clean_label(value),
        SourceTable = "reactives",
        EvidenceText = value,
        EvidenceURL = NA_character_,
        ExtractionRule = "legacy_reactive_group",
        Confidence = "medium",
        stringsAsFactors = FALSE
      )
    }
  }

  lotus = pubchem_profiles$LOTUSProfile
  if (is.data.frame(lotus) && nrow(lotus) > 0) {
    for (i in seq_len(nrow(lotus))) {
      row = lotus[i, , drop = FALSE]
      classes = .normalized_split_terms(row$NaturalProductClass)
      for (class_name in classes) {
        lookup = .normalized_term_lookup(class_name)
        rows[[length(rows) + 1]] = data.frame(
          Query = row$Query,
          CID = row$CID,
          ClassSystem = "LOTUS",
          ClassType = "natural_product_class",
          ClassID = row$LOTUS_ID,
          ClassName = class_name,
          ClassGroup = lookup$TermGroup,
          SourceTable = "LOTUSProfile",
          EvidenceText = row$RawValue,
          EvidenceURL = row$SourceURL,
          ExtractionRule = "natural_product_class_terms",
          Confidence = lookup$Confidence,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  classifications = pubchem_profiles$PubChemClassificationProfile
  if (is.data.frame(classifications) && nrow(classifications) > 0) {
    chemical_classes = classifications[
      classifications$TreeType == "chemical" &
        grepl("LOTUS", paste(classifications$Source,
                             classifications$TreeName,
                             classifications$ClassPath),
              ignore.case = TRUE),
      ,
      drop = FALSE
    ]
    for (i in seq_len(nrow(chemical_classes))) {
      row = chemical_classes[i, , drop = FALSE]
      path_terms = .normalized_classification_path_terms(row$ClassPath,
                                                         row$ClassName)
      if (length(path_terms) < 1) next
      class_group = .normalized_clean_label(path_terms[[1]])
      for (level in seq_along(path_terms)) {
        class_name = path_terms[[level]]
        class_type = if (level == 1L) {
          "natural_product_superclass"
        } else if (level == length(path_terms)) {
          if (length(path_terms) == 1L) {
            "natural_product_class"
          } else {
            "natural_product_subclass"
          }
        } else {
          "natural_product_class"
        }
        rows[[length(rows) + 1]] = data.frame(
          Query = row$Query,
          CID = row$CID,
          ClassSystem = "LOTUS",
          ClassType = class_type,
          ClassID = if (level == length(path_terms)) row$HNID else NA_character_,
          ClassName = class_name,
          ClassGroup = class_group,
          SourceTable = "PubChemClassificationProfile",
          EvidenceText = row$ClassPath,
          EvidenceURL = row$ClassificationURL,
          ExtractionRule = "lotus_chemical_tree",
          Confidence = "high",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  mesh = pubchem_profiles$MeSHProfile
  if (is.data.frame(mesh) && nrow(mesh) > 0) {
    for (i in seq_len(nrow(mesh))) {
      row = mesh[i, , drop = FALSE]
      tree_numbers = .normalized_split_codes(row$TreeCategory)
      if (length(tree_numbers) < 1) tree_numbers = NA_character_
      for (tree_number in tree_numbers) {
        rows[[length(rows) + 1]] = data.frame(
          Query = row$Query,
          CID = row$CID,
          ClassSystem = "MeSH",
          ClassType = "mesh_tree_category",
          ClassID = tree_number,
          ClassName = .normalized_mesh_label(row$PharmacologicAction,
                                             row$Descriptor),
          ClassGroup = .normalized_mesh_tree_group(tree_number),
          SourceTable = "MeSHProfile",
          EvidenceText = row$RawValue,
          EvidenceURL = row$SourceURL,
          ExtractionRule = "mesh_tree_number",
          Confidence = ifelse(is.na(tree_number), "medium", "high"),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (is.data.frame(kegg$pathways) && nrow(kegg$pathways) > 0) {
    for (i in seq_len(nrow(kegg$pathways))) {
      row = kegg$pathways[i, , drop = FALSE]
      group = .uaf_first_non_empty_text(row$PathwayGroup)
      if (is.na(group)) next
      rows[[length(rows) + 1]] = data.frame(
        Query = row$Query,
        CID = .normalized_lookup_cid(row$Query, cid_lookup),
        ClassSystem = "KEGG",
        ClassType = "pathway_group",
        ClassID = row$PathwayID,
        ClassName = row$PathwayName,
        ClassGroup = group,
        SourceTable = "KEGGPathways",
        EvidenceText = row$Evidence,
        EvidenceURL = row$EvidenceURL,
        ExtractionRule = "kegg_pathway_group",
        Confidence = "high",
        stringsAsFactors = FALSE
      )
    }
  }

  if (is.data.frame(kegg$enzymes) && nrow(kegg$enzymes) > 0) {
    for (i in seq_len(nrow(kegg$enzymes))) {
      row = kegg$enzymes[i, , drop = FALSE]
      group = .uaf_first_non_empty_text(row$EnzymeClass)
      if (is.na(group)) next
      rows[[length(rows) + 1]] = data.frame(
        Query = row$Query,
        CID = .normalized_lookup_cid(row$Query, cid_lookup),
        ClassSystem = "KEGG",
        ClassType = "enzyme_class",
        ClassID = row$ECNumber,
        ClassName = row$EnzymeName,
        ClassGroup = group,
        SourceTable = "KEGGEnzymes",
        EvidenceText = row$Evidence,
        EvidenceURL = row$EvidenceURL,
        ExtractionRule = "kegg_enzyme_class",
        Confidence = "high",
        stringsAsFactors = FALSE
      )
    }
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "CID", "ClassSystem", "ClassType", "ClassID",
                 "ClassName", "ClassGroup"),
    collapse_cols = c("SourceTable", "EvidenceText", "EvidenceURL",
                      "ExtractionRule")
  )
}

.normalized_chemical_measurements = function(pubchem, pubchem_profiles) {
  cols = .normalized_measurement_cols()
  rows = list()

  if (is.data.frame(pubchem$properties) && nrow(pubchem$properties) > 0) {
    property_units = c(MolecularWeight = "Da", ExactMass = "Da",
                       MonoisotopicMass = "Da", XLogP = NA_character_,
                       TPSA = "A^2", Complexity = NA_character_,
                       Charge = "count", HBondDonorCount = "count",
                       HBondAcceptorCount = "count",
                       RotatableBondCount = "count",
                       HeavyAtomCount = "count", IsotopeAtomCount = "count",
                       AtomStereoCount = "count",
                       DefinedAtomStereoCount = "count",
                       UndefinedAtomStereoCount = "count",
                       BondStereoCount = "count",
                       DefinedBondStereoCount = "count",
                       UndefinedBondStereoCount = "count",
                       CovalentUnitCount = "count")
    for (i in seq_len(nrow(pubchem$properties))) {
      row = pubchem$properties[i, , drop = FALSE]
      for (property in intersect(names(property_units), colnames(row))) {
        value = suppressWarnings(as.numeric(row[[property]][[1]]))
        if (is.na(value)) next
        rows[[length(rows) + 1]] = data.frame(
          Query = row$Query,
          CID = suppressWarnings(as.integer(row$CID)),
          Property = .normalized_property_name(property),
          Value = value,
          Unit = unname(property_units[[property]]),
          Condition = NA_character_,
          Method = "PubChem property",
          SourceTable = "PubChemProperties",
          Source = "PubChem",
          EvidenceText = property,
          EvidenceURL = row$SourceURL,
          PubChemURL = NA_character_,
          ExtractionRule = "pubchem_property_numeric",
          Confidence = "high",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  annotations = pubchem$annotations
  if (is.data.frame(annotations) && nrow(annotations) > 0) {
    keep = !is.na(annotations$ValueNumeric) &
      grepl("Boiling|Melting|Flash|Density|Solubility|Vapor|Pressure|pKa|LogP|Henry|Mass|Point|Temperature",
            paste(annotations$HeadingPath, annotations$Name, annotations$CleanValue),
            ignore.case = TRUE)
    measurement_rows = annotations[keep, , drop = FALSE]
    for (i in seq_len(nrow(measurement_rows))) {
      row = measurement_rows[i, , drop = FALSE]
      property = .normalized_property_name(.uaf_first_non_empty_text(row$Name,
                                                                     row$HeadingPath))
      if (is.na(property)) next
      rows[[length(rows) + 1]] = data.frame(
        Query = row$Query,
        CID = suppressWarnings(as.integer(row$CID)),
        Property = property,
        Value = row$ValueNumeric,
        Unit = row$UnitClean,
        Condition = NA_character_,
        Method = "PubChem annotation",
        SourceTable = "PubChemAnnotations",
        Source = row$Source,
        EvidenceText = row$CleanValue,
        EvidenceURL = row$SourceURL,
        PubChemURL = row$PubChemURL,
        ExtractionRule = "annotation_measurement",
        Confidence = "medium",
        stringsAsFactors = FALSE
      )
    }
  }

  safety = pubchem_profiles$SafetyProfile
  if (is.data.frame(safety) && nrow(safety) > 0) {
    for (i in seq_len(nrow(safety))) {
      row = safety[i, , drop = FALSE]
      tox = .normalized_toxicity_measurements(.pubchem_collapse(c(row$RawValue,
                                                                  row$HazardStatement)))
      if (nrow(tox) < 1) next
      for (j in seq_len(nrow(tox))) {
        rows[[length(rows) + 1]] = data.frame(
          Query = row$Query,
          CID = suppressWarnings(as.integer(row$CID)),
          Property = tolower(tox$Metric[[j]]),
          Value = tox$Value[[j]],
          Unit = tox$Unit[[j]],
          Condition = .pubchem_collapse(c(tox$Route[[j]], tox$Species[[j]])),
          Method = "toxicity extraction",
          SourceTable = "SafetyProfile",
          Source = row$Source,
          EvidenceText = row$RawValue,
          EvidenceURL = row$SourceURL,
          PubChemURL = row$PubChemURL,
          ExtractionRule = "toxicity_measurement_regex",
          Confidence = "medium",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out = .normalized_bind_rows(rows, cols)
  out = .normalized_enrich_measurements(out)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "CID", "Property", "Value", "Unit", "Condition"),
    collapse_cols = c("Method", "SourceTable", "Source", "EvidenceText",
                      "EvidenceURL", "PubChemURL", "ExtractionRule")
  )
}

.normalized_measurement_cols = function() {
  c("Query", "CID", "Property", "Value", "Unit", "ValueRelation",
    "StandardValue", "StandardUnit", "StandardValueLow",
    "StandardValueHigh", "MeasurementClass", "BehaviorBin", "Condition",
    "Method", "SourceTable", "Source", "EvidenceText", "EvidenceURL",
    "PubChemURL", "ExtractionRule", "Confidence")
}

.normalized_enrich_measurements = function(measurements) {
  cols = .normalized_measurement_cols()
  if (!is.data.frame(measurements) || nrow(measurements) < 1) {
    return(.uaf_empty_table(cols))
  }
  for (col in setdiff(cols, colnames(measurements))) measurements[[col]] = NA
  enriched = lapply(seq_len(nrow(measurements)), function(i) {
    row = measurements[i, , drop = FALSE]
    standard = .normalized_standard_measurement(
      property = row$Property[[1]],
      value = row$Value[[1]],
      unit = row$Unit[[1]],
      evidence_text = row$EvidenceText[[1]]
    )
    row$ValueRelation = standard$ValueRelation
    row$StandardValue = standard$StandardValue
    row$StandardUnit = standard$StandardUnit
    row$StandardValueLow = standard$StandardValueLow
    row$StandardValueHigh = standard$StandardValueHigh
    row$MeasurementClass = standard$MeasurementClass
    row$BehaviorBin = standard$BehaviorBin
    row
  })
  out = do.call(rbind, enriched)
  row.names(out) = NULL
  out[, cols, drop = FALSE]
}

.normalized_standard_measurement = function(property, value, unit,
                                            evidence_text = NA_character_) {
  property = .normalized_property_name(property)
  value = suppressWarnings(as.numeric(value))
  unit_clean = .normalized_measurement_unit(unit)
  relation = .normalized_measurement_relation(evidence_text)
  range = .normalized_measurement_range(evidence_text)
  class = .normalized_measurement_class(property)
  standard = .normalized_measurement_standard_value(property, value,
                                                    unit_clean)
  standard_value = standard$value
  standard_unit = standard$unit
  if (!is.na(range$low) && !is.na(range$high)) {
    low_standard = .normalized_measurement_standard_value(property,
                                                          range$low,
                                                          unit_clean)
    high_standard = .normalized_measurement_standard_value(property,
                                                           range$high,
                                                           unit_clean)
    standard_low = low_standard$value
    standard_high = high_standard$value
  } else if (relation %in% c("less_than", "less_or_equal")) {
    standard_low = NA_real_
    standard_high = standard_value
  } else if (relation %in% c("greater_than", "greater_or_equal")) {
    standard_low = standard_value
    standard_high = NA_real_
  } else {
    standard_low = standard_value
    standard_high = standard_value
  }
  if (!is.na(range$low) && !is.na(range$high)) relation = "range"
  list(
    ValueRelation = relation,
    StandardValue = standard_value,
    StandardUnit = standard_unit,
    StandardValueLow = standard_low,
    StandardValueHigh = standard_high,
    MeasurementClass = class,
    BehaviorBin = .normalized_measurement_behavior_bin(property,
                                                       standard_value,
                                                       standard_unit)
  )
}

.normalized_measurement_unit = function(unit) {
  unit = .uaf_first_non_empty_text(unit)
  if (is.na(unit)) return(NA_character_)
  unit = gsub("\u00b0", "", unit, fixed = TRUE)
  unit = gsub("\\s+", " ", unit)
  unit = trimws(unit)
  lower = tolower(unit)
  if (lower %in% c("c", "deg c", "degree c", "degrees c", "celsius")) return("C")
  if (lower %in% c("f", "deg f", "fahrenheit")) return("F")
  if (lower %in% c("k", "kelvin")) return("K")
  if (lower %in% c("da", "g/mol", "amu")) return("Da")
  if (lower %in% c("a^2", "angstrom^2", "angstrom2", "\u00c5^2")) return("A^2")
  if (lower %in% c("g/ml", "g/mL", "g cm-3", "g/cm3", "g/cc")) return("g/mL")
  if (lower %in% c("kg/m3", "kg m-3")) return("kg/m3")
  if (lower %in% c("mg/l", "mg/liter", "mg/litre")) return("mg/L")
  if (lower %in% c("g/l", "g/liter", "g/litre")) return("g/L")
  if (lower %in% c("mg/ml")) return("mg/mL")
  if (lower %in% c("ug/ml", "\u00b5g/ml", "mcg/ml")) return("ug/mL")
  if (lower %in% c("mg/kg", "mg kg-1")) return("mg/kg")
  if (lower %in% c("ug/kg", "\u00b5g/kg", "mcg/kg")) return("ug/kg")
  if (lower %in% c("g/kg")) return("g/kg")
  if (lower %in% c("pa")) return("Pa")
  if (lower %in% c("kpa")) return("kPa")
  if (lower %in% c("mmhg", "torr")) return("mmHg")
  if (lower %in% c("atm")) return("atm")
  if (lower %in% c("count")) return("count")
  unit
}

.normalized_measurement_relation = function(evidence_text) {
  text = tolower(.pubchem_collapse(evidence_text))
  if (is.na(text) || text == "") return("equal")
  range = .normalized_measurement_range(text)
  if (!is.na(range$low) && !is.na(range$high)) {
    return("range")
  }
  if (grepl("<=|\u2264", text, perl = TRUE)) return("less_or_equal")
  if (grepl(">=|\u2265", text, perl = TRUE)) return("greater_or_equal")
  if (grepl("<", text, fixed = TRUE)) return("less_than")
  if (grepl(">", text, fixed = TRUE)) return("greater_than")
  if (grepl("~|approx|approximately|about|ca\\.", text, perl = TRUE)) {
    return("approximate")
  }
  "equal"
}

.normalized_measurement_range = function(evidence_text) {
  empty = list(low = NA_real_, high = NA_real_)
  text = .uaf_first_non_empty_text(evidence_text)
  if (is.na(text)) return(empty)
  pattern = "([-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?)\\s*(?:-|to)\\s*([-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?)"
  match = regexec(pattern, text, perl = TRUE, ignore.case = TRUE)
  pieces = regmatches(text, match)[[1]]
  if (length(pieces) < 3) return(empty)
  low = suppressWarnings(as.numeric(pieces[[2]]))
  high = suppressWarnings(as.numeric(pieces[[3]]))
  if (is.na(low) || is.na(high)) return(empty)
  list(low = min(low, high), high = max(low, high))
}

.normalized_measurement_class = function(property) {
  property = .normalized_property_name(property)
  if (property %in% c("molecular_weight", "exact_mass",
                      "monoisotopic_mass")) return("mass")
  if (property %in% c("boiling_point", "melting_point",
                      "flash_point")) return("temperature")
  if (property %in% c("density")) return("density")
  if (property %in% c("solubility")) return("solubility")
  if (property %in% c("vapor_pressure")) return("vapor_pressure")
  if (property %in% c("logp", "pka", "tpsa", "complexity",
                      "charge")) return("descriptor")
  if (property %in% c("ld50", "lc50", "ec50", "ic50")) return("toxicity")
  if (grepl("count$", property)) return("count")
  "measurement"
}

.normalized_measurement_standard_value = function(property, value, unit) {
  if (is.na(value)) return(list(value = NA_real_, unit = NA_character_))
  class = .normalized_measurement_class(property)
  if (class == "temperature") {
    if (identical(unit, "F")) return(list(value = round((value - 32) * 5 / 9, 4),
                                         unit = "C"))
    if (identical(unit, "K")) return(list(value = round(value - 273.15, 4),
                                         unit = "C"))
    return(list(value = value, unit = "C"))
  }
  if (class == "density") {
    if (identical(unit, "kg/m3")) return(list(value = value / 1000,
                                             unit = "g/mL"))
    return(list(value = value, unit = "g/mL"))
  }
  if (class == "solubility") {
    if (identical(unit, "g/L")) return(list(value = value * 1000,
                                           unit = "mg/L"))
    if (identical(unit, "mg/mL")) return(list(value = value * 1000,
                                             unit = "mg/L"))
    if (identical(unit, "ug/mL")) return(list(value = value,
                                             unit = "mg/L"))
    return(list(value = value, unit = "mg/L"))
  }
  if (class == "vapor_pressure") {
    if (identical(unit, "kPa")) return(list(value = value * 1000,
                                           unit = "Pa"))
    if (identical(unit, "mmHg")) return(list(value = value * 133.322,
                                            unit = "Pa"))
    if (identical(unit, "atm")) return(list(value = value * 101325,
                                           unit = "Pa"))
    return(list(value = value, unit = "Pa"))
  }
  if (class == "toxicity") {
    if (identical(unit, "ug/kg")) return(list(value = value / 1000,
                                             unit = "mg/kg"))
    if (identical(unit, "g/kg")) return(list(value = value * 1000,
                                            unit = "mg/kg"))
    return(list(value = value, unit = "mg/kg"))
  }
  if (property %in% c("molecular_weight", "exact_mass",
                      "monoisotopic_mass")) {
    return(list(value = value, unit = "Da"))
  }
  if (property == "tpsa") return(list(value = value, unit = "A^2"))
  if (grepl("count$|charge", property)) return(list(value = value,
                                                    unit = "count"))
  list(value = value, unit = unit)
}

.normalized_measurement_behavior_bin = function(property, value, unit) {
  property = .normalized_property_name(property)
  if (is.na(value)) return(NA_character_)
  if (property %in% c("molecular_weight", "exact_mass",
                      "monoisotopic_mass")) {
    if (value < 200) return("small_molecule_lt_200_Da")
    if (value <= 500) return("drug_like_size_200_500_Da")
    return("large_molecule_gt_500_Da")
  }
  if (property == "boiling_point") {
    if (value < 100) return("highly_volatile_boiling_point_lt_100_C")
    if (value < 200) return("volatile_boiling_point_100_200_C")
    if (value < 300) return("semi_volatile_boiling_point_200_300_C")
    return("low_volatility_boiling_point_ge_300_C")
  }
  if (property == "melting_point") {
    if (value < 25) return("liquid_or_low_melting_lt_25_C")
    if (value < 100) return("low_melting_25_100_C")
    if (value < 200) return("moderate_melting_100_200_C")
    return("high_melting_ge_200_C")
  }
  if (property == "flash_point") {
    if (value < 23) return("highly_flammable_flash_point_lt_23_C")
    if (value < 60) return("flammable_flash_point_23_60_C")
    if (value < 93) return("combustible_flash_point_60_93_C")
    return("higher_flash_point_ge_93_C")
  }
  if (property == "density") {
    if (value < 0.8) return("low_density_lt_0_8_g_mL")
    if (value <= 1.2) return("water_like_density_0_8_1_2_g_mL")
    return("high_density_gt_1_2_g_mL")
  }
  if (property == "solubility") {
    if (value < 10) return("low_solubility_lt_10_mg_L")
    if (value < 1000) return("moderate_solubility_10_1000_mg_L")
    return("high_solubility_ge_1000_mg_L")
  }
  if (property == "vapor_pressure") {
    if (value < 1) return("low_vapor_pressure_lt_1_Pa")
    if (value < 1000) return("moderate_vapor_pressure_1_1000_Pa")
    return("high_vapor_pressure_ge_1000_Pa")
  }
  if (property == "ld50") {
    if (value <= 50) return("very_high_acute_toxicity_le_50_mg_kg")
    if (value <= 300) return("high_acute_toxicity_50_300_mg_kg")
    if (value <= 2000) return("moderate_acute_toxicity_300_2000_mg_kg")
    return("lower_acute_toxicity_gt_2000_mg_kg")
  }
  if (property == "logp") {
    if (value < 0) return("hydrophilic_logp_lt_0")
    if (value <= 3) return("balanced_logp_0_3")
    if (value <= 5) return("lipophilic_logp_3_5")
    return("highly_lipophilic_logp_gt_5")
  }
  NA_character_
}

#' Summarize standardized chemical measurements
#'
#' @description
#' `chemicalMeasurementSummary()` condenses `ChemicalMeasurements` into one row
#' per compound and property using standardized units where possible. It is
#' designed for filtering and plotting measurement behavior without manually
#' parsing source units or evidence text.
#'
#' @param x A categorate result containing `ChemicalMeasurements`, or a
#' `ChemicalMeasurements` data frame.
#'
#' @return A data frame with one row per compound-property combination.
#'
#' @examples
#' \dontrun{
#' result = categorate(compounds, library_data, detail = "research")
#' chemicalMeasurementSummary(result)
#' }
#'
#' @export
chemicalMeasurementSummary = function(x) {
  measurements = if (is.list(x) &&
                     !is.data.frame(x) &&
                     is.data.frame(x$ChemicalMeasurements)) {
    x$ChemicalMeasurements
  } else if (is.data.frame(x)) {
    x
  } else {
    .uaf_empty_table(.normalized_measurement_cols())
  }
  .normalized_chemical_measurement_summary(measurements)
}

.normalized_chemical_measurement_summary = function(measurements) {
  cols = .normalized_measurement_summary_cols()
  if (!is.data.frame(measurements) || nrow(measurements) < 1) {
    return(.uaf_empty_table(cols))
  }
  for (col in setdiff(.normalized_measurement_cols(), colnames(measurements))) {
    measurements[[col]] = NA
  }
  key = paste(.pubchem_key_value(measurements$Query),
              .pubchem_key_value(measurements$CID),
              .pubchem_key_value(measurements$Property),
              sep = "\r")
  rows = lapply(unique(key), function(group_key) {
    group = measurements[key == group_key, , drop = FALSE]
    values = suppressWarnings(as.numeric(group$StandardValue))
    values = values[!is.na(values)]
    best_row = .normalized_measurement_best_row(group)
    data.frame(
      Query = group$Query[[1]],
      CID = suppressWarnings(as.integer(group$CID[[1]])),
      Property = group$Property[[1]],
      MeasurementClass = .uaf_first_non_empty_text(group$MeasurementClass),
      MeasurementCount = nrow(group),
      StandardUnit = .uaf_first_non_empty_text(group$StandardUnit),
      MinStandardValue = if (length(values) > 0) min(values) else NA_real_,
      MedianStandardValue = if (length(values) > 0) stats::median(values) else NA_real_,
      MaxStandardValue = if (length(values) > 0) max(values) else NA_real_,
      BestValue = suppressWarnings(as.numeric(best_row$StandardValue[[1]])),
      BestValueRelation = best_row$ValueRelation[[1]],
      BehaviorBins = .pubchem_collapse(group$BehaviorBin),
      SourceTables = .pubchem_collapse(group$SourceTable),
      Sources = .pubchem_collapse(group$Source),
      EvidenceURLs = .pubchem_collapse(c(group$EvidenceURL, group$PubChemURL)),
      Confidence = .normalized_best_confidence(group$Confidence),
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out[order(out$Query, out$Property), cols, drop = FALSE]
}

.normalized_measurement_summary_cols = function() {
  c("Query", "CID", "Property", "MeasurementClass", "MeasurementCount",
    "StandardUnit", "MinStandardValue", "MedianStandardValue",
    "MaxStandardValue", "BestValue", "BestValueRelation", "BehaviorBins",
    "SourceTables", "Sources", "EvidenceURLs", "Confidence")
}

.normalized_measurement_best_row = function(group) {
  values = suppressWarnings(as.numeric(group$StandardValue))
  confidence = match(tolower(group$Confidence), c("low", "medium", "high"))
  confidence[is.na(confidence)] = 0
  has_value = !is.na(values)
  rank = order(!has_value, -confidence, group$SourceTable, group$Property)
  group[rank[[1]], , drop = FALSE]
}

.normalized_property_name = function(name) {
  name = .uaf_first_non_empty_text(name)
  if (is.na(name)) return(NA_character_)
  text = tolower(name)
  if (grepl("molecularweight|molecular weight", text)) return("molecular_weight")
  if (grepl("exactmass|exact mass", text)) return("exact_mass")
  if (grepl("monoisotopic", text)) return("monoisotopic_mass")
  if (grepl("xlogp|logp|log kow|logkow", text)) return("logp")
  if (grepl("tpsa|polar surface", text)) return("tpsa")
  if (grepl("complexity", text)) return("complexity")
  if (grepl("charge", text)) return("charge")
  if (grepl("hbonddonor|hydrogen bond donor", text)) return("h_bond_donor_count")
  if (grepl("hbondacceptor|hydrogen bond acceptor", text)) return("h_bond_acceptor_count")
  if (grepl("rotatable", text)) return("rotatable_bond_count")
  if (grepl("heavyatom|heavy atom", text)) return("heavy_atom_count")
  if (grepl("boiling", text)) return("boiling_point")
  if (grepl("melting", text)) return("melting_point")
  if (grepl("flash", text)) return("flash_point")
  if (grepl("vapor|vapour", text)) return("vapor_pressure")
  if (grepl("density", text)) return("density")
  if (grepl("solubility", text)) return("solubility")
  if (grepl("\\bpka\\b", text)) return("pka")
  if (grepl("henry", text)) return("henrys_law_constant")
  .normalized_clean_label(name)
}

.normalized_chemical_hazards = function(pubchem_profiles) {
  cols = c("Query", "CID", "HazardCode", "HazardCategory", "HazardGroup",
           "SignalWord", "HazardClass", "HazardStatement", "PrecautionCode",
           "ExposureRoute", "TargetOrgan", "ToxicityMetric", "ToxicityValue",
           "ToxicityUnit", "Species", "Source", "EvidenceText",
           "EvidenceURL", "PubChemURL", "ExtractionRule", "Confidence")
  safety = pubchem_profiles$SafetyProfile
  if (!is.data.frame(safety) || nrow(safety) < 1) return(.uaf_empty_table(cols))

  rows = list()
  for (i in seq_len(nrow(safety))) {
    row = safety[i, , drop = FALSE]
    evidence = .pubchem_collapse(c(row$RawValue, row$HazardStatement,
                                   row$HazardClass))
    hazard_codes = .normalized_split_terms(row$HazardCode)
    if (length(hazard_codes) < 1) {
      hazard_codes = .uaf_extract_pattern(evidence, "\\bH[0-9]{3}[A-Za-z]?\\b")
    }
    if (length(hazard_codes) < 1) hazard_codes = NA_character_
    precaution_codes = .pubchem_collapse(c(.normalized_split_terms(row$PrecautionCode),
                                           .uaf_extract_pattern(evidence, "\\bP[0-9]{3}[A-Za-z]?\\b")))
    routes = .pubchem_collapse(.pubchem_extract_routes(evidence))
    organs = .pubchem_collapse(.normalized_target_organs(evidence))
    toxicity = .normalized_toxicity_measurements(evidence)
    if (nrow(toxicity) < 1) {
      toxicity = data.frame(Metric = NA_character_,
                            Value = NA_real_,
                            Unit = NA_character_,
                            Route = NA_character_,
                            Species = NA_character_,
                            stringsAsFactors = FALSE)
    }
    for (hazard_code in hazard_codes) {
      for (j in seq_len(nrow(toxicity))) {
        rows[[length(rows) + 1]] = data.frame(
          Query = row$Query,
          CID = suppressWarnings(as.integer(row$CID)),
          HazardCode = hazard_code,
          HazardCategory = .normalized_hazard_category(hazard_code, evidence),
          HazardGroup = .normalized_hazard_group(hazard_code, evidence),
          SignalWord = row$SignalWord,
          HazardClass = row$HazardClass,
          HazardStatement = row$HazardStatement,
          PrecautionCode = precaution_codes,
          ExposureRoute = .uaf_first_non_empty_text(routes, toxicity$Route[[j]]),
          TargetOrgan = organs,
          ToxicityMetric = toxicity$Metric[[j]],
          ToxicityValue = toxicity$Value[[j]],
          ToxicityUnit = toxicity$Unit[[j]],
          Species = toxicity$Species[[j]],
          Source = row$Source,
          EvidenceText = evidence,
          EvidenceURL = row$SourceURL,
          PubChemURL = row$PubChemURL,
          ExtractionRule = "ghs_and_toxicity_regex",
          Confidence = ifelse(!is.na(hazard_code), "high", "medium"),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "CID", "HazardCode", "HazardCategory",
                 "HazardGroup", "SignalWord", "PrecautionCode",
                 "ExposureRoute", "TargetOrgan", "ToxicityMetric",
                 "ToxicityValue", "ToxicityUnit", "Species"),
    collapse_cols = c("HazardClass", "HazardStatement", "Source",
                      "EvidenceText", "EvidenceURL", "PubChemURL",
                      "ExtractionRule")
  )
}

.normalized_hazard_category = function(code, text = "") {
  if (!is.na(code) && grepl("^H2", code)) return("physical_hazard")
  if (!is.na(code) && grepl("^H3", code)) return("health_hazard")
  if (!is.na(code) && grepl("^H4", code)) return("environmental_hazard")
  text = tolower(.uaf_squish_text(text))
  if (grepl("flammable|explosive|oxidiz", text)) return("physical_hazard")
  if (grepl("aquatic|environment", text)) return("environmental_hazard")
  if (grepl("toxic|harmful|irritat|corros|carcin", text)) return("health_hazard")
  NA_character_
}

.normalized_hazard_group = function(code, text = "") {
  code = .uaf_squish_text(code)
  text_l = tolower(.uaf_squish_text(text))
  if (!is.na(code) && code %in% c("H318", "H319")) return("eye_damage_irritation")
  if (!is.na(code) && code %in% c("H315", "H317")) return("skin_irritation_sensitization")
  if (!is.na(code) && grepl("^H30|^H31|^H33", code)) return("acute_toxicity")
  if (!is.na(code) && grepl("^H22|^H23|^H24|^H25|^H26|^H27", code)) return("flammability_reactivity")
  if (!is.na(code) && grepl("^H4", code)) return("aquatic_environmental_toxicity")
  if (grepl("eye", text_l)) return("eye_damage_irritation")
  if (grepl("skin|dermal", text_l)) return("skin_irritation_sensitization")
  if (grepl("flammable", text_l)) return("flammability_reactivity")
  if (grepl("toxic|harmful|fatal", text_l)) return("acute_toxicity")
  if (grepl("aquatic", text_l)) return("aquatic_environmental_toxicity")
  NA_character_
}

.normalized_target_organs = function(text) {
  organs = c("eye", "skin", "respiratory tract", "liver", "kidney",
             "central nervous system", "blood", "heart", "lung",
             "reproductive system")
  text_l = tolower(.uaf_squish_text(text))
  unique(organs[vapply(organs, function(organ) {
    grepl(paste0("\\b", organ, "\\b"), text_l)
  }, logical(1))])
}

.normalized_toxicity_measurements = function(text) {
  cols = c("Metric", "Value", "Unit", "Route", "Species")
  text = .uaf_squish_text(text)
  if (is.na(text) || text == "") return(.uaf_empty_table(cols))
  pattern = "\\b(LD50|LC50|NOAEL|LOAEL|EC50|IC50)\\b[^0-9]{0,50}([-+]?\\d*\\.?\\d+)\\s*([A-Za-z0-9./% -]+)"
  matches = gregexpr(pattern, text, perl = TRUE, ignore.case = TRUE)
  values = regmatches(text, matches)[[1]]
  if (length(values) < 1 || identical(values, character(0)) ||
      (length(values) == 1 && values[[1]] == "")) {
    return(.uaf_empty_table(cols))
  }
  rows = list()
  for (value in values) {
    parts = regmatches(value, regexec(pattern, value, perl = TRUE,
                                      ignore.case = TRUE))[[1]]
    if (length(parts) < 4) next
    unit = .uaf_squish_text(parts[[4]])
    unit = sub("\\s+(rat|mouse|rabbit|human|fish|daphnia|oral|dermal|inhalation).*$",
               "", unit, ignore.case = TRUE)
    rows[[length(rows) + 1]] = data.frame(
      Metric = toupper(parts[[2]]),
      Value = suppressWarnings(as.numeric(parts[[3]])),
      Unit = unit,
      Route = .pubchem_collapse(.pubchem_extract_routes(value)),
      Species = .normalized_species(value),
      stringsAsFactors = FALSE
    )
  }
  .normalized_bind_rows(rows, cols)
}

.normalized_species = function(text) {
  species = c("rat", "mouse", "rabbit", "human", "fish", "daphnia", "dog")
  text_l = tolower(.uaf_squish_text(text))
  found = species[vapply(species, function(x) grepl(paste0("\\b", x, "\\b"),
                                                    text_l), logical(1))]
  .pubchem_collapse(found)
}

.normalized_chemical_uses = function(pubchem_profiles) {
  cols = c("Query", "CID", "UseDomain", "UseType", "UseTerm", "UseGroup",
           "SourceTable", "Source", "EvidenceText", "EvidenceURL",
           "PubChemURL", "ExtractionRule", "Confidence")
  rows = list()

  fema = pubchem_profiles$FEMAProfile
  if (is.data.frame(fema) && nrow(fema) > 0) {
    for (i in seq_len(nrow(fema))) {
      row = fema[i, , drop = FALSE]
      use_values = list(descriptor = row$DescriptorTerms,
                        flavor = row$FlavorTerms,
                        odor = row$OdorTerms,
                        regulatory_status = row$GRASStatus)
      for (use_type in names(use_values)) {
        terms = .normalized_split_terms(use_values[[use_type]])
        if (use_type == "pharmacologic_class") {
          terms = .normalized_pharmacologic_terms(terms)
        }
        for (term in terms) {
          lookup = .normalized_term_lookup(term)
          rows[[length(rows) + 1]] = data.frame(
            Query = row$Query,
            CID = row$CID,
            UseDomain = ifelse(use_type == "regulatory_status", "regulatory",
                               "sensory"),
            UseType = use_type,
            UseTerm = term,
            UseGroup = ifelse(use_type == "regulatory_status",
                              "flavor_regulatory", lookup$TermGroup),
            SourceTable = "FEMAProfile",
            Source = row$Source,
            EvidenceText = row$RawValue,
            EvidenceURL = row$SourceURL,
            PubChemURL = row$PubChemURL,
            ExtractionRule = "fema_discrete_terms",
            Confidence = lookup$Confidence,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  fda = pubchem_profiles$FDA_SPL_Profile
  if (is.data.frame(fda) && nrow(fda) > 0) {
    for (i in seq_len(nrow(fda))) {
      row = fda[i, , drop = FALSE]
      use_values = list(active_ingredient = row$ActiveIngredient,
                        pharmacologic_class = row$PharmacologicClass,
                        administration_route = row$Route,
                        dosage_form = row$DosageForm,
                        warning = row$WarningTerms)
      for (use_type in names(use_values)) {
        terms = .normalized_split_terms(use_values[[use_type]])
        if (use_type == "pharmacologic_class") {
          terms = .normalized_pharmacologic_terms(terms)
        }
        for (term in terms) {
          lookup = .normalized_term_lookup(term)
          rows[[length(rows) + 1]] = data.frame(
            Query = row$Query,
            CID = row$CID,
            UseDomain = ifelse(use_type == "warning", "safety",
                               "biomedical"),
            UseType = use_type,
            UseTerm = term,
            UseGroup = lookup$TermGroup,
            SourceTable = "FDA_SPL_Profile",
            Source = row$Source,
            EvidenceText = row$RawValue,
            EvidenceURL = row$SourceURL,
            PubChemURL = row$PubChemURL,
            ExtractionRule = "fda_spl_discrete_terms",
            Confidence = lookup$Confidence,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  mesh = pubchem_profiles$MeSHProfile
  if (is.data.frame(mesh) && nrow(mesh) > 0) {
    for (i in seq_len(nrow(mesh))) {
      row = mesh[i, , drop = FALSE]
      for (term in .normalized_split_terms(c(row$PharmacologicAction,
                                             row$Descriptor))) {
        lookup = .normalized_term_lookup(term)
        rows[[length(rows) + 1]] = data.frame(
          Query = row$Query,
          CID = row$CID,
          UseDomain = "biomedical",
          UseType = "pharmacologic_action",
          UseTerm = term,
          UseGroup = lookup$TermGroup,
          SourceTable = "MeSHProfile",
          Source = row$Source,
          EvidenceText = row$RawValue,
          EvidenceURL = row$SourceURL,
          PubChemURL = row$PubChemURL,
          ExtractionRule = "mesh_action_terms",
          Confidence = lookup$Confidence,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "CID", "UseDomain", "UseType", "UseTerm",
                 "UseGroup"),
    collapse_cols = c("SourceTable", "Source", "EvidenceText",
                      "EvidenceURL", "PubChemURL", "ExtractionRule")
  )
}

.normalized_pharmacologic_terms = function(terms) {
  terms = .uaf_non_empty(terms)
  if (length(terms) < 1) return(character())
  keep = vapply(terms, function(term) {
    route_or_form = length(.pubchem_extract_routes(term)) > 0 ||
      length(.pubchem_extract_dosage_forms(term)) > 0
    pharmacologic_signal = grepl("inhibitor|agonist|antagonist|agent|anti|blocker|modulator|class",
                                 term, ignore.case = TRUE)
    !route_or_form || pharmacologic_signal
  }, logical(1))
  terms[keep]
}

.normalized_chemical_bioassays = function(pubchem) {
  cols = c("Query", "CID", "AID", "SID", "PanelMemberID", "AssayName",
           "AssayType", "ActivityOutcome", "ActivityClass", "ActivityName",
           "ActivityValue", "ActivityUnit", "TargetAccession",
           "TargetGeneID", "TargetName", "TargetType", "TargetOtherID",
           "TargetTaxonomyID", "TargetOrganism", "TargetCommonName",
           "AssaySourceName", "AssaySourceID", "EndpointNames",
           "EndpointCount", "PubMedID", "RNAi", "SourceTable",
           "EvidenceText", "EvidenceURL", "ExtractionRule", "Confidence")
  bioactivity = if (is.list(pubchem) && is.data.frame(pubchem$bioactivity)) {
    pubchem$bioactivity
  } else {
    NULL
  }
  details = if (is.list(pubchem) &&
                is.data.frame(pubchem$bioassay_details)) {
    pubchem$bioassay_details
  } else {
    NULL
  }
  if (!is.data.frame(bioactivity) || nrow(bioactivity) < 1) {
    return(.uaf_empty_table(cols))
  }

  for (col in setdiff(.pubchem_bioactivity_cols(), colnames(bioactivity))) {
    bioactivity[[col]] = NA_character_
  }
  if (is.data.frame(details)) {
    for (col in setdiff(.pubchem_bioassay_detail_cols(), colnames(details))) {
      details[[col]] = NA_character_
    }
  }

  rows = list()
  for (i in seq_len(nrow(bioactivity))) {
    row = bioactivity[i, , drop = FALSE]
    detail = .normalized_bioassay_detail_for_row(row, details)
    has_detail = length(.uaf_non_empty(c(detail$AID, detail$TargetName,
                                         detail$EndpointNames,
                                         detail$AssaySourceName))) > 0
    outcome = .uaf_first_non_empty_text(row$ActivityOutcome, row$Value)
    activity_class = .uaf_first_non_empty_text(
      row$ActivityClass,
      .normalized_activity_class(outcome)
    )
    evidence = .pubchem_collapse(c(row$AssayName, row$ActivityName,
                                   outcome, row$TargetAccession,
                                   row$TargetGeneID, detail$TargetName,
                                   detail$TargetOrganism,
                                   detail$EndpointNames, row$PubMedID,
                                   detail$PubMedID))
    if (length(.uaf_non_empty(c(row$AID, evidence))) < 1) next

    rows[[length(rows) + 1]] = data.frame(
      Query = row$Query,
      CID = suppressWarnings(as.integer(row$CID)),
      AID = suppressWarnings(as.integer(row$AID)),
      SID = suppressWarnings(as.integer(row$SID)),
      PanelMemberID = row$PanelMemberID,
      AssayName = row$AssayName,
      AssayType = row$AssayType,
      ActivityOutcome = outcome,
      ActivityClass = activity_class,
      ActivityName = row$ActivityName,
      ActivityValue = suppressWarnings(as.numeric(row$ActivityValue)),
      ActivityUnit = row$ActivityUnit,
      TargetAccession = .uaf_first_non_empty_text(row$TargetAccession,
                                                  detail$TargetAccession),
      TargetGeneID = .uaf_first_non_empty_text(row$TargetGeneID,
                                               detail$TargetGeneID),
      TargetName = .normalized_assay_detail_target_name(
        detail$TargetName,
        detail$TargetDescription,
        row$AssayName,
        row$TargetAccession,
        row$TargetGeneID
      ),
      TargetType = detail$TargetType,
      TargetOtherID = detail$TargetOtherID,
      TargetTaxonomyID = detail$TargetTaxonomyID,
      TargetOrganism = detail$TargetOrganism,
      TargetCommonName = detail$TargetCommonName,
      AssaySourceName = detail$AssaySourceName,
      AssaySourceID = detail$AssaySourceID,
      EndpointNames = detail$EndpointNames,
      EndpointCount = suppressWarnings(as.integer(detail$EndpointCount)),
      PubMedID = .pubchem_collapse(c(row$PubMedID, detail$PubMedID)),
      RNAi = row$RNAi,
      SourceTable = .pubchem_collapse(c("PubChemBioactivity",
                                        if (has_detail) {
                                          "PubChemBioAssayDetails"
                                        })),
      EvidenceText = evidence,
      EvidenceURL = .pubchem_collapse(c(row$SourceURL,
                                        detail$AssayDetailURL)),
      ExtractionRule = .pubchem_collapse(c("pubchem_assay_summary",
                                           if (has_detail) {
                                             "pubchem_assay_detail"
                                           })),
      Confidence = ifelse(!is.na(row$AID) && !is.na(outcome) &&
                            length(.uaf_non_empty(c(detail$TargetName,
                                                    detail$EndpointNames))) > 0,
                          "high",
                          ifelse(!is.na(row$AID) && !is.na(outcome),
                                 "high", "medium")),
      stringsAsFactors = FALSE
    )
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "CID", "AID", "SID", "TargetAccession",
                 "TargetGeneID", "ActivityOutcome", "ActivityName",
                 "ActivityValue", "AssayName"),
    collapse_cols = c("PanelMemberID", "TargetName", "TargetType",
                      "TargetOtherID", "TargetTaxonomyID", "TargetOrganism",
                      "TargetCommonName", "AssaySourceName",
                      "AssaySourceID", "EndpointNames", "PubMedID", "RNAi",
                      "SourceTable", "EvidenceText", "EvidenceURL",
                      "ExtractionRule")
  )
}

.normalized_bioassay_detail_for_row = function(row, details) {
  cols = .pubchem_bioassay_detail_cols()
  empty = as.data.frame(stats::setNames(rep(list(NA_character_),
                                            length(cols)),
                                        cols),
                        stringsAsFactors = FALSE)
  if (!is.data.frame(details) || nrow(details) < 1) return(empty[1, , drop = FALSE])
  aid_rows = details[details$AID == row$AID[[1]], , drop = FALSE]
  if (nrow(aid_rows) < 1) return(empty[1, , drop = FALSE])

  query_match = .pubchem_key_value(aid_rows$Query) ==
    .pubchem_key_value(row$Query[[1]])
  cid_match = .pubchem_key_value(aid_rows$CID) ==
    .pubchem_key_value(row$CID[[1]])
  if (any(query_match & cid_match)) {
    aid_rows = aid_rows[query_match & cid_match, , drop = FALSE]
  }

  accession = .uaf_first_non_empty_text(row$TargetAccession)
  if (!is.na(accession) && "TargetAccession" %in% colnames(aid_rows)) {
    accession_match = aid_rows$TargetAccession == accession
    accession_match[is.na(accession_match)] = FALSE
    if (any(accession_match)) {
      return(aid_rows[which(accession_match)[[1]], , drop = FALSE])
    }
  }

  gene = .uaf_first_non_empty_text(row$TargetGeneID)
  if (!is.na(gene) && "TargetGeneID" %in% colnames(aid_rows)) {
    gene_match = aid_rows$TargetGeneID == gene
    gene_match[is.na(gene_match)] = FALSE
    if (any(gene_match)) {
      return(aid_rows[which(gene_match)[[1]], , drop = FALSE])
    }
  }

  aid_rows[1, , drop = FALSE]
}

.normalized_assay_detail_target_name = function(detail_name,
                                                detail_description,
                                                assay_name,
                                                target_accession,
                                                target_gene_id) {
  candidates = .uaf_non_empty(c(detail_description, detail_name))
  for (candidate in candidates) {
    acronym = .uaf_extract_pattern(
      candidate,
      "\\b[A-Z][A-Z0-9]{2,}(?:-[A-Z0-9]+)?\\b"
    )
    acronym = acronym[!acronym %in% c("HUMAN", "NAD", "CELL", "LINE",
                                      "CHEMBL", "CHAIN")]
    if (length(acronym) > 0) return(acronym[[1]])
  }

  clean = .normalized_clean_target_name(detail_name)
  if (!is.na(clean) &&
      !grepl("^Chain\\s+[A-Z],", clean, ignore.case = TRUE)) {
    return(clean)
  }

  .normalized_target_name(assay_name, target_accession, target_gene_id)
}

.normalized_chemical_bioactivities = function(chemical_bioassays) {
  cols = c("Query", "CID", "AID", "SID", "ActivityOutcome",
           "ActivityClass", "ActivityDirection", "BioactivityDomain",
           "AssayType", "ActivityName", "ActivityValue", "ActivityUnit",
           "TargetName", "TargetAccession", "TargetGeneID",
           "TargetTaxonomyID", "TargetOrganism", "EndpointNames",
           "PubMedID", "SourceTable", "EvidenceText", "EvidenceURL",
           "ExtractionRule", "Confidence")
  if (!is.data.frame(chemical_bioassays) ||
      nrow(chemical_bioassays) < 1) {
    return(.uaf_empty_table(cols))
  }

  rows = list()
  for (i in seq_len(nrow(chemical_bioassays))) {
    row = chemical_bioassays[i, , drop = FALSE]
    text = .pubchem_collapse(c(row$AssayName, row$ActivityName,
                               row$ActivityOutcome, row$TargetName,
                               row$TargetType, row$TargetOrganism,
                               row$EndpointNames))
    target_name = .uaf_first_non_empty_text(
      row$TargetName,
      .normalized_target_name(row$AssayName,
                              row$TargetAccession,
                              row$TargetGeneID)
    )
    rows[[length(rows) + 1]] = data.frame(
      Query = row$Query,
      CID = row$CID,
      AID = row$AID,
      SID = row$SID,
      ActivityOutcome = row$ActivityOutcome,
      ActivityClass = .normalized_activity_class(row$ActivityOutcome,
                                                 row$ActivityClass),
      ActivityDirection = .normalized_activity_direction(text),
      BioactivityDomain = .normalized_bioactivity_domain(text,
                                                         row$TargetAccession,
                                                         row$TargetGeneID),
      AssayType = row$AssayType,
      ActivityName = row$ActivityName,
      ActivityValue = row$ActivityValue,
      ActivityUnit = row$ActivityUnit,
      TargetName = target_name,
      TargetAccession = row$TargetAccession,
      TargetGeneID = row$TargetGeneID,
      TargetTaxonomyID = row$TargetTaxonomyID,
      TargetOrganism = row$TargetOrganism,
      EndpointNames = row$EndpointNames,
      PubMedID = row$PubMedID,
      SourceTable = "ChemicalBioassays",
      EvidenceText = row$EvidenceText,
      EvidenceURL = row$EvidenceURL,
      ExtractionRule = "pubchem_activity_call",
      Confidence = row$Confidence,
      stringsAsFactors = FALSE
    )
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "CID", "AID", "SID", "ActivityClass",
                 "ActivityDirection", "BioactivityDomain", "TargetName",
                 "TargetAccession", "TargetGeneID", "ActivityValue"),
    collapse_cols = c("ActivityOutcome", "AssayType", "ActivityName",
                      "TargetTaxonomyID", "TargetOrganism", "EndpointNames",
                      "PubMedID", "SourceTable", "EvidenceText",
                      "EvidenceURL", "ExtractionRule")
  )
}

.normalized_chemical_targets = function(chemical_bioassays) {
  cols = c("Query", "CID", "TargetName", "TargetType", "TargetAccession",
           "TargetOtherID", "TargetGeneID", "TargetTaxonomyID",
           "TargetOrganism", "TargetCommonName", "TargetGroup", "AIDCount",
           "ActiveAssayCount", "InactiveAssayCount", "InconclusiveAssayCount",
           "AssayTypes", "ActivityDirections", "BioactivityDomains",
           "EndpointNames", "AssaySources", "SourceTable", "EvidenceText",
           "EvidenceURL", "ExtractionRule", "Confidence")
  if (!is.data.frame(chemical_bioassays) ||
      nrow(chemical_bioassays) < 1) {
    return(.uaf_empty_table(cols))
  }

  target_name = .uaf_squish_text(chemical_bioassays$TargetName)
  inferred_name = mapply(.normalized_target_name,
                         chemical_bioassays$AssayName,
                         chemical_bioassays$TargetAccession,
                         chemical_bioassays$TargetGeneID,
                         USE.NAMES = FALSE)
  target_name[is.na(target_name) | target_name == ""] =
    inferred_name[is.na(target_name) | target_name == ""]
  keep = (!is.na(target_name) & target_name != "") |
    !is.na(chemical_bioassays$TargetAccession) |
    !is.na(chemical_bioassays$TargetGeneID)
  if (!any(keep)) return(.uaf_empty_table(cols))

  source = chemical_bioassays[keep, , drop = FALSE]
  source$TargetName = target_name[keep]
  domains = mapply(.normalized_bioactivity_domain,
                   source$AssayName,
                   source$TargetAccession,
                   source$TargetGeneID,
                   USE.NAMES = FALSE)
  directions = vapply(source$AssayName, .normalized_activity_direction,
                      character(1))
  group_keys = paste(.pubchem_key_value(source$Query),
                     .pubchem_key_value(source$CID),
                     .pubchem_key_value(source$TargetName),
                     .pubchem_key_value(source$TargetAccession),
                     .pubchem_key_value(source$TargetGeneID),
                     sep = "\r")

  rows = list()
  for (group_key in unique(group_keys)) {
    group = source[group_keys == group_key, , drop = FALSE]
    group_domains = domains[group_keys == group_key]
    group_directions = directions[group_keys == group_key]
    activity_classes = .uaf_non_empty(group$ActivityClass)
    rows[[length(rows) + 1]] = data.frame(
      Query = group$Query[[1]],
      CID = group$CID[[1]],
      TargetName = .uaf_first_non_empty_text(group$TargetName),
      TargetType = .uaf_first_non_empty_text(group$TargetType,
                                             .normalized_target_type(group)),
      TargetAccession = .uaf_first_non_empty_text(group$TargetAccession),
      TargetOtherID = .uaf_first_non_empty_text(group$TargetOtherID),
      TargetGeneID = .uaf_first_non_empty_text(group$TargetGeneID),
      TargetTaxonomyID = .uaf_first_non_empty_text(group$TargetTaxonomyID),
      TargetOrganism = .uaf_first_non_empty_text(
        group$TargetOrganism,
        .normalized_target_organism(group$AssayName)
      ),
      TargetCommonName = .uaf_first_non_empty_text(group$TargetCommonName),
      TargetGroup = .uaf_first_non_empty_text(group_domains),
      AIDCount = length(unique(.uaf_non_empty(group$AID))),
      ActiveAssayCount = sum(activity_classes == "active"),
      InactiveAssayCount = sum(activity_classes == "inactive"),
      InconclusiveAssayCount = sum(activity_classes == "inconclusive"),
      AssayTypes = .pubchem_collapse(group$AssayType),
      ActivityDirections = .pubchem_collapse(group_directions),
      BioactivityDomains = .pubchem_collapse(group_domains),
      EndpointNames = .pubchem_collapse(group$EndpointNames),
      AssaySources = .pubchem_collapse(group$AssaySourceName),
      SourceTable = "ChemicalBioassays",
      EvidenceText = .pubchem_collapse(group$EvidenceText),
      EvidenceURL = .pubchem_collapse(group$EvidenceURL),
      ExtractionRule = "pubchem_target_summary",
      Confidence = ifelse(length(.uaf_non_empty(c(group$TargetName,
                                                  group$TargetAccession,
                                                  group$TargetGeneID,
                                                  group$TargetTaxonomyID))) > 0,
                          "high", "medium"),
      stringsAsFactors = FALSE
    )
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "CID", "TargetName", "TargetAccession",
                 "TargetGeneID"),
    collapse_cols = c("AssayTypes", "ActivityDirections",
                      "BioactivityDomains", "EndpointNames", "AssaySources",
                      "SourceTable", "EvidenceText", "EvidenceURL",
                      "ExtractionRule")
  )
}

.normalized_chemical_potencies = function(chemical_bioassays) {
  cols = c("Query", "CID", "AID", "ActivityOutcome", "ActivityClass",
           "PotencyMetric", "PotencyValue", "PotencyUnit", "PotencyBucket",
           "TargetName", "TargetAccession", "TargetGeneID", "AssayName",
           "AssayType", "SourceTable", "EvidenceText", "EvidenceURL",
           "ExtractionRule", "Confidence")
  if (!is.data.frame(chemical_bioassays) ||
      nrow(chemical_bioassays) < 1 ||
      !"ActivityValue" %in% colnames(chemical_bioassays)) {
    return(.uaf_empty_table(cols))
  }

  rows = list()
  for (i in seq_len(nrow(chemical_bioassays))) {
    row = chemical_bioassays[i, , drop = FALSE]
    value = suppressWarnings(as.numeric(row$ActivityValue))
    if (is.na(value)) next
    metric = .normalized_potency_metric(row$ActivityName, row$AssayName)
    rows[[length(rows) + 1]] = data.frame(
      Query = row$Query,
      CID = row$CID,
      AID = row$AID,
      ActivityOutcome = row$ActivityOutcome,
      ActivityClass = .normalized_activity_class(row$ActivityOutcome,
                                                 row$ActivityClass),
      PotencyMetric = metric,
      PotencyValue = value,
      PotencyUnit = row$ActivityUnit,
      PotencyBucket = .normalized_potency_bucket(value, row$ActivityUnit),
      TargetName = .uaf_first_non_empty_text(
        row$TargetName,
        .normalized_target_name(row$AssayName,
                                row$TargetAccession,
                                row$TargetGeneID)
      ),
      TargetAccession = row$TargetAccession,
      TargetGeneID = row$TargetGeneID,
      AssayName = row$AssayName,
      AssayType = row$AssayType,
      SourceTable = "ChemicalBioassays",
      EvidenceText = row$EvidenceText,
      EvidenceURL = row$EvidenceURL,
      ExtractionRule = "pubchem_activity_value",
      Confidence = ifelse(!is.na(metric) && !is.na(row$ActivityUnit),
                          "high", "medium"),
      stringsAsFactors = FALSE
    )
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "CID", "AID", "PotencyMetric", "PotencyValue",
                 "PotencyUnit", "TargetName", "TargetAccession",
                 "TargetGeneID"),
    collapse_cols = c("ActivityOutcome", "AssayName", "AssayType",
                      "SourceTable", "EvidenceText", "EvidenceURL",
                      "ExtractionRule")
  )
}

.normalized_activity_class = function(outcome, fallback = NA_character_) {
  fallback = .uaf_first_non_empty_text(fallback)
  if (!is.na(fallback) && fallback != "unspecified") return(fallback)
  outcome = tolower(.uaf_squish_text(outcome))
  if (is.na(outcome) || outcome == "") return("unspecified")
  if (grepl("inactive", outcome)) return("inactive")
  if (grepl("active", outcome)) return("active")
  if (grepl("inconclusive|unspecified|not tested|probe", outcome)) {
    return("inconclusive")
  }
  .normalized_clean_label(outcome)
}

.normalized_activity_direction = function(text) {
  text = tolower(.pubchem_collapse(text))
  if (is.na(text) || text == "") return(NA_character_)
  if (grepl("antagonist|inhibit|blocker|blockade", text)) return("inhibitor")
  if (grepl("agonist", text)) return("agonist")
  if (grepl("activator|activation|stimulat", text)) return("activator")
  if (grepl("binding|binds?", text)) return("binding")
  if (grepl("cytotoxic|toxicity|dead cells|viability|growth inhibition",
            text)) {
    return("cell_viability")
  }
  if (grepl("substrate", text)) return("substrate")
  NA_character_
}

.normalized_bioactivity_domain = function(text, target_accession = NA_character_,
                                          target_gene_id = NA_character_) {
  text = tolower(.pubchem_collapse(c(text, target_accession, target_gene_id)))
  if (is.na(text) || text == "") return(NA_character_)
  if (grepl("tumou?r|cancer|cell line|cytotoxic|viability|growth inhibition|nci-",
            text)) {
    return("cell_viability_cytotoxicity")
  }
  if (grepl("antimicrobial|antibacterial|antifungal|parasite|malaria",
            text)) {
    return("antimicrobial_parasitic")
  }
  if (grepl("kinase|jnk|mapk|phosphatase", text)) return("kinase_signaling")
  if (grepl("receptor|agonist|antagonist|channel|transporter", text)) {
    return("receptor_transporter")
  }
  if (grepl("enzyme|hydrolase|dehydrogenase|protease|oxidase|inhibitor",
            text) ||
      length(.uaf_non_empty(c(target_accession, target_gene_id))) > 0) {
    return("enzyme_protein_target")
  }
  if (grepl("inflamm|cox|prostaglandin", text)) return("inflammation")
  "general_bioassay"
}

.normalized_target_name = function(assay_name,
                                   target_accession = NA_character_,
                                   target_gene_id = NA_character_) {
  assay_name = .uaf_squish_text(assay_name)
  if (is.na(assay_name) || assay_name == "") {
    return(.uaf_first_non_empty_text(target_accession, target_gene_id))
  }

  cell_line = regmatches(
    assay_name,
    regexec("Data for the\\s+([^.;:]+?\\s+cell line)", assay_name,
            perl = TRUE, ignore.case = TRUE)
  )[[1]]
  if (length(cell_line) >= 2) return(.uaf_squish_text(cell_line[[2]]))

  pattern = "(?:Inhibitors?|Agonists?|Antagonists?|Activators?|Modulators?|Binders?)\\s+of\\s+([^:;,]+)"
  hit = regmatches(assay_name, regexec(pattern, assay_name, perl = TRUE,
                                       ignore.case = TRUE))[[1]]
  if (length(hit) >= 2) {
    candidate = .normalized_clean_target_name(hit[[2]])
    if (!is.na(candidate)) return(candidate)
  }

  acronym = regmatches(assay_name,
                       regexec("\\b([A-Z0-9]{2,}(?:-[A-Z0-9]+)?)\\b",
                               assay_name, perl = TRUE))[[1]]
  if (length(acronym) >= 2 &&
      !acronym[[2]] %in% c("HTS", "QHTS", "RNAI", "NCI")) {
    return(acronym[[2]])
  }

  .uaf_first_non_empty_text(target_accession, target_gene_id)
}

.normalized_clean_target_name = function(x) {
  x = .uaf_squish_text(x)
  if (is.na(x) || x == "") return(NA_character_)
  x = sub("\\s*\\([^)]*\\).*$", "", x, perl = TRUE)
  x = sub("\\s*\\[[^]]*\\].*$", "", x, perl = TRUE)
  x = sub("\\s+Assay.*$", "", x, ignore.case = TRUE, perl = TRUE)
  x = sub("\\s+Measurement.*$", "", x, ignore.case = TRUE, perl = TRUE)
  x = sub("\\s+Response.*$", "", x, ignore.case = TRUE, perl = TRUE)
  x = gsub("^(the|human|recombinant)\\s+", "", x, ignore.case = TRUE,
           perl = TRUE)
  x = .uaf_squish_text(x)
  if (is.na(x) || x == "" || nchar(x) > 80) return(NA_character_)
  x
}

.normalized_target_type = function(group) {
  text = tolower(.pubchem_collapse(c(group$AssayName,
                                     group$TargetAccession,
                                     group$TargetGeneID,
                                     group$TargetName)))
  if (grepl("cell line|nci-|hek|cho|hela|tumou?r|cancer", text)) {
    return("cell_line")
  }
  if (length(.uaf_non_empty(c(group$TargetAccession,
                              group$TargetGeneID))) > 0) {
    return("protein_or_gene")
  }
  if (grepl("receptor|enzyme|kinase|transporter|channel", text)) {
    return("protein_or_gene")
  }
  "assay_target"
}

.normalized_target_organism = function(text) {
  text = tolower(.pubchem_collapse(text))
  if (is.na(text) || text == "") return(NA_character_)
  if (grepl("\\bhuman\\b|\\bhek\\b|\\bhela\\b|\\bnci-", text)) return("human")
  if (grepl("\\bmouse\\b|\\bmurine\\b", text)) return("mouse")
  if (grepl("\\brat\\b", text)) return("rat")
  if (grepl("bacteria|bacterial|escherichia|staphylococcus", text)) {
    return("bacteria")
  }
  NA_character_
}

.normalized_potency_metric = function(activity_name, assay_name) {
  text = .pubchem_collapse(c(activity_name, assay_name))
  if (is.na(text) || text == "") return(NA_character_)
  metric = .uaf_extract_pattern(text, "\\b(?:IC50|EC50|AC50|GI50|LC50|Ki|Kd)\\b")
  if (length(metric) > 0) return(toupper(metric[[1]]))
  if (grepl("potency", text, ignore.case = TRUE)) return("potency")
  .normalized_clean_label(activity_name)
}

.normalized_potency_bucket = function(value, unit) {
  value = suppressWarnings(as.numeric(value))
  unit = tolower(.uaf_squish_text(unit))
  if (is.na(value)) return(NA_character_)
  if (.normalized_is_micromolar_unit(unit)) {
    if (value <= 0.1) return("sub_100_nM")
    if (value <= 1) return("sub_1_uM")
    if (value <= 10) return("1_to_10_uM")
    if (value <= 50) return("10_to_50_uM")
    return("above_50_uM")
  }
  "numeric_activity_value"
}

.normalized_is_micromolar_unit = function(unit) {
  unit = tolower(.uaf_squish_text(unit))
  if (is.na(unit) || unit == "") return(FALSE)
  micro_pattern = paste0("[", intToUtf8(0x00b5), intToUtf8(0x03bc), "]")
  unit = gsub(micro_pattern, "u", unit, perl = TRUE)
  unit %in% c("um", "u m", "u-m", "micromolar", "micromol/l")
}

.normalized_taxonomy_organism_keys = function(rows) {
  if (length(rows) < 1) return(character())
  table = do.call(rbind, rows)
  if (!is.data.frame(table) || nrow(table) < 1 ||
      !"Organism" %in% colnames(table)) {
    return(character())
  }
  organism = .uaf_non_empty(table$Organism)
  if (length(organism) < 1) return(character())
  keep = !is.na(table$Organism) & table$Organism != ""
  unique(.normalized_taxonomy_organism_key(table$Query[keep],
                                           table$CID[keep],
                                           table$Organism[keep]))
}

.normalized_taxonomy_organism_key = function(query, cid, organism) {
  paste(.pubchem_key_value(query),
        .pubchem_key_value(cid),
        tolower(.pubchem_key_value(organism)),
        sep = "\r")
}

.normalized_chemical_taxonomy = function(pubchem_profiles) {
  cols = c("Query", "CID", "TaxonomySystem", "TaxonomyID", "Domain",
           "Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
           "Species", "Organism", "CommonName", "TaxonomyRank",
           "TaxonomyGroup", "TaxonomyTerms", "NaturalProductClass",
           "Occurrence", "Source", "EvidenceText", "EvidenceURL",
           "PubChemURL", "ExtractionRule", "Confidence")
  lotus = pubchem_profiles$LOTUSProfile
  classifications = pubchem_profiles$PubChemClassificationProfile
  has_lotus = is.data.frame(lotus) && nrow(lotus) > 0
  has_classifications = is.data.frame(classifications) &&
    nrow(classifications) > 0
  if (!has_lotus && !has_classifications) return(.uaf_empty_table(cols))
  rows = list()

  if (has_lotus) {
    for (col in c("TaxonomyID", "CommonName", "TaxonomyRank",
                  "TaxonomyLineage")) {
      if (!col %in% colnames(lotus)) lotus[[col]] = NA_character_
    }

    compound_keys = paste(.pubchem_key_value(lotus$Query),
                          .pubchem_key_value(lotus$CID),
                          sep = "\r")
    taxon_specific = !is.na(lotus$TaxonomyID) & lotus$TaxonomyID != "" &
      !grepl(";", lotus$TaxonomyID, fixed = TRUE) &
      (!is.na(lotus$Organism) | !is.na(lotus$TaxonomyLineage))
    compounds_with_taxa = unique(compound_keys[taxon_specific])
    group_source = lotus[!(compound_keys %in% compounds_with_taxa &
                             !taxon_specific), , drop = FALSE]
    group_source_compound_keys = compound_keys[
      !(compound_keys %in% compounds_with_taxa & !taxon_specific)
    ]
    semantic_key = ifelse(!is.na(group_source$TaxonomyID) &
                            group_source$TaxonomyID != "",
                          group_source$TaxonomyID, "compound")
    group_keys = paste(.pubchem_key_value(group_source$Query),
                       .pubchem_key_value(group_source$CID),
                       .pubchem_key_value(semantic_key),
                       sep = "\r")
    for (group_key in unique(group_keys)) {
      group = group_source[group_keys == group_key, , drop = FALSE]
      all_group = lotus[compound_keys == group_source_compound_keys[
        which(group_keys == group_key)[[1]]
      ], , drop = FALSE]
      evidence = .pubchem_collapse(c(group$RawValue, group$Organism,
                                     group$Taxonomy, all_group$Occurrence,
                                     group$TaxonomyLineage,
                                     all_group$NaturalProductClass))
      lineage = .uaf_first_non_empty_text(group$TaxonomyLineage,
                                          group$Taxonomy)
      taxonomy_terms = .lotus_extract_taxonomy_terms(lineage, group$Taxonomy)
      taxonomy_text = .pubchem_collapse(taxonomy_terms)
      organism_candidates = unique(c(
        .lotus_extract_organisms(.pubchem_collapse(group$Organism),
                                 "organism species source organism biological source"),
        .lotus_extract_organisms(evidence,
                                 "organism species source organism biological source")
      ))
      organism_candidates = organism_candidates[
        vapply(organism_candidates, .lotus_is_organism_name, logical(1))
      ]
      organism = .uaf_first_non_empty_text(organism_candidates)
      genus_species = .normalized_genus_species(organism)
      lineage_ranks = .lotus_lineage_ranks(lineage, organism)
      domain = lineage_ranks$Domain
      kingdom = .uaf_first_non_empty_text(
        lineage_ranks$Kingdom,
        .lotus_extract_rank(evidence, taxonomy_terms, "kingdom"),
        .normalized_kingdom(evidence, taxonomy_text)
      )
      phylum = .uaf_first_non_empty_text(lineage_ranks$Phylum,
                                         .lotus_extract_rank(evidence,
                                                             taxonomy_terms,
                                                             "phylum"))
      tax_class = .uaf_first_non_empty_text(lineage_ranks$Class,
                                            .lotus_extract_rank(evidence,
                                                                taxonomy_terms,
                                                                "class"))
      order = .uaf_first_non_empty_text(lineage_ranks$Order,
                                        .lotus_extract_rank(evidence,
                                                            taxonomy_terms,
                                                            "order"))
      family = .uaf_first_non_empty_text(
        lineage_ranks$Family,
        .lotus_extract_family(evidence, taxonomy_terms),
        .lotus_extract_rank(evidence, taxonomy_terms, "family")
      )
      genus = .uaf_first_non_empty_text(
        genus_species$Genus,
        lineage_ranks$Genus,
        .lotus_extract_rank(evidence, taxonomy_terms, "genus")
      )
      species = .uaf_first_non_empty_text(
        genus_species$Species,
        lineage_ranks$Species,
        .lotus_extract_rank(evidence, taxonomy_terms, "species")
      )
      taxonomy_terms = unique(.uaf_non_empty(c(lineage_ranks$Terms,
                                               taxonomy_terms, domain, kingdom,
                                               phylum, tax_class, order, family,
                                               genus, species)))
      taxonomy_text = .pubchem_collapse(taxonomy_terms)
      natural_product_class = .pubchem_collapse(c(group$NaturalProductClass,
                                                  all_group$NaturalProductClass))
      confidence = ifelse(
        !is.na(.uaf_first_non_empty_text(group$TaxonomyID)) &&
          !is.na(organism) && !is.na(taxonomy_text),
        "high",
        ifelse(length(.uaf_non_empty(c(kingdom, phylum, tax_class, order,
                                      family, genus, species, taxonomy_text,
                                      natural_product_class))) > 0,
               "medium", "low")
      )
      rows[[length(rows) + 1]] = data.frame(
        Query = group$Query[[1]],
        CID = group$CID[[1]],
        TaxonomySystem = "LOTUS",
        TaxonomyID = .uaf_first_non_empty_text(group$TaxonomyID),
        Domain = domain,
        Kingdom = kingdom,
        Phylum = phylum,
        Class = tax_class,
        Order = order,
        Family = family,
        Genus = genus,
        Species = species,
        Organism = organism,
        CommonName = .uaf_first_non_empty_text(group$CommonName),
        TaxonomyRank = .uaf_first_non_empty_text(group$TaxonomyRank),
        TaxonomyGroup = .normalized_taxonomy_group(kingdom, taxonomy_text),
        TaxonomyTerms = taxonomy_text,
        NaturalProductClass = natural_product_class,
        Occurrence = .pubchem_collapse(c(group$Occurrence, all_group$Occurrence)),
        Source = .pubchem_collapse(group$Source),
        EvidenceText = evidence,
        EvidenceURL = .pubchem_collapse(group$SourceURL),
        PubChemURL = .pubchem_collapse(group$PubChemURL),
        ExtractionRule = "lotus_taxonomy_regex",
        Confidence = confidence,
        stringsAsFactors = FALSE
      )
    }
  }

  if (has_classifications) {
    existing_keys = .normalized_taxonomy_organism_keys(rows)
    biological = classifications[
      classifications$TreeType == "biological" &
        grepl("LOTUS", paste(classifications$Source,
                             classifications$TreeName,
                             classifications$ClassPath),
              ignore.case = TRUE),
      ,
      drop = FALSE
    ]
    for (i in seq_len(nrow(biological))) {
      row = biological[i, , drop = FALSE]
      organism = .uaf_first_non_empty_text(row$ClassName)
      if (!.lotus_is_organism_name(organism)) next
      lineage = row$ClassPath
      lineage_ranks = .lotus_lineage_ranks(lineage, organism)
      taxonomy_terms = unique(.uaf_non_empty(c(lineage_ranks$Terms,
                                               lineage_ranks$Domain,
                                               lineage_ranks$Kingdom,
                                               lineage_ranks$Phylum,
                                               lineage_ranks$Class,
                                               lineage_ranks$Order,
                                               lineage_ranks$Family,
                                               lineage_ranks$Genus,
                                               lineage_ranks$Species)))
      taxonomy_text = .pubchem_collapse(taxonomy_terms)
      kingdom = .uaf_first_non_empty_text(lineage_ranks$Kingdom,
                                          .normalized_kingdom(lineage,
                                                              taxonomy_text))
      organism_key = .normalized_taxonomy_organism_key(row$Query, row$CID,
                                                       organism)
      if (organism_key %in% existing_keys) next
      existing_keys = c(existing_keys, organism_key)

      rows[[length(rows) + 1]] = data.frame(
        Query = row$Query,
        CID = row$CID,
        TaxonomySystem = "LOTUS",
        TaxonomyID = NA_character_,
        Domain = lineage_ranks$Domain,
        Kingdom = kingdom,
        Phylum = lineage_ranks$Phylum,
        Class = lineage_ranks$Class,
        Order = lineage_ranks$Order,
        Family = lineage_ranks$Family,
        Genus = lineage_ranks$Genus,
        Species = lineage_ranks$Species,
        Organism = organism,
        CommonName = NA_character_,
        TaxonomyRank = "species",
        TaxonomyGroup = .normalized_taxonomy_group(kingdom, taxonomy_text),
        TaxonomyTerms = taxonomy_text,
        NaturalProductClass = NA_character_,
        Occurrence = row$ClassPath,
        Source = .pubchem_collapse(c(row$Source, row$TreeName)),
        EvidenceText = row$ClassPath,
        EvidenceURL = row$ClassificationURL,
        PubChemURL = row$PubChemURL,
        ExtractionRule = "lotus_biological_tree",
        Confidence = "medium",
        stringsAsFactors = FALSE
      )
    }
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "CID", "TaxonomySystem", "TaxonomyID",
                 "Kingdom", "Phylum",
                 "Class", "Order", "Family", "Genus", "Species",
                 "Organism", "TaxonomyGroup", "TaxonomyTerms",
                 "NaturalProductClass"),
    collapse_cols = c("Occurrence", "Source", "EvidenceText",
                      "EvidenceURL", "PubChemURL", "ExtractionRule")
  )
}

.normalized_chemical_occurrences = function(chemical_taxonomy) {
  cols = c("Query", "CID", "SourceDatabase", "TaxonomyID", "Organism",
           "CommonName", "Domain", "Kingdom", "Phylum", "Class", "Order",
           "Family", "Genus", "Species", "OccurrenceType",
           "NaturalProductClass", "EvidenceText", "EvidenceURL",
           "PubChemURL", "ExtractionRule", "Confidence")
  if (!is.data.frame(chemical_taxonomy) || nrow(chemical_taxonomy) < 1) {
    return(.uaf_empty_table(cols))
  }

  for (col in setdiff(cols, colnames(chemical_taxonomy))) {
    chemical_taxonomy[[col]] = NA_character_
  }

  rows = list()
  for (i in seq_len(nrow(chemical_taxonomy))) {
    row = chemical_taxonomy[i, , drop = FALSE]
    if (length(.uaf_non_empty(c(row$TaxonomyID, row$Organism,
                                row$Species, row$Genus, row$Family))) < 1) {
      next
    }

    rows[[length(rows) + 1]] = data.frame(
      Query = row$Query,
      CID = row$CID,
      SourceDatabase = .uaf_first_non_empty_text(row$TaxonomySystem, "LOTUS"),
      TaxonomyID = row$TaxonomyID,
      Organism = row$Organism,
      CommonName = row$CommonName,
      Domain = row$Domain,
      Kingdom = row$Kingdom,
      Phylum = row$Phylum,
      Class = row$Class,
      Order = row$Order,
      Family = row$Family,
      Genus = row$Genus,
      Species = row$Species,
      OccurrenceType = .normalized_occurrence_type(row),
      NaturalProductClass = row$NaturalProductClass,
      EvidenceText = .uaf_first_non_empty_text(row$Occurrence,
                                               row$EvidenceText),
      EvidenceURL = row$EvidenceURL,
      PubChemURL = row$PubChemURL,
      ExtractionRule = .normalized_occurrence_extraction_rule(row),
      Confidence = row$Confidence,
      stringsAsFactors = FALSE
    )
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "CID", "SourceDatabase", "TaxonomyID",
                 "Organism", "Family", "Genus", "Species",
                 "OccurrenceType", "NaturalProductClass"),
    collapse_cols = c("EvidenceText", "EvidenceURL", "PubChemURL",
                      "ExtractionRule")
  )
}

.normalized_occurrence_extraction_rule = function(row) {
  rule = .uaf_first_non_empty_text(row$ExtractionRule)
  if (!is.na(rule) && rule == "lotus_biological_tree") {
    return("lotus_biological_tree")
  }
  "lotus_occurrence_taxon_link"
}

.normalized_occurrence_type = function(row) {
  text = tolower(.pubchem_collapse(c(row$TaxonomySystem,
                                     row$NaturalProductClass,
                                     row$Occurrence,
                                     row$EvidenceText)))
  if (is.na(text) || text == "") return("unknown")
  if (grepl("lotus|natural product|reported in|occurrence", text,
            ignore.case = TRUE)) {
    return("natural_product")
  }
  if (grepl("food|flavo[u]?r|odor|aroma", text, ignore.case = TRUE)) {
    return("food_or_sensory")
  }
  "unknown"
}

.normalized_label_value = function(text, label) {
  text = .uaf_squish_text(text)
  if (is.na(text) || text == "") return(NA_character_)
  pattern = paste0("\\b", label, "\\s*[:=]?\\s*([^;|,.]+)")
  match = regmatches(text, regexec(pattern, text, perl = TRUE,
                                   ignore.case = TRUE))[[1]]
  if (length(match) < 2) return(NA_character_)
  .uaf_squish_text(match[[2]])
}

.normalized_genus_species = function(organism) {
  organism = .uaf_squish_text(organism)
  out = list(Genus = NA_character_, Species = NA_character_)
  if (is.na(organism) || organism == "") return(out)
  parts = strsplit(organism, "\\s+", perl = TRUE)[[1]]
  if (length(parts) >= 1 && grepl("^[A-Z][a-z-]+$", parts[[1]])) {
    out$Genus = parts[[1]]
  }
  if (length(parts) >= 2 && grepl("^[a-z-]+$", parts[[2]])) {
    out$Species = paste(parts[1:2], collapse = " ")
  }
  out
}

.normalized_kingdom = function(text, taxonomy = NA_character_) {
  combined = tolower(.pubchem_collapse(c(text, taxonomy)))
  if (grepl("plantae|plant", combined)) return("Plantae")
  if (grepl("fungi|fungus", combined)) return("Fungi")
  if (grepl("bacteria|bacterium", combined)) return("Bacteria")
  if (grepl("animalia|animal", combined)) return("Animalia")
  if (grepl("archaea", combined)) return("Archaea")
  NA_character_
}

.normalized_taxonomy_group = function(kingdom, taxonomy) {
  if (!is.na(kingdom) && kingdom == "Plantae") return("plant")
  if (!is.na(kingdom) && kingdom == "Fungi") return("fungus")
  if (!is.na(kingdom) && kingdom == "Bacteria") return("bacteria")
  if (!is.na(kingdom) && kingdom == "Animalia") return("animal")
  .normalized_clean_label(taxonomy)
}

.normalized_mesh_tree_group = function(tree_number) {
  tree_number = .uaf_squish_text(tree_number)
  if (is.na(tree_number) || tree_number == "") return(NA_character_)
  prefix = substr(tree_number, 1, 1)
  groups = c(
    A = "Anatomy",
    B = "Organisms",
    C = "Diseases",
    D = "Chemicals and Drugs",
    E = "Analytical Diagnostic and Therapeutic Techniques",
    F = "Psychiatry and Psychology",
    G = "Phenomena and Processes",
    H = "Disciplines and Occupations",
    I = "Anthropology Education Sociology and Social Phenomena",
    J = "Technology Industry and Agriculture",
    K = "Humanities",
    L = "Information Science",
    M = "Named Groups",
    N = "Health Care",
    V = "Publication Characteristics",
    Z = "Geographicals"
  )
  out = unname(groups[[prefix]])
  if (is.null(out)) return(NA_character_)
  out
}

.normalized_mesh_label = function(...) {
  terms = .normalized_split_terms(unlist(list(...), use.names = FALSE))
  terms = terms[!grepl("^[A-Z][0-9]{2}(?:\\.[0-9]+)+$", terms)]
  if (length(terms) < 1) return(NA_character_)
  terms[[1]]
}

.normalized_chemical_pathway_roles = function(kegg, cid_lookup) {
  cols = c("Query", "CID", "KEGG_ID", "PathwayID", "PathwayName",
           "PathwayGroup", "ReactionID", "ReactionName", "ECNumber",
           "EnzymeName", "EnzymeClass", "RoleType", "SourceTable",
           "EvidenceURL", "ExtractionRule", "Confidence")
  rows = list()

  if (is.data.frame(kegg$pathways) && nrow(kegg$pathways) > 0) {
    for (i in seq_len(nrow(kegg$pathways))) {
      row = kegg$pathways[i, , drop = FALSE]
      rows[[length(rows) + 1]] = data.frame(
        Query = row$Query,
        CID = .normalized_lookup_cid(row$Query, cid_lookup),
        KEGG_ID = row$KEGG_ID,
        PathwayID = row$PathwayID,
        PathwayName = row$PathwayName,
        PathwayGroup = row$PathwayGroup,
        ReactionID = NA_character_,
        ReactionName = NA_character_,
        ECNumber = NA_character_,
        EnzymeName = NA_character_,
        EnzymeClass = NA_character_,
        RoleType = "pathway_member",
        SourceTable = "KEGGPathways",
        EvidenceURL = row$EvidenceURL,
        ExtractionRule = "kegg_pathway_role",
        Confidence = "high",
        stringsAsFactors = FALSE
      )
    }
  }

  if (is.data.frame(kegg$reactions) && nrow(kegg$reactions) > 0) {
    for (i in seq_len(nrow(kegg$reactions))) {
      row = kegg$reactions[i, , drop = FALSE]
      rows[[length(rows) + 1]] = data.frame(
        Query = row$Query,
        CID = .normalized_lookup_cid(row$Query, cid_lookup),
        KEGG_ID = row$KEGG_ID,
        PathwayID = NA_character_,
        PathwayName = NA_character_,
        PathwayGroup = NA_character_,
        ReactionID = row$ReactionID,
        ReactionName = row$ReactionName,
        ECNumber = NA_character_,
        EnzymeName = NA_character_,
        EnzymeClass = NA_character_,
        RoleType = "reaction_member",
        SourceTable = "KEGGReactions",
        EvidenceURL = row$EvidenceURL,
        ExtractionRule = "kegg_reaction_role",
        Confidence = "high",
        stringsAsFactors = FALSE
      )
    }
  }

  if (is.data.frame(kegg$enzymes) && nrow(kegg$enzymes) > 0) {
    for (i in seq_len(nrow(kegg$enzymes))) {
      row = kegg$enzymes[i, , drop = FALSE]
      rows[[length(rows) + 1]] = data.frame(
        Query = row$Query,
        CID = .normalized_lookup_cid(row$Query, cid_lookup),
        KEGG_ID = row$KEGG_ID,
        PathwayID = NA_character_,
        PathwayName = NA_character_,
        PathwayGroup = NA_character_,
        ReactionID = NA_character_,
        ReactionName = NA_character_,
        ECNumber = row$ECNumber,
        EnzymeName = row$EnzymeName,
        EnzymeClass = row$EnzymeClass,
        RoleType = "enzyme_associated",
        SourceTable = "KEGGEnzymes",
        EvidenceURL = row$EvidenceURL,
        ExtractionRule = "kegg_enzyme_role",
        Confidence = "high",
        stringsAsFactors = FALSE
      )
    }
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "CID", "KEGG_ID", "PathwayID", "ReactionID",
                 "ECNumber", "RoleType"),
    collapse_cols = c("PathwayName", "PathwayGroup", "ReactionName",
                      "EnzymeName", "EnzymeClass", "SourceTable",
                      "EvidenceURL", "ExtractionRule")
  )
}

.normalized_kegg_reaction_participants = function(kegg) {
  cols = c("Query", "KEGG_ID", "ReactionID", "Side", "ParticipantID",
           "ParticipantName", "ParticipantRole", "Equation", "Definition",
           "EvidenceURL", "ExtractionRule", "Confidence")
  reactions = kegg$reactions
  if (!is.data.frame(reactions) || nrow(reactions) < 1) return(.uaf_empty_table(cols))

  rows = list()
  for (i in seq_len(nrow(reactions))) {
    row = reactions[i, , drop = FALSE]
    parsed = .normalized_reaction_participants(row$Equation,
                                               row$ReactionDefinition)
    if (nrow(parsed) < 1) next
    for (j in seq_len(nrow(parsed))) {
      rows[[length(rows) + 1]] = data.frame(
        Query = row$Query,
        KEGG_ID = row$KEGG_ID,
        ReactionID = row$ReactionID,
        Side = parsed$Side[[j]],
        ParticipantID = parsed$ParticipantID[[j]],
        ParticipantName = parsed$ParticipantName[[j]],
        ParticipantRole = .normalized_participant_role(row$KEGG_ID,
                                                       parsed$ParticipantID[[j]],
                                                       parsed$Side[[j]]),
        Equation = row$Equation,
        Definition = row$ReactionDefinition,
        EvidenceURL = row$EvidenceURL,
        ExtractionRule = "kegg_equation_participant_parser",
        Confidence = ifelse(!is.na(parsed$ParticipantID[[j]]), "high", "medium"),
        stringsAsFactors = FALSE
      )
    }
  }

  out = .normalized_bind_rows(rows, cols)
  .normalized_dedupe_table(
    out,
    key_cols = c("Query", "KEGG_ID", "ReactionID", "Side", "ParticipantID",
                 "ParticipantName", "ParticipantRole", "Equation",
                 "Definition"),
    collapse_cols = c("EvidenceURL", "ExtractionRule")
  )
}

.normalized_reaction_participants = function(equation, definition = NA_character_) {
  cols = c("Side", "ParticipantID", "ParticipantName")
  equation = .uaf_squish_text(equation)
  if (is.na(equation) || equation == "") return(.uaf_empty_table(cols))
  delimiter = .normalized_reaction_delimiter(equation)
  sides = strsplit(equation, delimiter, perl = TRUE)[[1]]
  if (length(sides) < 2) return(.uaf_empty_table(cols))

  definition_sides = list(character(), character())
  definition = .uaf_squish_text(definition)
  if (!is.na(definition) && definition != "") {
    def_delimiter = .normalized_reaction_delimiter(definition)
    split_definition = strsplit(definition, def_delimiter, perl = TRUE)[[1]]
    if (length(split_definition) >= 2) {
      definition_sides = lapply(split_definition[1:2], .normalized_reaction_names)
    }
  }

  rows = list()
  side_names = c("substrate", "product")
  for (side_index in 1:2) {
    pieces = .normalized_reaction_pieces(sides[[side_index]])
    names = definition_sides[[side_index]]
    for (piece_index in seq_along(pieces)) {
      ids = .uaf_extract_pattern(pieces[[piece_index]], "\\b[CDG]\\d{5}\\b")
      if (length(ids) < 1) next
      for (id in ids) {
        participant_name = if (length(names) >= piece_index) {
          names[[piece_index]]
        } else {
          NA_character_
        }
        rows[[length(rows) + 1]] = data.frame(
          Side = side_names[[side_index]],
          ParticipantID = id,
          ParticipantName = participant_name,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  .normalized_bind_rows(rows, cols)
}

.normalized_reaction_delimiter = function(text) {
  if (grepl("<=>", text, fixed = TRUE)) return("\\s*<=>\\s*")
  if (grepl("=>", text, fixed = TRUE)) return("\\s*=>\\s*")
  if (grepl("<=", text, fixed = TRUE)) return("\\s*<=\\s*")
  "\\s*=\\s*"
}

.normalized_reaction_pieces = function(side) {
  pieces = unlist(strsplit(side, "\\s+\\+\\s+", perl = TRUE), use.names = FALSE)
  pieces = gsub("^\\s*[0-9]+\\s+", "", pieces)
  .uaf_non_empty(pieces)
}

.normalized_reaction_names = function(side) {
  pieces = unlist(strsplit(side, "\\s+\\+\\s+", perl = TRUE), use.names = FALSE)
  pieces = gsub("^\\s*[0-9]+\\s+", "", pieces)
  .uaf_non_empty(pieces)
}

.normalized_participant_role = function(kegg_id, participant_id, side) {
  if (is.na(participant_id) || is.na(kegg_id)) return("reaction_partner")
  if (participant_id == kegg_id && side == "substrate") return("query_substrate")
  if (participant_id == kegg_id && side == "product") return("query_product")
  "reaction_partner"
}
