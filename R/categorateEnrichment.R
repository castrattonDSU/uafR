.categorate_research_enrichment = function(compounds,
                                           data_list,
                                           detail,
                                           cache,
                                           cache_dir,
                                           throttle,
                                           assay_detail_limit = 50,
                                           trait_matrix_profile = "core",
                                           trait_matrix_mode = "binary",
                                           trait_matrix_min_confidence = 0,
                                           trait_matrix_max_traits = Inf,
                                           request_fun,
                                           kegg_request_fun,
                                           pubchem_query_overrides = NULL,
                                           kegg_throttle = NULL,
                                           strict_sources = FALSE) {
  compounds = .uaf_clean_compounds(compounds)
  pubchem_profile = ifelse(detail == "full", "full", "safety")
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

  pubchem_cache_dir = if (is.null(cache_dir)) NULL else file.path(cache_dir, "pubchem")
  kegg_cache_dir = if (is.null(cache_dir)) NULL else file.path(cache_dir, "kegg")

  pubchem = tryCatch(
    pubchemProfile(compounds = compounds,
                   profile = pubchem_profile,
                   sections = pubchem_sections,
                   sources = pubchem_sources,
                   cache = cache,
                   cache_dir = pubchem_cache_dir,
                   throttle = throttle,
                   assay_detail_limit = assay_detail_limit,
                   query_overrides = pubchem_query_overrides,
                   request_fun = request_fun),
    error = function(error) {
      if (isTRUE(strict_sources)) stop(error)
      warning("PubChem enrichment failed: ", conditionMessage(error),
              call. = FALSE)
      .categorate_empty_pubchem_profile(compounds, pubchem_profile)
    }
  )
  pubchem_profiles = .pubchem_extract_profiles(pubchem)

  kegg_ids = .categorate_kegg_ids(data_list$KEGG)
  kegg = tryCatch(
    keggProfile(compounds = compounds,
                kegg_ids = kegg_ids,
                pubchem_profile = pubchem,
                cache = cache,
                cache_dir = kegg_cache_dir,
                throttle = max(.uaf_first_numeric(kegg_throttle, throttle),
                               0.35),
                request_fun = kegg_request_fun),
    error = function(error) {
      if (isTRUE(strict_sources)) stop(error)
      warning("KEGG enrichment failed: ", conditionMessage(error),
              call. = FALSE)
      .categorate_empty_kegg_profile()
    }
  )

  normalized = .categorate_normalized_outputs(compounds = compounds,
                                              data_list = data_list,
                                              pubchem = pubchem,
                                              kegg = kegg,
                                              pubchem_profiles = pubchem_profiles,
                                              trait_matrix_profile = trait_matrix_profile,
                                              trait_matrix_mode = trait_matrix_mode,
                                              trait_matrix_min_confidence = trait_matrix_min_confidence,
                                              trait_matrix_max_traits = trait_matrix_max_traits)
  source_coverage = .categorate_source_coverage(compounds = compounds,
                                                data_list = data_list,
                                                pubchem = pubchem,
                                                kegg = kegg,
                                                pubchem_profiles = pubchem_profiles)
  derived_groups = .categorate_derived_groups(compounds = compounds,
                                              data_list = data_list,
                                              pubchem = pubchem,
                                              kegg = kegg,
                                              source_coverage = source_coverage,
                                              pubchem_profiles = pubchem_profiles,
                                              normalized = normalized)
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
    Provenance = provenance
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

.uaf_first_numeric = function(..., default = NA_real_) {
  values = suppressWarnings(as.numeric(unlist(list(...), use.names = FALSE)))
  values = values[is.finite(values)]
  if (length(values) < 1) default else values[[1]]
}

.categorate_kegg_ids = function(kegg_table) {
  if (!is.data.frame(kegg_table) || nrow(kegg_table) < 1) return(character())
  value_cols = setdiff(colnames(kegg_table), "Chemical")
  .kegg_normalize_ids(unlist(kegg_table[value_cols], use.names = FALSE))
}

.categorate_empty_pubchem_profile = function(compounds, profile) {
  identity = data.frame(Query = compounds,
                        CID = NA_integer_,
                        MatchStatus = "not_run",
                        SourceURL = NA_character_,
                        stringsAsFactors = FALSE)
  out = list(
    identity = identity,
    properties = .pubchem_empty_table(c("Query", "CID", "SourceURL")),
    synonyms = .pubchem_empty_table(c("Query", "CID", "Synonym", "SourceURL")),
    annotations = .pubchem_empty_annotation_table(),
    spectra = .pubchem_empty_annotation_table(),
    safety = .pubchem_empty_annotation_table(),
    experimental = .pubchem_empty_annotation_table(),
    source_annotations = .pubchem_empty_annotation_table(),
    taxonomy = .pubchem_empty_table(c("Query", "CID", "TaxonomyID",
                                      "Organism", "CommonName", "Rank",
                                      "Domain", "Kingdom", "Phylum",
                                      "Class", "Order", "Family", "Genus",
                                      "Species", "Lineage",
                                      "SourceAnnotation", "Source",
                                      "SourceURL", "TaxonomyURL",
                                      "PubChemURL")),
    classifications = .pubchem_empty_table(.pubchem_classification_cols()),
    bioactivity = .pubchem_empty_table(.pubchem_bioactivity_cols()),
    bioassay_details = .pubchem_empty_table(.pubchem_bioassay_detail_cols()),
    provenance = .pubchem_empty_table(c("Table", "SourceURL")),
    profile = profile,
    retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  class(out) = c("uaf_pubchem_profile", class(out))
  out
}

.categorate_empty_kegg_profile = function() {
  list(
    matches = .uaf_empty_table(c("Query", "KEGG_ID", "Database", "MatchName",
                                 "MatchStatus", "MatchScore", "MatchRank",
                                 "SourceURL", "RetrievedAt")),
    search_candidates = .kegg_empty_search_candidates(),
    records = .uaf_empty_table(c("Query", "KEGG_ID", "Database", "Field",
                                 "Value", "CleanValue", "ValueNumeric",
                                 "UnitClean", "EvidenceURL", "RetrievedAt")),
    identifiers = .uaf_empty_table(c("Query", "KEGG_ID", "Database",
                                     "IdentifierType", "Identifier",
                                     "SourceField", "EvidenceURL",
                                     "RetrievedAt")),
    pathways = .uaf_empty_table(c("Query", "KEGG_ID", "Database", "PathwayID",
                                  "PathwayName", "PathwayGroup", "Evidence",
                                  "EvidenceURL", "RetrievedAt")),
    reactions = .uaf_empty_table(c("Query", "KEGG_ID", "Database",
                                   "ReactionID", "ReactionName",
                                   "ReactionDefinition", "Equation",
                                   "Evidence", "EvidenceURL",
                                   "RetrievedAt")),
    enzymes = .uaf_empty_table(c("Query", "KEGG_ID", "Database",
                                 "ECNumber", "EnzymeName", "EnzymeClass",
                                 "Evidence", "EvidenceURL", "RetrievedAt")),
    modules = .uaf_empty_table(c("Query", "KEGG_ID", "Database", "ModuleID",
                                 "ModuleName", "ModuleDefinition",
                                 "Evidence", "EvidenceURL", "RetrievedAt")),
    links = .uaf_empty_table(c("Query", "KEGG_ID", "Database",
                               "TargetDatabase", "TargetID", "TargetPrefix",
                               "SourceEntry", "EvidenceURL", "RetrievedAt")),
    link_metadata = .uaf_empty_table(c("TargetDatabase", "TargetID",
                                       "TargetEntry", "Name", "Definition",
                                       "Equation", "Class", "PathwayGroup",
                                       "EvidenceURL", "RetrievedAt")),
    classifications = .uaf_empty_table(c("Query", "KEGG_ID", "Database",
                                         "ClassificationType",
                                         "Classification", "Evidence",
                                         "EvidenceURL", "RetrievedAt")),
    provenance = .uaf_empty_table(c("Table", "SourceURL")),
    retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
}

.categorate_pubchem_identifiers = function(pubchem) {
  cols = c("Query", "CID", "IdentifierType", "Identifier", "SourceField",
           "SourceURL")
  rows = list()

  if (is.data.frame(pubchem$annotations) && nrow(pubchem$annotations) > 0) {
    for (i in seq_len(nrow(pubchem$annotations))) {
      annotation = pubchem$annotations[i, , drop = FALSE]
      identifiers = .uaf_extract_identifiers(paste(annotation$Name,
                                                   annotation$Value,
                                                   annotation$Source))
      if (nrow(identifiers) < 1) next
      for (j in seq_len(nrow(identifiers))) {
        rows[[length(rows) + 1]] = data.frame(
          Query = annotation$Query,
          CID = annotation$CID,
          IdentifierType = identifiers$IdentifierType[[j]],
          Identifier = identifiers$Identifier[[j]],
          SourceField = annotation$HeadingPath,
          SourceURL = .uaf_first_non_empty_text(annotation$SourceURL,
                                                annotation$PubChemURL),
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
          SourceURL = synonym$SourceURL,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = unique(do.call(rbind, rows))
  row.names(out) = NULL
  out
}

.categorate_source_coverage = function(compounds, data_list, pubchem, kegg,
                                       pubchem_profiles = NULL) {
  cols = c("Query", "Source", "Present", "RecordCount", "RetrievedAt")
  retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  rows = list()

  legacy_sources = c("reactives", "LOTUS", "KEGG", "FEMA", "FDA_SPL")
  for (query in compounds) {
    for (source in legacy_sources) {
      count = .categorate_source_count(data_list[[source]], source, query)
      rows[[length(rows) + 1]] = data.frame(
        Query = query,
        Source = source,
        Present = count > 0,
        RecordCount = count,
        RetrievedAt = retrieved_at,
        stringsAsFactors = FALSE
      )
    }

    source_tables = list(
      PubChemSpectra = pubchem$spectra,
      PubChemSafety = pubchem$safety,
      PubChemExperimental = pubchem$experimental,
      PubChemBioactivity = pubchem$bioactivity,
      PubChemBioAssayDetails = pubchem$bioassay_details,
      SafetyProfile = pubchem_profiles$SafetyProfile,
      FEMAProfile = pubchem_profiles$FEMAProfile,
      FDA_SPL_Profile = pubchem_profiles$FDA_SPL_Profile,
      LOTUSProfile = pubchem_profiles$LOTUSProfile,
      PubChemClassificationProfile = pubchem_profiles$PubChemClassificationProfile,
      MeSHProfile = pubchem_profiles$MeSHProfile,
      LiteratureProfile = pubchem_profiles$LiteratureProfile,
      KEGGPathways = kegg$pathways,
      KEGGReactions = kegg$reactions,
      KEGGEnzymes = kegg$enzymes,
      KEGGModules = kegg$modules
    )
    for (source in names(source_tables)) {
      table = source_tables[[source]]
      count = if (is.data.frame(table) &&
                  "Query" %in% colnames(table)) {
        sum(table$Query == query, na.rm = TRUE)
      } else {
        0
      }
      rows[[length(rows) + 1]] = data.frame(
        Query = query,
        Source = source,
        Present = count > 0,
        RecordCount = count,
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

.categorate_source_count = function(table, value_column, query) {
  if (!is.data.frame(table) ||
      nrow(table) < 1 ||
      !"Chemical" %in% colnames(table) ||
      !value_column %in% colnames(table)) {
    return(0)
  }
  values = table[table$Chemical == query, value_column, drop = TRUE]
  values = .uaf_non_empty(values)
  sum(!values %in% c("None", "NA", ""))
}

.categorate_derived_groups = function(compounds, data_list, pubchem, kegg,
                                      source_coverage,
                                      pubchem_profiles = NULL,
                                      normalized = NULL) {
  cols = c("Query", "CID", "MolecularFormula", "MolecularWeight", "ExactMass",
           "XLogP", "TPSA", "HBondDonorCount", "HBondAcceptorCount",
           "RotatableBondCount", "HeavyAtomCount", "is_natural_product",
           "natural_product_superclasses", "natural_product_classes",
           "natural_product_subclasses", "occurrence_count",
           "organism_count", "kingdom_count", "family_count", "genus_count",
           "kingdoms_observed", "families_observed", "genera_observed",
           "dominant_kingdom", "dominant_family", "taxonomic_breadth",
           "is_plant_occurring", "is_fungal_occurring", "is_bacterial_occurring",
           "is_flavor_ingredient", "is_drug_or_label_ingredient",
           "has_kegg_pathway", "has_bioactivity", "has_active_bioactivity",
           "has_safety_hazard", "has_literature", "bioassay_count",
           "active_assay_count", "inactive_assay_count",
           "inconclusive_assay_count", "target_count", "active_target_count",
           "assay_detail_count", "protein_target_count",
           "cell_line_target_count", "human_target_count",
           "target_organisms", "bioactivity_domains",
           "active_bioactivity_domains",
           "min_potency_uM", "median_potency_uM", "potent_activity_count",
           "pathway_count", "reaction_count", "enzyme_count", "source_count",
           "molecular_size_bin",
           "polarity_bin", "lipophilicity_bin", "volatility_proxy",
           "oxygenated", "nitrogenous", "sulfur_containing", "halogenated",
           "metabolic_context", "biomedical_context", "ecological_context",
           "sensory_context", "safety_context", "analytical_context",
           "bioactivity_context", "confidence_score", "evidence_score")
  rows = list()

  for (query in compounds) {
    properties = .categorate_property_row(pubchem$properties, query)
    cid = .categorate_identity_value(pubchem$identity, query, "CID")
    formula = .categorate_property_value(properties, "MolecularFormula")
    molecular_weight = .categorate_property_number(properties, "MolecularWeight")
    exact_mass = .categorate_property_number(properties, "ExactMass")
    xlogp = .categorate_property_number(properties, "XLogP")
    tpsa = .categorate_property_number(properties, "TPSA")
    hbd = .categorate_property_number(properties, "HBondDonorCount")
    hba = .categorate_property_number(properties, "HBondAcceptorCount")
    rotatable = .categorate_property_number(properties, "RotatableBondCount")
    heavy_atoms = .categorate_property_number(properties, "HeavyAtomCount")

    coverage = source_coverage[source_coverage$Query == query, , drop = FALSE]
    source_count = sum(coverage$Present, na.rm = TRUE)
    pathway_count = .categorate_query_count(kegg$pathways, query)
    reaction_count = .categorate_query_count(kegg$reactions, query)
    enzyme_count = .categorate_query_count(kegg$enzymes, query)

    has_lotus_profile = .categorate_profile_count(pubchem_profiles$LOTUSProfile, query) > 0
    has_fema_profile = .categorate_profile_count(pubchem_profiles$FEMAProfile, query) > 0
    has_fda_profile = .categorate_profile_count(pubchem_profiles$FDA_SPL_Profile, query) > 0
    has_safety_profile = .categorate_profile_count(pubchem_profiles$SafetyProfile, query) > 0
    has_literature_profile = .categorate_profile_count(pubchem_profiles$LiteratureProfile, query) > 0
    has_mesh_profile = .categorate_profile_count(pubchem_profiles$MeSHProfile, query) > 0

    is_natural_product = .categorate_source_present(coverage, "LOTUS") ||
      has_lotus_profile
    is_flavor = .categorate_source_present(coverage, "FEMA") ||
      has_fema_profile
    is_drug = .categorate_source_present(coverage, "FDA_SPL") ||
      has_fda_profile ||
      .categorate_source_present(coverage, "PubChemBioactivity") ||
      any(kegg$records$Query == query & kegg$records$Database == "drug",
          na.rm = TRUE)
    has_safety = .categorate_source_present(coverage, "PubChemSafety") ||
      has_safety_profile
    bioactivity_summary = .categorate_bioactivity_summary(normalized, query)
    has_bioactivity = .categorate_source_present(coverage, "PubChemBioactivity") ||
      bioactivity_summary$bioassay_count > 0
    has_literature = .categorate_has_literature(query = query,
                                                pubchem = pubchem,
                                                kegg = kegg) ||
      has_literature_profile
    occurrence_summary = .categorate_occurrence_summary(normalized, query)
    np_class_summary = .categorate_natural_product_class_summary(normalized,
                                                                 query)
    metabolic_context = .categorate_metabolic_context(kegg, query)
    biomedical_context = .categorate_biomedical_context(pubchem_profiles,
                                                        is_drug,
                                                        has_mesh_profile,
                                                        query,
                                                        bioactivity_summary)
    ecological_context = .pubchem_collapse(c(
      .categorate_ecological_context(pubchem_profiles,
                                     is_natural_product,
                                     query),
      np_class_summary$natural_product_superclasses,
      np_class_summary$natural_product_classes,
      np_class_summary$natural_product_subclasses,
      occurrence_summary$dominant_kingdom,
      occurrence_summary$families_observed,
      occurrence_summary$taxonomic_breadth
    ))
    sensory_context = .categorate_sensory_context(pubchem_profiles,
                                                  is_flavor,
                                                  query)
    safety_context = .categorate_safety_context(pubchem_profiles,
                                                has_safety,
                                                query)
    analytical_context = .categorate_analytical_context(
      volatility_proxy = .categorate_volatility_proxy(molecular_weight,
                                                      tpsa, hbd),
      polarity_bin = .categorate_polarity_bin(tpsa),
      lipophilicity_bin = .categorate_lipophilicity_bin(xlogp)
    )
    confidence_score = .categorate_confidence_score(
      cid = cid,
      formula = formula,
      molecular_weight = molecular_weight,
      exact_mass = exact_mass,
      source_count = source_count,
      pathway_count = pathway_count,
      reaction_count = reaction_count,
      enzyme_count = enzyme_count,
      has_safety = has_safety,
      has_bioactivity = has_bioactivity
    )

    rows[[length(rows) + 1]] = data.frame(
      Query = query,
      CID = suppressWarnings(as.integer(cid)),
      MolecularFormula = formula,
      MolecularWeight = molecular_weight,
      ExactMass = exact_mass,
      XLogP = xlogp,
      TPSA = tpsa,
      HBondDonorCount = hbd,
      HBondAcceptorCount = hba,
      RotatableBondCount = rotatable,
      HeavyAtomCount = heavy_atoms,
      is_natural_product = is_natural_product,
      natural_product_superclasses = np_class_summary$natural_product_superclasses,
      natural_product_classes = np_class_summary$natural_product_classes,
      natural_product_subclasses = np_class_summary$natural_product_subclasses,
      occurrence_count = occurrence_summary$occurrence_count,
      organism_count = occurrence_summary$organism_count,
      kingdom_count = occurrence_summary$kingdom_count,
      family_count = occurrence_summary$family_count,
      genus_count = occurrence_summary$genus_count,
      kingdoms_observed = occurrence_summary$kingdoms_observed,
      families_observed = occurrence_summary$families_observed,
      genera_observed = occurrence_summary$genera_observed,
      dominant_kingdom = occurrence_summary$dominant_kingdom,
      dominant_family = occurrence_summary$dominant_family,
      taxonomic_breadth = occurrence_summary$taxonomic_breadth,
      is_plant_occurring = occurrence_summary$is_plant_occurring,
      is_fungal_occurring = occurrence_summary$is_fungal_occurring,
      is_bacterial_occurring = occurrence_summary$is_bacterial_occurring,
      is_flavor_ingredient = is_flavor,
      is_drug_or_label_ingredient = is_drug,
      has_kegg_pathway = pathway_count > 0,
      has_bioactivity = has_bioactivity,
      has_active_bioactivity = bioactivity_summary$active_assay_count > 0,
      has_safety_hazard = has_safety,
      has_literature = has_literature,
      bioassay_count = bioactivity_summary$bioassay_count,
      active_assay_count = bioactivity_summary$active_assay_count,
      inactive_assay_count = bioactivity_summary$inactive_assay_count,
      inconclusive_assay_count = bioactivity_summary$inconclusive_assay_count,
      target_count = bioactivity_summary$target_count,
      active_target_count = bioactivity_summary$active_target_count,
      assay_detail_count = bioactivity_summary$assay_detail_count,
      protein_target_count = bioactivity_summary$protein_target_count,
      cell_line_target_count = bioactivity_summary$cell_line_target_count,
      human_target_count = bioactivity_summary$human_target_count,
      target_organisms = bioactivity_summary$target_organisms,
      bioactivity_domains = bioactivity_summary$bioactivity_domains,
      active_bioactivity_domains = bioactivity_summary$active_bioactivity_domains,
      min_potency_uM = bioactivity_summary$min_potency_uM,
      median_potency_uM = bioactivity_summary$median_potency_uM,
      potent_activity_count = bioactivity_summary$potent_activity_count,
      pathway_count = pathway_count,
      reaction_count = reaction_count,
      enzyme_count = enzyme_count,
      source_count = source_count,
      molecular_size_bin = .categorate_mw_bin(molecular_weight),
      polarity_bin = .categorate_polarity_bin(tpsa),
      lipophilicity_bin = .categorate_lipophilicity_bin(xlogp),
      volatility_proxy = .categorate_volatility_proxy(molecular_weight,
                                                      tpsa, hbd),
      oxygenated = .categorate_formula_has(formula, "O"),
      nitrogenous = .categorate_formula_has(formula, "N"),
      sulfur_containing = .categorate_formula_has(formula, "S"),
      halogenated = .categorate_formula_has(formula, "F|Cl|Br|I"),
      metabolic_context = metabolic_context,
      biomedical_context = biomedical_context,
      ecological_context = ecological_context,
      sensory_context = sensory_context,
      safety_context = safety_context,
      analytical_context = analytical_context,
      bioactivity_context = bioactivity_summary$bioactivity_context,
      confidence_score = confidence_score,
      evidence_score = confidence_score,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.categorate_property_row = function(properties, query) {
  if (!is.data.frame(properties) ||
      nrow(properties) < 1 ||
      !"Query" %in% colnames(properties)) {
    return(data.frame())
  }
  properties[match(query, properties$Query), , drop = FALSE]
}

.categorate_identity_value = function(identity, query, column) {
  if (!is.data.frame(identity) ||
      nrow(identity) < 1 ||
      !"Query" %in% colnames(identity) ||
      !column %in% colnames(identity)) {
    return(NA_character_)
  }
  rows = identity[identity$Query == query, column, drop = TRUE]
  if (length(rows) < 1) return(NA_character_)
  rows[[1]]
}

.categorate_property_value = function(properties, column) {
  if (!is.data.frame(properties) ||
      nrow(properties) < 1 ||
      !column %in% colnames(properties)) {
    return(NA_character_)
  }
  value = properties[[column]][[1]]
  if (is.null(value) || is.na(value) || value == "") return(NA_character_)
  value
}

.categorate_property_number = function(properties, column) {
  value = .categorate_property_value(properties, column)
  if (is.na(value)) return(NA_real_)
  suppressWarnings(as.numeric(value))
}

.categorate_query_count = function(table, query) {
  if (!is.data.frame(table) ||
      nrow(table) < 1 ||
      !"Query" %in% colnames(table)) {
    return(0)
  }
  sum(table$Query == query, na.rm = TRUE)
}

.categorate_profile_count = function(table, query) {
  if (!is.data.frame(table) ||
      nrow(table) < 1 ||
      !"Query" %in% colnames(table)) {
    return(0)
  }
  sum(table$Query == query, na.rm = TRUE)
}

.categorate_occurrence_summary = function(normalized, query) {
  empty = list(
    occurrence_count = 0L,
    organism_count = 0L,
    kingdom_count = 0L,
    family_count = 0L,
    genus_count = 0L,
    kingdoms_observed = NA_character_,
    families_observed = NA_character_,
    genera_observed = NA_character_,
    dominant_kingdom = NA_character_,
    dominant_family = NA_character_,
    taxonomic_breadth = "none",
    is_plant_occurring = FALSE,
    is_fungal_occurring = FALSE,
    is_bacterial_occurring = FALSE
  )
  if (is.null(normalized) ||
      !is.list(normalized) ||
      is.null(normalized$ChemicalOccurrences) ||
      !is.data.frame(normalized$ChemicalOccurrences) ||
      nrow(normalized$ChemicalOccurrences) < 1) {
    return(empty)
  }

  occurrences = normalized$ChemicalOccurrences
  if (!"Query" %in% colnames(occurrences)) return(empty)
  rows = occurrences[occurrences$Query == query, , drop = FALSE]
  if (nrow(rows) < 1) return(empty)

  organisms = unique(.uaf_non_empty(c(rows$Organism, rows$Species)))
  kingdoms = unique(.uaf_non_empty(rows$Kingdom))
  families = unique(.uaf_non_empty(rows$Family))
  genera = unique(.uaf_non_empty(rows$Genus))

  empty$occurrence_count = nrow(rows)
  empty$organism_count = length(organisms)
  empty$kingdom_count = length(kingdoms)
  empty$family_count = length(families)
  empty$genus_count = length(genera)
  empty$kingdoms_observed = .pubchem_collapse(kingdoms)
  empty$families_observed = .pubchem_collapse(families)
  empty$genera_observed = .pubchem_collapse(genera)
  empty$dominant_kingdom = .categorate_mode(rows$Kingdom)
  empty$dominant_family = .categorate_mode(rows$Family)
  empty$taxonomic_breadth = .categorate_taxonomic_breadth(
    kingdom_count = empty$kingdom_count,
    family_count = empty$family_count,
    genus_count = empty$genus_count,
    organism_count = empty$organism_count
  )
  empty$is_plant_occurring = any(kingdoms %in% c("Plantae", "Viridiplantae"))
  empty$is_fungal_occurring = any(kingdoms %in% "Fungi")
  empty$is_bacterial_occurring = any(kingdoms %in% c("Bacteria", "Archaea"))
  empty
}

.categorate_natural_product_class_summary = function(normalized, query) {
  empty = list(
    natural_product_superclasses = NA_character_,
    natural_product_classes = NA_character_,
    natural_product_subclasses = NA_character_
  )
  if (is.null(normalized) ||
      !is.list(normalized) ||
      is.null(normalized$ChemicalClasses) ||
      !is.data.frame(normalized$ChemicalClasses) ||
      nrow(normalized$ChemicalClasses) < 1) {
    return(empty)
  }

  classes = normalized$ChemicalClasses
  if (!all(c("Query", "ClassSystem", "ClassType", "ClassName") %in%
           colnames(classes))) {
    return(empty)
  }
  rows = classes[classes$Query == query &
                   classes$ClassSystem == "LOTUS", , drop = FALSE]
  if (nrow(rows) < 1) return(empty)

  empty$natural_product_superclasses = .pubchem_collapse(
    unique(rows$ClassName[rows$ClassType == "natural_product_superclass"])
  )
  empty$natural_product_classes = .pubchem_collapse(
    unique(rows$ClassName[rows$ClassType == "natural_product_class"])
  )
  empty$natural_product_subclasses = .pubchem_collapse(
    unique(rows$ClassName[rows$ClassType == "natural_product_subclass"])
  )
  empty
}

.categorate_bioactivity_summary = function(normalized, query) {
  empty = list(
    bioassay_count = 0L,
    active_assay_count = 0L,
    inactive_assay_count = 0L,
    inconclusive_assay_count = 0L,
    target_count = 0L,
    active_target_count = 0L,
    assay_detail_count = 0L,
    protein_target_count = 0L,
    cell_line_target_count = 0L,
    human_target_count = 0L,
    target_organisms = NA_character_,
    bioactivity_domains = NA_character_,
    active_bioactivity_domains = NA_character_,
    min_potency_uM = NA_real_,
    median_potency_uM = NA_real_,
    potent_activity_count = 0L,
    bioactivity_context = NA_character_
  )
  if (is.null(normalized) || !is.list(normalized)) return(empty)

  bioactivities = normalized$ChemicalBioactivities
  targets = normalized$ChemicalTargets
  potencies = normalized$ChemicalPotencies
  bioassays = normalized$ChemicalBioassays

  if (is.data.frame(bioactivities) &&
      nrow(bioactivities) > 0 &&
      "Query" %in% colnames(bioactivities)) {
    bio_rows = bioactivities[bioactivities$Query == query, , drop = FALSE]
    classes = .uaf_non_empty(bio_rows$ActivityClass)
    domains = unique(.uaf_non_empty(bio_rows$BioactivityDomain))
    active_rows = bio_rows[bio_rows$ActivityClass == "active", , drop = FALSE]
    active_domains = unique(.uaf_non_empty(active_rows$BioactivityDomain))
    active_targets = if (nrow(active_rows) > 0) {
      target_keys = paste(.pubchem_key_value(active_rows$TargetName),
                          .pubchem_key_value(active_rows$TargetAccession),
                          .pubchem_key_value(active_rows$TargetGeneID),
                          sep = "\r")
      target_keys[target_keys != "<NA>\r<NA>\r<NA>"]
    } else {
      character()
    }

    empty$bioassay_count = nrow(bio_rows)
    empty$active_assay_count = sum(classes == "active")
    empty$inactive_assay_count = sum(classes == "inactive")
    empty$inconclusive_assay_count = sum(classes == "inconclusive")
    empty$active_target_count = length(unique(active_targets))
    empty$bioactivity_domains = .pubchem_collapse(domains)
    empty$active_bioactivity_domains = .pubchem_collapse(active_domains)
  }

  if (is.data.frame(targets) &&
      nrow(targets) > 0 &&
      "Query" %in% colnames(targets)) {
    target_rows = targets[targets$Query == query, , drop = FALSE]
    empty$target_count = nrow(target_rows)
    organisms = unique(.uaf_non_empty(target_rows$TargetOrganism))
    empty$target_organisms = .pubchem_collapse(organisms)
    empty$protein_target_count = sum(target_rows$TargetType == "protein_or_gene",
                                     na.rm = TRUE)
    empty$cell_line_target_count = sum(target_rows$TargetType == "cell_line",
                                       na.rm = TRUE)
    empty$human_target_count = sum(grepl("Homo sapiens|human",
                                         paste(target_rows$TargetOrganism,
                                               target_rows$TargetCommonName),
                                         ignore.case = TRUE),
                                   na.rm = TRUE)
    if (empty$active_target_count < 1) {
      empty$active_target_count = sum(target_rows$ActiveAssayCount > 0,
                                      na.rm = TRUE)
    }
  }

  if (is.data.frame(bioassays) &&
      nrow(bioassays) > 0 &&
      "Query" %in% colnames(bioassays)) {
    assay_rows = bioassays[bioassays$Query == query, , drop = FALSE]
    detailed = grepl("pubchem_assay_detail", assay_rows$ExtractionRule,
                     fixed = TRUE)
    detailed[is.na(detailed)] = FALSE
    empty$assay_detail_count = length(unique(.uaf_non_empty(
      assay_rows$AID[detailed]
    )))
  }

  if (is.data.frame(potencies) &&
      nrow(potencies) > 0 &&
      "Query" %in% colnames(potencies)) {
    potency_rows = potencies[potencies$Query == query, , drop = FALSE]
    uM = vapply(potency_rows$PotencyUnit,
                .normalized_is_micromolar_unit,
                logical(1))
    values = suppressWarnings(as.numeric(potency_rows$PotencyValue[uM]))
    values = values[!is.na(values)]
    if (length(values) > 0) {
      empty$min_potency_uM = min(values)
      empty$median_potency_uM = stats::median(values)
      empty$potent_activity_count = sum(values <= 10)
    }
  }

  empty$bioactivity_context = .pubchem_collapse(c(
    empty$bioactivity_domains,
    empty$active_bioactivity_domains,
    empty$target_organisms,
    if (empty$active_assay_count > 0) "active_assay",
    if (empty$potent_activity_count > 0) "potent_activity"
  ))
  empty
}

.categorate_mode = function(x) {
  x = .uaf_non_empty(x)
  if (length(x) < 1) return(NA_character_)
  tab = table(x)
  top_count = max(tab)
  top = names(tab)[tab == top_count]
  if (length(top) > 1) return("mixed")
  top[[1]]
}

.categorate_taxonomic_breadth = function(kingdom_count, family_count,
                                         genus_count, organism_count) {
  if (organism_count < 1) return("none")
  if (kingdom_count > 1) return("multi_kingdom")
  if (family_count > 1) return("multi_family")
  if (genus_count > 1) return("multi_genus")
  if (organism_count > 1) return("multi_species_same_genus")
  "single_organism"
}

.categorate_source_present = function(coverage, source) {
  rows = coverage[coverage$Source == source, , drop = FALSE]
  any(rows$Present, na.rm = TRUE)
}

.categorate_metabolic_context = function(kegg, query) {
  contexts = character()
  if (is.data.frame(kegg$classifications) && nrow(kegg$classifications) > 0) {
    contexts = c(contexts,
                 kegg$classifications$Classification[
                   kegg$classifications$Query == query &
                     kegg$classifications$ClassificationType == "KEGG pathway group"
                 ])
  }
  if (is.data.frame(kegg$pathways) && nrow(kegg$pathways) > 0) {
    contexts = c(contexts, kegg$pathways$PathwayGroup[kegg$pathways$Query == query])
  }
  .pubchem_collapse(contexts)
}

.categorate_biomedical_context = function(pubchem_profiles, is_drug,
                                          has_mesh_profile, query,
                                          bioactivity_summary = NULL) {
  contexts = character()
  if (isTRUE(is_drug)) contexts = c(contexts, "drug_or_label_ingredient")
  fda = pubchem_profiles$FDA_SPL_Profile
  if (is.data.frame(fda) && nrow(fda) > 0) {
    fda_query = fda[fda$Query == query, , drop = FALSE]
    contexts = c(contexts, fda_query$PharmacologicClass,
                 fda_query$Route, fda_query$DosageForm)
  }
  mesh = pubchem_profiles$MeSHProfile
  if (is.data.frame(mesh) && nrow(mesh) > 0) {
    mesh_query = mesh[mesh$Query == query, , drop = FALSE]
    contexts = c(contexts, mesh_query$PharmacologicAction)
  }
  if (isTRUE(has_mesh_profile)) contexts = c(contexts, "mesh_annotated")
  if (is.list(bioactivity_summary)) {
    contexts = c(contexts,
                 bioactivity_summary$bioactivity_domains,
                 bioactivity_summary$active_bioactivity_domains,
                 if (bioactivity_summary$active_assay_count > 0) {
                   "active_bioactivity"
                 },
                 if (bioactivity_summary$potent_activity_count > 0) {
                   "potent_bioactivity"
                 })
  }
  .pubchem_collapse(contexts)
}

.categorate_ecological_context = function(pubchem_profiles,
                                          is_natural_product,
                                          query) {
  contexts = character()
  if (isTRUE(is_natural_product)) contexts = c(contexts, "natural_product")
  lotus = pubchem_profiles$LOTUSProfile
  if (is.data.frame(lotus) && nrow(lotus) > 0) {
    lotus_query = lotus[lotus$Query == query, , drop = FALSE]
    contexts = c(contexts, lotus_query$NaturalProductClass,
                 lotus_query$Taxonomy, lotus_query$Organism)
  }
  .pubchem_collapse(contexts)
}

.categorate_sensory_context = function(pubchem_profiles, is_flavor, query) {
  contexts = character()
  if (isTRUE(is_flavor)) contexts = c(contexts, "flavor_ingredient")
  fema = pubchem_profiles$FEMAProfile
  if (is.data.frame(fema) && nrow(fema) > 0) {
    fema_query = fema[fema$Query == query, , drop = FALSE]
    contexts = c(contexts, fema_query$DescriptorTerms,
                 fema_query$FlavorTerms, fema_query$OdorTerms,
                 fema_query$GRASStatus)
  }
  .pubchem_collapse(contexts)
}

.categorate_safety_context = function(pubchem_profiles, has_safety, query) {
  contexts = character()
  if (isTRUE(has_safety)) contexts = c(contexts, "safety_annotated")
  safety = pubchem_profiles$SafetyProfile
  if (is.data.frame(safety) && nrow(safety) > 0) {
    safety_query = safety[safety$Query == query, , drop = FALSE]
    contexts = c(contexts, safety_query$SignalWord,
                 safety_query$HazardCode, safety_query$HazardClass)
  }
  .pubchem_collapse(contexts)
}

.categorate_analytical_context = function(volatility_proxy, polarity_bin,
                                          lipophilicity_bin) {
  .pubchem_collapse(c(volatility_proxy, polarity_bin, lipophilicity_bin))
}

.categorate_has_literature = function(query, pubchem, kegg) {
  pubchem_hit = is.data.frame(pubchem$annotations) &&
    nrow(pubchem$annotations) > 0 &&
    any(pubchem$annotations$Query == query &
          grepl("Literature|PubMed|PMID|Reference",
                paste(pubchem$annotations$Heading,
                      pubchem$annotations$HeadingPath,
                      pubchem$annotations$Name,
                      pubchem$annotations$Value,
                      pubchem$annotations$Source),
                ignore.case = TRUE),
        na.rm = TRUE)

  kegg_hit = is.data.frame(kegg$identifiers) &&
    nrow(kegg$identifiers) > 0 &&
    any(kegg$identifiers$Query == query &
          kegg$identifiers$IdentifierType %in% c("PubMed", "DOI"),
        na.rm = TRUE)
  pubchem_hit || kegg_hit
}

.categorate_formula_has = function(formula, element_pattern) {
  if (is.na(formula) || formula == "") return(FALSE)
  grepl(paste0("(", element_pattern, ")([0-9]|$)"), formula, perl = TRUE)
}

.categorate_mw_bin = function(molecular_weight) {
  if (is.na(molecular_weight)) return(NA_character_)
  if (molecular_weight < 150) return("small_<150_Da")
  if (molecular_weight < 350) return("medium_150_350_Da")
  if (molecular_weight < 750) return("large_350_750_Da")
  "very_large_>=750_Da"
}

.categorate_polarity_bin = function(tpsa) {
  if (is.na(tpsa)) return(NA_character_)
  if (tpsa < 40) return("low_TPSA_<40")
  if (tpsa < 90) return("moderate_TPSA_40_90")
  "high_TPSA_>=90"
}

.categorate_lipophilicity_bin = function(xlogp) {
  if (is.na(xlogp)) return(NA_character_)
  if (xlogp < 0) return("hydrophilic_XLogP_<0")
  if (xlogp < 3) return("balanced_XLogP_0_3")
  if (xlogp < 5) return("lipophilic_XLogP_3_5")
  "highly_lipophilic_XLogP_>=5"
}

.categorate_volatility_proxy = function(molecular_weight, tpsa, hbd) {
  if (is.na(molecular_weight) || is.na(tpsa)) return(NA_character_)
  if (molecular_weight < 200 && tpsa < 60 && (is.na(hbd) || hbd <= 1)) {
    return("more_volatile_proxy")
  }
  if (molecular_weight < 350 && tpsa < 90) return("moderate_volatile_proxy")
  "less_volatile_proxy"
}

.categorate_confidence_score = function(cid, formula, molecular_weight,
                                        exact_mass, source_count,
                                        pathway_count, reaction_count,
                                        enzyme_count, has_safety,
                                        has_bioactivity) {
  score = 0
  if (!is.na(suppressWarnings(as.integer(cid)))) score = score + 0.25
  if (!is.na(formula) && formula != "") score = score + 0.1
  if (!is.na(molecular_weight)) score = score + 0.1
  if (!is.na(exact_mass)) score = score + 0.1
  score = score + min(source_count, 5) * 0.05
  if (pathway_count > 0) score = score + 0.1
  if (reaction_count > 0) score = score + 0.05
  if (enzyme_count > 0) score = score + 0.05
  if (isTRUE(has_safety)) score = score + 0.05
  if (isTRUE(has_bioactivity)) score = score + 0.05
  round(min(score, 1), 2)
}

.categorate_research_provenance = function(pubchem, kegg) {
  cols = c("SourceSystem", "Table", "SourceURL")
  rows = list()
  if (is.data.frame(pubchem$provenance) && nrow(pubchem$provenance) > 0) {
    rows[[length(rows) + 1]] = data.frame(
      SourceSystem = "PubChem",
      Table = pubchem$provenance$Table,
      SourceURL = pubchem$provenance$SourceURL,
      stringsAsFactors = FALSE
    )
  }
  if (is.data.frame(kegg$provenance) && nrow(kegg$provenance) > 0) {
    rows[[length(rows) + 1]] = data.frame(
      SourceSystem = "KEGG",
      Table = kegg$provenance$Table,
      SourceURL = kegg$provenance$SourceURL,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = unique(do.call(rbind, rows))
  row.names(out) = NULL
  out
}
