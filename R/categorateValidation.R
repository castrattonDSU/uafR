#' Describe categorate output schemas
#'
#' @description
#' `chemicalDataDictionary()` returns the expected table and column contracts for
#' the main `categorate(detail = "research")` and `categorate(detail = "full")`
#' outputs. It is intended for users who need to understand, export, validate,
#' or join enriched chemical data without guessing what each table contains.
#'
#' @param tables Optional character vector of table names to return. If `NULL`,
#' all documented enriched output tables are returned.
#'
#' @return A data frame with table names, column names, expected types, required
#' flags, semantic roles, allowed values, and descriptions.
#'
#' @examples
#' \dontrun{
#' dictionary = chemicalDataDictionary()
#' chemicalDataDictionary("ChemicalTraitReport")
#' }
#'
#' @export
chemicalDataDictionary = function(tables = NULL) {
  dictionary = .categorate_data_dictionary()
  tables = .uaf_non_empty(tables)
  if (length(tables) > 0) {
    dictionary = dictionary[dictionary$Table %in% tables, , drop = FALSE]
  }
  row.names(dictionary) = NULL
  dictionary
}

#' Validate enriched categorate results
#'
#' @description
#' `validateCategorateResult()` audits a `categorate()` result against the uafR
#' data dictionary. It checks expected tables, required columns, basic column
#' types, controlled values, duplicate keys, table completeness, and source
#' coverage diagnostics. The return value is designed to be stored with results
#' and inspected before downstream analyses.
#'
#' @param x A list returned by `categorate()`, preferably with
#' `detail = "research"` or `detail = "full"`.
#' @param tables Optional character vector of table names to validate. If
#' `NULL`, all dictionary tables relevant to the result are checked.
#' @param strict Logical. If `TRUE`, missing optional documented columns are
#' reported as warnings.
#'
#' @return A list with `Summary`, `TableQuality`, `SourceDiagnostics`, `Issues`,
#' and `DataDictionary` data frames.
#'
#' @examples
#' \dontrun{
#' result = categorate(compounds, library_data, detail = "research")
#' audit = validateCategorateResult(result)
#' audit$Summary
#' audit$Issues
#' }
#'
#' @export
validateCategorateResult = function(x, tables = NULL, strict = FALSE) {
  dictionary = chemicalDataDictionary(tables)
  if (!is.list(x) || is.data.frame(x)) {
    issues = .categorate_validation_issue(
      severity = "error",
      table = NA_character_,
      column = NA_character_,
      issue = "Input is not a categorate result list",
      expected = "list",
      observed = class(x)[[1]],
      row_count = NA_integer_,
      examples = NA_character_
    )
    out = list(
      Summary = .categorate_validation_summary(
        table_quality = .uaf_empty_table(.categorate_table_quality_cols()),
        source_diagnostics = .uaf_empty_table(.categorate_source_diagnostic_cols()),
        issues = issues
      ),
      TableQuality = .uaf_empty_table(.categorate_table_quality_cols()),
      SourceDiagnostics = .uaf_empty_table(.categorate_source_diagnostic_cols()),
      Issues = issues,
      DataDictionary = dictionary
    )
    class(out) = c("uaf_categorate_validation", class(out))
    return(out)
  }

  table_names = unique(dictionary$Table)
  quality_rows = list()
  issue_rows = list()
  for (table_name in table_names) {
    table_dictionary = dictionary[dictionary$Table == table_name, ,
                                  drop = FALSE]
    table = x[[table_name]]
    quality_rows[[length(quality_rows) + 1]] =
      .categorate_table_quality(table_name, table, table_dictionary, strict)
    table_issues = .categorate_validate_table(table_name, table,
                                              table_dictionary, strict)
    if (nrow(table_issues) > 0) {
      issue_rows[[length(issue_rows) + 1]] = table_issues
    }
  }

  table_quality = .categorate_bind_quality(quality_rows,
                                           .categorate_table_quality_cols())
  source_diagnostics = .categorate_source_diagnostics(x)
  issues = .categorate_bind_quality(issue_rows,
                                    .categorate_validation_issue_cols())
  summary = .categorate_validation_summary(table_quality,
                                           source_diagnostics,
                                           issues)
  out = list(
    Summary = summary,
    TableQuality = table_quality,
    SourceDiagnostics = source_diagnostics,
    Issues = issues,
    DataDictionary = dictionary
  )
  class(out) = c("uaf_categorate_validation", class(out))
  out
}

#' @export
print.uaf_categorate_validation = function(x, ...) {
  summary = x$Summary
  if (is.data.frame(summary) && nrow(summary) > 0) {
    cat("uafR categorate validation\n")
    cat("Status:", summary$Status[[1]], "\n")
    cat("Tables present:", summary$TablesPresent[[1]], "/",
        summary$TablesExpected[[1]], "\n", sep = "")
    cat("Issues:", summary$IssueCount[[1]], " total; ",
        summary$ErrorCount[[1]], " error; ",
        summary$WarningCount[[1]], " warning\n", sep = "")
  } else {
    cat("uafR categorate validation\n")
  }
  invisible(x)
}

.categorate_data_dictionary = function() {
  specs = list(
    .categorate_schema("PubChemIdentity",
                       c("Query", "CID", "MatchStatus", "SourceURL"),
                       required = c("Query", "CID", "MatchStatus"),
                       role = "identity",
                       description = "PubChem compound resolution for each query."),
    .categorate_schema("PubChemProperties",
                       c("Query", "CID", "MolecularFormula",
                         "MolecularWeight", "IUPACName", "InChIKey",
                         "CanonicalSMILES", "IsomericSMILES", "ExactMass",
                         "XLogP", "TPSA", "HBondDonorCount",
                         "HBondAcceptorCount", "RotatableBondCount",
                         "HeavyAtomCount", "SourceURL"),
                       required = c("Query", "CID"),
                       role = "identity",
                       description = "Core PubChem structural and physicochemical properties."),
    .categorate_schema("PubChemSynonyms",
                       c("Query", "CID", "Synonym", "SourceURL"),
                       required = c("Query", "CID", "Synonym"),
                       role = "identity",
                       description = "PubChem synonyms and alternate names."),
    .categorate_schema("PubChemAnnotations", .pubchem_annotation_cols(),
                       required = c("Query", "CID", "Name", "CleanValue"),
                       role = "source_evidence",
                       description = "Parsed PubChem annotation rows."),
    .categorate_schema("PubChemSourceAnnotations", .pubchem_annotation_cols(),
                       required = c("Query", "CID", "Name", "CleanValue"),
                       role = "source_evidence",
                       description = "Parsed source-specific PubChem annotation rows."),
    .categorate_schema("PubChemSpectra", .pubchem_annotation_cols(),
                       required = c("Query", "CID", "Name", "CleanValue"),
                       role = "analytical",
                       description = "Mass-spectrometry and spectral annotations from PubChem."),
    .categorate_schema("PubChemSafety", .pubchem_annotation_cols(),
                       required = c("Query", "CID", "Name", "CleanValue"),
                       role = "safety",
                       description = "Safety and hazard annotations from PubChem."),
    .categorate_schema("PubChemExperimental", .pubchem_annotation_cols(),
                       required = c("Query", "CID", "Name", "CleanValue"),
                       role = "measurement",
                       description = "Experimental property annotations from PubChem."),
    .categorate_schema("PubChemBioactivity", .pubchem_bioactivity_cols(),
                       required = c("Query", "CID", "AID", "ActivityOutcome"),
                       role = "bioactivity",
                       description = "PubChem BioAssay summary rows."),
    .categorate_schema("PubChemBioAssayDetails", .pubchem_bioassay_detail_cols(),
                       required = c("AID", "AssayName"),
                       role = "bioactivity",
                       description = "Selected detailed PubChem BioAssay descriptions."),
    .categorate_schema("PubChemIdentifiers",
                       c("Query", "CID", "IdentifierType", "Identifier",
                         "SourceField", "RawValue", "Source", "SourceURL",
                         "PubChemURL"),
                       required = c("Query", "CID", "IdentifierType",
                                    "Identifier"),
                       role = "identifier",
                       description = "Identifiers extracted from PubChem annotations."),
    .categorate_schema("SafetyProfile",
                       c("Query", "CID", "SignalWord", "HazardCode",
                         "PrecautionCode", "HazardClass", "HazardStatement",
                         "RawValue", "Source", "SourceURL", "PubChemURL"),
                       required = c("Query", "CID"),
                       role = "safety",
                       description = "Cleaned safety, GHS, and hazard statements."),
    .categorate_schema("FEMAProfile",
                       c("Query", "CID", "FEMANumber", "GRASStatus",
                         "JECFANumber", "DescriptorTerms", "FlavorTerms",
                         "OdorTerms", "RawValue", "Source", "SourceURL",
                         "PubChemURL"),
                       required = c("Query", "CID"),
                       role = "sensory",
                       description = "FEMA flavor, odor, GRAS, and JECFA annotations."),
    .categorate_schema("FDA_SPL_Profile",
                       c("Query", "CID", "LabelSection", "ActiveIngredient",
                         "PharmacologicClass", "Route", "DosageForm",
                         "HasBoxedWarning", "WarningTerms", "RawValue",
                         "Source", "SourceURL", "PubChemURL"),
                       required = c("Query", "CID"),
                       role = "biomedical",
                       description = "FDA/SPL label-derived drug and route information."),
    .categorate_schema("LOTUSProfile",
                       c("Query", "CID", "LOTUS_ID", "TaxonomyID",
                         "Organism", "Taxonomy", "Kingdom", "Phylum",
                         "Class", "Order", "Family", "Genus", "Species",
                         "NaturalProductClass", "Occurrence", "ReferenceID",
                         "RawValue", "Source", "SourceURL", "PubChemURL"),
                       required = c("Query", "CID"),
                       role = "ecology",
                       description = "LOTUS occurrence, taxonomy, and natural-product annotations."),
    .categorate_schema("PubChemClassifications", .pubchem_classification_cols(),
                       required = c("Query", "CID", "ClassName"),
                       role = "classification",
                       description = "PubChem classification tree membership."),
    .categorate_schema("MeSHProfile",
                       c("Query", "CID", "Descriptor", "PharmacologicAction",
                         "TreeCategory", "TreeNumber", "RawValue", "Source",
                         "SourceURL", "PubChemURL"),
                       required = c("Query", "CID"),
                       role = "biomedical",
                       description = "MeSH descriptor, action, and tree-number annotations."),
    .categorate_schema("LiteratureProfile",
                       c("Query", "CID", "PubMedID", "DOI", "PatentID",
                         "ReferenceType", "TitleOrText", "RawValue", "Source",
                         "SourceURL", "PubChemURL"),
                       required = c("Query", "CID"),
                       role = "literature",
                       description = "Literature and patent identifiers from PubChem annotations."),
    .categorate_schema("ChemicalTerms",
                       c("Query", "CID", "Domain", "TermType", "TermRaw",
                         "TermClean", "TermGroup", "SourceTable", "Source",
                         "EvidenceText", "EvidenceURL", "PubChemURL",
                         "ExtractionRule", "Confidence"),
                       required = c("Query", "CID", "Domain", "TermType",
                                    "TermClean"),
                       role = "term",
                       description = "Long-form normalized source terms."),
    .categorate_schema("ChemicalTraits", .normalized_trait_cols(),
                       required = c("Query", "CID", "TraitType", "TraitGroup",
                                    "TraitValue", "TraitValueClean",
                                    "SourceDatabase", "Confidence",
                                    "ConfidenceScore"),
                       role = "trait",
                       description = "Cross-source discrete traits with evidence and confidence."),
    .categorate_schema("ChemicalTraitOntology", .normalized_ontology_cols(),
                       required = c("Query", "CID", "OntologyDomain",
                                    "OntologyGroup", "OntologyTerm",
                                    "OntologyKey", "SourceTraitKey",
                                    "Confidence", "ConfidenceScore"),
                       role = "ontology",
                       description = "Controlled-domain ontology terms mapped from source traits."),
    .categorate_schema("ChemicalTraitMatrix",
                       c("Query", "CID"),
                       required = c("Query", "CID"),
                       role = "matrix",
                       description = "Wide dynamic trait matrix; columns after Query/CID are trait keys."),
    .categorate_schema("ChemicalTraitOntologyMatrix",
                       c("Query", "CID"),
                       required = c("Query", "CID"),
                       role = "matrix",
                       description = "Wide dynamic ontology matrix; columns after Query/CID are ontology keys."),
    .categorate_schema("ChemicalTraitEvidence", .normalized_evidence_cols(),
                       required = c("Query", "CID", "EvidenceType",
                                    "AnalysisKey", "SourceTraitKey",
                                    "ConfidenceScore"),
                       role = "evidence",
                       description = "Audit trail from ontology or trait keys back to source evidence."),
    .categorate_schema("ChemicalTraitReport", .normalized_report_cols(),
                       required = c("Query", "CID", "TraitCount",
                                    "OntologyTermCount", "EvidenceCount",
                                    "SourceDatabases", "OntologyDomains"),
                       role = "report",
                       description = "One-row-per-compound researcher-facing trait report."),
    .categorate_schema("ChemicalTraitSummary",
                       c("Query", "CID", "SummaryLevel", "TraitType",
                         "SourceDatabase", "TraitCount", "MatrixEligibleCount",
                         "HighConfidenceCount", "MediumConfidenceCount",
                         "LowConfidenceCount", "MeanConfidenceScore",
                         "TopTraitGroups", "TopTraitValues", "SourceDatabases"),
                       required = c("Query", "CID", "SummaryLevel",
                                    "TraitCount"),
                       role = "summary",
                       description = "Compact trait breadth summaries by compound."),
    .categorate_schema("ChemicalTraitSimilarity",
                       c("QueryA", "CIDA", "QueryB", "CIDB",
                         "SharedTraitCount", "UnionTraitCount",
                         "OnlyA_TraitCount", "OnlyB_TraitCount",
                         "JaccardSimilarity", "OverlapCoefficient",
                         "SharedTraits", "OnlyA_Traits", "OnlyB_Traits"),
                       required = c("QueryA", "QueryB", "SharedTraitCount",
                                    "JaccardSimilarity"),
                       role = "similarity",
                       description = "Pairwise trait or ontology similarity between compounds."),
    .categorate_schema("ChemicalClasses",
                       c("Query", "CID", "ClassSystem", "ClassType",
                         "ClassID", "ClassName", "ClassPath", "SourceTable",
                         "EvidenceText", "EvidenceURL", "ExtractionRule",
                         "Confidence"),
                       required = c("Query", "CID", "ClassSystem",
                                    "ClassType", "ClassName"),
                       role = "classification",
                       description = "Discrete chemical class assignments."),
    .categorate_schema("ChemicalMeasurements",
                       .normalized_measurement_cols(),
                       required = c("Query", "CID", "Property", "Value"),
                       role = "measurement",
                       description = "Numeric and categorical measurements with standardized units and bins."),
    .categorate_schema("ChemicalMeasurementSummary",
                       .normalized_measurement_summary_cols(),
                       required = c("Query", "CID", "Property",
                                    "MeasurementClass", "MeasurementCount"),
                       role = "measurement",
                       description = "Per-compound measurement summaries using standardized units."),
    .categorate_schema("ChemicalHazards",
                       c("Query", "CID", "HazardCode", "HazardCategory",
                         "HazardGroup", "SignalWord", "HazardClass",
                         "HazardStatement", "PrecautionCode", "ExposureRoute",
                         "TargetOrgan", "ToxicityMetric", "ToxicityValue",
                         "ToxicityUnit", "Species", "Source", "EvidenceText",
                         "EvidenceURL", "PubChemURL", "ExtractionRule",
                         "Confidence"),
                       required = c("Query", "CID"),
                       role = "safety",
                       description = "Discrete hazard codes, routes, organs, and toxicity metrics."),
    .categorate_schema("ChemicalUses",
                       c("Query", "CID", "UseDomain", "UseType", "UseTerm",
                         "UseGroup", "SourceTable", "Source", "EvidenceText",
                         "EvidenceURL", "PubChemURL", "ExtractionRule",
                         "Confidence"),
                       required = c("Query", "CID", "UseDomain", "UseTerm"),
                       role = "use",
                       description = "Discrete chemical use, sensory, drug, and biomedical terms."),
    .categorate_schema("ChemicalBioassays",
                       c("Query", "CID", "AID", "SID", "PanelMemberID",
                         "AssayName", "AssayType", "ActivityOutcome",
                         "ActivityClass", "ActivityName", "ActivityDirection",
                         "ActivityValue", "ActivityUnit", "TargetName",
                         "TargetType", "TargetAccession", "TargetOtherID",
                         "TargetGeneID", "TargetTaxonomyID", "TargetOrganism",
                         "TargetCommonName", "EndpointNames",
                         "AssaySourceName", "PubMedID", "SourceTable",
                         "EvidenceText", "EvidenceURL", "ExtractionRule",
                         "Confidence"),
                       required = c("Query", "CID", "AID", "ActivityClass"),
                       role = "bioactivity",
                       description = "Normalized assay-level activity rows."),
    .categorate_schema("ChemicalBioactivities",
                       c("Query", "CID", "AID", "SID", "ActivityOutcome",
                         "ActivityClass", "AssayType", "ActivityName",
                         "ActivityDirection", "BioactivityDomain",
                         "TargetName", "TargetAccession", "TargetGeneID",
                         "TargetOrganism", "ActivityValue", "ActivityUnit",
                         "EvidenceText", "EvidenceURL", "ExtractionRule",
                         "Confidence"),
                       required = c("Query", "CID", "AID", "ActivityClass"),
                       role = "bioactivity",
                       description = "Discrete assay outcome and bioactivity-domain calls."),
    .categorate_schema("ChemicalTargets",
                       c("Query", "CID", "TargetName", "TargetType",
                         "TargetAccession", "TargetOtherID", "TargetGeneID",
                         "TargetTaxonomyID", "TargetOrganism",
                         "TargetCommonName", "TargetGroup", "AIDCount",
                         "ActiveAssayCount", "InactiveAssayCount",
                         "InconclusiveAssayCount", "AssayTypes",
                         "ActivityDirections", "BioactivityDomains",
                         "EndpointNames", "AssaySources", "SourceTable",
                         "EvidenceText", "EvidenceURL", "ExtractionRule",
                         "Confidence"),
                       required = c("Query", "CID", "TargetName"),
                       role = "bioactivity",
                       description = "Aggregated targets, genes, organisms, and assay counts."),
    .categorate_schema("ChemicalPotencies",
                       c("Query", "CID", "AID", "ActivityOutcome",
                         "ActivityClass", "AssayName", "AssayType",
                         "TargetName", "TargetAccession", "PotencyMetric",
                         "PotencyValue", "PotencyUnit", "ActivityDirection",
                         "EvidenceText", "EvidenceURL", "ExtractionRule",
                         "Confidence"),
                       required = c("Query", "CID", "AID", "PotencyMetric",
                                    "PotencyValue"),
                       role = "bioactivity",
                       description = "Numeric assay potency values."),
    .categorate_schema("ChemicalTaxonomy",
                       c("Query", "CID", "TaxonomySystem", "TaxonomyID",
                         "Domain", "Kingdom", "Phylum", "Class", "Order",
                         "Family", "Genus", "Species", "Organism",
                         "CommonName", "Rank", "Lineage", "TaxonomyTerm",
                         "TaxonomyRank", "NaturalProductClass", "SourceTable",
                         "Source", "EvidenceText", "EvidenceURL",
                         "PubChemURL", "ExtractionRule", "Confidence"),
                       required = c("Query", "CID"),
                       role = "taxonomy",
                       description = "Taxonomic terms and ranks from LOTUS/PubChem taxonomy sources."),
    .categorate_schema("ChemicalOccurrences",
                       c("Query", "CID", "SourceDatabase", "TaxonomyID",
                         "Organism", "Kingdom", "Phylum", "Class", "Order",
                         "Family", "Genus", "Species", "Occurrence",
                         "OccurrenceType", "EvidenceText", "EvidenceURL",
                         "PubChemURL", "ExtractionRule", "Confidence"),
                       required = c("Query", "CID", "SourceDatabase"),
                       role = "ecology",
                       description = "One row per chemical-organism occurrence signal."),
    .categorate_schema("ChemicalPathwayRoles",
                       c("Query", "CID", "KEGG_ID", "PathwayID",
                         "PathwayName", "PathwayGroup", "ReactionID",
                         "ReactionName", "ECNumber", "EnzymeName",
                         "EnzymeClass", "RoleType", "SourceTable",
                         "EvidenceURL", "ExtractionRule", "Confidence"),
                       required = c("Query", "KEGG_ID", "RoleType",
                                    "SourceTable"),
                       role = "metabolism",
                       description = "KEGG pathway, reaction, enzyme, and role terms."),
    .categorate_schema("KEGGReactionParticipants",
                       c("Query", "KEGG_ID", "ReactionID", "Side",
                         "ParticipantID", "ParticipantName", "EvidenceURL",
                         "ExtractionRule"),
                       required = c("Query", "KEGG_ID", "ReactionID",
                                    "Side", "ParticipantID"),
                       role = "metabolism",
                       description = "Parsed KEGG reaction substrates and products."),
    .categorate_schema("KEGGMatches",
                       c("Query", "KEGG_ID", "Database", "MatchName",
                         "MatchStatus", "MatchScore", "MatchRank",
                         "SourceURL", "RetrievedAt"),
                       required = c("Query", "KEGG_ID", "Database",
                                    "MatchStatus"),
                       role = "metabolism",
                       description = "KEGG name/ID resolution matches."),
    .categorate_schema("KEGGRecords",
                       c("Query", "KEGG_ID", "Database", "Field", "Value",
                         "CleanValue", "ValueNumeric", "UnitClean",
                         "EvidenceURL", "RetrievedAt"),
                       required = c("Query", "KEGG_ID", "Field", "Value"),
                       role = "metabolism",
                       description = "Parsed KEGG flat-file fields."),
    .categorate_schema("KEGGIdentifiers",
                       c("Query", "KEGG_ID", "Database", "IdentifierType",
                         "Identifier", "SourceField", "EvidenceURL",
                         "RetrievedAt"),
                       required = c("Query", "KEGG_ID", "IdentifierType",
                                    "Identifier"),
                       role = "identifier",
                       description = "Identifiers parsed from KEGG records."),
    .categorate_schema("KEGGPathways",
                       c("Query", "KEGG_ID", "Database", "PathwayID",
                         "PathwayName", "PathwayGroup", "Evidence",
                         "EvidenceURL", "RetrievedAt"),
                       required = c("Query", "KEGG_ID", "PathwayID"),
                       role = "metabolism",
                       description = "KEGG pathway membership."),
    .categorate_schema("KEGGReactions",
                       c("Query", "KEGG_ID", "Database", "ReactionID",
                         "ReactionName", "ReactionDefinition", "Equation",
                         "Evidence", "EvidenceURL", "RetrievedAt"),
                       required = c("Query", "KEGG_ID", "ReactionID"),
                       role = "metabolism",
                       description = "KEGG reaction membership and equations."),
    .categorate_schema("KEGGEnzymes",
                       c("Query", "KEGG_ID", "Database", "ECNumber",
                         "EnzymeName", "EnzymeClass", "Evidence",
                         "EvidenceURL", "RetrievedAt"),
                       required = c("Query", "KEGG_ID", "ECNumber"),
                       role = "metabolism",
                       description = "KEGG enzyme links and enzyme classes."),
    .categorate_schema("KEGGModules",
                       c("Query", "KEGG_ID", "Database", "ModuleID",
                         "ModuleName", "ModuleDefinition", "Evidence",
                         "EvidenceURL", "RetrievedAt"),
                       required = c("Query", "KEGG_ID", "ModuleID"),
                       role = "metabolism",
                       description = "KEGG module membership."),
    .categorate_schema("KEGGLinks",
                       c("Query", "KEGG_ID", "Database", "TargetDatabase",
                         "TargetID", "TargetPrefix", "SourceEntry",
                         "EvidenceURL", "RetrievedAt"),
                       required = c("Query", "KEGG_ID", "TargetDatabase",
                                    "TargetID"),
                       role = "metabolism",
                       description = "Raw KEGG linked database IDs."),
    .categorate_schema("KEGGLinkMetadata",
                       c("TargetDatabase", "TargetID", "TargetEntry", "Name",
                         "Definition", "Equation", "Class", "PathwayGroup",
                         "EvidenceURL", "RetrievedAt"),
                       required = c("TargetDatabase", "TargetID"),
                       role = "metabolism",
                       description = "Names, definitions, classes, and equations for linked KEGG IDs."),
    .categorate_schema("KEGGClassifications",
                       c("Query", "KEGG_ID", "Database", "ClassificationType",
                         "Classification", "Evidence", "EvidenceURL",
                         "RetrievedAt"),
                       required = c("Query", "KEGG_ID",
                                    "ClassificationType", "Classification"),
                       role = "metabolism",
                       description = "KEGG pathway, enzyme, and BRITE-style classification terms."),
    .categorate_schema("SourceCoverage",
                       c("Query", "Source", "Present", "RecordCount",
                         "RetrievedAt"),
                       required = c("Query", "Source", "Present",
                                    "RecordCount"),
                       role = "diagnostic",
                       description = "Per-query source presence and record counts."),
    .categorate_schema("DerivedGroups",
                       c("Query", "CID", "MolecularFormula",
                         "MolecularWeight", "ExactMass", "XLogP", "TPSA",
                         "source_count", "confidence_score",
                         "evidence_score"),
                       required = c("Query", "CID", "source_count",
                                    "confidence_score"),
                       role = "summary",
                       description = "Per-query derived grouping and context fields."),
    .categorate_schema("Provenance",
                       c("Table", "SourceURL"),
                       required = c("Table", "SourceURL"),
                       role = "provenance",
                       description = "Source URLs used for enrichment tables.")
  )
  out = do.call(rbind, specs)
  row.names(out) = NULL
  out
}

.categorate_schema = function(table, columns, required = character(),
                              role = "data", description = NA_character_) {
  columns = unique(columns)
  required = unique(required)
  data.frame(
    Table = table,
    Column = columns,
    Type = vapply(columns, .categorate_column_type, character(1)),
    Required = columns %in% required,
    Role = role,
    AllowedValues = vapply(columns, .categorate_column_allowed_values,
                           character(1)),
    Description = vapply(columns, .categorate_column_description,
                         character(1), table = table,
                         table_description = description),
    stringsAsFactors = FALSE
  )
}

.categorate_column_type = function(column) {
  if (column %in% c("MeanConfidence", "StandardValueLow",
                    "StandardValueHigh", "MinStandardValue",
                    "MedianStandardValue", "MaxStandardValue",
                    "BestValue")) {
    return("numeric")
  }
  if (grepl("Count$|^Count$|CID$|CIDA$|CIDB$|AID$|SID$|TaxonomyID$|GeneID$|RecordCount$|MatchRank$",
            column)) {
    return("integer")
  }
  if (grepl("Score$|Similarity$|Coefficient$|Mass$|Weight$|XLogP$|TPSA$|ValueNumeric$|Value$|uM$|PotencyValue$|ConfidenceScore$|evidence_score$|confidence_score$",
            column) &&
      !column %in% c("Value", "TraitValue", "SourceTraitValue",
                     "AnalysisLabel", "OntologyLabel", "CleanValue",
                     "RawValue")) {
    return("numeric")
  }
  if (grepl("^is_|^has_|Present$|Eligible$|Warning$|oxygenated$|nitrogenous$|halogenated$|sulfur_containing$",
            column)) {
    return("logical")
  }
  "character"
}

.categorate_column_allowed_values = function(column) {
  allowed = switch(
    column,
    Confidence = "low; medium; high",
    MappingConfidence = "low; medium; high",
    ActivityClass = "active; inactive; inconclusive; unspecified",
    EvidenceType = "ontology; trait",
    Present = "TRUE; FALSE",
    MatrixEligible = "TRUE; FALSE",
    NA_character_
  )
  allowed
}

.categorate_column_description = function(column, table,
                                          table_description) {
  descriptions = c(
    Query = "Original query chemical name.",
    CID = "PubChem compound identifier.",
    CIDA = "PubChem CID for the first compound in a pair.",
    CIDB = "PubChem CID for the second compound in a pair.",
    SourceDatabase = "Source database or provider for the row.",
    SourceTable = "uafR source table used to create the row.",
    EvidenceText = "Source-backed text used during extraction.",
    EvidenceURL = "URL for source evidence when available.",
    ExtractionRule = "Named extraction or mapping rule used by uafR.",
    Confidence = "Discrete confidence class assigned during extraction.",
    ConfidenceScore = "Numeric confidence score assigned during extraction.",
    OntologyKey = "Stable controlled ontology key for analysis.",
    MatrixKey = "Stable wide-matrix key for analysis.",
    TraitValueClean = "Normalized trait value used for grouping and matrices.",
    RecordCount = "Number of records observed for the query/source.",
    Present = "Whether source records were found for the query.",
    SourceURL = "URL used to retrieve or describe source data."
  )
  if (column %in% names(descriptions)) return(descriptions[[column]])
  if (!is.na(table_description)) {
    return(paste(table_description, "Column:", column))
  }
  paste(table, "column", column)
}

.categorate_validate_table = function(table_name, table, dictionary, strict) {
  issues = list()
  required_cols = dictionary$Column[dictionary$Required]
  optional_cols = dictionary$Column[!dictionary$Required]
  if (is.null(table)) {
    return(.categorate_validation_issue(
      severity = "error",
      table = table_name,
      column = NA_character_,
      issue = "Expected table is missing from result",
      expected = "data.frame",
      observed = "NULL",
      row_count = NA_integer_,
      examples = NA_character_
    ))
  }
  if (!is.data.frame(table)) {
    return(.categorate_validation_issue(
      severity = "error",
      table = table_name,
      column = NA_character_,
      issue = "Expected table is not a data frame",
      expected = "data.frame",
      observed = paste(class(table), collapse = "; "),
      row_count = NA_integer_,
      examples = NA_character_
    ))
  }

  missing_required = setdiff(required_cols, colnames(table))
  for (column in missing_required) {
    issues[[length(issues) + 1]] = .categorate_validation_issue(
      severity = "error",
      table = table_name,
      column = column,
      issue = "Required column is missing",
      expected = "present",
      observed = "missing",
      row_count = nrow(table),
      examples = NA_character_
    )
  }
  if (isTRUE(strict)) {
    missing_optional = setdiff(optional_cols, colnames(table))
    for (column in missing_optional) {
      issues[[length(issues) + 1]] = .categorate_validation_issue(
        severity = "warning",
        table = table_name,
        column = column,
        issue = "Documented optional column is missing",
        expected = "present when supported by the source",
        observed = "missing",
        row_count = nrow(table),
        examples = NA_character_
      )
    }
  }

  if (nrow(table) > 0) {
    present_dictionary = dictionary[dictionary$Column %in% colnames(table), ,
                                    drop = FALSE]
    for (i in seq_len(nrow(present_dictionary))) {
      column = present_dictionary$Column[[i]]
      type = present_dictionary$Type[[i]]
      type_issue = .categorate_column_type_issue(table[[column]], type)
      if (!is.na(type_issue)) {
        issues[[length(issues) + 1]] = .categorate_validation_issue(
          severity = "warning",
          table = table_name,
          column = column,
          issue = "Column values do not match expected type",
          expected = type,
          observed = type_issue,
          row_count = nrow(table),
          examples = .categorate_examples(table[[column]])
        )
      }
      allowed = .uaf_non_empty(strsplit(
        .uaf_first_non_empty_text(present_dictionary$AllowedValues[[i]]),
        "\\s*;\\s*",
        perl = TRUE
      )[[1]])
      if (length(allowed) > 0) {
        values = unique(.uaf_non_empty(table[[column]]))
        bad_values = setdiff(values, allowed)
        if (length(bad_values) > 0) {
          issues[[length(issues) + 1]] = .categorate_validation_issue(
            severity = "warning",
            table = table_name,
            column = column,
            issue = "Column contains values outside documented set",
            expected = .pubchem_collapse(allowed),
            observed = .pubchem_collapse(utils::head(bad_values, 8)),
            row_count = nrow(table),
            examples = .categorate_examples(bad_values)
          )
        }
      }
    }

    duplicate_info = .categorate_duplicate_info(table_name, table)
    if (duplicate_info$count > 0) {
      issues[[length(issues) + 1]] = .categorate_validation_issue(
        severity = "warning",
        table = table_name,
        column = .pubchem_collapse(duplicate_info$columns),
        issue = "Duplicate analysis keys detected",
        expected = "unique keys",
        observed = paste0(duplicate_info$count, " duplicate row(s)"),
        row_count = nrow(table),
        examples = duplicate_info$examples
      )
    }
  }

  .categorate_bind_quality(issues, .categorate_validation_issue_cols())
}

.categorate_column_type_issue = function(values, type) {
  values = .uaf_non_empty(values)
  if (length(values) < 1) return(NA_character_)
  if (type == "numeric") {
    parsed = suppressWarnings(as.numeric(values))
    if (any(is.na(parsed))) return("non-numeric values present")
  } else if (type == "integer") {
    parsed = suppressWarnings(as.numeric(values))
    if (any(is.na(parsed) | parsed != floor(parsed))) {
      return("non-integer values present")
    }
  } else if (type == "logical") {
    ok = tolower(values) %in% c("true", "false", "t", "f", "yes", "no",
                                "1", "0")
    if (!all(ok)) return("non-logical values present")
  }
  NA_character_
}

.categorate_table_quality = function(table_name, table, dictionary, strict) {
  cols = .categorate_table_quality_cols()
  if (is.null(table)) {
    return(data.frame(
      Table = table_name,
      Present = FALSE,
      RowCount = NA_integer_,
      ColumnCount = NA_integer_,
      RequiredColumnCount = sum(dictionary$Required),
      MissingRequiredColumns = .pubchem_collapse(
        dictionary$Column[dictionary$Required]
      ),
      MissingOptionalColumns = NA_character_,
      Completeness = NA_real_,
      DuplicateKeyColumns = NA_character_,
      DuplicateKeyCount = NA_integer_,
      QueryCount = NA_integer_,
      CIDCount = NA_integer_,
      SourceDatabaseCount = NA_integer_,
      HighConfidenceRows = NA_integer_,
      MediumConfidenceRows = NA_integer_,
      LowConfidenceRows = NA_integer_,
      Status = "missing",
      stringsAsFactors = FALSE
    )[, cols, drop = FALSE])
  }
  if (!is.data.frame(table)) {
    return(data.frame(
      Table = table_name,
      Present = FALSE,
      RowCount = NA_integer_,
      ColumnCount = NA_integer_,
      RequiredColumnCount = sum(dictionary$Required),
      MissingRequiredColumns = .pubchem_collapse(
        dictionary$Column[dictionary$Required]
      ),
      MissingOptionalColumns = NA_character_,
      Completeness = NA_real_,
      DuplicateKeyColumns = NA_character_,
      DuplicateKeyCount = NA_integer_,
      QueryCount = NA_integer_,
      CIDCount = NA_integer_,
      SourceDatabaseCount = NA_integer_,
      HighConfidenceRows = NA_integer_,
      MediumConfidenceRows = NA_integer_,
      LowConfidenceRows = NA_integer_,
      Status = "not_data_frame",
      stringsAsFactors = FALSE
    )[, cols, drop = FALSE])
  }

  required_cols = dictionary$Column[dictionary$Required]
  optional_cols = dictionary$Column[!dictionary$Required]
  missing_required = setdiff(required_cols, colnames(table))
  missing_optional = setdiff(optional_cols, colnames(table))
  required_present = intersect(required_cols, colnames(table))
  completeness = if (nrow(table) < 1 || length(required_present) < 1) {
    NA_real_
  } else {
    required_values = table[, required_present, drop = FALSE]
    complete = !is.na(required_values) & required_values != ""
    round(mean(as.matrix(complete)), 3)
  }
  duplicate_info = .categorate_duplicate_info(table_name, table)
  confidence = .categorate_confidence_counts(table)
  status = if (length(missing_required) > 0) {
    "schema_error"
  } else if (nrow(table) < 1) {
    "empty"
  } else if (duplicate_info$count > 0) {
    "duplicates"
  } else {
    "ok"
  }
  data.frame(
    Table = table_name,
    Present = TRUE,
    RowCount = nrow(table),
    ColumnCount = ncol(table),
    RequiredColumnCount = length(required_cols),
    MissingRequiredColumns = .pubchem_collapse(missing_required),
    MissingOptionalColumns = if (isTRUE(strict)) {
      .pubchem_collapse(missing_optional)
    } else {
      NA_character_
    },
    Completeness = completeness,
    DuplicateKeyColumns = .pubchem_collapse(duplicate_info$columns),
    DuplicateKeyCount = duplicate_info$count,
    QueryCount = .categorate_distinct_count(table, "Query"),
    CIDCount = .categorate_distinct_count(table, "CID"),
    SourceDatabaseCount = .categorate_distinct_count(table, "SourceDatabase"),
    HighConfidenceRows = confidence$high,
    MediumConfidenceRows = confidence$medium,
    LowConfidenceRows = confidence$low,
    Status = status,
    stringsAsFactors = FALSE
  )[, cols, drop = FALSE]
}

.categorate_duplicate_info = function(table_name, table) {
  key_cols = .categorate_duplicate_key_cols(table_name)
  key_cols = intersect(key_cols, colnames(table))
  empty = list(columns = key_cols, count = 0L, examples = NA_character_)
  if (!is.data.frame(table) || nrow(table) < 2 || length(key_cols) < 1) {
    return(empty)
  }
  semantic_cols = setdiff(key_cols, c("Query", "CID", "CIDA", "CIDB"))
  if (length(semantic_cols) > 0) {
    semantic_values = table[, semantic_cols, drop = FALSE]
    has_semantic_value = apply(semantic_values, 1, function(row) {
      length(.uaf_non_empty(row)) > 0
    })
    table = table[has_semantic_value, , drop = FALSE]
  }
  if (nrow(table) < 2) return(empty)
  key = .pubchem_row_key(table, key_cols, fallback_cols = key_cols)
  duplicated_key = duplicated(key)
  empty$count = sum(duplicated_key, na.rm = TRUE)
  if (empty$count > 0) {
    empty$examples = .pubchem_collapse(utils::head(key[duplicated_key], 5))
  }
  empty
}

.categorate_duplicate_key_cols = function(table_name) {
  keys = list(
    ChemicalTerms = c("Query", "CID", "Domain", "TermType", "TermClean"),
    ChemicalTraits = c("Query", "CID", "TraitType", "TraitGroup",
                       "TraitValueClean"),
    ChemicalTraitOntology = c("Query", "CID", "OntologyDomain",
                              "OntologyGroup", "OntologyTerm"),
    ChemicalTraitEvidence = c("Query", "CID", "EvidenceType", "AnalysisKey",
                              "SourceTraitKey"),
    ChemicalTraitReport = c("Query", "CID"),
    ChemicalClasses = c("Query", "CID", "ClassSystem", "ClassType",
                        "ClassID", "ClassName"),
    ChemicalMeasurements = c("Query", "CID", "Property", "Value", "Unit",
                             "Condition"),
    ChemicalMeasurementSummary = c("Query", "CID", "Property"),
    ChemicalHazards = c("Query", "CID", "HazardCode", "HazardCategory",
                        "ExposureRoute", "TargetOrgan", "ToxicityMetric",
                        "ToxicityValue", "HazardStatement", "Source",
                        "EvidenceText"),
    ChemicalUses = c("Query", "CID", "UseDomain", "UseType", "UseTerm",
                     "UseGroup"),
    ChemicalBioassays = c("Query", "CID", "AID", "SID", "ActivityClass",
                          "ActivityName", "TargetAccession"),
    ChemicalBioactivities = c("Query", "CID", "AID", "SID", "ActivityClass",
                              "BioactivityDomain", "TargetAccession"),
    ChemicalTargets = c("Query", "CID", "TargetName", "TargetAccession",
                        "TargetGeneID"),
    ChemicalPotencies = c("Query", "CID", "AID", "PotencyMetric",
                          "PotencyValue", "PotencyUnit"),
    ChemicalTaxonomy = c("Query", "CID", "TaxonomySystem", "TaxonomyID",
                         "Organism", "TaxonomyTerm", "TaxonomyRank"),
    ChemicalOccurrences = c("Query", "CID", "SourceDatabase", "TaxonomyID",
                            "Organism", "Family", "Genus", "Species"),
    ChemicalPathwayRoles = c("Query", "KEGG_ID", "PathwayID", "ReactionID",
                             "ECNumber", "RoleType", "RoleTerm"),
    KEGGReactionParticipants = c("Query", "KEGG_ID", "ReactionID", "Side",
                                 "ParticipantID"),
    SourceCoverage = c("Query", "Source")
  )
  if (table_name %in% names(keys)) return(keys[[table_name]])
  character()
}

.categorate_source_diagnostics = function(x) {
  cols = .categorate_source_diagnostic_cols()
  coverage = x$SourceCoverage
  if (!is.data.frame(coverage) || nrow(coverage) < 1) {
    return(.uaf_empty_table(cols))
  }
  for (col in setdiff(c("Query", "Source", "Present", "RecordCount",
                        "RetrievedAt"), colnames(coverage))) {
    coverage[[col]] = NA
  }
  rows = lapply(seq_len(nrow(coverage)), function(i) {
    row = coverage[i, , drop = FALSE]
    count = suppressWarnings(as.integer(row$RecordCount[[1]]))
    present = row$Present[[1]] %in% TRUE ||
      tolower(as.character(row$Present[[1]])) == "true"
    group = .categorate_source_group(row$Source[[1]])
    status = if (isTRUE(present) && !is.na(count) && count > 0) {
      "present"
    } else {
      "absent"
    }
    data.frame(
      Query = row$Query[[1]],
      Source = row$Source[[1]],
      SourceGroup = group,
      Present = present,
      RecordCount = count,
      QualityStatus = status,
      Diagnostic = if (status == "present") {
        "Source returned usable rows for this query."
      } else {
        "No usable rows were extracted for this query/source."
      },
      RetrievedAt = row$RetrievedAt[[1]],
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out[, cols, drop = FALSE]
}

.categorate_source_group = function(source) {
  source = .uaf_first_non_empty_text(source, "")
  if (grepl("KEGG", source, ignore.case = TRUE)) return("KEGG")
  if (grepl("LOTUS|Taxonomy|Occurrence", source, ignore.case = TRUE)) {
    return("ecology")
  }
  if (grepl("FEMA|Flavor|Spectra", source, ignore.case = TRUE)) {
    return("sensory_analytical")
  }
  if (grepl("FDA|SPL|MeSH|Bioactivity|BioAssay|Literature", source,
            ignore.case = TRUE)) {
    return("biomedical")
  }
  if (grepl("Safety|Hazard|Experimental", source, ignore.case = TRUE)) {
    return("safety_measurement")
  }
  if (grepl("reactives", source, ignore.case = TRUE)) return("reactivity")
  "general"
}

.categorate_validation_summary = function(table_quality, source_diagnostics,
                                          issues) {
  error_count = if (is.data.frame(issues) && nrow(issues) > 0) {
    sum(issues$Severity == "error", na.rm = TRUE)
  } else {
    0L
  }
  warning_count = if (is.data.frame(issues) && nrow(issues) > 0) {
    sum(issues$Severity == "warning", na.rm = TRUE)
  } else {
    0L
  }
  status = if (error_count > 0) {
    "fail"
  } else if (warning_count > 0) {
    "warn"
  } else {
    "pass"
  }
  data.frame(
    Status = status,
    TablesExpected = if (is.data.frame(table_quality)) nrow(table_quality) else 0L,
    TablesPresent = if (is.data.frame(table_quality) && nrow(table_quality) > 0) {
      sum(table_quality$Present, na.rm = TRUE)
    } else {
      0L
    },
    NonEmptyTables = if (is.data.frame(table_quality) && nrow(table_quality) > 0) {
      sum(table_quality$RowCount > 0, na.rm = TRUE)
    } else {
      0L
    },
    IssueCount = if (is.data.frame(issues)) nrow(issues) else 0L,
    ErrorCount = error_count,
    WarningCount = warning_count,
    SourceDiagnosticsRows = if (is.data.frame(source_diagnostics)) {
      nrow(source_diagnostics)
    } else {
      0L
    },
    RetrievedAt = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
}

.categorate_confidence_counts = function(table) {
  empty = list(high = NA_integer_, medium = NA_integer_, low = NA_integer_)
  if (!is.data.frame(table) ||
      nrow(table) < 1 ||
      !"Confidence" %in% colnames(table)) {
    return(empty)
  }
  confidence = tolower(.uaf_squish_text(table$Confidence))
  empty$high = sum(confidence == "high", na.rm = TRUE)
  empty$medium = sum(confidence == "medium", na.rm = TRUE)
  empty$low = sum(confidence == "low", na.rm = TRUE)
  empty
}

.categorate_distinct_count = function(table, column) {
  if (!is.data.frame(table) ||
      nrow(table) < 1 ||
      !column %in% colnames(table)) {
    return(NA_integer_)
  }
  length(unique(.uaf_non_empty(table[[column]])))
}

.categorate_examples = function(values, n = 5) {
  .pubchem_collapse(utils::head(unique(.uaf_non_empty(values)), n))
}

.categorate_validation_issue = function(severity, table, column, issue,
                                        expected, observed, row_count,
                                        examples) {
  data.frame(
    Severity = severity,
    Table = table,
    Column = column,
    Issue = issue,
    Expected = expected,
    Observed = observed,
    RowCount = suppressWarnings(as.integer(row_count)),
    Examples = examples,
    stringsAsFactors = FALSE
  )
}

.categorate_bind_quality = function(rows, cols) {
  rows = rows[!vapply(rows, is.null, logical(1))]
  rows = rows[vapply(rows, function(row) {
    is.data.frame(row) && nrow(row) > 0
  }, logical(1))]
  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = do.call(rbind, rows)
  for (col in setdiff(cols, colnames(out))) out[[col]] = NA
  row.names(out) = NULL
  out[, cols, drop = FALSE]
}

.categorate_validation_issue_cols = function() {
  c("Severity", "Table", "Column", "Issue", "Expected", "Observed",
    "RowCount", "Examples")
}

.categorate_table_quality_cols = function() {
  c("Table", "Present", "RowCount", "ColumnCount", "RequiredColumnCount",
    "MissingRequiredColumns", "MissingOptionalColumns", "Completeness",
    "DuplicateKeyColumns", "DuplicateKeyCount", "QueryCount", "CIDCount",
    "SourceDatabaseCount", "HighConfidenceRows", "MediumConfidenceRows",
    "LowConfidenceRows", "Status")
}

.categorate_source_diagnostic_cols = function() {
  c("Query", "Source", "SourceGroup", "Present", "RecordCount",
    "QualityStatus", "Diagnostic", "RetrievedAt")
}
