fixture_bioassay_detail_response = function(aids = c(9001, 9002, 9003)) {
  detail = function(aid, name, description, target_name, target_description,
                    accession, gene, taxid, organism, common,
                    endpoints = c("Potency", "Efficacy"),
                    include_target_organism = TRUE) {
    list(assay = list(descr = list(
      aid = list(id = aid, version = 1),
      aid_source = list(db = list(name = "Fixture BioAssay",
                                  source_id = list(str = paste0("FIX", aid)))),
      name = name,
      description = list(paste("Assay Overview:", description)),
      protocol = list("Fixture protocol with dose-response measurement."),
      comment = list("Keywords: qHTS, target metadata, bioactivity"),
      xref = list(
        list(xref = list(pmid = 123456)),
        list(xref = list(gene = gene)),
        list(xref = list(taxonomy = taxid), comment = organism),
        list(xref = list(dburl = "https://example.test/bioassay"))
      ),
      results = lapply(seq_along(endpoints), function(i) {
        list(tid = i, name = endpoints[[i]])
      }),
      target = list(c(
        list(
          name = target_name,
          descr = target_description,
          mol_id = list(protein_accession = accession)
        ),
        if (include_target_organism) {
          list(organism = list(org = list(
            taxname = organism,
            common = common,
            db = list(list(db = "taxon", tag = list(id = taxid)))
          )))
        } else {
          list()
        }
      )),
      activity_outcome_method = 2,
      project_category = 1
    )))
  }

  containers = list()
  if (9001 %in% aids) {
    containers[[length(containers) + 1]] = detail(
      9001,
      "qHTS Assay for Inhibitors of PTGS1 (Cyclooxygenase-1)",
      "PTGS1 inhibitor confirmation",
      "prostaglandin G/H synthase 1 [Homo sapiens]",
      "Human PTGS1 (COX-1)",
      "NP_000953.2",
      5743,
      9606,
      "Homo sapiens",
      "human",
      c("IC50", "Efficacy")
    )
  }
  if (9002 %in% aids) {
    containers[[length(containers) + 1]] = detail(
      9002,
      "qHTS Assay for Inhibitors of PTGS2 (Cyclooxygenase-2)",
      "PTGS2 inhibitor counter-screen",
      "prostaglandin G/H synthase 2 [Homo sapiens]",
      "Human PTGS2 (COX-2)",
      "NP_000954.1",
      5742,
      9606,
      "Homo sapiens",
      "human",
      c("Activity Outcome"),
      include_target_organism = FALSE
    )
  }
  if (9003 %in% aids) {
    containers[[length(containers) + 1]] = detail(
      9003,
      "qHTS Assay for Agonists of the Thyroid Stimulating Hormone Receptor",
      "TSHR agonist assay",
      "thyroid stimulating hormone receptor [Homo sapiens]",
      "Human TSHR",
      "AAR07906",
      7253,
      9606,
      "Homo sapiens",
      "human",
      c("Potency", "Response")
    )
  }
  list(PC_AssayContainer = containers)
}

fixture_pubchem_enrichment_request = function(url) {
  if (grepl("/pug/compound/name/aspirin/cids/JSON$", url)) {
    return(list(IdentifierList = list(CID = list(2244))))
  }

  if (grepl("/property/", url)) {
    return(list(PropertyTable = list(Properties = list(list(
      CID = 2244,
      Title = "Aspirin",
      MolecularFormula = "C9H8O4",
      MolecularWeight = "180.159",
      IUPACName = "2-acetyloxybenzoic acid",
      InChIKey = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
      CanonicalSMILES = "CC(=O)OC1=CC=CC=C1C(=O)O",
      IsomericSMILES = "CC(=O)OC1=CC=CC=C1C(=O)O",
      ExactMass = "180.04225873",
      XLogP = "1.2",
      TPSA = "63.6",
      HBondDonorCount = "1",
      HBondAcceptorCount = "4",
      RotatableBondCount = "3",
      HeavyAtomCount = "13"
    )))))
  }

  if (grepl("/synonyms/JSON$", url)) {
    return(list(InformationList = list(Information = list(list(
      CID = 2244,
      Synonym = list("aspirin", "50-78-2")
    )))))
  }

  if (grepl("heading=Names%20and%20Identifiers", url)) {
    return(list(Record = list(
      Reference = list(list(
        ReferenceNumber = 1,
        SourceName = "Fixture Names",
        URL = "https://example.test/names"
      )),
      Section = list(list(
        TOCHeading = "Names and Identifiers",
        Information = list(list(
          Name = "External IDs",
          ReferenceNumber = 1,
          Value = list(StringWithMarkup = list(list(
            String = "KEGG: C01405; PMID:123456"
          )))
        ))
      ))
    )))
  }

  if (grepl("heading=Safety%20and%20Hazards", url)) {
    return(list(Record = list(
      Reference = list(list(
        ReferenceNumber = 2,
        SourceName = "Fixture Safety",
        URL = "https://example.test/safety"
      )),
      Section = list(list(
        TOCHeading = "Safety and Hazards",
        Information = list(list(
          Name = "Hazard Statement",
          ReferenceNumber = 2,
          Value = list(StringWithMarkup = list(list(
            String = "H319 Warning Causes serious eye irritation; P264; LD50 oral rat 200 mg/kg"
          )))
        ))
      ))
    )))
  }

  if (grepl("source=.*Flavor|source=.*FEMA|Flavor|FEMA", url)) {
    return(list(Record = list(
      Reference = list(list(
        ReferenceNumber = 4,
        SourceName = "Flavor and Extract Manufacturers Association (FEMA)",
        URL = "https://example.test/fema"
      )),
      Section = list(list(
        TOCHeading = "FEMA Flavor Profile",
        Information = list(list(
          Name = "Flavor Profile",
          ReferenceNumber = 4,
          Value = list(StringWithMarkup = list(list(
            String = "FEMA No. 3009; GRAS; JECFA No. 75; citrus sweet odor and flavor"
          )))
        ))
      ))
    )))
  }

  if (grepl("source=.*FDA|source=.*SPL|FDA|SPL", url)) {
    return(list(Record = list(
      Reference = list(list(
        ReferenceNumber = 5,
        SourceName = "FDA/SPL Indexing Data",
        URL = "https://example.test/fda-spl"
      )),
      Section = list(list(
        TOCHeading = "FDA SPL",
        Information = list(
          list(
            Name = "Active Ingredient",
            ReferenceNumber = 5,
            Value = list(StringWithMarkup = list(list(String = "Aspirin")))
          ),
          list(
            Name = "Pharmacologic Class",
            ReferenceNumber = 5,
            Value = list(StringWithMarkup = list(list(String = "Cyclooxygenase Inhibitors; oral tablet")))
          )
        )
      ))
    )))
  }

  if (grepl("source=.*LOTUS|LOTUS", url)) {
    return(list(Record = list(
      Reference = list(list(
        ReferenceNumber = 6,
        SourceName = "LOTUS - the natural products occurrence database",
        URL = "https://example.test/lotus"
      )),
      Section = list(list(
        TOCHeading = "LOTUS Occurrence",
        Information = list(
          list(
            Name = "LOTUS ID",
            ReferenceNumber = 6,
            Value = list(StringWithMarkup = list(list(String = "LTS000001")))
          ),
          list(
            Name = "Organism",
            ReferenceNumber = 6,
            Value = list(StringWithMarkup = list(list(String = "Salvia officinalis")))
          ),
          list(
            Name = "Taxonomy",
            ReferenceNumber = 6,
            Value = list(StringWithMarkup = list(list(
              String = "Eukaryota; Plantae; Streptophyta; Magnoliopsida; Lamiales; Lamiaceae"
            )))
          ),
          list(
            Name = "Occurrence",
            ReferenceNumber = 6,
            Value = list(StringWithMarkup = list(list(
              String = "Reported in leaves and flowers; occurrence evidence is not a taxonomic lineage"
            )))
          ),
          list(
            Name = "Natural Product Class",
            ReferenceNumber = 6,
            Value = list(StringWithMarkup = list(list(String = "phenylpropanoid")))
          )
        )
      ))
    )))
  }

  if (grepl("source=.*Medical|source=.*MeSH|MeSH|Medical", url)) {
    return(list(Record = list(
      Reference = list(list(
        ReferenceNumber = 7,
        SourceName = "Medical Subject Headings (MeSH)",
        URL = "https://example.test/mesh"
      )),
      Section = list(list(
        TOCHeading = "MeSH Pharmacologic Action",
        Information = list(list(
          Name = "Pharmacologic Action",
          ReferenceNumber = 7,
          Value = list(StringWithMarkup = list(list(
            String = "Anti-Inflammatory Agents, Non-Steroidal; D27.505.519"
          )))
        ))
      ))
    )))
  }

  if (grepl("heading=Experimental%20Properties", url)) {
    return(list(Record = list(
      Reference = list(list(
        ReferenceNumber = 3,
        SourceName = "Fixture Experimental",
        URL = "https://example.test/experimental"
      )),
      Section = list(list(
        TOCHeading = "Experimental Properties",
        Information = list(
          list(
            Name = "Boiling Point",
            ReferenceNumber = 3,
            Value = list(StringWithMarkup = list(list(String = "140 C")))
          ),
          list(
            Name = "Density",
            ReferenceNumber = 3,
            Value = list(StringWithMarkup = list(list(String = "1.40 g/mL")))
          )
        )
      ))
    )))
  }

  if (grepl("/assaysummary/JSON$", url)) {
    cols = c("AID", "Panel Member ID", "SID", "CID",
             "Activity Outcome", "Target Accession", "Target GeneID",
             "Activity Value [uM]", "Activity Name", "Assay Name",
             "Assay Type", "PubMed ID", "RNAi")
    return(list(Table = list(
      Columns = list(Column = as.list(cols)),
      Row = list(
        list(Cell = as.list(c(
          "9001", "", "70001", "2244", "Active", "NP_000953", "5743",
          "0.35", "IC50",
          "qHTS Assay for Inhibitors of PTGS1 (Cyclooxygenase-1)",
          "Confirmatory", "123456", ""
        ))),
        list(Cell = as.list(c(
          "9002", "", "70002", "2244", "Inactive", "NP_000954", "5742",
          "", "",
          "qHTS Assay for Inhibitors of PTGS2 (Cyclooxygenase-2)",
          "Confirmatory", "", ""
        ))),
        list(Cell = as.list(c(
          "9003", "", "70003", "2244", "Inconclusive", "AAR07906", "",
          "12.5", "Potency",
          "qHTS Assay for Agonists of the Thyroid Stimulating Hormone Receptor",
          "Confirmatory", "", ""
        )))
      )
    )))
  }

  if (grepl("/assay/aid/", url) && grepl("/description/JSON$", url)) {
    aid_text = sub("^.*/assay/aid/([^/]+)/description/JSON$", "\\1", url)
    aids = suppressWarnings(as.integer(strsplit(aid_text, ",", fixed = TRUE)[[1]]))
    return(fixture_bioassay_detail_response(aids))
  }

  list(Record = list(Section = list()))
}

fixture_kegg_enrichment_request = function(url) {
  compound_record = paste(c(
    "ENTRY       C01405                      Compound",
    "NAME        Aspirin;",
    "            Acetylsalicylic acid",
    "FORMULA     C9H8O4",
    "EXACT_MASS  180.0423",
    "MOL_WEIGHT  180.159",
    "PATHWAY     map00590  Arachidonic acid metabolism",
    "REACTION    R01335",
    "ENZYME      3.1.1.55",
    "DBLINKS     PubChem: 2244",
    "            ChEBI: 15365"
  ), collapse = "\n")

  if (grepl("/find/compound/aspirin$", url)) {
    return("cpd:C01405\tAspirin; Acetylsalicylic acid")
  }
  if (grepl("/find/drug/aspirin$", url)) {
    return("")
  }
  if (grepl("/get/C01405$", url)) {
    return(paste(compound_record, "///", sep = "\n"))
  }
  if (grepl("/get/path:map00590$", url)) {
    return(paste(c(
      "ENTRY       map00590                    Pathway",
      "NAME        Arachidonic acid metabolism",
      "///"
    ), collapse = "\n"))
  }
  if (grepl("/get/rn:R01335$", url)) {
    return(paste(c(
      "ENTRY       R01335                      Reaction",
      "NAME        acetylsalicylate deacetylase reaction",
      "DEFINITION  Aspirin + H2O <=> Salicylate + Acetate",
      "EQUATION    C01405 + C00001 <=> C00805 + C00033",
      "///"
    ), collapse = "\n"))
  }
  if (grepl("/get/ec:3.1.1.55$", url)) {
    return(paste(c(
      "ENTRY       EC 3.1.1.55                 Enzyme",
      "NAME        acetylsalicylate deacetylase",
      "///"
    ), collapse = "\n"))
  }
  if (grepl("/link/pathway/cpd:C01405$", url)) {
    return("cpd:C01405\tpath:map00590")
  }
  if (grepl("/link/reaction/cpd:C01405$", url)) {
    return("cpd:C01405\trn:R01335")
  }
  if (grepl("/link/enzyme/cpd:C01405$", url)) {
    return("cpd:C01405\tec:3.1.1.55")
  }
  ""
}

fixture_taxonomy_record = function(taxid, scientific_name, common_name,
                                   lineage) {
  common_info = if (is.na(common_name)) {
    list()
  } else {
    list(list(
      ReferenceNumber = 1,
      Value = list(StringWithMarkup = list(list(String = common_name)))
    ))
  }
  list(Record = list(
    RecordType = "Taxonomy",
    RecordNumber = taxid,
    RecordTitle = scientific_name,
    Section = list(list(
      TOCHeading = "Taxonomy Information",
      Section = list(
        list(
          TOCHeading = "Scientific Name",
          Information = list(list(
            ReferenceNumber = 1,
            Value = list(StringWithMarkup = list(list(String = scientific_name)))
          ))
        ),
        list(
          TOCHeading = "Common Name",
          Information = common_info
        ),
        list(
          TOCHeading = "Rank",
          Information = list(list(
            ReferenceNumber = 1,
            Value = list(StringWithMarkup = list(list(String = "species")))
          ))
        ),
        list(
          TOCHeading = "Domain",
          Information = list(list(
            ReferenceNumber = 1,
            Value = list(StringWithMarkup = list(list(String = "Eukaryota")))
          ))
        ),
        list(
          TOCHeading = "Lineage",
          Information = list(list(
            ReferenceNumber = 1,
            Value = list(StringWithMarkup = list(list(String = lineage)))
          ))
        )
      )
    ))
  ))
}

fixture_classification_response = function(hid, tree_name, source_id, paths) {
  list(Hierarchies = list(Hierarchy = lapply(seq_along(paths), function(path_index) {
    terms = paths[[path_index]]
    leaf_to_root = rev(terms)
    nodes = lapply(seq_along(leaf_to_root), function(node_index) {
      term = leaf_to_root[[node_index]]
      node_id = paste0("N", hid, "_", path_index, "_", node_index)
      parent_id = if (node_index < length(leaf_to_root)) {
        paste0("N", hid, "_", path_index, "_", node_index + 1L)
      } else {
        NA_character_
      }
      list(
        NodeID = node_id,
        ParentID = parent_id,
        Information = list(
          HNID = paste0("HNID", hid, "_", path_index, "_", node_index),
          Name = list(StringWithMarkup = list(list(String = term))),
          Match = node_index == 1L,
          Counts = list(
            list(Type = "Compound", Count = 10L + path_index),
            list(Type = "Taxonomy", Count = 2L + path_index),
            list(Type = "DOI", Count = path_index),
            list(Type = "PubMed", Count = path_index + 1L)
          )
        )
      )
    })
    list(
      SourceName = "LOTUS - the natural products occurrence database",
      SourceID = source_id,
      HID = hid,
      Information = list(
        HID = hid,
        HNID = paste0("ROOT", hid),
        Name = list(StringWithMarkup = list(list(String = tree_name)))
      ),
      Node = nodes
    )
  })))
}

fixture_lotus_taxonomy_request = function(url) {
  if (grepl("/pug/compound/name/methyl%20salicylate/cids/JSON$", url)) {
    return(list(IdentifierList = list(CID = list(4133))))
  }
  if (grepl("/property/", url)) {
    return(list(PropertyTable = list(Properties = list(list(
      CID = 4133,
      Title = "Methyl Salicylate"
    )))))
  }
  if (grepl("/synonyms/JSON$", url)) {
    return(list(InformationList = list(Information = list(list(
      CID = 4133,
      Synonym = list("methyl salicylate")
    )))))
  }
  if (grepl("/pug_view/data/compound/4133/JSON\\?source=.*LOTUS", url)) {
    description = "Methyl Salicylate has been reported in Camellia sinensis, Phellinus tremulae, and other organisms with data available."
    return(list(Record = list(
      Reference = list(list(
        ReferenceNumber = 72,
        SourceName = "LOTUS - the natural products occurrence database",
        URL = "https://lotus.nprod.net/"
      )),
      Section = list(
        list(
          TOCHeading = "Names and Identifiers",
          Section = list(list(
            TOCHeading = "Record Description",
            Information = list(list(
              ReferenceNumber = 72,
              Value = list(StringWithMarkup = list(list(
                String = description,
                Markup = list(
                  list(Start = 39, Length = 17,
                       URL = "https://pubchem.ncbi.nlm.nih.gov/taxonomy/4442#section=Natural-Products"),
                  list(Start = 58, Length = 18,
                       URL = "https://pubchem.ncbi.nlm.nih.gov/taxonomy/108899#section=Natural-Products"),
                  list(Start = 82, Length = 15,
                       URL = "#section=Taxonomy")
                )
              )))
            ))
          ))
        ),
        list(
          TOCHeading = "Taxonomy",
          Information = list(list(
            ReferenceNumber = 72,
            Value = list(ExternalTableName = "consolidatedcompoundtaxonomy")
          ))
        ),
        list(
          TOCHeading = "Classification",
          Section = list(
            list(
              TOCHeading = "LOTUS: Biological Tree",
              Information = list(list(
                ReferenceNumber = 72,
                Name = "HID",
                Value = list(Number = list(115))
              ))
            ),
            list(
              TOCHeading = "LOTUS: Chemical Tree",
              Information = list(list(
                ReferenceNumber = 72,
                Name = "HID",
                Value = list(Number = list(142))
              ))
            )
          )
        )
      )
    )))
  }
  if (grepl("/pug_view/data/taxonomy/4442/JSON$", url)) {
    return(fixture_taxonomy_record(
      4442,
      "Camellia sinensis",
      "black tea",
      "Eukaryota; Viridiplantae; Streptophyta; Magnoliopsida; Ericales; Theaceae; Camellia"
    ))
  }
  if (grepl("/pug_view/data/taxonomy/108899/JSON$", url)) {
    return(fixture_taxonomy_record(
      108899,
      "Phellinus tremulae",
      NA_character_,
      "Eukaryota; Fungi; Basidiomycota; Agaricomycetes; Hymenochaetales; Hymenochaetaceae; Phellinus"
    ))
  }
  if (grepl("classifications\\.fcgi", url) &&
      grepl("hid=115", url) &&
      grepl("search_uid=4133", url)) {
    return(fixture_classification_response(
      hid = 115,
      tree_name = "LOTUS: Biological Tree",
      source_id = "LOTUS",
      paths = list(
        c("Cytota", "Eukaryota", "Archaeplastida", "Plantae",
          "Viridiplantae", "Streptophyta", "Magnoliopsida", "Ericales",
          "Theaceae", "Camellia", "Camellia sinensis"),
        c("Cytota", "Eukaryota", "Opisthokonta", "Holomycota", "Fungi",
          "Basidiomycota", "Agaricomycetes", "Hymenochaetales",
          "Hymenochaetaceae", "Phellinus", "Phellinus tremulae"),
        c("Cytota", "Eukaryota", "Archaeplastida", "Plantae",
          "Viridiplantae", "Streptophyta", "Magnoliopsida", "Lamiales",
          "Lamiaceae", "Salvia", "Salvia officinalis")
      )
    ))
  }
  if (grepl("classifications\\.fcgi", url) &&
      grepl("hid=142", url) &&
      grepl("search_uid=4133", url)) {
    return(fixture_classification_response(
      hid = 142,
      tree_name = "LOTUS: Chemical Tree",
      source_id = "LOTUS",
      paths = list(
        c("Shikimates and Phenylpropanoids",
          "Phenolic acids (C6-C1)",
          "Simple phenolic acids")
      )
    ))
  }
  list(Record = list(Section = list()))
}

expect_no_duplicate_keys = function(table, cols) {
  if (!is.data.frame(table) || nrow(table) < 2) return(invisible(TRUE))
  keys = table[, cols, drop = FALSE]
  expect_false(any(duplicated(keys)))
}

test_that("categorate research enrichment builds analysis-ready tables", {
  data_list = list(
    reactives = data.frame(reactives = "Carboxylic acids",
                           Chemical = "aspirin"),
    LOTUS = data.frame(LOTUS = "None", Chemical = "aspirin"),
    KEGG = data.frame(KEGG = "KEGG: C01405", Chemical = "aspirin"),
    FEMA = data.frame(FEMA = "None", Chemical = "aspirin"),
    FDA_SPL = data.frame(FDA_SPL = "Aspirin tablet", Chemical = "aspirin")
  )

  enrichment = .categorate_research_enrichment(
    compounds = "aspirin",
    data_list = data_list,
    detail = "research",
    cache = FALSE,
    cache_dir = NULL,
    throttle = 0,
    request_fun = fixture_pubchem_enrichment_request,
    kegg_request_fun = fixture_kegg_enrichment_request
  )

  expect_true(all(c("PubChemProperties", "PubChemSourceAnnotations",
                    "SafetyProfile", "FEMAProfile", "FDA_SPL_Profile",
                    "LOTUSProfile", "PubChemClassifications",
                    "MeSHProfile", "LiteratureProfile",
                    "ChemicalTerms", "ChemicalTraits",
                    "ChemicalTraitOntology", "ChemicalTraitMatrix",
                    "ChemicalTraitOntologyMatrix", "ChemicalTraitEvidence",
                    "ChemicalTraitReport", "ChemicalTraitSummary",
                    "ChemicalTraitSimilarity", "ChemicalClasses",
                    "ChemicalMeasurements", "ChemicalMeasurementSummary",
                    "ChemicalHazards", "ChemicalUses", "ChemicalBioassays",
                    "ChemicalBioactivities", "ChemicalTargets",
                    "ChemicalPotencies", "PubChemBioAssayDetails",
                    "ChemicalTaxonomy",
                    "ChemicalOccurrences",
                    "ChemicalPathwayRoles", "KEGGReactionParticipants",
                    "KEGGPathways", "KEGGReactions", "KEGGEnzymes",
                    "SourceCoverage", "DerivedGroups", "DataDictionary",
                    "TableQuality", "SourceDiagnostics",
                    "ValidationIssues", "ValidationSummary") %in%
                    names(enrichment)))
  expect_true(any(enrichment$PubChemIdentifiers$IdentifierType == "PubMed"))
  expect_true(any(enrichment$SafetyProfile$HazardCode == "H319",
                  na.rm = TRUE))
  expect_true(any(enrichment$FEMAProfile$FEMANumber == "3009",
                  na.rm = TRUE))
  expect_true(any(grepl("citrus", enrichment$FEMAProfile$DescriptorTerms),
                  na.rm = TRUE))
  expect_true(any(enrichment$FDA_SPL_Profile$ActiveIngredient == "Aspirin",
                  na.rm = TRUE))
  expect_true(any(grepl("oral", enrichment$FDA_SPL_Profile$Route),
                  na.rm = TRUE))
  expect_true(any(grepl("LTS000001", enrichment$LOTUSProfile$LOTUS_ID),
                  na.rm = TRUE))
  expect_true(any(grepl("phenylpropanoid", enrichment$LOTUSProfile$NaturalProductClass),
                  na.rm = TRUE))
  expect_true(any(grepl("Anti-Inflammatory", enrichment$MeSHProfile$PharmacologicAction),
                  na.rm = TRUE))
  expect_true(any(enrichment$KEGGPathways$PathwayGroup == "Lipid metabolism",
                  na.rm = TRUE))
  expect_true(any(enrichment$KEGGReactions$ReactionName == "acetylsalicylate deacetylase reaction",
                  na.rm = TRUE))
  expect_true(any(enrichment$KEGGReactions$Equation == "C01405 + C00001 <=> C00805 + C00033",
                  na.rm = TRUE))
  expect_true(any(enrichment$KEGGEnzymes$EnzymeClass == "Hydrolases",
                  na.rm = TRUE))
  expect_true(any(enrichment$KEGGLinkMetadata$TargetDatabase == "reaction"))
  expect_no_duplicate_keys(enrichment$PubChemAnnotations,
                           c("Query", "CID", "Name", "CleanValue",
                             "ValueNumeric", "UnitClean", "Source"))
  expect_no_duplicate_keys(enrichment$FEMAProfile,
                           c("Query", "CID", "FEMANumber", "GRASStatus",
                             "JECFANumber", "DescriptorTerms", "FlavorTerms",
                             "OdorTerms"))
  expect_no_duplicate_keys(enrichment$FDA_SPL_Profile,
                           c("Query", "CID", "ActiveIngredient",
                             "PharmacologicClass", "Route", "DosageForm",
                             "HasBoxedWarning", "WarningTerms"))
  expect_no_duplicate_keys(enrichment$LOTUSProfile,
                           c("Query", "CID", "LOTUS_ID", "Organism",
                             "Taxonomy", "NaturalProductClass", "Occurrence",
                             "ReferenceID"))
  expect_true(any(enrichment$ChemicalTerms$TermClean == "citrus" &
                    enrichment$ChemicalTerms$TermGroup == "citrus"))
  expect_true(any(enrichment$ChemicalTerms$TermClean == "anti-inflammatory" &
                    enrichment$ChemicalTerms$TermGroup == "anti-inflammatory"))
  expect_true(any(enrichment$ChemicalTerms$TermClean == "cyclooxygenase inhibitor" &
                    enrichment$ChemicalTerms$TermGroup == "enzyme inhibitor"))
  expect_true(any(enrichment$ChemicalTraits$TraitType == "hazard" &
                    enrichment$ChemicalTraits$TraitGroup == "hazard_code" &
                    enrichment$ChemicalTraits$TraitValue == "H319"))
  expect_true(any(enrichment$ChemicalTraits$TraitType == "pathway" &
                    enrichment$ChemicalTraits$TraitGroup == "pathway_group" &
                    enrichment$ChemicalTraits$TraitValue == "Lipid metabolism"))
  expect_true(any(enrichment$ChemicalTraits$TraitType == "occurrence" &
                    enrichment$ChemicalTraits$TraitGroup == "family" &
                    enrichment$ChemicalTraits$TraitValue == "Lamiaceae"))
  expect_true(any(enrichment$ChemicalTraits$TraitType == "property" &
                    enrichment$ChemicalTraits$TraitGroup == "xlogp_bin" &
                    enrichment$ChemicalTraits$TraitValue == "balanced_XLogP_0_3"))
  expect_no_duplicate_keys(enrichment$ChemicalTraits,
                           c("Query", "CID", "TraitType", "TraitGroup",
                             "TraitValueClean"))
  expect_true(any(enrichment$ChemicalTraitOntology$OntologyDomain == "safety" &
                    enrichment$ChemicalTraitOntology$OntologyGroup == "ghs_hazard_code" &
                    enrichment$ChemicalTraitOntology$OntologyTerm == "h319" &
                    enrichment$ChemicalTraitOntology$OntologyID == "H319" &
                    enrichment$ChemicalTraitOntology$OntologyIDSource == "GHS" &
                    enrichment$ChemicalTraitOntology$SourceTraitKey ==
                      "hazard__hazard_code__h319"))
  expect_true(any(enrichment$ChemicalTraitOntology$OntologyDomain == "safety" &
                    enrichment$ChemicalTraitOntology$OntologyGroup == "hazard_endpoint" &
                    enrichment$ChemicalTraitOntology$OntologyTerm == "eye damage irritation"))
  expect_true(any(enrichment$ChemicalTraitOntology$OntologyDomain == "metabolism" &
                    enrichment$ChemicalTraitOntology$OntologyGroup == "pathway_group" &
                    enrichment$ChemicalTraitOntology$OntologyTerm == "lipid metabolism" &
                    enrichment$ChemicalTraitOntology$OntologyID == "map00590" &
                    enrichment$ChemicalTraitOntology$OntologyIDSource == "KEGG PATHWAY" &
                    enrichment$ChemicalTraitOntology$OntologyIDURL ==
                      "https://www.kegg.jp/entry/map00590"))
  expect_true(any(enrichment$ChemicalTraitOntology$OntologyDomain == "ecology" &
                    enrichment$ChemicalTraitOntology$OntologyGroup == "taxonomic_family" &
                    enrichment$ChemicalTraitOntology$OntologyTerm == "lamiaceae"))
  expect_true(any(enrichment$ChemicalTraitOntology$OntologyDomain == "physicochemical" &
                    enrichment$ChemicalTraitOntology$OntologyGroup == "lipophilicity" &
                    enrichment$ChemicalTraitOntology$OntologyTerm == "balanced xlogp 0 3"))
  expect_true(any(enrichment$ChemicalTraitOntology$OntologyDomain == "biomedical" &
                    enrichment$ChemicalTraitOntology$OntologyGroup == "mesh_tree_category" &
                    enrichment$ChemicalTraitOntology$OntologyID == "D27.505.519" &
                    enrichment$ChemicalTraitOntology$OntologyIDSource ==
                      "MeSH tree number"))
  expect_true(any(enrichment$ChemicalTraitOntology$OntologyDomain == "metabolism" &
                    enrichment$ChemicalTraitOntology$OntologyGroup == "ec_number" &
                    enrichment$ChemicalTraitOntology$OntologyID == "ec:3.1.1.55" &
                    enrichment$ChemicalTraitOntology$OntologyIDURL ==
                      "https://www.kegg.jp/entry/ec:3.1.1.55"))
  expect_true(any(enrichment$ChemicalTraitOntology$OntologyDomain == "metabolism" &
                    enrichment$ChemicalTraitOntology$OntologyGroup ==
                      "reaction_participant_id" &
                    enrichment$ChemicalTraitOntology$OntologyTerm == "c01405" &
                    enrichment$ChemicalTraitOntology$OntologyID == "C01405" &
                    enrichment$ChemicalTraitOntology$OntologyIDSource ==
                      "KEGG COMPOUND"))
  expect_no_duplicate_keys(enrichment$ChemicalTraitOntology,
                           c("Query", "CID", "OntologyDomain",
                             "OntologyGroup", "OntologyTerm"))
  expect_true(any(enrichment$ChemicalTraitEvidence$EvidenceType == "ontology" &
                    enrichment$ChemicalTraitEvidence$AnalysisKey ==
                      "safety__ghs_hazard_code__h319" &
                    enrichment$ChemicalTraitEvidence$ExternalID == "H319" &
                    enrichment$ChemicalTraitEvidence$SourceTraitKey ==
                      "hazard__hazard_code__h319" &
                    enrichment$ChemicalTraitEvidence$EvidenceURL ==
                      "https://example.test/safety" &
                    enrichment$ChemicalTraitEvidence$ExtractionRule ==
                      "ghs_and_toxicity_regex"))
  h319_evidence = chemicalTraitEvidence(
    enrichment,
    keys = "safety__ghs_hazard_code__h319"
  )
  expect_equal(nrow(h319_evidence), 1)
  expect_equal(h319_evidence$ExternalID, "H319")
  expect_equal(h319_evidence$SourceTraitKey, "hazard__hazard_code__h319")
  expect_equal(h319_evidence$EvidenceURL, "https://example.test/safety")
  kegg_evidence = chemicalTraitEvidence(enrichment, keys = "map00590")
  expect_true(any(kegg_evidence$AnalysisKey ==
                    "metabolism__pathway_group__lipid_metabolism" &
                    kegg_evidence$ExternalID == "map00590" &
                    grepl("rest.kegg.jp", kegg_evidence$EvidenceURL,
                          fixed = TRUE)))
  trait_evidence = chemicalTraitEvidence(
    enrichment,
    keys = "hazard__hazard_code__h319",
    type = "trait"
  )
  expect_equal(nrow(trait_evidence), 1)
  expect_equal(trait_evidence$EvidenceType, "trait")
  expect_equal(trait_evidence$AnalysisKey, "hazard__hazard_code__h319")
  expect_equal(trait_evidence$EvidenceID, "H319")
  report = enrichment$ChemicalTraitReport
  expect_equal(nrow(report), 1)
  expect_equal(report$Query, "aspirin")
  expect_gt(report$TraitCount, 0)
  expect_gt(report$OntologyTermCount, 0)
  expect_gt(report$EvidenceCount, 0)
  expect_true(grepl("safety", report$OntologyDomains))
  expect_true(grepl("ghs_hazard_code=H319", report$SafetySignals,
                    fixed = TRUE))
  expect_true(grepl("lipophilicity=balanced_XLogP_0_3",
                    report$PhysicochemicalSignals,
                    fixed = TRUE))
  expect_true(grepl("pathway_group=Lipid metabolism",
                    report$MetabolismSignals,
                    fixed = TRUE))
  expect_true(grepl("taxonomic_family=Lamiaceae", report$EcologySignals,
                    fixed = TRUE))
  expect_true(grepl("GHS:H319", report$ExternalIDs, fixed = TRUE))
  expect_true(grepl("https://example.test/safety", report$EvidenceURLs,
                    fixed = TRUE))
  expect_true(is.na(report$NearestNeighbors))
  report_from_helper = chemicalTraitReport(enrichment)
  expect_equal(report_from_helper$SafetySignals, report$SafetySignals)
  dictionary = chemicalDataDictionary("ChemicalTraitReport")
  expect_true(all(c("SafetySignals", "MetabolismSignals",
                    "NearestNeighbors") %in% dictionary$Column))
  expect_true(any(dictionary$Column == "TraitCount" &
                    dictionary$Type == "integer" &
                    dictionary$Required))
  measurement_dictionary = chemicalDataDictionary("ChemicalMeasurements")
  expect_true(any(measurement_dictionary$Column == "StandardValue" &
                    measurement_dictionary$Type == "numeric"))
  expect_true(any(measurement_dictionary$Column == "BehaviorBin"))
  validation = validateCategorateResult(enrichment)
  expect_s3_class(validation, "uaf_categorate_validation")
  expect_equal(validation$Summary$Status, "pass")
  expect_equal(nrow(validation$Issues), 0)
  expect_equal(enrichment$ValidationSummary$Status, "pass")
  expect_equal(nrow(enrichment$ValidationIssues), 0)
  quality = enrichment$TableQuality
  report_quality = quality[quality$Table == "ChemicalTraitReport", ,
                           drop = FALSE]
  expect_equal(report_quality$RowCount, 1)
  expect_equal(report_quality$Completeness, 1)
  expect_equal(report_quality$DuplicateKeyCount, 0)
  expect_equal(report_quality$Status, "ok")
  source_diagnostics = enrichment$SourceDiagnostics
  expect_true(any(source_diagnostics$Query == "aspirin" &
                    source_diagnostics$Source == "KEGGReactions" &
                    source_diagnostics$QualityStatus == "present"))
  broken = enrichment
  broken$ChemicalTraitReport$Query = NULL
  broken_validation = validateCategorateResult(broken)
  expect_equal(broken_validation$Summary$Status, "fail")
  expect_true(any(broken_validation$Issues$Table == "ChemicalTraitReport" &
                    broken_validation$Issues$Column == "Query" &
                    broken_validation$Issues$Severity == "error"))
  export_dir = file.path(tempdir(), "uafR_categorate_export_fixture")
  manifest = exportCategorateWorkbook(enrichment,
                                      export_dir,
                                      format = "csv",
                                      overwrite = TRUE)
  expect_true(dir.exists(export_dir))
  expect_equal(manifest$Table[[1]], "ExportManifest")
  expect_true(all(c("ChemicalTraitReport", "ChemicalTraitEvidence",
                    "ChemicalMeasurementSummary", "DataDictionary",
                    "ValidationSummary") %in%
                    manifest$Table))
  report_file = file.path(export_dir,
                          manifest$FileName[
                            manifest$Table == "ChemicalTraitReport"
                          ])
  evidence_file = file.path(export_dir,
                            manifest$FileName[
                              manifest$Table == "ChemicalTraitEvidence"
                            ])
  expect_true(file.exists(report_file))
  expect_true(file.exists(evidence_file))
  exported_report = utils::read.csv(report_file,
                                    stringsAsFactors = FALSE,
                                    check.names = FALSE)
  expect_equal(exported_report$Query, "aspirin")
  expect_true(grepl("ghs_hazard_code=H319",
                    exported_report$SafetySignals,
                    fixed = TRUE))
  subset_dir = file.path(tempdir(), "uafR_categorate_export_subset")
  subset_manifest = NULL
  expect_warning(
    subset_manifest <- exportCategorateWorkbook(
      enrichment,
      subset_dir,
      tables = c("ChemicalTraitReport", "MissingTable"),
      format = "csv",
      overwrite = TRUE
    ),
    "Requested table"
  )
  expect_equal(subset_manifest$Table,
               c("ExportManifest", "ChemicalTraitReport"))
  expect_true(any(enrichment$ChemicalTraitSummary$TraitType == "hazard" &
                    enrichment$ChemicalTraitSummary$TraitCount >= 1))
  expect_true(any(enrichment$ChemicalTraitSummary$TraitType == "pathway" &
                    enrichment$ChemicalTraitSummary$SourceDatabases == "KEGG"))
  expect_equal(nrow(enrichment$ChemicalTraitSimilarity), 0)
  trait_matrix = enrichment$ChemicalTraitMatrix
  expect_equal(nrow(trait_matrix), 1)
  expect_true(all(c("hazard__hazard_code__h319",
                    "pathway__pathway_group__lipid_metabolism",
                    "occurrence__family__lamiaceae",
                    "property__xlogp_bin__balanced_xlogp_0_3") %in%
                    colnames(trait_matrix)))
  expect_equal(trait_matrix$hazard__hazard_code__h319[[1]], 1)
  expect_equal(trait_matrix$pathway__pathway_group__lipid_metabolism[[1]], 1)
  expect_equal(trait_matrix$occurrence__family__lamiaceae[[1]], 1)
  expect_equal(trait_matrix$property__xlogp_bin__balanced_xlogp_0_3[[1]], 1)
  ontology_matrix = enrichment$ChemicalTraitOntologyMatrix
  expect_equal(nrow(ontology_matrix), 1)
  expect_true(all(c("safety__ghs_hazard_code__h319",
                    "metabolism__pathway_group__lipid_metabolism",
                    "ecology__taxonomic_family__lamiaceae",
                    "physicochemical__lipophilicity__balanced_xlogp_0_3") %in%
                    colnames(ontology_matrix)))
  expect_equal(ontology_matrix$safety__ghs_hazard_code__h319[[1]], 1)
  expect_equal(ontology_matrix$metabolism__pathway_group__lipid_metabolism[[1]], 1)
  expect_equal(ontology_matrix$ecology__taxonomic_family__lamiaceae[[1]], 1)
  expect_equal(ontology_matrix$physicochemical__lipophilicity__balanced_xlogp_0_3[[1]], 1)
  full_trait_matrix = chemicalTraitMatrix(enrichment, profile = "full")
  expect_gt(ncol(full_trait_matrix), ncol(trait_matrix))
  expect_true("use__active_ingredient__aspirin" %in%
                colnames(full_trait_matrix))
  limited_trait_matrix = chemicalTraitMatrix(enrichment$ChemicalTraits,
                                             profile = "full",
                                             max_traits = 3)
  expect_equal(ncol(limited_trait_matrix), 5)
  confidence_trait_matrix = chemicalTraitMatrix(enrichment,
                                                profile = "safety",
                                                mode = "confidence",
                                                min_confidence = "medium")
  expect_true("hazard__hazard_code__h319" %in%
                colnames(confidence_trait_matrix))
  expect_true(confidence_trait_matrix$hazard__hazard_code__h319[[1]] >= 0.65)
  sensory_trait_matrix = chemicalTraitMatrix(enrichment, profile = "sensory")
  expect_true("use__flavor__sweet" %in% colnames(sensory_trait_matrix))
  expect_false(any(grepl("target__", colnames(sensory_trait_matrix))))
  expect_true(any(enrichment$ChemicalClasses$ClassSystem == "MeSH" &
                    enrichment$ChemicalClasses$ClassGroup == "Chemicals and Drugs"))
  expect_true(any(enrichment$ChemicalMeasurements$Property == "boiling_point" &
                    enrichment$ChemicalMeasurements$Value == 140 &
                    enrichment$ChemicalMeasurements$ValueRelation == "equal" &
                    enrichment$ChemicalMeasurements$StandardValue == 140 &
                    enrichment$ChemicalMeasurements$StandardUnit == "C" &
                    enrichment$ChemicalMeasurements$MeasurementClass == "temperature" &
                    enrichment$ChemicalMeasurements$BehaviorBin ==
                      "volatile_boiling_point_100_200_C"))
  expect_true(any(enrichment$ChemicalMeasurements$Property == "density" &
                    enrichment$ChemicalMeasurements$Value == 1.40 &
                    enrichment$ChemicalMeasurements$Unit == "g/mL" &
                    enrichment$ChemicalMeasurements$StandardValue == 1.40 &
                    enrichment$ChemicalMeasurements$StandardUnit == "g/mL" &
                    enrichment$ChemicalMeasurements$MeasurementClass == "density" &
                    enrichment$ChemicalMeasurements$BehaviorBin ==
                      "high_density_gt_1_2_g_mL"))
  expect_true(any(enrichment$ChemicalMeasurements$Property == "ld50" &
                    enrichment$ChemicalMeasurements$Value == 200 &
                    enrichment$ChemicalMeasurements$Unit == "mg/kg" &
                    enrichment$ChemicalMeasurements$ValueRelation == "equal" &
                    enrichment$ChemicalMeasurements$StandardValue == 200 &
                    enrichment$ChemicalMeasurements$StandardUnit == "mg/kg" &
                    enrichment$ChemicalMeasurements$MeasurementClass == "toxicity" &
                    enrichment$ChemicalMeasurements$BehaviorBin ==
                      "high_acute_toxicity_50_300_mg_kg"))
  expect_equal(.normalized_standard_measurement(
    property = "vapor_pressure",
    value = 1e-5,
    unit = "atm",
    evidence_text = "1e-5 atm"
  )$ValueRelation, "equal")
  measurement_summary = enrichment$ChemicalMeasurementSummary
  expect_true(any(measurement_summary$Property == "boiling_point" &
                    measurement_summary$MeasurementClass == "temperature" &
                    measurement_summary$MeasurementCount >= 1 &
                    measurement_summary$StandardUnit == "C" &
                    measurement_summary$BestValue == 140 &
                    grepl("volatile_boiling_point_100_200_C",
                          measurement_summary$BehaviorBins,
                          fixed = TRUE)))
  expect_true(any(measurement_summary$Property == "density" &
                    measurement_summary$MeasurementClass == "density" &
                    measurement_summary$StandardUnit == "g/mL" &
                    measurement_summary$BestValue == 1.40 &
                    grepl("high_density_gt_1_2_g_mL",
                          measurement_summary$BehaviorBins,
                          fixed = TRUE)))
  expect_true(any(measurement_summary$Property == "ld50" &
                    measurement_summary$MeasurementClass == "toxicity" &
                    measurement_summary$StandardUnit == "mg/kg" &
                    measurement_summary$BestValue == 200 &
                    measurement_summary$BestValueRelation == "equal" &
                    grepl("high_acute_toxicity_50_300_mg_kg",
                          measurement_summary$BehaviorBins,
                          fixed = TRUE)))
  helper_summary = chemicalMeasurementSummary(enrichment)
  expect_equal(nrow(helper_summary), nrow(measurement_summary))
  expect_equal(sort(helper_summary$Property), sort(measurement_summary$Property))
  expect_true(any(enrichment$ChemicalHazards$HazardCode == "H319" &
                    enrichment$ChemicalHazards$HazardGroup == "eye_damage_irritation" &
                    enrichment$ChemicalHazards$PrecautionCode == "P264" &
                    enrichment$ChemicalHazards$ExposureRoute == "oral" &
                    enrichment$ChemicalHazards$TargetOrgan == "eye" &
                    enrichment$ChemicalHazards$ToxicityMetric == "LD50" &
                    enrichment$ChemicalHazards$ToxicityValue == 200 &
                    enrichment$ChemicalHazards$Species == "rat",
                  na.rm = TRUE))
  expect_true(any(enrichment$ChemicalUses$UseType == "flavor" &
                    enrichment$ChemicalUses$UseTerm == "sweet" &
                    enrichment$ChemicalUses$UseGroup == "sweet"))
  expect_true(any(enrichment$ChemicalUses$UseType == "administration_route" &
                    enrichment$ChemicalUses$UseTerm == "oral" &
                    enrichment$ChemicalUses$UseGroup == "enteral route"))
  expect_false(any(enrichment$ChemicalUses$UseType == "pharmacologic_class" &
                     enrichment$ChemicalUses$UseTerm == "oral tablet"))
  expect_true(any(enrichment$ChemicalTaxonomy$Kingdom == "Plantae" &
                    enrichment$ChemicalTaxonomy$Phylum == "Streptophyta" &
                    enrichment$ChemicalTaxonomy$Class == "Magnoliopsida" &
                    enrichment$ChemicalTaxonomy$Order == "Lamiales" &
                    enrichment$ChemicalTaxonomy$Family == "Lamiaceae" &
                    enrichment$ChemicalTaxonomy$Genus == "Salvia" &
                    enrichment$ChemicalTaxonomy$Species == "Salvia officinalis"))
  expect_true(any(enrichment$ChemicalOccurrences$Organism == "Salvia officinalis" &
                    enrichment$ChemicalOccurrences$OccurrenceType == "natural_product" &
                    enrichment$ChemicalOccurrences$Kingdom == "Plantae" &
                    enrichment$ChemicalOccurrences$Family == "Lamiaceae"))
  expect_false(any(grepl("reported|occurrence evidence|database|lotus",
                         enrichment$ChemicalTaxonomy$TaxonomyTerms,
                         ignore.case = TRUE),
                   na.rm = TRUE))
  expect_false(any(grepl("reported|occurrence evidence|database|lotus",
                         enrichment$ChemicalTerms$TermClean[
                           enrichment$ChemicalTerms$TermType == "taxonomy"
                         ],
                         ignore.case = TRUE),
                   na.rm = TRUE))
  expect_true(any(enrichment$ChemicalPathwayRoles$RoleType == "enzyme_associated" &
                    enrichment$ChemicalPathwayRoles$ECNumber == "3.1.1.55" &
                    enrichment$ChemicalPathwayRoles$EnzymeClass == "Hydrolases"))
  expect_true(any(enrichment$KEGGReactionParticipants$ParticipantID == "C01405" &
                    enrichment$KEGGReactionParticipants$Side == "substrate" &
                    enrichment$KEGGReactionParticipants$ParticipantRole == "query_substrate"))
  expect_equal(nrow(enrichment$KEGGReactionParticipants), 4)
  expect_no_duplicate_keys(enrichment$KEGGReactionParticipants,
                           c("Query", "KEGG_ID", "ReactionID", "Side",
                             "ParticipantID", "ParticipantRole"))

  derived = enrichment$DerivedGroups
  expect_equal(derived$MolecularFormula, "C9H8O4")
  expect_true(derived$is_drug_or_label_ingredient)
  expect_true(derived$has_kegg_pathway)
  expect_true(derived$has_safety_hazard)
  expect_true(derived$is_natural_product)
  expect_equal(derived$occurrence_count, 1)
  expect_equal(derived$organism_count, 1)
  expect_equal(derived$dominant_kingdom, "Plantae")
  expect_equal(derived$dominant_family, "Lamiaceae")
  expect_equal(derived$taxonomic_breadth, "single_organism")
  expect_true(derived$is_plant_occurring)
  expect_false(derived$is_fungal_occurring)
  expect_true(derived$is_flavor_ingredient)
  expect_equal(derived$lipophilicity_bin, "balanced_XLogP_0_3")
  expect_equal(derived$polarity_bin, "moderate_TPSA_40_90")
  expect_true(grepl("Lipid metabolism", derived$metabolic_context))
  expect_true(grepl("drug_or_label_ingredient", derived$biomedical_context))
  expect_true(grepl("natural_product", derived$ecological_context))
  expect_true(grepl("flavor_ingredient", derived$sensory_context))
  expect_true(grepl("H319", derived$safety_context))
  expect_true(grepl("moderate_TPSA_40_90", derived$analytical_context))
  expect_true(derived$oxygenated)
  expect_gt(derived$confidence_score, 0.6)
  expect_equal(derived$evidence_score, derived$confidence_score)
})

test_that("categorate full enrichment normalizes PubChem bioactivity", {
  data_list = list(
    reactives = data.frame(reactives = "Carboxylic acids",
                           Chemical = "aspirin"),
    LOTUS = data.frame(LOTUS = "None", Chemical = "aspirin"),
    KEGG = data.frame(KEGG = "KEGG: C01405", Chemical = "aspirin"),
    FEMA = data.frame(FEMA = "None", Chemical = "aspirin"),
    FDA_SPL = data.frame(FDA_SPL = "Aspirin tablet", Chemical = "aspirin")
  )

  enrichment = .categorate_research_enrichment(
    compounds = "aspirin",
    data_list = data_list,
    detail = "full",
    cache = FALSE,
    cache_dir = NULL,
    throttle = 0,
    request_fun = fixture_pubchem_enrichment_request,
    kegg_request_fun = fixture_kegg_enrichment_request
  )

  expect_equal(nrow(enrichment$PubChemBioactivity), 3)
  expect_equal(nrow(enrichment$PubChemBioAssayDetails), 3)
  expect_true(any(enrichment$PubChemBioAssayDetails$AID == 9001 &
                    enrichment$PubChemBioAssayDetails$TargetName ==
                      "prostaglandin G/H synthase 1 [Homo sapiens]" &
                    enrichment$PubChemBioAssayDetails$EndpointNames ==
                      "IC50; Efficacy"))
  expect_true(any(enrichment$PubChemBioAssayDetails$AID == 9002 &
                    enrichment$PubChemBioAssayDetails$TargetOrganism ==
                      "Homo sapiens"))
  expect_true(any(enrichment$PubChemBioactivity$Value == "Active"))
  expect_true(any(enrichment$PubChemBioactivity$ActivityValue == 0.35 &
                    enrichment$PubChemBioactivity$ActivityUnit == "uM"))
  expect_true(any(enrichment$ChemicalBioassays$AID == 9001 &
                    enrichment$ChemicalBioassays$ActivityClass == "active" &
                    enrichment$ChemicalBioassays$ActivityValue == 0.35 &
                    enrichment$ChemicalBioassays$TargetOrganism == "Homo sapiens" &
                    enrichment$ChemicalBioassays$EndpointNames == "IC50; Efficacy"))
  expect_true(any(enrichment$ChemicalBioactivities$ActivityDirection == "inhibitor" &
                    enrichment$ChemicalBioactivities$BioactivityDomain ==
                      "enzyme_protein_target" &
                    enrichment$ChemicalBioactivities$TargetName == "PTGS1"))
  expect_true(any(enrichment$ChemicalTargets$TargetName == "PTGS1" &
                    enrichment$ChemicalTargets$TargetGeneID == "5743" &
                    enrichment$ChemicalTargets$TargetTaxonomyID == "9606" &
                    enrichment$ChemicalTargets$TargetOrganism == "Homo sapiens" &
                    enrichment$ChemicalTargets$ActiveAssayCount == 1))
  expect_true(any(enrichment$ChemicalPotencies$PotencyMetric == "IC50" &
                    enrichment$ChemicalPotencies$PotencyBucket == "sub_1_uM" &
                    enrichment$ChemicalPotencies$TargetName == "PTGS1"))
  expect_true(any(enrichment$ChemicalTerms$Domain == "bioactivity" &
                    enrichment$ChemicalTerms$TermClean == "active"))
  expect_true(any(enrichment$ChemicalTerms$Domain == "bioactivity" &
                    enrichment$ChemicalTerms$TermClean == "ptgs1"))
  expect_true(any(enrichment$ChemicalTraits$TraitType == "target" &
                    enrichment$ChemicalTraits$TraitGroup == "target_name" &
                    enrichment$ChemicalTraits$TraitValue == "PTGS1"))
  expect_true(any(enrichment$ChemicalTraits$TraitType == "target" &
                    enrichment$ChemicalTraits$TraitGroup == "target_gene" &
                    enrichment$ChemicalTraits$TraitValue == "5743"))
  expect_true(any(enrichment$ChemicalTraits$TraitType == "target" &
                    enrichment$ChemicalTraits$TraitGroup == "target_organism" &
                    enrichment$ChemicalTraits$TraitValue == "Homo sapiens"))
  expect_true(any(enrichment$ChemicalTraits$TraitType == "bioactivity" &
                    enrichment$ChemicalTraits$TraitGroup == "endpoint_metric" &
                    enrichment$ChemicalTraits$TraitValue == "IC50"))
  expect_true(any(enrichment$ChemicalTraits$TraitType == "potency" &
                    enrichment$ChemicalTraits$TraitGroup == "potency_bucket" &
                    enrichment$ChemicalTraits$TraitValue == "sub_1_uM"))
  expect_true(any(enrichment$ChemicalTraitOntology$OntologyDomain == "bioactivity" &
                    enrichment$ChemicalTraitOntology$OntologyGroup == "target_gene" &
                    enrichment$ChemicalTraitOntology$OntologyTerm == "5743" &
                    enrichment$ChemicalTraitOntology$OntologyID == "5743" &
                    enrichment$ChemicalTraitOntology$OntologyIDSource ==
                      "NCBI Gene" &
                    enrichment$ChemicalTraitOntology$OntologyIDURL ==
                      "https://www.ncbi.nlm.nih.gov/gene/5743"))
  expect_true(any(enrichment$ChemicalTraitOntology$OntologyDomain == "bioactivity" &
                    enrichment$ChemicalTraitOntology$OntologyGroup == "target_taxonomy" &
                    enrichment$ChemicalTraitOntology$OntologyTerm == "9606" &
                    enrichment$ChemicalTraitOntology$OntologyIDSource ==
                      "NCBI Taxonomy"))
  expect_true(any(enrichment$ChemicalTraitOntology$OntologyDomain == "bioactivity" &
                    enrichment$ChemicalTraitOntology$OntologyGroup == "endpoint_metric" &
                    enrichment$ChemicalTraitOntology$OntologyTerm == "ic50" &
                    enrichment$ChemicalTraitOntology$OntologyID == "AID:9001" &
                    enrichment$ChemicalTraitOntology$OntologyIDURL ==
                      "https://pubchem.ncbi.nlm.nih.gov/bioassay/9001"))
  expect_false("target__target_name__ptgs1" %in%
                 colnames(enrichment$ChemicalTraitMatrix))
  expect_true(all(c("target__target_organism__homo_sapiens",
                    "bioactivity__endpoint_metric__ic50",
                    "potency__potency_bucket__sub_1_um") %in%
                    colnames(enrichment$ChemicalTraitMatrix)))
  bioactivity_trait_matrix = chemicalTraitMatrix(enrichment,
                                                 profile = "bioactivity")
  expect_true(all(c("target__target_name__ptgs1",
                    "target__target_gene__5743",
                    "target__target_organism__homo_sapiens",
                    "bioactivity__endpoint_metric__ic50",
                    "potency__potency_bucket__sub_1_um") %in%
                    colnames(bioactivity_trait_matrix)))
  expect_equal(bioactivity_trait_matrix$target__target_name__ptgs1[[1]], 1)

  derived = enrichment$DerivedGroups
  expect_true(derived$has_bioactivity)
  expect_true(derived$has_active_bioactivity)
  expect_equal(derived$bioassay_count, 3)
  expect_equal(derived$active_assay_count, 1)
  expect_equal(derived$inactive_assay_count, 1)
  expect_equal(derived$inconclusive_assay_count, 1)
  expect_equal(derived$target_count, 3)
  expect_equal(derived$active_target_count, 1)
  expect_equal(derived$assay_detail_count, 3)
  expect_equal(derived$protein_target_count, 3)
  expect_equal(derived$human_target_count, 3)
  expect_true(grepl("Homo sapiens", derived$target_organisms))
  expect_equal(derived$min_potency_uM, 0.35)
  expect_equal(derived$potent_activity_count, 1)
  expect_true(grepl("enzyme_protein_target", derived$bioactivity_domains))
  expect_true(grepl("active_bioactivity", derived$biomedical_context))
  expect_true(grepl("potent_activity", derived$bioactivity_context))
})

test_that("chemical trait helpers summarize and compare multiple compounds", {
  traits = data.frame(
    Query = c("compound_a", "compound_a", "compound_a",
              "compound_b", "compound_b", "compound_b",
              "compound_c"),
    CID = c(1L, 1L, 1L, 2L, 2L, 2L, 3L),
    TraitType = c("hazard", "pathway", "bioactivity",
                  "hazard", "pathway", "bioactivity",
                  "sensory"),
    TraitGroup = c("hazard_code", "pathway_group", "activity_class",
                   "hazard_code", "pathway_group", "activity_class",
                   "flavor"),
    TraitValue = c("H319", "Lipid metabolism", "active",
                   "H319", "Amino acid metabolism", "inactive",
                   "sweet"),
    SourceDatabase = c("PubChem", "KEGG", "PubChem BioAssay",
                       "PubChem", "KEGG", "PubChem BioAssay",
                       "FEMA"),
    Confidence = c("high", "high", "high", "high", "medium", "medium",
                   "high"),
    stringsAsFactors = FALSE
  )

  full_matrix = chemicalTraitMatrix(traits, profile = "full")
  expect_true(all(c("hazard__hazard_code__h319",
                    "pathway__pathway_group__lipid_metabolism",
                    "bioactivity__activity_class__active") %in%
                    colnames(full_matrix)))

  ontology = chemicalTraitOntology(traits)
  expect_true(any(ontology$Query == "compound_a" &
                    ontology$OntologyDomain == "safety" &
                    ontology$OntologyGroup == "ghs_hazard_code" &
                    ontology$OntologyTerm == "h319" &
                    ontology$OntologyID == "H319" &
                    ontology$OntologyIDSource == "GHS" &
                    ontology$SourceTraitKey == "hazard__hazard_code__h319"))
  expect_true(any(ontology$Query == "compound_a" &
                    ontology$OntologyDomain == "metabolism" &
                    ontology$OntologyGroup == "pathway_group" &
                    ontology$OntologyTerm == "lipid metabolism"))
  expect_true(any(ontology$Query == "compound_a" &
                    ontology$OntologyDomain == "bioactivity" &
                    ontology$OntologyGroup == "assay_outcome" &
                    ontology$OntologyTerm == "active"))
  expect_true(any(ontology$Query == "compound_c" &
                    ontology$OntologyDomain == "sensory" &
                    ontology$OntologyGroup == "flavor" &
                    ontology$OntologyTerm == "sweet"))
  expect_no_duplicate_keys(ontology,
                           c("Query", "CID", "OntologyDomain",
                             "OntologyGroup", "OntologyTerm"))
  ontology_evidence = chemicalTraitEvidence(traits,
                                            keys = "H319",
                                            type = "ontology")
  expect_equal(nrow(ontology_evidence), 2)
  expect_true(all(ontology_evidence$EvidenceType == "ontology"))
  expect_true(all(ontology_evidence$ExternalID == "H319"))
  trait_key_evidence = chemicalTraitEvidence(
    traits,
    keys = "hazard__hazard_code__h319",
    type = "trait"
  )
  expect_equal(nrow(trait_key_evidence), 2)
  expect_true(all(trait_key_evidence$EvidenceType == "trait"))
  expect_true(all(trait_key_evidence$AnalysisKey ==
                    "hazard__hazard_code__h319"))
  report = chemicalTraitReport(traits)
  report_a = report[report$Query == "compound_a", , drop = FALSE]
  expect_equal(nrow(report), 3)
  expect_equal(report_a$TraitCount, 3)
  expect_equal(report_a$OntologyTermCount, 3)
  expect_true(grepl("ghs_hazard_code=H319", report_a$SafetySignals,
                    fixed = TRUE))
  expect_true(grepl("pathway_group=Lipid metabolism",
                    report_a$MetabolismSignals,
                    fixed = TRUE))
  expect_true(grepl("assay_outcome=active", report_a$BioactivitySignals,
                    fixed = TRUE))
  expect_true(grepl("GHS:H319", report_a$ExternalIDs, fixed = TRUE))
  expect_true(grepl("compound_b (0.2)", report_a$NearestNeighbors,
                    fixed = TRUE))
  expect_equal(report_a$NearestNeighborSimilarity, 0.2)
  expect_true(grepl("safety__ghs_hazard_code__h319",
                    report_a$SharedOntologyKeys,
                    fixed = TRUE))

  ontology_matrix = chemicalTraitOntologyMatrix(traits)
  expect_true(all(c("safety__ghs_hazard_code__h319",
                    "metabolism__pathway_group__lipid_metabolism",
                    "bioactivity__assay_outcome__active",
                    "sensory__flavor__sweet") %in%
                    colnames(ontology_matrix)))
  expect_equal(ontology_matrix$safety__ghs_hazard_code__h319[
    ontology_matrix$Query == "compound_a"
  ], 1)
  limited_ontology_matrix = chemicalTraitOntologyMatrix(ontology,
                                                        max_terms = 2)
  expect_equal(ncol(limited_ontology_matrix), 4)

  taxonomy_ontology = chemicalTraitOntology(data.frame(
    Query = "camellia",
    CID = 4133L,
    TraitType = "taxonomy",
    TraitGroup = "family",
    TraitValue = "Theaceae",
    SourceDatabase = "LOTUS",
    EvidenceID = "4442",
    Confidence = "high",
    stringsAsFactors = FALSE
  ))
  expect_true(any(taxonomy_ontology$OntologyDomain == "ecology" &
                    taxonomy_ontology$OntologyGroup == "taxonomic_family" &
                    taxonomy_ontology$OntologyID == "4442" &
                    taxonomy_ontology$OntologyIDSource == "NCBI Taxonomy" &
                    taxonomy_ontology$OntologyIDURL ==
                      "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4442"))

  summary_by_type = chemicalTraitSummary(traits, by = "type")
  expect_true(any(summary_by_type$Query == "compound_a" &
                    summary_by_type$TraitType == "hazard" &
                    summary_by_type$TraitCount == 1 &
                    summary_by_type$HighConfidenceCount == 1))

  summary_by_source = chemicalTraitSummary(traits,
                                           by = "source",
                                           min_confidence = "high")
  expect_false(any(summary_by_source$Query == "compound_b" &
                     summary_by_source$SourceDatabase == "KEGG"))

  similarity = chemicalTraitSimilarity(traits, profile = "full")
  ab = similarity[similarity$QueryA == "compound_a" &
                    similarity$QueryB == "compound_b", , drop = FALSE]
  expect_equal(nrow(ab), 1)
  expect_equal(ab$SharedTraitCount, 1)
  expect_equal(ab$UnionTraitCount, 5)
  expect_equal(ab$JaccardSimilarity, 0.2)
  expect_true(grepl("hazard__hazard_code__h319", ab$SharedTraits,
                    fixed = TRUE))

  similarity_from_matrix = chemicalTraitSimilarity(full_matrix)
  expect_equal(nrow(similarity_from_matrix), 3)
  expect_true(all(c("compound_a", "compound_b") %in%
                    c(similarity_from_matrix$QueryA,
                      similarity_from_matrix$QueryB)))
})

test_that("LOTUS taxonomy follows PubChem taxon links and ignores table placeholders", {
  profile = pubchemProfile("methyl salicylate",
                           profile = "minimal",
                           sources = "LOTUS - the natural products occurrence database",
                           cache = FALSE,
                           throttle = 0,
                           request_fun = fixture_lotus_taxonomy_request)
  profiles = .pubchem_extract_profiles(profile)
  classes = .normalized_chemical_classes(
    compounds = "methyl salicylate",
    data_list = list(reactives = data.frame()),
    pubchem_profiles = profiles,
    kegg = .categorate_empty_kegg_profile(),
    cid_lookup = c("methyl salicylate" = 4133)
  )
  taxonomy = .normalized_chemical_taxonomy(profiles)
  occurrences = .normalized_chemical_occurrences(taxonomy)

  expect_equal(nrow(profiles$LOTUSProfile), 2)
  expect_equal(nrow(profiles$PubChemClassificationProfile), 4)
  expect_equal(nrow(occurrences), 3)
  expect_true(any(profiles$PubChemClassificationProfile$TreeType == "chemical" &
                    profiles$PubChemClassificationProfile$ClassName ==
                      "Simple phenolic acids" &
                    grepl("Shikimates and Phenylpropanoids",
                          profiles$PubChemClassificationProfile$ClassPath)))
  expect_true(any(classes$ClassType == "natural_product_superclass" &
                    classes$ClassName == "Shikimates and Phenylpropanoids"))
  expect_true(any(classes$ClassType == "natural_product_subclass" &
                    classes$ClassName == "Simple phenolic acids"))
  expect_true(any(profiles$LOTUSProfile$TaxonomyID == "4442" &
                    profiles$LOTUSProfile$Organism == "Camellia sinensis" &
                    grepl("Theaceae", profiles$LOTUSProfile$Taxonomy)))
  expect_true(any(occurrences$TaxonomyID == "4442" &
                    occurrences$OccurrenceType == "natural_product" &
                    occurrences$Organism == "Camellia sinensis" &
                    occurrences$Family == "Theaceae"))
  summary = .categorate_occurrence_summary(
    list(ChemicalOccurrences = occurrences),
    "methyl salicylate"
  )
  expect_equal(summary$occurrence_count, 3)
  expect_equal(summary$organism_count, 3)
  expect_equal(summary$kingdom_count, 2)
  expect_equal(summary$family_count, 3)
  expect_equal(summary$dominant_kingdom, "Plantae")
  expect_equal(summary$dominant_family, "mixed")
  expect_equal(summary$taxonomic_breadth, "multi_kingdom")
  expect_true(summary$is_plant_occurring)
  expect_true(summary$is_fungal_occurring)
  expect_true(any(taxonomy$TaxonomyID == "4442" &
                    taxonomy$Kingdom == "Plantae" &
                    taxonomy$Family == "Theaceae" &
                    taxonomy$Genus == "Camellia" &
                    taxonomy$Species == "Camellia sinensis"))
  expect_true(any(taxonomy$TaxonomyID == "108899" &
                    taxonomy$Kingdom == "Fungi" &
                    taxonomy$Family == "Hymenochaetaceae" &
                    taxonomy$Genus == "Phellinus" &
                    taxonomy$Species == "Phellinus tremulae"))
  expect_true(any(is.na(taxonomy$TaxonomyID) &
                    taxonomy$ExtractionRule == "lotus_biological_tree" &
                    taxonomy$Kingdom == "Plantae" &
                    taxonomy$Family == "Lamiaceae" &
                    taxonomy$Genus == "Salvia" &
                    taxonomy$Species == "Salvia officinalis"))
  expect_true(any(occurrences$Organism == "Salvia officinalis" &
                    occurrences$ExtractionRule == "lotus_biological_tree"))
  bad_text = paste(c(profiles$LOTUSProfile$Organism,
                     profiles$LOTUSProfile$Taxonomy,
                     taxonomy$TaxonomyTerms,
                     classes$ClassName),
                   collapse = "; ")
  expect_false(grepl("Salicylate has|consolidatedcompoundtaxonomy|(^|; )115(;|$)|(^|; )142(;|$)",
                     bad_text,
                     ignore.case = TRUE))
})

test_that("PubChem duplicate facts are collapsed without losing provenance", {
  annotations = data.frame(
    Query = c("aspirin", "aspirin", "aspirin"),
    CID = c(2244L, 2244L, 2244L),
    Heading = c("Safety and Hazards", "Source: FDA/SPL Indexing Data",
                "Safety and Hazards"),
    HeadingPath = c("Safety and Hazards", "FDA SPL > Safety",
                    "Safety and Hazards"),
    Name = c("Hazard Statement", "Hazard Statement", "Signal Word"),
    Value = c("H319 Warning", "H319 Warning", "Warning"),
    CleanValue = c("H319 Warning", "H319 Warning", "Warning"),
    ValueNumeric = c(NA_real_, NA_real_, NA_real_),
    Unit = c(NA_character_, NA_character_, NA_character_),
    UnitClean = c(NA_character_, NA_character_, NA_character_),
    Source = c("Fixture Safety", "Fixture Safety", "Fixture Safety"),
    SourceURL = c("https://example.test/a", "https://example.test/b",
                  "https://example.test/a"),
    PubChemURL = c("https://pubchem.test/heading",
                   "https://pubchem.test/source",
                   "https://pubchem.test/heading"),
    stringsAsFactors = FALSE
  )

  deduped_annotations = .pubchem_dedupe_annotations(annotations)

  expect_equal(nrow(deduped_annotations), 2)
  hazard_row = deduped_annotations[deduped_annotations$Name == "Hazard Statement",
                                   , drop = FALSE]
  expect_equal(nrow(hazard_row), 1)
  expect_true(grepl("Safety and Hazards", hazard_row$Heading, fixed = TRUE))
  expect_true(grepl("Source: FDA/SPL Indexing Data", hazard_row$Heading,
                    fixed = TRUE))
  expect_true(grepl("https://example.test/a", hazard_row$SourceURL,
                    fixed = TRUE))
  expect_true(grepl("https://example.test/b", hazard_row$SourceURL,
                    fixed = TRUE))

  safety_profile = data.frame(
    Query = c("aspirin", "aspirin", "aspirin"),
    CID = c(2244L, 2244L, 2244L),
    SignalWord = c("Warning", "Warning", "Danger"),
    HazardCode = c("H319", "H319", "H301"),
    PrecautionCode = c(NA_character_, NA_character_, NA_character_),
    HazardClass = c("Eye irritation", "Safety and Hazards",
                    "Acute toxicity"),
    HazardStatement = c("Causes serious eye irritation",
                        "Causes serious eye irritation",
                        "Toxic if swallowed"),
    RawValue = c("H319 Causes serious eye irritation",
                 "H319 Causes serious eye irritation",
                 "H301 Toxic if swallowed"),
    Source = c("GHS", "GHS", "GHS"),
    SourceURL = c("https://example.test/ghs-a",
                  "https://example.test/ghs-b",
                  "https://example.test/ghs-c"),
    PubChemURL = c("https://pubchem.test/heading",
                   "https://pubchem.test/source",
                   "https://pubchem.test/heading"),
    stringsAsFactors = FALSE
  )

  deduped_profile = .pubchem_dedupe_profile(
    table = safety_profile,
    key_cols = c("Query", "CID", "SignalWord", "HazardCode",
                 "PrecautionCode", "HazardStatement"),
    collapse_cols = c("HazardClass", "RawValue", "Source", "SourceURL",
                      "PubChemURL"),
    fallback_cols = c("RawValue", "HazardClass")
  )

  expect_equal(nrow(deduped_profile), 2)
  h319 = deduped_profile[deduped_profile$HazardCode == "H319", , drop = FALSE]
  expect_equal(nrow(h319), 1)
  expect_true(grepl("Eye irritation", h319$HazardClass, fixed = TRUE))
  expect_true(grepl("Safety and Hazards", h319$HazardClass, fixed = TRUE))
  expect_true(grepl("https://example.test/ghs-a", h319$SourceURL,
                    fixed = TRUE))
  expect_true(grepl("https://example.test/ghs-b", h319$SourceURL,
                    fixed = TRUE))

  sparse_profile = data.frame(
    Query = c("aspirin", "aspirin", "aspirin"),
    CID = c(2244L, 2244L, 2244L),
    ExtractedFact = c(NA_character_, NA_character_, NA_character_),
    RawValue = c("first source-only note", "second source-only note",
                 "first source-only note"),
    SourceURL = c("https://example.test/sparse-a",
                  "https://example.test/sparse-b",
                  "https://example.test/sparse-c"),
    stringsAsFactors = FALSE
  )

  deduped_sparse = .pubchem_dedupe_profile(
    table = sparse_profile,
    key_cols = c("Query", "CID", "ExtractedFact"),
    collapse_cols = c("RawValue", "SourceURL"),
    fallback_cols = "RawValue"
  )

  expect_equal(nrow(deduped_sparse), 2)
})
