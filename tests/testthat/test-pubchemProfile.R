fixture_pubchem_request = function(url) {
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
      TPSA = "63.6",
      Complexity = "212"
    )))))
  }

  if (grepl("/synonyms/JSON$", url)) {
    return(list(InformationList = list(Information = list(list(
      CID = 2244,
      Synonym = list("aspirin", "acetylsalicylic acid")
    )))))
  }

  if (grepl("heading=Mass%20Spectrometry|heading=GC-MS|heading=MS-MS", url)) {
    return(list(Record = list(
      Reference = list(list(
        ReferenceNumber = 1,
        SourceName = "Fixture Spectra",
        URL = "https://example.test/spectra"
      )),
      Section = list(list(
        TOCHeading = "Mass Spectrometry",
        Section = list(list(
          TOCHeading = "GC-MS",
          Information = list(list(
            Name = "Top m/z peaks",
            ReferenceNumber = 1,
            Value = list(StringWithMarkup = list(
              list(String = "120 100"),
              list(String = "92 80")
            ))
          ))
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
        Section = list(list(
          TOCHeading = "GHS Classification",
          Information = list(
            list(
              Name = "Signal",
              ReferenceNumber = 2,
              Value = list(StringWithMarkup = list(list(String = "Warning")))
            ),
            list(
              Name = "Hazard Statement",
              ReferenceNumber = 2,
              Value = list(StringWithMarkup = list(list(String = "Causes serious eye irritation")))
            )
          )
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
        Information = list(list(
          Name = "Boiling Point",
          ReferenceNumber = 3,
          Value = list(StringWithMarkup = list(list(String = "140 C")))
        ))
      ))
    )))
  }

  if (grepl("/assaysummary/JSON$", url)) {
    return(list(Table = list(Row = list(list(AID = 123, ActivityOutcome = "Active")))))
  }

  if (grepl("/assay/aid/123/description/JSON$", url)) {
    return(list(PC_AssayContainer = list(list(assay = list(descr = list(
      aid = list(id = 123, version = 1),
      aid_source = list(db = list(name = "Fixture BioAssay",
                                  source_id = list(str = "FIX123"))),
      name = "qHTS Assay for Inhibitors of PTGS1",
      description = list("Assay Overview: PTGS1 inhibition assay"),
      protocol = list("Measure prostaglandin synthesis inhibition."),
      comment = list("Keywords: inhibitor, cyclooxygenase"),
      xref = list(
        list(xref = list(pmid = 123456)),
        list(xref = list(gene = 5743)),
        list(xref = list(taxonomy = 9606)),
        list(xref = list(dburl = "https://example.test/bioassay"))
      ),
      results = list(
        list(tid = 1, name = "Potency"),
        list(tid = 2, name = "Efficacy")
      ),
      target = list(list(
        name = "prostaglandin G/H synthase 1 [Homo sapiens]",
        descr = "Human PTGS1 (COX-1)",
        mol_id = list(protein_accession = "NP_000953.2"),
        organism = list(org = list(
          taxname = "Homo sapiens",
          common = "human",
          db = list(list(db = "taxon", tag = list(id = 9606)))
        ))
      )),
      activity_outcome_method = 2,
      project_category = 1
    ))))))
  }

  list(Record = list(Section = list()))
}

test_that("pubchemProfile returns structured profile tables from PubChem responses", {
  profile = pubchemProfile("aspirin",
                           profile = "full",
                           cache = FALSE,
                           throttle = 0,
                           request_fun = fixture_pubchem_request)

  expect_s3_class(profile, "uaf_pubchem_profile")
  expect_equal(profile$identity$CID, 2244)
  expect_equal(profile$properties$MolecularFormula, "C9H8O4")
  expect_true("TPSA" %in% colnames(profile$properties))
  expect_true("acetylsalicylic acid" %in% profile$synonyms$Synonym)
  expect_true(any(profile$spectra$Value == "120 100"))
  expect_true(any(profile$safety$Value == "Warning"))
  expect_true(any(profile$experimental$Name == "Boiling Point"))
  expect_true("CleanValue" %in% colnames(profile$annotations))
  expect_equal(profile$experimental$ValueNumeric[profile$experimental$Name == "Boiling Point"][[1]],
               140)
  expect_true(any(profile$bioactivity$Value == "Active"))
  expect_true(all(c("AID", "ActivityOutcome", "ActivityClass",
                    "ActivityValue", "AssayName", "TargetAccession") %in%
                    colnames(profile$bioactivity)))
  expect_true(any(profile$bioactivity$ActivityClass == "active"))
  expect_equal(nrow(profile$bioassay_details), 1)
  expect_true(any(profile$bioassay_details$TargetName ==
                    "prostaglandin G/H synthase 1 [Homo sapiens]" &
                    profile$bioassay_details$TargetGeneID == "5743" &
                    profile$bioassay_details$TargetTaxonomyID == "9606"))
  expect_true(nrow(profile$provenance) > 0)
})

test_that("pubchemProfile resolves InChIKey queries through the InChIKey endpoint", {
  inchikey = "GWCVYMLSGTZAHQ-HNNXBMFYSA-N"
  request_fun = function(url) {
    if (grepl("/pug/compound/name/", url)) {
      stop("InChIKey query should not use the name endpoint")
    }
    if (grepl("/pug/compound/inchikey/GWCVYMLSGTZAHQ-HNNXBMFYSA-N/cids/JSON$",
              url)) {
      return(list(IdentifierList = list(CID = list(162905822))))
    }
    if (grepl("/property/", url)) {
      return(list(PropertyTable = list(Properties = list(list(
        CID = 162905822,
        Title = "Fixture InChIKey compound",
        MolecularFormula = "C22H23NO6",
        InChIKey = inchikey,
        CanonicalSMILES = "COc1ccccc1",
        IsomericSMILES = "COc1ccccc1"
      )))))
    }
    if (grepl("/synonyms/JSON$", url)) {
      return(list(InformationList = list(Information = list(list(
        CID = 162905822,
        Synonym = list("fixture compound")
      )))))
    }
    list()
  }

  profile = pubchemProfile(
    inchikey,
    profile = "minimal",
    include_annotations = FALSE,
    request_fun = request_fun
  )

  expect_equal(profile$identity$CID, 162905822)
  expect_equal(profile$identity$MatchStatus, "resolved_inchikey")
  expect_match(profile$identity$SourceURL, "/compound/inchikey/")
  expect_equal(profile$properties$InChIKey, inchikey)
})

test_that("pubchemProfile validates empty compound input", {
  expect_error(pubchemProfile(c("", NA), request_fun = fixture_pubchem_request),
               "at least one non-empty")
})

test_that("pubchemProfile can skip PUG-View annotations for property-only workflows", {
  property_only_request = function(url) {
    if (grepl("/pug_view/", url)) {
      stop("PUG-View should not be requested when include_annotations = FALSE")
    }
    fixture_pubchem_request(url)
  }

  profile = pubchemProfile("aspirin",
                           profile = "ms",
                           cache = FALSE,
                           throttle = 0,
                           include_annotations = FALSE,
                           request_fun = property_only_request)

  expect_equal(profile$identity$CID, 2244)
  expect_equal(profile$properties$MolecularFormula, "C9H8O4")
  expect_equal(nrow(profile$annotations), 0)
  expect_equal(nrow(profile$source_annotations), 0)
})

test_that("pubchemProfile retries conservative name aliases for PubChem identity", {
  requested = character()
  alias_request = function(url) {
    requested <<- c(requested, utils::URLdecode(url))
    if (grepl("/pug/compound/name/beta-pinene/cids/JSON$", utils::URLdecode(url))) {
      return(list(IdentifierList = list(CID = list(14896))))
    }
    if (grepl("/property/", url)) {
      return(list(PropertyTable = list(Properties = list(list(
        CID = 14896,
        Title = "beta-Pinene",
        MolecularFormula = "C10H16",
        MolecularWeight = "136.238",
        IUPACName = "6,6-dimethyl-2-methylidenebicyclo[3.1.1]heptane",
        InChIKey = "WTARULDDTDQWMU-UHFFFAOYSA-N",
        SMILES = "CC1(C2CCC(=C)C1C2)C"
      )))))
    }
    if (grepl("/synonyms/JSON$", url)) {
      return(list(InformationList = list(Information = list(list(
        CID = 14896,
        Synonym = list("beta-pinene")
      )))))
    }
    list()
  }
  compound = paste0(intToUtf8(0x03b2), "-pinene")

  profile = pubchemProfile(compound,
                           profile = "minimal",
                           cache = FALSE,
                           throttle = 0,
                           include_annotations = FALSE,
                           request_fun = alias_request)

  expect_equal(profile$identity$CID, 14896)
  expect_equal(profile$identity$MatchStatus, "resolved_alias")
  expect_equal(profile$identity$QueriedName, "beta-pinene")
  expect_true(any(grepl("/name/beta-pinene/cids/JSON", requested,
                        fixed = TRUE)))
  expect_equal(profile$properties$SMILES, "CC1(C2CCC(=C)C1C2)C")
})
