
# uafR - A new standard for mass spectrometry data processing

<!-- badges: start -->
<!-- badges: end -->

## Objective

An R package that automates GC-MS processing.

## Installation

### Student install from a private bundle

Student machines should install uafR from a private student bundle distributed
by the DSU dsDNA Core program. This keeps the repository private and gives
students a direct RStudio workflow.

The bundle contains a local uafR source archive, checksum verification,
preflight checks, installer and update scripts, verification scripts, an
offline acceptance test, a quick reference, and the training manual. The
installer still uses CRAN and Bioconductor for dependencies.

From the unzipped bundle folder, students should open `START_HERE.md`, then run
these scripts in RStudio:

``` r
source("verify_bundle_integrity.R")
source("preflight_check.R")
source("install_uafR_from_bundle.R")
source("run_student_acceptance_test.R")
```

After installation, confirm the package loads:

``` r
library(uafR)
data("library_data", package = "uafR")
```

### Build a student bundle

From the repository root:

``` sh
Rscript tools/build_student_bundle.R
```

The script creates:

``` text
student_bundle/uafR_student_bundle_<version>/
student_bundle/uafR_student_bundle_<version>.zip
```

The bundle folder includes:

``` text
START_HERE.md
README_STUDENT_INSTALL.md
uafR_QUICK_REFERENCE.md
CHECKSUMS.csv
verify_bundle_integrity.R
preflight_check.R
install_uafR_from_bundle.R
update_uafR_from_bundle.R
verify_uafR_install.R
run_student_acceptance_test.R
packages/uafR_<version>.tar.gz
training/uafR_training_manual.pdf
training/scripts/
training/data/
examples/test_install.R
```

### Developer install from a local checkout

Developers working from the private repository can install from the local source
checkout:

``` r
install.packages(c("remotes", "BiocManager"))
BiocManager::install(c("ChemmineR", "fmcsR"), ask = FALSE, update = FALSE)
remotes::install_local(
  ".",
  dependencies = c("Depends", "Imports"),
  upgrade = "never",
  build_vignettes = FALSE
)
```

### Optional live smoke test

After installation, a small live database smoke test can be run when internet
access is available:

``` r
library(uafR)
data("library_data", package = "uafR")

quick_result = categorate(
  compounds = c("aspirin", "caffeine"),
  chemical_library = library_data,
  input_format = "wide",
  detail = "research",
  cache = FALSE,
  throttle = 0.1,
  assay_detail_limit = 0
)

validateCategorateResult(quick_result)$Summary
quick_result$ChemicalTraitReport
quick_result$ChemicalMeasurementSummary
```

## Species-first plant phytochemistry workflows

uafR can now start from plant species names instead of a known compound list.
The plant phytochemistry resolver attempts to normalize plant names, collect
reported species-compound evidence from normalized provider adapters or curated
intake tables, include conservative literature/PubTator candidate evidence when
available, and then reuse uafR compound enrichment for resolved compounds.
When a `chemical_library` is supplied, enrichment can use the full
`categorate()` workflow, including library-based FMCS matching. When no
chemical library is supplied, uafR now falls back to `pubchemProfile()` and
still returns PubChem-derived properties, traits, matrices, validation, and
provenance.
Current live-capable public adapters include PubMed literature search,
PubTator candidate co-mentions, KNApSAcK organism-metabolite lookup, conservative
LOTUS API parsing when taxon evidence is present, and PubChem taxonomy
annotations after NCBI taxonomy resolution. NPASS species-source rows should be
supplied through curated intake or `provider_results` until a small stable
species-query endpoint is added.

This workflow is designed for plant, ecology, remediation, metabolomics,
chemical ecology, natural-products, and environmental chemistry projects. It
does not produce complete metabolomes. Absence of public records is not absence
of compounds. Literature co-mentions are candidate evidence unless curated.
PubChem taxonomy annotations are source-backed chemical associations, not proof
that a student's sample contains the compound. Direct species evidence is
stronger than genus or family fallback evidence. Always inspect validation and
provenance before interpretation.

``` r
library(uafR)
data("library_data", package = "uafR")

plants = c("Salix nigra", "Camellia sinensis", "Zea mays")

phyto = resolvePlantPhytochemistry(
  plants = plants,
  sources = c("lotus", "pubmed", "pubtator"),
  taxon_fallback = c("species", "genus"),
  enrich_compounds = TRUE,
  chemical_library = library_data,
  detail = "research",
  cache = TRUE,
  cache_dir = "uafR_plant_cache",
  max_pubmed_records = 25,
  max_provider_records = 100
)

phyto$SpeciesChemistrySummary
phyto$PlantCompoundOccurrences
validatePlantPhytochemistryResult(phyto)$Summary

# Direct species database records and curated rows are analysis-ready by
# default. Literature co-mentions and genus/family fallbacks remain available
# for review but are not treated as analysis-ready occurrence evidence.
analysis_ready = filterPlantPhytochemistryEvidence(phyto)

# Restrict analyses to compounds reported from a known biological context when
# that metadata is available. Supported groups include root_belowground, leaf,
# stem_shoot, bark_wood, flower, fruit_seed, aerial, whole_plant, and
# exudate_rhizosphere.
leaf_records = filterPlantPhytochemistryEvidence(
  phyto,
  plant_part_group = "leaf"
)

exportPlantPhytochemistryWorkbook(
  phyto,
  path = "plant_phytochemistry_export",
  format = "csv",
  overwrite = TRUE
)

exportPlantPhytochemistryWorkbook(
  phyto,
  path = "plant_phytochemistry_analysis_ready_export",
  format = "csv",
  preset = "analysis_ready",
  overwrite = TRUE
)
```

Use the full export for audit/provenance and the analysis-ready preset for
downstream matrices or figures where candidate co-mentions and taxon fallbacks
should be excluded.

For quick species-first discovery without a library, omit `chemical_library`.
This produces PubChem-only enrichment and skips FMCS library matching:

``` r
phyto = resolvePlantPhytochemistry(
  plants = c("Camellia sinensis", "Salix nigra"),
  sources = c("lotus", "knapsack", "pubmed"),
  enrich_compounds = TRUE,
  detail = "research",
  cache = TRUE,
  cache_dir = "uafR_plant_cache",
  max_provider_records = 50
)

names(phyto$SpeciesChemistryMatrix)
phyto$CompoundResolution
```

For projects that already have local curation, use the curated intake fallback:

``` r
curated = data.frame(
  species = "Salix nigra",
  compound_name = "salicin",
  source_database = "manual",
  citation_or_url = "https://example.org/source",
  evidence_tier = "manual_curated",
  plant_part = "bark",
  method = "LC-MS"
)

occurrences = standardizePlantCompoundIntake(curated)
phyto = resolvePlantPhytochemistry(
  plants = unique(curated$species),
  sources = character(),
  curated_data = occurrences,
  enrich_compounds = TRUE,
  chemical_library = library_data,
  detail = "research"
)
```

An offline training example is available at
`training/scripts/10_species_phytochemistry.R`. It uses simulated curated
species-compound rows and does not call live web services unless the
`UAFR_LIVE_PLANT_DISCOVERY` environment variable is set to `"true"`.

## Example Mass Spectrometry Workflows

These are basic examples of how to use core functions. The input .CSV file has strict column name/input data requirements. The column names MUST include: 'Component.RT', 'Component.Area', 'Base.Peak.MZ', 'File.Name', 'Compound.Name', and 'Match.Factor' in no particular order.

``` r
library(uafR)

input_dat = read.csv("your/gcms/dataset.csv")
```
 Component.RT  |  Base.Peak.MZ    |  Component.Area  |       Compound.Name        |  Match.Factor  |  File.Name  
:-------------:|:----------------:|:----------------:|:---------------------------|:--------------:|:------------:
8.229034       |84.00             |906.4701          |Pipradrol                   |62.62271        |Std_soln_07    
8.286703       |120.00            |209705.1878       |Methyl salicylate           |98.16152        |Std_soln_00a    
8.296408       |119.99            |30332.9022        |Methyl salicylate           |95.79911        |Std_soln_00    
8.303958       |120.00            |6476.4785         |Methyl salicylate           |86.29569        |Std_soln_07    
8.348031       |105.00            |420.8119          |3-Hexen-1-ol, benzoate, (Z)-|68.78156        |Std_soln_00    
**...**        |**...**           |**...**           |**...**                     |**...**         |**...**         

### In this example, the user knows what chemicals they are interested in:
``` r
input_spread = spreadOut(input_dat)
query_chemicals = c("Linalool", "Methyl Salicylate", "Limonene", "alpha-Thujene")

### extract the query_chemicals from the "spread out" input:
input_exacto = mzExacto(input_spread, query_chemicals)
```
### In this example, the user just wants to keep the top hits:
``` r
query_chemicals = input_dat$Compound.Name[input_dat$Match.Factor > 80]

input_exacto = mzExacto(input_spread, query_chemicals)
```

## Example Cheminformatics Workflow
``` r
## example usage for chemical informatics:
query_chemicals = c("Linalool", "Methyl Salicylate", "Limonene", "alpha-Thujene")
GroupA = c("Guaiacol",	"Tridecane",	"Ethyl heptanoate", "Caffeine")
GroupB = c("2-Aminothiazole", "Aspirin", "Octanoic acid", "alpha-Pinene", "Toluene")
chem_library = data.frame(cbind(GroupA, GroupB))

query_categorated = categorate(query_chemicals, chem_library, input_format = "wide")
```

## Detailed PubChem Enrichment

`pubchemProfile()` pulls richer PubChem data into reusable tables before downstream filtering or mass spectrometry extraction. This keeps web enrichment separate from `spreadOut()` and `mzExacto()`, so results can be cached, inspected, and tested independently.

``` r
query_chemicals = c("Methyl salicylate", "Octanal", "Undecane")

chem_profile = pubchemProfile(query_chemicals, profile = "ms")

chem_profile$identity
chem_profile$properties
chem_profile$synonyms
chem_profile$spectra
```

Available profiles are:

- `"minimal"`: PubChem CID, names, formula, identifiers, exact mass, and molecular weight.
- `"ms"`: minimal data plus chemical descriptors and mass spectrometry annotations.
- `"safety"`: minimal/descriptive data plus safety, hazard, and experimental property annotations.
- `"bioactivity"`: minimal/descriptive data plus PubChem assay summary data.
- `"full"`: all supported profile sections.

## aiNsect Molecular-Olfaction Export

`exportAiNsectMolOlfInputs()` converts paired wide essential-oil GC-MS profile
tables into aiNsect molecular-olfaction input files. The identity table should
contain treatment columns with compound names by ranked row. The abundance table
should contain the same treatments and rows with abundance values. Compound
structures are resolved through `pubchemProfile()`; uafR does not fabricate
SMILES, InChIKeys, CIDs, formulas, or abundance values.

``` sh
Rscript tools/export_ainsect_mololf_inputs.R \
  --chem-id-csv /Users/chasestratton/src/github/castrattonDSU/aiNsect_tracker/EO_PCA_2026/20240612-EO-gcms-data_all.csv \
  --chem-quant-csv /Users/chasestratton/src/github/castrattonDSU/aiNsect_tracker/EO_PCA_2026/20240612-EO-gcms-quant_all.csv \
  --out-dir /Users/chasestratton/src/github/castrattonDSU/aiNsect_tracker/EO_PCA_2026/mololf_export \
  --cache-dir /Users/chasestratton/src/github/castrattonDSU/aiNsect_tracker/EO_PCA_2026/pubchem_cache \
  --profile ms \
  --throttle 0.2
```

The export directory contains:

- `uafR_compounds.csv`: compound IDs, names, SMILES, InChIKeys, PubChem CIDs,
  molecular formulas, source labels, and notes.
- `treatment_compound_abundance.csv`: treatment/compound abundance rows with
  raw parsed abundance and within-treatment relative abundance.
- `uafR_compounds_unresolved.csv`: compounds that could not be exported because
  PubChem did not provide a usable SMILES.
- `uafR_mololf_export_summary.json`: run metadata, counts, output paths, and
  quality checks.
- `uafR_pubchem_identity_audit.csv` and `uafR_pubchem_properties_audit.csv`:
  optional PubChem audit tables.

## Research-Grade Database Enrichment

`categorate()` keeps the original eight-table output by default. Use
`detail = "research"` or `detail = "full"` when you want tidy PubChem and KEGG
tables appended to the categorated result. These tables are designed for
filtering, grouping, and downstream analyses rather than one-off text lookup.

``` r
query_categorated = categorate(
  query_chemicals,
  chem_library,
  input_format = "wide",
  detail = "research"
)

query_categorated$PubChemProperties
query_categorated$SafetyProfile
query_categorated$FEMAProfile
query_categorated$FDA_SPL_Profile
query_categorated$LOTUSProfile
query_categorated$PubChemClassifications
query_categorated$MeSHProfile
query_categorated$LiteratureProfile
query_categorated$ChemicalTerms
query_categorated$ChemicalTraits
query_categorated$ChemicalTraitOntology
query_categorated$ChemicalTraitMatrix
query_categorated$ChemicalTraitOntologyMatrix
query_categorated$ChemicalTraitEvidence
query_categorated$ChemicalTraitReport
query_categorated$ChemicalTraitSummary
query_categorated$ChemicalTraitSimilarity
query_categorated$ChemicalClasses
query_categorated$ChemicalMeasurements
query_categorated$ChemicalMeasurementSummary
query_categorated$ChemicalHazards
query_categorated$ChemicalUses
query_categorated$ChemicalBioassays
query_categorated$ChemicalBioactivities
query_categorated$ChemicalTargets
query_categorated$ChemicalPotencies
query_categorated$PubChemBioAssayDetails
query_categorated$ChemicalTaxonomy
query_categorated$ChemicalOccurrences
query_categorated$ChemicalPathwayRoles
query_categorated$KEGGReactionParticipants
query_categorated$KEGGPathways
query_categorated$KEGGClassifications
query_categorated$DerivedGroups
query_categorated$SourceCoverage
query_categorated$DataDictionary
query_categorated$TableQuality
query_categorated$SourceDiagnostics
query_categorated$ValidationIssues
query_categorated$ValidationSummary
```

The source-specific PubChem profile tables preserve the raw annotation text
while adding cleaned fields such as GHS hazard codes, FEMA/JECFA identifiers,
flavor or odor terms, FDA/SPL route and dosage-form terms, LOTUS occurrence
signals, PubChem classification-tree paths, MeSH pharmacologic actions, PubMed
IDs, DOI values, PubChem BioAssay activity summaries, and bounded BioAssay
description metadata when `detail = "full"` is used. The `assay_detail_limit`
argument controls how many BioAssay descriptions are fetched per query so broad
screens stay inspectable and polite to PubChem. `DerivedGroups` summarizes these
into analysis columns for metabolic, biomedical, ecological, sensory, safety,
bioactivity, and analytical context.

The normalized `Chemical*` tables go one step further and split evidence text
into discrete values for analysis. `ChemicalTerms` is a long-form vocabulary
table for grouping chemicals by sensory, biomedical, ecological, safety,
classification, and metabolic terms. `ChemicalTraits` unifies those cleaned
signals into one cross-source long table with `TraitType`, `TraitGroup`,
`TraitValue`, `SourceDatabase`, evidence, confidence, and a stable matrix key.
`ChemicalTraitOntology` maps whitelisted discrete traits into controlled domains
such as safety, physicochemical behavior, metabolism, ecology, sensory,
bioactivity, biomedical use, regulatory status, and reactivity. Every ontology
row keeps the source trait key, evidence text/URL, source database, extraction
rule, confidence, and source-backed identifiers where available: GHS hazard
codes, KEGG pathway/reaction/compound/EC IDs, MeSH tree numbers, NCBI Taxonomy
IDs, NCBI Gene IDs, PubChem BioAssay AIDs, and target accessions. `ChemicalTraitMatrix`
and `ChemicalTraitOntologyMatrix` convert curated high-value traits into binary
columns for filtering, clustering, ordination, heatmaps, and model inputs. Use
the helper functions to rebuild specialized matrices without rerunning web
requests:

``` r
core_matrix = query_categorated$ChemicalTraitMatrix
ontology = query_categorated$ChemicalTraitOntology
ontology_matrix = query_categorated$ChemicalTraitOntologyMatrix
ontology_evidence = query_categorated$ChemicalTraitEvidence
trait_report = query_categorated$ChemicalTraitReport
bioactivity_matrix = chemicalTraitMatrix(query_categorated, profile = "bioactivity")
kegg_matrix = chemicalTraitMatrix(query_categorated$ChemicalTraits, profile = "kegg")
confidence_matrix = chemicalTraitMatrix(
  query_categorated,
  profile = "full",
  mode = "confidence",
  min_confidence = "high",
  max_traits = 250
)
ontology_confidence_matrix = chemicalTraitOntologyMatrix(
  query_categorated,
  mode = "confidence",
  min_confidence = "medium",
  max_terms = 250
)

trait_summary = chemicalTraitSummary(query_categorated)
source_summary = chemicalTraitSummary(query_categorated, by = "source")
trait_similarity = chemicalTraitSimilarity(query_categorated, profile = "core")
report = chemicalTraitReport(query_categorated, min_confidence = "medium")
dictionary = chemicalDataDictionary("ChemicalTraitReport")
audit = validateCategorateResult(query_categorated)
audit$Summary
audit$TableQuality
audit$SourceDiagnostics
audit$Issues
export_manifest = exportCategorateWorkbook(
  query_categorated,
  "categorate_export",
  format = "csv",
  overwrite = TRUE
)
export_manifest
h319_evidence = chemicalTraitEvidence(
  query_categorated,
  keys = "safety__ghs_hazard_code__h319"
)
trait_evidence = chemicalTraitEvidence(
  query_categorated,
  keys = "hazard__hazard_code__h319",
  type = "trait"
)
```

`ChemicalMeasurements` stores numeric
properties and extracted experimental/toxicity values with raw units plus
standardized values, canonical units, measurement classes, relations, and
behavior bins. `ChemicalMeasurementSummary` condenses those values into one row
per compound-property combination for plotting and filtering. `ChemicalHazards`
separates GHS codes, hazard groups, routes, target organs, precaution codes, and
toxicity metrics. `ChemicalBioassays`, `ChemicalBioactivities`, `ChemicalTargets`,
and `ChemicalPotencies` split PubChem assay summaries into assay outcomes,
activity domains, target identifiers/names, target organism/taxonomy fields,
assay endpoint names, activity directions, and numeric potency values.
`PubChemBioAssayDetails` preserves the raw selected assay-description metadata
behind those normalized bioactivity rows. `ChemicalTaxonomy` extracts organism,
genus, species, family,
order, class, phylum, kingdom, cleaned taxonomy terms, and natural-product
class. `ChemicalClasses` includes discrete natural-product superclass, class,
and subclass rows from LOTUS chemical classification trees. `ChemicalOccurrences`
gives one row per chemical-organism occurrence with taxonomic ranks, source
evidence, occurrence type, and confidence, while `DerivedGroups` summarizes
natural-product classes, bioactivity breadth, occurrence breadth, dominant
kingdom/family, and plant/fungal/bacterial occurrence flags. `ChemicalPathwayRoles` and
`KEGGReactionParticipants` split KEGG pathways, enzymes, reactions, substrates,
and products into analysis-ready rows.

Every detailed `categorate()` result also carries a data-quality layer.
`DataDictionary` describes expected tables, columns, types, required fields,
roles, and allowed controlled values. `TableQuality` reports row counts,
required-column completeness, duplicate analysis-key counts, confidence counts,
and table status. `SourceDiagnostics` shows which source/query combinations
returned usable rows. `ValidationSummary` and `ValidationIssues` provide the
same audit generated by `validateCategorateResult()`, so downstream scripts can
stop early when a schema or extraction problem appears.

Use `exportCategorateWorkbook()` to share results outside R. It writes a
manifest plus clean CSV files by default, and can write an `.xlsx` workbook when
`openxlsx` or `writexl` is installed.

`keggProfile()` can also be used directly when the goal is KEGG-specific
annotation. It resolves names or supplied KEGG IDs, parses KEGG flat-file
records, follows KEGG links, and returns pathways, reactions, enzymes, modules,
identifiers, and reproducible pathway/enzyme classifications.

``` r
kegg_profile = keggProfile(
  c("aspirin", "glucose"),
  max_matches_per_query = 3,
  link_targets = c("pathway", "reaction", "enzyme")
)

kegg_profile$pathways
kegg_profile$reactions
kegg_profile$enzymes
kegg_profile$classifications
```

Linked KEGG IDs are resolved into names and definitions by default, so tables
such as `kegg_profile$reactions` include reaction names, definitions, and
equations when KEGG exposes them.

For package checks, live PubChem/NCI integration tests are opt-in. Set `UAFR_RUN_LIVE_TESTS=true` before running tests when you want to exercise live web calls.

## Combined Mass Spectrometry + Cheminformatics Workflow

``` r
query_chemicals = input_dat$Compound.Name[input_dat$Match.Factor > 70]
query_categorated = categorate(query_chemicals, chem_library, input_format = "wide")

## example of using the info from categorate() to get a user-defined set of chemicals with exactoThese():
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "All")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "reactives")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "LOTUS")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "KEGG")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "FEMA")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "FDA_SPL")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = c("reactives", "FEMA"))
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Greater Than", subset_input = 125)
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Less Than", subset_input = 205)
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Between", subset_input = c(125, 200))
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "Rings", subsetArgs2 = "Greater Than", subset_input = 1)
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "Groups", subsetArgs2 = "Greater Than", subset_input = 2)
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "Atoms", subsetArgs2 = "Greater Than", subset_input = 6)
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "NCharges", subsetArgs2 = "Greater Than", subset_input = 2)
these_chems = exactoThese(query_categorated, subsetBy = "Library", subsetArgs = "GroupB")

input_exacto = mzExacto(input_spread, these_chems)
```
