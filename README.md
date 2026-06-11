
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

For LOTUS specifically, use a local index for serious plant panels. The LOTUS
simple web API is useful for small smoke tests, but common species can return
very large unpaged payloads. uafR therefore supports `lotus_index`, a flat CSV,
TSV, JSON, JSONL, NDJSON, RDS, or data-frame index produced from LOTUS MongoDB,
Wikidata, SDF metadata, or another LOTUS export. `standardizeLotusIndex()`
normalizes flexible column names such as `allTaxa`, `traditional_name`,
`lotus_id`, `inchikey`, `smiles`, `molecular_formula`, `doi`, and `pmid`.
`queryLotusIndex()` returns the same `PlantCompoundOccurrences` schema used by
the full resolver. When `lotus_index` is supplied, `resolvePlantPhytochemistry()`
uses the local index instead of the live LOTUS API.

Build the compact index once, then reuse it across projects:

``` r
lotus_build = buildLotusIndex(
  input = "lotus_flat_export.jsonl",
  out_file = "lotus_compact_index.csv",
  overwrite = TRUE
)

lotus_build$BuildSummary
```

The same step can be run from Terminal:

``` sh
Rscript tools/build_lotus_index.R \
  --input lotus_flat_export.jsonl \
  --out-file lotus_compact_index.csv \
  --overwrite
```

The official LOTUS SMILES and SDF downloads are useful for structure work, but
they do not carry the source-backed species occurrence records needed for
species-first phytochemistry. For plant runs, use the official MongoDB ZIP
download and flatten the `lotusUniqueNaturalProduct.bson` collection into a
taxon-compound CSV first. The repository includes a Python standard-library
helper for this step, so MongoDB command-line tools are not required:

``` sh
mkdir -p lotus_cache/downloads lotus_cache/exports

curl -L \
  -o lotus_cache/downloads/LOTUSlatest.zip \
  https://lotus.naturalproducts.net/download/mongo

PYTHONDONTWRITEBYTECODE=1 python3 tools/flatten_lotus_mongo_dump.py \
  --input lotus_cache/downloads/LOTUSlatest.zip \
  --out-file lotus_cache/exports/LOTUS_mongo_flat.csv \
  --compact-index-file lotus_cache/exports/LOTUS_compact_index.csv \
  --lookup-dir lotus_cache/exports/LOTUS_lookup_index \
  --overwrite
```

`LOTUS_mongo_flat.csv` is the auditable one-row-per-compound/taxon/reference
export. `LOTUS_compact_index.csv` is the smaller, analysis-ready index to pass
to uafR for small panels or ad hoc inspection. `LOTUS_lookup_index/` is the
production path for larger plant panels because uafR can read only the shards
needed for the submitted species, genus, and family keys instead of scanning
the complete CSV. Build products under `lotus_cache/` are ignored by package
builds and should be treated as local data cache files, not source files.

If the compact CSV already exists, build only the lookup directory without
rewriting the flat export:

``` sh
PYTHONDONTWRITEBYTECODE=1 python3 tools/flatten_lotus_mongo_dump.py \
  --from-compact-index lotus_cache/exports/LOTUS_compact_index.csv \
  --lookup-dir lotus_cache/exports/LOTUS_lookup_index \
  --overwrite
```

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
  lotus_index = "lotus_cache/exports/LOTUS_lookup_index",
  enrich_compounds = TRUE,
  chemical_library = library_data,
  detail = "research",
  cache = TRUE,
  cache_dir = "uafR_plant_cache",
  max_pubmed_records = 25,
  max_provider_records = 100,
  enrichment_batch_size = 25,
  resume_enrichment = TRUE
)

lotus_occurrences = queryLotusIndex(
  plants,
  lotus_index = "lotus_cache/exports/LOTUS_lookup_index",
  taxon_fallback = c("species", "genus")
)

phyto$SpeciesChemistrySummary
phyto$PlantCompoundOccurrences
phyto$PlantContextEvidence
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

# Candidate, fallback, literature, and unresolved rows should be reviewed
# before they are promoted into analysis-ready occurrence evidence.
review = plantPhytochemistryReviewTable(phyto)
if (nrow(review) > 0) {
  review$review_decision[1] = "promote_curated"
  review$reviewed_by[1] = "researcher name"
  review$review_note[1] = "Source reports the compound from leaf tissue."
  review$proposed_citation_or_url[1] = "https://doi.org/example"
  review$proposed_plant_part[1] = "leaf"
  review$proposed_method[1] = "LC-MS"
  phyto_reviewed = applyPlantPhytochemistryReview(phyto, review)
}

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

Use the comparability layer before clustering, ordination, scoring, or group
comparisons. `metabolite` is a broad biological term, while `volatile`
describes an analytical or physicochemical fraction; those should not be
compared as if they were the same axis. uafR therefore records biological
domain, biosynthetic family, analytical behavior, comparison scope, confidence,
source-backed classification fields, and caveats in `ChemistryComparability`.
The classifier prioritizes normalized source fields from `ChemicalClasses`,
`LOTUSProfile`, `PubChemClassifications`, `ChemicalTraitOntology`,
`ChemicalTerms`, `KEGGClassifications`, `KEGGPathways`, and `DerivedGroups`
before it falls back to compound-name patterns.
Biological-context evidence is tracked separately in `PlantContextEvidence` so
plant-part, tissue, and method filters can be audited instead of being hidden
inside a matrix. Provider-level coverage and review burden are summarized in
`ProviderContextAudit`, which helps identify whether context gaps come from a
specific source, candidate-only literature evidence, or sparse provider records.

``` r
# Inspect how each plant-compound row was assigned to a comparison scope.
phyto$ChemistryComparability[, c(
  "species", "compound_name", "metabolism_domain",
  "biosynthetic_family", "chemical_behavior",
  "comparison_scope", "comparison_group",
  "comparability_confidence", "comparability_basis",
  "classification_source_table", "classification_source_value"
)]

# Inspect provider-level biological-context coverage before trusting a
# context-specific matrix.
phyto$ProviderContextAudit[, c(
  "source_database", "occurrence_count", "context_known_fraction",
  "review_required_context_fraction", "audit_status",
  "recommended_action"
)]

# Add source-backed context from local literature/source text. Rows must match
# occurrence evidence by PMID or DOI; context is not filled unless the source
# text contains plant-part, tissue, or method terms linked to the species/genus
# or compound in a chemical-source context.
source_text = data.frame(
  pmid = "12345678",
  doi = "10.1000/example",
  title = "Phytochemical constituents from Salix nigra leaves",
  abstract = paste(
    "Salicin was identified from leaves of Salix nigra.",
    "The extract was analyzed by LC-MS."
  )
)

phyto = enrichPlantContextEvidence(
  phyto,
  context_sources = source_text
)

# Optional live PubMed context fetches are explicit, cached, capped, and
# prioritized toward source records with useful context text or broad
# species/compound coverage.
phyto = enrichPlantContextEvidence(
  phyto,
  fetch_pubmed = TRUE,
  cache = TRUE,
  cache_dir = "uafR_plant_cache/context",
  max_sources = 50
)

# Compare specialized plant chemistry separately from primary metabolism.
specialized_matrix = plantComparableChemistryMatrix(
  phyto,
  comparison_scope = "specialized_metabolites",
  feature = "comparison_group",
  mode = "binary"
)

# Compare volatile specialized chemistry only within the volatile fraction.
volatile_matrix = plantComparableChemistryMatrix(
  phyto,
  comparison_scope = "volatile_specialized_metabolites",
  feature = "comparison_group",
  mode = "binary"
)

# Compare leaf volatile chemistry only when plant-part or method context is
# known. This prevents root exudate, whole-plant, and unspecified-extract rows
# from being mixed into a leaf-focused matrix.
leaf_volatile_matrix = plantComparableChemistryMatrix(
  phyto,
  comparison_scope = "volatile_specialized_metabolites",
  feature = "comparison_group",
  plant_part_group = "leaf",
  require_context = TRUE,
  mode = "binary"
)

# Primary metabolites are a different comparison scope.
primary_matrix = plantComparableChemistryMatrix(
  phyto,
  comparison_scope = "primary_metabolites",
  feature = "comparison_group",
  mode = "binary"
)
```

Unknown or broad rows stay in the audit table but are excluded from comparable
matrices by default. This is intentional: a missing class is not evidence that a
compound belongs in a broad mixed analysis.

Before scaling to hundreds of plants, run a small provider pilot and inspect
the output. The pilot workflow writes the same batch audit tables plus a compact
species summary, a QA report, a review-needed table, and context-aware
comparable matrices for specialized, volatile-specialized, leaf-associated,
root/exudate-associated, and primary-metabolism chemistry.

``` r
plants = plantPhytochemistryPilotPanel()$species

pilot = runPlantPhytochemistryPilot(
  plants = plants,
  sources = c("lotus", "knapsack", "npass", "pubchem", "pubmed", "pubtator"),
  out_dir = "plant_phytochemistry_pilot",
  cache_dir = "uafR_plant_cache",
  species_chunk_size = 10,
  compound_resolution_profile = "identity",
  max_pubmed_records = 25,
  max_provider_records = 100,
  request_timeout = 30,
  overwrite = TRUE
)

pilot$PilotSummary
pilot$PilotQAReport
pilot$ProviderContextAudit
pilot$PilotMatrices$matrix_volatile_specialized_metabolites
pilot$PilotMatrices$matrix_root_exudate_associated_chemistry
```

The command-line wrapper accepts a CSV plant list:

``` sh
Rscript tools/run_plant_phytochemistry_pilot.R \
  --plant-csv plants.csv \
  --species-col species \
  --out-dir plant_phytochemistry_pilot \
  --cache-dir uafR_plant_cache \
  --lotus-index lotus_compact_index.csv \
  --sources lotus,knapsack,npass,pubchem,pubmed,pubtator \
  --compound-resolution-profile identity \
  --max-pubmed-records 25 \
  --max-provider-records 100 \
  --request-timeout 30 \
  --overwrite
```

The wrapper can also use the built-in 15-species panel:

``` sh
Rscript tools/run_plant_phytochemistry_pilot.R \
  --default-panel \
  --out-dir plant_phytochemistry_pilot \
  --compound-resolution-profile identity \
  --overwrite
```

Start with `compound_resolution_profile = "identity"` for a provider/data-depth
pilot. Move to `"research"` only after the occurrence evidence, context
coverage, provider context audit, and review burden look reasonable.

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

For larger plant lists, use the staged batch workflow. It discovers
plant-compound evidence first, writes resumable checkpoints, filters to
analysis-ready direct or curated records by default, and then performs a fast
PubChem identity-only resolution pass. Richer `detail = "research"` enrichment
should be run later on a reviewed subset of compounds.
The identity pass preserves chemically meaningful alpha/beta/gamma and
plus/minus prefixes in compound keys, and reports deterministic PubChem alias
matches with `MatchStatus = "resolved_alias"` in the identity audit table.
When `lotus_index` is supplied, the identity pass first recovers source-backed
LOTUS SMILES, InChIKeys, formulas, and CIDs by LOTUS record ID. PubChem name
lookup is then used only for remaining gaps. A separate
`CompoundIdentityReview` table flags source-ambiguous structures and rows where
the source provides a structure but the displayed label looks like a class,
mixture, plant product, or other non-discrete compound name.

``` r
plants = c("Salix nigra", "Camellia sinensis", "Zea mays")

phyto_batch = runPlantPhytochemistryBatch(
  plants = plants,
  sources = c("lotus", "knapsack", "pubmed", "pubtator"),
  out_dir = "plant_phytochemistry_batch",
  cache_dir = "uafR_plant_cache",
  species_chunk_size = 25,
  compound_resolution_profile = "identity",
  max_pubmed_records = 25,
  max_provider_records = 100,
  request_timeout = 30,
  compound_batch_size = 100,
  resume = TRUE,
  overwrite = TRUE
)

phyto_batch$BatchRunManifest
phyto_batch$SpeciesChemistrySummary
phyto_batch$CompoundResolution
phyto_batch$CompoundIdentityReview
```

Completed identity review worksheets can be reapplied to the result object so
manual structure decisions are reproducible:

``` r
identity_review = phyto_batch$CompoundIdentityReview

# Example: only fill these fields after checking the LOTUS source record,
# PubChem record, DOI, or another source-backed structure record.
identity_review$review_decision[1] = "update_identity"
identity_review$reviewed_by[1] = "researcher name"
identity_review$review_note[1] = "Source record supports this structure."
identity_review$proposed_compound_name[1] = "reviewed compound name"
identity_review$proposed_smiles[1] = "source-backed SMILES"
identity_review$proposed_inchikey[1] = "source-backed InChIKey"
identity_review$proposed_molecular_formula[1] = "source-backed formula"
identity_review$proposed_resolution_source[1] = "manual_identity_review"

phyto_batch_reviewed = applyPlantCompoundIdentityReview(
  phyto_batch,
  identity_review,
  reviewer = "researcher name"
)
```

The batch output directory includes `all_occurrences.csv`,
`analysis_ready_occurrences.csv`, `review_required_occurrences.csv`,
`compound_identity_resolution.csv`, `compound_identity_review.csv`,
`species_chemistry_summary.csv`, `species_chemistry_matrix.csv`,
`chemistry_comparability.csv`,
`comparable_chemistry_matrix.csv`, `provider_diagnostics.csv`,
`context_coverage_report.csv`, `batch_chunk_manifest.csv`, and
`run_manifest.json`.

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

### Package release checks

Use the release-check wrapper before sharing a source archive, building a
student bundle, or opening a pull request:

``` sh
Rscript tools/check_package_release.R
```

The script builds a source tarball outside the repository and runs
`R CMD check --no-manual` against the tarball. That is the production check
path. Running `R CMD check .` directly on the live checkout can report local
artifacts such as `.git`, `.DS_Store`, `.Rhistory`, or generated check
directories that are not included in the source package.

Live PubChem/NCI integration tests are opt-in. Set `UAFR_RUN_LIVE_TESTS=true`
before running tests, or pass `--run-live-tests` to
`tools/check_package_release.R`, when you want to exercise live web calls.

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
