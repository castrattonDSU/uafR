# uafR <a href="https://castrattonDSU.github.io/uafR/"><img src="man/figures/logo.svg" align="right" height="139" alt="uafR evidence-prism hex logo" /></a>

[![R-CMD-check](https://github.com/castrattonDSU/uafR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/castrattonDSU/uafR/actions/workflows/R-CMD-check.yaml)
[![PLOS ONE](https://img.shields.io/badge/PLOS%20ONE-10.1371%2Fjournal.pone.0306202-0A7BBB)](https://doi.org/10.1371/journal.pone.0306202)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**uafR** is an R package for processing tentative compound identifications from
GC–MS workflows. It prepares vendor-exported peak tables, aggregates detections
across samples, extracts user-selected compounds, standardizes abundance, and
adds structured chemical context from PubChem and KEGG.

## What uafR does

The package supports three connected workflows:

1. **GC–MS data processing**
   Prepare raw peak tables with `spreadOut()`, extract target compounds with
   `mzExacto()`, and standardize abundance with `standardifyIt()`.

2. **Chemical interpretation**
   Group and annotate compound names with `categorate()`, then select compounds
   for downstream extraction with `exactoThese()`.

3. **Research-grade enrichment**
   Retrieve reusable PubChem and KEGG tables with `pubchemProfile()`,
   `keggProfile()`, and the enriched `categorate()` workflow; validate and export
   results with `validateCategorateResult()` and
   `exportCategorateWorkbook()`.

## Production workflow guide

uafR now includes package-level helpers for choosing workflows, checking API
stability, planning large runs, and documenting scientific limits:

``` r
uafRWorkflowGuide()
uafRApiStability()
uafRProviderContracts()
uafRClaimGuidance()
```

Stable APIs are intended for downstream scripts and training materials.
Experimental APIs are usable but may still gain columns, diagnostics, or
provider-specific hardening as large plant and database workflows mature.
Project-specific scripts under `tools/` remain examples or wrappers and should
not be treated as general package APIs unless they are promoted into exported
functions.

Large species-first projects should be planned before live queries are run:

``` r
plan = planPlantChemistryRun(
  plants = "project_species.csv",
  compounds = "resolved_compounds.csv",
  sources = c("lotus", "pubmed", "pubtator"),
  cache_dir = "uafR_plant_cache",
  lotus_index = "lotus_cache/exports/LOTUS_lookup_index"
)

plan$Summary
plan$ProviderPlan
plan$OutputEstimates
plan$Recommendations
```

Use `inspectUafRCache()` and `summarizeUafRCache()` to audit local cache
coverage before repeating expensive PubChem, KEGG, PubMed, PubTator, or plant
provider workflows. These planning helpers do not query live web services.

## Installation

uafR requires R 4.2 or newer.

```r
install.packages(c("remotes", "BiocManager"))

BiocManager::install(
  c("ChemmineR", "fmcsR"),
  ask = FALSE,
  update = FALSE
)

remotes::install_github(
  "castrattonDSU/uafR",
  dependencies = c("Depends", "Imports"),
  upgrade = "never",
  build_vignettes = FALSE
)
```

Load the package:

```r
library(uafR)
```

Windows users may need
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) when a dependency must
be compiled from source. macOS users may need the Xcode Command Line Tools.

## GC–MS input

`spreadOut()` expects a data frame containing these columns:

| Column | Meaning |
|---|---|
| `Component.RT` | Component retention time |
| `Component.Area` | Integrated component area |
| `Base.Peak.MZ` | Base-peak mass-to-charge ratio |
| `File.Name` | Sample or source file |
| `Compound.Name` | Tentative compound identity |
| `Match.Factor` | Library match score |

Read vendor output without changing its column names:

```r
raw_gcms <- read.csv(
  "path/to/gcms-export.csv",
  check.names = FALSE
)

required_columns <- c(
  "Component.RT",
  "Component.Area",
  "Base.Peak.MZ",
  "File.Name",
  "Compound.Name",
  "Match.Factor"
)

stopifnot(all(required_columns %in% names(raw_gcms)))
```

## Core workflow: prepare and extract GC–MS data

Prepare the peak table:

```r
spread <- spreadOut(raw_gcms)
```

Choose compounds directly or use a match-factor threshold:

```r
targets <- unique(
  raw_gcms$Compound.Name[
    !is.na(raw_gcms$Match.Factor) &
    raw_gcms$Match.Factor >= 80
  ]
)
```

Extract target compounds across samples:

```r
extracted <- mzExacto(
  data_in = spread,
  chemicals = targets,
  decontaminate = TRUE
)

head(extracted)
```

`mzExacto()` returns compound identities, optimal retention time, exact mass,
best match factor, and aggregated component area across samples.

## Core workflow: annotate and select compounds

uafR includes an example chemical library:

```r
data("library_data", package = "uafR")

query_compounds <- c(
  "Linalool",
  "Methyl salicylate",
  "Limonene",
  "alpha-Pinene"
)
```

Run the standard annotation workflow:

```r
annotations <- categorate(
  compounds = query_compounds,
  chemical_library = library_data,
  input_format = "wide"
)

names(annotations)
```

Standard mode returns eight tables:

```text
reactives
LOTUS
KEGG
FEMA
FDA_SPL
FMCS
FunctionalGroups
BestChemMatch
```

Select compounds using database, structure, or library-group evidence:

```r
selected <- exactoThese(
  annotations,
  subsetBy = "Database",
  subsetArgs = c("LOTUS", "KEGG")
)
```

Use the selected names in the GC–MS extraction workflow:

```r
selected_data <- mzExacto(
  data_in = spread,
  chemicals = selected
)
```

## Research-grade enrichment

Use `detail = "research"` to append normalized PubChem and KEGG tables while
preserving the original eight standard tables:

```r
research_result <- categorate(
  compounds = query_compounds,
  chemical_library = library_data,
  input_format = "wide",
  detail = "research",
  cache = TRUE,
  throttle = 0.2,
  assay_detail_limit = 0
)
```

<<<<<<< HEAD
Frequently used outputs include:
||||||| 6833f2c
## Example Mass Spectrometry Workflows
=======
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

Use structural Tanimoto similarity when the question is molecular resemblance,
not shared database traits. `chemicalTraitSimilarity()` compares normalized
traits from categorate enrichment. `chemicalTanimotoSimilarity()` and
`plantChemicalTanimotoSimilarity()` compare PubChem Fingerprint2D bit vectors
and report true compound-pair structural similarity. The plant helper keeps
species names, compound IDs, evidence tiers, and plant context attached so the
result can be joined directly to phylogenetic, ecological, remediation, or
sample-similarity analyses.

``` r
# Pairwise molecular similarity for a plant phytochemistry result. For small
# projects, the pairwise tables can be returned in memory.
plant_tanimoto = plantChemicalTanimotoSimilarity(
  phyto,
  cache = TRUE,
  cache_dir = "uafR_plant_cache/tanimoto",
  return_group_compound_pairs = FALSE
)

plant_tanimoto$PlantPairTanimotoSummary[, c(
  "species_a", "species_b", "compound_pair_count",
  "shared_compound_count", "mean_tanimoto", "median_tanimoto",
  "p95_tanimoto", "max_tanimoto",
  "compound_pair_count_ge_0_85"
)]

# For larger plant panels, stream the large pairwise tables to compressed CSV
# files and keep the compact species-pair summary in the returned object.
plant_tanimoto = plantChemicalTanimotoSimilarity(
  phyto,
  cache = TRUE,
  cache_dir = "uafR_plant_cache/tanimoto",
  out_dir = "plant_tanimoto_export",
  return_compound_pairs = TRUE,
  return_group_compound_pairs = TRUE
)

plant_tanimoto$ExportManifest
plant_tanimoto$PlantPairTanimotoSummary
```

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
>>>>>>> origin/plant-phytochemistry-hardening

```r
research_result$PubChemProperties
research_result$ChemicalTraits
research_result$ChemicalTraitOntology
research_result$ChemicalTraitMatrix
research_result$ChemicalTraitReport
research_result$ChemicalMeasurements
research_result$ChemicalMeasurementSummary
research_result$ChemicalHazards
research_result$ChemicalOccurrences
research_result$KEGGPathways
research_result$SourceCoverage
```

Validate the result before downstream analysis:

```r
audit <- validateCategorateResult(research_result)

<<<<<<< HEAD
||||||| 6833f2c
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
=======
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
>>>>>>> origin/plant-phytochemistry-hardening
audit$Summary
audit$TableQuality
audit$SourceDiagnostics
audit$Issues
```

Export a curated CSV bundle:

```r
manifest <- exportCategorateWorkbook(
  research_result,
  path = "categorate-export",
  format = "csv",
  overwrite = TRUE
)

manifest
```

To write an Excel workbook, install `openxlsx` or `writexl` and provide an
`.xlsx` path.

## Source-specific profiles

### PubChem

<<<<<<< HEAD
```r
pubchem <- pubchemProfile(
  query_compounds,
  profile = "ms",
  cache = TRUE
)
||||||| 6833f2c
`keggProfile()` can also be used directly when the goal is KEGG-specific
annotation. It resolves names or supplied KEGG IDs, parses KEGG flat-file
records, follows KEGG links, and returns pathways, reactions, enzymes, modules,
identifiers, and reproducible pathway/enzyme classifications.
=======
Large plant-chemistry projects often run `categorate()` enrichment in resumable
batches. Use the batch-analysis helpers to audit those outputs before handing
them to a downstream statistics or figure-making workspace:

``` r
batches = readCategorateBatchDirectory(
  "dsi_categorate_research_pubchem_cid_20260611/categorate_batches"
)

batch_summary = summarizeCategorateBatches(batches)
batch_summary[, c("BatchIndex", "Status", "ValidationStatus",
                  "ResolvedCIDCount", "PubChemPropertyRatio")]

combined = combineCategorateTables(
  batches,
  tables = c("ChemicalTraits", "PubChemProperties", "ChemicalTraitEvidence")
)

manifest = exportPlantChemistryAnalysisBundle(
  categorate_batches = batches,
  path = "plant_chemistry_analysis_bundle",
  plant_membership = "dsi_species_compound_membership.csv",
  species_pair_tanimoto = "dsi_species_pair_tanimoto_summary.csv",
  resolved_compounds = "dsi_resolved_compounds_for_categorate.csv",
  pubchem_fingerprints = "dsi_pubchem_fingerprints.csv",
  plant_compound_pair_tanimoto = "dsi_plant_compound_pair_tanimoto.csv.gz",
  file_references = c(
    plant_compound_pairs = "dsi_plant_compound_pair_tanimoto.csv.gz",
    compound_pairs = "dsi_compound_pair_tanimoto.csv.gz",
    chemical_traits_full = "combined_chemical_traits_if_exported_separately.csv.gz",
    trait_evidence_full = "combined_trait_evidence_if_exported_separately.csv.gz"
  ),
  plant_list = "project_species.csv",
  project_id = "plant_project",
  tables = c("ChemicalTraitSummary", "DerivedGroups", "PubChemProperties",
             "SourceCoverage", "ValidationSummary", "ValidationIssues"),
  format = "csv",
  include_comparable_tanimoto = TRUE,
  overwrite = TRUE
)
manifest

validation = validatePlantChemistryAnalysisBundle(
  "plant_chemistry_analysis_bundle",
  use_python = TRUE,
  use_pandas = TRUE
)
validation$Summary
```

The exported `BatchSummary` should be checked first. Batches with `Status =
"error"` should not be interpreted; batches with `Status = "incomplete"` need
review because they have fewer PubChem property rows than resolved CIDs at the
chosen threshold. Large compressed pairwise files can be recorded in
`FileReferences` instead of copied into every bundle. For very large runs, keep
row-level tables such as `ChemicalTraits` and `ChemicalTraitEvidence` as
separate targeted exports unless the downstream analysis needs a full combined
copy. Very wide sparse tables such as `ChemicalTraitMatrix` are useful for
focused modeling, but they can be slow and awkward to move for thousands of
compounds; keep the batch directory in `FileReferences` and extract a focused
matrix profile when the analysis plan is settled.

For CSV bundles, `exportPlantChemistryAnalysisBundle()` finalizes the handoff by
default. The finalized bundle includes `00_DataDictionary.csv`,
`03b_PlantCompoundMembershipEnriched.csv`, `11b_SourceCoverageSummary.csv`,
`12b_ValidationOverview.csv`, `14_PlantChemistrySummary.csv`,
`15_PlantChemistryMissingSpecies.csv`, `16_EvidenceGradeSummary.csv`,
`17_ReviewRequiredOccurrences.csv`, model-ready feature matrices
(`18_` through `21_` plus evidence/context matrices as `24_` through `28_`),
`README.md`, and `METHODS_TEXT.md`. If
`include_comparable_tanimoto = TRUE` and a plant-compound pair Tanimoto table is
supplied, the bundle also includes scope- and group-filtered comparable
Tanimoto summaries as `22_ComparableScopeTanimotoSummary.csv` and
`23_ComparableGroupTanimotoSummary.csv`.
The enriched membership table joins occurrence evidence to PubChem fingerprints,
PubChem/categorate properties, natural-product context, source coverage, and the
uafR comparable-chemistry fields (`comparison_scope`, `comparison_group`,
`comparison_subgroup`, `comparability_confidence`, `comparability_basis`,
`comparable_for_matrix`, and `comparison_caveat`). These fields are intended to
help downstream projects compare like with like; they do not turn public
database records into confirmed sample measurements.

For reusable project handoffs, use `runPlantChemistryProject()` or the CLI
wrapper. The runner accepts a plant list plus either curated plant-compound
records, a saved plant phytochemistry result, cached categorate batches, or an
explicit live-discovery request. It writes a project manifest and a finalized
bundle without hardcoding project-specific paths:

``` r
project = runPlantChemistryProject(
  plant_list = "project_species.csv",
  metadata = "project_plant_metadata.csv",
  plant_compounds = "curated_or_discovered_plant_compounds.csv",
  output_dir = "plant_chemistry_project_run",
  plant_compound_pair_tanimoto = "plant_compound_pair_tanimoto.csv.gz",
  project_id = "plant_project",
  overwrite = TRUE
)

project$Project
validatePlantChemistryAnalysisBundle(project$Project$bundle_dir)
```

``` sh
Rscript tools/run_plant_chemistry_project.R \
  --plant-list project_species.csv \
  --metadata project_plant_metadata.csv \
  --plant-compounds curated_or_discovered_plant_compounds.csv \
  --plant-compound-pair-tanimoto plant_compound_pair_tanimoto.csv.gz \
  --out-dir plant_chemistry_project_run \
  --project-id plant_project \
  --overwrite true
```

Use `finalizePlantChemistryAnalysisBundle()` when a bundle has already been
written and only needs the analysis-ready handoff files refreshed. Use
`validatePlantChemistryAnalysisBundle()` before handing a bundle to another
project. It checks manifest row/column counts, required columns, accidental row
index columns, and CSV parser consistency. Optional Python and pandas checks can
be enabled on machines where those tools are available.

For model-ready exports outside a full bundle, `exportPlantChemistryFeatureSet()`
can write count, binary, fraction, and confidence-weighted matrices:

``` r
features = exportPlantChemistryFeatureSet(
  membership = enriched_membership,
  modes = c("count", "binary", "fraction", "confidence"),
  path = "plant_feature_set",
  overwrite = TRUE
)

features$Manifest
features$SpeciesMetadata
```

Evidence and classification helpers support reproducible subsetting and
review:

``` r
plantOccurrenceEvidenceDictionary()
chemistryComparisonDictionary()

direct = filterPlantEvidenceDirect(enriched_membership)
comparable = filterPlantEvidenceComparable(enriched_membership)
review = filterPlantEvidenceReviewRequired(enriched_membership)
```

Use `standardizeChemistryClassificationOverrides()` and
`applyChemistryClassificationOverrides()` when a project needs source-backed or
human-reviewed chemistry-scope corrections. Unknown chemistry is retained for
audit but excluded from comparable matrices by default.

`keggProfile()` can also be used directly when the goal is KEGG-specific
annotation. It resolves names or supplied KEGG IDs, parses KEGG flat-file
records, follows KEGG links, and returns pathways, reactions, enzymes, modules,
identifiers, and reproducible pathway/enzyme classifications.
>>>>>>> origin/plant-phytochemistry-hardening

pubchem$identity
pubchem$properties
pubchem$spectra
pubchem$provenance
```

Available profiles are:

- `"minimal"` — identifiers, formula, exact mass, and molecular weight
- `"ms"` — mass-spectrometry-oriented descriptors and annotations
- `"safety"` — safety, hazard, and experimental-property annotations
- `"bioactivity"` — PubChem BioAssay summaries
- `"full"` — all supported sections

### KEGG

```r
kegg <- keggProfile(
  compounds = query_compounds,
  pubchem_profile = pubchem,
  link_targets = c("pathway", "reaction", "enzyme")
)

kegg$matches
kegg$pathways
kegg$reactions
kegg$enzymes
kegg$classifications
kegg$provenance
```

## Standardization

<<<<<<< HEAD
Standardize `mzExacto()` output with an internal standard:
||||||| 6833f2c
For package checks, live PubChem/NCI integration tests are opt-in. Set `UAFR_RUN_LIVE_TESTS=true` before running tests when you want to exercise live web calls.
=======
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
>>>>>>> origin/plant-phytochemistry-hardening

```r
standardized <- standardifyIt(
  data_in = extracted,
  standard_type = "Internal",
  standard_used = "Tetradecane",
  IS_ng = 190.5,
  IS_uL = 1,
  collect_time = 1,
  sample_amt = 1
)
```

For external standardization, provide a calibration matrix:

```r
data("ExternalStandard_data", package = "uafR")

standardized <- standardifyIt(
  data_in = extracted,
  standard_type = "External",
  ES_calibration = ExternalStandard_data
)
```

## Main user-facing functions

| Function | Purpose |
|---|---|
| `spreadOut()` | Prepare raw GC–MS peak tables for downstream processing |
| `mzExacto()` | Extract and aggregate user-selected compounds across samples |
| `standardifyIt()` | Standardize abundance using internal or external standards |
| `categorate()` | Annotate compounds and compare them with chemical-library groups |
| `exactoThese()` | Select compounds from `categorate()` results |
| `pubchemProfile()` | Retrieve structured PubChem profiles |
| `keggProfile()` | Retrieve structured KEGG profiles |
| `chemicalTraitMatrix()` | Build analysis-ready trait matrices |
| `chemicalTraitOntologyMatrix()` | Build controlled ontology matrices |
| `chemicalTraitEvidence()` | Trace matrix or ontology terms back to evidence |
| `chemicalTraitReport()` | Build compact per-compound research summaries |
| `chemicalTraitSummary()` | Summarize trait breadth by compound or source |
| `chemicalTraitSimilarity()` | Compare compounds by shared and distinct traits |
| `chemicalMeasurementSummary()` | Summarize normalized numeric measurements |
| `validateCategorateResult()` | Audit schemas, completeness, duplicates, and source coverage |
| `exportCategorateWorkbook()` | Export curated CSV or Excel result bundles |

## Reproducibility and network access

Several workflows query public web services. Results can change as external
databases evolve.

- Keep `cache = TRUE` for reusable PubChem and KEGG responses.
- Store result objects and validation tables with each analysis.
- Record the uafR version and session information.
- Use bounded request limits and the default throttling settings.
- Inspect `SourceCoverage`, `SourceDiagnostics`, and validation output before
  interpreting missing values.

Record the environment used for an analysis:

```r
packageVersion("uafR")
sessionInfo()
```

## Documentation and support

- Function reference: <https://castrattonDSU.github.io/uafR/>
- Issues and bug reports: <https://github.com/castrattonDSU/uafR/issues>
- Source code: <https://github.com/castrattonDSU/uafR>

## Citation

When using uafR, cite:

> Stratton CA, Thompson Y, Zio K, Morrison WR III, Murrell EG (2024).
> uafR: An R package that automates mass spectrometry data processing.
> *PLOS ONE* 19(7): e0306202.
> <https://doi.org/10.1371/journal.pone.0306202>

## License

uafR is released under the MIT License.
