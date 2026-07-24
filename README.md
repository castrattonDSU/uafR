# uafR <a href="https://castrattonDSU.github.io/uafR/"><img src="man/figures/logo.svg" align="right" height="139" alt="uafR evidence-prism hex logo" /></a>

[![R-CMD-check](https://github.com/castrattonDSU/uafR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/castrattonDSU/uafR/actions/workflows/R-CMD-check.yaml)
[![PLOS ONE](https://img.shields.io/badge/PLOS%20ONE-10.1371%2Fjournal.pone.0306202-0A7BBB)](https://doi.org/10.1371/journal.pone.0306202)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**uafR** provides reproducible workflows for tentative mass-spectrometry
annotations, public-database chemical enrichment, species-first plant
phytochemistry evidence, comparable chemistry, structural similarity, and
validated research handoff bundles. Source diagnostics, evidence grades,
identity review, and scientific guardrails preserve uncertainty and provenance
for downstream analysis.

## What uafR does

The package supports five connected workflows:

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

4. **Species-first plant phytochemistry**
   Search supported public sources for reported plant-compound relationships,
   retain direct and fallback evidence separately, resolve compound identities,
   and build species summaries with `resolvePlantPhytochemistry()`.

5. **Comparable analysis and project handoff**
   Grade evidence, separate biologically comparable chemistry, calculate
   PubChem Fingerprint2D Tanimoto summaries, and export validated feature and
   analysis bundles.

## Choose a workflow

| Starting point or goal | Recommended entry point |
|---|---|
| GC-MS hit table | `spreadOut()`, `mzExacto()`, `standardifyIt()` |
| Compound names or PubChem CIDs | `categorate(detail = "research")` |
| Small plant species list | `resolvePlantPhytochemistry()` |
| Curated plant-compound table | `standardizePlantCompoundIntake()` and `runPlantChemistryProject()` |
| Large plant panel | `planPlantChemistryRun()` and `runPlantChemistryPanel()` |
| Compound or plant structural similarity | `chemicalTanimotoSimilarity()` or `plantChemicalTanimotoSimilarity()` |
| Model-ready plant chemistry features | `exportPlantChemistryFeatureSet()` |
| Portable analysis handoff | `finalizePlantChemistryAnalysisBundle()` |

The same guidance is available from R:

```r
uafRWorkflowGuide()
uafRApiStability()
uafRProviderContracts()
uafRClaimGuidance()
```

Functions labeled `stable` are intended for reusable downstream scripts.
Functions labeled `experimental` are available for research use but may gain
provider-specific fields or diagnostics before the next stable release.

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

Frequently used outputs include:

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

## Species-first plant phytochemistry

`resolvePlantPhytochemistry()` starts with plant names and attempts to discover
reported plant-compound relationships from enabled public providers. It returns
normalized evidence, provider diagnostics, compound-resolution tables, species
summaries, matrices, validation, a data dictionary, and provenance.

A bounded exploratory query looks like this:

```r
plants <- c("Salix nigra", "Camellia sinensis", "Zea mays")

phyto <- resolvePlantPhytochemistry(
  plants = plants,
  sources = c("lotus", "pubmed", "pubtator"),
  taxon_fallback = c("species", "genus"),
  enrich_compounds = FALSE,
  cache = TRUE,
  cache_dir = "uafR_plant_cache",
  max_pubmed_records = 25,
  max_provider_records = 100,
  throttle = 0.5
)

phyto$PlantCompoundOccurrences
phyto$LiteratureCandidates
phyto$ProviderDiagnostics
phyto$SpeciesChemistrySummary

plant_audit <- validatePlantPhytochemistryResult(phyto)
plant_audit$Summary
plant_audit$Issues
```

This example uses live services and may return different coverage as public
databases change. Start with a few plants, keep caching enabled, and inspect
`ProviderDiagnostics` before increasing the query size.

### Run the offline example

The package includes a small **simulated teaching fixture**. It tests the
workflow without network access; it is not evidence that the listed compounds
occur in the listed plants.

```r
example_dir <- system.file(
  "extdata",
  "offline_plant_chemistry",
  package = "uafR"
)

plants <- read.csv(
  file.path(example_dir, "plant_list.csv"),
  stringsAsFactors = FALSE
)
metadata <- read.csv(
  file.path(example_dir, "plant_metadata.csv"),
  stringsAsFactors = FALSE
)
plant_compounds <- read.csv(
  file.path(example_dir, "plant_compounds.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
compound_pairs <- read.csv(
  file.path(example_dir, "plant_compound_pair_tanimoto.csv"),
  stringsAsFactors = FALSE
)

offline_project <- runPlantChemistryProject(
  plant_list = plants,
  metadata = metadata,
  plant_compounds = plant_compounds,
  plant_compound_pair_tanimoto = compound_pairs,
  output_dir = "offline-plant-example",
  project_id = "simulated_offline_example",
  overwrite = TRUE
)

bundle_dir <- offline_project$Project$bundle_dir[[1]]
bundle_audit <- validatePlantChemistryAnalysisBundle(bundle_dir)
bundle_audit$Summary
```

The same example can be regenerated from a repository checkout:

```sh
Rscript tools/build_offline_plant_chemistry_example.R \
  --out-dir offline-plant-example
```

### Use a curated plant-compound table

Public-source discovery is necessarily incomplete. A locally reviewed table can
enter the same data model:

```r
curated <- read.csv("plant_compounds_reviewed.csv", stringsAsFactors = FALSE)
occurrences <- standardizePlantCompoundIntake(curated)

project <- runPlantChemistryProject(
  plant_list = unique(curated$species),
  plant_compounds = curated,
  output_dir = "plant-chemistry-project",
  project_id = "reviewed_curated_intake",
  overwrite = TRUE
)
```

At minimum, curated intake requires `species` and `compound_name`. Include
`source_database`, `source_record_id`, `citation_or_url`, `evidence_tier`,
`plant_part`, `tissue`, `method`, and review notes whenever those fields are
available.

### Use a local LOTUS index

For medium or large panels, use an official LOTUS export locally rather than
repeating broad live searches:

```r
lotus_build <- buildLotusIndex(
  input = "path/to/flat-lotus-export",
  out_file = "uafR_indexes/lotus_index.rds",
  overwrite = FALSE
)

lotus_hits <- queryLotusIndex(
  plants = c("Salix nigra", "Zea mays"),
  lotus_index = lotus_build$LotusIndex
)

phyto <- resolvePlantPhytochemistry(
  plants = plants,
  sources = c("lotus", "pubmed", "pubtator"),
  lotus_index = "uafR_indexes/lotus_index.rds",
  enrich_compounds = FALSE,
  cache = TRUE,
  cache_dir = "uafR_plant_cache"
)
```

Raw LOTUS MongoDB downloads can be converted with
`tools/flatten_lotus_mongo_dump.py`. The source URL, checksums, build manifest,
and retrieval date should be retained with the resulting index.

### Review evidence and identity

Plant evidence remains separable by source and strength:

```r
graded <- plantOccurrenceEvidenceGrade(
  phyto$PlantCompoundOccurrences
)

direct <- filterPlantEvidenceDirect(graded)
review_required <- filterPlantEvidenceReviewRequired(graded)

identity_review <- plantCompoundIdentityReviewTable(
  phyto
)
```

Literature co-mentions are candidates, not confirmed occurrence records.
Genus/family fallback remains separate from direct species evidence. Review
decisions can be exported and replayed with the plant-evidence and
compound-identity review helpers.

### Separate comparable chemistry

Do not mix primary metabolites, specialized metabolites, volatiles, lipids, and
unknown chemistry without an explicit scientific reason:

```r
comparability <- plantChemistryComparability(plant_compounds)

comparable_matrix <- plantComparableChemistryMatrix(
  comparability,
  level = "species",
  mode = "count"
)

chemistryComparisonDictionary()
```

Unknown, unresolved, and non-comparable compounds remain available for audit
but are excluded from comparable matrices by default.

### Calculate plant-labeled structural similarity

Tanimoto workflows use source-backed, structure-resolved compounds:

```r
similarity <- plantChemicalTanimotoSimilarity(
  phyto,
  out_dir = "plant-tanimoto",
  cache_dir = "uafR_plant_cache/pubchem",
  return_group_compound_pairs = FALSE
)

similarity$PlantPairTanimotoSummary
similarity$ComparableScopeTanimotoSummary
similarity$ComparableGroupTanimotoSummary
```

Structural similarity is not biological equivalence, shared function, pathway
activity, or evidence of occurrence in a project sample.

### Plan and run a large panel

Preflight does not make live requests:

```r
run_plan <- planPlantChemistryRun(
  plants = "project_species.csv",
  sources = c("lotus", "npass", "knapsack", "pubchem", "pubmed", "pubtator"),
  cache_dir = "uafR_plant_cache",
  lotus_index = "uafR_indexes/LOTUS_lookup_index",
  species_chunk_size = 25,
  compound_batch_size = 25,
  max_pubmed_records = 25
)

run_plan$InputNameAudit
run_plan$ProviderPlan
run_plan$OutputEstimates
run_plan$ReadinessChecks
run_plan$Recommendations
```

Use `runPlantChemistryPanel()` or its command-line wrapper for an audited,
resumable production run:

```sh
Rscript tools/run_plant_chemistry_panel.R \
  --mode preflight \
  --plant-csv project_species.csv \
  --out-dir plant_phytochemistry_panel \
  --cache-dir uafR_plant_cache \
  --lotus-index uafR_indexes/LOTUS_lookup_index \
  --sources lotus,npass,knapsack,pubchem,pubmed,pubtator
```

Run `pilot` before `discovery`. Exit status 75 means a public service requested
a pause; rerun the same command later to resume from validated checkpoints and
caches. Full pairwise compound output is opt-in because it can become very
large.

### Interpret plant results conservatively

- A database record supports a reported association, not measured chemistry in
  the current sample.
- No returned record is a coverage limitation, not evidence of absence.
- Direct species records are stronger than genus or family fallback.
- PubMed/PubTator co-mentions remain candidate evidence unless reviewed.
- Pathway annotations provide context, not evidence of pathway activity.
- Natural-product occurrence in another organism does not establish occurrence
  in the queried plant.
- Tanimoto similarity describes structure, not efficacy or mechanism.

## Source-specific profiles

### PubChem

```r
pubchem <- pubchemProfile(
  query_compounds,
  profile = "ms",
  cache = TRUE
)

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

Standardize `mzExacto()` output with an internal standard:

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
| `resolvePlantPhytochemistry()` | Discover reported species-first plant chemistry evidence |
| `buildLotusIndex()` / `queryLotusIndex()` | Build and query scalable local LOTUS resources |
| `runPlantPhytochemistryBatch()` | Run resumable species discovery in chunks |
| `runPlantChemistryPanel()` | Coordinate a staged production plant panel |
| `plantOccurrenceEvidenceGrade()` | Assign conservative, reviewable evidence grades |
| `plantChemistryComparability()` | Separate biologically comparable chemistry scopes |
| `plantChemicalTanimotoSimilarity()` | Calculate plant-labeled structure similarity |
| `runPlantChemistryProject()` | Create a reusable local plant chemistry project |
| `exportPlantChemistryFeatureSet()` | Write analysis-neutral species feature matrices |
| `validatePlantChemistryAnalysisBundle()` | Audit final bundle schemas and references |
| `exportAiNsectMolOlfInputs()` | Export structure-resolved molecular-olfaction inputs |

## Reproducibility and network access

Several workflows query public web services. Results can change as external
databases evolve.

- Keep `cache = TRUE` for reusable PubChem and KEGG responses.
- Use local LOTUS/NPASS indexes for medium or large species panels.
- Store result objects and validation tables with each analysis.
- Record the uafR version and session information.
- Use bounded request limits and the default throttling settings.
- Inspect `SourceCoverage`, `SourceDiagnostics`, and validation output before
  interpreting missing values.
- Treat HTTP 429/503 responses as a pause signal; resume from the same cache
  and output paths rather than restarting.
- Never write NCBI API keys or other credentials into scripts, command-line
  arguments, manifests, or exported data.

Record the environment used for an analysis:

```r
packageVersion("uafR")
sessionInfo()
```

## Documentation and support

- Function reference: <https://castrattonDSU.github.io/uafR/>
- Public training materials: [`training/`](training/)
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
