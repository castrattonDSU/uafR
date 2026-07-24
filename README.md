# uafR

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
