# uafR Student Quick Reference

Keep this file open while setting up uafR and running the first training
workflows.

## First-Time Setup

1. Unzip the whole bundle folder.
2. Open RStudio.
3. Open `START_HERE.md`.
4. Source the preflight check:

```r
source("preflight_check.R")
```

5. If preflight has no failed checks, install uafR:

```r
source("install_uafR_from_bundle.R")
```

6. Run the student acceptance test:

```r
source("run_student_acceptance_test.R")
```

Success means the Console ends with:

```text
Acceptance test result: PASS.
```

## Check the Bundle Files

To confirm the bundle was not damaged while copying or unzipping:

```r
source("verify_bundle_integrity.R")
```

Success means the Console ends with:

```text
Bundle integrity result: PASS.
```

## Load uafR

```r
library(uafR)
data("library_data", package = "uafR")
```

## Offline Training Workflow

This does not query live web services.

```r
source("training/scripts/04_core_workflow.R")

result <- run_core_workflow(live_lookup = FALSE)
result$exact
```

## Small Live Database Test

Run this only after the offline checks pass and internet access is available.

```r
library(uafR)
data("library_data", package = "uafR")

quick_result <- categorate(
  compounds = c("aspirin", "caffeine"),
  chemical_library = library_data,
  input_format = "wide",
  detail = "research",
  cache = TRUE,
  throttle = 0.2,
  assay_detail_limit = 0
)

validateCategorateResult(quick_result)$Summary
quick_result$ChemicalTraitReport
quick_result$ChemicalMeasurementSummary
```

## Offline Species-First Plant Extension

This example uses simulated curated rows. It demonstrates the package workflow;
it is not evidence about a research sample.

```r
source("training/scripts/10_species_phytochemistry.R")

plant_example <- run_species_phytochemistry_smoke(
  output_dir = "results/plant_training_example"
)

plant_example$SpeciesChemistrySummary
plant_example$PlantCompoundOccurrences
```

## Update From a New Bundle

When a new bundle is distributed, unzip the new bundle and run:

```r
source("update_uafR_from_bundle.R")
source("run_student_acceptance_test.R")
```

## Common Problems

| Problem | What to do |
|---|---|
| A script cannot find `packages/uafR_<version>.tar.gz`. | Unzip the whole bundle again and keep the folder structure unchanged. |
| `preflight_check.R` says the R library is not writable. | Restart RStudio. If it still fails, ask for help setting a user library. |
| Bioconductor dependencies fail. | Confirm internet access, update R if it is old, and rerun `preflight_check.R`. |
| `library(uafR)` fails right after install. | Restart RStudio, then run `library(uafR)` again. |
| `categorate()` says no library was detected. | Load `library_data` and pass `chemical_library = library_data`. |
| A live query is slow. | Start with one or two compounds, use `cache = TRUE`, and keep `assay_detail_limit = 0`. |

## Files to Look For

| File | Purpose |
|---|---|
| `START_HERE.md` | First setup checklist. |
| `README_STUDENT_INSTALL.md` | Detailed installation and troubleshooting directions. |
| `uafR_QUICK_REFERENCE.md` | This short command reference. |
| `CHECKSUMS.csv` | File sizes and MD5 checksums for the generated bundle. |
| `training/uafR_training_manual.pdf` | Full 10-week training manual. |
| `uafR_*_log.txt` | Logs written after setup scripts run. |
