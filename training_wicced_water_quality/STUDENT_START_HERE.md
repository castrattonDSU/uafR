# Student Start Here: WiCCED Water-Quality uafR/ML Training Bundle

This bundle contains the 8-week Project WiCCED water-quality training manual,
simulated teaching data, R scripts, worksheets, rubrics, machine-learning
resources, and optional reference outputs. The simulated data are for training
only. They are not real WiCCED measurements and must not be used as
environmental evidence.

## 1. Unzip and open the folder

1. Save the zipped training bundle to a local folder on your computer.
2. Unzip it.
3. Move the unzipped folder somewhere easy to find, such as Desktop or
   Documents. Avoid editing the training folder from inside the zip preview.
4. Open RStudio.
5. In RStudio, choose `Session > Set Working Directory > Choose Directory`.
6. Select the unzipped training folder, not one of its subfolders.

Run this in the RStudio Console to confirm that RStudio is pointed at the
right folder:

```r
getwd()
file.exists("main.pdf")
file.exists("data/wicced_teaching_water_quality.csv")
file.exists("scripts/run_wicced_smoke_tests.R")
```

All three `file.exists()` checks should return `TRUE`.

## 2. Open the manual first

Open `main.pdf`. This is the full 8-week training manual. Keep it open while
you run the scripts. The manual explains what each week produces, how to
interpret the outputs, and what claims should not be made from simulated or
tentative evidence.

## 3. Check the R session

Run the setup check:

```r
source("scripts/00_setup_check.R")
```

This writes `setup_package_status.csv` in the training folder. If `uafR` is not
available, that is not an emergency. The core WiCCED teaching scripts are
offline and use simulated data. The uafR package is needed for live or expanded
uafR work.

## 4. Install uafR if a local package archive is included

Some bundles include a `uafR_package/` folder containing a local uafR source
archive. If it is present, install it with:

```r
pkg <- list.files(
  "uafR_package",
  pattern = "^uafR_.*[.]tar[.]gz$",
  full.names = TRUE
)

if (length(pkg) == 0) {
  message("No bundled uafR package archive found.")
} else {
  install.packages(pkg[[1]], repos = NULL, type = "source")
}
```

After installation, run the setup check again:

```r
source("scripts/00_setup_check.R")
```

If installation fails, save the full error message in your research log and
continue with the offline training scripts until the issue is reviewed.

## 5. Run the full smoke test

The smoke test builds a temporary project and runs the setup, water-quality
QA/QC, simulated GC-MS workflow, beginner machine-learning ladder, and advanced
machine-learning extension. It does not use live web services.

```r
source("scripts/run_wicced_smoke_tests.R")
```

The final line should say:

```text
WiCCED training smoke test passed.
```

If the smoke test fails, copy the error message into `logs/daily_research_log.md`
and stop before changing scripts.

## 6. Create your working folders

Run:

```r
source("scripts/01_project_setup.R")
```

This creates folders such as `results/`, `models/`, `figures/`, `logs/`,
`reports/`, and `references/`. The file `logs/daily_research_log.md` is the
main place to record commands, outputs, decisions, errors, and interpretation
limits.

## 7. Run the beginner workflow in order

Run each script from the RStudio Console while the training folder is the
working directory.

```r
source("scripts/02_water_quality_qaqc.R")
source("scripts/03_uafr_chemical_screening_template.R")
source("scripts/04_ml_01_build_feature_table.R")
source("scripts/04_ml_02_explore_features.R")
source("scripts/04_ml_03_train_models.R")
source("scripts/04_ml_04_diagnostics_and_model_card.R")
```

Expected output locations:

- `results/` contains water-quality QA/QC outputs and simulated GC-MS outputs.
- `results/ml/` contains feature tables, train/test split files, metrics,
  predictions, plots, diagnostics, and a model card.
- `logs/daily_research_log.md` should be updated after each work session.

The beginner machine-learning scripts use simple baseline models so that the
workflow is understandable. Do not skip the beginner steps.

## 8. Run the advanced extension only after the beginner steps make sense

The advanced extension adds blocked validation, feature-set sensitivity,
permutation importance, decision gates, and an advanced model review. It should
be run only after the beginner feature table, leakage exclusions, train/test
split, metrics, prediction errors, and model card are understood.

```r
source("scripts/05_ml_01_advanced_validation.R")
source("scripts/05_ml_02_permutation_importance_and_sensitivity.R")
source("scripts/05_ml_03_advanced_model_report.R")
```

Expected output location:

- `results/ml_advanced/`

## 9. Compare against reference outputs if included

If the bundle contains an `expected_outputs/` folder, use it only as a
reference. Your own outputs should be created by running the scripts. Do not
edit the reference outputs.

## 10. What to bring to the first review

Bring these items:

- `setup_package_status.csv`
- `logs/daily_research_log.md`
- `results/water_quality_qaqc_flags.csv`
- `results/wicced_simulated_gcms_pipeline_summary.txt`
- `results/ml/ml_feature_dictionary.csv`
- `results/ml/ml_model_comparison_summary.md`
- `results/ml/ml_model_card.md`
- Any error messages that appeared while running scripts

Be prepared to explain what was simulated, what was measured, what was modeled,
and what claims are not supported.
