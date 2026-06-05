# Project WiCCED Water Quality uafR Training Manual

This folder contains a separate 8-week undergraduate training manual focused on
Project WiCCED, water quality, uafR, reproducible R workflows, and meaningful
machine-learning practice.

Build from this folder:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The rendered PDF is `main.pdf`.

Run the offline script smoke test from the repository root:

```sh
Rscript training_wicced_water_quality/scripts/run_wicced_smoke_tests.R
```

The smoke test creates a temporary project, runs setup, QA/QC, simulated GC-MS,
the full beginner ML ladder, and verifies expected outputs without live database
access.

Run the simulated GC-MS/uafR-compatible teaching pipeline from this folder:

```sh
Rscript scripts/03_uafr_chemical_screening_template.R
```

The simulated pipeline uses clearly labeled teaching data in `data/` and writes
ignored outputs to `results/`. These files demonstrate how a GC-MS hit table can
be validated, spread, aggregated into exact-match abundance tables, converted to
relative abundance, and summarized by water-quality context. They are not real
WiCCED measurements.

Run the beginner machine-learning ladder from this folder:

```sh
Rscript scripts/04_ml_01_build_feature_table.R
Rscript scripts/04_ml_02_explore_features.R
Rscript scripts/04_ml_03_train_models.R
Rscript scripts/04_ml_04_diagnostics_and_model_card.R
```

Those scripts write beginner-friendly tables, plots, metrics, diagnostics, and a
model card to `results/ml/`. Guiding resources are in `ml_resources/`.

After the beginner ladder is understood, run the advanced ML extension:

```sh
Rscript scripts/05_ml_01_advanced_validation.R
Rscript scripts/05_ml_02_permutation_importance_and_sensitivity.R
Rscript scripts/05_ml_03_advanced_model_report.R
```

Those scripts write blocked-validation, feature-set sensitivity, permutation
importance, decision-gate, and advanced model-review outputs to
`results/ml_advanced/`. The advanced extension should not be interpreted until a
student can explain the beginner feature table, leakage exclusions, baseline
comparison, prediction errors, and model card.

The manual is intentionally separate from the main 10-week dsDNA Core uafR
training program in `training/`.
