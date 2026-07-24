# Advanced ML Extension Guide

## Use This Only After the Beginner Ladder

The advanced extension should be used only after the student can explain:

- What the target is.
- Which variables are predictors.
- Which variables were excluded and why.
- What leakage means.
- What the train/test split evaluates.
- Whether the beginner model beat the baseline.
- Which rows were predicted poorly.
- What the beginner model card says the model cannot do.

If those points are not clear, rerun and discuss scripts `04_ml_01` through
`04_ml_04` before opening the advanced outputs.

## What the Advanced Extension Adds

The advanced scripts do not chase a more impressive model. They add stricter
checks:

1. Blocked validation.
   The model is tested with a time holdout and leave-one-site-out folds.

2. Feature-set sensitivity.
   The beginner feature set is compared with conductivity and context
   sensitivity feature sets. Conductivity is flagged because it is closely
   related to salinity.

3. Permutation importance.
   Each feature is shuffled in the held-out rows to test how much the fitted
   model depends on that feature under the chosen split.

4. Decision gate.
   The final report separates workflow success from scientific readiness.

## Script Order

Run the beginner scripts first:

```sh
Rscript scripts/04_ml_01_build_feature_table.R
Rscript scripts/04_ml_02_explore_features.R
Rscript scripts/04_ml_03_train_models.R
Rscript scripts/04_ml_04_diagnostics_and_model_card.R
```

Then run the advanced extension:

```sh
Rscript scripts/05_ml_01_advanced_validation.R
Rscript scripts/05_ml_02_permutation_importance_and_sensitivity.R
Rscript scripts/05_ml_03_advanced_model_report.R
```

## Main Outputs

- `results/ml_advanced/advanced_ml_validation_folds.csv`
- `results/ml_advanced/advanced_ml_split_comparison.csv`
- `results/ml_advanced/advanced_ml_feature_set_metrics.csv`
- `results/ml_advanced/advanced_ml_permutation_importance.csv`
- `results/ml_advanced/advanced_ml_sensitivity_warnings.csv`
- `results/ml_advanced/advanced_ml_sensitivity_plots.pdf`
- `results/ml_advanced/advanced_ml_decision_gate.csv`
- `results/ml_advanced/advanced_ml_final_recommendations.csv`
- `results/ml_advanced/advanced_ml_model_review.md`

## How to Read the Outputs

Start with `advanced_ml_model_review.md`. It gives the shortest interpretation
of what changed after the beginner model.

Then open `advanced_ml_decision_gate.csv`. Any row marked
`requires_human_review`, `needs_work`, or `not_ready_for_environmental_claims`
must be discussed before the model is used in a report.

Open `advanced_ml_split_comparison.csv` next. A model that only looks good under
one split is not stable enough for strong claims.

Finally, open `advanced_ml_permutation_importance.csv`. Large importance values
show model dependence under the selected split. They do not prove cause,
mechanism, source, exposure, toxicity, risk, or management relevance.

## Required Interpretation Rule

An advanced model is acceptable only if it remains reproducible, source-backed,
leakage-reviewed, and scientifically limited. Higher accuracy is not enough.
