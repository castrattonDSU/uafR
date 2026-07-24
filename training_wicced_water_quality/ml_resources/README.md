# Machine-Learning Resources

This folder supports the machine-learning part of the WiCCED water-quality
training manual. The goal is to help students understand each modeling step
before they try more advanced tools.

Run the beginner ML ladder from the `training_wicced_water_quality/` folder:

```sh
Rscript scripts/04_ml_01_build_feature_table.R
Rscript scripts/04_ml_02_explore_features.R
Rscript scripts/04_ml_03_train_models.R
Rscript scripts/04_ml_04_diagnostics_and_model_card.R
```

The scripts use only base R and the synthetic teaching water-quality dataset by
default. Outputs are written to `results/ml/`, which is ignored by git. The
teaching data are synthetic and cannot support environmental claims.

## Beginner Path

1. Build the feature table.
   Students learn what the target is, what the features are, what is excluded,
   and how the train/test split is made.

2. Explore before modeling.
   Students inspect ranges, missingness, simple correlations, class balance, and
   plots before fitting any model.

3. Train baseline models.
   Students fit a mean-only baseline, a linear regression model, a majority-class
   classifier, and a logistic classifier. The point is comparison, not maximum
   accuracy.

4. Diagnose and document.
   Students inspect prediction plots, residuals, confusion-matrix counts, model
   coefficients, and a model card.

## Advanced Extension

Use the advanced extension only after the beginner path is understood and the
beginner model card is complete:

```sh
Rscript scripts/05_ml_01_advanced_validation.R
Rscript scripts/05_ml_02_permutation_importance_and_sensitivity.R
Rscript scripts/05_ml_03_advanced_model_report.R
```

The advanced extension adds blocked validation, feature-set sensitivity,
permutation importance, and decision gates. Start with
`advanced_ml_extension_guide.md` before interpreting the outputs.

## Required Interpretation Rule

Model output is evidence about a pattern under a recorded split. It is not proof
of causation, regulatory status, site safety, or field deployment readiness.
