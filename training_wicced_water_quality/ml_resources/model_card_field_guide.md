# Model Card Field Guide

## Project Question

Write the specific prediction question. Example: "Can basic water-quality
measurements predict high-salinity class in the synthetic teaching data?"

## Target

Name the target column and task type. Regression predicts a number.
Classification predicts a category.

## Rows

Record how many rows were used for training and testing. If rows were excluded,
state why.

## Features

List the features used by the model. For each feature, record source, unit,
transformation, and leakage status.

## Split Design

State whether the split was time-based, random, site-holdout, or another design.
Explain what kind of generalization it tests.

## Baseline

State the simple model used for comparison. For regression, this is often the
training mean. For classification, this is often the majority class.

## Metrics

Use RMSE and MAE for regression. Use confusion-matrix counts, accuracy,
sensitivity, and specificity for classification.

## Diagnostics

Inspect observed-versus-predicted plots, residuals, and class predictions. State
where the model fails.

## Limits

Say what the model cannot support. For this training program, all teaching-data
models are for workflow practice and hypothesis generation only.

