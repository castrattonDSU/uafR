# Machine-Learning Troubleshooting

## The model runs but the result looks too good

Check for leakage. Make sure the feature table does not include the target
itself, a direct copy of the target, a post-outcome interpretation, or a site/date
identifier that memorizes the answer.

## The model performs worse than the baseline

That is a valid result. Record it. A weak model can mean the dataset is too small,
the features are not informative, the split is difficult, or the target is not
well defined.

## The logistic model gives extreme probabilities

Small teaching datasets can separate classes too cleanly. Report the behavior
and avoid strong claims. Use the confusion matrix and model-card limitations
rather than claiming deployment readiness.

## The train/test split feels unfair

A hard split is often honest. If the goal is future prediction, a time split can
be more realistic than a random split. If the goal is new-site prediction, a
site-holdout split may be appropriate. The split must match the claim.

## A feature has missing values

Do not silently delete data. Record missingness by variable and group. Decide
whether to exclude the feature, impute with a documented rule, or keep the model
limited to complete rows.

## A figure or metric does not support the story

Revise the story, not the evidence. The final report should reflect what the
model actually did under the recorded split.

