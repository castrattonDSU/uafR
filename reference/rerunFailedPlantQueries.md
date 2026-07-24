# Rerun failed plant chemistry query batches

Provides a generic retry runner for queue rows produced by
\`writePlantChemistryRetryQueue()\`. By default \`dry_run = TRUE\`, so
the function only returns an execution plan. To actually rerun work,
provide a \`runner_fun\` that accepts one retry-queue row as its first
argument and set \`dry_run = FALSE\`.

## Usage

``` r
rerunFailedPlantQueries(
  retry_queue,
  runner_fun = NULL,
  dry_run = TRUE,
  out_dir = NULL,
  overwrite = FALSE,
  ...
)
```

## Arguments

- retry_queue:

  Retry queue data frame, CSV path, or manifest-like object.

- runner_fun:

  Optional function called for each retry row when \`dry_run = FALSE\`.

- dry_run:

  Logical. If \`TRUE\`, do not execute retries.

- out_dir:

  Optional directory for \`retry_run_manifest.csv\`.

- overwrite:

  Logical. If \`FALSE\`, existing retry manifests are preserved.

- ...:

  Additional arguments passed to \`runner_fun\`.

## Value

List with \`Plan\`, \`RetryQueue\`, and \`RetryRunManifest\`.
