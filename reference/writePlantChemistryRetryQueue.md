# Write a plant chemistry retry queue

Extracts failed, incomplete, timed-out, rate-limited, or not-started
batches from a large-run manifest into a reproducible retry queue. The
queue can be inspected, edited, archived, or passed to
\`rerunFailedPlantQueries()\`.

## Usage

``` r
writePlantChemistryRetryQueue(
  manifest,
  path = NULL,
  include_status = c("failed", "error", "timeout", "timed_out", "rate_limited",
    "incomplete", "planned", "not_started", "running", "started", "stopped", "retry"),
  overwrite = FALSE
)
```

## Arguments

- manifest:

  Data frame, CSV path, JSON path, or list containing a manifest table.

- path:

  Optional CSV path to write.

- include_status:

  Batch statuses to include in the retry queue.

- overwrite:

  Logical. If \`FALSE\`, an existing \`path\` is not replaced.

## Value

Retry queue data frame.
