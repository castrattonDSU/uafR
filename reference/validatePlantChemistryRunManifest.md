# Validate a plant chemistry run manifest

Standardizes and validates a batch/run manifest from large plant
chemistry projects. The function is offline-only: it checks table shape,
status values, and output-file existence when paths are supplied, but it
never reruns failed queries or contacts providers.

## Usage

``` r
validatePlantChemistryRunManifest(manifest, base_dir = NULL)
```

## Arguments

- manifest:

  Data frame, CSV path, JSON path, or list containing a manifest table.

- base_dir:

  Optional directory used to resolve relative output paths.

## Value

A list with \`Summary\`, \`TableQuality\`, \`Issues\`, \`RetryQueue\`,
and standardized \`Manifest\` tables.
