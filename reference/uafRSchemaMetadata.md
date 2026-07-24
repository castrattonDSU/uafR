# uafR schema metadata

Creates package/schema metadata rows for result objects, manifests, and
downstream bundles.

## Usage

``` r
uafRSchemaMetadata(
  workflow_name = "unspecified",
  workflow_parameters = NULL,
  schema_version = .uaf_schema_version()
)
```

## Arguments

- workflow_name:

  Name of the workflow creating the object.

- workflow_parameters:

  Optional named list of important parameters.

- schema_version:

  uafR schema version. Defaults to the current internal schema version.

## Value

One-row data frame with schema and package metadata.
