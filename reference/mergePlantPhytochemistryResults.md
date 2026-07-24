# Merge plant phytochemistry discovery results

Deterministically combines separately completed provider or species
batches. Occurrence evidence is preserved, exact duplicate evidence keys
are removed, and per-species provider accounting prefers completed
record/no-record states over earlier retry-required states.

## Usage

``` r
mergePlantPhytochemistryResults(..., strict = FALSE)
```

## Arguments

- ...:

  Plant phytochemistry result objects, or one list containing them.

- strict:

  Logical passed to result validation.

## Value

A merged \`uaf_plant_phytochemistry\` result.
