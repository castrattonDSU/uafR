# Score plant chemistry candidates

Score plant chemistry candidates

## Usage

``` r
scorePlantChemistryCandidates(
  x,
  weights = c(direct_species = 0.25, resolved_compounds = 0.2, source_coverage = 0.2,
    literature = 0.15, traits = 0.2)
)
```

## Arguments

- x:

  Plant phytochemistry result or species summary table.

- weights:

  Named numeric vector for transparent score components.

## Value

Data frame with component scores and total prioritization score.
