# Extract biological context evidence from plant occurrence rows

\`plantContextEvidence()\` extracts auditable plant-part, tissue, and
method context evidence from normalized plant-compound occurrence rows.
It preserves the raw matched text, normalized context group, extraction
rule, source field, source record identifiers, confidence, and review
flag. The result is used by plant phytochemistry workflows to separate
context-known rows from context-missing rows without fabricating
biological location.

## Usage

``` r
plantContextEvidence(x, min_confidence = "low")
```

## Arguments

- x:

  Plant phytochemistry result, normalized occurrence table, or named
  list containing \`PlantCompoundOccurrences\`.

- min_confidence:

  Minimum context confidence to retain.

## Value

A \`PlantContextEvidence\` data frame.
