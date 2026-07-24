# Audit provider-specific biological-context coverage

\`plantProviderContextAudit()\` summarizes how much source-backed
plant-part, tissue, and method context each provider contributes. It is
designed for pilot runs and scale-up decisions: sparse or low-confidence
providers are flagged for source inspection, parser hardening, or manual
curation before context-aware matrices are interpreted.

## Usage

``` r
plantProviderContextAudit(x, context_evidence = NULL)
```

## Arguments

- x:

  Plant phytochemistry result, normalized occurrence table, or named
  list containing \`PlantCompoundOccurrences\`.

- context_evidence:

  Optional \`PlantContextEvidence\` table. If omitted, context evidence
  is taken from \`x\` when available or extracted from occurrences.

## Value

A \`ProviderContextAudit\` data frame with one row per source database.
