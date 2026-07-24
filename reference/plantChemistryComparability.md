# Classify plant chemistry into comparable analysis scopes

\`plantChemistryComparability()\` adds a conservative comparison layer
to plant-compound evidence. The output separates biological role
(\`metabolism_domain\`), biosynthetic or chemical family
(\`biosynthetic_family\`), and analytical/physicochemical behavior
(\`chemical_behavior\`). This prevents downstream analyses from mixing
unmatched concepts such as primary metabolites, specialized metabolites,
lipids, hormone signals, and volatile fractions.

The classifier is deterministic and source-transparent. It prioritizes
source-backed uafR class fields from \`ChemicalClasses\`,
\`LOTUSProfile\`, \`PubChemClassifications\`, \`ChemicalTraitOntology\`,
\`ChemicalTerms\`, \`KEGGClassifications\`, \`KEGGPathways\`, and
\`DerivedGroups\`, then falls back to
\`ChemicalTraits\`/\`ChemicalClasses\` aggregates, compound-name
patterns, and occurrence evidence text. Unknown compounds remain
\`unknown\` and are excluded from comparable matrices unless explicitly
requested.

## Usage

``` r
plantChemistryComparability(
  x,
  categorate_result = NULL,
  min_confidence = "low"
)
```

## Arguments

- x:

  Plant phytochemistry result, normalized occurrence table, or named
  list containing \`PlantCompoundOccurrences\` and optional
  \`CategorateResult\`.

- categorate_result:

  Optional categorate-like enrichment result when \`x\` is an occurrence
  table.

- min_confidence:

  Minimum occurrence confidence to include.

## Value

A \`ChemistryComparability\` data frame with one row per retained
plant-compound occurrence.
