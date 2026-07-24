# Build matrices from comparable plant chemistry scopes

Builds wide matrices only after selecting a defensible chemistry
comparison scope. Use this for clustering, ordination, candidate
scoring, and figures when chemistry should be compared like-for-like.
For example, \`"volatile_specialized_metabolites"\` compares volatile
specialized chemistry, while \`"primary_metabolites"\` compares
primary-metabolism compounds.

## Usage

``` r
plantComparableChemistryMatrix(
  x,
  level = c("species", "genus", "family"),
  comparison_scope = "specialized_metabolites",
  feature = c("comparison_group", "biosynthetic_family", "chemical_behavior", "compound",
    "comparison_scope"),
  mode = c("binary", "count", "confidence"),
  min_comparability_confidence = "medium",
  plant_part_group = NULL,
  tissue_group = NULL,
  method_group = NULL,
  biological_context_status = NULL,
  require_context = FALSE,
  include_unknown = FALSE,
  max_features = Inf
)
```

## Arguments

- x:

  Plant phytochemistry result, \`ChemistryComparability\` table, or
  named list containing \`ChemistryComparability\`.

- level:

  One of \`"species"\`, \`"genus"\`, or \`"family"\`.

- comparison_scope:

  Scope(s) to retain. Use \`"all_classified"\` to retain every
  classified non-unknown scope.

- feature:

  Feature family to turn into columns: \`"comparison_group"\`,
  \`"biosynthetic_family"\`, \`"chemical_behavior"\`, \`"compound"\`, or
  \`"comparison_scope"\`.

- mode:

  One of \`"binary"\`, \`"count"\`, or \`"confidence"\`.

- min_comparability_confidence:

  Minimum classification confidence.

- plant_part_group:

  Optional plant-part groups to retain.

- tissue_group:

  Optional tissue groups to retain.

- method_group:

  Optional method groups to retain.

- biological_context_status:

  Optional context-status values to retain.

- require_context:

  Logical. If \`TRUE\`, retain only rows with known plant part/tissue or
  analytical method context.

- include_unknown:

  Logical. If \`TRUE\`, unknown comparison groups can be retained.

- max_features:

  Maximum number of feature columns to retain.

## Value

Wide matrix-like data frame.
