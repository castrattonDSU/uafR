# Pairwise PubChem fingerprint Tanimoto similarity for chemicals

\`chemicalTanimotoSimilarity()\` resolves chemicals to PubChem CIDs,
fetches PubChem Fingerprint2D bit strings, decodes them with
\`ChemmineR\`, and computes pairwise Tanimoto similarity. It accepts a
plain character vector or a data frame with compound identifiers. When a
grouping column is supplied, it also returns group-level summaries that
are suitable for joining to other sample, species, treatment, or site
similarity tables.

## Usage

``` r
chemicalTanimotoSimilarity(
  compounds,
  compound_col = NULL,
  cid_col = NULL,
  inchikey_col = NULL,
  smiles_col = NULL,
  compound_id_col = NULL,
  group_cols = NULL,
  thresholds = c(0.5, 0.7, 0.85, 0.95),
  top_n_pairs = 5,
  return_compound_pairs = TRUE,
  return_group_compound_pairs = FALSE,
  out_dir = NULL,
  max_in_memory_pairs = 1e+06,
  pair_block_size = 250,
  pair_shard_rows = Inf,
  cache = TRUE,
  cache_dir = NULL,
  throttle = 0.2,
  refresh = FALSE,
  request_fun = NULL,
  name_fallback = TRUE,
  progress_fun = NULL,
  progress_every = 25,
  service_busy_limit = Inf
)
```

## Arguments

- compounds:

  Character vector or data frame. Character vectors are treated as
  compound names. Data frames can contain compound names, PubChem CIDs,
  InChIKeys, SMILES, and optional grouping columns.

- compound_col:

  Column containing compound names when \`compounds\` is a data frame.
  If \`NULL\`, common names such as \`compound_name\`, \`Query\`, or
  \`Chemical\` are detected.

- cid_col:

  Column containing PubChem CIDs. If \`NULL\`, common names such as
  \`CID\`, \`pubchem_cid\`, or \`cid\` are detected.

- inchikey_col:

  Column containing InChIKeys. If \`NULL\`, common names such as
  \`InChIKey\` or \`inchikey\` are detected.

- smiles_col:

  Optional SMILES column retained as metadata. PubChem Fingerprint2D
  values are still fetched from PubChem so the Tanimoto definition is
  consistent across workflows.

- compound_id_col:

  Optional stable compound identifier column. If absent, identifiers are
  derived from InChIKey, PubChem CID, or compound name.

- group_cols:

  Optional grouping columns. For plants, use \`"species"\` or call
  \`plantChemicalTanimotoSimilarity()\`.

- thresholds:

  Numeric Tanimoto thresholds to count in group summaries.

- top_n_pairs:

  Number of top compound-pair explanations to include in group
  summaries.

- return_compound_pairs:

  Logical. If \`TRUE\`, include or write the compound-compound pair
  table.

- return_group_compound_pairs:

  Logical. If \`TRUE\` and \`group_cols\` are supplied, include or write
  all cross-group compound-pair rows.

- out_dir:

  Optional output directory. When supplied, large pairwise tables are
  written as compressed CSV files and the result records their paths.

- max_in_memory_pairs:

  Maximum pair rows allowed in memory when \`out_dir\` is \`NULL\`.

- pair_block_size:

  Number of focal rows per pairwise computation block.

- pair_shard_rows:

  Maximum rows per compressed pairwise CSV shard. Use \`Inf\` (the
  default) to preserve one file per requested pair table.

- cache:

  Logical. If \`TRUE\`, PubChem JSON responses are cached.

- cache_dir:

  Cache directory. Defaults to the uafR user cache.

- throttle:

  Seconds to wait after uncached PubChem requests.

- refresh:

  Logical. If \`TRUE\`, ignore cached PubChem responses.

- request_fun:

  Optional request function for tests. It receives a URL and returns
  parsed JSON or JSON text.

- name_fallback:

  Logical. If \`TRUE\`, unresolved compounds may be queried by compound
  name after CID and InChIKey attempts fail.

- progress_fun:

  Optional callback receiving named progress-event lists. It is intended
  for durable project runners and is ignored by default.

- progress_every:

  Emit identity and pair-block progress after this many completed items
  or blocks.

- service_busy_limit:

  Number of consecutive PubChem HTTP 429/503 responses allowed before a
  \`uaf_pubchem_service_busy\` condition stops the run. The default
  \`Inf\` preserves ordinary interactive behavior.

## Value

A list with class \`"uaf_tanimoto_similarity"\` containing
\`CompoundInput\`, \`CompoundResolution\`, \`PubChemFingerprints\`,
\`CompoundIdentityMap\`, \`ExcludedCompoundIdentities\`,
\`CanonicalCompounds\`, \`GroupCompoundMembershipEvidence\`,
\`CompoundTanimoto\`, \`GroupCompoundMembership\`,
\`GroupPairTanimotoSummary\`, \`GroupCompoundTanimoto\`,
\`ComparableScopeTanimotoSummary\`, \`ComparableGroupTanimotoSummary\`,
\`ExportManifest\`, and \`Provenance\`. The identity tables retain every
input mapping while analysis membership is collapsed to verified
canonical PubChem structures. Comparable summaries require source-backed
comparison labels and exclude unknown or explicitly non-comparable
membership by default.

## Examples

``` r
if (FALSE) { # \dontrun{
compounds = data.frame(
  species = c("Plant A", "Plant A", "Plant B"),
  compound_name = c("caffeine", "theobromine", "aspirin")
)
sim = chemicalTanimotoSimilarity(compounds, group_cols = "species")
sim$GroupPairTanimotoSummary
} # }
```
