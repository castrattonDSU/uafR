# Prepare plant chemistry for structure-based Tanimoto analysis

Builds a conservative, offline handoff between species-first plant
chemistry discovery and PubChem Fingerprint2D Tanimoto analysis.
Evidence rows remain available for provenance, while the analysis
membership is collapsed to one row per plant and source-backed chemical
identity. Repeated publications, source records, and aliases are
summarized rather than treated as distinct compounds. Exact
source-record identities are preferred when \`SourceCompoundIdentity\`
is available. Any inherited fingerprint columns are deliberately
removed. No fingerprints or pairwise similarities are computed.

Valid full InChIKeys are preferred as identity keys. A positive PubChem
CID is used only when an InChIKey is unavailable. When an InChIKey is
present, any existing CID is retained only as reported audit metadata
and is not passed as trusted input to the next PubChem stage. This
forces the downstream workflow to resolve and verify the InChIKey with
the collision-resistant uafR cache.

\`DuplicateAudit\` reports repeated evidence collapsed within one
plant-structure membership. \`NameStructureAudit\` separately reports
normalized compound labels associated with multiple exact structures.
Those structures remain distinct and require label-level review; they
are never collapsed merely because their names match.

## Usage

``` r
preparePlantTanimotoInput(
  plant_chemistry,
  level = "species",
  occurrence_status = c("direct_reported", "curated_reported"),
  analysis_ready = TRUE,
  min_confidence = "medium",
  include_review_required = FALSE,
  out_dir = NULL,
  overwrite = FALSE,
  strict = TRUE
)
```

## Arguments

- plant_chemistry:

  A plant phytochemistry result containing \`PlantCompoundOccurrences\`
  and \`CompoundResolution\`, or a data frame that already combines
  plant labels with source-backed identity fields.

- level:

  Grouping column, normally \`"species"\`.

- occurrence_status:

  Allowed occurrence statuses when that column is present. Defaults to
  direct or curated reported evidence.

- analysis_ready:

  If \`TRUE\`, require analysis-ready evidence when the input contains
  an \`analysis_ready\` column.

- min_confidence:

  Minimum evidence confidence when a \`confidence\` column is present.

- include_review_required:

  If \`FALSE\`, identity rows flagged for review are excluded from the
  Tanimoto-ready membership and retained in the excluded table.

- out_dir:

  Optional directory for the CSV/JSON server handoff bundle.

- overwrite:

  Logical. If \`FALSE\`, existing handoff files are protected.

- strict:

  Logical. If \`TRUE\`, stop when a validation check fails.

## Value

A list with \`PlantCompoundMembership\`, \`CompoundInput\`,
\`EvidenceRows\`, \`DuplicateAudit\`, \`NameStructureAudit\`,
\`ExcludedRows\`, \`ValidationSummary\`, \`Summary\`,
\`ExportManifest\`, and \`Provenance\`.

## Examples

``` r
if (FALSE) { # \dontrun{
prepared = preparePlantTanimotoInput(
  phyto,
  out_dir = "plant_tanimoto_handoff"
)
prepared$PlantCompoundMembership
} # }
```
