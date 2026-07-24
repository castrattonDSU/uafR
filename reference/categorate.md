# categorate

Searches a vector of chemical names on PubChem for published
information. Uses regular expressions to extract information from 6
data- bases: LOTUS the natural products occurrence database
\[LOTUS\](https://lotus.naturalproducts.net/), Flavor and Extract
Manufacturers Association \[FEMA\](https://www.femaflavor.org/), Kyoto
Encyclopedia of Genes and Genomes \[KEGG\](https://www.genome.jp/kegg/),
Food and Drug Administration \[FDA/SPL\](https://www.fda.gov/), and
Reactive Groups from \[PubChem\](https://pubchem.ncbi.nlm.nih.gov/).
Also downloads structural data in SDF format to summarize the atomic
structures and common molecular groups of each chemical. Finally, the
function downloads structural data for the chemicals within each group
of the input library (\`chemical_library\`). Each chemical of interest
(\`compounds\`) is tested against every group of chemicals for
structural overlaps that exceed a threshold of 0.85.

## Usage

``` r
categorate(
  compounds,
  chemical_library,
  input_format = "wide",
  detail = c("standard", "research", "full"),
  cache = TRUE,
  cache_dir = NULL,
  throttle = 0.2,
  assay_detail_limit = 50,
  trait_matrix_profile = "core",
  trait_matrix_mode = "binary",
  trait_matrix_min_confidence = 0,
  trait_matrix_max_traits = Inf,
  request_fun = NULL,
  kegg_request_fun = NULL
)
```

## Arguments

- compounds:

  A vector containing chemical names in IUPAC notation (preferred).

- chemical_library:

  A data frame containing columns with chemical groups. Column names
  label the group while the rows contain chemicals that are described by
  the label. Use \`data("library_data", package = "uafR")\` to load the
  bundled example library.

- input_format:

  Designates the data structure for input library
  (\`chemical_library\`). Default setting is "wide" can be changed to
  "long."

- detail:

  One of \`"standard"\`, \`"research"\`, or \`"full"\`. \`"standard"\`
  returns the original eight output tables. \`"research"\` appends tidy
  PubChem and KEGG enrichment tables, source coverage, and derived
  grouping columns. \`"full"\` additionally requests broader PubChem
  annotations and assay summary data.

- cache:

  Logical. If \`TRUE\`, enrichment requests are cached. Only used when
  \`detail\` is \`"research"\` or \`"full"\`.

- cache_dir:

  Directory for cached enrichment responses. Defaults to a user-cache
  location when available.

- throttle:

  Seconds to wait between uncached enrichment requests.

- assay_detail_limit:

  Maximum number of unique PubChem BioAssay AIDs per query for which
  assay-description metadata is fetched when \`detail = "full"\`.
  Active, numeric, and target-bearing assays are prioritized. Set to
  \`0\` to skip assay-description requests.

- trait_matrix_profile:

  Trait profile used to build \`ChemicalTraitMatrix\`. Defaults to
  \`"core"\` for a compact cross-domain matrix. Use
  \`chemicalTraitMatrix(result, profile = "full")\` to rebuild wider
  matrices from \`ChemicalTraits\` without repeating database requests.

- trait_matrix_mode:

  Matrix value mode for \`ChemicalTraitMatrix\`: \`"binary"\` for 0/1
  presence, \`"count"\` for trait-row counts, or \`"confidence"\` for
  the maximum trait confidence score.

- trait_matrix_min_confidence:

  Minimum confidence required for traits in \`ChemicalTraitMatrix\`.
  Accepts a numeric score or \`"low"\`, \`"medium"\`, or \`"high"\`.

- trait_matrix_max_traits:

  Maximum number of trait columns to include in \`ChemicalTraitMatrix\`,
  ranked by prevalence and confidence. Defaults to \`Inf\`.

- request_fun:

  Optional PubChem request function for tests and advanced users. Only
  used when \`detail\` is \`"research"\` or \`"full"\`.

- kegg_request_fun:

  Optional KEGG request function for tests and advanced users. Only used
  when \`detail\` is \`"research"\` or \`"full"\`.

## Value

List with at least 8 data frames: \`reactives\`, \`LOTUS\`, \`KEGG\`,
\`FEMA\`, \`FDA_SPL\`, \`FMCS\`, \`FunctionalGroups\`, and
\`BestChemMatch\`. Database tables store extracted source annotations by
chemical; \`FMCS\` stores atomic and molecular sub-group summaries;
\`FunctionalGroups\` stores strong (similarity \> 0.95) or moderate
(similarity \> 0.85) matches with input library groups; and
\`BestChemMatch\` stores top chemicals from groups a chemical shared a
strong match with. With \`detail = "research"\` or \`detail = "full"\`,
the list also includes PubChem and KEGG enrichment tables:
\`PubChemIdentity\`, \`PubChemProperties\`, \`PubChemSynonyms\`,
\`PubChemAnnotations\`, \`PubChemSourceAnnotations\`,
\`PubChemSpectra\`, \`PubChemSafety\`, \`PubChemExperimental\`,
\`PubChemBioactivity\`, \`PubChemIdentifiers\`, \`SafetyProfile\`,
\`FEMAProfile\`, \`FDA_SPL_Profile\`, \`LOTUSProfile\`,
\`PubChemClassifications\`, \`MeSHProfile\`, \`LiteratureProfile\`,
\`ChemicalTerms\`, \`ChemicalTraits\`, \`ChemicalTraitOntology\`,
\`ChemicalTraitMatrix\`, \`ChemicalTraitOntologyMatrix\`,
\`ChemicalTraitEvidence\`, \`ChemicalTraitReport\`,
\`ChemicalTraitSummary\`, \`ChemicalTraitSimilarity\`,
\`ChemicalClasses\`, \`ChemicalMeasurements\`,
\`ChemicalMeasurementSummary\`, \`ChemicalHazards\`, \`ChemicalUses\`,
\`ChemicalBioassays\`, \`ChemicalBioactivities\`, \`ChemicalTargets\`,
\`ChemicalPotencies\`, \`ChemicalTaxonomy\`, \`ChemicalOccurrences\`,
\`ChemicalPathwayRoles\`, \`KEGGReactionParticipants\`, \`KEGGMatches\`,
\`KEGGSearchCandidates\`, \`KEGGRecords\`, \`KEGGIdentifiers\`,
\`KEGGPathways\`, \`KEGGReactions\`, \`KEGGEnzymes\`, \`KEGGModules\`,
\`KEGGLinks\`, \`KEGGLinkMetadata\`, \`KEGGClassifications\`,
\`SourceCoverage\`, \`DerivedGroups\`, \`Provenance\`,
\`DataDictionary\`, \`TableQuality\`, \`SourceDiagnostics\`,
\`ValidationIssues\`, and \`ValidationSummary\`. The normalized
\`Chemical\*\` tables split source evidence into discrete terms, traits,
classes, measurements, standardized measurement summaries, hazards,
uses, bioactivity calls, assay targets, potency values, taxonomy fields,
organism occurrences, pathway roles, and reaction participants for
grouping and comparative analysis. \`ChemicalTraits\` provides a
cross-source long-form trait vocabulary. \`ChemicalTraitOntology\` maps
whitelisted discrete traits to controlled ontology domains and groups
while preserving source evidence and source-backed identifiers where
available. \`ChemicalTraitMatrix\` and \`ChemicalTraitOntologyMatrix\`
provide compact binary wide matrices for filtering, clustering,
ordination, heatmaps, and model inputs. Use \`chemicalTraitMatrix()\`
and \`chemicalTraitOntologyMatrix()\` to rebuild full,
confidence-weighted, or capped matrices without rerunning web requests.
\`ChemicalTraitEvidence\` and \`chemicalTraitEvidence()\` trace matrix
and ontology keys back to source evidence, IDs, URLs, extraction rules,
and confidence values. \`ChemicalTraitReport\` and
\`chemicalTraitReport()\` provide one compact, researcher-facing row per
compound with domain signals, evidence coverage, identifiers, and
nearest ontology neighbors. \`ChemicalTraitSummary\` summarizes trait
breadth by compound and trait type, while \`ChemicalTraitSimilarity\`
compares compounds by shared and distinct traits. \`DataDictionary\`,
\`TableQuality\`, \`SourceDiagnostics\`, \`ValidationIssues\`, and
\`ValidationSummary\` document expected schemas and audit each result
for missing tables, missing required columns, type mismatches, duplicate
analysis keys, completeness, and source-specific coverage.

## Details

Provides a detailed overview of categorical and structural information
for every chemical of interest. Functional matches are also generated by
matching with groups from a chemical library.

## Examples

``` r
if (FALSE) { # \dontrun{
compounds = c("3-Octanone","Decane","Mesitylene","1,2,4-trimethyl-benzene",
"D-Limonene","beta-ethyl benzeneethanol","1,4-diethyl benzene",
"1,2-diethyl benzene","1,3,8-p-Menthatriene","(2-methyl-1-propenyl)-Benzene",
"1-Phenyl-1-butene","Linalool","Nonanal","5-nonyl-2-Thiophenecarboxylic acid",
"Dichloroacetaldehyde","Linalyl acetate","Beta-Ocimene")
categorate(compounds, library_data, input_format = "wide")
} # }
```
