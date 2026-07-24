# Enrich plant context evidence from source-backed literature text

\`enrichPlantContextEvidence()\` extracts plant-part, tissue, and method
context from source text linked to occurrence rows by PMID or DOI. The
function is intentionally conservative: source text must contain a
context term and must mention the occurrence species/genus, the
compound, or a chemical-profile source phrase before context is emitted.
Rows preserve the source PMID/DOI, source field, matched text,
extraction rule, confidence, and review flag.

When \`x\` is a plant phytochemistry result, the returned result has
updated \`PlantContextEvidence\`, \`PlantCompoundOccurrences\`, and
\`ProviderContextAudit\` tables. When \`x\` is an occurrence table and
\`apply = FALSE\`, the function returns only the source-backed
\`PlantContextEvidence\` rows.

## Usage

``` r
enrichPlantContextEvidence(
  x,
  context_sources = NULL,
  fetch_pubmed = FALSE,
  cache = TRUE,
  cache_dir = NULL,
  throttle = 0.34,
  ncbi_email = Sys.getenv("NCBI_EMAIL", ""),
  ncbi_tool = Sys.getenv("NCBI_TOOL", "uafR"),
  ncbi_api_key = Sys.getenv("NCBI_API_KEY", ""),
  max_sources = 100,
  request_fun = NULL,
  request_timeout = 30,
  min_confidence = "low",
  apply = inherits(x, "uaf_plant_phytochemistry") || (is.list(x) &&
    is.data.frame(x$PlantCompoundOccurrences))
)
```

## Arguments

- x:

  Plant phytochemistry result, normalized occurrence table, or named
  list containing \`PlantCompoundOccurrences\`.

- context_sources:

  Optional source text table. Columns are matched to the
  \`LiteratureCandidates\` schema; useful fields include \`pmid\`,
  \`doi\`, \`title\`, \`abstract\`, \`evidence_text\`, and
  \`evidence_url\`.

- fetch_pubmed:

  Logical. If \`TRUE\`, fetch PubMed abstracts for occurrence PMIDs and
  DOI-to-PMID matches. Network access is never attempted when this is
  \`FALSE\`.

- cache:

  Logical. If \`TRUE\`, cache PubMed requests.

- cache_dir:

  Cache directory for PubMed context requests.

- throttle:

  Seconds to wait between uncached PubMed requests.

- ncbi_email:

  Optional NCBI email.

- ncbi_tool:

  Optional NCBI tool name.

- ncbi_api_key:

  Optional NCBI API key. It is not written to outputs.

- max_sources:

  Maximum unique PubMed source records to fetch. Candidate PMID/DOI
  records are prioritized by available source text, direct occurrence
  evidence, and shared species/compound coverage before fetching.

- request_fun:

  Optional request function for tests or controlled HTTP.

- request_timeout:

  Maximum seconds allowed for uncached PubMed requests.

- min_confidence:

  Minimum context confidence to retain.

- apply:

  Logical. If \`TRUE\`, apply the source-backed context evidence to
  occurrence rows when returning a result/list.

## Value

A source-backed \`PlantContextEvidence\` table, or an updated plant
phytochemistry result/list when \`x\` is a result/list and \`apply =
TRUE\`.
