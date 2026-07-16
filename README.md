
# uafR: reproducible chemical enrichment and plant chemistry workflows

<!-- badges: start -->
<!-- badges: end -->

## Objective

uafR organizes tentative GC-MS annotations, public-database chemical
enrichment, reported plant-compound evidence, comparable chemistry,
structure-based similarity, and analysis handoff bundles. It preserves source
provenance and uncertainty; it does not convert database records or tentative
library hits into confirmed sample chemistry.

The current development line is an internal release candidate. Stable and
experimental interfaces are listed by `uafRApiStability()`. The package remains
private during hardening, and source ZIP/student-bundle installation is the
default sharing route.

## Production workflow guide

uafR now includes package-level helpers for choosing workflows, checking API
stability, planning large runs, and documenting scientific limits:

``` r
uafRWorkflowGuide()
uafRApiStability()
uafRProviderContracts()
uafRClaimGuidance()
```

Stable APIs are intended for downstream scripts and training materials.
Experimental APIs are usable but may still gain columns, diagnostics, or
provider-specific hardening as large plant and database workflows mature.
Project-specific scripts under `tools/` remain examples or wrappers and should
not be treated as general package APIs unless they are promoted into exported
functions.

Large species-first projects should be planned before live queries are run:

``` r
plan = planPlantChemistryRun(
  plants = "project_species.csv",
  compounds = "resolved_compounds.csv",
  sources = c("lotus", "pubmed", "pubtator"),
  cache_dir = "uafR_plant_cache",
  lotus_index = "lotus_cache/exports/LOTUS_lookup_index"
)

plan$Summary
plan$ProviderPlan
plan$OutputEstimates
plan$Recommendations
```

Use `inspectUafRCache()` and `summarizeUafRCache()` to audit local cache
coverage before repeating expensive PubChem, KEGG, PubMed, PubTator, or plant
provider workflows. These planning helpers do not query live web services.
Provider interpretation contracts are available from `uafRProviderContracts()`;
the package also ships `PROVIDER_SOURCES.md` with official source entry points
and third-party data-use boundaries. Users remain responsible for checking
current provider terms before redistributing downloaded source records.

### Which workflow should I use?

| Starting point | Recommended path | Main outputs |
| --- | --- | --- |
| GC-MS hit tables | `spreadOut()`, `mzExacto()`, `exactoThese()` | cleaned hit tables and exact-mass/library subsets |
| compound names or CIDs | `categorate(detail = "research")`, `pubchemProfile()`, `keggProfile()` | source-backed enrichment, validation, trait tables |
| plant species names | `resolvePlantPhytochemistry()` for discovery, then `runPlantChemistryProject()` for bundle handoff | plant-compound occurrences, evidence grades, matrices |
| curated plant-compound rows | `standardizePlantCompoundIntake()`, `runPlantChemistryProject()` | finalized analysis bundle without live discovery |
| large plant panel | `planPlantChemistryRun()`, cached/local LOTUS workflows, retry queues | run plan, cache summary, retry queue, manifest |
| Tanimoto similarity | `preparePlantTanimotoInput()`, `chemicalTanimotoSimilarity()`, `plantChemicalTanimotoSimilarity()`, `plantComparableTanimotoSummary()` | identity-audited server handoffs plus compound and plant-pair similarity summaries |
| model-ready matrices | `exportPlantChemistryFeatureSet()` or finalized bundle feature files | species x chemistry/source/evidence/context matrices |
| private student install | `tools/build_student_bundle.R` and the bundle scripts | source ZIP, install scripts, offline acceptance test |

Large-run recovery helpers are intentionally offline and manifest-based:

``` r
validation = validatePlantChemistryRunManifest("batch_manifest.csv")
retry_queue = writePlantChemistryRetryQueue(
  validation$Manifest,
  path = "retry_queue.csv",
  overwrite = TRUE
)

# No provider calls occur until dry_run is turned off and a runner is supplied.
rerunFailedPlantQueries(retry_queue, dry_run = TRUE)
```

Compound identity decisions should be reviewed through explicit audit tables,
not hidden name substitutions:

``` r
audit = standardizeCompoundIdentityAudit(phyto$CompoundResolution)
validateCompoundIdentityAudit(audit)$Summary
template = exportCompoundIdentityReviewTemplate(
  phyto$CompoundResolution,
  path = "compound_identity_review.csv",
  overwrite = TRUE
)
```

## Installation

### Student install from a private bundle

Student machines should install uafR from a private student bundle distributed
by the DSU dsDNA Core program. This keeps the repository private and gives
students a direct RStudio workflow.

The bundle contains a local uafR source archive, checksum verification,
preflight checks, installer and update scripts, verification scripts, an
offline acceptance test, a quick reference, and the training manual. The
installer still uses CRAN and Bioconductor for dependencies.

From the unzipped bundle folder, students should open `START_HERE.md`, then run
these scripts in RStudio:

``` r
source("verify_bundle_integrity.R")
source("preflight_check.R")
source("install_uafR_from_bundle.R")
source("run_student_acceptance_test.R")
```

After installation, confirm the package loads:

``` r
library(uafR)
data("library_data", package = "uafR")
```

### Build a student bundle

From the repository root:

``` sh
Rscript tools/build_student_bundle.R
```

The script creates:

``` text
student_bundle/uafR_student_bundle_<version>/
student_bundle/uafR_student_bundle_<version>.zip
```

The bundle folder includes:

``` text
START_HERE.md
README_STUDENT_INSTALL.md
uafR_QUICK_REFERENCE.md
CHECKSUMS.csv
verify_bundle_integrity.R
preflight_check.R
install_uafR_from_bundle.R
update_uafR_from_bundle.R
verify_uafR_install.R
run_student_acceptance_test.R
packages/uafR_<version>.tar.gz
training/uafR_training_manual.pdf
training/scripts/
training/data/
examples/test_install.R
```

### Developer install from a local checkout

Developers working from the private repository can install from the local source
checkout:

``` r
install.packages(c("remotes", "BiocManager"))
BiocManager::install(c("ChemmineR", "fmcsR"), ask = FALSE, update = FALSE)
remotes::install_local(
  ".",
  dependencies = c("Depends", "Imports"),
  upgrade = "never",
  build_vignettes = FALSE
)
```

### Optional live smoke test

After installation, a small live database smoke test can be run when internet
access is available:

``` r
library(uafR)
data("library_data", package = "uafR")

quick_result = categorate(
  compounds = c("aspirin", "caffeine"),
  chemical_library = library_data,
  input_format = "wide",
  detail = "research",
  cache = FALSE,
  throttle = 0.1,
  assay_detail_limit = 0
)

validateCategorateResult(quick_result)$Summary
quick_result$ChemicalTraitReport
quick_result$ChemicalMeasurementSummary
```

## Species-first plant phytochemistry workflows

uafR can now start from plant species names instead of a known compound list.
The plant phytochemistry resolver attempts to normalize plant names, collect
reported species-compound evidence from normalized provider adapters or curated
intake tables, include conservative literature/PubTator candidate evidence when
available, and then reuse uafR compound enrichment for resolved compounds.
When a `chemical_library` is supplied, enrichment can use the full
`categorate()` workflow, including library-based FMCS matching. When no
chemical library is supplied, uafR now falls back to `pubchemProfile()` and
still returns PubChem-derived properties, traits, matrices, validation, and
provenance.
Current live-capable public adapters include PubMed literature search,
PubTator candidate co-mentions, KNApSAcK organism-metabolite lookup, conservative
LOTUS API parsing when taxon evidence is present, and PubChem taxonomy
annotations after NCBI taxonomy resolution. NPASS is supported through a local,
manifest-backed index built from the official NPASS 3.0/NPASS-2026 general,
structure, species-source, and taxonomy downloads. Large projects should use
local LOTUS and NPASS indexes rather than repeatedly scraping provider pages.

For LOTUS specifically, use a local index for serious plant panels. The LOTUS
simple web API is useful for small smoke tests, but common species can return
very large unpaged payloads. uafR therefore supports `lotus_index`, a flat CSV,
TSV, JSON, JSONL, NDJSON, RDS, or data-frame index produced from LOTUS MongoDB,
Wikidata, SDF metadata, or another LOTUS export. `standardizeLotusIndex()`
normalizes flexible column names such as `allTaxa`, `traditional_name`,
`lotus_id`, `inchikey`, `smiles`, `molecular_formula`, `doi`, and `pmid`.
`queryLotusIndex()` returns the same `PlantCompoundOccurrences` schema used by
the full resolver. When `lotus_index` is supplied, `resolvePlantPhytochemistry()`
uses the local index instead of the live LOTUS API.

Build the compact index once, then reuse it across projects:

``` r
lotus_build = buildLotusIndex(
  input = "lotus_flat_export.jsonl",
  out_file = "lotus_compact_index.csv",
  overwrite = TRUE
)

lotus_build$BuildSummary
```

The same step can be run from Terminal:

``` sh
Rscript tools/build_lotus_index.R \
  --input lotus_flat_export.jsonl \
  --out-file lotus_compact_index.csv \
  --overwrite
```

The official LOTUS SMILES and SDF downloads are useful for structure work, but
they do not carry the source-backed species occurrence records needed for
species-first phytochemistry. For plant runs, use the official MongoDB ZIP
download and flatten the `lotusUniqueNaturalProduct.bson` collection into a
taxon-compound CSV first. The repository includes a Python standard-library
helper for this step, so MongoDB command-line tools are not required:

``` sh
mkdir -p lotus_cache/downloads lotus_cache/exports

curl -L \
  -o lotus_cache/downloads/LOTUSlatest.zip \
  https://lotus.naturalproducts.net/download/mongo

PYTHONDONTWRITEBYTECODE=1 python3 tools/flatten_lotus_mongo_dump.py \
  --input lotus_cache/downloads/LOTUSlatest.zip \
  --out-file lotus_cache/exports/LOTUS_mongo_flat.csv \
  --compact-index-file lotus_cache/exports/LOTUS_compact_index.csv \
  --lookup-dir lotus_cache/exports/LOTUS_lookup_index \
  --overwrite
```

`LOTUS_mongo_flat.csv` is the auditable one-row-per-compound/taxon/reference
export. `LOTUS_compact_index.csv` is the smaller, analysis-ready index to pass
to uafR for small panels or ad hoc inspection. `LOTUS_lookup_index/` is the
production path for larger plant panels because uafR can read only the shards
needed for the submitted species, genus, and family keys instead of scanning
the complete CSV. Build products under `lotus_cache/` are ignored by package
builds and should be treated as local data cache files, not source files.

If the compact CSV already exists, build only the lookup directory without
rewriting the flat export:

``` sh
PYTHONDONTWRITEBYTECODE=1 python3 tools/flatten_lotus_mongo_dump.py \
  --from-compact-index lotus_cache/exports/LOTUS_compact_index.csv \
  --lookup-dir lotus_cache/exports/LOTUS_lookup_index \
  --overwrite
```

This workflow is designed for plant, ecology, remediation, metabolomics,
chemical ecology, natural-products, and environmental chemistry projects. It
does not produce complete metabolomes. Absence of public records is not absence
of compounds. Literature co-mentions are candidate evidence unless curated.
PubChem taxonomy annotations are source-backed chemical associations, not proof
that a student's sample contains the compound. Direct species evidence is
stronger than genus or family fallback evidence. Always inspect validation and
provenance before interpretation.

Evidence filters distinguish exact source records from aggregate or inferred
associations. Exact species records from a validated local LOTUS index, the
official NPASS index, or an exact-taxon KNApSAcK result can qualify as direct
database evidence. PubChem taxonomy collection associations remain
review-required unless a row-specific source explicitly links the species and
compound. Broad collection citation lists are retained as source diagnostics;
they are not converted into PMIDs, biological context, or direct occurrence
claims. Genus/family fallbacks and PubMed/PubTator candidates remain separate
from direct evidence in every analysis-ready export.

``` r
library(uafR)
data("library_data", package = "uafR")

plants = c("Salix nigra", "Camellia sinensis", "Zea mays")

phyto = resolvePlantPhytochemistry(
  plants = plants,
  sources = c("lotus", "pubmed", "pubtator"),
  taxon_fallback = c("species", "genus"),
  lotus_index = "lotus_cache/exports/LOTUS_lookup_index",
  enrich_compounds = TRUE,
  chemical_library = library_data,
  detail = "research",
  cache = TRUE,
  cache_dir = "uafR_plant_cache",
  max_pubmed_records = 25,
  max_provider_records = 100,
  enrichment_batch_size = 25,
  resume_enrichment = TRUE
)

lotus_occurrences = queryLotusIndex(
  plants,
  lotus_index = "lotus_cache/exports/LOTUS_lookup_index",
  taxon_fallback = c("species", "genus")
)

phyto$SpeciesChemistrySummary
phyto$PlantCompoundOccurrences
phyto$PlantContextEvidence
validatePlantPhytochemistryResult(phyto)$Summary

# Direct species database records and curated rows are analysis-ready by
# default. Literature co-mentions and genus/family fallbacks remain available
# for review but are not treated as analysis-ready occurrence evidence.
analysis_ready = filterPlantPhytochemistryEvidence(phyto)

# Restrict analyses to compounds reported from a known biological context when
# that metadata is available. Supported groups include root_belowground, leaf,
# stem_shoot, bark_wood, flower, fruit_seed, aerial, whole_plant, and
# exudate_rhizosphere.
leaf_records = filterPlantPhytochemistryEvidence(
  phyto,
  plant_part_group = "leaf"
)

# Candidate, fallback, literature, and unresolved rows should be reviewed
# before they are promoted into analysis-ready occurrence evidence.
review = plantPhytochemistryReviewTable(phyto)
if (nrow(review) > 0) {
  review$review_decision[1] = "promote_curated"
  review$reviewed_by[1] = "researcher name"
  review$review_note[1] = "Source reports the compound from leaf tissue."
  review$proposed_citation_or_url[1] = "https://doi.org/example"
  review$proposed_plant_part[1] = "leaf"
  review$proposed_method[1] = "LC-MS"
  phyto_reviewed = applyPlantPhytochemistryReview(phyto, review)
}

exportPlantPhytochemistryWorkbook(
  phyto,
  path = "plant_phytochemistry_export",
  format = "csv",
  overwrite = TRUE
)

exportPlantPhytochemistryWorkbook(
  phyto,
  path = "plant_phytochemistry_analysis_ready_export",
  format = "csv",
  preset = "analysis_ready",
  overwrite = TRUE
)
```

Use the full export for audit/provenance and the analysis-ready preset for
downstream matrices or figures where candidate co-mentions and taxon fallbacks
should be excluded.

Use structural Tanimoto similarity when the question is molecular resemblance,
not shared database traits. `chemicalTraitSimilarity()` compares normalized
traits from categorate enrichment. `chemicalTanimotoSimilarity()` and
`plantChemicalTanimotoSimilarity()` compare PubChem Fingerprint2D bit vectors
and report true compound-pair structural similarity. The plant helper keeps
species names, compound IDs, evidence tiers, and plant context attached so the
result can be joined directly to phylogenetic, ecological, remediation, or
sample-similarity analyses.

``` r
# Pairwise molecular similarity for a plant phytochemistry result. For small
# projects, the pairwise tables can be returned in memory.
plant_tanimoto = plantChemicalTanimotoSimilarity(
  phyto,
  cache = TRUE,
  cache_dir = "uafR_plant_cache/tanimoto",
  return_group_compound_pairs = FALSE
)

plant_tanimoto$PlantPairTanimotoSummary[, c(
  "species_a", "species_b", "compound_pair_count",
  "shared_compound_count", "mean_tanimoto", "median_tanimoto",
  "p95_tanimoto", "max_tanimoto",
  "compound_pair_count_ge_0_85"
)]

# For larger plant panels, first create an offline, identity-audited handoff.
# This preserves all evidence rows but emits only one membership per species
# and source-backed structure. It does not contact PubChem or calculate pairs.
prepared = preparePlantTanimotoInput(
  phyto,
  occurrence_status = c("direct_reported", "curated_reported"),
  analysis_ready = TRUE,
  min_confidence = "medium",
  include_review_required = FALSE,
  out_dir = "plant_tanimoto_handoff",
  overwrite = FALSE,
  strict = TRUE
)

prepared$Summary
prepared$ValidationSummary
prepared$DuplicateAudit
prepared$NameStructureAudit
prepared$ExcludedRows
```

`preparePlantTanimotoInput()` treats repeated publications and source records
as evidence, not extra compounds. It collapses only the species-by-structure
membership, prefers exact source-record InChIKeys, removes any inherited
fingerprint columns, and leaves CID blank whenever an InChIKey is available so
the server must resolve and verify that InChIKey against PubChem. The generated
`NameStructureAudit` preserves cases where one normalized compound label maps
to multiple exact structures; these are not treated as duplicates or silently
merged. The generated CSV and JSON manifests contain row counts and checksums.
Copy the complete handoff directory and the validated uafR source version to
the server. Production panel handoffs also contain
`plant_species_universe.csv`; pass it to the server so plants with fewer than
two usable structures remain present as explicit `insufficient_support` pairs
instead of disappearing from the result.

Run the server gates in order. Explicit `preflight` mode performs no network
requests:

``` sh
Rscript tools/run_plant_tanimoto_server.R \
  --input plant_tanimoto_handoff/plant_compound_membership_tanimoto_ready.csv \
  --species-universe plant_tanimoto_handoff/plant_species_universe.csv \
  --manifest plant_tanimoto_handoff/tanimoto_input_export_manifest.csv \
  --release-manifest uafR_release_manifest.json \
  --out-dir plant_tanimoto_server_output \
  --cache-dir pubchem_tanimoto_cache \
  --mode preflight --full-pairs true \
  --throttle 1.1
```

The preflight validates checksums, identity formats, unique
species-by-structure keys, review exclusions, installed package version,
dependencies, write access, disk space, and estimated pair counts. Continue
with the same handoff and cache:

``` sh
UAFR_CONFIRM_SERVER_TANIMOTO=YES Rscript tools/run_plant_tanimoto_server.R \
  --input plant_tanimoto_handoff/plant_compound_membership_tanimoto_ready.csv \
  --species-universe plant_tanimoto_handoff/plant_species_universe.csv \
  --manifest plant_tanimoto_handoff/tanimoto_input_export_manifest.csv \
  --release-manifest uafR_release_manifest.json \
  --out-dir plant_tanimoto_server_output \
  --cache-dir pubchem_tanimoto_cache \
  --mode smoke --smoke-structures 25 --throttle 1.1

UAFR_CONFIRM_SERVER_TANIMOTO=YES Rscript tools/run_plant_tanimoto_server.R \
  --input plant_tanimoto_handoff/plant_compound_membership_tanimoto_ready.csv \
  --species-universe plant_tanimoto_handoff/plant_species_universe.csv \
  --manifest plant_tanimoto_handoff/tanimoto_input_export_manifest.csv \
  --release-manifest uafR_release_manifest.json \
  --out-dir plant_tanimoto_server_output \
  --cache-dir pubchem_tanimoto_cache \
  --mode summary --throttle 1.1

UAFR_CONFIRM_SERVER_TANIMOTO=YES Rscript tools/run_plant_tanimoto_server.R \
  --input plant_tanimoto_handoff/plant_compound_membership_tanimoto_ready.csv \
  --species-universe plant_tanimoto_handoff/plant_species_universe.csv \
  --manifest plant_tanimoto_handoff/tanimoto_input_export_manifest.csv \
  --release-manifest uafR_release_manifest.json \
  --out-dir plant_tanimoto_server_output \
  --cache-dir pubchem_tanimoto_cache \
  --mode full --full-pairs true --throttle 1.1
```

The smoke subset is deterministic and selected by round-robin plant coverage.
The summary mode resolves all fingerprints and writes compact plant-pair
summaries without full row-level pair tables. Full mode reuses that cache and
publishes requested pair files only after row-count and identity validation.
Runs pause with exit status 75 after repeated PubChem `429`/`503` responses;
successful cache entries remain reusable. Inspect `server_run_status.json`,
`server_progress.csv`, and the completion markers before advancing. Outputs
larger than 25 million rows are written in compressed shards; outputs above 250
million rows require a completed summary gate and explicit extreme-output
confirmation. Compound-name fallback is disabled throughout. The server also
writes exact comparable-scope and comparable-group summaries from the same
decoded fingerprints. Unknown and non-comparable classifications are retained
for audit but excluded from these comparable summaries by default.

For a direct large run from a trusted in-memory result, stream the large
pairwise tables to compressed CSV files and keep the compact species-pair
summary in the returned object:

``` r
plant_tanimoto = plantChemicalTanimotoSimilarity(
  phyto,
  cache = TRUE,
  cache_dir = "uafR_plant_cache/tanimoto",
  out_dir = "plant_tanimoto_export",
  return_compound_pairs = TRUE,
  return_group_compound_pairs = TRUE
)

plant_tanimoto$ExportManifest
plant_tanimoto$PlantPairTanimotoSummary
```

Use the comparability layer before clustering, ordination, scoring, or group
comparisons. `metabolite` is a broad biological term, while `volatile`
describes an analytical or physicochemical fraction; those should not be
compared as if they were the same axis. uafR therefore records biological
domain, biosynthetic family, analytical behavior, comparison scope, confidence,
source-backed classification fields, and caveats in `ChemistryComparability`.
The classifier prioritizes normalized source fields from `ChemicalClasses`,
`LOTUSProfile`, `PubChemClassifications`, `ChemicalTraitOntology`,
`ChemicalTerms`, `KEGGClassifications`, `KEGGPathways`, and `DerivedGroups`
before it falls back to compound-name patterns.
Biological-context evidence is tracked separately in `PlantContextEvidence` so
plant-part, tissue, and method filters can be audited instead of being hidden
inside a matrix. Provider-level coverage and review burden are summarized in
`ProviderContextAudit`, which helps identify whether context gaps come from a
specific source, candidate-only literature evidence, or sparse provider records.

``` r
# Inspect how each plant-compound row was assigned to a comparison scope.
phyto$ChemistryComparability[, c(
  "species", "compound_name", "metabolism_domain",
  "biosynthetic_family", "chemical_behavior",
  "comparison_scope", "comparison_group",
  "comparability_confidence", "comparability_basis",
  "classification_source_table", "classification_source_value"
)]

# Inspect provider-level biological-context coverage before trusting a
# context-specific matrix.
phyto$ProviderContextAudit[, c(
  "source_database", "occurrence_count", "context_known_fraction",
  "review_required_context_fraction", "audit_status",
  "recommended_action"
)]

# Add source-backed context from local literature/source text. Rows must match
# occurrence evidence by PMID or DOI; context is not filled unless the source
# text contains plant-part, tissue, or method terms linked to the species/genus
# or compound in a chemical-source context.
source_text = data.frame(
  pmid = "12345678",
  doi = "10.1000/example",
  title = "Phytochemical constituents from Salix nigra leaves",
  abstract = paste(
    "Salicin was identified from leaves of Salix nigra.",
    "The extract was analyzed by LC-MS."
  )
)

phyto = enrichPlantContextEvidence(
  phyto,
  context_sources = source_text
)

# Optional live PubMed context fetches are explicit, cached, capped, and
# prioritized toward source records with useful context text or broad
# species/compound coverage.
phyto = enrichPlantContextEvidence(
  phyto,
  fetch_pubmed = TRUE,
  cache = TRUE,
  cache_dir = "uafR_plant_cache/context",
  max_sources = 50
)

# Compare specialized plant chemistry separately from primary metabolism.
specialized_matrix = plantComparableChemistryMatrix(
  phyto,
  comparison_scope = "specialized_metabolites",
  feature = "comparison_group",
  mode = "binary"
)

# Compare volatile specialized chemistry only within the volatile fraction.
volatile_matrix = plantComparableChemistryMatrix(
  phyto,
  comparison_scope = "volatile_specialized_metabolites",
  feature = "comparison_group",
  mode = "binary"
)

# Compare leaf volatile chemistry only when plant-part or method context is
# known. This prevents root exudate, whole-plant, and unspecified-extract rows
# from being mixed into a leaf-focused matrix.
leaf_volatile_matrix = plantComparableChemistryMatrix(
  phyto,
  comparison_scope = "volatile_specialized_metabolites",
  feature = "comparison_group",
  plant_part_group = "leaf",
  require_context = TRUE,
  mode = "binary"
)

# Primary metabolites are a different comparison scope.
primary_matrix = plantComparableChemistryMatrix(
  phyto,
  comparison_scope = "primary_metabolites",
  feature = "comparison_group",
  mode = "binary"
)
```

Unknown or broad rows stay in the audit table but are excluded from comparable
matrices by default. This is intentional: a missing class is not evidence that a
compound belongs in a broad mixed analysis.

Before scaling to hundreds of plants, run a small provider pilot and inspect
the output. The pilot workflow writes the same batch audit tables plus a compact
species summary, a QA report, a review-needed table, and context-aware
comparable matrices for specialized, volatile-specialized, leaf-associated,
root/exudate-associated, and primary-metabolism chemistry.

``` r
plants = plantPhytochemistryPilotPanel()$species

pilot = runPlantPhytochemistryPilot(
  plants = plants,
  sources = c("lotus", "knapsack", "npass", "pubchem", "pubmed", "pubtator"),
  out_dir = "plant_phytochemistry_pilot",
  cache_dir = "uafR_plant_cache",
  species_chunk_size = 10,
  compound_resolution_profile = "identity",
  max_pubmed_records = 25,
  max_provider_records = 100,
  request_timeout = 30,
  overwrite = TRUE
)

pilot$PilotSummary
pilot$PilotQAReport
pilot$ProviderContextAudit
pilot$PilotMatrices$matrix_volatile_specialized_metabolites
pilot$PilotMatrices$matrix_root_exudate_associated_chemistry
```

The command-line wrapper accepts a CSV plant list:

``` sh
Rscript tools/run_plant_phytochemistry_pilot.R \
  --plant-csv plants.csv \
  --species-col species \
  --out-dir plant_phytochemistry_pilot \
  --cache-dir uafR_plant_cache \
  --lotus-index lotus_compact_index.csv \
  --sources lotus,knapsack,npass,pubchem,pubmed,pubtator \
  --compound-resolution-profile identity \
  --max-pubmed-records 25 \
  --max-provider-records 100 \
  --request-timeout 30 \
  --overwrite
```

The wrapper can also use the built-in 15-species panel:

``` sh
Rscript tools/run_plant_phytochemistry_pilot.R \
  --default-panel \
  --out-dir plant_phytochemistry_pilot \
  --compound-resolution-profile identity \
  --overwrite
```

Start with `compound_resolution_profile = "identity"` for a provider/data-depth
pilot. Move to `"research"` only after the occurrence evidence, context
coverage, provider context audit, and review burden look reasonable.

For a large local-LOTUS panel, use a stratified pilot before touching the full
run directory. The representative wrapper selects deterministic, alphabetically
spread species from three precomputed lookup strata: exact species key, genus
key without an exact species key, and no species/genus key. The default 25-row
panel uses 10, 10, and 5 species, respectively. It performs species-level local
LOTUS discovery only and leaves compound resolution disabled:

``` sh
Rscript tools/run_representative_lotus_pilot.R \
  --plant-csv plant_species_run_input.csv \
  --lotus-index LOTUS_lookup_index \
  --out-dir representative_lotus_pilot \
  --cache-dir representative_lotus_pilot_cache \
  --exact-count 10 \
  --genus-only-count 10 \
  --no-key-count 5 \
  --overwrite true
```

Inspect `representative_pilot_query_accounting.csv`,
`representative_pilot_validation.csv`, `provider_diagnostics.csv`, and
`representative_pilot_manifest.csv`. Exact-key species should exercise source
record extraction; genus-only and no-key species test explicit no-hit behavior
when genus fallback is disabled. The pilot never interprets a no-hit as
biological absence and never launches PubChem or literature providers.

For quick species-first discovery without a library, omit `chemical_library`.
This produces PubChem-only enrichment and skips FMCS library matching:

``` r
phyto = resolvePlantPhytochemistry(
  plants = c("Camellia sinensis", "Salix nigra"),
  sources = c("lotus", "knapsack", "pubmed"),
  enrich_compounds = TRUE,
  detail = "research",
  cache = TRUE,
  cache_dir = "uafR_plant_cache",
  max_provider_records = 50
)

names(phyto$SpeciesChemistryMatrix)
phyto$CompoundResolution
```

For 100 or more plants, use a staged run. Do not submit the complete panel to
every live provider in one call. A large run requires `cache = TRUE`, a
persistent `out_dir`, and local or explicitly authorized live discovery.
`planPlantChemistryRun()` performs an offline preflight and reports discovery
chunks, provider request lower bounds, identity batches, pairwise output size,
and blocking readiness checks.

``` r
plants = read.csv("plant_species.csv", stringsAsFactors = FALSE)$species
lotus_index = "lotus_cache/exports/LOTUS_lookup_index"
run_dir = "plant_phytochemistry_700"
cache_dir = "uafR_plant_cache"

plan = planPlantChemistryRun(
  plants = plants,
  sources = c("lotus", "pubmed"),
  cache_dir = cache_dir,
  lotus_index = lotus_index,
  species_chunk_size = 25,
  compound_batch_size = 25,
  max_pubmed_records = 10
)

plan$Summary
plan$InputNameAudit[plan$InputNameAudit$review_required |
                      plan$InputNameAudit$query_status == "blank_input", ]
plan$ReadinessChecks
plan$OutputEstimates
plan$RunConfiguration
```

For a multi-provider production panel, use `runPlantChemistryPanel()` or its
installed CLI instead of assembling the stages by hand. The runner freezes the
input and review ledgers, verifies local LOTUS and NPASS resources, executes a
small pilot, performs resumable discovery and compound enrichment, writes the
Tanimoto server handoff, and reconciles the returned server products into a
validated analysis bundle. Public-service calls are cached. Re-running the
same command reuses only checkpoints whose input, parameter, provider-resource,
and artifact signatures still match.

An audited production panel should run from a checked source artifact, not from
`devtools::load_all()`. After creating the approved clean Git commit, generate
the release manifest and tarball together:

``` sh
Rscript tools/check_package_release.R \
  --out-dir /private/tmp/uafR_production_release \
  --release-manifest /private/tmp/uafR_production_release/uafR_release_manifest.json
```

Install that tarball into a dedicated project library and run the installed
CLI. Supply the same manifest and tarball to every mode, preferably through one
JSON configuration file, with `require_release_artifact` set to `true`.
Preflight validates the clean Git commit, package version, tarball filename,
byte count, MD5, and SHA-256. These artifacts become part of the run signature
and are copied into the Tanimoto server handoff.

``` sh
PRODUCTION_LIB="project_software/R_library"
mkdir -p "$PRODUCTION_LIB"
R CMD INSTALL --library="$PRODUCTION_LIB" \
  /private/tmp/uafR_production_release/build/uafR_0.4.0.9000.tar.gz

RUNNER=$(R_LIBS_USER="$PRODUCTION_LIB" Rscript -e \
  'cat(system.file("scripts", "run_plant_chemistry_panel.R", package = "uafR"))')

R_LIBS_USER="$PRODUCTION_LIB" Rscript "$RUNNER" \
  --config plant_chemistry_run_config.json \
  --mode preflight
```

The repository `tools/run_plant_chemistry_panel.R` wrapper is for development;
it deliberately loads the checkout. Do not use it for a commit-bound
production data run.

Set NCBI contact information in the environment; never put an API key in a
script, JSON config, or command line:

``` sh
export NCBI_EMAIL="researcher@example.edu"
export NCBI_TOOL="uafR"
# Optional: export NCBI_API_KEY in the local shell only.
```

Run the gates separately while validating a new project. The abbreviated
commands below show the required structure; keep the same paths and options on
every resume:

``` sh
Rscript tools/run_plant_chemistry_panel.R \
  --mode preflight \
  --plant-csv plant_species_run_input.csv \
  --out-dir plant_phytochemistry_full_panel \
  --lotus-index LOTUS_lookup_index \
  --npass-raw-dir resources/NPASS_3.0_2026/raw \
  --sources lotus,npass,knapsack,pubchem,pubmed,pubtator \
  --expected-species-count 701 \
  --resume true

Rscript tools/run_plant_chemistry_panel.R \
  --mode build-indexes \
  --plant-csv plant_species_run_input.csv \
  --out-dir plant_phytochemistry_full_panel \
  --lotus-index LOTUS_lookup_index \
  --npass-raw-dir resources/NPASS_3.0_2026/raw \
  --download-npass true \
  --sources lotus,npass,knapsack,pubchem,pubmed,pubtator \
  --resume true

Rscript tools/run_plant_chemistry_panel.R \
  --mode pilot \
  --plant-csv plant_species_run_input.csv \
  --out-dir plant_phytochemistry_full_panel \
  --sources lotus,npass,knapsack,pubchem,pubmed,pubtator \
  --pilot-count 12 \
  --pilot-enrichment-compounds 25 \
  --resume true
```

The pilot bounds each provider to at most
`max(25, pilot_enrichment_compounds)` records per species (or an explicitly
smaller `max_provider_records` value). Discovery does not perform unbounded
name-based compound resolution. Its enrichment gate deterministically selects
only records carrying a source-backed CID or full InChIKey, and the second
discovery pass blocks live requests so `cache_only_rerun_identical` is a real
checkpoint/cache test.

Production discovery checkpoints use schema `3.3.0`. Each species chunk stores
only normalized provider tables; context linking, evidence summaries,
comparability, matrices, identity review, and validation are generated once
after the chunks are combined. A checkpoint created by an older schema is
rebuilt automatically. This makes resume behavior deterministic while avoiding
repeated derivation work in every five- or twenty-five-species chunk.

Continue through `discovery`, `identity`, `research-enrichment`,
`full-enrichment`, and `tanimoto-handoff` only after the preceding gate has
passed. `pipeline_status.json`, `pipeline_progress.csv`, stage manifests,
`failed_queries.csv`, and `retry_queue.csv` are the operational source of
truth. Exit status 75 means a provider returned repeated service-busy responses:
leave caches and completed checkpoints in place, wait for the service to
recover, and rerun the identical stage. Exit status 2 at finalization means the
validated Tanimoto server output has not yet been supplied. It is not a request
to restart discovery.

The NPASS downloader verifies remote byte counts, keeps interrupted transfers
as `.partial` files, and resumes them when the server supports byte ranges.
`buildNpassIndex()` streams the species-source table and writes lossless RDS
shards plus source and shard checksums. Raw provider downloads and KNApSAcK
responses remain local project resources; they are not package data and should
not be copied into redistributable result bundles without reviewing current
provider terms.

Resolve every `fail` readiness check before starting. In particular, genus-only
and `sp.`/`spp.` entries cannot support direct species-level interpretation;
replace them with accepted species names or preserve them in a separately
labeled lower-rank analysis. Stage 1 should use only the local LOTUS index and
should not contact PubChem. This creates versioned, atomic discovery
checkpoints plus an incremental manifest and retry queue.

``` r
phyto_discovery = runPlantPhytochemistryBatch(
  plants = plants,
  sources = "lotus",
  lotus_index = lotus_index,
  out_dir = run_dir,
  cache_dir = cache_dir,
  species_chunk_size = 25,
  compound_resolution_profile = "none",
  cache = TRUE,
  resume = TRUE,
  progress = TRUE,
  overwrite = TRUE
)

manifest_check = validatePlantChemistryRunManifest(
  phyto_discovery$BatchChunkManifest,
  base_dir = run_dir
)
manifest_check$Summary
manifest_check$RetryQueue
```

If R is interrupted, run the same call again. Completed checkpoints are reused;
failed, stopped, corrupt, or not-started chunks are retried. Compound resolution
is deferred automatically until every discovery chunk is complete. Once the
occurrence and evidence-review tables have been inspected, first recover exact
source-record structures and prepare the identity-safe server handoff. PubChem
fingerprints and pairwise similarity then follow the server preflight, smoke,
summary, and full gates described above. Rich categorate enrichment should use
the verified canonical CIDs from the successful server identity map rather than
restart ambiguous name resolution.

``` r
phyto_batch = runPlantPhytochemistryBatch(
  plants = plants,
  sources = "lotus",
  lotus_index = lotus_index,
  out_dir = run_dir,
  cache_dir = cache_dir,
  species_chunk_size = 25,
  compound_resolution_profile = "identity",
  compound_batch_size = 25,
  throttle = 0.5,
  cache = TRUE,
  resume = TRUE,
  progress = TRUE,
  overwrite = TRUE
)

phyto_batch$BatchRunManifest
phyto_batch$CompoundResolution
phyto_batch$CompoundIdentityReview
```

Run PubMed, PubTator, KNApSAcK, or PubChem occurrence searches only after a
small all-provider pilot passes and an identical cache-only replay is
reproducible. For a live follow-up of 100 or more species, use the production
panel runner or set `allow_large_live_run = TRUE` explicitly in lower-level
batch calls. NPASS discovery uses the documented local index; it is not queried
through an inferred live species endpoint.

The identity pass preserves chemically meaningful alpha/beta/gamma and
plus/minus prefixes in compound keys, and reports deterministic PubChem alias
matches with `MatchStatus = "resolved_alias"` in the identity audit table.
When `lotus_index` is supplied, the identity pass first recovers source-backed
LOTUS SMILES, InChIKeys, formulas, and CIDs by LOTUS record ID. PubChem name
lookup is then used only for remaining gaps. A separate
`CompoundIdentityReview` table flags source-ambiguous structures and rows where
the source provides a structure but the displayed label looks like a class,
mixture, plant product, or other non-discrete compound name.

Completed identity review worksheets can be reapplied to the result object so
manual structure decisions are reproducible:

``` r
identity_review = phyto_batch$CompoundIdentityReview

# Example: only fill these fields after checking the LOTUS source record,
# PubChem record, DOI, or another source-backed structure record.
identity_review$review_decision[1] = "update_identity"
identity_review$reviewed_by[1] = "researcher name"
identity_review$review_note[1] = "Source record supports this structure."
identity_review$proposed_compound_name[1] = "reviewed compound name"
identity_review$proposed_smiles[1] = "source-backed SMILES"
identity_review$proposed_inchikey[1] = "source-backed InChIKey"
identity_review$proposed_molecular_formula[1] = "source-backed formula"
identity_review$proposed_resolution_source[1] = "manual_identity_review"

phyto_batch_reviewed = applyPlantCompoundIdentityReview(
  phyto_batch,
  identity_review,
  reviewer = "researcher name"
)
```

The batch output directory includes `all_occurrences.csv`,
`analysis_ready_occurrences.csv`, `review_required_occurrences.csv`,
`compound_identity_resolution.csv`, `compound_identity_review.csv`,
`species_chemistry_summary.csv`, `species_chemistry_matrix.csv`,
`chemistry_comparability.csv`,
`comparable_chemistry_matrix.csv`, `provider_diagnostics.csv`,
`context_coverage_report.csv`, `batch_chunk_manifest.csv`,
`discovery_chunk_manifest.csv`, `failed_queries.csv`, `retry_queue.csv`, and
`run_manifest.json`. The offline test suite exercises interrupted/resumed
discovery across 705 synthetic species; this verifies orchestration and does
not imply complete public chemistry coverage for any real plant panel.

For projects that already have local curation, use the curated intake fallback:

``` r
curated = data.frame(
  species = "Salix nigra",
  compound_name = "salicin",
  source_database = "manual",
  citation_or_url = "https://example.org/source",
  evidence_tier = "manual_curated",
  plant_part = "bark",
  method = "LC-MS"
)

occurrences = standardizePlantCompoundIntake(curated)
phyto = resolvePlantPhytochemistry(
  plants = unique(curated$species),
  sources = character(),
  curated_data = occurrences,
  enrich_compounds = TRUE,
  chemical_library = library_data,
  detail = "research"
)
```

An offline training example is available at
`training/scripts/10_species_phytochemistry.R`. It uses simulated curated
species-compound rows and does not call live web services unless the
`UAFR_LIVE_PLANT_DISCOVERY` environment variable is set to `"true"`.

## Example Mass Spectrometry Workflows

These are basic examples of how to use core functions. The input .CSV file has strict column name/input data requirements. The column names MUST include: 'Component.RT', 'Component.Area', 'Base.Peak.MZ', 'File.Name', 'Compound.Name', and 'Match.Factor' in no particular order.

``` r
library(uafR)

input_dat = read.csv("your/gcms/dataset.csv")
```
 Component.RT  |  Base.Peak.MZ    |  Component.Area  |       Compound.Name        |  Match.Factor  |  File.Name  
:-------------:|:----------------:|:----------------:|:---------------------------|:--------------:|:------------:
8.229034       |84.00             |906.4701          |Pipradrol                   |62.62271        |Std_soln_07    
8.286703       |120.00            |209705.1878       |Methyl salicylate           |98.16152        |Std_soln_00a    
8.296408       |119.99            |30332.9022        |Methyl salicylate           |95.79911        |Std_soln_00    
8.303958       |120.00            |6476.4785         |Methyl salicylate           |86.29569        |Std_soln_07    
8.348031       |105.00            |420.8119          |3-Hexen-1-ol, benzoate, (Z)-|68.78156        |Std_soln_00    
**...**        |**...**           |**...**           |**...**                     |**...**         |**...**         

### In this example, the user knows what chemicals they are interested in:
``` r
input_spread = spreadOut(input_dat)
query_chemicals = c("Linalool", "Methyl Salicylate", "Limonene", "alpha-Thujene")

### extract the query_chemicals from the "spread out" input:
input_exacto = mzExacto(input_spread, query_chemicals)
```
### In this example, the user just wants to keep the top hits:
``` r
query_chemicals = input_dat$Compound.Name[input_dat$Match.Factor > 80]

input_exacto = mzExacto(input_spread, query_chemicals)
```

## Example Cheminformatics Workflow
``` r
## example usage for chemical informatics:
query_chemicals = c("Linalool", "Methyl Salicylate", "Limonene", "alpha-Thujene")
GroupA = c("Guaiacol",	"Tridecane",	"Ethyl heptanoate", "Caffeine")
GroupB = c("2-Aminothiazole", "Aspirin", "Octanoic acid", "alpha-Pinene", "Toluene")
chem_library = data.frame(cbind(GroupA, GroupB))

query_categorated = categorate(query_chemicals, chem_library, input_format = "wide")
```

## Detailed PubChem Enrichment

`pubchemProfile()` pulls richer PubChem data into reusable tables before downstream filtering or mass spectrometry extraction. This keeps web enrichment separate from `spreadOut()` and `mzExacto()`, so results can be cached, inspected, and tested independently.

``` r
query_chemicals = c("Methyl salicylate", "Octanal", "Undecane")

chem_profile = pubchemProfile(query_chemicals, profile = "ms")

chem_profile$identity
chem_profile$properties
chem_profile$synonyms
chem_profile$spectra
```

Available profiles are:

- `"minimal"`: PubChem CID, names, formula, identifiers, exact mass, and molecular weight.
- `"ms"`: minimal data plus chemical descriptors and mass spectrometry annotations.
- `"safety"`: minimal/descriptive data plus safety, hazard, and experimental property annotations.
- `"bioactivity"`: minimal/descriptive data plus PubChem assay summary data.
- `"full"`: all supported profile sections.

## aiNsect Molecular-Olfaction Export

`exportAiNsectMolOlfInputs()` converts paired wide essential-oil GC-MS profile
tables into aiNsect molecular-olfaction input files. The identity table should
contain treatment columns with compound names by ranked row. The abundance table
should contain the same treatments and rows with abundance values. Compound
structures are resolved through `pubchemProfile()`; uafR does not fabricate
SMILES, InChIKeys, CIDs, formulas, or abundance values.

``` sh
Rscript tools/export_ainsect_mololf_inputs.R \
  --chem-id-csv project_data/gcms_identity.csv \
  --chem-quant-csv project_data/gcms_abundance.csv \
  --out-dir project_results/mololf_export \
  --cache-dir project_cache/pubchem \
  --profile ms \
  --throttle 0.2
```

The export directory contains:

- `uafR_compounds.csv`: compound IDs, names, SMILES, InChIKeys, PubChem CIDs,
  molecular formulas, source labels, and notes.
- `treatment_compound_abundance.csv`: treatment/compound abundance rows with
  raw parsed abundance and within-treatment relative abundance.
- `uafR_compounds_unresolved.csv`: compounds that could not be exported because
  PubChem did not provide a usable SMILES.
- `uafR_mololf_export_summary.json`: run metadata, counts, output paths, and
  quality checks.
- `uafR_pubchem_identity_audit.csv` and `uafR_pubchem_properties_audit.csv`:
  optional PubChem audit tables.

## Research-Grade Database Enrichment

`categorate()` keeps the original eight-table output by default. Use
`detail = "research"` or `detail = "full"` when you want tidy PubChem and KEGG
tables appended to the categorated result. These tables are designed for
filtering, grouping, and downstream analyses rather than one-off text lookup.

``` r
query_categorated = categorate(
  query_chemicals,
  chem_library,
  input_format = "wide",
  detail = "research"
)

query_categorated$PubChemProperties
query_categorated$SafetyProfile
query_categorated$FEMAProfile
query_categorated$FDA_SPL_Profile
query_categorated$LOTUSProfile
query_categorated$PubChemClassifications
query_categorated$MeSHProfile
query_categorated$LiteratureProfile
query_categorated$ChemicalTerms
query_categorated$ChemicalTraits
query_categorated$ChemicalTraitOntology
query_categorated$ChemicalTraitMatrix
query_categorated$ChemicalTraitOntologyMatrix
query_categorated$ChemicalTraitEvidence
query_categorated$ChemicalTraitReport
query_categorated$ChemicalTraitSummary
query_categorated$ChemicalTraitSimilarity
query_categorated$ChemicalClasses
query_categorated$ChemicalMeasurements
query_categorated$ChemicalMeasurementSummary
query_categorated$ChemicalHazards
query_categorated$ChemicalUses
query_categorated$ChemicalBioassays
query_categorated$ChemicalBioactivities
query_categorated$ChemicalTargets
query_categorated$ChemicalPotencies
query_categorated$PubChemBioAssayDetails
query_categorated$ChemicalTaxonomy
query_categorated$ChemicalOccurrences
query_categorated$ChemicalPathwayRoles
query_categorated$KEGGReactionParticipants
query_categorated$KEGGPathways
query_categorated$KEGGClassifications
query_categorated$DerivedGroups
query_categorated$SourceCoverage
query_categorated$DataDictionary
query_categorated$TableQuality
query_categorated$SourceDiagnostics
query_categorated$ValidationIssues
query_categorated$ValidationSummary
```

The source-specific PubChem profile tables preserve the raw annotation text
while adding cleaned fields such as GHS hazard codes, FEMA/JECFA identifiers,
flavor or odor terms, FDA/SPL route and dosage-form terms, LOTUS occurrence
signals, PubChem classification-tree paths, MeSH pharmacologic actions, PubMed
IDs, DOI values, PubChem BioAssay activity summaries, and bounded BioAssay
description metadata when `detail = "full"` is used. The `assay_detail_limit`
argument controls how many BioAssay descriptions are fetched per query so broad
screens stay inspectable and polite to PubChem. `DerivedGroups` summarizes these
into analysis columns for metabolic, biomedical, ecological, sensory, safety,
bioactivity, and analytical context.

The normalized `Chemical*` tables go one step further and split evidence text
into discrete values for analysis. `ChemicalTerms` is a long-form vocabulary
table for grouping chemicals by sensory, biomedical, ecological, safety,
classification, and metabolic terms. `ChemicalTraits` unifies those cleaned
signals into one cross-source long table with `TraitType`, `TraitGroup`,
`TraitValue`, `SourceDatabase`, evidence, confidence, and a stable matrix key.
`ChemicalTraitOntology` maps whitelisted discrete traits into controlled domains
such as safety, physicochemical behavior, metabolism, ecology, sensory,
bioactivity, biomedical use, regulatory status, and reactivity. Every ontology
row keeps the source trait key, evidence text/URL, source database, extraction
rule, confidence, and source-backed identifiers where available: GHS hazard
codes, KEGG pathway/reaction/compound/EC IDs, MeSH tree numbers, NCBI Taxonomy
IDs, NCBI Gene IDs, PubChem BioAssay AIDs, and target accessions. `ChemicalTraitMatrix`
and `ChemicalTraitOntologyMatrix` convert curated high-value traits into binary
columns for filtering, clustering, ordination, heatmaps, and model inputs. Use
the helper functions to rebuild specialized matrices without rerunning web
requests:

``` r
core_matrix = query_categorated$ChemicalTraitMatrix
ontology = query_categorated$ChemicalTraitOntology
ontology_matrix = query_categorated$ChemicalTraitOntologyMatrix
ontology_evidence = query_categorated$ChemicalTraitEvidence
trait_report = query_categorated$ChemicalTraitReport
bioactivity_matrix = chemicalTraitMatrix(query_categorated, profile = "bioactivity")
kegg_matrix = chemicalTraitMatrix(query_categorated$ChemicalTraits, profile = "kegg")
confidence_matrix = chemicalTraitMatrix(
  query_categorated,
  profile = "full",
  mode = "confidence",
  min_confidence = "high",
  max_traits = 250
)
ontology_confidence_matrix = chemicalTraitOntologyMatrix(
  query_categorated,
  mode = "confidence",
  min_confidence = "medium",
  max_terms = 250
)

trait_summary = chemicalTraitSummary(query_categorated)
source_summary = chemicalTraitSummary(query_categorated, by = "source")
trait_similarity = chemicalTraitSimilarity(query_categorated, profile = "core")
report = chemicalTraitReport(query_categorated, min_confidence = "medium")
dictionary = chemicalDataDictionary("ChemicalTraitReport")
audit = validateCategorateResult(query_categorated)
audit$Summary
audit$TableQuality
audit$SourceDiagnostics
audit$Issues
export_manifest = exportCategorateWorkbook(
  query_categorated,
  "categorate_export",
  format = "csv",
  overwrite = TRUE
)
export_manifest
h319_evidence = chemicalTraitEvidence(
  query_categorated,
  keys = "safety__ghs_hazard_code__h319"
)
trait_evidence = chemicalTraitEvidence(
  query_categorated,
  keys = "hazard__hazard_code__h319",
  type = "trait"
)
```

`ChemicalMeasurements` stores numeric
properties and extracted experimental/toxicity values with raw units plus
standardized values, canonical units, measurement classes, relations, and
behavior bins. `ChemicalMeasurementSummary` condenses those values into one row
per compound-property combination for plotting and filtering. `ChemicalHazards`
separates GHS codes, hazard groups, routes, target organs, precaution codes, and
toxicity metrics. `ChemicalBioassays`, `ChemicalBioactivities`, `ChemicalTargets`,
and `ChemicalPotencies` split PubChem assay summaries into assay outcomes,
activity domains, target identifiers/names, target organism/taxonomy fields,
assay endpoint names, activity directions, and numeric potency values.
`PubChemBioAssayDetails` preserves the raw selected assay-description metadata
behind those normalized bioactivity rows. `ChemicalTaxonomy` extracts organism,
genus, species, family,
order, class, phylum, kingdom, cleaned taxonomy terms, and natural-product
class. `ChemicalClasses` includes discrete natural-product superclass, class,
and subclass rows from LOTUS chemical classification trees. `ChemicalOccurrences`
gives one row per chemical-organism occurrence with taxonomic ranks, source
evidence, occurrence type, and confidence, while `DerivedGroups` summarizes
natural-product classes, bioactivity breadth, occurrence breadth, dominant
kingdom/family, and plant/fungal/bacterial occurrence flags. `ChemicalPathwayRoles` and
`KEGGReactionParticipants` split KEGG pathways, enzymes, reactions, substrates,
and products into analysis-ready rows.

Every detailed `categorate()` result also carries a data-quality layer.
`DataDictionary` describes expected tables, columns, types, required fields,
roles, and allowed controlled values. `TableQuality` reports row counts,
required-column completeness, duplicate analysis-key counts, confidence counts,
and table status. `SourceDiagnostics` shows which source/query combinations
returned usable rows. `ValidationSummary` and `ValidationIssues` provide the
same audit generated by `validateCategorateResult()`, so downstream scripts can
stop early when a schema or extraction problem appears.

Use `exportCategorateWorkbook()` to share results outside R. It writes a
manifest plus clean CSV files by default, and can write an `.xlsx` workbook when
`openxlsx` or `writexl` is installed.

Large plant-chemistry projects often run `categorate()` enrichment in resumable
batches. Use the batch-analysis helpers to audit those outputs before handing
them to a downstream statistics or figure-making workspace:

``` r
batches = readCategorateBatchDirectory(
  "dsi_categorate_research_pubchem_cid_20260611/categorate_batches"
)

batch_summary = summarizeCategorateBatches(batches)
batch_summary[, c("BatchIndex", "Status", "ValidationStatus",
                  "ResolvedCIDCount", "PubChemPropertyRatio")]

combined = combineCategorateTables(
  batches,
  tables = c("ChemicalTraits", "PubChemProperties", "ChemicalTraitEvidence")
)

manifest = exportPlantChemistryAnalysisBundle(
  categorate_batches = batches,
  path = "plant_chemistry_analysis_bundle",
  plant_membership = "dsi_species_compound_membership.csv",
  species_pair_tanimoto = "dsi_species_pair_tanimoto_summary.csv",
  resolved_compounds = "dsi_resolved_compounds_for_categorate.csv",
  pubchem_fingerprints = "dsi_pubchem_fingerprints.csv",
  plant_compound_pair_tanimoto = "dsi_plant_compound_pair_tanimoto.csv.gz",
  file_references = c(
    plant_compound_pairs = "dsi_plant_compound_pair_tanimoto.csv.gz",
    compound_pairs = "dsi_compound_pair_tanimoto.csv.gz",
    chemical_traits_full = "combined_chemical_traits_if_exported_separately.csv.gz",
    trait_evidence_full = "combined_trait_evidence_if_exported_separately.csv.gz"
  ),
  plant_list = "project_species.csv",
  project_id = "plant_project",
  tables = c("ChemicalTraitSummary", "DerivedGroups", "PubChemProperties",
             "SourceCoverage", "ValidationSummary", "ValidationIssues"),
  format = "csv",
  include_comparable_tanimoto = TRUE,
  overwrite = TRUE
)
manifest

validation = validatePlantChemistryAnalysisBundle(
  "plant_chemistry_analysis_bundle",
  use_python = TRUE,
  use_pandas = TRUE
)
validation$Summary
```

The exported `BatchSummary` should be checked first. Batches with `Status =
"error"` should not be interpreted; batches with `Status = "incomplete"` need
review because they have fewer PubChem property rows than resolved CIDs at the
chosen threshold. Large compressed pairwise files can be recorded in
`FileReferences` instead of copied into every bundle. For very large runs, keep
row-level tables such as `ChemicalTraits` and `ChemicalTraitEvidence` as
separate targeted exports unless the downstream analysis needs a full combined
copy. Very wide sparse tables such as `ChemicalTraitMatrix` are useful for
focused modeling, but they can be slow and awkward to move for thousands of
compounds; keep the batch directory in `FileReferences` and extract a focused
matrix profile when the analysis plan is settled.

For CSV bundles, `exportPlantChemistryAnalysisBundle()` finalizes the handoff by
default. The finalized bundle includes `00_DataDictionary.csv`,
`03b_PlantCompoundMembershipEnriched.csv`, `11b_SourceCoverageSummary.csv`,
`12b_ValidationOverview.csv`, `14_PlantChemistrySummary.csv`,
`15_PlantChemistryMissingSpecies.csv`, `16_EvidenceGradeSummary.csv`,
`17_ReviewRequiredOccurrences.csv`, model-ready feature matrices
(`18_` through `21_` plus evidence/context matrices as `24_` through `28_`),
`README.md`, and `METHODS_TEXT.md`. If
`include_comparable_tanimoto = TRUE` and a plant-compound pair Tanimoto table is
supplied, the bundle also includes scope- and group-filtered comparable
Tanimoto summaries as `22_ComparableScopeTanimotoSummary.csv` and
`23_ComparableGroupTanimotoSummary.csv`.
The enriched membership table joins occurrence evidence to PubChem fingerprints,
PubChem/categorate properties, natural-product context, source coverage, and the
uafR comparable-chemistry fields (`comparison_scope`, `comparison_group`,
`comparison_subgroup`, `comparability_confidence`, `comparability_basis`,
`comparable_for_matrix`, and `comparison_caveat`). These fields are intended to
help downstream projects compare like with like; they do not turn public
database records into confirmed sample measurements.

For reusable project handoffs, use `runPlantChemistryProject()` or the CLI
wrapper. The runner accepts a plant list plus either curated plant-compound
records, a saved plant phytochemistry result, cached categorate batches, or an
explicit live-discovery request. It writes a project manifest and a finalized
bundle without hardcoding project-specific paths:

``` r
project = runPlantChemistryProject(
  plant_list = "project_species.csv",
  metadata = "project_plant_metadata.csv",
  plant_compounds = "curated_or_discovered_plant_compounds.csv",
  output_dir = "plant_chemistry_project_run",
  plant_compound_pair_tanimoto = "plant_compound_pair_tanimoto.csv.gz",
  project_id = "plant_project",
  overwrite = TRUE
)

project$Project
validatePlantChemistryAnalysisBundle(project$Project$bundle_dir)
```

``` sh
Rscript tools/run_plant_chemistry_project.R \
  --plant-list project_species.csv \
  --metadata project_plant_metadata.csv \
  --plant-compounds curated_or_discovered_plant_compounds.csv \
  --plant-compound-pair-tanimoto plant_compound_pair_tanimoto.csv.gz \
  --out-dir plant_chemistry_project_run \
  --project-id plant_project \
  --overwrite true
```

Use `finalizePlantChemistryAnalysisBundle()` when a bundle has already been
written and only needs the analysis-ready handoff files refreshed. Use
`validatePlantChemistryAnalysisBundle()` before handing a bundle to another
project. It checks manifest row/column counts, required columns, accidental row
index columns, portable artifact paths, file sizes, MD5 checksums, nested
feature-manifest file references, identical species IDs and row order across
feature matrices, and CSV parser consistency. The manifest intentionally omits
only its own checksum to avoid a circular hash. Optional Python and pandas
checks can be enabled on machines where those tools are available.

Schema `1.0.0` introduces portable manifest paths and checksums. Older bundles
that lack those fields should be treated as legacy artifacts: rerun
`finalizePlantChemistryAnalysisBundle()` with the current package before
handoff, then require a passing validation result.

For model-ready exports outside a full bundle, `exportPlantChemistryFeatureSet()`
can write count, binary, fraction, and confidence-weighted matrices:

``` r
project_species = read.csv("project_species.csv")$species
project_plant_metadata = read.csv("project_plant_metadata.csv")

features = exportPlantChemistryFeatureSet(
  membership = enriched_membership,
  species_universe = project_species,
  plant_metadata = project_plant_metadata,
  modes = c("count", "binary", "fraction", "confidence"),
  path = "plant_feature_set",
  overwrite = TRUE
)

features$Manifest
features$SpeciesMetadata
```

Every feature matrix uses the same `species_id` rows and order. Species in the
declared universe with no supplied chemistry records are retained with zero
feature values and `chemistry_record_status = "no_records_in_membership"`.
Those zeros mean that no records were supplied to the feature builder; they do
not demonstrate biological absence. Context fields are similarly conservative:
`biological_context_known` requires a reported plant part or tissue,
`analytical_method_known` requires an analytical method, and
`source_provenance_record` records generic database/literature provenance
separately.

Evidence and classification helpers support reproducible subsetting and
review:

``` r
plantOccurrenceEvidenceDictionary()
chemistryComparisonDictionary()

direct = filterPlantEvidenceDirect(enriched_membership)
comparable = filterPlantEvidenceComparable(enriched_membership)
review = filterPlantEvidenceReviewRequired(enriched_membership)
```

Use `standardizeChemistryClassificationOverrides()` and
`applyChemistryClassificationOverrides()` when a project needs source-backed or
human-reviewed chemistry-scope corrections. Overrides use the canonical values
listed by `chemistryComparisonDictionary()`; common legacy aliases are
normalized, but invalid scope/group combinations fail explicitly. Unknown,
broad/uncertain, and xenobiotic/contaminant chemistry is retained for audit but
excluded from comparable matrices by default.

`keggProfile()` can also be used directly when the goal is KEGG-specific
annotation. It resolves names or supplied KEGG IDs, parses KEGG flat-file
records, follows KEGG links, and returns pathways, reactions, enzymes, modules,
identifiers, and reproducible pathway/enzyme classifications. Name searches
accept only exact normalized KEGG synonyms for enrichment. Broader substring
results are preserved with `Accepted = "No"` and an explicit reason in
`search_candidates` (and `KEGGSearchCandidates` in categorate results), so a
result such as `1-Hexanol` matching `2-Ethylhexan-1-ol` cannot silently acquire
the wrong biochemical context.

``` r
kegg_profile = keggProfile(
  c("aspirin", "glucose"),
  max_matches_per_query = 3,
  link_targets = c("pathway", "reaction", "enzyme")
)

kegg_profile$pathways
kegg_profile$search_candidates
kegg_profile$reactions
kegg_profile$enzymes
kegg_profile$classifications
```

Linked KEGG IDs are resolved into names and definitions by default, so tables
such as `kegg_profile$reactions` include reaction names, definitions, and
equations when KEGG exposes them.

### Package release checks

Use the release-check wrapper before sharing a source archive, building a
student bundle, or opening a pull request:

``` sh
Rscript tools/check_package_release.R
```

The script builds a source tarball outside the repository and runs
`R CMD check --no-manual` against the tarball. That is the production check
path. Running `R CMD check .` directly on the live checkout can report local
artifacts such as `.git`, `.DS_Store`, `.Rhistory`, or generated check
directories that are not included in the source package.

By default the wrapper also installs the built tarball into a clean temporary
library and checks that core production helpers load. Useful release options
include:

``` sh
Rscript tools/check_package_release.R --out-dir uafR_release_check
Rscript tools/check_package_release.R --skip-install-check
Rscript tools/check_package_release.R --student-bundle student_bundle/uafR_student_bundle_VERSION
```

After the worktree is clean, create a release-candidate manifest beside the
built tarball:

``` sh
Rscript tools/check_package_release.R \
  --out-dir /private/tmp/uafR_release_check \
  --release-manifest /private/tmp/uafR_release_check/uafR_release_manifest.json
```

Manifest generation refuses a dirty worktree. It records the exact Git commit,
package and R versions, tarball byte size, MD5 and SHA-256 hashes, package-test
result, R CMD check result, and clean-library installation result.

Live PubChem/NCI integration tests are opt-in. Set `UAFR_RUN_LIVE_TESTS=true`
before running tests, or pass `--run-live-tests` to
`tools/check_package_release.R`, when you want to exercise live web calls.

### Offline plant chemistry example

The package includes a simulated offline fixture under
`inst/extdata/offline_plant_chemistry`. It is for workflow testing and teaching,
not real phytochemical evidence. It exercises direct database-style evidence,
fallback context, candidate-only literature rows, unresolved compounds,
comparable chemistry, and small Tanimoto summaries without querying live
services.

From the repository root:

``` sh
Rscript tools/build_offline_plant_chemistry_example.R uafR_offline_example
```

The output is a finalized plant chemistry analysis bundle that can be opened in
RStudio or inspected as CSV files.

## Combined Mass Spectrometry + Cheminformatics Workflow

``` r
query_chemicals = input_dat$Compound.Name[input_dat$Match.Factor > 70]
query_categorated = categorate(query_chemicals, chem_library, input_format = "wide")

## example of using the info from categorate() to get a user-defined set of chemicals with exactoThese():
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "All")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "reactives")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "LOTUS")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "KEGG")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "FEMA")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "FDA_SPL")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = c("reactives", "FEMA"))
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Greater Than", subset_input = 125)
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Less Than", subset_input = 205)
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Between", subset_input = c(125, 200))
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "Rings", subsetArgs2 = "Greater Than", subset_input = 1)
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "Groups", subsetArgs2 = "Greater Than", subset_input = 2)
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "Atoms", subsetArgs2 = "Greater Than", subset_input = 6)
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "NCharges", subsetArgs2 = "Greater Than", subset_input = 2)
these_chems = exactoThese(query_categorated, subsetBy = "Library", subsetArgs = "GroupB")

input_exacto = mzExacto(input_spread, these_chems)
```
