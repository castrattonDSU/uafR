# Offline Plant Phytochemistry Workflow

## Purpose

This article runs a complete plant-chemistry project without web access.
The included rows are a **simulated teaching fixture**. Their source
labels, citations, occurrence relationships, and project context are
synthetic and must not be treated as evidence about the listed plants.

The fixture exists to demonstrate the package contracts:

- one stable plant universe;
- direct, fallback, candidate, and review-required evidence;
- resolved and unresolved compounds;
- comparable and non-comparable chemistry;
- plant-labeled Tanimoto rows;
- model-ready feature matrices; and
- a validated analysis bundle.

## Load the fixture

``` r

example_dir = system.file(
  "extdata",
  "offline_plant_chemistry",
  package = "uafR"
)
stopifnot(nzchar(example_dir))

plants = read.csv(
  file.path(example_dir, "plant_list.csv"),
  stringsAsFactors = FALSE
)
metadata = read.csv(
  file.path(example_dir, "plant_metadata.csv"),
  stringsAsFactors = FALSE
)
plant_compounds = read.csv(
  file.path(example_dir, "plant_compounds.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
compound_pairs = read.csv(
  file.path(example_dir, "plant_compound_pair_tanimoto.csv"),
  stringsAsFactors = FALSE
)

data.frame(
  plants = nrow(plants),
  metadata_rows = nrow(metadata),
  occurrence_rows = nrow(plant_compounds),
  pair_rows = nrow(compound_pairs)
)
#>   plants metadata_rows occurrence_rows pair_rows
#> 1      6             6              12         5
```

Inspect the evidence and identity fields before building outputs:

``` r

plant_compounds[, c(
  "species",
  "compound_name",
  "source_database",
  "evidence_tier",
  "matched_rank",
  "confidence",
  "InChIKey",
  "comparison_scope",
  "comparable_for_matrix",
  "curation_flag"
)]
#>              species            compound_name source_database
#> 1        Salix nigra                  salicin       LOTUS_SIM
#> 2        Salix nigra                quercetin      PubMed_SIM
#> 3           Zea mays             ferulic acid       LOTUS_SIM
#> 4           Zea mays                  glucose     PubChem_SIM
#> 5  Camellia sinensis                 caffeine       LOTUS_SIM
#> 6  Camellia sinensis epigallocatechin gallate       LOTUS_SIM
#> 7    Mentha piperita                  menthol       LOTUS_SIM
#> 8    Mentha piperita                 limonene       LOTUS_SIM
#> 9    Brassica juncea                 sinigrin       LOTUS_SIM
#> 10   Brassica juncea             lead complex    PubTator_SIM
#> 11       Salix nigra           willow extract    PubTator_SIM
#> 12          Zea mays benzoxazinoid derivative    KNApSAcK_SIM
#>                        evidence_tier matched_rank confidence
#> 1            direct_species_database      species       high
#> 2          direct_species_literature      species     medium
#> 3            direct_species_database      species       high
#> 4            direct_species_database      species     medium
#> 5            direct_species_database      species       high
#> 6            direct_species_database      species       high
#> 7            direct_species_database      species       high
#> 8            direct_species_database      species       high
#> 9            direct_species_database      species       high
#> 10 direct_species_pubtator_candidate      species        low
#> 11 direct_species_pubtator_candidate      species        low
#> 12           genus_database_fallback        genus     medium
#>                       InChIKey                 comparison_scope
#> 1  NGFMICBWJRZIBI-UHFFFAOYSA-N          specialized_metabolites
#> 2  REFJWTPEDVJJIY-UHFFFAOYSA-N          specialized_metabolites
#> 3  KSEBMYQBYZTDHS-UHFFFAOYSA-N          specialized_metabolites
#> 4  WQZGKKKJIJFFOK-UHFFFAOYSA-N              primary_metabolites
#> 5  RYYVLZVUVIJVGH-UHFFFAOYSA-N          specialized_metabolites
#> 6  WMBWREPUVVBILR-UHFFFAOYSA-N          specialized_metabolites
#> 7  NOOLISFMXDJSKH-UHFFFAOYSA-N volatile_specialized_metabolites
#> 8  XMGQYMWWDOXHJM-UHFFFAOYSA-N volatile_specialized_metabolites
#> 9  XTQHRLGUOYLDGD-UHFFFAOYSA-N          specialized_metabolites
#> 10                                    xenobiotic_or_contaminant
#> 11                                                      unknown
#> 12                                      specialized_metabolites
#>    comparable_for_matrix                       curation_flag
#> 1                    Yes                   simulated_fixture
#> 2                    Yes                   simulated_fixture
#> 3                    Yes                   simulated_fixture
#> 4                    Yes                   simulated_fixture
#> 5                    Yes                   simulated_fixture
#> 6                    Yes                   simulated_fixture
#> 7                    Yes                   simulated_fixture
#> 8                    Yes                   simulated_fixture
#> 9                    Yes                   simulated_fixture
#> 10                    No simulated_candidate_review_required
#> 11                    No simulated_candidate_review_required
#> 12                    No          simulated_fallback_context
```

## Build and validate a project

[`runPlantChemistryProject()`](https://castrattonDSU.github.io/uafR/reference/runPlantChemistryProject.md)
is the offline/cached project coordinator. Live discovery is disabled
unless `run_discovery = TRUE` is supplied explicitly.

``` r

project_dir = tempfile("uafr_offline_project_")

project = runPlantChemistryProject(
  plant_list = plants,
  metadata = metadata,
  plant_compounds = plant_compounds,
  plant_compound_pair_tanimoto = compound_pairs,
  output_dir = project_dir,
  project_id = "simulated_offline_example",
  overwrite = TRUE
)

bundle_dir = project$Project$bundle_dir[[1]]
validation = validatePlantChemistryAnalysisBundle(bundle_dir)
stopifnot(identical(
  validation$Summary$ExportReadyStatus[[1]],
  "pass"
))

validation$Summary
#>   CSVFileCount CSVPassCount CSVWarnCount CSVFailCount RequiredColumnFailCount
#> 1           23           23            0            0                       0
#>   ArtifactWarnCount ArtifactFailCount ManifestReferenceFailCount
#> 1                 0                 0                          0
#>   PythonCsvChecked PandasChecked ExportReadyStatus
#> 1               No            No              pass
```

The validator checks CSV parsing, required columns, artifact
relationships, manifest references, and the bundle’s export-ready
status. A passing bundle is structurally ready for handoff; it does not
make the simulated evidence real.

## Inspect evidence and review queues

``` r

membership = read.csv(
  file.path(bundle_dir, "03b_PlantCompoundMembershipEnriched.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

grade_summary = read.csv(
  file.path(bundle_dir, "16_EvidenceGradeSummary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

direct = filterPlantEvidenceDirect(membership)
comparable = filterPlantEvidenceComparable(membership)
review_required = filterPlantEvidenceReviewRequired(membership)

grade_summary
#>                               evidence_grade evidence_grade_rank
#> 1             direct_species_database_record                   1
#> 2 direct_species_literature_supported_record                   2
#> 3             pubtator_pubmed_candidate_only                   4
#>                         evidence_grade_label occurrence_count species_count
#> 1             Direct species database record                9             5
#> 2 Direct species literature-supported record                1             1
#> 3             PubTator/PubMed candidate only                2             2
#>   compound_count source_backed_count structure_resolved_count comparable_count
#> 1              8                   9                        0                7
#> 2              1                   1                        0                1
#> 3              0                   0                        0                0
#>   review_required_count
#> 1                     0
#> 2                     0
#> 3                     2
#>                                                               recommended_use
#> 1  Useful for occurrence summaries but not structure-based Tanimoto analyses.
#> 2  Useful for occurrence summaries but not structure-based Tanimoto analyses.
#> 3 Use only for literature triage until a source-backed occurrence is curated.
data.frame(
  all_rows = nrow(membership),
  direct_structure_resolved = nrow(direct),
  direct_comparable = nrow(comparable),
  review_required = nrow(review_required)
)
#>   all_rows direct_structure_resolved direct_comparable review_required
#> 1       12                         0                 0               2
```

These filters do not delete the original evidence. They make the
analysis rule explicit while the bundle preserves fallback, candidate,
unresolved, and excluded rows for audit.

## Inspect model-ready features

``` r

feature_manifest = read.csv(
  file.path(bundle_dir, "28_FeatureMatrixManifest.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

feature_manifest[, c(
  "Table",
  "FileName",
  "Mode",
  "RowCount",
  "ColumnCount"
)]
#>                        Table                                 FileName     Mode
#> 1 ComparisonGroupCountMatrix 18_FeatureComparisonGroupCountMatrix.csv    count
#> 2 ComparisonScopeCountMatrix 19_FeatureComparisonScopeCountMatrix.csv    count
#> 3  SourceCoverageCountMatrix       20_FeatureSourceCoverageMatrix.csv    count
#> 4   EvidenceGradeCountMatrix   24_FeatureEvidenceGradeCountMatrix.csv    count
#> 5       PlantPartCountMatrix       25_FeaturePlantPartCountMatrix.csv    count
#> 6          TissueCountMatrix          26_FeatureTissueCountMatrix.csv    count
#> 7          MethodCountMatrix          27_FeatureMethodCountMatrix.csv    count
#> 8            SpeciesMetadata            21_FeatureSpeciesMetadata.csv metadata
#>   RowCount ColumnCount
#> 1        6           8
#> 2        6           5
#> 3        6           7
#> 4        6           5
#> 5        6           2
#> 6        6           2
#> 7        6           3
#> 8        6          16

feature_files = file.path(bundle_dir, feature_manifest$FileName)
species_order = lapply(feature_files, function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )$species
})

stopifnot(all(vapply(
  species_order[-1],
  identical,
  logical(1),
  species_order[[1]]
)))
```

Feature matrices are analysis inputs. They are not evidence of
mechanism, efficacy, remediation performance, exposure, pathway
activity, or chemistry measured in a project sample.

## Reproduce outside the vignette

From a repository checkout:

``` sh
Rscript tools/build_offline_plant_chemistry_example.R \
  --out-dir offline-plant-example
```

Record `packageVersion("uafR")`,
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html), the input
checksums, and the generated manifest with any downstream analysis.
