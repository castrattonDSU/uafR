# uafR Student Bundle Install

This bundle installs uafR from a local package archive distributed by the DSU
dsDNA Core program. uafR is also publicly available from its GitHub repository;
this bundle provides a fixed, auditable version that can be installed without
downloading uafR itself during setup.

The installer uses CRAN and Bioconductor only for required dependencies. The
uafR package itself is installed from the local archive in `packages/`.

## Quick Start

1. Unzip the whole bundle folder.
2. Open RStudio.
3. Open `START_HERE.md` for the short checklist.
4. Source `verify_bundle_integrity.R`.
5. Source `preflight_check.R`.
6. If preflight has no failed checks, source `install_uafR_from_bundle.R`.
7. Source `run_student_acceptance_test.R`.
8. Keep `uafR_QUICK_REFERENCE.md` nearby and open `training/uafR_training_manual.pdf`.

Do not run scripts from inside a zip preview window. The folder must be fully
unzipped so R can see the `packages/`, `training/`, and `examples/` folders.

## Confirm the Bundle Was Copied Correctly

The generated bundle includes `CHECKSUMS.csv`, which lists the file size and
MD5 checksum for every file in the bundle at build time. After unzipping, run:

```r
source("verify_bundle_integrity.R")
```

When it finishes, the Console should say:

```text
Bundle integrity result: PASS.
```

## Install

In RStudio, open this file from the unzipped bundle folder:

```text
install_uafR_from_bundle.R
```

Click **Source**. The installer will print five steps:

```text
Step 1 of 5: Checking the bundle and R library
Step 2 of 5: Installing CRAN dependencies
Step 3 of 5: Installing Bioconductor dependencies
Step 4 of 5: Installing uafR from the local bundle archive
Step 5 of 5: Verifying the installation
```

When it finishes, the Console should say:

```text
Install complete.
```

The installer writes `uafR_install_log.txt` in this folder.

## Verify

To verify the installation again later, source:

```r
source("verify_uafR_install.R")
```

The verification script does not query live web services. It checks package
loading, required dependencies, bundled data, and the saved standard
categorate result.

## Run the Student Acceptance Test

After installation, source:

```r
source("run_student_acceptance_test.R")
```

This runs the offline workflow used in training. It should end with:

```text
Acceptance test result: PASS.
```

To include the optional live database smoke test:

```r
Sys.setenv(UAFR_STUDENT_LIVE = "1")
source("run_student_acceptance_test.R")
```

## Quick Reference

Open `uafR_QUICK_REFERENCE.md` for a compact command list covering setup,
offline training, the live smoke test, updates, and common troubleshooting.

## Update From a New Bundle

When a newer bundle is distributed, unzip the new bundle and source:

```r
source("update_uafR_from_bundle.R")
```

The update script reinstalls uafR from the new local archive and verifies the
result.

## Basic Use After Install

```r
library(uafR)
data("library_data", package = "uafR")
```

For a small live database enrichment test:

```r
quick_result <- categorate(
  compounds = c("aspirin", "caffeine"),
  chemical_library = library_data,
  input_format = "wide",
  detail = "research",
  cache = TRUE,
  throttle = 0.2,
  assay_detail_limit = 0
)

validateCategorateResult(quick_result)$Summary
quick_result$ChemicalTraitReport
```

## Troubleshooting

| Problem | Likely cause | What to do |
|---|---|---|
| Bundle integrity check fails. | One or more files changed, did not copy, or did not unzip correctly. | Delete the unzipped folder and unzip the original bundle again. If it still fails, get a fresh copy of the zip. |
| The script cannot find the package archive. | The bundle was not fully unzipped, or the `packages/` folder was moved. | Unzip the whole bundle again and keep the folder structure unchanged. |
| Preflight says the R library is not writable. | R cannot install packages into the first library path. | Restart RStudio. If it still fails, ask for help setting a user library. |
| Bioconductor packages fail to install. | Internet access, R version, or Bioconductor access is blocked. | Run `preflight_check.R`, update R if it is old, and try again on a network that can reach Bioconductor. |
| R asks about updating many packages. | RStudio or Bioconductor is offering broad updates. | Choose the option that does not update unrelated packages. |
| `library(uafR)` fails right after install. | RStudio has not refreshed the library. | Restart RStudio, then run `library(uafR)`. |
| A live database query is slow. | Live PubChem, KEGG, LOTUS, or other services are responding slowly. | Start with one or two compounds, keep `assay_detail_limit = 0`, use caching, and scale up after the smoke test works. |
