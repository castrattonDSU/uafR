# Start Here: uafR Student Bundle

Use this folder to install and test uafR on a student machine.

## What to Do First

1. Install the current version of R from <https://cran.r-project.org/>.
2. Install the current version of RStudio Desktop from <https://posit.co/download/rstudio-desktop/>.
3. Unzip this whole bundle folder. Do not run scripts from inside the zip preview.
4. Open RStudio.
5. In RStudio, open `verify_bundle_integrity.R` from this folder and click **Source**.
6. Open `preflight_check.R` from this folder and click **Source**.
7. If preflight has no failed checks, open `install_uafR_from_bundle.R` and click **Source**.
8. When installation finishes, open `run_student_acceptance_test.R` and click **Source**.
9. Keep `uafR_QUICK_REFERENCE.md` nearby and open `training/uafR_training_manual.pdf`.

## What Each Script Does

| Script | When to run it | What it checks or changes |
|---|---|---|
| `verify_bundle_integrity.R` | After unzipping | Confirms copied files match the generated `CHECKSUMS.csv` file. |
| `preflight_check.R` | Before installing | Checks R, RStudio, internet access, the bundle files, the R library, and required dependencies. |
| `install_uafR_from_bundle.R` | First install | Installs dependencies, installs uafR from `packages/`, and verifies the result. |
| `verify_uafR_install.R` | Any time | Confirms that uafR loads and package data are available. |
| `run_student_acceptance_test.R` | After install | Runs the offline workflow used in training and confirms the package is ready for class. |
| `update_uafR_from_bundle.R` | When a newer bundle is distributed | Reinstalls uafR from the newer local archive and verifies the result. |

## Expected Success Message

The final check should end with:

```text
Acceptance test result: PASS.
```

If a script fails, read the message printed at the bottom of the Console. Each
script also writes a log file in this folder.

## Quick Reference

Open `uafR_QUICK_REFERENCE.md` for a compact list of setup commands, the
offline training workflow, the live database smoke test, and common fixes.

## Optional Live Database Test

The acceptance test avoids live database queries by default. To include the
live `categorate()` smoke test, run this in the RStudio Console from the bundle
folder:

```r
Sys.setenv(UAFR_STUDENT_LIVE = "1")
source("run_student_acceptance_test.R")
```

The live test needs internet access and can take longer than the offline checks.
