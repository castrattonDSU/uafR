# Manual Upgrade Report

## Baseline Audit

Date: 2026-05-27

Working directory: local working copy of the `uafR` repository.

Baseline Git state at start of upgrade:

- Existing local edits were present from the student bundle hardening pass.
- `training/assets/dsdna_core_logo.png` was already present and untracked at the start of this upgrade.
- Generated LaTeX scratch files were present before editing and were excluded from content review.

Main LaTeX entry point:

- `training/main.tex`

Included LaTeX source groups:

- Front matter: `training/chapters/preface.tex`
- Program architecture: `training/chapters/program_architecture.tex`
- Daily schedule: `training/chapters/daily_plans.tex`
- Weekly chapters: `training/chapters/week01.tex` through `training/chapters/week10.tex`
- Labs: `training/chapters/lab_compendium.tex`
- Readings/exercises: `training/chapters/readings_exercises.tex`
- Project tracks: `training/chapters/project_tracks.tex`
- Worksheets: `training/worksheets/worksheets.tex`
- Rubrics: `training/rubrics/rubrics.tex`
- Instructor/program notes: `training/program_notes/program_implementation_notes.tex`
- Review guidance: `training/program_notes/dsdna_core_review_notes.tex`
- Installation appendix: `training/chapters/appendix_installation.tex`
- Table/script guide: `training/chapters/appendix_tables.tex`

Style and build support:

- `training/preamble.tex`
- `training/bibliography/references.bib`
- `training/assets/dsdna_core_logo.png`
- `training/scripts/*.R`
- `training/data/*.csv`

Baseline build command:

```sh
cd training
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Baseline build status:

- Build succeeded.
- `latexmk` reported the PDF was already up to date.
- Baseline PDF: 89 pages.
- `training/main.log` showed only font-substitution warnings; no unresolved references or citations were detected in the baseline log check.

Baseline R/training smoke status:

```sh
Rscript training/scripts/run_training_smoke_tests.R
```

- Completed successfully.
- The default library did not have installed `uafR`, so package-dependent offline workflow checks were skipped in that smoke-test run.
- `devtools::load_all(".")` successfully loaded the package and confirmed exported function availability.

## Current Manual Inventory

Usable student-facing detail already present:

- Title page and program promise.
- Program outcomes and final package overview.
- High-level weekly chapters for Weeks 1-10.
- A 10-week weekday table.
- Initial lab compendium with 10 labs.
- Initial worksheets and rubrics.
- Installation appendix and bundle scripts.
- Project tracks for six broad interest areas.

Skeletal or insufficient sections at baseline:

- The day-by-day schedule covered only weekdays and used table summaries rather than full daily modules.
- Weekend required work was not fully represented.
- Daily objectives, prerequisites, time budgets, materials, log prompts, deliverables, quality gates, troubleshooting, and advanced challenges were not present for all 70 days.
- Lab compendium lacked the requested 16-lab production bank and detailed answer guidance.
- Worksheets were printable but too thin for actual student artifacts across all required stages.
- Rubrics did not yet cover all requested categories with point weights.
- Instructor notes needed a more complete facilitator guide, decision-gate scripts, and support for students with different backgrounds.
- Project tracks needed stronger examples of good/bad questions, outputs, caveats, minimum viable products, and advanced paths.
- Final program completion standard was implicit rather than auditable.

Search findings:

- No standard unfinished-work markers were found in source text.
- Existing use of `mentor` and `program lead` is intentional for instructor-facing and review sections.
- Logo fallback language is limited to build instructions and does not appear as student-facing unfinished content.

## Verified uafR Reference List

Verified exported functions from `NAMESPACE`, source files, Rd docs, and training scripts:

- `spreadOut()`
- `mzExacto()`
- `exactoThese()`
- `categorate()`
- `validateCategorateResult()`
- `chemicalDataDictionary()`
- `chemicalTraitMatrix()`
- `chemicalTraitEvidence()`
- `chemicalTraitOntology()`
- `chemicalTraitOntologyMatrix()`
- `chemicalTraitReport()`
- `chemicalTraitSummary()`
- `chemicalTraitSimilarity()`
- `chemicalMeasurementSummary()`
- `exportCategorateWorkbook()`
- `pubchemProfile()`
- `keggProfile()`
- `personalLib()`
- `standardifyIt()`

Verified package datasets referenced in training materials:

- `library_data`
- `standard_data`
- `standard_spread`
- `standard_exacto`
- `standard_categorated`

Verified training scripts:

- `training/scripts/00_install_check.R`
- `training/scripts/01_project_setup.R`
- `training/scripts/04_core_workflow.R`
- `training/scripts/05_categorate_research.R`
- `training/scripts/06_trait_matrix_analysis.R`
- `training/scripts/07_visualization.R`
- `training/scripts/09_reproducibility_check.R`
- `training/scripts/run_training_smoke_tests.R`

Verified student bundle scripts:

- `verify_bundle_integrity.R`
- `preflight_check.R`
- `install_uafR_from_bundle.R`
- `verify_uafR_install.R`
- `run_student_acceptance_test.R`
- `update_uafR_from_bundle.R`
- `tools/build_student_bundle.R`

Verified `categorate(detail = "research")` output and validation table names from `chemicalDataDictionary()`, Rd docs, source files, and tests:

- `PubChemIdentity`
- `PubChemProperties`
- `PubChemSynonyms`
- `PubChemAnnotations`
- `PubChemSourceAnnotations`
- `PubChemSpectra`
- `PubChemSafety`
- `PubChemExperimental`
- `PubChemBioactivity`
- `PubChemBioAssayDetails`
- `PubChemIdentifiers`
- `SafetyProfile`
- `FEMAProfile`
- `FDA_SPL_Profile`
- `LOTUSProfile`
- `PubChemClassifications`
- `MeSHProfile`
- `LiteratureProfile`
- `ChemicalTerms`
- `ChemicalTraits`
- `ChemicalTraitOntology`
- `ChemicalTraitMatrix`
- `ChemicalTraitOntologyMatrix`
- `ChemicalTraitEvidence`
- `ChemicalTraitReport`
- `ChemicalTraitSummary`
- `ChemicalTraitSimilarity`
- `ChemicalClasses`
- `ChemicalMeasurements`
- `ChemicalMeasurementSummary`
- `ChemicalHazards`
- `ChemicalUses`
- `ChemicalBioassays`
- `ChemicalBioactivities`
- `ChemicalTargets`
- `ChemicalPotencies`
- `ChemicalTaxonomy`
- `ChemicalOccurrences`
- `ChemicalPathwayRoles`
- `KEGGReactionParticipants`
- `KEGGMatches`
- `KEGGRecords`
- `KEGGIdentifiers`
- `KEGGPathways`
- `KEGGReactions`
- `KEGGEnzymes`
- `KEGGModules`
- `KEGGLinks`
- `KEGGLinkMetadata`
- `KEGGClassifications`
- `SourceCoverage`
- `DerivedGroups`
- `Provenance`
- `DataDictionary`
- `TableQuality`
- `SourceDiagnostics`
- `ValidationIssues`
- `ValidationSummary`

Validation naming note:

- Saved `categorate()` results include `ValidationSummary`, `TableQuality`, `SourceDiagnostics`, `ValidationIssues`, and `DataDictionary`.
- Direct `validateCategorateResult(result)` calls return `Summary`, `TableQuality`, `SourceDiagnostics`, `Issues`, and `DataDictionary`.

Conservative interpretation rules to carry into the manual:

- Tentative GC-MS hits are not confirmed chemical identities.
- Database presence is not sample presence.
- Source annotation is not automatically an analysis variable.
- KEGG pathway context is not pathway activity.
- Hazard screening is not risk assessment.
- Natural-product occurrence in another organism is not occurrence in the student sample.
- Bioactivity annotation is not mechanism, efficacy, or clinical relevance.
- Sparse trait matrices need coverage and sensitivity checks before interpretation.

## Upgrade Progress

Files changed during the manual upgrade:

- `training/main.tex`
- `training/preamble.tex`
- `training/README.md`
- `training/assets/README.md`
- `training/chapters/preface.tex`
- `training/chapters/program_architecture.tex`
- `training/chapters/daily_plans.tex`
- `training/chapters/lab_compendium.tex`
- `training/chapters/readings_exercises.tex`
- `training/chapters/project_tracks.tex`
- `training/chapters/completion_standard.tex`
- `training/chapters/appendix_tables.tex`
- `training/worksheets/worksheets.tex`
- `training/rubrics/rubrics.tex`
- `training/program_notes/program_implementation_notes.tex`
- `training/main.pdf`
- `manual_upgrade_report.md`

Major content additions:

- Expanded front matter with audience, prerequisites, student commitments, student products, mentor review targets, calendar-day expectations, and a 10-week arc.
- Expanded program architecture with full-day weekday structure, weekend asynchronous expectations, work block definitions, product ladder, decision gates, student briefing, preparation checklist, and artifact accumulation logic.
- Replaced the daily schedule with a full 70-day spine covering all 10 weeks, including required weekend work. Each day includes purpose, objectives, prerequisites, materials, time budget, concept work, guided lab or reading discussion, independent project work, evidence/writing work, log prompt, deliverable, quality gate, common problems, and advanced challenge.
- Rebuilt the lab compendium into a 16-lab production bank covering question triage, reproducible project build, source audit, offline GC-MS workflow, research enrichment smoke test, Evidence Court, trait matrix design, figure clinic, literature bridge, handoff audit, identifier ambiguity, wrong-column input repair, source coverage missingness, matrix sparsity, reviewer response, and final package recovery.
- Rebuilt worksheets into 20 student artifacts with purpose, timing, instructions, required fields, common mistakes, and completion standards.
- Rebuilt rubrics into 10 weighted assessment categories covering weekly deliverables, research log, code/reproducibility, chemical evidence, trait matrix/grouping dictionary, figures, literature, report, presentation, and final package/handoff.
- Expanded the implementation guide with pre-program setup, daily and weekly review scripts, decision-gate questions, expected answers, coaching responses, pacing adjustments, data-mode choices, privacy guidance, final evaluation, and archive guidance.
- Expanded project tracks for plant/agriculture/natural products, environmental exposure/toxicology, food/flavor/aroma, biomedical/public health, methods/data science, and chemistry education/communication.
- Added a Program Completion Standard chapter with final package requirements, non-claim boundaries, final review questions, minimum passing package, and high-impact package criteria.
- Added scientific claim-boundary examples for tentative annotation, database presence, source annotation, pathway context, hazard screening, natural product occurrence, bioactivity, and sparse trait matrices.
- Added a validation naming note clarifying `ValidationSummary` in saved `categorate()` results versus `Summary` from direct `validateCategorateResult()` audits.

Build validation:

```sh
cd training
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

- Final build status: passed.
- Final PDF: `training/main.pdf`
- Final PDF page count: 162 pages.
- Log review after final build found no overfull boxes, unresolved references, unresolved citations, LaTeX errors, or pdfTeX destination warnings.
- Remaining LaTeX notices were font-substitution warnings for bold typewriter text; they do not block rendering.
- Transient LaTeX scratch files were removed with `latexmk -c main.tex` after the successful build.

R and source validation:

```sh
Rscript -e 'devtools::load_all(".", quiet=TRUE); source("training/scripts/run_training_smoke_tests.R")'
```

- Passed.
- Offline project setup and offline core workflow smoke tests completed.
- Service-dependent live tests were intentionally skipped because `UAFR_TRAINING_LIVE` was not set.

```sh
Rscript -e 'devtools::load_all(".", quiet=TRUE); ...'
```

- Verified 16 manual-critical exported functions.
- Verified 17 manual-critical dictionary tables.
- Verified 5 direct validation outputs from `validateCategorateResult()`.
- Verified `exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "reactives")` returns records with the current package data.

Structural validation:

- Confirmed exactly 70 day modules are present in `training/chapters/daily_plans.tex`.
- Confirmed every weekend day in Weeks 1 through 10 has required asynchronous work.
- Confirmed Lab 16, Worksheet 20, Rubric 10, and Program Completion Standard are present.
- Confirmed source text has no standard unfinished-work markers, tool/agent references, local user paths, obvious secrets, or private-key text.

Known remaining issues:

- Human proofreading is still required for tone, institutional preference, and final logo placement.
- The manual intentionally avoids live database validation by default. Live query behavior should be tested only in a controlled session with small query sizes, caching, and throttling.
- The PDF is long by design. If it becomes too large for student distribution, a later pass could split student and mentor packets.

Assumptions made:

- The package source, Rd files, tests, and training scripts are the source of truth for function names and output names.
- Student installation will use the private zip/bundle workflow rather than public GitHub installation.
- Weekdays are full internship days; weekends are lighter but required asynchronous days.
