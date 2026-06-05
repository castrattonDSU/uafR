# WiCCED Water Quality Manual Upgrade Report

## Source audit

- Created a separate manual folder so the existing 10-week dsDNA Core uafR
  training manual remains intact.
- Online context reviewed: Project WiCCED Salt home/about/research/education
  pages, DSU faculty profiles for Dr. Gulnihal Ozbay and Dr. Chase Stratton,
  DNREC Delaware water-quality monitoring and marsh migration context, and
  recent water-quality machine-learning review literature.
- uafR function references were checked against repository source, NAMESPACE,
  man pages, tests, and existing training scripts before inclusion:
  `spreadOut()`, `mzExacto()`, `exactoThese()`,
  `categorate(detail = "research")`, `validateCategorateResult()`,
  `chemicalTraitMatrix()`, `chemicalTraitEvidence()`,
  `chemicalTraitOntology()`, `chemicalTraitOntologyMatrix()`, and
  `pubchemProfile()`.

## Files added

- `main.tex`, `preamble.tex`, `README.md`, and this report.
- Bibliography: `bibliography/references.bib`.
- Chapters: program overview, WiCCED context, full 8-week spine, uafR workflow
  guidance, machine-learning guidance, lab compendium, completion standard, and
  appendix script inventory.
- Mentor guide: `program_notes/mentor_guide.tex`.
- Worksheets and rubrics: `worksheets/worksheets.tex`,
  `rubrics/rubrics.tex`.
- Teaching data: `data/wicced_teaching_water_quality.csv`,
  `data/chemical_screening_template.csv`,
  `data/wicced_simulated_gcms_hits.csv`,
  `data/wicced_simulated_gcms_queries.csv`,
  `data/wicced_simulated_gcms_compound_reference.csv`, and
  `data/wicced_simulated_gcms_sample_metadata.csv`.
- Training scripts:
  `00_setup_check.R`, `01_project_setup.R`, `02_water_quality_qaqc.R`,
  `03_uafr_chemical_screening_template.R`,
  `04_ml_01_build_feature_table.R`, `04_ml_02_explore_features.R`,
  `04_ml_03_train_models.R`,
  `04_ml_04_diagnostics_and_model_card.R`,
  `05_ml_01_advanced_validation.R`,
  `05_ml_02_permutation_importance_and_sensitivity.R`,
  `05_ml_03_advanced_model_report.R`, and
  `run_wicced_smoke_tests.R`.
- Machine-learning resources: `ml_resources/README.md`,
  `ml_resources/ml_glossary.csv`,
  `ml_resources/ml_workflow_checklist.csv`,
  `ml_resources/ml_troubleshooting.md`,
  `ml_resources/model_card_field_guide.md`, and
  `ml_resources/advanced_ml_extension_guide.md`.
- Local `.gitignore` for LaTeX auxiliary files.

## Major content additions

- Complete 8-week, 56-day active training schedule with full weekday modules
  and lighter required weekend work.
- WiCCED-specific framing around saltwater intrusion, Delaware water quality,
  aquaculture/ecosystem health, marsh migration, and workforce development.
- Explicit co-advising alignment between uafR/cheminformatics/data science and
  water resources/aquaculture/aquatic ecology.
- Production lab bank covering project triage, reproducible setup, QA/QC,
  source audit, simulated uafR-compatible GC-MS evidence, trait matrices, ML
  baselines, leakage checks, figure critique, reviewer response, and final
  recovery.
- Beginner ML ladder covering feature-table construction, feature dictionaries,
  train/test split plans, pre-model exploration, baseline comparison,
  regression/classification metrics, diagnostics, interpretation prompts, and a
  populated model card.
- Advanced ML extension covering time-holdout and leave-one-site-out
  validation, conductivity/context sensitivity, skipped overfit feature sets,
  held-out permutation importance, decision gates, and a written advanced model
  review. This extension is explicitly positioned after beginner mastery.
- Worksheets for project brief, variable dictionary, QA/QC, source evidence,
  uafR traits, ML features/splits/model cards, figures, literature synthesis,
  research logs, and final audits.
- Rubrics for weekly deliverables, logs, reproducibility, water-quality QA/QC,
  uafR evidence traceability, machine learning, figures, reports, and
  presentations.
- Conservative scientific language distinguishing measurements, tentative
  annotations, database records, source-backed traits, ML predictions, and
  supported claims.

## Build status

- Command: `latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex`
  from `training_wicced_water_quality/`.
- Status: passed.
- Output: `training_wicced_water_quality/main.pdf`.
- Notes: final build produced no hard errors, no unresolved citations, and no
  overfull boxes. Remaining LaTeX messages are underfull table-wrapping
  warnings in dense technical tables.

## R/script status

- Command: `Rscript training_wicced_water_quality/scripts/run_wicced_smoke_tests.R`
  from the repository root.
- Status: passed.
- Coverage: setup check, project skeleton, water-quality QA/QC, simulated
  GC-MS hit-table validation/spread/exact-match aggregation/relative abundance
  outputs, four-step beginner ML ladder, regression/classification metrics,
  diagnostic plots, populated model card, three-step advanced ML extension,
  blocked validation, feature-set sensitivity, permutation importance,
  decision-gate review, and expected output verification in a temporary project.
- Live web/database queries are intentionally not part of the smoke test.

## Source checks

- 56 `\DayPlan{}` modules found.
- 8 `\WeekHeader{}` modules found.
- No unresolved marker hits in authored source.
- No authored-source absolute local paths or credential strings found. The only
  credential-pattern false positives were LaTeX macro names such as
  `\detokenize`; generated build logs contain local build paths and are ignored
  by the manual folder `.gitignore`.

## Known remaining issues

- The manual uses the existing dsDNA Core logo path as a fallback. If a
  WiCCED-specific or final DSU/dsDNA Core logo should appear on the title page,
  place it in the expected asset path and adjust `\CoreLogo` if needed.
- Teaching water-quality and GC-MS data are synthetic and support workflow
  training only; real WiCCED data must be approved and documented before
  environmental claims are made.
- Live `categorate()`/`pubchemProfile()` examples are described as optional
  templates because service availability and package dependencies vary by
  machine.

## Assumptions made

- The program should remain separate from the existing 10-week manual.
- Weekends should be active but lighter than weekdays.
- uafR is included as a chemical-screening and source-evidence layer, not as a
  replacement for water-quality field/lab measurement.
- Machine learning should emphasize baseline models, leakage control,
  validation, model cards, and cautious interpretation before complex
  algorithms.
