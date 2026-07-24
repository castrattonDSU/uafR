# DSU dsDNA Core uafR Undergraduate Research Training Program

This directory contains a public, production-oriented 10-week curriculum for
undergraduate researchers using uafR and complementary reproducible-research
practices in chemical and mass-spectrometry projects.

The curriculum has two companion publications:

- [Student manual](uafR_training_manual.pdf): the complete 70-day program,
  labs, readings, exercises, project tracks, worksheets, rubrics, setup
  instructions, and completion standard.
- [Instructor and mentor guide](uafR_instructor_guide.pdf): preparation,
  pacing, review questions, coaching guidance, assessment evidence, and
  curriculum/bundle release procedures.

The student manual does not contain instructor answer guidance or maintainer
procedures.

## Program Scope

The program supports full internship days on weekdays and lighter, required
asynchronous work on weekends. Every calendar day has a purpose, deliverable,
research-log prompt, and quality gate. Local programs must separately document
their approved schedule, compensation or academic-credit terms, holidays,
accommodations, safety procedures, supervision, and data-access rules.

The curriculum teaches students to:

1. Define a feasible chemical research question.
2. Build a reproducible R/uafR project.
3. Interpret tentative GC-MS annotations conservatively.
4. Audit evidence from PubChem, KEGG, LOTUS, FEMA, FDA/SPL, MeSH, and
   literature sources.
5. Build source-backed traits, matrices, figures, tables, and written claims.
6. Use the optional species-first plant phytochemistry workflow after mastering
   the core offline pathway.
7. Deliver a report, presentation, evidence trail, reproducible folder, and
   continuation note.

## Build and Check Both Manuals

Requirements: a current R installation, `latexmk`, BibTeX, and the LaTeX
packages used in `preamble.tex`.

From the repository root:

```sh
Rscript tools/build_training_manuals.R
Rscript tools/check_training_publication.R
```

The stable outputs are:

```text
training/uafR_training_manual.pdf
training/uafR_instructor_guide.pdf
```

To build only one edition:

```sh
Rscript tools/build_training_manuals.R --student-only
Rscript tools/build_training_manuals.R --instructor-only
```

The publication check verifies the 70-day sequence, student/instructor
separation, logo dimensions, required scripts, documented uafR functions,
public-language rules, PDF page counts, and PDF text boundaries.

## Run the Training Smoke Tests

Install or load the current uafR package, then run from the repository root:

```sh
Rscript training/scripts/run_training_smoke_tests.R
```

The default run is offline. Live public-service checks are opt-in and are not
required for CI:

```sh
UAFR_TRAINING_LIVE=1 Rscript training/scripts/run_training_smoke_tests.R
```

## Installation Routes for Students

Public GitHub installation is documented in the package
[README](../README.md#installation) and in the student manual. A versioned
offline bundle is available for classrooms that need a fixed package archive
or cannot access GitHub during setup.

Build the offline bundle from the repository root:

```sh
Rscript tools/build_student_bundle.R
```

The generated zip contains the local uafR source archive, student manual,
training scripts, integrity metadata, preflight check, installer, update
script, verification script, quick reference, and offline acceptance test. The
builder regenerates the student manual before packaging it.

## Directory Guide

- `chapters/`: student-facing curriculum and appendices.
- `scripts/`: verified R support scripts used in training.
- `worksheets/`: printable planning and analysis artifacts.
- `rubrics/`: public assessment criteria.
- `program_notes/`: instructor-facing implementation, review, and release
  guidance; these files are excluded from the student PDF.
- `assets/`: the approved dsDNA Core wordmark and asset documentation.
- `data/`: clearly labeled templates or simulated training inputs.
- `bibliography/`: references used by the student manual.

## Version and Support

- Curriculum version: 1.1
- Revision: July 2026
- Package and curriculum repository:
  <https://github.com/castrattonDSU/uafR>
- Maintainer: Chase Stratton, DSU dsDNA Core
- Package citation: run `citation("uafR")`
- Problems with package code or public instructions:
  <https://github.com/castrattonDSU/uafR/issues>

Public database records are evidence, not proof that a compound occurs in a
project sample. Simulated examples are labeled and must not be reported as
observed or published chemistry.

## Reuse and Institutional Marks

See [PUBLICATION_NOTICE.md](PUBLICATION_NOTICE.md). The repository's MIT
software license does not, by itself, grant permission to reuse DSU or dsDNA
Core names, logos, or wordmarks. No separate open-content license has yet been
declared for the curriculum prose, worksheets, or rubrics.
