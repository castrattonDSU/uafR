# DSU dsDNA Core Undergraduate Research Training Program

This directory contains a full-day, 10-week undergraduate research training
program for using uafR and complementary reproducible research practices in
chemical and mass spectrometry projects. It includes student chapters, daily
plans, lab protocols, reading lists, real-world exercises, worksheets, rubrics,
and program implementation notes.

## Build the manual

From the repository root:

```sh
cd training
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The compiled manual is `training/main.pdf`.

If `latexmk` is unavailable, run `pdflatex main.tex`, `bibtex main`,
`pdflatex main.tex`, and `pdflatex main.tex`.

## Build the student install bundle

From the repository root:

```sh
Rscript tools/build_student_bundle.R
```

The generated zip in `student_bundle/` contains the local uafR package archive,
`START_HERE.md`, the preflight check, installer, update script, verification
script, student acceptance test, training manual PDF, classroom support scripts,
and data templates. Students install from that bundle rather than from an
online source repository.

## Logo

Place the dsDNA Core logo at:

```text
training/assets/dsdna_core_logo.png
```

The manual compiles without the logo and replaces the placeholder automatically
when the file is present.

## Directory guide

- `chapters/`: student-facing weekly curriculum and appendices.
- `scripts/`: R scripts used in the weekly labs.
- `worksheets/`: printable planning and analysis worksheets.
- `rubrics/`: assessment rubrics for students and program review.
- `program_notes/`: program implementation notes, facilitation notes, and pacing guidance.
- `assets/`: logo and future visual assets.
- `data/`: optional local training datasets.
- `bibliography/`: references used by the manual.
