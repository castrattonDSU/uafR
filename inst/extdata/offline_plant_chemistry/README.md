# Offline Plant Chemistry Example

This directory contains a small simulated teaching fixture for offline uafR
workflow checks. It is not a real phytochemistry database extract and must not
be interpreted as evidence that these compounds occur in the listed plants.

The fixture is designed to exercise production workflows without live web
queries:

- direct species database-style evidence
- genus/family fallback context
- PubMed/PubTator candidate-only rows
- resolved and unresolved compounds
- comparable and non-comparable chemistry
- plant-pair compound Tanimoto rows

Regenerate the example bundle from the repository root with:

```sh
Rscript tools/build_offline_plant_chemistry_example.R
```
