# Communications Biology public-code release skeleton

This directory is the safe public release for the frozen manuscript package. It contains no raw GWAS summary statistics, individual-level data, LD reference matrices, credentials or absolute local paths.

Release identifiers: local release label `v1.0-submission-rc1`; public code repository [GitHub](https://github.com/seefreewind/immune_cm_psych_multimorbidity); archived code and derived-data package [Zenodo DOI 10.5281/zenodo.22689138](https://doi.org/10.5281/zenodo.22689138). The Zenodo record is version 1 and contains the deposited release ZIP.

## Intended structure

- `figure_scripts/`: copies of the project figure-generation entry points.
- `figure_source_data/`: frozen TSV inputs used by the displayed figures.
- `frozen_tables/`: canonical derived tables for the manuscript summaries.
- `manifests/`: input/output and exclusion manifests.
- `qc/`: figure and manuscript number/source audits.
- `environment/`: software versions and platform notes.
- `validation/`: commands for figure and number audits.

## Reproduction order

1. Obtain the public source GWAS files from the accession and provider links in the manuscript Data Availability statement.
2. Obtain any controlled-access resources directly from the original providers under their terms.
3. Use the frozen derived tables and figure source-data TSV files deposited with the manuscript.
4. Run the descriptive figure scripts and the validation commands after setting relative project paths.

The statistical analyses are frozen in the accompanying project and are not re-run by this release skeleton.
