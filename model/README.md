# Ecoli-GEM model and annotation files

This directory contains the Ecoli-GEM model and annotation files.

## Model

The model is available as `.xml`, `.xlsx`, `.txt`, `.yml`, and `.mat`. Note that only the `.yml` version is available on branches other than `main` (e.g., `develop`), to facilitate tracking of model changes.

## FROG report

A FROG report generated using https://runfrog.de/ on the SBML model file is available in the `FROG` folder. This folder contains `tsv` files detailing the results of the FROG analysis as well as a `json` file describing the metadata. The `frog.tsv` is skipped (as in, not versioned) because its content is already duplicated in the other files.
