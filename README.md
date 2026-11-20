# Insect mock-community nanopore sequencing analysis

Code supporting the manuscript:

**Evaluating nanopore adaptive sampling (NAS) for recovering the mitochondrial genomes of low-abundance pest and beneficial species in mixed insect samples**

## Overview

This repository contains the workflow used to compare NAS with metagenomic sequencing in a defined insect mock community including sixteen species.

The workflow performs adapter trimming, read filtering, read mapping and statistical analyses to reproduce the manuscript results.

## Data

- **Raw sequencing data:** Seven Oxford Nanopore FASTQ datasets are available through [NCBI BioProject PRJNA1562222](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1562222).
- **Supporting data:** Six sequencing summaries and three MinKNOW reports are available through [Zenodo](https://doi.org/10.5281/zenodo.24239885).

Raw SUP sequencing summaries are for the whole-run metrics table. Analytical summaries support read classification and other downstream analyses. Both are required.

## Installation

Install the unified environment using environment.yml:

    micromamba create -f environment.yml
    micromamba activate paper1-analysis

Dorado v1.4.0 must be installed separately.

See docs/INSTALLATION.md for details.

## Workflow

1. Prepare inputs.
2. Trim, filter and map reads.
3. Generate canonical dataset tables.
4. Statistic, table and figure generation.

Repository contents:

- workflow/: execution scripts.
- scripts/: preprocessing, mapping and analyses.
- metadata/: input mappings and reference metadata.
- references/: mitochondrial genome reference sequences for each of the sixteen insect species included in the mock.
- docs/: installation and execution instructions.

See docs/WORKFLOW.md for complete instructions.

## Citation and licence

Citation information is provided in CITATION.cff.

Released under the MIT License.
Copyright (c) 2026 Agriculture Victoria.

Contact: luc.shepherd@agriculture.vic.gov.au
