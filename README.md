# Insect mock-community nanopore sequencing analysis

Code supporting the manuscript:

**Evaluating nanopore adaptive sampling (NAS) for recovering the mitochondrial genomes of low-abundance pest and beneficial species in mixed insect samples**

## Overview

This repository contains the workflow used to compare NAS with metagenomic 
sequencing in a defined 16-species insect mock community.

The workflow performs adapter trimming, read filtering, mitochondrial
mapping and statistical analyses to reproduce the manuscript results.

## Data

- NCBI SRA: Seven FASTQs (submission SUB16635091; BioSample SAMN65318136).
- Zenodo: Six sequencing summaries and three MinKNOW reports.
  Reserved DOI: 10.5281/zenodo.24239885.

Raw SUP summaries support complete-run metrics. Analytical summaries
support read classification and downstream analyses. Both are required.

## Installation

Install the unified environment using environment.yml:

    micromamba create -f environment.yml
    micromamba activate paper1-analysis

Dorado v1.4.0 must be installed separately.

See docs/INSTALLATION.md for details.

## Workflow

1. Prepare inputs.
2. Trim, filter and map reads.
3. Generate canonical datasets.
4. Reproduce manuscript statistics, tables and figures.

Repository contents:

- workflow/: execution scripts.
- scripts/: preprocessing, mapping and analyses.
- metadata/: input mappings and reference metadata.
- references/: mitochondrial reference sequences.
- docs/: installation and execution instructions.

See docs/WORKFLOW.md for complete instructions.

## Citation and licence

Citation information is provided in CITATION.cff.

Released under the MIT License.
Copyright (c) 2026 Agriculture Victoria.

Contact: luc.shepherd@agriculture.vic.gov.au
