# Insect mock-community nanopore sequencing analysis

Code supporting the manuscript:

**Evaluating nanopore adaptive sampling (NAS) for recovering the mitochondrial genomes of low-abundance pest and beneficial species in mixed insect samples**

## Scope

The workflow starts from seven deposited FASTQs and supporting
sequencing metadata. It performs adapter trimming, length filtering,
read-subset extraction, mitochondrial mapping, canonical data
generation and manuscript analyses.

Offline rebasecalling is upstream provenance, not an execution
requirement. Deposited reads were generated using Dorado v1.4.0
with model dna_r10.4.1_e8.2_400bps_sup@v5.2.0.
POD5 files are not required by this workflow.

## Data sources

- NCBI SRA: seven sequencing FASTQs, submission SUB16635091.
  BioSample accession: SAMN65318136. BioProject and SRA run
  accessions will be added when assigned.
- Zenodo: six compressed sequencing summaries and three MinKNOW
  reports. Reserved DOI: https://doi.org/10.5281/zenodo.24239885.
  The DOI becomes active when the dataset is published.
- GitHub: https://github.com/lab-rat-luc/NAS_mock_community_manuscript
- This repository includes the mitochondrial reference FASTA,
  analytical reference metadata, corrected community composition
  and analytical code.

The raw SUP and analytical sequencing summaries serve different
workflow roles. Both are required by the current implementation.
Their filenames are mapped in metadata/input_mapping.tsv.

## Repository contents

- workflow/: input preparation and stage drivers.
- scripts/preprocessing/: adapter trimming and read filtering.
- scripts/mapping/: mitochondrial mapping.
- scripts/canonical/: canonical data generators.
- scripts/analyses/: statistical analyses and plotting.
- references/: exact mitochondrial reference FASTA.
- metadata/: input mappings, reference metadata and manifest templates.
- docs/: detailed workflow instructions.

The analytical reference manifest preserves the reference identifiers
used during analysis. The corrected community composition table
describes the specimens, including Apis mellifera identified to
species level. A reference identifier containing a subspecies name
does not establish the specimen's subspecies identity.

## Software and installation

Requires a Linux-compatible command-line environment with Bash,
Python, R, Dorado, SeqKit, pigz, minimap2, samtools, flock and standard
Linux utilities.

Use the unified environment specified in environment.yml.
Installation instructions are in [docs/INSTALLATION.md](docs/INSTALLATION.md).
Activate paper1-analysis before running the workflow.
Dorado 1.4.0 is configured separately as described in that guide.

## Input organisation

Download the seven FASTQs into one directory:

- mock_meta_rep1.fastq.gz
- mock_meta_rep2.fastq.gz
- mock_meta_rep3.fastq.gz
- mock_nas_rep1.fastq.gz
- mock_nas_rep2.fastq.gz
- mock_nas_rep3.fastq.gz
- mock_half_flowcell.fastq.gz

Download the ten Zenodo supporting-data files into a directory
named `supporting_data/`.

Keep the six compressed sequencing summaries directly inside
`supporting_data/`.

Create three subdirectories and move the corresponding MinKNOW
reports into them:

- `supporting_data/full-meta-dataset/`
- `supporting_data/full-nas-dataset/`
- `supporting_data/half-flowcell-dataset/`

The importer takes `supporting_data/` as its `--supporting-dir`.
The canonical driver takes the same directory as its
supporting-data root.

## Execution order

1. workflow/prepare_inputs.py
   Restore the analytical filenames from the downloaded gzip files.
   Defaults to a dry run; --execute performs decompression.

2. workflow/run_preprocessing.sh
   Merge untrimmed technical replicates, prepare the reference index,
   trim, filter and map reads. Defaults to plan mode.

3. workflow/run_canonical.sh
   Prepare canonical manifests and generate canonical datasets.
   Defaults to plan mode.

4. workflow/run_manuscript_analyses.sh
   Generate the computational manuscript tables, statistics and figures.
   Defaults to plan mode.

Detailed preprocessing instructions are in
[docs/PREPROCESSING.md](docs/PREPROCESSING.md).

Use a separate reproduction workspace and an appropriately resourced
workstation or compute allocation. Do not run substantial processing
on an HPC login node.

Drivers are first-run workflows, not automatic restart managers.
Inspect failures before rerunning. Intermediate FASTQs and BAMs can
require substantial storage; measured requirements will be documented
after internal testing.

## Analytical definitions

- Length filtering: 200-20,000 bp; no additional read-quality threshold.
- Accepted NAS reads: end_reason == signal_positive.
- Rejected NAS reads:
  end_reason == data_service_unblock_mux_change.
- Half-flowcell NAS: channels 1-256.
- Half-flowcell Metagenomic: channels 257-512.
- Mapping: minimap2 map-ont with secondary alignments disabled.
- Retained BAM alignments: MAPQ >=20; exclusion mask 2308.

Individual and merged whole-flowcell FASTQs are trimmed separately.
The general filtering and mapping scripts also process merged inputs;
the merged-only alternatives are not additional driver stages.

## Citation, licensing and contact

This software is released under the MIT License.
Copyright (c) 2026 Agriculture Victoria.

Citation information is provided in `CITATION.cff`.

Correspondence: luc.shepherd@agriculture.vic.gov.au
