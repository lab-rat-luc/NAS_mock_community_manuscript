# Preprocessing and mitochondrial mapping

## Status
Workflow assembly is in progress. End-to-end reproduction testing is pending.

## Inputs
- Seven compressed FASTQs and six compressed sequencing summaries,
  listed in metadata/input_mapping.tsv.
- The exact mitochondrial reference FASTA in references/.

Run workflow/prepare_inputs.py first to reconstruct the input workspace.

## Configuration
Set SRA_FASTQ_DIR, SUPPORTING_DATA_DIR and OUTPUT_DIR.
Optionally set DORADO to the Dorado executable path.
OUTPUT_DIR must be a separate reproduction workspace.

## Commands
Print the plan:
    bash workflow/run_preprocessing.sh --plan

Execute after installation, configuration and input import:
    bash workflow/run_preprocessing.sh --execute

Use an appropriately resourced workstation or compute allocation,
not an HPC login node.

## Processing order
1. Copy the reference FASTA into workspace/refs/.
2. Merge untrimmed Metagenomic barcodes 01, 02 and 03.
3. Merge untrimmed NAS barcodes 09, 10 and 11.
4. Build the minimap2 reference index.
5. Adapter-trim individual and merged FASTQs separately.
6. Length-filter reads and extract analytical subsets.
7. Map retained subsets and generate indexed BAMs.

The proposed merge order will be checked against the historical
merged datasets during internal reproduction testing.

## Analytical settings
- Dorado kits: SQK-NBD114-24 for whole-flowcell libraries;
  SQK-LSK114 for the half-flowcell library.
- Trimming uses --no-trim-primers.
- Length filter: 200-20,000 bp; no additional quality threshold.
- Accepted: end_reason == signal_positive.
- Rejected: end_reason == data_service_unblock_mux_change.
- Half-flowcell NAS: channels 1-256.
- Half-flowcell Metagenomic: channels 257-512.
- Mapping: map-ont; secondary alignments disabled.
- BAM retention: MAPQ >=20; exclusion mask 2308.

The general filtering and mapping scripts process merged datasets too.
Do not additionally run the merged-only scripts in this pathway.

## Software
Requires Bash, Python, Dorado, SeqKit, pigz, minimap2, samtools,
flock and standard Linux command-line utilities.
See INSTALLATION.md for the unified environment and Dorado setup.

Current settings retain 16 Dorado threads and mapping settings of
12 minimap2 threads, 4 samtools view threads and 4 samtools sort threads.
Mapping pipeline components can run concurrently.

## Outputs and failures
Outputs and stage logs remain within the reproduction workspace.
The driver requires a fresh imported workspace and does not provide
automatic restart or rollback. Inspect failures before rerunning.

Canonical generation and manuscript analyses are separate stages.
