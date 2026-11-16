# Reproducing the manuscript analyses

## Status

The workflow is assembled. Full analytical reproduction is pending.
Do not interpret syntax checks or plan-mode success as result validation.

## 1. Install software

Follow INSTALLATION.md, then activate the unified environment:

    micromamba activate paper1-analysis

Configure the separately installed Dorado 1.4.0 executable:

    export DORADO="/path/to/dorado-1.4.0/bin/dorado"

## 2. Download inputs and configure paths

Obtain the seven FASTQs from SRA and the supporting dataset from Zenodo.
Accessions and download links will be added to the repository README.

Set absolute paths appropriate to your installation:

    export REPO_DIR="/path/to/downloaded/code"
    export SRA_FASTQ_DIR="/path/to/downloaded/fastqs"
    export SUPPORTING_ROOT="/path/to/extracted/zenodo/dataset"
    export SUPPORTING_DATA_DIR="$SUPPORTING_ROOT/sequencing_metadata"
    export OUTPUT_DIR="/path/to/separate/reproduction_workspace"

The workspace must be outside the repository and input directories.
Use a new workspace. Do not point it at existing analytical results.

The supporting-data root must contain sequencing_metadata/ and the
three dataset directories holding the MinKNOW report JSON files.

## 3. Preview input import

    python3 "$REPO_DIR/workflow/prepare_inputs.py" \
      --sra-dir "$SRA_FASTQ_DIR" \
      --supporting-dir "$SUPPORTING_DATA_DIR" \
      --output-dir "$OUTPUT_DIR"

This checks mappings and file availability without decompressing data.

## 4. Import inputs

Run substantial processing on an appropriately resourced workstation
or compute allocation, not an HPC login node.

    python3 "$REPO_DIR/workflow/prepare_inputs.py" \
      --sra-dir "$SRA_FASTQ_DIR" \
      --supporting-dir "$SUPPORTING_DATA_DIR" \
      --output-dir "$OUTPUT_DIR" \
      --execute

This restores the analytical filenames using metadata/input_mapping.tsv.
Source files remain unchanged.

## 5. Preprocess and map

Preview:

    bash "$REPO_DIR/workflow/run_preprocessing.sh" --plan

Execute:

    bash "$REPO_DIR/workflow/run_preprocessing.sh" --execute

This merges untrimmed whole-flowcell technical replicates, builds the
reference index, trims individual and merged FASTQs separately,
length-filters reads, extracts subsets and maps mitochondrial reads.

Do not run the merged-only filtering or mapping scripts afterwards.
See PREPROCESSING.md for analytical settings.

## 6. Generate canonical data

Preview:

    bash "$REPO_DIR/workflow/run_canonical.sh" \
      "$OUTPUT_DIR" "$SUPPORTING_ROOT" --plan

Execute:

    bash "$REPO_DIR/workflow/run_canonical.sh" \
      "$OUTPUT_DIR" "$SUPPORTING_ROOT" --execute

This prepares canonical manifests and MinKNOW report inputs, then
generates main and barcode mapping summaries, read metadata and
the canonical run summary.

## 7. Generate manuscript analyses

Preview:

    bash "$REPO_DIR/workflow/run_manuscript_analyses.sh" \
      "$OUTPUT_DIR" --plan

Execute:

    bash "$REPO_DIR/workflow/run_manuscript_analyses.sh" \
      "$OUTPUT_DIR" --execute

The driver runs Python generators and R analyses in dependency order.
The whole-versus-half-flowcell analysis supplies the species-colour
map used by the expected-versus-observed regression.

## Outputs and failures

Generated data, figures and logs remain in the reproduction workspace.
Some analysis outputs are under canonical_data/analyses/; others are
under analysis/. These locations preserve the analytical structure.

Drivers stop on errors. They are not automatic restart managers.
Inspect the relevant logs before rerunning; do not delete partial
results or disable checks without understanding the failure.

Storage and runtime requirements will be documented after testing.

## Supporting-summary scope

The current implementation imports both sequencing-summary types.
Internal testing will determine whether one type per physical run
can support all reported analyses without changing the results.
The final deposition and instructions will reflect that decision.

## Validation

The authors will compare reconstructed inputs, canonical tables,
statistics and figures against the authoritative analysis before release.
Readers do not need internal audit records to execute this workflow.
