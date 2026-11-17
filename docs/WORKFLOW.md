# Reproducing the manuscript analyses

The workflow starts from the seven FASTQs deposited on SRA and the sequencing summaries. Run the four stages in order: input preparation, preprocessing and mapping, canonical data table generation, and analyses.

Use a Linux-compatible workstation or compute allocation. Keep downloaded inputs, the code repository and the reproduction workspace in separate directories.

## 1. Install software

Follow [INSTALLATION.md](INSTALLATION.md), then activate the environment:

```bash
micromamba activate paper1-analysis
```

Dorado v1.4.0 must be installed separately.

## 2. Download and organise inputs

See the repository [README](../README.md) for data identifiers and availability.

Download the seven compressed FASTQs into one directory:

```text
fastqs/
├── mock_meta_rep1.fastq.gz
├── mock_meta_rep2.fastq.gz
├── mock_meta_rep3.fastq.gz
├── mock_nas_rep1.fastq.gz
├── mock_nas_rep2.fastq.gz
├── mock_nas_rep3.fastq.gz
└── mock_half_flowcell.fastq.gz
```

The Zenodo deposit contains six compressed sequencing summaries, three MinKNOW reports and a README. Download these files into a separate directory named `supporting_data/`.

Keep the summaries compressed. Create three subdirectories and move the JSON reports into the following locations:

```text
supporting_data/
├── README.md
├── full_meta_raw_sup_summary.tsv.gz
├── full_nas_raw_sup_summary.tsv.gz
├── half_flowcell_raw_sup_summary.tsv.gz
├── full_meta_analytical_summary.tsv.gz
├── full_nas_analytical_summary.tsv.gz
├── half_flowcell_analytical_summary.tsv.gz
├── full-meta-dataset/
│   └── report_FBF72598_20260219_1223_4f326e35.json
├── full-nas-dataset/
│   └── report_FBF77024_20260220_1337_8b313231.json
└── half-flowcell-dataset/
    └── report_FBF71559_20260305_1539_2c752f6f.json
```

Both summary types are required: raw SUP summaries support whole-run table metrics, while analytical summaries support read classification and timing analyses.

## 3. Configure paths

In the terminal used to run the workflow, replace these examples with absolute paths:

```bash
export REPO_DIR="/path/to/NAS_mock_community_manuscript"
export SRA_FASTQ_DIR="/path/to/fastqs"
export SUPPORTING_ROOT="/path/to/supporting_data"
export SUPPORTING_DATA_DIR="$SUPPORTING_ROOT"
export OUTPUT_DIR="/path/to/reproduction_workspace"
export REFERENCE_FASTA="$REPO_DIR/references/mock_mitogenomes_clean.fasta"
export DORADO="/path/to/dorado-1.4.0/bin/dorado"
```

`OUTPUT_DIR` must be absent or empty before input import. It must not overlap the repository or either input directory.

Keep these settings available for subsequent stages.

## 4. Prepare inputs

Preview the import:

```bash
python3 "$REPO_DIR/workflow/prepare_inputs.py" \
  --sra-dir "$SRA_FASTQ_DIR" \
  --supporting-dir "$SUPPORTING_DATA_DIR" \
  --output-dir "$OUTPUT_DIR"
```

This checks all 13 input mappings without writing files.

Import the inputs:

```bash
python3 "$REPO_DIR/workflow/prepare_inputs.py" \
  --sra-dir "$SRA_FASTQ_DIR" \
  --supporting-dir "$SUPPORTING_DATA_DIR" \
  --output-dir "$OUTPUT_DIR" \
  --execute
```

The importer decompresses the seven FASTQs and six sequencing summaries into the filenames specified in [input_mapping.tsv](../metadata/input_mapping.tsv). 

## 5. Preprocess and map

Preview:

```bash
bash "$REPO_DIR/workflow/run_preprocessing.sh" --plan
```

Execute:

```bash
bash "$REPO_DIR/workflow/run_preprocessing.sh" --execute
```

This stage:

- Merges untrimmed whole-flowcell technical replicates.
- Builds the mitochondrial reference index.
- Trims individual and merged FASTQs separately.
- Length-filters reads and extracts read subsets.
- Maps retained reads and indexes the resulting BAMs.

Do not run the merged-only filtering or mapping scripts afterwards; the driver already processes the merged inputs.

See [PREPROCESSING.md](PREPROCESSING.md) for analytical settings.

## 6. Generate canonical datasets

Canonical dataset tables are the standardised intermediate files used by downstream analyses.

Preview:

```bash
bash "$REPO_DIR/workflow/run_canonical.sh" \
  "$OUTPUT_DIR" "$SUPPORTING_ROOT" --plan
```

Execute:

```bash
bash "$REPO_DIR/workflow/run_canonical.sh" \
  "$OUTPUT_DIR" "$SUPPORTING_ROOT" --execute
```

This stage prepares manifests and MinKNOW report inputs, then generates:

- Dataset-level and barcode-level mitochondrial mapping summaries.
- Read-level alignment tables.
- Read metadata, including timing and NAS classifications.
- The sequencing run summary.

## 7. Generate manuscript analyses

Preview:

```bash
bash "$REPO_DIR/workflow/run_manuscript_analyses.sh" \
  "$OUTPUT_DIR" --plan
```

Execute:

```bash
bash "$REPO_DIR/workflow/run_manuscript_analyses.sh" \
  "$OUTPUT_DIR" --execute
```

The driver runs the Python generators and R analyses in dependency order. This includes generating the species-colour map before the expected-versus-observed regression analysis.

## Find the outputs

All generated files remain inside `OUTPUT_DIR`.

| Location | Contents |
|---|---|
| `full-meta-dataset/` | Whole-flowcell metagenomic preprocessing and mapping |
| `full-nas-dataset/` | Whole-flowcell NAS preprocessing and mapping |
| `half-flowcell-dataset/` | Half-flowcell preprocessing and mapping |
| `canonical_data/outputs/` | Canonical mapping, read-metadata and run-summary tables |
| `canonical_data/analyses/` | Accepted/rejected read-length, coverage-over-time and species-mapping-table analyses |
| `analysis/` | Run-performance, enrichment and relative-abundance analyses |

Analysis directories contain output tables, figures and logs as applicable. 
