#!/usr/bin/env bash
set -euo pipefail

# Run after prepare_inputs.py has reconstructed the input files.
# Default: print the execution plan without writing files.

MODE="${1:---plan}"
if [[ "$#" -gt 1 ||
      ( "$MODE" != "--plan" && "$MODE" != "--execute" ) ]]; then
    echo "Usage: $0 [--plan|--execute]" >&2
    exit 1
fi

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$REPO_DIR/workflow/config.sh"

OUTPUT_DIR="$(realpath -m "$OUTPUT_DIR")"
REFERENCE_FASTA="$(realpath -m "$REFERENCE_FASTA")"

DATASETS=(
    full-meta-dataset
    full-nas-dataset
    half-flowcell-dataset
)

# Ignore scheduler variables inherited from an HPC session.
unset SLURM_SUBMIT_DIR SLURM_ARRAY_TASK_ID PAPER1_LOCAL_TASK_ID

echo "=== Paper 1 preprocessing and mapping ==="
echo "Repository: $REPO_DIR"
echo "Workspace: $OUTPUT_DIR"
echo "Reference: $REFERENCE_FASTA"
echo
echo "1. Copy reference FASTA into workspace/refs/"
echo "2. Merge untrimmed whole-flowcell barcode FASTQs"
echo "3. Build mitochondrial reference index"
echo "4. Adapter-trim individual and merged FASTQs"
echo "5. Length-filter and extract read subsets"
echo "6. Map retained subsets and index BAMs"
echo
echo "Merged-only filtering and mapping scripts will not be run."

if [[ "$MODE" == "--plan" ]]; then
    echo "PLAN ONLY: No files written or analytical commands executed."
    exit 0
fi

# Restrict execution to a separate reproduction workspace.
SRA_REAL="$(realpath -e "$SRA_FASTQ_DIR")"
SUPPORT_REAL="$(realpath -e "$SUPPORTING_DATA_DIR")"

for input_dir in "$SRA_REAL" "$SUPPORT_REAL"; do
    if [[ "$OUTPUT_DIR" == "$input_dir" ||
          "$OUTPUT_DIR" == "$input_dir/"* ||
          "$input_dir" == "$OUTPUT_DIR/"* ]]; then
        echo "ERROR: Workspace overlaps an input directory." >&2
        exit 1
    fi
done

if [[ "$OUTPUT_DIR" == "/" ]]; then
    echo "ERROR: Root directory cannot be an output workspace." >&2
    exit 1
fi

if [[ "$OUTPUT_DIR" == "$REPO_DIR" ||
      "$OUTPUT_DIR" == "$REPO_DIR/"* ||
      "$REPO_DIR" == "$OUTPUT_DIR/"* ]]; then
    echo "ERROR: Workspace cannot contain the repository." >&2
    exit 1
fi

[[ -d "$OUTPUT_DIR" && -f "$REFERENCE_FASTA" ]] || {
    echo "ERROR: Imported workspace or reference FASTA missing." >&2
    exit 1
}

for tool in "${DORADO:-dorado}" seqkit pigz minimap2 samtools flock; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: Required executable unavailable: $tool" >&2
        exit 1
    }
done

for dataset in "${DATASETS[@]}"; do
    [[ -f "$OUTPUT_DIR/$dataset/input/sequencing_summary.tsv" ]] || {
        echo "ERROR: Missing sequencing summary for $dataset" >&2
        exit 1
    }

    for directory in trimmed filtered logs; do
        if [[ -e "$OUTPUT_DIR/$dataset/$directory" ]]; then
            echo "ERROR: Existing stage directory: $dataset/$directory" >&2
            echo "Use a fresh imported workspace; do not overwrite results." >&2
            exit 1
        fi
    done
done

for relative in \
    full-meta-dataset/input/barcode01.fastq \
    full-meta-dataset/input/barcode02.fastq \
    full-meta-dataset/input/barcode03.fastq \
    full-nas-dataset/input/barcode09.fastq \
    full-nas-dataset/input/barcode10.fastq \
    full-nas-dataset/input/barcode11.fastq \
    half-flowcell-dataset/input/unclassified.fastq
do
    [[ -f "$OUTPUT_DIR/$relative" ]] || {
        echo "ERROR: Missing imported FASTQ: $relative" >&2
        exit 1
    }
done

if [[ -e "$OUTPUT_DIR/refs" ]]; then
    echo "ERROR: Workspace refs directory already exists." >&2
    exit 1
fi

mkdir "$OUTPUT_DIR/refs"
cp -- "$REFERENCE_FASTA" \
    "$OUTPUT_DIR/refs/mock_mitogenomes_clean.fasta"

bash "$REPO_DIR/workflow/prepare_mapping_inputs.sh" \
    "$OUTPUT_DIR" \
    "$OUTPUT_DIR/refs/mock_mitogenomes_clean.fasta" \
    --execute

cd "$OUTPUT_DIR"

for dataset in "${DATASETS[@]}"; do
    mkdir -p "$dataset/logs"

    bash "$REPO_DIR/scripts/preprocessing/adapter_trim_fastqs.sh" \
        "$dataset" > "$dataset/logs/portable_trimming.log" 2>&1

    bash "$REPO_DIR/scripts/preprocessing/filter_fastqs.sh" \
        "$dataset" > "$dataset/logs/portable_filtering.log" 2>&1

    bash "$REPO_DIR/scripts/mapping/map_filtered_fastqs.sh" \
        "$dataset" > "$dataset/logs/portable_mapping.log" 2>&1
done

echo "PASS: Preprocessing and mapping completed."
echo "Canonical generation and manuscript analyses are separate stages."
