#!/bin/bash

set -euo pipefail

# Paths are resolved from the reproduction workspace.
unset SLURM_SUBMIT_DIR

if [[ -z "${1-}" ]]; then
    echo "Usage:"
    echo "  bash scripts/preprocessing/adapter_trim_fastqs.sh <dataset_dir>"
    exit 1
fi

DATASET="$1"

DORADO="${DORADO:-dorado}"
command -v "$DORADO" >/dev/null 2>&1 || {
    echo "ERROR: Dorado unavailable: $DORADO" >&2
    exit 1
}

INPUT_DIR="$DATASET/input"
TRIM_DIR="$DATASET/trimmed"
LOG_DIR="$DATASET/logs"

mkdir -p "$TRIM_DIR" "$LOG_DIR"

# ------------------------------------------------------------
# Determine sequencing kit
# ------------------------------------------------------------

if [[ "$DATASET" == *half-flowcell* ]]; then
    KIT="SQK-LSK114"
else
    KIT="SQK-NBD114-24"
fi

echo "=================================================="
echo "DATE=$(date)"
echo "HOST=$(hostname)"
echo "DATASET=$DATASET"
echo "KIT=$KIT"
echo "INPUT_DIR=$INPUT_DIR"
echo "TRIM_DIR=$TRIM_DIR"
echo "DORADO_VERSION=$("$DORADO" --version 2>&1 | tail -n1)"
echo "=================================================="

# ------------------------------------------------------------
# Copy sequencing summary
# ------------------------------------------------------------

cp -av \
    "$INPUT_DIR/sequencing_summary.tsv" \
    "$TRIM_DIR/sequencing_summary.tsv"

# ------------------------------------------------------------
# Trim all FASTQs
# ------------------------------------------------------------

SUMMARY_TSV="$LOG_DIR/trim_summary.tsv"

echo -e "file\tinput_bytes\toutput_bytes" > "$SUMMARY_TSV"

find "$INPUT_DIR" \
    -maxdepth 1 \
    -name "*.fastq" \
    | sort |
while read -r FASTQ
do
    BASENAME=$(basename "$FASTQ")

    OUTPUT="$TRIM_DIR/$BASENAME"

    echo
    echo "=================================================="
    echo "INPUT  : $FASTQ"
    echo "OUTPUT : $OUTPUT"
    echo "=================================================="

    INPUT_SIZE=$(stat -c%s "$FASTQ")

    "$DORADO" trim \
        --threads 16 \
        --sequencing-kit "$KIT" \
        --no-trim-primers \
        --emit-fastq \
        "$FASTQ" \
        > "$OUTPUT" \
        2> "$LOG_DIR/${BASENAME}.trim.log"

    OUTPUT_SIZE=$(stat -c%s "$OUTPUT")

    printf "%s\t%s\t%s\n" \
        "$BASENAME" \
        "$INPUT_SIZE" \
        "$OUTPUT_SIZE" \
        >> "$SUMMARY_TSV"

    echo
    echo "INPUT_BYTES=$INPUT_SIZE"
    echo "OUTPUT_BYTES=$OUTPUT_SIZE"
done

echo
echo "========== COMPLETE =========="
date
