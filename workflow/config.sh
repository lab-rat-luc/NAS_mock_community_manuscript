#!/usr/bin/env bash

# Paper 1 reproducibility workflow configuration.
# Set SRA_FASTQ_DIR and SUPPORTING_DATA_DIR before execution.

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

SRA_FASTQ_DIR="${SRA_FASTQ_DIR:-}"
SUPPORTING_DATA_DIR="${SUPPORTING_DATA_DIR:-}"

OUTPUT_DIR="${OUTPUT_DIR:-}"

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: Set OUTPUT_DIR to a separate reproduction workspace." >&2
    return 1 2>/dev/null || exit 1
fi

REFERENCE_FASTA="$REPO_DIR/references/mock_mitogenomes_clean.fasta"
REFERENCE_INDEX="$OUTPUT_DIR/refs/mock_mitogenomes_clean.mmi"

if [[ -z "$SRA_FASTQ_DIR" ]]; then
    echo "ERROR: SRA_FASTQ_DIR is not configured." >&2
    return 1 2>/dev/null || exit 1
fi

if [[ -z "$SUPPORTING_DATA_DIR" ]]; then
    echo "ERROR: SUPPORTING_DATA_DIR is not configured." >&2
    return 1 2>/dev/null || exit 1
fi
