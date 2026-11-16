#!/usr/bin/env bash
set -euo pipefail

# Paper 1: reconstruct merged FASTQs and prepare reference index.
# Usage:
#   bash prepare_mapping_inputs.sh OUTPUT_DIR REFERENCE_FASTA --dry-run
#   bash prepare_mapping_inputs.sh OUTPUT_DIR REFERENCE_FASTA --execute

if [[ "$#" -ne 3 ]]; then
    echo "Usage: $0 OUTPUT_DIR REFERENCE_FASTA --dry-run|--execute" >&2
    exit 1
fi

OUTPUT_DIR="$(realpath -m "$1")"
REFERENCE_FASTA="$(realpath -e "$2")"
MODE="$3"

if [[ "$MODE" != "--dry-run" && "$MODE" != "--execute" ]]; then
    echo "ERROR: Invalid mode." >&2
    exit 1
fi

META_DIR="$OUTPUT_DIR/full-meta-dataset/input"
NAS_DIR="$OUTPUT_DIR/full-nas-dataset/input"

META_MERGED="$META_DIR/full-meta-merged.fastq"
NAS_MERGED="$NAS_DIR/full-nas-merged.fastq"

INDEX="$OUTPUT_DIR/refs/mock_mitogenomes_clean.mmi"

META_INPUTS=(
    "$META_DIR/barcode01.fastq"
    "$META_DIR/barcode02.fastq"
    "$META_DIR/barcode03.fastq"
)

NAS_INPUTS=(
    "$NAS_DIR/barcode09.fastq"
    "$NAS_DIR/barcode10.fastq"
    "$NAS_DIR/barcode11.fastq"
)

for file in "${META_INPUTS[@]}" "${NAS_INPUTS[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Missing input: $file" >&2
        exit 1
    fi
done

if [[ -e "$META_MERGED" || -e "$NAS_MERGED" || -e "$INDEX" ]]; then
    echo "ERROR: One or more outputs already exist." >&2
    exit 1
fi

if [[ "$OUTPUT_DIR" == "/" ]]; then
    echo "ERROR: Root directory cannot be an output workspace." >&2
    exit 1
fi

echo "=== Paper 1 mapping input preparation ==="
echo "Output directory: $OUTPUT_DIR"
echo "Reference FASTA: $REFERENCE_FASTA"

echo
echo "Metagenomic merge:"
printf '%s\n' "${META_INPUTS[@]}"
echo "Destination: $META_MERGED"

echo
echo "NAS merge:"
printf '%s\n' "${NAS_INPUTS[@]}"
echo "Destination: $NAS_MERGED"

echo
echo "Reference index: $INDEX"

if [[ "$MODE" == "--dry-run" ]]; then
    echo
    echo "PASS: Dry run completed. No files created."
    exit 0
fi

command -v minimap2 >/dev/null 2>&1 || {
    echo "ERROR: minimap2 unavailable." >&2
    exit 1
}

mkdir -p "$OUTPUT_DIR/refs"

META_TMP="${META_MERGED}.partial"
NAS_TMP="${NAS_MERGED}.partial"
INDEX_TMP="${INDEX}.partial"

for file in "$META_TMP" "$NAS_TMP" "$INDEX_TMP"; do
    if [[ -e "$file" ]]; then
        echo "ERROR: Partial output already exists: $file" >&2
        exit 1
    fi
done

cat "${META_INPUTS[@]}" > "$META_TMP"
mv -- "$META_TMP" "$META_MERGED"

cat "${NAS_INPUTS[@]}" > "$NAS_TMP"
mv -- "$NAS_TMP" "$NAS_MERGED"

minimap2 -d "$INDEX_TMP" "$REFERENCE_FASTA"
mv -- "$INDEX_TMP" "$INDEX"

echo
echo "PASS: Merged FASTQs and reference index created."
