#!/usr/bin/env bash
set -euo pipefail

# Usage:
# bash run_canonical.sh WORKSPACE SUPPORTING_ROOT [--plan|--execute]

if [[ "$#" -lt 2 || "$#" -gt 3 ]]; then
    echo "Usage: $0 WORKSPACE SUPPORTING_ROOT [--plan|--execute]" >&2
    exit 1
fi

WORKSPACE="$(realpath -m "$1")"
SUPPORTING_ROOT="$(realpath -m "$2")"
MODE="${3:---plan}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ "$MODE" != "--plan" && "$MODE" != "--execute" ]]; then
    echo "ERROR: Expected --plan or --execute." >&2
    exit 1
fi

GENERATORS=(
    make_canonical_mito_mapping_summary.py
    make_canonical_barcode_mito_mapping_summary.py
    make_canonical_read_metadata.py
    make_canonical_run_summary_reanalysis.py
)

for name in "${GENERATORS[@]}"; do
    [[ -f "$REPO_DIR/scripts/canonical/$name" ]] || {
        echo "ERROR: Missing generator: $name" >&2
        exit 1
    }
done

[[ -f "$REPO_DIR/workflow/prepare_canonical_inputs.py" ]] || {
    echo "ERROR: Canonical-input helper missing." >&2
    exit 1
}

echo "=== Paper 1 canonical data generation ==="
echo "Repository: $REPO_DIR"
echo "Workspace: $WORKSPACE"
echo "Supporting-data root: $SUPPORTING_ROOT"
echo
echo "1. Prepare five manifests, reference metadata and three reports."
echo "2. Generate main mitochondrial alignments and summary."
echo "3. Generate barcode mitochondrial alignments and summary."
echo "4. Generate canonical read metadata."
echo "5. Generate canonical run summary."
echo

if [[ "$MODE" == "--plan" ]]; then
    echo "PLAN ONLY: No files written or generators executed."
    exit 0
fi

if [[ "$WORKSPACE" == "/" ]]; then
    echo "ERROR: Root directory cannot be an output workspace." >&2
    exit 1
fi

for protected in "$REPO_DIR" "$SUPPORTING_ROOT"; do
    if [[ "$WORKSPACE" == "$protected" ||
          "$WORKSPACE" == "$protected/"* ||
          "$protected" == "$WORKSPACE/"* ]]; then
        echo "ERROR: Workspace overlaps repository or supporting data." >&2
        exit 1
    fi
done

[[ -d "$WORKSPACE" && -d "$SUPPORTING_ROOT" ]] || {
    echo "ERROR: Workspace or supporting-data root missing." >&2
    exit 1
}

# This is a first-run driver, not an overwrite/restart procedure.
if [[ -e "$WORKSPACE/canonical_data" ||
      -L "$WORKSPACE/canonical_data" ]]; then
    echo "ERROR: canonical_data already exists; inspect before proceeding." >&2
    exit 1
fi

command -v python3 >/dev/null
command -v samtools >/dev/null
command -v seqkit >/dev/null
python3 -c 'import numpy; import pandas'

python3 "$REPO_DIR/workflow/prepare_canonical_inputs.py" \
    --workspace "$WORKSPACE" \
    --supporting-root "$SUPPORTING_ROOT" \
    --execute

mkdir -p "$WORKSPACE/canonical_data/logs"

for name in "${GENERATORS[@]}"; do
    echo "RUNNING: $name"

    python3 "$REPO_DIR/scripts/canonical/$name" \
        --project-dir "$WORKSPACE" \
        > "$WORKSPACE/canonical_data/logs/${name%.py}.console.log" 2>&1

    echo "COMPLETED: $name"
done

echo "PASS: All four canonical generators completed."
echo "Manuscript analyses are a separate stage."
