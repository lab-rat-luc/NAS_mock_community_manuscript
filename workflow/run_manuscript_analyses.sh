#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -lt 1 || "$#" -gt 2 ]]; then
    echo "Usage: $0 WORKSPACE [--plan|--execute]" >&2
    exit 1
fi

WORKSPACE="$(realpath -m "$1")"
MODE="${2:---plan}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ANALYSES="$REPO_DIR/scripts/analyses"
# Python and R are supplied by the active paper1-analysis environment.

if [[ "$MODE" != "--plan" && "$MODE" != "--execute" ]]; then
    echo "ERROR: Expected --plan or --execute." >&2
    exit 1
fi

SCRIPTS=(
    run_performance/make_run_performance_tables.py
    run_performance/extract_minknow_mux_scan_pores.py
    run_performance/plot_run_performance_publication.R
    run_performance/plot_minknow_mux_scan_pores.R
    accept_reject_length/make_accept_reject_length_tables.py
    accept_reject_length/plot_accept_reject_length.R
    coverage_over_time/make_coverage_over_time.py
    coverage_over_time/plot_coverage_over_time.R
    enrichment/analyse_enrichment.R
    ra_full_vs_half/plot_ra_full_vs_half.R
    ra_deviation/analyse_ra_deviation_linear_regression.R
    species_mito_mapping_table/make_species_mapping_table.py
)

echo "=== Paper 1 manuscript analyses ==="
echo "Workspace: $WORKSPACE"
echo "Execution: Python and R from the active environment"
echo

for index in "${!SCRIPTS[@]}"; do
    relative="${SCRIPTS[$index]}"
    [[ -f "$ANALYSES/$relative" ]] || {
        echo "ERROR: Missing script: $relative" >&2
        exit 1
    }
    printf '%02d. %s\n' "$((index + 1))" "$relative"
done

if [[ "$MODE" == "--plan" ]]; then
    echo
    echo "PLAN ONLY: No files written or analyses executed."
    exit 0
fi

if [[ "$WORKSPACE" == "/" ]]; then
    echo "ERROR: Root directory cannot be an output workspace." >&2
    exit 1
fi

if [[ "$WORKSPACE" == "$REPO_DIR" ||
      "$WORKSPACE" == "$REPO_DIR/"* ||
      "$REPO_DIR" == "$WORKSPACE/"* ]]; then
    echo "ERROR: Workspace overlaps the public repository." >&2
    exit 1
fi

[[ -d "$WORKSPACE" ]] || {
    echo "ERROR: Workspace missing." >&2
    exit 1
}

for relative in \
    canonical_data/inputs/reference_manifest.tsv \
    canonical_data/outputs/canonical_mito_mapping_summary.tsv \
    canonical_data/outputs/canonical_mapped_alignments.tsv.gz \
    canonical_data/outputs/canonical_read_metadata.tsv.gz
do
    [[ -f "$WORKSPACE/$relative" ]] || {
        echo "ERROR: Missing canonical input: $relative" >&2
        exit 1
    }
done

command -v python3 >/dev/null
command -v Rscript >/dev/null
python3 -c 'import numpy; import pandas'
Rscript --version

# An existing run directory requires inspection, not blind rerunning.
LOG_DIR="$WORKSPACE/workflow_logs/manuscript_analyses"
if [[ -e "$LOG_DIR" || -L "$LOG_DIR" ]]; then
    echo "ERROR: Analysis run logs already exist; inspect before rerunning." >&2
    exit 1
fi
mkdir -p "$LOG_DIR"

cd "$WORKSPACE"

for index in "${!SCRIPTS[@]}"; do
    relative="${SCRIPTS[$index]}"
    script="$ANALYSES/$relative"
    name="$(basename "$relative")"
    printf -v number '%02d' "$((index + 1))"
    log="$LOG_DIR/${number}_${name}.log"

    echo "RUNNING: $relative"

    case "$name" in
        *.py)
            python3 "$script" --project-dir "$WORKSPACE" \
                > "$log" 2>&1
            ;;
        *.R)
            Rscript "$script" "$WORKSPACE" > "$log" 2>&1
            ;;
    esac

    echo "COMPLETED: $relative"
done

echo "PASS: All 12 manuscript-analysis scripts completed."
echo "Numerical and figure comparisons are separate internal checks."
