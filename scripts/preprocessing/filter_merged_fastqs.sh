#!/bin/bash

set -euo pipefail

# Paths are resolved from the reproduction workspace.
unset SLURM_SUBMIT_DIR
mkdir -p logs

if [[ -z "${1-}" ]]; then
  echo "Usage: bash scripts/preprocessing/filter_merged_fastqs.sh <dataset_dir>"
  echo "Example:"
  echo "  bash scripts/preprocessing/filter_merged_fastqs.sh full-meta-dataset"
  echo "  bash scripts/preprocessing/filter_merged_fastqs.sh full-nas-dataset"
  exit 1
fi

DATASET="$1"


command -v seqkit >/dev/null 2>&1 || { echo "ERROR: seqkit not found"; exit 1; }
command -v pigz >/dev/null 2>&1 || { echo "ERROR: pigz not found"; exit 1; }

RUN_DIR="$DATASET/trimmed"
SUMMARY="$RUN_DIR/sequencing_summary.tsv"

if [[ ! -f "$SUMMARY" ]]; then
  echo "ERROR: missing $SUMMARY"
  exit 1
fi

SUMMARY_TABLE="$DATASET/filtering_summary_merged.tsv"
LOCK_FILE="$SUMMARY_TABLE.lock"
FASTQ_LIST="$DATASET/logs/input_merged_fastqs.txt"

mkdir -p "$DATASET/logs"

echo -e "dataset\tfile\tstep\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > "$SUMMARY_TABLE"

echo "===== MERGED FILTERING START ====="
echo "DATE=$(date)"
echo "HOST=$(hostname)"
echo "DATASET=$DATASET"
echo "SUMMARY=$SUMMARY"
echo "RUN_DIR=$RUN_DIR"
echo "SEQKIT=$(seqkit version)"
echo "PIGZ=$(pigz --version 2>&1 | head -n 1)"
echo "LENGTH_FILTER=200-20000"
echo "ACCEPTED=signal_positive"
echo "REJECTED=data_service_unblock_mux_change"

# ---------- dataset type ----------
if [[ "$DATASET" == *nas* ]]; then
  TYPE="NAS"
elif [[ "$DATASET" == *meta* ]]; then
  TYPE="META"
else
  echo "ERROR: this merged-only script is intended for full-meta-dataset or full-nas-dataset"
  exit 1
fi

echo "TYPE=$TYPE"

# ---------- merged FASTQ manifest only ----------
: > "$FASTQ_LIST"

if [[ "$TYPE" == "NAS" ]]; then
  mapfile -t MERGED_FASTQS < <(find "$RUN_DIR" -maxdepth 1 -name "full-nas-merged.fastq" | sort)
elif [[ "$TYPE" == "META" ]]; then
  mapfile -t MERGED_FASTQS < <(find "$RUN_DIR" -maxdepth 1 -name "full-meta-merged.fastq" | sort)
fi

FASTQ_COUNT="${#MERGED_FASTQS[@]}"
echo "INPUT_MERGED_FASTQ_COUNT=$FASTQ_COUNT"

if [[ "$FASTQ_COUNT" -ne 1 ]]; then
  echo "ERROR: expected exactly one merged FASTQ for $TYPE, found $FASTQ_COUNT"
  printf '%s\n' "${MERGED_FASTQS[@]}"
  exit 1
fi

printf '%s\n' "${MERGED_FASTQS[@]}" > "$FASTQ_LIST"
cat "$FASTQ_LIST"

# ---------- detect columns ----------
header="$(head -n 1 "$SUMMARY" | tr '\t' '\n')"

rid_col="$(echo "$header" | nl -ba | awk '$2=="read_id"{print $1}')"
er_col="$(echo "$header" | nl -ba | awk '$2=="end_reason"{print $1}')"

if [[ -z "$rid_col" || -z "$er_col" ]]; then
  echo "ERROR: read_id or end_reason column not found in $SUMMARY"
  exit 1
fi

echo "COLUMN_READ_ID=$rid_col"
echo "COLUMN_END_REASON=$er_col"

# ---------- build ID lists ----------
POS="$DATASET/positive.ids"
UNB="$DATASET/unblock.ids"

echo "CMD: build positive.ids"
awk -F'\t' -v r="$rid_col" -v e="$er_col" \
  'NR>1 && $e=="signal_positive"{print $r}' \
  "$SUMMARY" | sort -u > "$POS"

echo "CMD: build unblock.ids"
awk -F'\t' -v r="$rid_col" -v e="$er_col" \
  'NR>1 && $e=="data_service_unblock_mux_change"{print $r}' \
  "$SUMMARY" | sort -u > "$UNB"

echo "ID_COUNTS:"
echo "  POS=$(wc -l < "$POS")"
echo "  UNB=$(wc -l < "$UNB")"

# ---------- export for parallel ----------
export DATASET TYPE POS UNB SUMMARY_TABLE LOCK_FILE

process_fastq() {
  set -euo pipefail

  fq="$1"

  base="$(basename "$fq" .fastq)"
  outdir="$DATASET/filtered/$base"
  mkdir -p "$outdir"

  filtered="$outdir/len200_20000.fastq.gz"

  echo "===== PROCESSING $base ====="
  echo "FASTQ=$fq"
  echo "OUTDIR=$outdir"

  # ---------- length filter ----------
  echo "CMD: seqkit seq -m 200 -M 20000 $fq | pigz -p 2 > $filtered"
  seqkit seq -m 200 -M 20000 "$fq" | pigz -p 2 > "$filtered"

  stats_line="$(seqkit stats -T "$filtered" | tail -n 1)"
  stats="$(echo "$stats_line" | cut -f4-)"

  (
    flock -x 200
    echo -e "$DATASET\t$base\tlen_filtered\t$stats" >> "$SUMMARY_TABLE"
  ) 200>>"$LOCK_FILE"

  if [[ "$TYPE" == "META" ]]; then

    echo "CMD: seqkit grep -f $POS $filtered -o $outdir/positive.fastq.gz"
    seqkit grep -f "$POS" "$filtered" -o "$outdir/positive.fastq.gz"

    stats_line="$(seqkit stats -T "$outdir/positive.fastq.gz" | tail -n 1)"
    stats="$(echo "$stats_line" | cut -f4-)"

    (
      flock -x 200
      echo -e "$DATASET\t$base\tpositive\t$stats" >> "$SUMMARY_TABLE"
    ) 200>>"$LOCK_FILE"

  elif [[ "$TYPE" == "NAS" ]]; then

    echo "CMD: seqkit grep -f $POS $filtered -o $outdir/positive.fastq.gz"
    seqkit grep -f "$POS" "$filtered" -o "$outdir/positive.fastq.gz"

    echo "CMD: seqkit grep -f $UNB $filtered -o $outdir/unblock.fastq.gz"
    seqkit grep -f "$UNB" "$filtered" -o "$outdir/unblock.fastq.gz"

    stats_line="$(seqkit stats -T "$outdir/positive.fastq.gz" | tail -n 1)"
    stats="$(echo "$stats_line" | cut -f4-)"

    (
      flock -x 200
      echo -e "$DATASET\t$base\tpositive\t$stats" >> "$SUMMARY_TABLE"
    ) 200>>"$LOCK_FILE"

    stats_line="$(seqkit stats -T "$outdir/unblock.fastq.gz" | tail -n 1)"
    stats="$(echo "$stats_line" | cut -f4-)"

    (
      flock -x 200
      echo -e "$DATASET\t$base\tunblock\t$stats" >> "$SUMMARY_TABLE"
    ) 200>>"$LOCK_FILE"
  fi

  echo "DONE $base"
}

export -f process_fastq

# ---------- run the single merged FASTQ ----------
xargs -P 1 -I {} bash -c 'process_fastq "$@"' _ {} < "$FASTQ_LIST"

echo "===== MERGED FILTERING COMPLETE ====="
echo "DATE=$(date)"

rm -f "$LOCK_FILE"
