#!/bin/bash

set -euo pipefail

# Paths are resolved from the reproduction workspace.
unset SLURM_SUBMIT_DIR
mkdir -p logs

if [[ -z "${1-}" ]]; then
  echo "Usage: bash scripts/preprocessing/filter_fastqs.sh <dataset_dir>"
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

SUMMARY_TABLE="$DATASET/filtering_summary.tsv"
LOCK_FILE="$SUMMARY_TABLE.lock"
FASTQ_LIST="$DATASET/logs/input_fastqs.txt"

echo -e "dataset\tfile\tstep\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > "$SUMMARY_TABLE"

echo "===== FILTERING START ====="
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
echo "HALF_FLOWCELL_NAS_CHANNELS=1-256"
echo "HALF_FLOWCELL_META_CHANNELS=257-512"

# ---------- dataset type ----------
if [[ "$DATASET" == *half-flowcell* ]]; then
  TYPE="HALF"
elif [[ "$DATASET" == *nas* ]]; then
  TYPE="NAS"
else
  TYPE="META"
fi

echo "TYPE=$TYPE"

# ---------- FASTQ manifest ----------
find "$RUN_DIR"     -maxdepth 1     -name "*.fastq"     | sort > "$FASTQ_LIST"

FASTQ_COUNT="$(wc -l < "$FASTQ_LIST")"
echo "INPUT_FASTQ_COUNT=$FASTQ_COUNT"

if [[ "$FASTQ_COUNT" -eq 0 ]]; then
  echo "ERROR: no FASTQ files found under $RUN_DIR"
  exit 1
fi

# ---------- detect columns ----------
header="$(head -n 1 "$SUMMARY" | tr '\t' '\n')"

rid_col="$(echo "$header" | nl -ba | awk '$2=="read_id"{print $1}')"
er_col="$(echo "$header" | nl -ba | awk '$2=="end_reason"{print $1}')"
ch_col="$(echo "$header" | nl -ba | awk '$2=="channel"{print $1}')"

if [[ -z "$rid_col" || -z "$er_col" ]]; then
  echo "ERROR: read_id or end_reason column not found in $SUMMARY"
  exit 1
fi

if [[ "$TYPE" == "HALF" && -z "$ch_col" ]]; then
  echo "ERROR: channel column not found in $SUMMARY but required for half-flowcell dataset"
  exit 1
fi

echo "COLUMN_READ_ID=$rid_col"
echo "COLUMN_END_REASON=$er_col"
echo "COLUMN_CHANNEL=${ch_col:-NA}"

# ---------- build ID lists ----------
POS="$DATASET/positive.ids"
UNB="$DATASET/unblock.ids"
NAS_IDS="$DATASET/nas.ids"
META_IDS="$DATASET/meta.ids"

echo "CMD: build positive.ids"
awk -F'\t' -v r="$rid_col" -v e="$er_col" \
  'NR>1 && $e=="signal_positive"{print $r}' \
  "$SUMMARY" | sort -u > "$POS"

echo "CMD: build unblock.ids"
awk -F'\t' -v r="$rid_col" -v e="$er_col" \
  'NR>1 && $e=="data_service_unblock_mux_change"{print $r}' \
  "$SUMMARY" | sort -u > "$UNB"

if [[ "$TYPE" == "HALF" ]]; then
  echo "CMD: build nas.ids from channels 1-256"
  awk -F'\t' -v r="$rid_col" -v c="$ch_col" \
    'NR>1 && $c>=1 && $c<=256{print $r}' \
    "$SUMMARY" | sort -u > "$NAS_IDS"

  echo "CMD: build meta.ids from channels 257-512"
  awk -F'\t' -v r="$rid_col" -v c="$ch_col" \
    'NR>1 && $c>=257 && $c<=512{print $r}' \
    "$SUMMARY" | sort -u > "$META_IDS"
else
  : > "$NAS_IDS"
  : > "$META_IDS"
fi

echo "ID_COUNTS:"
echo "  POS=$(wc -l < "$POS")"
echo "  UNB=$(wc -l < "$UNB")"
echo "  NAS=$(wc -l < "$NAS_IDS")"
echo "  META=$(wc -l < "$META_IDS")"

# ---------- export for parallel ----------
export DATASET TYPE POS UNB NAS_IDS META_IDS SUMMARY_TABLE LOCK_FILE

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

  elif [[ "$TYPE" == "HALF" ]]; then

    echo "CMD: seqkit grep -f $NAS_IDS $filtered -o $outdir/NAS.fastq.gz"
    seqkit grep -f "$NAS_IDS" "$filtered" -o "$outdir/NAS.fastq.gz"

    echo "CMD: seqkit grep -f $META_IDS $filtered -o $outdir/META.fastq.gz"
    seqkit grep -f "$META_IDS" "$filtered" -o "$outdir/META.fastq.gz"

    stats_line="$(seqkit stats -T "$outdir/NAS.fastq.gz" | tail -n 1)"
    stats="$(echo "$stats_line" | cut -f4-)"

    (
      flock -x 200
      echo -e "$DATASET\t$base\tNAS_channel_split\t$stats" >> "$SUMMARY_TABLE"
    ) 200>>"$LOCK_FILE"

    stats_line="$(seqkit stats -T "$outdir/META.fastq.gz" | tail -n 1)"
    stats="$(echo "$stats_line" | cut -f4-)"

    (
      flock -x 200
      echo -e "$DATASET\t$base\tMETA_channel_split\t$stats" >> "$SUMMARY_TABLE"
    ) 200>>"$LOCK_FILE"

    echo "CMD: seqkit grep -f $POS $outdir/NAS.fastq.gz -o $outdir/NAS.positive.fastq.gz"
    seqkit grep -f "$POS" "$outdir/NAS.fastq.gz" -o "$outdir/NAS.positive.fastq.gz"

    echo "CMD: seqkit grep -f $UNB $outdir/NAS.fastq.gz -o $outdir/NAS.unblock.fastq.gz"
    seqkit grep -f "$UNB" "$outdir/NAS.fastq.gz" -o "$outdir/NAS.unblock.fastq.gz"

    stats_line="$(seqkit stats -T "$outdir/NAS.positive.fastq.gz" | tail -n 1)"
    stats="$(echo "$stats_line" | cut -f4-)"

    (
      flock -x 200
      echo -e "$DATASET\t$base\tNAS_positive\t$stats" >> "$SUMMARY_TABLE"
    ) 200>>"$LOCK_FILE"

    stats_line="$(seqkit stats -T "$outdir/NAS.unblock.fastq.gz" | tail -n 1)"
    stats="$(echo "$stats_line" | cut -f4-)"

    (
      flock -x 200
      echo -e "$DATASET\t$base\tNAS_unblock\t$stats" >> "$SUMMARY_TABLE"
    ) 200>>"$LOCK_FILE"
  fi

  echo "DONE $base"
}

export -f process_fastq

# ---------- parallel execution ----------
xargs -P 2 -I {} bash -c 'process_fastq "$@"' _ {} < "$FASTQ_LIST"

echo "===== FILTERING COMPLETE ====="
echo "DATE=$(date)"

rm -f "$LOCK_FILE"
