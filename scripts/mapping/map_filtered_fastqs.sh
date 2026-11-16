#!/bin/bash

set -euo pipefail

# Paths are resolved from the reproduction workspace.
unset SLURM_SUBMIT_DIR
mkdir -p logs


if [[ -z "${1-}" ]]; then
  echo "Usage: bash scripts/mapping/map_filtered_fastqs.sh <dataset_dir>"
  exit 1
fi

DATASET="$1"

REF_MMI="refs/mock_mitogenomes_clean.mmi"
REF_FASTA="refs/mock_mitogenomes_clean.fasta"

MAPQ=20
EXCLUDE_FLAGS=2308

MINIMAP_THREADS=12
VIEW_THREADS=4
SORT_THREADS=4

MAX_CONCURRENT="${MAX_CONCURRENT:-3}"

FILTERED_DIR="$DATASET/filtered"
LOG_DIR="$DATASET/logs"
TASKFILE="$LOG_DIR/map_filtered_tasks.tsv"
INPUT_MANIFEST="$LOG_DIR/mapping_inputs.txt"
PROVENANCE_LOG="$LOG_DIR/mapping_provenance.txt"

mkdir -p "$LOG_DIR"

add_task() {
  fq="$1"

  dir="$(dirname "$fq")"
  base="$(basename "$fq" .fastq.gz)"
  outdir="${dir}/bams"

  mkdir -p "$outdir"

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$fq" \
    "${outdir}/${base}.strict.bam" \
    "${outdir}/${base}.flagstat.txt" \
    "${outdir}/${base}.idxstats.txt" \
    "${outdir}/${base}.minimap2.log" \
    "${outdir}/${base}.mapping_command.txt" \
    >> "$TASKFILE"
}

############################
# DRIVER MODE
############################
if [[ -z "${PAPER1_LOCAL_TASK_ID:-}" ]]; then

  [[ -d "$FILTERED_DIR" ]] || { echo "ERROR: filtered dir not found: $FILTERED_DIR"; exit 1; }
  [[ -f "$REF_MMI" ]] || { echo "ERROR: missing reference index: $REF_MMI"; exit 1; }
  [[ -f "$REF_FASTA" ]] || { echo "ERROR: missing reference FASTA: $REF_FASTA"; exit 1; }

  command -v minimap2 >/dev/null 2>&1 || { echo "ERROR: minimap2 missing"; exit 1; }
  command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools missing"; exit 1; }

  rm -f "$TASKFILE" "$INPUT_MANIFEST" "$PROVENANCE_LOG"

  echo "===== MAPPING PROVENANCE =====" | tee "$PROVENANCE_LOG"
  echo "DATE=$(date)" | tee -a "$PROVENANCE_LOG"
  echo "HOST=$(hostname)" | tee -a "$PROVENANCE_LOG"
  echo "DATASET=$DATASET" | tee -a "$PROVENANCE_LOG"
  echo "FILTERED_DIR=$FILTERED_DIR" | tee -a "$PROVENANCE_LOG"
  echo "REF_MMI=$REF_MMI" | tee -a "$PROVENANCE_LOG"
  echo "REF_FASTA=$REF_FASTA" | tee -a "$PROVENANCE_LOG"
  echo "INPUT_FORMAT=FASTQ.GZ" | tee -a "$PROVENANCE_LOG"
  echo "MINIMAP2_VERSION=$(minimap2 --version)" | tee -a "$PROVENANCE_LOG"
  echo "SAMTOOLS_VERSION=$(samtools --version | head -n1)" | tee -a "$PROVENANCE_LOG"
  echo "MAP_PRESET=map-ont" | tee -a "$PROVENANCE_LOG"
  echo "MAPQ=$MAPQ" | tee -a "$PROVENANCE_LOG"
  echo "EXCLUDE_FLAGS=$EXCLUDE_FLAGS" | tee -a "$PROVENANCE_LOG"
  echo "MINIMAP_THREADS=$MINIMAP_THREADS" | tee -a "$PROVENANCE_LOG"
  echo "VIEW_THREADS=$VIEW_THREADS" | tee -a "$PROVENANCE_LOG"
  echo "SORT_THREADS=$SORT_THREADS" | tee -a "$PROVENANCE_LOG"
  echo "MAX_CONCURRENT=$MAX_CONCURRENT" | tee -a "$PROVENANCE_LOG"

  echo "===== INPUT FASTQS =====" | tee -a "$PROVENANCE_LOG"

  if [[ "$DATASET" == *full-nas-dataset* ]]; then

    find "$FILTERED_DIR" -type f \
      \( -name "positive.fastq.gz" -o -name "unblock.fastq.gz" \) \
      | sort > "$INPUT_MANIFEST"

  elif [[ "$DATASET" == *half-flowcell-dataset* ]]; then

    find "$FILTERED_DIR" -type f \
      \( \
        -name "META.fastq.gz" -o \
        -name "NAS.positive.fastq.gz" -o \
        -name "NAS.unblock.fastq.gz" \
      \) | sort > "$INPUT_MANIFEST"

  elif [[ "$DATASET" == *full-meta-dataset* ]]; then

    find "$FILTERED_DIR" -type f \
      -name "positive.fastq.gz" \
      | sort > "$INPUT_MANIFEST"

  else
    echo "ERROR: unknown dataset type: $DATASET"
    exit 1
  fi

  cat "$INPUT_MANIFEST" | tee -a "$PROVENANCE_LOG"

  while IFS= read -r fq; do
    add_task "$fq"
  done < "$INPUT_MANIFEST"

  ntasks=$(wc -l < "$TASKFILE")
  [[ "$ntasks" -gt 0 ]] || { echo "ERROR: no mapping tasks created"; exit 1; }

  echo "TASK_COUNT=$ntasks" | tee -a "$PROVENANCE_LOG"
  echo "TASKFILE=$TASKFILE" | tee -a "$PROVENANCE_LOG"

  echo "===== TASKFILE =====" | tee -a "$PROVENANCE_LOG"
  column -t -s $'\t' "$TASKFILE" 2>/dev/null | tee -a "$PROVENANCE_LOG" || cat "$TASKFILE" | tee -a "$PROVENANCE_LOG"

  for ((task_id=0; task_id<ntasks; task_id++)); do
    PAPER1_LOCAL_TASK_ID="$task_id" bash "$0" "$DATASET"
  done

  exit 0
fi

############################
# ARRAY MODE
############################

line=$(sed -n "$((PAPER1_LOCAL_TASK_ID + 1))p" "$TASKFILE")
IFS=$'\t' read -r fastq bam flagstat idxstats minimap_log cmd_log <<< "$line"

echo "===== TASK PROVENANCE ====="
echo "TASK_ID=$PAPER1_LOCAL_TASK_ID"
echo "FASTQ=$fastq"
echo "BAM=$bam"
echo "FLAGSTAT=$flagstat"
echo "IDXSTATS=$idxstats"
echo "MINIMAP_LOG=$minimap_log"
echo "CMD_LOG=$cmd_log"
echo "HOST=$(hostname)"
echo "START=$(date)"

[[ -f "$fastq" ]] || { echo "ERROR: FASTQ missing: $fastq"; exit 2; }
[[ -f "$REF_MMI" ]] || { echo "ERROR: reference index missing: $REF_MMI"; exit 3; }

echo "FASTQ_SIZE:"
du -sh "$fastq"

mkdir -p "$(dirname "$bam")"

cat > "$cmd_log" <<EOF
minimap2 -t $MINIMAP_THREADS -ax map-ont --secondary=no $REF_MMI $fastq |
samtools view -@ $VIEW_THREADS -b -q $MAPQ -F $EXCLUDE_FLAGS |
samtools sort -@ $SORT_THREADS -o $bam -
samtools index $bam
samtools flagstat $bam > $flagstat
samtools idxstats $bam > $idxstats
samtools view -c $bam
EOF

echo "===== COMMAND ====="
cat "$cmd_log"

set +e
set -o pipefail

minimap2 -t "$MINIMAP_THREADS" -ax map-ont --secondary=no "$REF_MMI" "$fastq" 2> "$minimap_log" \
  | samtools view -@ "$VIEW_THREADS" -b -q "$MAPQ" -F "$EXCLUDE_FLAGS" \
  | samtools sort -@ "$SORT_THREADS" -o "$bam" -

exit_code=$?

set -e
set +o pipefail

echo "PIPELINE_EXIT_CODE=$exit_code"

if [[ "$exit_code" -ne 0 ]]; then
  echo "ERROR: mapping failed for $fastq"
  exit "$exit_code"
fi

samtools index "$bam"
samtools flagstat "$bam" > "$flagstat"
samtools idxstats "$bam" > "$idxstats"

echo "BAM_SIZE:"
du -sh "$bam"

echo "MAPPED_READS_STRICT_BAM:"
samtools view -c "$bam"

echo "END=$(date)"
