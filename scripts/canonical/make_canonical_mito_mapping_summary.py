#!/usr/bin/env python3
"""Create canonical mitochondrial alignment and species-summary tables.

Definitions
-----------
Complete SAM SEQ length describes the sequenced read.
CIGAR M + = + X describes the query sequence aligned to the mitochondrial reference.
Mitochondrial relative abundance is calculated from M + = + X bases.

Inputs
------
canonical_data/inputs/dataset_manifest.tsv with columns:
dataset, read_subset, experiment, bam

Outputs
-------
canonical_data/outputs/canonical_mapped_alignments.tsv.gz
canonical_data/outputs/canonical_mito_mapping_summary.tsv
canonical_data/logs/canonical_mito_mapping_summary_run.log
canonical_data/logs/canonical_mito_mapping_summary_provenance.json
canonical_data/scripts/make_canonical_mito_mapping_summary.py
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

CIGAR_TOKEN = re.compile(r"(\d+)([MIDNSHP=X])")
EXPECTED_REFERENCE_COUNT = 16


def now_iso():
    return datetime.now(timezone.utc).astimezone().isoformat()


def sha256_file(path, chunk_size=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def run_text(command):
    completed = subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return completed.stdout


def parse_cigar(cigar):
    tokens = CIGAR_TOKEN.findall(cigar)
    reconstructed = "".join(length + op for length, op in tokens)
    if not tokens or reconstructed != cigar:
        raise ValueError(f"Malformed or unsupported CIGAR: {cigar}")

    counts = defaultdict(int)
    for length_text, operation in tokens:
        counts[operation] += int(length_text)

    aligned_query = counts["M"] + counts["="] + counts["X"]
    query_consumed = aligned_query + counts["I"] + counts["S"]
    reference_consumed = aligned_query + counts["D"] + counts["N"]

    return {
        "aligned_query": aligned_query,
        "soft": counts["S"],
        "inserted": counts["I"],
        "deleted": counts["D"],
        "hard": counts["H"],
        "query_consumed": query_consumed,
        "reference_consumed": reference_consumed,
    }


def read_manifest(path):
    required = {"dataset", "read_subset", "experiment", "bam"}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                "Manifest must contain: dataset, read_subset, experiment, bam"
            )
        rows = list(reader)

    if not rows:
        raise ValueError("Manifest contains no datasets")

    labels = [row["dataset"] for row in rows]
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise ValueError(f"Duplicate dataset labels: {duplicates}")

    for row in rows:
        bam = Path(row["bam"])
        if not bam.is_file():
            raise FileNotFoundError(f"Missing BAM: {bam}")
        if not Path(str(bam) + ".bai").is_file():
            raise FileNotFoundError(f"Missing BAM index: {bam}.bai")

    return rows


def reference_lengths(bam):
    header = run_text(["samtools", "view", "-H", str(bam)])
    refs = {}
    order = []
    for line in header.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        fields = line.split("\t")
        name = next((x[3:] for x in fields if x.startswith("SN:")), None)
        length = next((x[3:] for x in fields if x.startswith("LN:")), None)
        if name is None or length is None:
            raise ValueError(f"Incomplete @SQ line in {bam}: {line}")
        refs[name] = int(length)
        order.append(name)
    if len(refs) != EXPECTED_REFERENCE_COUNT:
        raise ValueError(
            f"Expected {EXPECTED_REFERENCE_COUNT} references in {bam}, found {len(refs)}"
        )
    return refs, order


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-dir",
        required=True,
    )
    args = parser.parse_args()

    project = Path(args.project_dir).resolve()
    canonical = project / "canonical_data"
    manifest = canonical / "inputs" / "dataset_manifest.tsv"
    outputs = canonical / "outputs"
    logs = canonical / "logs"
    scripts = canonical / "scripts"
    outputs.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    scripts.mkdir(parents=True, exist_ok=True)

    read_output = outputs / "canonical_mapped_alignments.tsv.gz"
    summary_output = outputs / "canonical_mito_mapping_summary.tsv"
    log_path = logs / "canonical_mito_mapping_summary_run.log"
    provenance_path = logs / "canonical_mito_mapping_summary_provenance.json"
    frozen_script = Path(__file__).resolve()

    datasets = read_manifest(manifest)
    samtools_version = run_text(["samtools", "--version"]).splitlines()[0]

    read_fields = [
        "dataset", "read_subset", "experiment", "read_id", "reference",
        "reference_length", "mapq", "cigar", "alignment_start_1based",
        "alignment_end_1based", "mapped_read_sequence_bases",
        "cigar_aligned_query_bases", "soft_clipped_query_bases",
        "inserted_query_bases", "deleted_reference_bases",
        "hard_clipped_bases",
    ]

    dataset_stats = {}
    reference_signature = None
    started = now_iso()

    with open(log_path, "w") as log, gzip.open(
        read_output, "wt", newline=""
    ) as read_handle:
        read_writer = csv.DictWriter(
            read_handle, delimiter="\t", fieldnames=read_fields
        )
        read_writer.writeheader()

        log.write(f"Started: {started}\n")
        log.write(f"Samtools: {samtools_version}\n")
        log.write("Definition: SEQ length describes the sequenced read.\n")
        log.write("Definition: M + = + X describes aligned query bases.\n")
        log.write("RA denominator: total M + = + X bases within each dataset.\n")

        for item in datasets:
            dataset = item["dataset"]
            bam = Path(item["bam"])
            refs, ref_order = reference_lengths(bam)
            signature = tuple((name, refs[name]) for name in ref_order)
            if reference_signature is None:
                reference_signature = signature
            elif signature != reference_signature:
                raise ValueError(f"Reference list differs for dataset {dataset}")

            log.write(f"[DATASET] {dataset}\n")
            log.write(f"[BAM] {bam}\n")
            log.write(f"[CMD] samtools view {bam}\n")

            by_ref = {
                ref: {
                    "mapped_reads": 0,
                    "seq": 0,
                    "aligned": 0,
                    "soft": 0,
                    "inserted": 0,
                    "deleted": 0,
                    "hard": 0,
                }
                for ref in ref_order
            }

            seen = set()
            duplicate_read_ids = 0
            accounting_failures = 0
            total_alignments = 0

            process = subprocess.Popen(
                ["samtools", "view", str(bam)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            assert process.stdout is not None
            for line in process.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    raise ValueError(f"SAM record has fewer than 11 fields: {line[:200]}")

                read_id = fields[0]
                flag = int(fields[1])
                ref = fields[2]
                pos = int(fields[3])
                mapq = int(fields[4])
                cigar = fields[5]
                seq = fields[9]

                if flag & 2308:
                    raise ValueError(
                        f"Dataset {dataset} contains alignment excluded by flag 2308: {read_id}"
                    )
                if mapq < 20:
                    raise ValueError(
                        f"Dataset {dataset} contains MAPQ below 20: {read_id}"
                    )
                if ref not in refs:
                    raise ValueError(f"Unexpected reference in {dataset}: {ref}")
                if cigar == "*" or seq == "*":
                    raise ValueError(f"Missing CIGAR or SEQ in {dataset}: {read_id}")

                if read_id in seen:
                    duplicate_read_ids += 1
                else:
                    seen.add(read_id)

                cigar_counts = parse_cigar(cigar)
                seq_length = len(seq)
                if seq_length != cigar_counts["query_consumed"]:
                    accounting_failures += 1
                    raise ValueError(
                        f"SEQ/CIGAR accounting failure for {read_id}: "
                        f"SEQ={seq_length}, CIGAR query={cigar_counts['query_consumed']}"
                    )

                reference_end = pos + cigar_counts["reference_consumed"] - 1
                values = by_ref[ref]
                values["mapped_reads"] += 1
                values["seq"] += seq_length
                values["aligned"] += cigar_counts["aligned_query"]
                values["soft"] += cigar_counts["soft"]
                values["inserted"] += cigar_counts["inserted"]
                values["deleted"] += cigar_counts["deleted"]
                values["hard"] += cigar_counts["hard"]
                total_alignments += 1

                read_writer.writerow({
                    "dataset": dataset,
                    "read_subset": item["read_subset"],
                    "experiment": item["experiment"],
                    "read_id": read_id,
                    "reference": ref,
                    "reference_length": refs[ref],
                    "mapq": mapq,
                    "cigar": cigar,
                    "alignment_start_1based": pos,
                    "alignment_end_1based": reference_end,
                    "mapped_read_sequence_bases": seq_length,
                    "cigar_aligned_query_bases": cigar_counts["aligned_query"],
                    "soft_clipped_query_bases": cigar_counts["soft"],
                    "inserted_query_bases": cigar_counts["inserted"],
                    "deleted_reference_bases": cigar_counts["deleted"],
                    "hard_clipped_bases": cigar_counts["hard"],
                })

            stderr = process.stderr.read() if process.stderr is not None else ""
            return_code = process.wait()
            if stderr:
                log.write(stderr)
            if return_code != 0:
                raise RuntimeError(f"samtools view failed for {bam}")
            if duplicate_read_ids:
                raise ValueError(
                    f"Dataset {dataset} contains {duplicate_read_ids} duplicate primary read IDs"
                )

            total_aligned = sum(x["aligned"] for x in by_ref.values())
            total_seq = sum(x["seq"] for x in by_ref.values())
            if total_aligned == 0:
                raise ValueError(f"Dataset {dataset} has zero aligned query bases")

            dataset_stats[dataset] = {
                "read_subset": item["read_subset"],
                "experiment": item["experiment"],
                "bam": str(bam),
                "bam_size_bytes": bam.stat().st_size,
                "bam_mtime": datetime.fromtimestamp(
                    bam.stat().st_mtime, timezone.utc
                ).isoformat(),
                "reference_count": len(refs),
                "total_alignments": total_alignments,
                "total_mapped_read_sequence_bases": total_seq,
                "total_cigar_aligned_query_bases": total_aligned,
                "total_soft_clipped_query_bases": sum(x["soft"] for x in by_ref.values()),
                "total_inserted_query_bases": sum(x["inserted"] for x in by_ref.values()),
                "duplicate_primary_read_ids": duplicate_read_ids,
                "cigar_seq_accounting_failures": accounting_failures,
                "by_reference": by_ref,
            }

            log.write(
                f"[OK] {dataset}: reads={total_alignments}, "
                f"SEQ={total_seq}, M_eq_X={total_aligned}\n"
            )

    summary_fields = [
        "dataset", "read_subset", "experiment", "reference",
        "reference_length", "mapped_reads", "mapped_read_sequence_bases",
        "cigar_aligned_query_bases", "soft_clipped_query_bases",
        "inserted_query_bases", "deleted_reference_bases",
        "hard_clipped_bases", "mean_mapped_read_length",
        "mitochondrial_RA_fraction", "mitochondrial_RA_percent",
    ]
    summary_rows = []

    reference_map = dict(reference_signature)
    for item in datasets:
        dataset = item["dataset"]
        stats = dataset_stats[dataset]
        total_aligned = stats["total_cigar_aligned_query_bases"]
        for ref, ref_length in reference_signature:
            values = stats["by_reference"][ref]
            mapped_reads = values["mapped_reads"]
            fraction = values["aligned"] / total_aligned
            mean_length = values["seq"] / mapped_reads if mapped_reads else 0.0
            summary_rows.append({
                "dataset": dataset,
                "read_subset": item["read_subset"],
                "experiment": item["experiment"],
                "reference": ref,
                "reference_length": ref_length,
                "mapped_reads": mapped_reads,
                "mapped_read_sequence_bases": values["seq"],
                "cigar_aligned_query_bases": values["aligned"],
                "soft_clipped_query_bases": values["soft"],
                "inserted_query_bases": values["inserted"],
                "deleted_reference_bases": values["deleted"],
                "hard_clipped_bases": values["hard"],
                "mean_mapped_read_length": f"{mean_length:.10f}",
                "mitochondrial_RA_fraction": f"{fraction:.12f}",
                "mitochondrial_RA_percent": f"{100.0 * fraction:.10f}",
            })

        ra_sum = sum(
            float(row["mitochondrial_RA_percent"])
            for row in summary_rows
            if row["dataset"] == dataset
        )
        dataset_stats[dataset]["mitochondrial_RA_percent_sum"] = ra_sum
        if abs(ra_sum - 100.0) > 1e-6:
            raise ValueError(f"RA does not sum to 100 for {dataset}: {ra_sum}")

    write_tsv(summary_output, summary_rows, summary_fields)

    expected_benchmarks = {
        "full_meta": {"seq": 45777683, "aligned": 33912560},
        "full_nas_accepted": {"seq": 25941340, "aligned": 24344637},
        "full_nas_rejected": {"seq": 29280653, "aligned": 27890945},
    }
    benchmark_results = {}
    for dataset, expected in expected_benchmarks.items():
        if dataset not in dataset_stats:
            continue
        observed_seq = dataset_stats[dataset]["total_mapped_read_sequence_bases"]
        observed_aligned = dataset_stats[dataset]["total_cigar_aligned_query_bases"]
        passed = observed_seq == expected["seq"] and observed_aligned == expected["aligned"]
        benchmark_results[dataset] = {
            "expected": expected,
            "observed": {"seq": observed_seq, "aligned": observed_aligned},
            "passed": passed,
        }
        if not passed:
            print(
                f"[HISTORICAL BENCHMARK DIFFERENCE] {dataset}: "
                f"expected_seq={expected['seq']} "
                f"observed_seq={observed_seq} "
                f"expected_aligned={expected['aligned']} "
                f"observed_aligned={observed_aligned}"
            )

    provenance = {
        "analysis_name": "canonical_mito_mapping_summary",
        "started": started,
        "finished": now_iso(),
        "project_dir": str(project),
        "manifest": str(manifest),
        "manifest_sha256": sha256_file(manifest),
        "samtools_version": samtools_version,
        "python_version": sys.version,
        "hostname": os.uname().nodename,
        "base_definitions": {
            "mapped_read_sequence_bases": "Complete SAM SEQ length",
            "cigar_aligned_query_bases": "Sum of M, = and X CIGAR operations",
            "mitochondrial_RA_denominator": (
                "Total M, = and X query bases across all 16 references within each dataset"
            ),
        },
        "input_filter_expectation": {
            "minimum_mapq": 20,
            "excluded_flag_mask": 2308,
            "filtering_applied_upstream": True,
        },
        "datasets": dataset_stats,
        "independent_benchmark_results": benchmark_results,
        "outputs": {
            str(read_output): sha256_file(read_output),
            str(summary_output): sha256_file(summary_output),
            str(frozen_script): sha256_file(frozen_script),
        },
        "log": str(log_path),
    }

    with open(provenance_path, "w") as handle:
        json.dump(provenance, handle, indent=2)

    with open(log_path, "a") as log:
        log.write(f"Finished: {provenance['finished']}\n")
        log.write(f"[OUTPUT] {read_output}\n")
        log.write(f"[OUTPUT] {summary_output}\n")
        log.write(f"[PROVENANCE] {provenance_path}\n")

    print(f"Wrote: {read_output}")
    print(f"Wrote: {summary_output}")
    print(f"Wrote: {provenance_path}")
    print(f"Wrote: {log_path}")


if __name__ == "__main__":
    main()
