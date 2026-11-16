#!/usr/bin/env python3
"""Create canonical_run_summary.tsv/csv for the adapter-trimmed reanalysis.

Raw-output metrics are calculated from pre-demultiplexing
offline Dorado SUP-rebasecalled pass-read summaries:
  * pass-read count
  * pass-base yield
  * read N50

Filtered metrics are obtained from adapter-trimmed,
length-filtered FASTQs with `seqkit stats -a -T`:
  * average read length
  * AvgQual

Mapped metrics are aggregated from filtered-read BAMs mapped
to the mitochondrial reference sequences.
Provenance is written separately under canonical_data/logs/.
"""

import argparse
import csv
import hashlib
import json
import platform
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

TARGETS = {
    "NC_022727_Diuraphis_noxia",
    "NC_029727_Myzus_persicae",
    "NC_069039_Trichogramma_pretiosum",
}

OUTPUT_FIELDS = [
    "canonical_dataset", "experiment", "sample", "replicate",
    "actual_barcode", "replicate_type", "initial_pore_count",
    "initial_pore_count_type", "raw_pass_read_count_millions",
    "raw_pass_base_yield_gb", "raw_read_N50_bp",
    "filtered_read_count_millions", "filtered_mean_read_length_bp", "filtered_mean_read_quality",
    "on_target_reads", "on_target_cigar_aligned_query_bases",
    "on_target_percent_of_mapped_mitochondrial_aligned_bases",
    "on_target_read_scope",
]

# Order must match run_statistics_manifest.tsv.
ROW_SPECS = {
    "half_meta": {
        "filtered": "half-flowcell-dataset/filtered/unclassified/META.fastq.gz",
        "run": "half_flowcell", "half": "meta",
        "raw_sup_scope": "half_meta",
        "mito_main": ["half_meta"], "scope": "all_metagenomic_reads",
    },
    "half_nas_combined": {
        "filtered": "half-flowcell-dataset/filtered/unclassified/NAS.fastq.gz",
        "run": "half_flowcell", "half": "nas",
        "raw_sup_scope": "half_nas",
        "mito_main": ["half_nas_accepted"],
        "mito_denominator_main": [
            "half_nas_accepted",
            "half_nas_rejected",
        ],
        "scope": (
            "NAS_accepted_targets_over_all_NAS_"
            "mitochondrial_bases"
        ),
    },
    "full_meta": {
        "filtered": "full-meta-dataset/filtered/full-meta-merged/len200_20000.fastq.gz",
        "run": "full_meta",
        "raw_sup_scope": "full_meta",
        "mito_main": ["full_meta"],
        "scope": "all_metagenomic_reads",
    },
    "full_meta_rep1": {
        "filtered": "full-meta-dataset/filtered/barcode01/len200_20000.fastq.gz",
        "run": "full_meta", "mito_barcode": ["full_meta_rep1"],
        "scope": "all_metagenomic_reads",
    },
    "full_meta_rep2": {
        "filtered": "full-meta-dataset/filtered/barcode02/len200_20000.fastq.gz",
        "run": "full_meta", "mito_barcode": ["full_meta_rep2"],
        "scope": "all_metagenomic_reads",
    },
    "full_meta_rep3": {
        "filtered": "full-meta-dataset/filtered/barcode03/len200_20000.fastq.gz",
        "run": "full_meta", "mito_barcode": ["full_meta_rep3"],
        "scope": "all_metagenomic_reads",
    },
    "full_nas_combined": {
        "filtered": "full-nas-dataset/filtered/full-nas-merged/len200_20000.fastq.gz",
        "run": "full_nas",
        "raw_sup_scope": "full_nas",
        "mito_main": ["full_nas_accepted"],
        "mito_denominator_main": [
            "full_nas_accepted",
            "full_nas_rejected",
        ],
        "scope": (
            "NAS_accepted_targets_over_all_NAS_"
            "mitochondrial_bases"
        ),
    },
    "full_nas_rep1_combined": {
        "filtered": "full-nas-dataset/filtered/barcode09/len200_20000.fastq.gz",
        "run": "full_nas",
        "mito_barcode": ["full_nas_rep1_accepted"],
        "mito_denominator_barcode": [
            "full_nas_rep1_accepted",
            "full_nas_rep1_rejected",
        ],
        "scope": (
            "NAS_accepted_targets_over_all_NAS_"
            "mitochondrial_bases"
        ),
    },
    "full_nas_rep2_combined": {
        "filtered": "full-nas-dataset/filtered/barcode10/len200_20000.fastq.gz",
        "run": "full_nas",
        "mito_barcode": ["full_nas_rep2_accepted"],
        "mito_denominator_barcode": [
            "full_nas_rep2_accepted",
            "full_nas_rep2_rejected",
        ],
        "scope": (
            "NAS_accepted_targets_over_all_NAS_"
            "mitochondrial_bases"
        ),
    },
    "full_nas_rep3_combined": {
        "filtered": "full-nas-dataset/filtered/barcode11/len200_20000.fastq.gz",
        "run": "full_nas",
        "mito_barcode": ["full_nas_rep3_accepted"],
        "mito_denominator_barcode": [
            "full_nas_rep3_accepted",
            "full_nas_rep3_rejected",
        ],
        "scope": (
            "NAS_accepted_targets_over_all_NAS_"
            "mitochondrial_bases"
        ),
    },
}

def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_command(command):
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            "Command failed:\n" + " ".join(map(str, command)) +
            "\nSTDOUT:\n" + result.stdout + "\nSTDERR:\n" + result.stderr
        )
    return result.stdout


def seqkit_stats(path):
    text = run_command(["seqkit", "stats", "-a", "-T", str(path)])
    records = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if len(records) != 1:
        raise ValueError(f"Expected one seqkit stats row for {path}; found {len(records)}")
    row = records[0]
    return {
        "num_seqs": int(row["num_seqs"]),
        "sum_len": int(row["sum_len"]),
        "N50": int(row["N50"]),
        "avg_len": float(row["avg_len"]),
        "AvgQual": float(row["AvgQual"]),
    }


def calculate_n50(lengths):
    if lengths.size == 0:
        raise ValueError(
            "Cannot calculate N50 from an empty length array"
        )

    ordered = np.sort(lengths)[::-1]

    cumulative = np.cumsum(
        ordered,
        dtype=np.int64,
    )

    threshold = int(ordered.sum()) / 2

    index = int(
        np.searchsorted(
            cumulative,
            threshold,
            side="left",
        )
    )

    return int(ordered[index])


def raw_sup_stats(path, scope):
    print(
        f"[progress] Calculating raw SUP metrics "
        f"for {scope} from {path}",
        flush=True,
    )

    length_chunks = []
    read_count = 0
    base_count = 0

    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=[
            "channel",
            "sequence_length_template",
        ],
        dtype={
            "channel": "int64",
            "sequence_length_template": "int64",
        },
        chunksize=500_000,
    ):
        if scope == "half_nas":
            chunk = chunk.loc[
                chunk["channel"].between(1, 256)
            ]

        elif scope == "half_meta":
            chunk = chunk.loc[
                chunk["channel"].between(257, 512)
            ]

        elif scope not in {
            "full_meta",
            "full_nas",
        }:
            raise ValueError(
                f"Unsupported raw SUP scope: {scope}"
            )

        lengths = chunk[
            "sequence_length_template"
        ].to_numpy(
            dtype=np.int64,
            copy=True,
        )

        if len(lengths) == 0:
            continue

        if (lengths <= 0).any():
            raise ValueError(
                "Non-positive sequence_length_template "
                f"detected in {path} for {scope}"
            )

        length_chunks.append(lengths)
        read_count += int(len(lengths))
        base_count += int(lengths.sum())

    if not length_chunks:
        raise ValueError(
            f"No raw SUP reads found for {scope} in {path}"
        )

    all_lengths = np.concatenate(length_chunks)

    if len(all_lengths) != read_count:
        raise ValueError(
            f"Raw SUP read-count mismatch for {scope}"
        )

    if int(all_lengths.sum()) != base_count:
        raise ValueError(
            f"Raw SUP base-count mismatch for {scope}"
        )

    return {
        "num_seqs": read_count,
        "sum_len": base_count,
        "N50": calculate_n50(all_lengths),
    }


def find_total_pores(path):
    """Return the initial MinKNOW total_pores value.

    ONT report JSON stores total_pores as strings in chronological
    user_messages entries. The first valid value is the initial pore count.
    """
    with open(path) as handle:
        obj = json.load(handle)

    for message in obj.get("user_messages", []):
        if not isinstance(message, dict):
            continue

        extra_data = message.get("extra_data", {})

        if not isinstance(extra_data, dict):
            continue

        value = extra_data.get("total_pores")

        if value is None:
            continue

        try:
            pore_count = int(value)
        except (TypeError, ValueError):
            continue

        if pore_count > 0:
            return pore_count

    raise ValueError(
        "No valid user_messages[].extra_data.total_pores "
        f"value found in {path}"
    )

def allocate_split_pores(total_pores, sequencing_summary):
    observations = {"nas": set(), "meta": set()}
    with open(sequencing_summary, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"channel", "mux", "start_time"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{sequencing_summary} missing columns: {sorted(missing)}")
        for row in reader:
            if float(row["start_time"]) <= 600:
                channel = int(row["channel"])
                group = "nas" if channel <= 256 else "meta"
                observations[group].add((channel, int(row["mux"])))
    denominator = sum(len(values) for values in observations.values())
    if denominator == 0:
        raise ValueError("No channel/mux observations in first 10 minutes")
    allocated = {
        group: round(total_pores * len(values) / denominator)
        for group, values in observations.items()
    }
    return allocated, {group: len(values) for group, values in observations.items()}


def aggregate_mito(path):
    result = {}
    for row in read_tsv(path):
        dataset = row["dataset"]
        current = result.setdefault(dataset, {"target_reads": 0, "target_bases": 0, "all_bases": 0})
        bases = int(row["cigar_aligned_query_bases"])
        current["all_bases"] += bases
        if row["reference"] in TARGETS:
            current["target_reads"] += int(row["mapped_reads"])
            current["target_bases"] += bases
    return result


def combine_mito(
    source,
    numerator_dataset_names,
    denominator_dataset_names=None,
):
    if denominator_dataset_names is None:
        denominator_dataset_names = numerator_dataset_names

    combined = {
        "target_reads": 0,
        "target_bases": 0,
        "all_bases": 0,
    }

    for dataset in numerator_dataset_names:
        if dataset not in source:
            raise KeyError(
                f"Missing canonical mito numerator dataset: {dataset}"
            )

        combined["target_reads"] += source[dataset]["target_reads"]
        combined["target_bases"] += source[dataset]["target_bases"]

    for dataset in denominator_dataset_names:
        if dataset not in source:
            raise KeyError(
                f"Missing canonical mito denominator dataset: {dataset}"
            )

        combined["all_bases"] += source[dataset]["all_bases"]

    if combined["all_bases"] == 0:
        raise ValueError(
            "Zero mitochondrial aligned-base denominator for "
            f"{denominator_dataset_names}"
        )

    combined["percent"] = (
        100
        * combined["target_bases"]
        / combined["all_bases"]
    )

    return combined


def write_table(path, rows, delimiter):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-dir",
        required=True,
    )
    args = parser.parse_args()
    started = datetime.now().astimezone().isoformat()
    root = Path(args.project_dir).resolve()
    inputs_dir = root / "canonical_data" / "inputs"
    outputs_dir = root / "canonical_data" / "outputs"
    logs_dir = root / "canonical_data" / "logs"
    outputs_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    if shutil.which("seqkit") is None:
        raise RuntimeError("seqkit is not available on PATH")
    seqkit_version = run_command(["seqkit", "version"]).strip()

    run_manifest_path = inputs_dir / "run_statistics_manifest.tsv"
    ont_manifest_path = inputs_dir / "ont_report_manifest.tsv"
    main_mito_path = outputs_dir / "canonical_mito_mapping_summary.tsv"
    barcode_mito_path = outputs_dir / "canonical_barcode_mito_mapping_summary.tsv"

    raw_sup_dir = (
        root
        / "analysis"
        / "run_performance"
        / "inputs"
    )

    raw_sup_paths = {
        "full_meta": (
            raw_sup_dir
            / "full_meta_sequencing_summary_sup.tsv"
        ),
        "full_nas": (
            raw_sup_dir
            / "full_nas_sequencing_summary_sup.tsv"
        ),
        "half_flowcell": (
            raw_sup_dir
            / "half_flowcell_sequencing_summary_sup.tsv"
        ),
    }

    for raw_sup_path in raw_sup_paths.values():
        if not raw_sup_path.is_file():
            raise FileNotFoundError(raw_sup_path)

        if raw_sup_path.stat().st_size <= 0:
            raise ValueError(
                f"Empty raw SUP input: {raw_sup_path}"
            )

    manifest = read_tsv(run_manifest_path)
    observed_order = [row["canonical_dataset"] for row in manifest]
    expected_order = list(ROW_SPECS)
    if observed_order != expected_order:
        raise ValueError(
            "run_statistics_manifest row order differs from specification:\n"
            f"Observed: {observed_order}\nExpected: {expected_order}"
        )

    ont_manifest = {row["run"]: Path(row["report_json"]) for row in read_tsv(ont_manifest_path)}
    total_pores = {run: find_total_pores(path) for run, path in ont_manifest.items()}
    split_pores, early_counts = allocate_split_pores(
        total_pores["half_flowcell"],
        root / "half-flowcell-dataset" / "trimmed" / "sequencing_summary.tsv",
    )

    main_mito = aggregate_mito(main_mito_path)
    barcode_mito = aggregate_mito(barcode_mito_path)

    raw_sup_metrics = {
        "full_meta": raw_sup_stats(
            raw_sup_paths["full_meta"],
            "full_meta",
        ),
        "full_nas": raw_sup_stats(
            raw_sup_paths["full_nas"],
            "full_nas",
        ),
        "half_meta": raw_sup_stats(
            raw_sup_paths["half_flowcell"],
            "half_meta",
        ),
        "half_nas": raw_sup_stats(
            raw_sup_paths["half_flowcell"],
            "half_nas",
        ),
    }

    expected_raw_totals = {
        "full_meta": {
            "num_seqs": 8_914_639,
            "sum_len": 9_389_637_136,
        },
        "full_nas": {
            "num_seqs": 27_514_001,
            "sum_len": 11_191_661_832,
        },
        "half_meta": {
            "num_seqs": 2_792_100,
            "sum_len": 3_657_561_050,
        },
        "half_nas": {
            "num_seqs": 10_304_584,
            "sum_len": 4_030_959_029,
        },
    }

    for scope, expected in expected_raw_totals.items():
        observed = raw_sup_metrics[scope]

        if observed["num_seqs"] != expected["num_seqs"]:
            raise ValueError(
                f"Unexpected raw SUP read count for "
                f"{scope}: {observed['num_seqs']}"
            )

        if observed["sum_len"] != expected["sum_len"]:
            raise ValueError(
                f"Unexpected raw SUP base count for "
                f"{scope}: {observed['sum_len']}"
            )

    output_rows = []

    input_files = {
        run_manifest_path,
        ont_manifest_path,
        main_mito_path,
        barcode_mito_path,
        *raw_sup_paths.values(),
    }

    for metadata in manifest:
        dataset = metadata["canonical_dataset"]
        spec = ROW_SPECS[dataset]

        filtered = root / spec["filtered"]

        if not filtered.is_file():
            raise FileNotFoundError(filtered)

        input_files.add(filtered)

        raw_scope = spec.get("raw_sup_scope")

        if raw_scope is None:
            raw = None
        else:
            raw = raw_sup_metrics[raw_scope]

        filtered_stats = seqkit_stats(filtered)

        if "mito_main" in spec:
            mito = combine_mito(
                main_mito,
                spec["mito_main"],
                spec.get("mito_denominator_main"),
            )
        else:
            mito = combine_mito(
                barcode_mito,
                spec["mito_barcode"],
                spec.get("mito_denominator_barcode"),
            )

        if spec["run"] == "half_flowcell":
            pore_count = split_pores[spec["half"]]
        else:
            pore_count = total_pores[spec["run"]]

        output_rows.append({
            "canonical_dataset": dataset,
            "experiment": metadata["analysis_experiment"],
            "sample": metadata["analysis_sample"],
            "replicate": metadata["replicate"],
            "actual_barcode": metadata["actual_barcode"],
            "replicate_type": metadata["replicate_type"],
            "initial_pore_count": pore_count,
            "initial_pore_count_type": metadata["initial_pore_count_type"],
            "raw_pass_read_count_millions": (
                f"{raw['num_seqs'] / 1_000_000:.2f}"
                if raw is not None
                else ""
            ),
            "raw_pass_base_yield_gb": (
                f"{raw['sum_len'] / 1_000_000_000:.2f}"
                if raw is not None
                else ""
            ),
            "raw_read_N50_bp": (
                raw["N50"]
                if raw is not None
                else ""
            ),
            "filtered_read_count_millions": f"{filtered_stats['num_seqs'] / 1_000_000:.2f}",
            "filtered_mean_read_length_bp": f"{filtered_stats['avg_len']:.1f}",
            "filtered_mean_read_quality": f"{filtered_stats['AvgQual']:.2f}",
            "on_target_reads": mito["target_reads"],
            "on_target_cigar_aligned_query_bases": mito["target_bases"],
            "on_target_percent_of_mapped_mitochondrial_aligned_bases": f"{mito['percent']:.3f}",
            "on_target_read_scope": spec["scope"],
        })

    tsv_output = outputs_dir / "canonical_run_summary.tsv"
    csv_output = outputs_dir / "canonical_run_summary.csv"
    write_table(tsv_output, output_rows, "\t")
    write_table(csv_output, output_rows, ",")

    finished = datetime.now().astimezone().isoformat()
    provenance = {
        "analysis_name": "canonical_run_summary_reanalysis",
        "started": started,
        "finished": finished,
        "project_directory": str(root),
        "script": str(Path(__file__).resolve()),
        "rows": len(output_rows),
        "outputs": {
            str(tsv_output): sha256_file(tsv_output),
            str(csv_output): sha256_file(csv_output),
        },
        "metric_definitions": {
            "raw_pass_read_count_millions": (
                "read count from pre-demultiplexing offline "
                "SUP-rebasecalled pass-read summary divided by "
                "1e6; populated only for experiment-level rows"
            ),
            "raw_pass_base_yield_gb": (
                "sum of sequence_length_template from "
                "pre-demultiplexing offline SUP-rebasecalled "
                "pass-read summary divided by 1e9; populated "
                "only for experiment-level rows"
            ),
            "raw_read_N50_bp": (
                "N50 calculated from sequence_length_template "
                "in the pre-demultiplexing offline "
                "SUP-rebasecalled pass-read summary; populated "
                "only for experiment-level rows"
            ),
            "filtered_read_count_millions": (
                "seqkit num_seqs from adapter-trimmed, "
                "length-filtered FASTQ divided by 1e6"
            ),
            "filtered_mean_read_length_bp": (
                "seqkit avg_len from adapter-trimmed, "
                "length-filtered FASTQ"
            ),
            "filtered_mean_read_quality": (
                "seqkit AvgQual from adapter-trimmed, "
                "length-filtered FASTQ"
            ),
            "on_target_reads": (
                "sum of mapped_reads from filtered-read BAMs "
                "mapped to the three target mitochondrial "
                "reference sequences"
            ),
            "on_target_cigar_aligned_query_bases": (
                "sum of M+=+X query bases from filtered-read "
                "BAMs mapped to the three target mitochondrial "
                "reference sequences"
            ),
            "on_target_percent": (
                "META: 100 * target aligned query bases / all "
                "16-reference aligned query bases; NAS: 100 * "
                "accepted target aligned query bases / combined "
                "accepted and rejected aligned query bases across "
                "all 16 references"
            ),
            "whole_flowcell_initial_pores": "MinKNOW report JSON total_pores",
            "split_flowcell_initial_pores": "MinKNOW total allocated by first-10-minute unique channel/mux proportions",
        },
        "half_flowcell_channel_rules": {
            "nas": "channels 1--256",
            "meta": "channels 257--512",
        },
        "half_flowcell_early_channel_mux_counts": early_counts,
        "input_sha256": {
            str(path): sha256_file(path)
            for path in sorted(input_files, key=str)
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "seqkit": seqkit_version,
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }

    provenance_path = logs_dir / "canonical_run_summary_reanalysis_provenance.json"
    log_path = logs_dir / "canonical_run_summary_reanalysis_run.log"
    with open(provenance_path, "w") as handle:
        json.dump(provenance, handle, indent=2)
    with open(log_path, "w") as handle:
        handle.write(f"Started: {started}\n")
        handle.write(f"Finished: {finished}\n")
        handle.write(f"Rows: {len(output_rows)}\n")
        handle.write(f"TSV: {tsv_output}\n")
        handle.write(f"CSV: {csv_output}\n")
        handle.write(f"Provenance: {provenance_path}\n")
        handle.write(f"Seqkit: {seqkit_version}\n")
        handle.write(
            "Raw output: pre-demultiplexing offline "
            "SUP-rebasecalled pass reads\n"
        )
        handle.write(
            "Raw-output rows: half_meta, half_nas_combined, "
            "full_meta, full_nas_combined\n"
        )
        handle.write(
            "Barcode-row raw-output fields: blank\n"
        )
        handle.write(
            "Half-flowcell raw allocation: NAS channels "
            "1-256; META channels 257-512\n"
        )
        handle.write(
            "Filtered: adapter-trimmed, "
            "length-filtered FASTQs; no additional quality filter\n"
        )
        handle.write(
            "Mapped: filtered-read BAMs mapped to the "
            "mitochondrial reference sequences\n"
        )
        handle.write(
            "Filtered mean read quality: seqkit AvgQual\n"
        )

    print(f"Wrote: {tsv_output}")
    print(f"Wrote: {csv_output}")
    print(f"Wrote: {provenance_path}")
    print(f"Wrote: {log_path}")
    print(f"Rows: {len(output_rows)}")


if __name__ == "__main__":
    main()
