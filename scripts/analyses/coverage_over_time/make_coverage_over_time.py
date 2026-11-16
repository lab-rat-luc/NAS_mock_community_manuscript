#!/usr/bin/env python3

import json
import math
import platform
import re
import sys
from datetime import datetime
from pathlib import Path
import argparse

import numpy as np
import pandas as pd


_parser = argparse.ArgumentParser()
_parser.add_argument("--project-dir", required=True)
_args = _parser.parse_args()
PROJECT_DIR = Path(_args.project_dir).resolve()

ANALYSIS_DIR = (
    PROJECT_DIR
    / "canonical_data"
    / "analyses"
    / "coverage_over_time"
)

ALIGNMENT_FILE = (
    PROJECT_DIR
    / "canonical_data"
    / "outputs"
    / "canonical_mapped_alignments.tsv.gz"
)

READ_METADATA_FILE = (
    PROJECT_DIR
    / "canonical_data"
    / "outputs"
    / "canonical_read_metadata.tsv.gz"
)

REFERENCE_FILE = (
    PROJECT_DIR
    / "canonical_data"
    / "inputs"
    / "reference_manifest.tsv"
)

OUTPUT_DIR = ANALYSIS_DIR / "outputs"
LOG_DIR = ANALYSIS_DIR / "logs"

TIME_SERIES_OUTPUT = OUTPUT_DIR / "target_coverage_over_time.tsv"
THRESHOLD_OUTPUT = OUTPUT_DIR / "target_coverage_threshold_times.tsv"
FINAL_SUMMARY_OUTPUT = OUTPUT_DIR / "target_coverage_final_summary.tsv"
JOIN_VALIDATION_OUTPUT = OUTPUT_DIR / "target_timing_join_validation.tsv"
PROVENANCE_OUTPUT = LOG_DIR / "coverage_over_time_provenance.json"
LOG_OUTPUT = LOG_DIR / "coverage_over_time_run.log"

MAX_HOURS = 24

TARGET_REFS = {
    "NC_022727_Diuraphis_noxia",
    "NC_029727_Myzus_persicae",
    "NC_069039_Trichogramma_pretiosum",
}

DATASET_CONFIG = {
    "full_meta": {
        "method": "Metagenomic",
        "metadata_run": "full_meta",
        "metadata_decision_class": "not_applicable",
    },
    "full_nas_accepted": {
        "method": "NAS",
        "metadata_run": "full_nas",
        "metadata_decision_class": "accepted",
    },
}

EXPECTED_TARGET_ALIGNMENTS = {
    "full_meta": {
        "NC_022727_Diuraphis_noxia": 36,
        "NC_029727_Myzus_persicae": 297,
        "NC_069039_Trichogramma_pretiosum": 95,
    },
    "full_nas_accepted": {
        "NC_022727_Diuraphis_noxia": 111,
        "NC_029727_Myzus_persicae": 758,
        "NC_069039_Trichogramma_pretiosum": 281,
    },
}

EXPECTED_TOTAL_TARGET_ALIGNMENTS = 1578
BREADTH_THRESHOLDS = [50.0, 75.0, 90.0, 95.0, 99.0]
CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_reference_blocks(start_1based, cigar):
    """Return zero-based, half-open blocks covered by M, =, and X."""
    if not cigar or cigar == "*":
        return []

    reference_position = int(start_1based) - 1
    blocks = []
    tokens = CIGAR_PATTERN.findall(str(cigar))
    reconstructed = "".join(
        f"{length}{operation}" for length, operation in tokens
    )

    if reconstructed != cigar:
        raise ValueError(f"CIGAR could not be parsed completely: {cigar}")

    for length_text, operation in tokens:
        length = int(length_text)

        if operation in {"M", "=", "X"}:
            block_start = reference_position
            block_end = reference_position + length
            blocks.append((block_start, block_end))
            reference_position = block_end
        elif operation in {"D", "N"}:
            reference_position += length
        elif operation in {"I", "S", "H", "P"}:
            continue
        else:
            raise ValueError(f"Unexpected CIGAR operation: {operation}")

    return blocks


def load_target_alignments():
    required_columns = [
        "dataset",
        "read_id",
        "reference",
        "reference_length",
        "cigar",
        "alignment_start_1based",
        "alignment_end_1based",
        "cigar_aligned_query_bases",
    ]

    alignments = pd.read_csv(
        ALIGNMENT_FILE,
        sep="\t",
        compression="gzip",
        usecols=required_columns,
        dtype={
            "dataset": "string",
            "read_id": "string",
            "reference": "string",
            "reference_length": "int64",
            "cigar": "string",
            "alignment_start_1based": "int64",
            "alignment_end_1based": "int64",
            "cigar_aligned_query_bases": "int64",
        },
    )

    alignments = alignments.loc[
        alignments["dataset"].isin(DATASET_CONFIG)
        & alignments["reference"].isin(TARGET_REFS)
    ].copy()

    if len(alignments) != EXPECTED_TOTAL_TARGET_ALIGNMENTS:
        raise ValueError(
            f"Expected {EXPECTED_TOTAL_TARGET_ALIGNMENTS} target alignments, "
            f"observed {len(alignments)}"
        )

    observed_counts = (
        alignments.groupby(["dataset", "reference"], observed=True)
        .size()
        .to_dict()
    )

    for dataset, reference_counts in EXPECTED_TARGET_ALIGNMENTS.items():
        for reference, expected_count in reference_counts.items():
            observed_count = observed_counts.get((dataset, reference), 0)
            if observed_count != expected_count:
                raise ValueError(
                    f"Unexpected count for {dataset}, {reference}: "
                    f"expected {expected_count}, observed {observed_count}"
                )

    if alignments.duplicated(["dataset", "read_id"]).any():
        raise ValueError("Duplicate dataset and read_id records detected")

    if (alignments["cigar_aligned_query_bases"] <= 0).any():
        raise ValueError("Non-positive aligned-query-base values detected")

    return alignments


def load_reference_metadata():
    references = pd.read_csv(
        REFERENCE_FILE,
        sep="\t",
        dtype="string",
        usecols=[
            "reference",
            "species",
            "reference_length",
            "target_status",
        ],
    )

    references = references.loc[
        references["reference"].isin(TARGET_REFS)
    ].copy()

    if len(references) != 3:
        raise ValueError(
            f"Expected three target reference records, observed {len(references)}"
        )

    if references["reference"].duplicated().any():
        raise ValueError("Duplicate target references in manifest")

    references["reference_length"] = references["reference_length"].astype(int)

    if set(references["reference"]) != TARGET_REFS:
        raise ValueError("Target references do not match the required target set")

    if set(references["target_status"]) != {"NAS_target"}:
        raise ValueError("One or more focal references are not labelled NAS_target")

    return references


def load_target_timing(alignments):
    target_ids_by_run = {
        "full_meta": set(
            alignments.loc[
                alignments["dataset"] == "full_meta", "read_id"
            ].astype(str)
        ),
        "full_nas": set(
            alignments.loc[
                alignments["dataset"] == "full_nas_accepted", "read_id"
            ].astype(str)
        ),
    }

    target_start_times = {"full_meta": {}, "full_nas": {}}
    minimum_start_times = {"full_meta": math.inf, "full_nas": math.inf}

    for chunk in pd.read_csv(
        READ_METADATA_FILE,
        sep="\t",
        compression="gzip",
        usecols=["run", "read_id", "start_time", "decision_class"],
        dtype={
            "run": "string",
            "read_id": "string",
            "start_time": "float64",
            "decision_class": "string",
        },
        chunksize=500_000,
    ):
        chunk = chunk.loc[chunk["run"].isin(minimum_start_times)].copy()
        if chunk.empty:
            continue

        for run_name in minimum_start_times:
            run_chunk = chunk.loc[chunk["run"] == run_name]
            if run_chunk.empty:
                continue

            if run_chunk["start_time"].isna().any():
                raise ValueError(f"Missing start_time values in {run_name}")

            minimum_start_times[run_name] = min(
                minimum_start_times[run_name],
                float(run_chunk["start_time"].min()),
            )

            target_rows = run_chunk.loc[
                run_chunk["read_id"].isin(target_ids_by_run[run_name]),
                ["read_id", "start_time", "decision_class"],
            ]

            expected_decision = (
                "accepted" if run_name == "full_nas" else "not_applicable"
            )
            invalid = target_rows.loc[
                target_rows["decision_class"] != expected_decision
            ]
            if not invalid.empty:
                raise ValueError(
                    f"Target alignments in {run_name} joined to an unexpected "
                    "decision class"
                )

            for row in target_rows.itertuples(index=False):
                read_id = str(row.read_id)
                start_time = float(row.start_time)

                if read_id in target_start_times[run_name]:
                    if target_start_times[run_name][read_id] != start_time:
                        raise ValueError(
                            f"Conflicting start times for {run_name}, {read_id}"
                        )
                else:
                    target_start_times[run_name][read_id] = start_time

    for run_name in minimum_start_times:
        if not math.isfinite(minimum_start_times[run_name]):
            raise ValueError(f"No minimum start time found for {run_name}")

        expected_ids = target_ids_by_run[run_name]
        observed_ids = set(target_start_times[run_name])
        missing_ids = expected_ids - observed_ids
        unexpected_ids = observed_ids - expected_ids

        if missing_ids:
            raise ValueError(
                f"{len(missing_ids)} focal read IDs lack timing metadata in "
                f"{run_name}"
            )
        if unexpected_ids:
            raise ValueError(f"Unexpected timing IDs found in {run_name}")

    timing_records = []

    for dataset, config in DATASET_CONFIG.items():
        run_name = config["metadata_run"]
        run_minimum = minimum_start_times[run_name]
        dataset_alignments = alignments.loc[alignments["dataset"] == dataset]

        for row in dataset_alignments.itertuples(index=False):
            read_id = str(row.read_id)
            start_time = target_start_times[run_name][read_id]
            elapsed_hours = (start_time - run_minimum) / 3600.0

            if elapsed_hours < 0:
                raise ValueError(
                    f"Negative elapsed time for {dataset}, {read_id}"
                )

            hour_bin = max(1, int(math.ceil(elapsed_hours)))

            timing_records.append(
                {
                    "dataset": dataset,
                    "read_id": read_id,
                    "run": run_name,
                    "run_minimum_start_time": run_minimum,
                    "start_time": start_time,
                    "elapsed_hours": elapsed_hours,
                    "hour_bin": hour_bin,
                    "within_24_hours": hour_bin <= MAX_HOURS,
                }
            )

    timing_table = pd.DataFrame(timing_records)

    if len(timing_table) != len(alignments):
        raise ValueError(
            "Timing-table row count differs from alignment-table row count"
        )

    return timing_table, minimum_start_times


def build_outputs(alignments, timing_table, references):
    timed_alignments = alignments.merge(
        timing_table[
            [
                "dataset",
                "read_id",
                "start_time",
                "elapsed_hours",
                "hour_bin",
                "within_24_hours",
            ]
        ],
        on=["dataset", "read_id"],
        how="left",
        validate="one_to_one",
    )

    if timed_alignments[
        ["start_time", "elapsed_hours", "hour_bin", "within_24_hours"]
    ].isna().any().any():
        raise ValueError("Missing timing values after alignment join")

    reference_lookup = references.set_index("reference").to_dict(orient="index")
    time_series_rows = []
    threshold_rows = []
    final_rows = []

    for dataset, config in DATASET_CONFIG.items():
        dataset_rows = timed_alignments.loc[
            timed_alignments["dataset"] == dataset
        ].copy()

        for reference in sorted(TARGET_REFS):
            metadata = reference_lookup[reference]
            reference_length = int(metadata["reference_length"])
            reference_rows = dataset_rows.loc[
                dataset_rows["reference"] == reference
            ].copy()
            reference_rows = reference_rows.sort_values(
                ["hour_bin", "start_time", "read_id"]
            )

            depth = np.zeros(reference_length, dtype=np.int32)
            cumulative_reads = 0
            cumulative_aligned_bases = 0
            reference_time_rows = []

            for hour in range(1, MAX_HOURS + 1):
                hour_rows = reference_rows.loc[
                    reference_rows["hour_bin"] == hour
                ]

                for alignment in hour_rows.itertuples(index=False):
                    blocks = cigar_reference_blocks(
                        alignment.alignment_start_1based,
                        alignment.cigar,
                    )
                    block_aligned_bases = 0

                    for block_start, block_end in blocks:
                        if block_start < 0 or block_end > reference_length:
                            raise ValueError(
                                f"Reference block outside bounds for "
                                f"{alignment.read_id}"
                            )
                        if block_end <= block_start:
                            raise ValueError(
                                f"Invalid reference block for {alignment.read_id}"
                            )

                        depth[block_start:block_end] += 1
                        block_aligned_bases += block_end - block_start

                    if block_aligned_bases != alignment.cigar_aligned_query_bases:
                        raise ValueError(
                            "CIGAR block sum differs from canonical aligned-query "
                            f"bases for {alignment.read_id}"
                        )

                    cumulative_reads += 1
                    cumulative_aligned_bases += int(
                        alignment.cigar_aligned_query_bases
                    )

                covered_bases = int(np.count_nonzero(depth > 0))
                total_depth = int(depth.sum())

                if total_depth != cumulative_aligned_bases:
                    raise ValueError(
                        "Per-position depth sum differs from cumulative aligned "
                        f"bases for {dataset}, {reference}, hour {hour}"
                    )

                output_row = {
                    "dataset": dataset,
                    "method": config["method"],
                    "reference": reference,
                    "species": metadata["species"],
                    "reference_length": reference_length,
                    "time_hours": hour,
                    "cumulative_reads": cumulative_reads,
                    "cumulative_cigar_aligned_query_bases": (
                        cumulative_aligned_bases
                    ),
                    "covered_reference_bases_1x": covered_bases,
                    "coverage_breadth_percent_1x": (
                        100.0 * covered_bases / reference_length
                    ),
                    "cumulative_mean_depth": total_depth / reference_length,
                }
                reference_time_rows.append(output_row)
                time_series_rows.append(output_row)

            final_row = reference_time_rows[-1]
            within_24 = int((reference_rows["hour_bin"] <= MAX_HOURS).sum())
            after_24 = int((reference_rows["hour_bin"] > MAX_HOURS).sum())

            if final_row["cumulative_reads"] != within_24:
                raise ValueError(
                    f"Final cumulative read count mismatch for {dataset}, "
                    f"{reference}"
                )

            final_rows.append(
                {
                    "dataset": dataset,
                    "method": config["method"],
                    "reference": reference,
                    "species": metadata["species"],
                    "reference_length": reference_length,
                    "final_hour": MAX_HOURS,
                    "final_cumulative_reads": final_row["cumulative_reads"],
                    "final_cumulative_cigar_aligned_query_bases": final_row[
                        "cumulative_cigar_aligned_query_bases"
                    ],
                    "final_covered_reference_bases_1x": final_row[
                        "covered_reference_bases_1x"
                    ],
                    "final_coverage_breadth_percent_1x": final_row[
                        "coverage_breadth_percent_1x"
                    ],
                    "final_cumulative_mean_depth": final_row[
                        "cumulative_mean_depth"
                    ],
                    "total_target_alignments_all_times": len(reference_rows),
                    "target_alignments_within_24_hours": within_24,
                    "target_alignments_after_24_hours": after_24,
                }
            )

            for threshold in BREADTH_THRESHOLDS:
                first_reached = next(
                    (
                        row["time_hours"]
                        for row in reference_time_rows
                        if row["coverage_breadth_percent_1x"] >= threshold
                    ),
                    None,
                )

                threshold_rows.append(
                    {
                        "dataset": dataset,
                        "method": config["method"],
                        "reference": reference,
                        "species": metadata["species"],
                        "breadth_threshold_percent": threshold,
                        "first_hour_reached": first_reached,
                        "reached_within_24_hours": first_reached is not None,
                    }
                )

    return (
        pd.DataFrame(time_series_rows),
        pd.DataFrame(threshold_rows),
        pd.DataFrame(final_rows),
        timed_alignments,
    )


def main():
    started = datetime.now().astimezone().isoformat()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    for required_path in [
        ALIGNMENT_FILE,
        READ_METADATA_FILE,
        REFERENCE_FILE,
    ]:
        if not required_path.is_file():
            raise FileNotFoundError(required_path)

    alignments = load_target_alignments()
    references = load_reference_metadata()
    timing_table, minimum_start_times = load_target_timing(alignments)

    (
        time_series,
        threshold_table,
        final_summary,
        timed_alignments,
    ) = build_outputs(alignments, timing_table, references)

    if len(time_series) != 144:
        raise ValueError(
            f"Expected 144 time-series rows, observed {len(time_series)}"
        )
    if len(threshold_table) != 30:
        raise ValueError(
            f"Expected 30 threshold rows, observed {len(threshold_table)}"
        )
    if len(final_summary) != 6:
        raise ValueError(
            f"Expected 6 final-summary rows, observed {len(final_summary)}"
        )

    join_validation = (
        timed_alignments.groupby("dataset", observed=True)
        .agg(
            target_alignments=("read_id", "size"),
            target_alignments_with_timing=("start_time", "count"),
            target_alignments_within_24_hours=("within_24_hours", "sum"),
        )
        .reset_index()
    )
    join_validation["target_alignments_after_24_hours"] = (
        join_validation["target_alignments"]
        - join_validation["target_alignments_within_24_hours"]
    )
    join_validation["metadata_run"] = join_validation["dataset"].map(
        lambda dataset: DATASET_CONFIG[dataset]["metadata_run"]
    )
    join_validation["run_minimum_start_time"] = join_validation[
        "metadata_run"
    ].map(minimum_start_times)
    join_validation["timing_join_passed"] = (
        join_validation["target_alignments"]
        == join_validation["target_alignments_with_timing"]
    )

    if not join_validation["timing_join_passed"].all():
        raise ValueError("One or more timing joins failed")

    time_series.to_csv(TIME_SERIES_OUTPUT, sep="\t", index=False)
    threshold_table.to_csv(THRESHOLD_OUTPUT, sep="\t", index=False)
    final_summary.to_csv(FINAL_SUMMARY_OUTPUT, sep="\t", index=False)
    join_validation.to_csv(JOIN_VALIDATION_OUTPUT, sep="\t", index=False)

    finished = datetime.now().astimezone().isoformat()
    provenance = {
        "analysis_name": "target_coverage_over_time",
        "started": started,
        "finished": finished,
        "inputs": {
            "canonical_alignments": str(ALIGNMENT_FILE),
            "canonical_read_metadata": str(READ_METADATA_FILE),
            "reference_manifest": str(REFERENCE_FILE),
        },
        "datasets": {
            "Metagenomic": "full_meta",
            "NAS": "full_nas_accepted",
        },
        "targets": sorted(TARGET_REFS),
        "maximum_hours": MAX_HOURS,
        "hour_bin_definition": (
            "max(1, ceiling((start_time - minimum run start_time) / 3600))"
        ),
        "coverage_definition": {
            "covered_operations": ["M", "=", "X"],
            "reference_advance_without_coverage": ["D", "N"],
            "no_reference_advance": ["I", "S", "H", "P"],
            "breadth": "reference positions with depth >= 1 / reference length",
            "mean_depth": "sum of per-position depth / reference length",
            "reference_topology": "linear",
        },
        "breadth_thresholds_percent": BREADTH_THRESHOLDS,
        "expected_target_alignments": EXPECTED_TOTAL_TARGET_ALIGNMENTS,
        "time_series_rows": len(time_series),
        "threshold_rows": len(threshold_table),
        "final_summary_rows": len(final_summary),
        "validation_passed": True,
        "outputs": {
            "time_series": str(TIME_SERIES_OUTPUT),
            "threshold_times": str(THRESHOLD_OUTPUT),
            "final_summary": str(FINAL_SUMMARY_OUTPUT),
            "join_validation": str(JOIN_VALIDATION_OUTPUT),
        },
        "software": {
            "python": sys.version,
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
    }

    with open(PROVENANCE_OUTPUT, "w") as handle:
        json.dump(provenance, handle, indent=2)

    with open(LOG_OUTPUT, "w") as handle:
        handle.write(f"Started: {started}\n")
        handle.write(f"Target alignments: {len(alignments)}\n")
        handle.write(f"Targets: {', '.join(sorted(TARGET_REFS))}\n")
        handle.write(f"Maximum hours: {MAX_HOURS}\n")
        handle.write("Timing column: start_time\n")
        handle.write(
            "Hour bin: max(1, ceiling((start_time - minimum run start_time) / 3600))\n"
        )
        handle.write("Covered CIGAR operations: M, =, X\n")
        handle.write("Reference topology: linear\n")
        handle.write(f"Time-series rows: {len(time_series)}\n")
        handle.write(f"Threshold rows: {len(threshold_table)}\n")
        handle.write(f"Final-summary rows: {len(final_summary)}\n")
        handle.write(f"Time series: {TIME_SERIES_OUTPUT}\n")
        handle.write(f"Threshold times: {THRESHOLD_OUTPUT}\n")
        handle.write(f"Final summary: {FINAL_SUMMARY_OUTPUT}\n")
        handle.write(f"Join validation: {JOIN_VALIDATION_OUTPUT}\n")
        handle.write(f"Finished: {finished}\n")

    print(f"Wrote: {TIME_SERIES_OUTPUT}")
    print(f"Wrote: {THRESHOLD_OUTPUT}")
    print(f"Wrote: {FINAL_SUMMARY_OUTPUT}")
    print(f"Wrote: {JOIN_VALIDATION_OUTPUT}")
    print(f"Wrote: {PROVENANCE_OUTPUT}")
    print(f"Target alignments: {len(alignments)}")
    print(f"Time-series rows: {len(time_series)}")
    print(f"Threshold rows: {len(threshold_table)}")
    print(f"Final-summary rows: {len(final_summary)}")


if __name__ == "__main__":
    main()
