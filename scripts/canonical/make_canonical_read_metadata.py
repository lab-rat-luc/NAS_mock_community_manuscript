#!/usr/bin/env python3

import argparse
import csv
import gzip
import hashlib
import json
import platform
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path


parser = argparse.ArgumentParser(
    description="Generate canonical read metadata."
)
parser.add_argument("--project-dir", required=True)
args = parser.parse_args()
PROJECT_DIR = Path(args.project_dir).resolve()

CANONICAL_DIR = PROJECT_DIR / "canonical_data"

MANIFEST_PATH = (
    CANONICAL_DIR
    / "inputs"
    / "sequencing_summary_manifest.tsv"
)

OUTPUT_PATH = (
    CANONICAL_DIR
    / "outputs"
    / "canonical_read_metadata.tsv.gz"
)

LOG_PATH = (
    CANONICAL_DIR
    / "logs"
    / "canonical_read_metadata_run.log"
)

PROVENANCE_PATH = (
    CANONICAL_DIR
    / "logs"
    / "canonical_read_metadata_provenance.json"
)


EXPECTED_SOURCE_FIELDS = {
    "read_id",
    "run_id",
    "channel",
    "mux",
    "start_time",
    "duration",
    "passes_filtering",
    "template_start",
    "template_duration",
    "sequence_length_template",
    "mean_qscore_template",
    "pore_type",
    "end_reason",
}


OUTPUT_FIELDS = [
    "run",
    "experiment",
    "read_id",
    "run_id",
    "channel",
    "mux",
    "start_time",
    "duration",
    "template_start",
    "template_duration",
    "sequence_length_template",
    "mean_qscore_template",
    "pore_type",
    "end_reason",
    "sequencing_method",
    "flowcell_region",
    "decision_class",
]


EXPECTED_ROWS_BY_RUN = {
    "full_meta": 8914639,
    "full_nas": 27514001,
    "half_flowcell": 13096684,
}


EXPECTED_TOTAL_ROWS = sum(EXPECTED_ROWS_BY_RUN.values())


def sha256_file(path):
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        while True:
            block = handle.read(1024 * 1024)

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


def classify_read(run_name, channel, end_reason):
    if run_name == "full_meta":
        return (
            "meta",
            "whole_flowcell",
            "not_applicable",
        )

    if run_name == "full_nas":
        if end_reason == "signal_positive":
            decision_class = "accepted"
        elif end_reason == "data_service_unblock_mux_change":
            decision_class = "rejected"
        else:
            decision_class = "other"

        return (
            "nas",
            "whole_flowcell",
            decision_class,
        )

    if run_name == "half_flowcell":
        if 1 <= channel <= 256:
            sequencing_method = "nas"
            flowcell_region = "channels_1_to_256"

            if end_reason == "signal_positive":
                decision_class = "accepted"
            elif end_reason == "data_service_unblock_mux_change":
                decision_class = "rejected"
            else:
                decision_class = "other"

            return (
                sequencing_method,
                flowcell_region,
                decision_class,
            )

        if 257 <= channel <= 512:
            return (
                "meta",
                "channels_257_to_512",
                "not_applicable",
            )

        raise ValueError(
            f"Channel outside expected range: {channel}"
        )

    raise ValueError(
        f"Unexpected run name in manifest: {run_name}"
    )


def read_manifest():
    required_fields = {
        "run",
        "experiment",
        "sequencing_summary",
    }

    with open(MANIFEST_PATH, newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        observed_fields = set(reader.fieldnames or [])

        missing_fields = required_fields - observed_fields

        if missing_fields:
            raise ValueError(
                "Manifest is missing required fields: "
                f"{sorted(missing_fields)}"
            )

        rows = list(reader)

    if len(rows) != 3:
        raise ValueError(
            "Expected exactly three sequencing-summary "
            f"sources, found {len(rows)}"
        )

    run_names = [row["run"] for row in rows]

    if len(set(run_names)) != len(run_names):
        raise ValueError(
            "Duplicate run labels found in sequencing "
            "summary manifest"
        )

    expected_runs = set(EXPECTED_ROWS_BY_RUN)
    observed_runs = set(run_names)

    if observed_runs != expected_runs:
        raise ValueError(
            "Unexpected run labels. "
            f"Expected {sorted(expected_runs)}, "
            f"observed {sorted(observed_runs)}"
        )

    return rows


def main():
    started = datetime.now().astimezone().isoformat()

    OUTPUT_PATH.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    LOG_PATH.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    manifest_rows = read_manifest()

    source_provenance = []
    total_rows = 0

    with open(LOG_PATH, "w") as log_handle:
        log_handle.write(
            f"Started: {started}\n"
        )

        log_handle.write(
            f"Manifest: {MANIFEST_PATH}\n"
        )

        log_handle.write(
            "Half-flowcell NAS region: "
            "channels 1 to 256\n"
        )

        log_handle.write(
            "Half-flowcell Meta region: "
            "channels 257 to 512\n"
        )

        log_handle.write(
            "Accepted NAS end reason: "
            "signal_positive\n"
        )

        log_handle.write(
            "Rejected NAS end reason: "
            "data_service_unblock_mux_change\n"
        )

        with gzip.open(
            OUTPUT_PATH,
            "wt",
            newline="",
        ) as output_handle:
            writer = csv.DictWriter(
                output_handle,
                delimiter="\t",
                fieldnames=OUTPUT_FIELDS,
                lineterminator="\n",
            )

            writer.writeheader()

            for manifest_row in manifest_rows:
                run_name = manifest_row["run"]
                experiment = manifest_row["experiment"]

                source_path = Path(
                    manifest_row["sequencing_summary"]
                )

                if not source_path.is_file():
                    raise FileNotFoundError(
                        f"Missing source file: {source_path}"
                    )

                log_handle.write(
                    f"[SOURCE] {source_path}\n"
                )

                row_count = 0
                pass_status_counts = Counter()
                end_reason_counts = Counter()
                decision_class_counts = Counter()
                sequencing_method_counts = Counter()
                flowcell_region_counts = Counter()
                channel_counts = Counter()

                with open(
                    source_path,
                    newline="",
                ) as source_handle:
                    reader = csv.DictReader(
                        source_handle,
                        delimiter="\t",
                    )

                    observed_fields = set(
                        reader.fieldnames or []
                    )

                    missing_fields = (
                        EXPECTED_SOURCE_FIELDS
                        - observed_fields
                    )

                    if missing_fields:
                        raise ValueError(
                            f"Missing fields in {source_path}: "
                            f"{sorted(missing_fields)}"
                        )

                    for row in reader:
                        read_id = row["read_id"]
                        pass_status = row[
                            "passes_filtering"
                        ]
                        end_reason = row["end_reason"]

                        if not read_id:
                            raise ValueError(
                                "Encountered a row with an "
                                f"empty read_id in {source_path}"
                            )

                        if pass_status != "TRUE":
                            raise ValueError(
                                "Non-pass read encountered in "
                                f"{source_path}: {read_id}"
                            )

                        try:
                            channel = int(row["channel"])
                        except ValueError as error:
                            raise ValueError(
                                "Invalid channel value for "
                                f"{read_id}: {row['channel']}"
                            ) from error

                        (
                            sequencing_method,
                            flowcell_region,
                            decision_class,
                        ) = classify_read(
                            run_name,
                            channel,
                            end_reason,
                        )

                        writer.writerow(
                            {
                                "run": run_name,
                                "experiment": experiment,
                                "read_id": read_id,
                                "run_id": row["run_id"],
                                "channel": row["channel"],
                                "mux": row["mux"],
                                "start_time": row[
                                    "start_time"
                                ],
                                "duration": row["duration"],
                                "template_start": row[
                                    "template_start"
                                ],
                                "template_duration": row[
                                    "template_duration"
                                ],
                                "sequence_length_template": (
                                    row[
                                        "sequence_length_template"
                                    ]
                                ),
                                "mean_qscore_template": row[
                                    "mean_qscore_template"
                                ],
                                "pore_type": row["pore_type"],
                                "end_reason": end_reason,
                                "sequencing_method": (
                                    sequencing_method
                                ),
                                "flowcell_region": (
                                    flowcell_region
                                ),
                                "decision_class": (
                                    decision_class
                                ),
                            }
                        )

                        row_count += 1
                        total_rows += 1

                        pass_status_counts[
                            pass_status
                        ] += 1

                        end_reason_counts[
                            end_reason
                        ] += 1

                        decision_class_counts[
                            decision_class
                        ] += 1

                        sequencing_method_counts[
                            sequencing_method
                        ] += 1

                        flowcell_region_counts[
                            flowcell_region
                        ] += 1

                        channel_counts[
                            str(channel)
                        ] += 1

                expected_rows = EXPECTED_ROWS_BY_RUN[
                    run_name
                ]

                if row_count != expected_rows:
                    raise ValueError(
                        f"Row-count validation failed for "
                        f"{run_name}: expected "
                        f"{expected_rows}, observed "
                        f"{row_count}"
                    )

                observed_channels = [
                    int(channel)
                    for channel in channel_counts
                ]

                source_provenance.append(
                    {
                        "run": run_name,
                        "experiment": experiment,
                        "source": str(source_path),
                        "source_sha256": sha256_file(
                            source_path
                        ),
                        "rows": row_count,
                        "expected_rows": expected_rows,
                        "row_count_passed": True,
                        "minimum_channel": min(
                            observed_channels
                        ),
                        "maximum_channel": max(
                            observed_channels
                        ),
                        "pass_status_counts": dict(
                            pass_status_counts
                        ),
                        "end_reason_counts": dict(
                            end_reason_counts
                        ),
                        "decision_class_counts": dict(
                            decision_class_counts
                        ),
                        "sequencing_method_counts": dict(
                            sequencing_method_counts
                        ),
                        "flowcell_region_counts": dict(
                            flowcell_region_counts
                        ),
                    }
                )

                log_handle.write(
                    f"[OK] {run_name}: "
                    f"rows={row_count}\n"
                )

    if total_rows != EXPECTED_TOTAL_ROWS:
        raise ValueError(
            "Total row-count validation failed: "
            f"expected {EXPECTED_TOTAL_ROWS}, "
            f"observed {total_rows}"
        )

    finished = datetime.now().astimezone().isoformat()

    provenance = {
        "analysis_name": (
            "canonical_read_metadata"
        ),
        "started": started,
        "finished": finished,
        "project_directory": str(PROJECT_DIR),
        "manifest": str(MANIFEST_PATH),
        "manifest_sha256": sha256_file(
            MANIFEST_PATH
        ),
        "output": str(OUTPUT_PATH),
        "output_sha256": sha256_file(
            OUTPUT_PATH
        ),
        "total_rows": total_rows,
        "expected_total_rows": (
            EXPECTED_TOTAL_ROWS
        ),
        "validation_passed": True,
        "classification_rules": {
            "half_nas_channels": "1-256",
            "half_meta_channels": "257-512",
            "accepted_end_reason": (
                "signal_positive"
            ),
            "rejected_end_reason": (
                "data_service_unblock_mux_change"
            ),
            "other_end_reasons": (
                "Preserved verbatim and classified "
                "as other for NAS reads"
            ),
            "meta_decision_class": (
                "not_applicable"
            ),
        },
        "sources": source_provenance,
        "python_version": sys.version,
        "platform": platform.platform(),
        "log": str(LOG_PATH),
    }

    with open(
        PROVENANCE_PATH,
        "w",
    ) as provenance_handle:
        json.dump(
            provenance,
            provenance_handle,
            indent=2,
        )

    with open(LOG_PATH, "a") as log_handle:
        log_handle.write(
            f"Finished: {finished}\n"
        )

        log_handle.write(
            f"[OUTPUT] {OUTPUT_PATH}\n"
        )

        log_handle.write(
            f"[ROWS] {total_rows}\n"
        )

        log_handle.write(
            f"[PROVENANCE] {PROVENANCE_PATH}\n"
        )

    print(f"Wrote: {OUTPUT_PATH}")
    print(f"Wrote: {PROVENANCE_PATH}")
    print(f"Wrote: {LOG_PATH}")
    print(f"Rows: {total_rows}")


if __name__ == "__main__":
    main()
