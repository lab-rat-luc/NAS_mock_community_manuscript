#!/usr/bin/env python3

import gzip
import json
import platform
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
    / "accept_reject_length"
)

ALIGNMENT_FILE = (
    PROJECT_DIR
    / "canonical_data"
    / "outputs"
    / "canonical_mapped_alignments.tsv.gz"
)

REFERENCE_FILE = (
    PROJECT_DIR
    / "canonical_data"
    / "inputs"
    / "reference_manifest.tsv"
)

OUTPUT_DIR = ANALYSIS_DIR / "outputs"
LOG_DIR = ANALYSIS_DIR / "logs"

READ_OUTPUT = (
    OUTPUT_DIR
    / "full_nas_accept_reject_read_lengths.tsv.gz"
)

DECISION_SUMMARY_OUTPUT = (
    OUTPUT_DIR
    / "full_nas_accept_reject_summary.tsv"
)

LENGTH_SUMMARY_OUTPUT = (
    OUTPUT_DIR
    / "full_nas_accept_reject_length_summary.tsv"
)

PROVENANCE_OUTPUT = (
    LOG_DIR
    / "accept_reject_length_provenance.json"
)

LOG_OUTPUT = (
    LOG_DIR
    / "accept_reject_length_run.log"
)

EXPECTED_DATASETS = {
    "full_nas_accepted": "Accepted",
    "full_nas_rejected": "Rejected",
}

# Reanalysis adapter-trimmed canonical dataset counts
EXPECTED_ROWS = {
    "full_nas_accepted": 76833,
    "full_nas_rejected": 87054,
}

EXPECTED_TOTAL_ROWS = sum(EXPECTED_ROWS.values())


def summarise_lengths(values):
    values = np.asarray(values, dtype=float)

    return {
        "n": int(values.size),
        "minimum_read_length_bases": float(
            np.min(values)
        ),
        "first_quartile_read_length_bases": float(
            np.quantile(values, 0.25)
        ),
        "median_read_length_bases": float(
            np.median(values)
        ),
        "mean_read_length_bases": float(
            np.mean(values)
        ),
        "third_quartile_read_length_bases": float(
            np.quantile(values, 0.75)
        ),
        "maximum_read_length_bases": float(
            np.max(values)
        ),
        "standard_deviation_read_length_bases": float(
            np.std(values, ddof=1)
        )
        if values.size > 1
        else 0.0,
        "interquartile_range_bases": float(
            np.quantile(values, 0.75)
            - np.quantile(values, 0.25)
        ),
    }


def main():
    started = datetime.now().astimezone().isoformat()

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    LOG_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not ALIGNMENT_FILE.is_file():
        raise FileNotFoundError(
            f"Missing canonical alignment file: "
            f"{ALIGNMENT_FILE}"
        )

    if not REFERENCE_FILE.is_file():
        raise FileNotFoundError(
            f"Missing reference manifest: "
            f"{REFERENCE_FILE}"
        )

    required_alignment_columns = [
        "dataset",
        "read_id",
        "reference",
        "mapped_read_sequence_bases",
        "cigar_aligned_query_bases",
        "soft_clipped_query_bases",
        "inserted_query_bases",
    ]

    alignments = pd.read_csv(
        ALIGNMENT_FILE,
        sep="\t",
        compression="gzip",
        usecols=required_alignment_columns,
        dtype={
            "dataset": "string",
            "read_id": "string",
            "reference": "string",
            "mapped_read_sequence_bases": "int64",
            "cigar_aligned_query_bases": "int64",
            "soft_clipped_query_bases": "int64",
            "inserted_query_bases": "int64",
        },
    )

    alignments = alignments.loc[
        alignments["dataset"].isin(
            EXPECTED_DATASETS
        )
    ].copy()

    observed_counts = (
        alignments.groupby("dataset")
        .size()
        .to_dict()
    )

    if observed_counts != EXPECTED_ROWS:
        raise ValueError(
            "Canonical dataset row counts differ from "
            f"expectation. Expected {EXPECTED_ROWS}, "
            f"observed {observed_counts}"
        )

    if len(alignments) != EXPECTED_TOTAL_ROWS:
        raise ValueError(
            f"Expected {EXPECTED_TOTAL_ROWS} records, "
            f"observed {len(alignments)}"
        )

    if (
        alignments["mapped_read_sequence_bases"]
        <= 0
    ).any():
        raise ValueError(
            "Non-positive complete SEQ lengths detected"
        )

    if (
        alignments["cigar_aligned_query_bases"]
        <= 0
    ).any():
        raise ValueError(
            "Non-positive CIGAR-aligned query lengths "
            "detected"
        )

    duplicate_keys = alignments.duplicated(
        subset=[
            "dataset",
            "read_id",
        ],
        keep=False,
    )

    if duplicate_keys.any():
        raise ValueError(
            "Duplicate dataset and read_id records "
            "detected"
        )

    references = pd.read_csv(
        REFERENCE_FILE,
        sep="\t",
        dtype="string",
        usecols=[
            "reference",
            "species",
            "refseq_accession",
            "target_status",
        ],
    )

    if len(references) != 16:
        raise ValueError(
            f"Expected 16 reference records, "
            f"observed {len(references)}"
        )

    if references["reference"].duplicated().any():
        raise ValueError(
            "Duplicate references detected in manifest"
        )

    alignments = alignments.merge(
        references,
        on="reference",
        how="left",
        validate="many_to_one",
    )

    if alignments[
        [
            "species",
            "refseq_accession",
            "target_status",
        ]
    ].isna().any().any():
        raise ValueError(
            "One or more alignment records did not join "
            "to the reference manifest"
        )

    alignments["decision_class"] = (
        alignments["dataset"]
        .map(EXPECTED_DATASETS)
    )

    alignments["read_length_bases"] = (
        alignments[
            "mapped_read_sequence_bases"
        ]
    )

    alignments["aligned_fraction_of_read"] = (
        alignments[
            "cigar_aligned_query_bases"
        ]
        / alignments[
            "mapped_read_sequence_bases"
        ]
    )

    if (
        alignments["aligned_fraction_of_read"]
        <= 0
    ).any():
        raise ValueError(
            "Non-positive aligned fractions detected"
        )

    if (
        alignments["aligned_fraction_of_read"]
        > 1
    ).any():
        raise ValueError(
            "Aligned fraction greater than one detected"
        )

    output_columns = [
        "dataset",
        "decision_class",
        "read_id",
        "reference",
        "refseq_accession",
        "species",
        "target_status",
        "read_length_bases",
        "cigar_aligned_query_bases",
        "soft_clipped_query_bases",
        "inserted_query_bases",
        "aligned_fraction_of_read",
    ]

    alignments = alignments[
        output_columns
    ].sort_values(
        [
            "reference",
            "decision_class",
            "read_id",
        ]
    )

    alignments.to_csv(
        READ_OUTPUT,
        sep="\t",
        index=False,
        compression={
            "method": "gzip",
            "compresslevel": 6,
            "mtime": 0,
        },
    )

    read_counts = (
        alignments.groupby(
            [
                "reference",
                "refseq_accession",
                "species",
                "target_status",
                "decision_class",
            ],
            observed=True,
        )
        .size()
        .unstack(
            fill_value=0
        )
        .reset_index()
    )

    for decision in [
        "Accepted",
        "Rejected",
    ]:
        if decision not in read_counts.columns:
            read_counts[decision] = 0

    read_counts = read_counts.rename(
        columns={
            "Accepted": "accepted_reads",
            "Rejected": "rejected_reads",
        }
    )

    read_counts["total_reads"] = (
        read_counts["accepted_reads"]
        + read_counts["rejected_reads"]
    )

    read_counts["accepted_percent"] = (
        100.0
        * read_counts["accepted_reads"]
        / read_counts["total_reads"]
    )

    read_counts["rejected_percent"] = (
        100.0
        * read_counts["rejected_reads"]
        / read_counts["total_reads"]
    )

    if not np.allclose(
        read_counts["accepted_percent"]
        + read_counts["rejected_percent"],
        100.0,
        rtol=0,
        atol=1e-9,
    ):
        raise ValueError(
            "Accepted and rejected percentages do not "
            "sum to 100"
        )

    read_counts = read_counts.sort_values(
        [
            "accepted_percent",
            "total_reads",
        ],
        ascending=[
            False,
            False,
        ],
    )

    read_counts.to_csv(
        DECISION_SUMMARY_OUTPUT,
        sep="\t",
        index=False,
    )

    summary_rows = []

    grouped = alignments.groupby(
        [
            "reference",
            "refseq_accession",
            "species",
            "target_status",
            "decision_class",
        ],
        observed=True,
        sort=True,
    )

    for group_keys, group_data in grouped:
        (
            reference,
            accession,
            species,
            target_status,
            decision_class,
        ) = group_keys

        statistics = summarise_lengths(
            group_data[
                "read_length_bases"
            ]
        )

        summary_rows.append(
            {
                "reference": reference,
                "refseq_accession": accession,
                "species": species,
                "target_status": target_status,
                "decision_class": decision_class,
                **statistics,
            }
        )

    length_summary = pd.DataFrame(
        summary_rows
    )

    if len(length_summary) != 32:
        raise ValueError(
            "Expected 32 species and decision-class "
            f"summary records, observed "
            f"{len(length_summary)}"
        )

    length_summary = length_summary.sort_values(
        [
            "reference",
            "decision_class",
        ]
    )

    length_summary.to_csv(
        LENGTH_SUMMARY_OUTPUT,
        sep="\t",
        index=False,
    )

    finished = datetime.now().astimezone().isoformat()

    provenance = {
        "analysis_name": (
            "full_nas_accept_reject_read_length"
        ),
        "started": started,
        "finished": finished,
        "canonical_alignment_input": str(
            ALIGNMENT_FILE
        ),
        "reference_manifest": str(
           REFERENCE_FILE
        ),
        "datasets": {
            "accepted": "full_nas_accepted",
            "rejected": "full_nas_rejected",
        },
        "records": {
            "accepted": int(
                observed_counts[
                    "full_nas_accepted"
                ]
            ),
            "rejected": int(
                observed_counts[
                    "full_nas_rejected"
                ]
            ),
            "total": int(
                len(alignments)
            ),
        },
        "definitions": {
            "read_length": (
                "Complete SAM SEQ length, stored as "
                "mapped_read_sequence_bases"
            ),
            "aligned_query_bases": (
                "CIGAR M plus = plus X query bases"
            ),
            "aligned_fraction_of_read": (
                "CIGAR-aligned query bases divided by "
                "complete SAM SEQ length"
            ),
            "sampling": (
                "No sampling. All canonical records used."
            ),
        },
        "outputs": {
            "read_level_table": str(
                READ_OUTPUT
            ),
            "decision_summary": str(
                DECISION_SUMMARY_OUTPUT
            ),
            "length_summary": str(
                LENGTH_SUMMARY_OUTPUT
            ),
        },
        "software": {
            "python": sys.version,
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "validation_passed": True,
    }

    with open(
        PROVENANCE_OUTPUT,
        "w",
    ) as handle:
        json.dump(
            provenance,
            handle,
            indent=2,
        )

    with open(
        LOG_OUTPUT,
        "w",
    ) as handle:
        handle.write(
            f"Started: {started}\n"
        )
        handle.write(
            f"Accepted reads: "
            f"{observed_counts['full_nas_accepted']}\n"
        )
        handle.write(
            f"Rejected reads: "
            f"{observed_counts['full_nas_rejected']}\n"
        )
        handle.write(
            f"Total reads: {len(alignments)}\n"
        )
        handle.write(
            "Read-length definition: complete SAM SEQ "
            "length\n"
        )
        handle.write(
            "Sampling: none; all records retained\n"
        )
        handle.write(
            f"Read-level output: {READ_OUTPUT}\n"
        )
        handle.write(
            f"Decision summary: "
            f"{DECISION_SUMMARY_OUTPUT}\n"
        )
        handle.write(
            f"Length summary: "
            f"{LENGTH_SUMMARY_OUTPUT}\n"
        )
        handle.write(
            f"Finished: {finished}\n"
        )

    print(
        f"Wrote: {READ_OUTPUT}"
    )
    print(
        f"Wrote: {DECISION_SUMMARY_OUTPUT}"
    )
    print(
        f"Wrote: {LENGTH_SUMMARY_OUTPUT}"
    )
    print(
        f"Wrote: {PROVENANCE_OUTPUT}"
    )
    print(
        f"Accepted reads: "
        f"{observed_counts['full_nas_accepted']}"
    )
    print(
        f"Rejected reads: "
        f"{observed_counts['full_nas_rejected']}"
    )
    print(
        f"Total reads: {len(alignments)}"
    )


if __name__ == "__main__":
    main()
