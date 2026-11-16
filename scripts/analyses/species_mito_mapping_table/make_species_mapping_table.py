#!/usr/bin/env python3

import csv
import hashlib
import json
import platform
import sys
from datetime import datetime
from pathlib import Path
import argparse


_parser = argparse.ArgumentParser()
_parser.add_argument("--project-dir", required=True)
_args = _parser.parse_args()
PROJECT_DIR = Path(_args.project_dir).resolve()

ANALYSIS_DIR = (
    PROJECT_DIR
    / "canonical_data"
    / "analyses"
    / "species_mito_mapping_table"
)

CANONICAL_SUMMARY = (
    PROJECT_DIR
    / "canonical_data"
    / "outputs"
    / "canonical_mito_mapping_summary.tsv"
)

REFERENCE_MANIFEST = (
    PROJECT_DIR
    / "canonical_data"
    / "inputs"
    / "reference_manifest.tsv"
)

OUTPUT_PATH = (
    ANALYSIS_DIR
    / "outputs"
    / "species_mito_mapping_cov_ra_nas_meta_table.tsv"
)

PROVENANCE_PATH = (
    ANALYSIS_DIR
    / "logs"
    / "species_mito_mapping_table_provenance.json"
)

LOG_PATH = (
    ANALYSIS_DIR
    / "logs"
    / "species_mito_mapping_table_run.log"
)


OUTPUT_FIELDS = [
    "reference",
    "species",
    "refseq_accession",
    "order",
    "target_status",
    "expected_DNA_RA_percent",
    "reference_length",
    "meta_mitochondrial_RA_percent",
    "meta_mapped_reads",
    "meta_cigar_aligned_query_bases",
    "meta_mean_reference_coverage_depth",
    "meta_mean_mapped_read_length",
    "nas_mitochondrial_RA_percent",
    "nas_mapped_reads",
    "nas_cigar_aligned_query_bases",
    "nas_mean_reference_coverage_depth",
    "nas_mean_mapped_read_length",
]


REQUIRED_CANONICAL_FIELDS = {
    "dataset",
    "reference",
    "reference_length",
    "mapped_reads",
    "mapped_read_sequence_bases",
    "cigar_aligned_query_bases",
    "mean_mapped_read_length",
    "mitochondrial_RA_fraction",
    "mitochondrial_RA_percent",
}


REQUIRED_REFERENCE_FIELDS = {
    "reference",
    "refseq_accession",
    "species",
    "order",
    "expected_DNA_RA_percent",
    "target_status",
    "reference_length",
}


def sha256_file(path):
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        while True:
            block = handle.read(1024 * 1024)

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


def read_tsv(path, required_fields):
    with open(path, newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        observed_fields = set(reader.fieldnames or [])

        missing_fields = required_fields - observed_fields

        if missing_fields:
            raise ValueError(
                f"Missing required fields in {path}: "
                f"{sorted(missing_fields)}"
            )

        rows = list(reader)

    return rows


def calculate_mean_depth(
    aligned_query_bases,
    reference_length,
):
    aligned_query_bases = int(aligned_query_bases)
    reference_length = int(reference_length)

    if reference_length <= 0:
        raise ValueError(
            "Reference length must be greater than zero"
        )

    return aligned_query_bases / reference_length


def main():
    started = datetime.now().astimezone().isoformat()

    if not CANONICAL_SUMMARY.is_file():
        raise FileNotFoundError(
            f"Missing canonical summary: "
            f"{CANONICAL_SUMMARY}"
        )

    if not REFERENCE_MANIFEST.is_file():
        raise FileNotFoundError(
            f"Missing reference manifest: "
            f"{REFERENCE_MANIFEST}"
        )

    canonical_rows = read_tsv(
        CANONICAL_SUMMARY,
        REQUIRED_CANONICAL_FIELDS,
    )

    reference_rows = read_tsv(
        REFERENCE_MANIFEST,
        REQUIRED_REFERENCE_FIELDS,
    )

    if len(reference_rows) != 16:
        raise ValueError(
            "Expected 16 references in the reference "
            f"manifest, found {len(reference_rows)}"
        )

    selected_datasets = {
        "full_meta",
        "full_nas_accepted",
    }

    selected_rows = [
        row
        for row in canonical_rows
        if row["dataset"] in selected_datasets
    ]

    if len(selected_rows) != 32:
        raise ValueError(
            "Expected 32 canonical rows for full_meta "
            "and full_nas_accepted, found "
            f"{len(selected_rows)}"
        )

    canonical_lookup = {}

    for row in selected_rows:
        key = (
            row["dataset"],
            row["reference"],
        )

        if key in canonical_lookup:
            raise ValueError(
                f"Duplicate canonical row: {key}"
            )

        canonical_lookup[key] = row

    reference_lookup = {}

    for row in reference_rows:
        reference = row["reference"]

        if reference in reference_lookup:
            raise ValueError(
                "Duplicate reference in manifest: "
                f"{reference}"
            )

        reference_lookup[reference] = row

    canonical_references = {
        row["reference"]
        for row in selected_rows
    }

    manifest_references = set(reference_lookup)

    if canonical_references != manifest_references:
        missing_from_canonical = sorted(
            manifest_references - canonical_references
        )

        missing_from_manifest = sorted(
            canonical_references - manifest_references
        )

        raise ValueError(
            "Reference mismatch. "
            f"Missing from canonical: "
            f"{missing_from_canonical}. "
            f"Missing from manifest: "
            f"{missing_from_manifest}"
        )

    output_rows = []

    for reference, metadata in reference_lookup.items():
        meta_key = (
            "full_meta",
            reference,
        )

        nas_key = (
            "full_nas_accepted",
            reference,
        )

        if meta_key not in canonical_lookup:
            raise ValueError(
                f"Missing full_meta row for {reference}"
            )

        if nas_key not in canonical_lookup:
            raise ValueError(
                "Missing full_nas_accepted row for "
                f"{reference}"
            )

        meta = canonical_lookup[meta_key]
        nas = canonical_lookup[nas_key]

        manifest_length = int(
            metadata["reference_length"]
        )

        meta_length = int(
            meta["reference_length"]
        )

        nas_length = int(
            nas["reference_length"]
        )

        if manifest_length != meta_length:
            raise ValueError(
                "Meta reference-length mismatch for "
                f"{reference}"
            )

        if manifest_length != nas_length:
            raise ValueError(
                "NAS reference-length mismatch for "
                f"{reference}"
            )

        meta_mean_depth = calculate_mean_depth(
            meta["cigar_aligned_query_bases"],
            manifest_length,
        )

        nas_mean_depth = calculate_mean_depth(
            nas["cigar_aligned_query_bases"],
            manifest_length,
        )

        output_rows.append(
            {
                "reference": reference,
                "species": metadata["species"],
                "refseq_accession": metadata[
                    "refseq_accession"
                ],
                "order": metadata["order"],
                "target_status": metadata[
                    "target_status"
                ],
                "expected_DNA_RA_percent": metadata[
                    "expected_DNA_RA_percent"
                ],
                "reference_length": manifest_length,
                "meta_mitochondrial_RA_percent": (
                    meta["mitochondrial_RA_percent"]
                ),
                "meta_mapped_reads": meta[
                    "mapped_reads"
                ],
                "meta_cigar_aligned_query_bases": (
                    meta["cigar_aligned_query_bases"]
                ),
                "meta_mean_reference_coverage_depth": (
                    f"{meta_mean_depth:.10f}"
                ),
                "meta_mean_mapped_read_length": (
                    meta["mean_mapped_read_length"]
                ),
                "nas_mitochondrial_RA_percent": (
                    nas["mitochondrial_RA_percent"]
                ),
                "nas_mapped_reads": nas[
                    "mapped_reads"
                ],
                "nas_cigar_aligned_query_bases": (
                    nas["cigar_aligned_query_bases"]
                ),
                "nas_mean_reference_coverage_depth": (
                    f"{nas_mean_depth:.10f}"
                ),
                "nas_mean_mapped_read_length": (
                    nas["mean_mapped_read_length"]
                ),
            }
        )

    output_rows.sort(
        key=lambda row: float(
            row["meta_mitochondrial_RA_percent"]
        ),
        reverse=True,
    )

    meta_ra_sum = sum(
        float(
            row["meta_mitochondrial_RA_percent"]
        )
        for row in output_rows
    )

    nas_ra_sum = sum(
        float(
            row["nas_mitochondrial_RA_percent"]
        )
        for row in output_rows
    )

    tolerance = 0.000001

    if abs(meta_ra_sum - 100.0) > tolerance:
        raise ValueError(
            f"Meta RA sum validation failed: "
            f"{meta_ra_sum}"
        )

    if abs(nas_ra_sum - 100.0) > tolerance:
        raise ValueError(
            f"NAS RA sum validation failed: "
            f"{nas_ra_sum}"
        )

    target_count = sum(
        row["target_status"] == "NAS_target"
        for row in output_rows
    )

    related_non_target_count = sum(
        row["target_status"]
        == "related_non_target"
        for row in output_rows
    )

    other_non_target_count = sum(
        row["target_status"]
        == "other_non_target"
        for row in output_rows
    )

    if target_count != 3:
        raise ValueError(
            "Expected 3 NAS targets, found "
            f"{target_count}"
        )

    if related_non_target_count != 1:
        raise ValueError(
            "Expected 1 related non-target, found "
            f"{related_non_target_count}"
        )

    if other_non_target_count != 12:
        raise ValueError(
            "Expected 12 other non-targets, found "
            f"{other_non_target_count}"
        )

    OUTPUT_PATH.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    LOG_PATH.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with open(
        OUTPUT_PATH,
        "w",
        newline="",
    ) as output_handle:
        writer = csv.DictWriter(
            output_handle,
            delimiter="\t",
            fieldnames=OUTPUT_FIELDS,
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(output_rows)

    finished = datetime.now().astimezone().isoformat()

    provenance = {
        "analysis_name": (
            "species_mito_mapping_table"
        ),
        "started": started,
        "finished": finished,
        "analysis_directory": str(
            ANALYSIS_DIR
        ),
        "canonical_summary": str(
            CANONICAL_SUMMARY
        ),
        "canonical_summary_sha256": sha256_file(
            CANONICAL_SUMMARY
        ),
        "reference_manifest": str(
            REFERENCE_MANIFEST
        ),
        "reference_manifest_sha256": sha256_file(
            REFERENCE_MANIFEST
        ),
        "output": str(OUTPUT_PATH),
        "output_sha256": sha256_file(
            OUTPUT_PATH
        ),
        "datasets": {
            "meta": "full_meta",
            "nas": "full_nas_accepted",
        },
        "dataset_scope": {
            "meta": (
                "All retained metagenomic alignments"
            ),
            "nas": (
                "Accepted NAS alignments only"
            ),
        },
        "rows": len(output_rows),
        "target_count": target_count,
        "related_non_target_count": (
            related_non_target_count
        ),
        "other_non_target_count": (
            other_non_target_count
        ),
        "meta_RA_sum_percent": meta_ra_sum,
        "nas_RA_sum_percent": nas_ra_sum,
        "validation_passed": True,
        "metric_definitions": {
            "mitochondrial_RA": (
                "Species CIGAR-aligned query bases "
                "participating in M, = or X operations "
                "divided by total M, = and X query bases "
                "across all 16 mitochondrial references "
                "within the dataset"
            ),
            "mapped_reads": (
                "Retained primary mapped alignments"
            ),
            "mean_reference_coverage_depth": (
                "Total M, = and X aligned query bases "
                "divided by mitochondrial reference length"
            ),
            "mean_mapped_read_length": (
                "Complete SAM SEQ bases divided by "
                "mapped-read count"
            ),
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
        },
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

    with open(
        LOG_PATH,
        "w",
    ) as log_handle:
        log_handle.write(
            f"Started: {started}\n"
        )

        log_handle.write(
            f"Analysis directory: {ANALYSIS_DIR}\n"
        )

        log_handle.write(
            f"Canonical summary: {CANONICAL_SUMMARY}\n"
        )

        log_handle.write(
            f"Reference manifest: {REFERENCE_MANIFEST}\n"
        )

        log_handle.write(
            "Meta dataset: full_meta\n"
        )

        log_handle.write(
            "NAS dataset: full_nas_accepted\n"
        )

        log_handle.write(
            "Mitochondrial RA base definition: "
            "M + = + X aligned query bases\n"
        )

        log_handle.write(
            "Mean mapped-read length definition: "
            "complete SAM SEQ length\n"
        )

        log_handle.write(
            "Mean coverage-depth definition: "
            "M + = + X bases divided by reference length\n"
        )

        log_handle.write(
            f"Rows: {len(output_rows)}\n"
        )

        log_handle.write(
            f"Meta RA sum: {meta_ra_sum:.12f}\n"
        )

        log_handle.write(
            f"NAS RA sum: {nas_ra_sum:.12f}\n"
        )

        log_handle.write(
            f"NAS targets: {target_count}\n"
        )

        log_handle.write(
            "Related non-targets: "
            f"{related_non_target_count}\n"
        )

        log_handle.write(
            "Other non-targets: "
            f"{other_non_target_count}\n"
        )

        log_handle.write(
            f"Output: {OUTPUT_PATH}\n"
        )

        log_handle.write(
            f"Provenance: {PROVENANCE_PATH}\n"
        )

        log_handle.write(
            f"Finished: {finished}\n"
        )

    print(f"Wrote: {OUTPUT_PATH}")
    print(f"Wrote: {PROVENANCE_PATH}")
    print(f"Wrote: {LOG_PATH}")
    print(f"Rows: {len(output_rows)}")
    print(f"Meta RA sum: {meta_ra_sum:.10f}")
    print(f"NAS RA sum: {nas_ra_sum:.10f}")


if __name__ == "__main__":
    main()
