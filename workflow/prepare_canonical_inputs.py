#!/usr/bin/env python3
"""Prepare canonical manifests and report inputs; dry run by default."""
import argparse
import csv
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--workspace", required=True, type=Path)
parser.add_argument("--supporting-root", required=True, type=Path)
parser.add_argument("--execute", action="store_true")
args = parser.parse_args()

repo = Path(__file__).resolve().parents[1]
workspace = args.workspace.resolve()
supporting = args.supporting_root.resolve()

if not workspace.is_dir():
    parser.error("Workspace must already exist.")
if not supporting.is_dir():
    parser.error("Supporting-data root does not exist.")
for protected in (repo, supporting):
    if (
        workspace == protected
        or protected in workspace.parents
        or workspace in protected.parents
    ):
        parser.error("Workspace overlaps the repository or supporting data.")

names = [
    "dataset_manifest.tsv",
    "barcode_dataset_manifest.tsv",
    "sequencing_summary_manifest.tsv",
    "ont_report_manifest.tsv",
    "run_statistics_manifest.tsv",
]
tasks = []

for name in names:
    source = repo / "metadata/templates" / name
    text = source.read_text().replace("{WORKSPACE}", str(workspace))
    if "{WORKSPACE}" in text:
        parser.error("Unresolved workspace placeholder.")
    destination = workspace / "canonical_data/inputs" / name
    tasks.append((destination, text.encode()))

reference = repo / "metadata/reference_manifest.tsv"
tasks.append((
    workspace / "canonical_data/inputs/reference_manifest.tsv",
    reference.read_bytes()
))

report_template = repo / "metadata/templates/ont_report_manifest.tsv"
with report_template.open(newline="") as handle:
    reports = list(csv.DictReader(handle, delimiter="\t"))

if len(reports) != 3:
    parser.error("Expected three MinKNOW reports.")

for row in reports:
    value = row["report_json"]
    prefix = "{WORKSPACE}/"
    if not value.startswith(prefix):
        parser.error("Unexpected report path: " + value)
    relative = Path(value[len(prefix):])
    if relative.is_absolute() or ".." in relative.parts:
        parser.error("Unsafe report path.")
    source = supporting / relative
    if not source.is_file():
        parser.error("Missing supporting report: " + str(source))
    tasks.append((workspace / relative, source.read_bytes()))

destinations = set()
for destination, content in tasks:
    if not destination.resolve().is_relative_to(workspace):
        parser.error("Destination escapes workspace.")
    if destination in destinations:
        parser.error("Duplicate destination.")
    if destination.exists() or destination.is_symlink():
        parser.error("Destination already exists: " + str(destination))
    destinations.add(destination)
    print(destination.relative_to(workspace))

if not args.execute:
    print("PLAN ONLY: Nine files prepared; nothing written.")
else:
    for destination, content in tasks:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("xb") as handle:
            handle.write(content)
    print("PASS: Five manifests, reference metadata and three reports prepared.")
