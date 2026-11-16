#!/usr/bin/env python3
"""Restore Paper 1 input filenames from downloaded gzip files.

Dry run by default. Add --execute to decompress into an empty output directory.
"""
import argparse
import csv
import gzip
import shutil
from pathlib import Path


def safe_relative(value, single=False):
    p = Path(value)
    if not value or p.is_absolute() or '..' in p.parts or (single and len(p.parts) != 1):
        raise ValueError(f'Unsafe path: {value}')
    return p


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sra-dir', type=Path, required=True)
    parser.add_argument('--supporting-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    manifest = repo / 'metadata/input_mapping.tsv'
    sra = args.sra_dir.resolve()
    supporting = args.supporting_dir.resolve()
    output = args.output_dir.resolve()

    for directory in (sra, supporting):
        if not directory.is_dir():
            parser.error(f'Input directory missing: {directory}')
        if output == directory or directory in output.parents or output in directory.parents:
            parser.error('Output directory overlaps an input directory')
    if output == repo or output in repo.parents or repo in output.parents:
        parser.error('Output must not be the repository root or an ancestor')
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error('Output directory must be absent or empty')

    with manifest.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len(rows) != 13:
        parser.error('Expected 13 input mappings')

    tasks = []
    seen = set()
    counts = {'sra_fastq': 0, 'supporting_data': 0}
    for row in rows:
        kind = row['source_type']
        if kind not in counts:
            parser.error(f'Unknown source type: {kind}')
        filename = safe_relative(row['source_filename'], single=True)
        relative = safe_relative(row['destination_relative_path'])
        if filename.suffix != '.gz':
            parser.error(f'Expected gzip source: {filename}')
        source = (sra if kind == 'sra_fastq' else supporting) / filename
        destination = output / relative
        if destination in seen:
            parser.error(f'Duplicate destination: {relative}')
        if not source.is_file() or source.is_symlink():
            parser.error(f'Missing or symlinked input: {source}')
        seen.add(destination)
        counts[kind] += 1
        tasks.append((source, destination))
    if counts != {'sra_fastq': 7, 'supporting_data': 6}:
        parser.error(f'Unexpected input counts: {counts}')

    print('PASS: 13 input mappings checked')
    for source, destination in tasks:
        print(f'{source.name} -> {destination.relative_to(output)}')
    if not args.execute:
        print('DRY RUN: No files written')
        return

    output.mkdir(parents=True, exist_ok=True)
    for source, destination in tasks:
        destination.parent.mkdir(parents=True, exist_ok=True)
        partial = destination.with_name(destination.name + '.partial')
        if destination.exists() or destination.is_symlink() or partial.exists() or partial.is_symlink():
            raise RuntimeError(f'Output or partial already exists: {destination}')
        try:
            with gzip.open(source, 'rb') as incoming, partial.open('xb') as outgoing:
                shutil.copyfileobj(incoming, outgoing, length=8 * 1024 * 1024)
            partial.rename(destination)
        except Exception:
            print(f'ERROR: Inspect partial output, if present: {partial}')
            raise
        print(f'CREATED: {destination.relative_to(output)}')
    print('PASS: All 13 inputs restored')


if __name__ == '__main__':
    main()
