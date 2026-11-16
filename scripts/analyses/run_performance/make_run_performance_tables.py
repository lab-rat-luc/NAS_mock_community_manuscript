#!/usr/bin/env python3

import json
import math
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
ANALYSIS_DIR = PROJECT_DIR / 'analysis' / 'run_performance'
INPUT_DIR = ANALYSIS_DIR / 'inputs'
OUTPUT_DIR = ANALYSIS_DIR / 'outputs'
LOG_DIR = ANALYSIS_DIR / 'logs'

INPUT_FILES = {
    'full_meta': (
        INPUT_DIR
        / 'full_meta_sequencing_summary_sup.tsv'
    ),
    'full_nas': (
        INPUT_DIR
        / 'full_nas_sequencing_summary_sup.tsv'
    ),
    'half_flowcell': (
        INPUT_DIR
        / 'half_flowcell_sequencing_summary_sup.tsv'
    ),
}



BIN_MINUTES = 10
BIN_SECONDS = BIN_MINUTES * 60
RUN_HOURS = 24
N_BINS = RUN_HOURS * 60 // BIN_MINUTES

BINNED_OUTPUT = OUTPUT_DIR / 'run_performance_binned.tsv'
SUMMARY_OUTPUT = OUTPUT_DIR / 'run_performance_summary.tsv'
VALIDATION_OUTPUT = OUTPUT_DIR / 'run_performance_validation.tsv'
PROVENANCE_OUTPUT = LOG_DIR / 'run_performance_provenance.json'
LOG_OUTPUT = LOG_DIR / 'run_performance_run.log'

DATASET_ORDER = ['full_meta', 'full_nas', 'half_meta', 'half_nas']
DATASET_METADATA = {
    'full_meta': {
        'run_type': 'Whole-flowcell',
        'library_type': 'meta',
        'method': 'Metagenomic',
        'physical_run': 'full_meta',
    },
    'full_nas': {
        'run_type': 'Whole-flowcell',
        'library_type': 'nas',
        'method': 'NAS',
        'physical_run': 'full_nas',
    },
    'half_meta': {
        'run_type': 'Half-flowcell',
        'library_type': 'meta',
        'method': 'Metagenomic',
        'physical_run': 'half_flowcell',
    },
    'half_nas': {
        'run_type': 'Half-flowcell',
        'library_type': 'nas',
        'method': 'NAS',
        'physical_run': 'half_flowcell',
    },
}


RAW_SUP_USECOLS = [
    'read_id',
    'channel',
    'template_start',
    'start_time',
    'sequence_length_template',
]

RAW_SUP_DTYPES = {
    'read_id': 'string',
    'channel': 'int64',
    'template_start': 'float64',
    'start_time': 'float64',
    'sequence_length_template': 'int64',
}


def assign_dataset(chunk, physical_run):
    if physical_run == 'full_meta':
        return pd.Series(
            'full_meta',
            index=chunk.index,
            dtype='string',
        )

    if physical_run == 'full_nas':
        return pd.Series(
            'full_nas',
            index=chunk.index,
            dtype='string',
        )

    if physical_run != 'half_flowcell':
        raise ValueError(
            f'Unsupported physical run: {physical_run}'
        )

    invalid_channel = ~chunk['channel'].between(
        1,
        512,
    )

    if invalid_channel.any():
        bad_channels = sorted(
            int(value)
            for value in chunk.loc[
                invalid_channel,
                'channel',
            ].unique()
        )

        raise ValueError(
            'Half-flowcell channels outside 1-512: '
            f'{bad_channels[:20]}'
        )

    return pd.Series(
        np.where(
            chunk['channel'].between(1, 256),
            'half_nas',
            'half_meta',
        ),
        index=chunk.index,
        dtype='string',
    )


def iter_raw_sup_chunks():
    for physical_run, input_file in INPUT_FILES.items():
        for chunk in pd.read_csv(
            input_file,
            sep='\t',
            usecols=RAW_SUP_USECOLS,
            dtype=RAW_SUP_DTYPES,
            chunksize=500_000,
        ):
            chunk['run'] = physical_run
            chunk['dataset'] = assign_dataset(
                chunk,
                physical_run,
            )
            yield chunk


def main():
    started = datetime.now().astimezone().isoformat()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    missing_inputs = [
        input_file
        for input_file in INPUT_FILES.values()
        if not input_file.is_file()
    ]

    if missing_inputs:
        raise FileNotFoundError(
            'Missing raw SUP input files: '
            + ', '.join(
                str(input_file)
                for input_file in missing_inputs
            )
        )

    empty_inputs = [
        input_file
        for input_file in INPUT_FILES.values()
        if input_file.stat().st_size <= 0
    ]

    if empty_inputs:
        raise ValueError(
            'Empty raw SUP input files: '
            + ', '.join(
                str(input_file)
                for input_file in empty_inputs
            )
        )

    physical_run_minima = {
        'full_meta': math.inf,
        'full_nas': math.inf,
        'half_flowcell': math.inf,
    }

    # Pass 1: physical-run time origins. Preserve prior priority: template_start.
    for chunk in iter_raw_sup_chunks():
        if chunk['template_start'].isna().any():
            raise ValueError('Missing template_start values detected')
        for run_name in physical_run_minima:
            run_chunk = chunk.loc[chunk['run'] == run_name, 'template_start']
            if not run_chunk.empty:
                physical_run_minima[run_name] = min(
                    physical_run_minima[run_name],
                    float(run_chunk.min()),
                )

    for run_name, minimum in physical_run_minima.items():
        if not math.isfinite(minimum):
            raise ValueError(f'No template_start minimum found for {run_name}')

    aggregate = {
        dataset: {
            'pass_reads': np.zeros(N_BINS, dtype=np.int64),
            'pass_bases': np.zeros(N_BINS, dtype=np.int64),
            'read_length_sum': np.zeros(N_BINS, dtype=np.float64),
            'read_lengths': [[] for _ in range(N_BINS)],
            'channels': [set() for _ in range(N_BINS)],
            'all_rows': 0,
            'rows_within_24h': 0,
            'rows_after_24h': 0,
            'bases_all_rows': 0,
            'bases_within_24h': 0,
        }
        for dataset in DATASET_ORDER
    }

    # Pass 2: aggregate all raw SUP pass reads into 10-minute bins.
    for chunk in iter_raw_sup_chunks():
        if (chunk['sequence_length_template'] <= 0).any():
            raise ValueError('Non-positive sequence_length_template detected')

        if chunk['dataset'].isna().any():
            raise ValueError(
                'Unassigned raw SUP dataset rows detected'
            )

        chunk['run_minimum'] = chunk['run'].map(physical_run_minima)
        chunk['elapsed_seconds'] = chunk['template_start'] - chunk['run_minimum']
        if (chunk['elapsed_seconds'] < 0).any():
            raise ValueError('Negative elapsed time detected')

        chunk['bin_index'] = np.floor(
            chunk['elapsed_seconds'] / BIN_SECONDS
        ).astype(np.int64)

        for dataset in DATASET_ORDER:
            dataset_chunk = chunk.loc[chunk['dataset'] == dataset]
            if dataset_chunk.empty:
                continue

            stats = aggregate[dataset]
            stats['all_rows'] += int(len(dataset_chunk))
            stats['bases_all_rows'] += int(
                dataset_chunk['sequence_length_template'].sum()
            )

            within = dataset_chunk.loc[
                dataset_chunk['bin_index'].between(0, N_BINS - 1)
            ].copy()
            after = len(dataset_chunk) - len(within)
            stats['rows_within_24h'] += int(len(within))
            stats['rows_after_24h'] += int(after)
            stats['bases_within_24h'] += int(
                within['sequence_length_template'].sum()
            )

            if within.empty:
                continue

            for bin_index, bin_data in within.groupby('bin_index', sort=False):
                idx = int(bin_index)
                lengths = bin_data['sequence_length_template'].to_numpy(
                    dtype=np.int64
                )
                stats['pass_reads'][idx] += len(lengths)
                stats['pass_bases'][idx] += int(lengths.sum())
                stats['read_length_sum'][idx] += float(lengths.sum())
                stats['read_lengths'][idx].extend(lengths.tolist())
                stats['channels'][idx].update(
                    int(value) for value in bin_data['channel'].unique()
                )

    binned_rows = []
    summary_rows = []
    validation_rows = []

    for dataset in DATASET_ORDER:
        metadata = DATASET_METADATA[dataset]
        stats = aggregate[dataset]
        active_channels = np.array(
            [len(values) for values in stats['channels']],
            dtype=np.int64,
        )
        nonzero = active_channels[active_channels > 0]
        if len(nonzero) == 0:
            raise ValueError(f'No productive channels found for {dataset}')
        starting_productive_channels = int(nonzero[0])

        cumulative_reads = np.cumsum(stats['pass_reads'])
        cumulative_bases = np.cumsum(stats['pass_bases'])

        for idx in range(N_BINS):
            lengths = stats['read_lengths'][idx]
            mean_length = float(np.mean(lengths)) if lengths else 0.0
            median_length = float(np.median(lengths)) if lengths else 0.0
            productive_count = int(active_channels[idx])

            binned_rows.append(
                {
                    'dataset': dataset,
                    'dataset_label': (
                        f"{metadata['run_type']} {metadata['method']}"
                    ),
                    'library_type': metadata['library_type'],
                    'method': metadata['method'],
                    'run_type': metadata['run_type'],
                    'physical_run': metadata['physical_run'],
                    'time_column_used': 'template_start',
                    'bin_minutes': BIN_MINUTES,
                    'time_bin_sec': idx * BIN_SECONDS,
                    'time_bin_hours': idx * BIN_SECONDS / 3600.0,
                    'pass_reads': int(stats['pass_reads'][idx]),
                    'pass_bases': int(stats['pass_bases'][idx]),
                    'productive_channels': productive_count,
                    'starting_productive_channels': (
                        starting_productive_channels
                    ),
                    'productive_channels_pct_of_first_nonzero_bin': (
                        100.0
                        * productive_count
                        / starting_productive_channels
                    ),
                    'mean_pass_read_length': mean_length,
                    'median_pass_read_length': median_length,
                    'cumulative_pass_reads': int(cumulative_reads[idx]),
                    'cumulative_pass_bases': int(cumulative_bases[idx]),
                    'cumulative_pass_bases_Gb': (
                        float(cumulative_bases[idx]) / 1e9
                    ),
                }
            )

        summary_rows.append(
            {
                'dataset': dataset,
                'dataset_label': f"{metadata['run_type']} {metadata['method']}",
                'run_type': metadata['run_type'],
                'method': metadata['method'],
                'n_bins': N_BINS,
                'bin_minutes': BIN_MINUTES,
                'starting_productive_channels': starting_productive_channels,
                'maximum_productive_channels': int(active_channels.max()),
                'final_productive_channels': int(active_channels[-1]),
                'all_pass_reads': int(stats['all_rows']),
                'pass_reads_within_24h': int(stats['rows_within_24h']),
                'pass_reads_after_24h': int(stats['rows_after_24h']),
                'all_pass_bases': int(stats['bases_all_rows']),
                'pass_bases_within_24h': int(stats['bases_within_24h']),
                'final_cumulative_pass_reads': int(cumulative_reads[-1]),
                'final_cumulative_pass_bases': int(cumulative_bases[-1]),
                'final_cumulative_pass_bases_Gb': float(cumulative_bases[-1]) / 1e9,
            }
        )

        validation_rows.append(
            {
                'dataset': dataset,
                'physical_run': metadata['physical_run'],
                'physical_run_minimum_template_start': (
                    physical_run_minima[metadata['physical_run']]
                ),
                'all_pass_reads': int(stats['all_rows']),
                'within_24h_plus_after_24h': int(
                    stats['rows_within_24h'] + stats['rows_after_24h']
                ),
                'row_partition_passed': (
                    stats['all_rows']
                    == stats['rows_within_24h'] + stats['rows_after_24h']
                ),
                'binned_read_sum': int(stats['pass_reads'].sum()),
                'binned_base_sum': int(stats['pass_bases'].sum()),
                'read_sum_passed': (
                    int(stats['pass_reads'].sum())
                    == stats['rows_within_24h']
                ),
                'base_sum_passed': (
                    int(stats['pass_bases'].sum())
                    == stats['bases_within_24h']
                ),
                'bin_count_passed': N_BINS == 144,
            }
        )

    binned = pd.DataFrame(binned_rows)
    summary = pd.DataFrame(summary_rows)
    validation = pd.DataFrame(validation_rows)

    if len(binned) != 576:
        raise ValueError(f'Expected 576 binned rows, observed {len(binned)}')
    if len(summary) != 4:
        raise ValueError(f'Expected 4 summary rows, observed {len(summary)}')
    if len(validation) != 4:
        raise ValueError(
            f'Expected 4 validation rows, observed {len(validation)}'
        )
    if not validation[
        [
            'row_partition_passed',
            'read_sum_passed',
            'base_sum_passed',
            'bin_count_passed',
        ]
    ].all().all():
        raise ValueError('One or more run-performance validations failed')

    binned.to_csv(BINNED_OUTPUT, sep='\t', index=False)
    summary.to_csv(SUMMARY_OUTPUT, sep='\t', index=False)
    validation.to_csv(VALIDATION_OUTPUT, sep='\t', index=False)

    finished = datetime.now().astimezone().isoformat()
    provenance = {
        'analysis_name': 'run_performance',
        'started': started,
        'finished': finished,
        'inputs': {
            scope: str(input_file)
            for scope, input_file in INPUT_FILES.items()
        },
        'supporting_data_source': (
            'Sequencing summaries supplied with the Zenodo '
            'supporting dataset linked from the repository README'
        ),
        'input_definition': (
            'Offline Dorado SUP pass-read summaries '
            'generated before demultiplexing'
        ),
        'datasets': DATASET_METADATA,
        'bin_minutes': BIN_MINUTES,
        'run_hours': RUN_HOURS,
        'time_column': 'template_start',
        'time_origin': 'minimum template_start within each physical run',
        'half_flowcell_clock': (
            'half_meta and half_nas share the half_flowcell physical-run origin'
        ),
        'pass_read_definition': (
            'All rows in pre-demultiplexing offline '
            'Dorado SUP summaries; all source rows '
            'represent SUP pass reads'
        ),
        'yield_definition': (
            'sum of sequence_length_template from '
            'pre-demultiplexing offline SUP summaries'
        ),
        'half_flowcell_allocation': {
            'half_nas': 'channels 1-256',
            'half_meta': 'channels 257-512',
        },
        'productive_channel_definition': (
            'number of unique channel values producing at least one pass read '
            'within the 10-minute bin'
        ),
        'starting_productive_channel_definition': (
            'productive channels in the first non-zero 10-minute bin for each '
            'analytical dataset'
        ),
        'percentage_definition': (
            'productive channels in bin divided by productive channels in the '
            'first non-zero bin, multiplied by 100'
        ),
        'interpretation_warning': (
            'Productive channels are not direct measurements of active pores. '
            'Direct pore analysis will use MinKNOW report JSON separately.'
        ),
        'outputs': {
            'binned': str(BINNED_OUTPUT),
            'summary': str(SUMMARY_OUTPUT),
            'validation': str(VALIDATION_OUTPUT),
        },
        'validation_passed': True,
        'software': {
            'python': sys.version,
            'pandas': pd.__version__,
            'numpy': np.__version__,
            'platform': platform.platform(),
        },
    }
    with open(PROVENANCE_OUTPUT, 'w') as handle:
        json.dump(provenance, handle, indent=2)

    with open(LOG_OUTPUT, 'w') as handle:
        handle.write(f'Started: {started}\n')
        for scope, input_file in INPUT_FILES.items():
            handle.write(
                f'Input [{scope}]: {input_file}\n'
            )

        handle.write(
            'Input definition: pre-demultiplexing '
            'offline Dorado SUP pass-read summaries\n'
        )
        handle.write(
            'Half-flowcell NAS allocation: '
            'channels 1-256\n'
        )
        handle.write(
            'Half-flowcell metagenomic allocation: '
            'channels 257-512\n'
        )
        handle.write(f'Bin size: {BIN_MINUTES} minutes\n')
        handle.write(f'Analysis window: {RUN_HOURS} hours\n')
        handle.write('Time column: template_start\n')
        handle.write(
            'Productive channel: unique channel producing at least one pass '
            'read in a bin\n'
        )
        handle.write(
            'Percentage denominator: productive channels in first non-zero '
            'bin, not active pores\n'
        )
        handle.write(f'Binned rows: {len(binned)}\n')
        handle.write(f'Summary rows: {len(summary)}\n')
        handle.write(f'Validation rows: {len(validation)}\n')
        handle.write(f'Binned output: {BINNED_OUTPUT}\n')
        handle.write(f'Summary output: {SUMMARY_OUTPUT}\n')
        handle.write(f'Validation output: {VALIDATION_OUTPUT}\n')
        handle.write(f'Finished: {finished}\n')

    print(f'Wrote: {BINNED_OUTPUT}')
    print(f'Wrote: {SUMMARY_OUTPUT}')
    print(f'Wrote: {VALIDATION_OUTPUT}')
    print(f'Wrote: {PROVENANCE_OUTPUT}')
    print(f'Binned rows: {len(binned)}')
    print(f'Summary rows: {len(summary)}')


if __name__ == '__main__':
    main()
