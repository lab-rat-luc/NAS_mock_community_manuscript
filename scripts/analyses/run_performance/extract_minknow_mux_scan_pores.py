#!/usr/bin/env python3

import csv
import json
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
import argparse

_parser = argparse.ArgumentParser()
_parser.add_argument("--project-dir", required=True)
_args = _parser.parse_args()
PROJECT_DIR = Path(_args.project_dir).resolve()
ANALYSIS_DIR = PROJECT_DIR / 'analysis' / 'run_performance'
MANIFEST = PROJECT_DIR / 'canonical_data' / 'inputs' / 'ont_report_manifest.tsv'
OUTPUT_DIR = ANALYSIS_DIR / 'outputs'
LOG_DIR = ANALYSIS_DIR / 'logs'

OBSERVATIONS_OUT = OUTPUT_DIR / 'minknow_mux_scan_pore_observations.tsv'
SUMMARY_OUT = OUTPUT_DIR / 'minknow_mux_scan_pore_summary.tsv'
INTERVALS_OUT = OUTPUT_DIR / 'minknow_mux_scan_intervals.tsv'
LOG_OUT = LOG_DIR / 'extract_minknow_mux_scan_pores.log'
PROVENANCE_OUT = LOG_DIR / 'minknow_mux_scan_pore_provenance.json'

RUN_LABELS = {
    'full_meta': 'Whole-flowcell Metagenomic',
    'full_nas': 'Whole-flowcell NAS',
    'half_flowcell': 'Half-flowcell physical run',
}

EXPECTED_INITIAL_TOTAL = {
    'full_meta': 1574,
    'full_nas': 1669,
    'half_flowcell': 1312,
}

EXPECTED_OBSERVATIONS_PER_RUN = 16
MESSAGE_PATTERN = re.compile(
    r'has found a total of\s+(\d+)\s+pores\.\s*'
    r'(\d+)\s+pores available for immediate sequencing',
    re.IGNORECASE,
)


def parse_time(value):
    if not value:
        raise ValueError('Missing ISO timestamp')
    value = str(value)
    if value.endswith('Z'):
        value = value[:-1] + '+00:00'
    parsed = datetime.fromisoformat(value)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def write_tsv(path, rows, fieldnames):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(
            handle,
            delimiter='\t',
            fieldnames=fieldnames,
            lineterminator='\n',
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    started = datetime.now().astimezone().isoformat()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    if not MANIFEST.is_file():
        raise FileNotFoundError(f'Missing manifest: {MANIFEST}')

    with MANIFEST.open(newline='') as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter='\t'))

    if len(manifest_rows) != 3:
        raise ValueError(
            f'Expected 3 manifest rows, observed {len(manifest_rows)}'
        )

    observations = []
    interval_rows = []
    summary_rows = []

    for manifest_row in manifest_rows:
        run = manifest_row['run']
        if run not in RUN_LABELS:
            raise ValueError(f'Unexpected run in manifest: {run}')

        report_path = Path(manifest_row['report_json'])
        if not report_path.is_file():
            raise FileNotFoundError(f'Missing report JSON: {report_path}')

        with report_path.open() as handle:
            report = json.load(handle)

        messages = report.get('user_messages', [])
        if not isinstance(messages, list):
            raise ValueError(f'user_messages is not a list in {report_path}')

        sequencing_start_messages = [
            message
            for message in messages
            if isinstance(message, dict)
            and message.get('identifier') == 'sequencing.start'
        ]
        if len(sequencing_start_messages) != 1:
            raise ValueError(
                f'Expected one sequencing.start message for {run}, '
                f'observed {len(sequencing_start_messages)}'
            )

        sequencing_start_time = parse_time(
            sequencing_start_messages[0].get('time')
        )

        run_records = []

        for message_index, message in enumerate(messages):
            if not isinstance(message, dict):
                continue
            if message.get('identifier') != 'mux_scan_result':
                continue

            extra = message.get('extra_data', {})
            if not isinstance(extra, dict):
                raise ValueError(
                    f'extra_data is not an object for {run}, message {message_index}'
                )

            if 'total_pores' not in extra or 'num_pores' not in extra:
                raise ValueError(
                    f'Missing total_pores or num_pores for {run}, '
                    f'message {message_index}'
                )

            total_pores = int(extra['total_pores'])
            available_pores = int(extra['num_pores'])
            if total_pores <= 0 or available_pores <= 0:
                raise ValueError(
                    f'Non-positive pore count for {run}, message {message_index}'
                )
            if available_pores > total_pores:
                raise ValueError(
                    f'Available pores exceed total pores for {run}, '
                    f'message {message_index}'
                )

            observation_time = parse_time(message.get('time'))
            message_text = str(message.get('user_message', ''))
            match = MESSAGE_PATTERN.search(message_text)
            if match is None:
                raise ValueError(
                    f'Could not parse pore counts from user_message for '
                    f'{run}, message {message_index}: {message_text}'
                )

            message_total = int(match.group(1))
            message_available = int(match.group(2))
            if message_total != total_pores or message_available != available_pores:
                raise ValueError(
                    f'Message counts do not match extra_data for {run}, '
                    f'message {message_index}'
                )

            elapsed_hours = (
                observation_time - sequencing_start_time
            ).total_seconds() / 3600.0
            if elapsed_hours < 0:
                raise ValueError(
                    f'Mux scan precedes sequencing.start for {run}, '
                    f'message {message_index}'
                )

            run_records.append(
                {
                    'run': run,
                    'run_label': RUN_LABELS[run],
                    'source_json': str(report_path),
                    'message_index': message_index,
                    'flow_cell_id': str(extra.get('flow_cell_id', '')),
                    'sequencing_start_time_utc': sequencing_start_time.isoformat(),
                    'mux_scan_time_utc': observation_time.isoformat(),
                    'elapsed_hours_from_sequencing_start': elapsed_hours,
                    'total_pores_found': total_pores,
                    'pores_available_for_immediate_sequencing': available_pores,
                    'available_fraction_of_total_pores': (
                        available_pores / total_pores
                    ),
                    'user_message': message_text,
                }
            )

        run_records.sort(key=lambda row: row['mux_scan_time_utc'])

        if len(run_records) != EXPECTED_OBSERVATIONS_PER_RUN:
            raise ValueError(
                f'Expected {EXPECTED_OBSERVATIONS_PER_RUN} mux scans for {run}, '
                f'observed {len(run_records)}'
            )

        times = [parse_time(row['mux_scan_time_utc']) for row in run_records]
        if any(later <= earlier for earlier, later in zip(times, times[1:])):
            raise ValueError(f'Mux-scan times are not strictly increasing for {run}')

        initial_total = run_records[0]['total_pores_found']
        initial_available = run_records[0][
            'pores_available_for_immediate_sequencing'
        ]
        if initial_total != EXPECTED_INITIAL_TOTAL[run]:
            raise ValueError(
                f'Initial total_pores mismatch for {run}: expected '
                f'{EXPECTED_INITIAL_TOTAL[run]}, observed {initial_total}'
            )

        for observation_index, row in enumerate(run_records, start=1):
            row['observation_index'] = observation_index
            row['initial_total_pores_found'] = initial_total
            row['initial_pores_available_for_immediate_sequencing'] = (
                initial_available
            )
            row['total_pores_percent_of_initial'] = (
                100.0 * row['total_pores_found'] / initial_total
            )
            row[
                'available_pores_percent_of_initial_available'
            ] = (
                100.0
                * row['pores_available_for_immediate_sequencing']
                / initial_available
            )

        intervals = [
            (later - earlier).total_seconds() / 60.0
            for earlier, later in zip(times, times[1:])
        ]
        for interval_index, interval_minutes in enumerate(intervals, start=1):
            interval_rows.append(
                {
                    'run': run,
                    'run_label': RUN_LABELS[run],
                    'interval_index': interval_index,
                    'from_observation_index': interval_index,
                    'to_observation_index': interval_index + 1,
                    'interval_minutes': interval_minutes,
                }
            )

        total_series = [row['total_pores_found'] for row in run_records]
        available_series = [
            row['pores_available_for_immediate_sequencing']
            for row in run_records
        ]

        summary_rows.append(
            {
                'run': run,
                'run_label': RUN_LABELS[run],
                'flow_cell_id': run_records[0]['flow_cell_id'],
                'observations': len(run_records),
                'initial_total_pores_found': initial_total,
                'final_total_pores_found': total_series[-1],
                'final_total_pores_percent_of_initial': (
                    100.0 * total_series[-1] / initial_total
                ),
                'initial_pores_available_for_immediate_sequencing': (
                    initial_available
                ),
                'final_pores_available_for_immediate_sequencing': (
                    available_series[-1]
                ),
                'final_available_pores_percent_of_initial_available': (
                    100.0 * available_series[-1] / initial_available
                ),
                'minimum_available_pores': min(available_series),
                'maximum_available_pores': max(available_series),
                'available_pore_increase_events': sum(
                    later > earlier
                    for earlier, later in zip(
                        available_series, available_series[1:]
                    )
                ),
                'total_pore_increase_events': sum(
                    later > earlier
                    for earlier, later in zip(total_series, total_series[1:])
                ),
                'first_mux_scan_elapsed_hours': run_records[0][
                    'elapsed_hours_from_sequencing_start'
                ],
                'last_mux_scan_elapsed_hours': run_records[-1][
                    'elapsed_hours_from_sequencing_start'
                ],
                'mean_mux_scan_interval_minutes': (
                    sum(intervals) / len(intervals)
                ),
                'minimum_mux_scan_interval_minutes': min(intervals),
                'maximum_mux_scan_interval_minutes': max(intervals),
            }
        )

        observations.extend(run_records)

    if len(observations) != 48:
        raise ValueError(
            f'Expected 48 mux-scan observations, observed {len(observations)}'
        )
    if len(summary_rows) != 3:
        raise ValueError(
            f'Expected 3 summary rows, observed {len(summary_rows)}'
        )
    if len(interval_rows) != 45:
        raise ValueError(
            f'Expected 45 interval rows, observed {len(interval_rows)}'
        )

    observation_fields = [
        'run',
        'run_label',
        'source_json',
        'flow_cell_id',
        'observation_index',
        'message_index',
        'sequencing_start_time_utc',
        'mux_scan_time_utc',
        'elapsed_hours_from_sequencing_start',
        'total_pores_found',
        'initial_total_pores_found',
        'total_pores_percent_of_initial',
        'pores_available_for_immediate_sequencing',
        'initial_pores_available_for_immediate_sequencing',
        'available_pores_percent_of_initial_available',
        'available_fraction_of_total_pores',
        'user_message',
    ]

    summary_fields = list(summary_rows[0].keys())
    interval_fields = list(interval_rows[0].keys())

    write_tsv(OBSERVATIONS_OUT, observations, observation_fields)
    write_tsv(SUMMARY_OUT, summary_rows, summary_fields)
    write_tsv(INTERVALS_OUT, interval_rows, interval_fields)

    finished = datetime.now().astimezone().isoformat()
    provenance = {
        'analysis_name': 'minknow_mux_scan_pore_series',
        'started': started,
        'finished': finished,
        'manifest': str(MANIFEST),
        'physical_runs': list(RUN_LABELS),
        'observations_per_run': EXPECTED_OBSERVATIONS_PER_RUN,
        'total_observations': len(observations),
        'definitions': {
            'total_pores_found': (
                'extra_data.total_pores; message text: total pores found '
                'during the pore scan'
            ),
            'pores_available_for_immediate_sequencing': (
                'extra_data.num_pores; message text: pores available for '
                'immediate sequencing'
            ),
            'time_origin': 'sequencing.start message time for each physical run',
            'half_flowcell_scope': (
                'whole half-flowcell physical run; no regional allocation'
            ),
        },
        'validation': {
            'message_text_matches_extra_data': True,
            'strictly_increasing_timestamps': True,
            'available_pores_not_greater_than_total_pores': True,
            'expected_initial_total_pores_match': True,
            'expected_record_counts_match': True,
        },
        'outputs': {
            'observations': str(OBSERVATIONS_OUT),
            'summary': str(SUMMARY_OUT),
            'intervals': str(INTERVALS_OUT),
        },
        'software': {
            'python': sys.version,
            'platform': platform.platform(),
        },
    }

    with PROVENANCE_OUT.open('w') as handle:
        json.dump(provenance, handle, indent=2)

    with LOG_OUT.open('w') as handle:
        handle.write(f'Started: {started}\n')
        handle.write(f'Manifest: {MANIFEST}\n')
        handle.write('Metric 1: total pores found during each mux scan\n')
        handle.write(
            'Metric 2: pores available for immediate sequencing at each '
            'mux scan\n'
        )
        handle.write(
            'Half-flowcell scope: whole physical split-flowcell run; no '
            'regional allocation\n'
        )
        handle.write(f'Observations: {len(observations)}\n')
        handle.write(f'Summary rows: {len(summary_rows)}\n')
        handle.write(f'Interval rows: {len(interval_rows)}\n')
        handle.write(f'Observations output: {OBSERVATIONS_OUT}\n')
        handle.write(f'Summary output: {SUMMARY_OUT}\n')
        handle.write(f'Intervals output: {INTERVALS_OUT}\n')
        handle.write(f'Finished: {finished}\n')

    print(f'Wrote: {OBSERVATIONS_OUT}')
    print(f'Wrote: {SUMMARY_OUT}')
    print(f'Wrote: {INTERVALS_OUT}')
    print(f'Wrote: {PROVENANCE_OUT}')
    print(f'Observations: {len(observations)}')
    print(f'Summary rows: {len(summary_rows)}')


if __name__ == '__main__':
    main()
