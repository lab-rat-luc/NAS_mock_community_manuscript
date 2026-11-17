# Preprocessing and mitochondrial mapping

This stage prepares the sequencing reads for downstream analyses. It merges whole-flowcell technical replicates, trims adapters, filters reads, extracts analytical subsets and maps reads to the provided mitochondrial references.

For installation, input organisation and path configuration, follow [WORKFLOW.md](WORKFLOW.md).

## Before starting

- Activate the `paper1-analysis` environment.
- Configure `REPO_DIR`, `SRA_FASTQ_DIR`, `SUPPORTING_DATA_DIR`, `OUTPUT_DIR`, `REFERENCE_FASTA` and `DORADO` as described in workflow.md.
- Run `workflow/prepare_inputs.py` to import the seven FASTQs and six sequencing summaries.

## Run preprocessing and mapping

Preview the processing plan:

```bash
bash "$REPO_DIR/workflow/run_preprocessing.sh" --plan
```

Execute:

```bash
bash "$REPO_DIR/workflow/run_preprocessing.sh" --execute
```

Plan mode does not execute analytical commands.

## Processing order

1. Copy the supplied mitochondrial reference FASTA into `OUTPUT_DIR/refs/`.
2. Merge untrimmed metagenomic barcodes in the order `barcode01`, `barcode02`, `barcode03`.
3. Merge untrimmed NAS barcodes in the order `barcode09`, `barcode10`, `barcode11`.
4. Build the minimap2 reference index.
5. Adapter-trim individual and merged FASTQs separately.
6. Length-filter reads and extract read subsets.
7. Map retained read subsets, sort the resulting BAMs and generate indexes and mapping statistics.

The general filtering and mapping scripts process both individual and merged inputs. Do not additionally run `filter_merged_fastqs.sh` or `map_merged_filtered_fastqs.sh`.

## Analytical settings

| Step | Setting |
|---|---|
| Whole-flowcell trimming kit | `SQK-NBD114-24` |
| Half-flowcell trimming kit | `SQK-LSK114` |
| Primer trimming | Disabled using `--no-trim-primers` |
| Read-length filter | 200–20,000 bp, inclusive |
| Additional read-quality filter | None |
| Accepted NAS reads | `end_reason == signal_positive` |
| Rejected NAS reads | `end_reason == data_service_unblock_mux_change` |
| Half-flowcell NAS region | Channels 1–256 |
| Half-flowcell metagenomic region | Channels 257–512 |
| Mapping | minimap2 `-ax map-ont --secondary=no` |
| Alignment retention | SAMtools `-q 20 -F 2308` |

The alignment filter retains mapped primary alignments with MAPQ ≥20 and excludes unmapped, secondary and supplementary alignments.

## Software and resources

Required tools are provided by the provided environment, except Dorado v1.4.0, which is installed separately. See [INSTALLATION.md](INSTALLATION.md).

The current scripts use:

- 16 Dorado threads.
- 12 minimap2 threads.
- 4 SAMtools view threads.
- 4 SAMtools sort threads.

The mapping commands form a pipeline, so their resource demands can overlap. Do not size a compute allocation using the minimap2 thread count alone.

Intermediate FASTQs and BAMs require substantial storage. Ensure sufficient space is available before starting.

## Outputs

All generated files remain inside `OUTPUT_DIR`:

| Directory | Contents |
|---|---|
| `refs/` | Mitochondrial reference FASTA and minimap2 index |
| `full-meta-dataset/` | Whole-flowcell metagenomic inputs, trimmed and filtered reads, BAMs and logs |
| `full-nas-dataset/` | Whole-flowcell NAS inputs, trimmed and filtered reads, BAMs and logs |
| `half-flowcell-dataset/` | Half-flowcell inputs, trimmed and filtered read subsets, BAMs and logs |

Canonical data table generation and manuscript analyses are separate stages. Continue with [WORKFLOW.md](WORKFLOW.md) after preprocessing completes successfully.

## If processing fails

The driver stops on errors and does not provide automatic restart or rollback.

Inspect the terminal error and relevant stage logs before rerunning. Completed or partial outputs may prevent a rerun. Do not delete results or disable safeguards without identifying the cause.
