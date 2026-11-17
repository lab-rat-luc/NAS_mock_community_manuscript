# Installation

Use the same `paper1-analysis` environment for input preparation, preprocessing, canonical data table generation and analyses.

## Requirements

- A Linux-compatible system.
- Micromamba installed and configured for shell activation.
- Standard linux command-line utilities, including Bash, `realpath` and `flock`.
- Sufficient storage for downloaded data, decompressed reads and generated files.
- Dorado v1.4.0, installed separately.

Run substantial processing on an appropriately resourced workstation or compute allocation.

## 1. Create the analysis environment

From the repository root, run:

```bash
micromamba create --file environment.yml --strict-channel-priority --yes
micromamba activate paper1-analysis
```

The environment specifications are provided in ../environment.yml.

## 2. Install and configure Dorado

Download Dorado v1.4.0 (Linux) from the official Oxford Nanopore Dorado releases.

Set the path to the installed executable:

```bash
export DORADO="/path/to/dorado-1.4.0/bin/dorado"
"$DORADO" --version
```

Replace the example path with your installation path. The executable used for this workflow is `1.4.0+ba44a013`.

Dorado is used for adapter trimming, not basecalling. No basecalling model download or POD5 files are required.

Set `DORADO` again when opening a new terminal or compute session, unless it is already configured in your shell environment.

## 3. Check the installation

With `paper1-analysis` active, run:

```bash
python3 -c 'import numpy, pandas; print(numpy.__version__, pandas.__version__)'
Rscript -e 'library(data.table); library(dplyr); library(ggplot2)'
seqkit version
minimap2 --version
samtools --version
pigz --version
"$DORADO" --version
```

Check the required system utilities:

```bash
command -v bash realpath flock
```

The Python and R commands should complete without import errors. Each tool should report its version, and each system utility should return an executable path.

### Principal software versions

| Software | Version |
|---|---|
| Python | 3.10.19 |
| pandas | 2.3.3 |
| NumPy | 2.2.6 |
| SeqKit | 2.11.0 |
| minimap2 | 2.30-r1287 |
| SAMtools | 1.22.1 |
| R | 4.5.3 |
| pigz | 2.8 |
| Dorado | 1.4.0 |

See ../environment.yml for the environment specification, including R packages.

## Next steps

Follow [WORKFLOW.md](WORKFLOW.md) to download and organise the inputs, configure paths and run the analyses.

Set `OUTPUT_DIR` to a new directory outside the code repository and the directories containing your downloaded data.

See [PREPROCESSING.md](PREPROCESSING.md) for preprocessing settings and thread usage.
