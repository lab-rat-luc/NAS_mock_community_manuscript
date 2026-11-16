# Installation

## Status

The unified environment has been installed successfully on the authors'
Linux HPC system. Specified Python and R packages load successfully.
End-to-end analytical reproduction remains pending.

## Requirements

- Linux-compatible system.
- Micromamba.
- Standard command-line utilities, including Bash, realpath and flock.
- Sufficient storage for decompressed reads and generated intermediates.
- Dorado 1.4.0, installed separately.

Resource requirements will be finalised after reproduction testing.

## Create the analysis environment

From the repository root:

    micromamba create --file environment.yml --strict-channel-priority --yes
    micromamba activate paper1-analysis

Use this same environment for input preparation, preprocessing,
canonical generation and manuscript analyses.

## Configure Dorado

Obtain Dorado 1.4.0 from the official Oxford Nanopore Dorado releases.
The final release documentation will include the verified download link.

Set the executable path for the current shell:

    export DORADO="/path/to/dorado-1.4.0/bin/dorado"
    "$DORADO" --version

The authors' executable reports 1.4.0+ba44a013.
Dorado is used for adapter trimming in this workflow, not rebasecalling.
No basecalling model download or POD5 input is required for trimming.

## Principal environment versions

- Python 3.10.19
- pandas 2.3.3
- NumPy 2.2.6
- SeqKit 2.11.0
- minimap2 2.30-r1287
- Samtools 1.22.1
- R 4.5.3

Exact requested R package versions are listed in environment.yml.
The reconstructed environment currently reports HTSlib 1.24 and pigz 2.8.

## Historical versus release versions

Historical filtering used SeqKit 2.8.2. The unified release environment
uses SeqKit 2.11.0, which was already used for canonical run-summary
statistics. Equivalence of the filtering results will be assessed during
internal reproduction testing.

Some original Python analyses ran under different Python/NumPy versions.
The unified environment is a release reconstruction, not an assertion
that every historical execution used the same environment.

## Check the installation

    python3 -c 'import numpy, pandas; print(numpy.__version__, pandas.__version__)'
    Rscript -e 'library(data.table); library(dplyr); library(ggplot2)'
    seqkit version
    minimap2 --version
    samtools --version
    pigz --version
    "$DORADO" --version

## Next steps

Download the SRA FASTQs and Zenodo supporting dataset as described in
README.md. Configure a separate reproduction workspace.

See PREPROCESSING.md for preprocessing instructions.
Do not run substantial processing on an HPC login node.
