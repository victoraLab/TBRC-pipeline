# TBRC Backend

This directory contains the Snakemake backend for TBRC.

It is responsible for:

- read trimming
- paired-read assembly
- barcode-based demultiplexing
- final FASTA generation
- QC generation
- optional IGBlast annotation
- optional clonality analysis
- packaging and optional result transfer

## Main Entry Points

- [`workflow/Snakefile`](backend/workflow/Snakefile)
- [`workflow/rules/`](backend/workflow/rules)
- [`workflow/scripts/`](backend/workflow/scripts)
- [`igblast/`](backend/igblast)

## Runtime Expectations

The backend is usually launched by the Shiny frontend and expects:

- `DATA_TABLE`: path to the generated `sample.json`
- `INPUTPATH`: path to the FASTQ input directory

The backend then reads run metadata from the JSON handoff file and executes the workflow.

## Environment

Create the Conda environment from the repository root [`environment.yml`](environment.yml):

```bash
conda env create -f environment.yml
conda activate snakemake
Rscript install_r_pins.R   # CRAN versions that differ from the conda builds (ggplate 0.2.0 etc.)
```

Useful verification:

```bash
which snakemake
which pandaseq
which cutadapt
which fastx_trimmer
which fastx_collapser
```

If you want QC, IGBlast, or clonality, you also need:

- `Rscript`
- `ggplot2`
- `dplyr`
- `readr`
- `ggplate`
- `clonality` for clonality analysis

## Local Run Pattern

```bash
cd /path/to/TBRC/backend
conda activate snakemake
export DATA_TABLE=/path/to/sample.json
export INPUTPATH=/path/to/input_fastqs
snakemake --cores 4
```

## Transfer Modes

Transfer behavior is controlled upstream in [`frontend/server_config.json`](frontend/server_config.json).

- `push`: backend pushes the final archive to the storage server
- `pull`: backend keeps the archive locally and another machine should pull it later

## IGBlast

See the full instructions in [`igblast_setup.md`](backend/igblast_setup.md).

The short version is:

```bash
cd /path/to/TBRC/backend/igblast
bash setup_igblast_databases.sh
./check_igblast_setup.sh mouse ig_all
./smoke_test_igblast.sh /path/to/sample.imgt.ready.fasta mouse ig_all
```
