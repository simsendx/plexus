# multiplexDesigner

`multiplexDesigner` is a Python-based bioinformatics tool designed to automate the creation of multiplex PCR panels. It specifically targets workflows like cfDNA and simsen, integrating genomic data processing, primer design (via `primer3`), and specificity checking (via BLAST) to generate optimized primer sets for multiple targets simultaneously.

## Features

- **Automated Primer Design**: Uses `primer3-py` to generate primer candidates for provided genomic junctions.
- **Multiplex Optimization**: tailored for multiplex PCR panels.
- **Specificity Checking**: Integrates BLAST to check for off-target amplification and primer specificity.
- **Configuration Presets**: Includes default and lenient configuration presets for different design stringency.
- **CLI Interface**: Easy-to-use command line interface for running the design pipeline.

## Installation

### Prerequisites

- Python 3.10+
- NCBI BLAST+ (if running specificity checks locally)

### Using uv (Recommended)

This project uses `uv` for package management.

```bash
git clone https://github.com/simsendx/multiplexDesigner
cd multiplexDesigner
uv pip install -e .
```

### Using Conda

You can also set up the environment using Conda:

```bash
git clone https://github.com/simsendx/multiplexDesigner
cd multiplexDesigner
conda env create -f config/environment.yml
conda activate multiplexdesigner-run
pip install -e .
```

## Usage

The primary interface is the `multiplexdesigner` CLI.

### Basic Command

To run the complete design pipeline:

```bash
multiplexdesigner run \
  --input data/junctions.csv \
  --fasta data/genome.fa \
  --output results/ \
  --name my_panel
```

### Input file format

The input file should be a CSV containing the target junctions with the following columns:

```csv
Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate
EGFR_T790M,chr7,55181378,55181378
KRAS_G12D,chr12,25245350,25245350
...
```

See `data/junctions.csv` for a complete example.

### Key Options

- `-i, --input`: Path to the input CSV file (Required).
- `-f, --fasta`: Path to the reference genome FASTA file (Required).
- `-o, --output`: Output directory (Default: `./output`).
- `-p, --preset`: Configuration preset (`default` or `lenient`).
- `--skip-blast`: Skip the BLAST specificity check (faster, but less validation).

Run `multiplexdesigner --help` for a full list of commands and options.

## Configuration

The design parameters (melting temperature, primer length, penalties, etc.) are controlled via configuration files in the `config/` directory. You can supply a custom JSON config file using the `--config` option.

## Dependencies & References

- [primer3-py](https://github.com/libnano/primer3-py): Python interface to Primer3.
- [Biopython](https://biopython.org/): For BLAST integration and sequence handling.
- [Typer](https://typer.tiangolo.com/): CLI framework.

### Acknowledgements

Inspired by and utilizes concepts from:

- [Primer3](https://primer3.org/)
- [Multiply](https://github.com/JasonAHendry/multiply)
