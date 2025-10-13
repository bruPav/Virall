# Virall

A comprehensive tool for viral genome analysis including assembly, classification, gene prediction, and annotation.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioinformatics](https://img.shields.io/badge/bioinformatics-viral%20genomics-green.svg)](https://github.com/bruPav/virall)

## Features

- **Multi-read type support**: Short reads (Illumina), long reads (Oxford Nanopore, PacBio), and hybrid assembly. DNA and RNA Seq
- **Reference-guided assembly**: Optimized strategies for different read types
- **Viral classification**: Taxonomic classification using Kaiju
- **Gene prediction**: Automated gene finding with Prodigal
- **Functional annotation**: VOG database integration for viral protein annotation
- **Quality assessment**: CheckV integration for viral genome completeness

## Installation

### Quick Start

```bash
# Clone the repository
git clone https://github.com/bruPav/virall.git
cd virall

# Run the installation script
bash install.sh
```

The installation script will:
- Install all required dependencies (SPAdes, Flye, Kaiju, CheckV, etc.)
- Download and set up viral databases
- Install the Python package

## Usage

### Basic analysis

```bash
# Short reads only (single or paired-end)
virall analyse --single-reads-2 reads.fastq -o output_dir

virall analyse --short-reads-1 reads_1.fastq --short-reads-2 reads_2.fastq -o output_dir

# Long reads only
virall analyse --long-reads reads.fastq -o output_dir

# Hybrid assembly (short + long reads)
virall analyse --short-reads-1 reads_1.fastq --short-reads-2 reads_2.fastq --long-reads reads.fastq -o output_dir
```

### Reference-guided analysis
Add a '--reference' flag. For example

```bash

# Short reads With reference genome (single and paired end)
virall analyse --single-reads reads.fastq --reference ref.fasta -o output_dir

```

### Memory-efficient Mode
Add a '--mem-efficient' flag. For example
```bash
# For large datasets
virall analyse --long-reads large_dataset.fastq --mem-efficient -o output_dir
```

## All Commands

### `analyse` - Complete Pipeline
Run the full viral genome analysis pipeline from reads to annotated contigs.

```bash
virall analyse --short-reads-1 reads_1.fastq --short-reads-2 reads_2.fastq -o output_dir
```

### `preprocess` - Read Preprocessing
Quality control, trimming, and filtering of sequencing reads.

```bash
virall preprocess --short-reads-1 reads_1.fastq --short-reads-2 reads_2.fastq -o preprocessed/
```

### `assemble` - Assembly Only
Assemble reads and identify viral contigs without full analysis.

```bash
virall assemble --short-reads-1 reads_1.fastq --short-reads-2 reads_2.fastq -o assembly/
```

### `identify` - Viral Identification
Identify viral sequences in sequencing reads using mapping-based approaches.

```bash
virall identify --short-reads-1 reads_1.fastq --short-reads-2 reads_2.fastq -o identification/
```

### `classify` - Taxonomic Classification
Classify viral contigs using Kaiju nucleotide-based classification.

```bash
virall classify --contigs contigs.fasta -o classification/
```

### `validate` - Quality Assessment
Validate and assess quality of assembled viral genomes using CheckV.

```bash
virall validate --contigs viral_contigs.fasta -o validation/
```

### `annotate` - Gene Annotation
Annotate viral contigs with gene prediction and functional annotation.

```bash
virall annotate --contigs viral_contigs.fasta -o annotation/
```

### `quantify` - Viral Quantification
Quantify viral contigs by mapping reads back to them.

```bash
virall quantify --contigs viral_contigs.fasta --short-reads-1 reads_1.fastq --short-reads-2 reads_2.fastq -o quantification/
```

### `train-model` - Model Training
Train a viral sequence identification model (work in progress, not available yet).

## Command Line Options

```
Usage: virall analyse [OPTIONS]

Options:
  -1, --short-reads-1 PATH     Forward short reads (FASTQ)
  -2, --short-reads-2 PATH     Reverse short reads (FASTQ)
  -l, --long-reads PATH        Long reads (FASTQ)
  -s, --single-reads PATH      Single-end reads (FASTQ)
  -r, --reference PATH         Reference genome (FASTA)
  -o, --output-dir PATH        Output directory [required]
  -t, --threads INTEGER        Number of threads [default: 8]
  -m, --memory TEXT            Memory limit [default: 16G]
  -c, --config PATH            Configuration file
  --min-contig-length INTEGER  Minimum contig length [default: 1000]
  --viral-confidence FLOAT     Viral classification confidence [default: 0.7]
  --assembly-strategy TEXT     Assembly strategy [default: auto]
  --rna-mode                   Enable RNA mode
  --mem-efficient              Enable memory-efficient mode
  --help                       Show this message and exit.
```

## Output Structure

```
output_dir/
├── 00_qc/                     # Quality control reports
│   ├── fastqc/               # FastQC reports
│   └── trimmomatic/          # Trimming logs
├── 01_assembly/              # Assembly results
│   ├── spades/               # SPAdes output
│   ├── flye/                 # Flye output
│   └── reference_guided/     # Reference-guided assembly
├── 02_identification/        # Viral identification
│   ├── kaiju/                # Kaiju classification
│   └── viral_contigs.fasta   # Identified viral contigs
├── 03_validation/            # Quality validation
│   └── checkv/               # CheckV results
├── 04_annotation/            # Gene prediction and annotation
│   ├── prodigal/             # Gene predictions
│   └── vog/                  # VOG annotations
└── final_report.html         # Summary report
```

## Configuration

Create a custom configuration file to modify default parameters:

```yaml
# config.yaml
min_contig_length: 1000
min_coverage: 5
viral_confidence: 0.7
max_long_reads: 50000
long_read_subsample_size: 20000

databases:
  kaiju_db_path: "/path/to/kaiju_db"
  checkv_db_path: "/path/to/checkv_db"
  vog_db_path: "/path/to/vog_db"
```

## Dependencies

### Required Tools
- **SPAdes**: Short-read assembly
- **Flye**: Long-read assembly
- **Kaiju**: Taxonomic classification
- **CheckV**: Viral genome validation
- **Prodigal**: Gene prediction
- **HMMER**: Protein annotation
- **FastQC**: Quality control
- **Trimmomatic**: Read trimming
- **Porechop**: Long-read trimming

### Python Packages
- click
- pyyaml
- biopython
- pandas
- numpy

## Examples

### Example 1: Short-read Assembly

```bash
virall analyse \
  --short-reads-1 sample_R1.fastq.gz \
  --short-reads-2 sample_R2.fastq.gz \
  --threads 16 \
  --memory 32G \
  -o short_read_analysis
```

### Example 2: Long-read Reference-guided Assembly

```bash
virall analyse \
  --long-reads nanopore_reads.fastq \
  --reference viral_reference.fasta \
  --mem-efficient \
  -o long_read_analysis
```

### Example 3: Hybrid Assembly with RNA Mode

```bash
virall analyse \
  --short-reads-1 rna_R1.fastq \
  --short-reads-2 rna_R2.fastq \
  --long-reads rna_long.fastq \
  --rna-mode \
  -o rna_analysis
```

## Troubleshooting

### Common Issues

1. **Memory errors**: Use `--mem-efficient` flag for large datasets
2. **Assembly failures**: Check read quality and try different assembly strategies

### Getting Help

- Check the [Issues](https://github.com/bruPav/virall/issues) page
- Review the logs in the output directory
- Ensure all dependencies are properly installed



## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Acknowledgments

- SPAdes team for the assembly algorithm
- Flye team for long-read assembly
- Kaiju team for taxonomic classification
- CheckV team for viral genome validation
- VOG team for viral protein annotation
- This software was developed by Dr. Aurora Britania Diaz Fernandez and Bruno Pavletic, Msc

## Changelog

### v0.1.0
- Initial release
- Support for short and long reads
- Reference-guided assembly
- Viral classification and annotation
- Memory-efficient mode
