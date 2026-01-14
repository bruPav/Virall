# Virall

<div align="center">
  <img src="logo.png" alt="Virall Logo" width="300">
</div>

A comprehensive tool for viral genome analysis including assembly, classification, gene prediction, and annotation.

[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioinformatics](https://img.shields.io/badge/bioinformatics-viral%20genomics-green.svg)](https://github.com/bruPav/virall)

## Features

- **Multi-read type support**: Short reads (Illumina), long reads (Oxford Nanopore, PacBio), and hybrid assembly. DNA and RNA Seq
- **Reference-guided assembly**: Optimized strategies for different read types
- **Viral classification**: Taxonomic classification using Kaiju
- **Gene prediction**: Automated gene finding with Prodigal
- **Functional annotation**: VOG database integration for viral protein annotation
- **Quality assessment**: CheckV integration for viral genome completeness
- **Interactive Plotting**: Dynamic visualizations for viral genome analysis
- **Single-cell RNA-seq (pooled)**: Supports pooled scRNA-seq by assembling cDNA (R2) with SPAdes rna-viral mode

## Workflow

See the [complete workflow diagram](viral_assembly_workflow.md) for a visual overview of the pipeline, including preprocessing, assembly, viral identification, and annotation steps.

## Installation

### Quick Start

```bash
# Clone the repository
git clone https://github.com/bruPav/Virall.git
cd Virall

# Run the installation script
bash install.sh
```

The installation script will:
- Install all required dependencies (SPAdes, Flye, Kaiju, CheckV, BWA, minimap2, samtools, fastp, fastplong, etc.)
- Download and set up viral databases (Kaiju viral database, VOG database)
- Install the Python package

## Usage

### Basic analysis

```bash
# Short reads only (single or paired-end)
virall analyse --single-reads reads.fastq -o output_dir

virall analyse --short-reads-1 reads_1.fastq --short-reads-2 reads_2.fastq -o output_dir

# Long reads only (ONT Nanopore)
virall analyse --nanopore nanopore_reads.fastq -o output_dir

# Long reads only (PacBio)
virall analyse --pacbio pacbio_reads.fastq -o output_dir

# Hybrid assembly (short + long reads)
virall analyse --short-reads-1 reads_1.fastq --short-reads-2 reads_2.fastq --nanopore long_reads.fastq -o output_dir
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
virall analyse --nanopore large_dataset.fastq --mem-efficient -o output_dir
```

### Host Filtering
Add a `--filter` flag to provide a host genome or sequences to filter out. For example:
```bash
# Filter out host reads (e.g., human) from analysis
virall analyse --single-reads reads.fastq --filter host_genome.fasta -o output_dir
```

### Single-cell RNA-seq (pooled)

```bash
# Provide paired-end inputs; R2 (cDNA) will be pooled and assembled as single-end
virall analyse \
  --single-cell \
  --short-reads-1 sample_R1.fastq.gz \
  --short-reads-2 sample_R2.fastq.gz \
  -o sc_rna_pooled_out
```

Notes for single-cell mode:
- Only RNA-seq is supported at the moment. DNA single-cell support will be added in future versions.
- The pipeline does not track individual cells yet. Reads are pooled across cells for assembly. Per-cell tracking and outputs will be added in future versions.
- Implementation details: when `--single-cell` is set, the pipeline treats R2 (cDNA) as single-end input, enables RNA mode automatically, and assembles using SPAdes `--rnaviral` without the genomic single-cell flag.
- Preprocessing: short-read trimming is performed with fastp (automatic adapter detection) according to configuration (`trim_adapters: true`).

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
virall classify --viral-contigs contigs.fasta -o classification/
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
  --short-reads-1 TEXT         Path to first mate of paired-end reads
  --short-reads-2 TEXT         Path to second mate of paired-end reads
  --nanopore TEXT              Path to ONT Nanopore long reads
  --pacbio TEXT                Path to PacBio long reads
  --single-reads TEXT          Path to single-end reads
  --reference TEXT             Optional reference genome for guided assembly
  -o, --output-dir TEXT        Output directory [required]
  -t, --threads INTEGER        Number of threads [default: 8]
  --memory TEXT                Memory allocation [default: 16G]
  --config TEXT                Configuration file path
  --min-contig-length INTEGER  Minimum contig length [default: 1000]
  --viral-confidence FLOAT     Viral sequence confidence threshold [default: 0.8]
  --assembly-strategy [hybrid|short_only|long_only]  Assembly strategy [default: hybrid]
  --rna-mode                   Enable RNA-specific assembly parameters
  -m, --mem-efficient          Enable memory-efficient mode with read subsampling for large datasets
  --single-cell                Enable single-cell sequencing mode (pooled scRNA-seq)
  --filter TEXT                Path to host genome for filtering (e.g. human.fna)
  --help                       Show this message and exit.
```

## Output Structure

```
output_dir/
├── 00_qc/
│   └── fastqc/                         # FastQC HTML/zip reports
├── 01_assemblies/                      # Assembly outputs (SPAdes/Flye)
├── 02_viral_contigs/                   # Viral contigs and scaffolds
│   ├── viral_contigs.fasta
│   └── viral_scaffolds.fasta
├── 03_classifications/                 # Kaiju results
│   └── kaiju_summary.tsv
├── 04_quality_assessment/              # CheckV results
├── 05_gene_predictions/                # Prodigal + VOG annotations
├── 06_quantification/                  # Read mapping and abundance
└── logs/                               # Pipeline logs (if enabled)
```

## Configuration

Create a custom configuration file to modify default parameters:

```yaml
# config.yaml
min_contig_length: 1000
min_coverage: 5
viral_confidence: 0.8
max_long_reads: 50000
long_read_subsample_size: 20000

databases:
  kaiju_db_path: "/path/to/kaiju_db"
  checkv_db_path: "/path/to/checkv_db"
  vog_db_path: "/path/to/vog_db"
```

## Dependencies

### Required Tools
- **SPAdes**: Short-read assembly (including RNA viral mode)
- **Flye**: Long-read assembly
- **Kaiju**: Taxonomic classification
- **CheckV**: Viral genome validation
- **Prodigal**: Gene prediction
- **HMMER**: Protein annotation (VOG database searches)
- **BWA**: Short-read mapping for quantification
- **minimap2**: Long-read mapping for quantification and polishing
- **samtools**: BAM file processing
- **bcftools**: Variant calling and consensus generation
- **pilon**: Assembly polishing
- **FastQC**: Quality control
- **fastp**: Short-read trimming (automatic adapter detection)
- **fastplong**: Long-read trimming (automatic adapter detection)

### Python Packages
- click (command-line interface)
- pyyaml (configuration files)
- biopython (sequence handling)
- pandas (data manipulation)
- numpy (numerical computing)
- scikit-learn (machine learning)
- loguru (logging)
- tqdm (progress bars)
- psutil (system utilities)
- matplotlib (plotting - for future features)
- seaborn (plotting - for future features)
- plotly (plotting - for future features)

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
  --nanopore nanopore_reads.fastq \
  --reference viral_reference.fasta \
  --mem-efficient \
  -o long_read_analysis
```

### Example 3: Hybrid Assembly with RNA Mode

```bash
virall analyse \
  --short-reads-1 rna_R1.fastq \
  --short-reads-2 rna_R2.fastq \
  --nanopore rna_long.fastq \
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
- This software was developed by Dr. Aurora Britania Diaz Fernandez and Bruno Pavletic, Msc together in collaboration with Nidia Trovao, PhD and Prof. Windy McNerny

## Changelog

### v0.2.1
- Added `--version` flag to CLI
- Added `--filter` option for host genome filtering
- Documentation updates

### v0.2.0
- Added interactive plotting capabilities (Sunburst charts, etc.)
- Improved database path resolution for containerized environments
- Enhanced documentation

### v0.1.1
- Support for short, long and Single Cell Seq reads
- Reference-guided assembly
- Viral assembly, classification, annotation, validation and quantification
- Memory-efficient mode
- Very short reads (<75 bp) also processed but only with reference-guided assembly
