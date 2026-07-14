# Virall

<div align="center">
  <img src="logo.png" alt="Virall Logo" width="300">
</div>

A comprehensive Nextflow pipeline for viral metagenome analysis including assembly, classification, gene prediction, annotation, quality assessment, and quantification.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioinformatics](https://img.shields.io/badge/bioinformatics-viral%20genomics-green.svg)](https://github.com/bruPav/Virall)

## Features

- **Multi-read type support**: Short reads (Illumina / Ion Torrent), long reads (Oxford Nanopore, PacBio), and hybrid assembly
- **RNA virus support**: SPAdes `--rnaviral` mode for RNA-seq data
- **Reference-guided assembly**: Reference-guided consensus (medaka, samtools consensus) alongside de novo assembly
- **Reference detection**: Maps reads to a reference genome and reports mapping rate, coverage breadth, and mean depth
- **Host filtering**: Remove host reads before assembly using minimap2
- **Viral classification**: Taxonomic classification using Kaiju (virus-subset database)
- **Viral filtering**: Extracts viral contigs based on Kaiju taxonomy
- **Quality assessment**: CheckV (bacteriophages) + geNomad (RNA viruses, eukaryotic viruses) with merged output
- **Gene prediction & annotation**: prodigal-gv + VOG HMM database; results organised by taxonomy
- **Abundance quantification**: BWA / minimap2 read mapping with MAPQ-filtered depth output
- **Interactive plots**: Taxonomy sunburst, contig length distribution, coverage heatmap, quality summary
- **Single-cell 10x Genomics**: Barcode / UMI extraction, per-cell viral UMI counting, MTX matrix output
- **Containerised**: All tools and databases bundled in Docker / Singularity images — nothing else to install

## Workflow

See the [complete workflow diagram](viral_assembly_workflow.md) for a visual overview of the pipeline steps.

---

## Quick Start

### Prerequisites

| Requirement | Version |
|-------------|---------|
| [Nextflow](https://www.nextflow.io/) | ≥ 24.10.0 |
| [Docker](https://docs.docker.com/get-docker/) **or** [Singularity / Apptainer](https://apptainer.org/) | any recent |

**macOS users:** Install Docker Desktop (`brew install --cask docker`) and Nextflow (`brew install nextflow`). The pipeline runs identically via Docker on macOS — no native tool installation needed.

**No container?** Set the `VIRALL_DATABASE_DIR` environment variable, or place databases in `$HOME/.virall/databases/` (see [configuration](#3-configure-run-parameters)).

### 1. Get the container image

**Docker** — pull the pre-built image from Docker Hub:

```bash
docker pull pavle17/virall
```

Or build locally from the repository:

```bash
docker build -t pavle17/virall .
```

**Singularity / Apptainer** — pull from Sylabs Cloud:

```bash
singularity pull virall.sif library://brupav/virall/virall:latest
```

Or build locally:

```bash
singularity build --fakeroot virall.sif virall.def
```

> The container includes all bioinformatics tools and pre-indexed databases (~6 GB compressed, ~21 GB on disk).

---

### 2. Prepare your sample sheet

Create (or edit) `nextflow/samples.csv`. Leave columns empty when a read type is absent.

| Column | Description |
|--------|-------------|
| `sample_id` | Unique sample name (required) |
| `illumina_r1` / `illumina_r2` | Paired-end Illumina short reads |
| `illumina_single` | Single-end Illumina short reads |
| `iontorrent_r1` / `iontorrent_r2` | Paired-end Ion Torrent reads |
| `iontorrent_single` | Single-end Ion Torrent reads |
| `nanopore_reads` | Oxford Nanopore long reads |
| `pacbio_reads` | PacBio long reads |
| `sc_read1` / `sc_read2` | 10x Genomics single-cell reads |

Legacy columns `read1`, `read2`, `single`, and `long` are still accepted for backward compatibility.

Example:

```csv
sample_id,illumina_r1,illumina_r2,nanopore_reads
sample1,/data/s1_R1.fastq.gz,/data/s1_R2.fastq.gz,
sample2,,,/data/s2_long.fastq.gz
```

---

### 3. Configure run parameters

Copy and edit the template parameter file:

```bash
cp nextflow/run_params.yaml my_run.yaml
```

Key settings:

```yaml
samples: "nextflow/samples.csv"
outdir: "/full/path/to/results"

# Assembly strategy: "auto" | "hybrid" | "short_only" | "long_only"
assembly_strategy: "auto"

# Long-read technology: "nanopore" | "pacbio"
long_read_tech: "nanopore"

# Optional: reference genome for guided assembly + detection report
reference: null

# Optional: host genome FASTA to filter out host reads
host_genome: null

# Set true for RNA-seq data (SPAdes --rnaviral)
rna_mode: false

# Resources
threads: 8
memory: "16G"
```

See [run_params.yaml](nextflow/run_params.yaml) for the full parameter reference including single-cell, Ion Torrent, metaviral mode, Flye tuning, and quantification thresholds.

---

### 4. Run the pipeline

**Docker:**

```bash
nextflow run nextflow/main.nf \
  -params-file my_run.yaml \
  -profile docker \
  -resume
```

**Singularity:**

```bash
nextflow run nextflow/main.nf \
  -params-file my_run.yaml \
  -profile singularity \
  -resume
```

**SLURM cluster (with Singularity):**

```bash
nextflow run nextflow/main.nf \
  -params-file my_run.yaml \
  -profile singularity,slurm \
  -resume
```

> `-resume` skips steps that already completed successfully — useful after failures or parameter changes.

---

## Output Structure

Results are written per sample under the output directory:

```
results/
└── <sample_id>/
    ├── 00_preprocess/              # Trimmed reads (fastp / fastplong) and QC reports
    │   └── host_filtered/          # (if host_genome set) host-removed reads
    ├── 01_assembly/                # Assembled contigs and scaffolds (SPAdes / Flye)
    ├── 02_viral_contigs/           # Kaiju-filtered viral contigs + taxonomy-based renaming
    ├── 03_classifications/         # Kaiju results (raw + with taxon names)
    ├── 04_quality_assessment/      # CheckV (phages) + geNomad (RNA/eukaryotic viruses)
    │   └── merged_quality/         # Merged quality table across both tools
    ├── 05_gene_predictions/        # prodigal-gv proteins + VOG HMM annotation, by taxonomy
    ├── 06_quantification/          # BWA/minimap2 depth, MAPQ-filtered BAMs
    ├── 07_plots/                   # Taxonomy, contig length, coverage, quality plots
    └── 08_single_cell/             # (single-cell mode) barcode counts, MTX matrix
    └── 08_reference_check/         # (if reference set) detection report + coverage plot
```

---

## Pipeline Parameters

Run `nextflow run nextflow/main.nf --help` for a complete, grouped parameter reference — or see [nextflow_schema.json](nextflow/nextflow_schema.json).

Key parameters are covered in the [Quick Start](#3-configure-run-parameters) above. Database paths (`kaiju_db`, `checkv_db`, `genomad_db`, `vog_db`) are set automatically by `-profile docker` or `-profile singularity`. For local runs, databases are resolved from `$VIRALL_DATABASE_DIR/` or `$HOME/.virall/databases/`.

---

## Troubleshooting

**Pipeline exits immediately with "MalformedInputException" or encoding error**
Your sample sheet CSV may contain non-UTF-8 characters. The pipeline automatically falls back to ISO-8859-1 encoding with a logged warning. If the error persists, save your CSV as plain UTF-8 from your editor.

**"Sample sheet missing required column 'sample_id'" or no sequencing columns detected**
The pipeline validates your sample sheet at startup. Verify column names match [the sample sheet format](#2-prepare-your-sample-sheet). These checks appear in the startup banner.

**No viral contigs found**
Kaiju classifies against the virus-subset database. If no contigs are classified as viral, the pipeline continues and downstream steps produce empty-but-valid output with informative messages.

**Assembly empty / all reads filtered as host**
Check `00_preprocess/host_filtered/` for read counts. If all reads map to the host genome, assembly is intentionally skipped with a warning.

**CheckV / geNomad database not found**
When running outside a container, set `kaiju_db`, `checkv_db`, `genomad_db`, and `vog_db` in your params file, or export `VIRALL_DATABASE_DIR` pointing to a directory containing `kaiju_db/`, `checkv_db/`, `genomad_db/`, and `vog_db/` subdirectories.

**Out of memory during assembly**
Increase `memory` in your params file and ensure your `nextflow.config` process block for `ASSEMBLE` reflects the new value. For very large metagenomes, consider `flye_genome_size` to help Flye pre-allocate correctly.

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Acknowledgments

- SPAdes team for the assembly algorithm
- Flye team for long-read assembly
- Kaiju team for taxonomic classification
- CheckV team for viral genome quality assessment
- geNomad team for RNA virus and eukaryotic virus identification
- VOG team for viral protein annotation
- This software was developed by Dr. Aurora Britania Diaz Fernandez and Bruno Pavletic, MSc, in collaboration with Nidia Trovao, PhD and Prof. Windy McNerny

---

## Changelog

### v1.0.0
- **Architecture:** Modularized pipeline into 8 DSL2 module files (down from 1666-line monolith)
- **UX:** Added pre-flight checks (sample sheet validation, startup banner, completion summary)
- **UX:** Added `nextflow_schema.json` with all 30 params documented; `--help` now shows grouped reference
- **UX:** Added `manifest` block (version, description, homepage) to `nextflow.config`
- **Robustness:** Replaced implicit `|| true` error suppression with logged warnings on critical tool failures
- **Robustness:** Fixed CSV encoding crash on non-UTF-8 sample sheets (ISO-8859-1 fallback)
- **Correctness:** Fixed fuzzy contig matching producing false species labels in plots
- **Correctness:** Fixed genome fraction correction silently disabled by NaN completeness values
- **Correctness:** Fixed last-write-wins bug destroying Kaiju taxonomy in `organize_genes_by_taxonomy.py`
- **Correctness:** Fixed `sys.exit(0)` masking plot generation failures in `reference_check_plot.py`
- **Correctness:** Removed duplicated single-end/long-read processing in REFERENCE_CHECK
- **Maintenance:** Aligned numpy version pins between Dockerfile and `install.sh`
- **Maintenance:** Unified Kaiju parsing across 5 scripts into shared `kaiju_utils.py`
- **Maintenance:** Unified name sanitization between `rename_contigs.py` and `organize_genes_by_taxonomy.py`
- **Maintenance:** Pinned bioinformatics tool versions in Dockerfile for reproducible builds
- **Portability:** Added macOS compatibility fixes (BSD `sed`, BSD `grep`, Bash 3.2 `declare -A`)
- **Portability:** Added `$HOME/.virall/databases` database path fallback for local runs
- **Portability:** Removed hardcoded user paths from example `asgsr_run_params.yaml`
- **Performance:** Switched intermediate `publishDir` modes from `copy` to `symlink`
- **Performance:** Replaced pointless `cp contigs scaffolds` with `touch` in assembly processes
- **Correctness:** Made `rename_contigs.py` idempotent on re-run
- **Maintenance:** Removed dead `run_plots.py` and stale references throughout codebase
- **Maintenance:** Collapsed 4 repeated auto-detection filter closures into a single helper
- Migrated to pure Nextflow pipeline; Python CLI removed
- Fix HPC scheduling: replace `params.threads` with `task.cpus` in all processes so SLURM/SGE allocate the correct number of CPUs
- Add resource directives (cpus/memory/time) to processes that were missing them (HOST_FILTER, GENOMAD, REFERENCE_CHECK, SC_*, etc.)
- Fix CheckV running on single thread regardless of configuration
- Extract MERGE_QUALITY inline Python into standalone `bin/merge_quality.py`
- Fix Medaka polishing command (`medaka_consensus` replaces non-existent `medaka_polish`)
- Add Medaka polishing for Nanopore long-only assemblies (PacBio HiFi skipped as already high-accuracy)
- Add geNomad support for RNA viruses and eukaryotic viruses alongside CheckV
- Add `genomad_db` path to Docker/Singularity profiles
- Singularity image now available on Sylabs Cloud (`library://brupav/virall/virall:latest`)

### v0.2.2
- Added configurable assembly parameters for Flye (min-overlap, genome-size hint)
- Exposed advanced assembly options via pipeline parameters

### v0.2.1
- Added `--filter` option for host genome filtering
- Documentation updates

### v0.2.0
- Added interactive plotting capabilities (sunburst charts, coverage heatmaps)
- Improved database path resolution for containerised environments

### v0.1.1
- Support for short, long, and single-cell sequencing reads
- Reference-guided assembly
- Viral assembly, classification, annotation, validation, and quantification
