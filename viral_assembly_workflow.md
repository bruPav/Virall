# Viral Genome Assembly Tool - Process Flow Diagram

```mermaid
flowchart TD
    %% Input Data
    A[Sequencing Reads] --> B{Read Type Detection}
    B --> C[Paired-End Reads R1/R2]
    B --> D[Single-End Reads]
    B --> E[Long Reads ONT/PacBio]
    B --> F[Single-Cell Data 10X Genomics]

    %% Single-Cell Preprocessing
    F --> SC1[SC_EXTRACT_BARCODES\numi_tools extract]
    SC1 --> SC2[SC_POOL_READS\nPool cDNA as single-end]
    SC2 --> G

    %% Standard Preprocessing
    C --> G[PREPROCESS\nfastp / fastplong trimming]
    D --> G
    E --> G

    %% Host Filtering
    G --> HF{Host Genome\nProvided?}
    HF -->|Yes| HFF[HOST_FILTER\nminimap2 remove host reads]
    HF -->|No| U
    HFF --> U

    %% Assembly
    U[ASSEMBLE\nSPAdes / Flye] --> KJ

    %% Viral Identification
    KJ[KAIJU\nTaxonomic classification] --> FV
    FV[FILTER_VIRAL\nExtract viral contigs]

    %% Quality Assessment - Parallel
    FV --> VAL[VALIDATE\nCheckV - phages]
    FV --> GEN[GENOMAD\ngeNomad - RNA/eukaryotic viruses]
    VAL --> MQ[MERGE_QUALITY\nCombined quality report]
    GEN --> MQ

    %% Annotation
    FV --> ANN[ANNOTATE\nVOG annotation]
    ANN --> OG[ORGANIZE_GENES\nSort by taxonomy]

    %% Quantification and Plotting
    FV --> QT[QUANTIFY\nRead mapping and depth]
    QT --> PL[PLOT\nVisualization]
    MQ --> PL

    %% Optional Reference Check
    HFF -.->|Reference set| RC[REFERENCE_CHECK\nDetection and coverage]

    %% Single-Cell Trace-Back
    SC1 -->|Tagged R2| SCM[SC_MAP_VIRAL\nbwa mem to viral contigs]
    FV -->|Viral contigs| SCM
    SCM --> SCC[SC_COUNT_CELLS\numi_tools count per cell]
    SCC --> SCB[SC_BUILD_MATRIX\nMTX / barcodes / features]
    FV -->|Contigs and Kaiju| SCB

    %% Outputs
    PL --> OUT[Output Files]
    SCB --> OUT
    OUT --> O1[Viral Genomes FASTA]
    OUT --> O2[Quality Reports]
    OUT --> O3[Gene Annotations]
    OUT --> O4[Abundance Plots]
    OUT --> O5[Cell x Virus Matrix MTX]

    %% Styling
    classDef input fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    classDef process fill:#f3e5f5,stroke:#4a148c,stroke-width:2px
    classDef decision fill:#fff3e0,stroke:#e65100,stroke-width:2px
    classDef output fill:#e8f5e8,stroke:#2e7d32,stroke-width:2px
    classDef viral fill:#ffebee,stroke:#c62828,stroke-width:2px
    classDef singlecell fill:#fff8e1,stroke:#f57f17,stroke-width:2px

    class A,C,D,E,F input
    class G,HFF,U,KJ,QT,PL process
    class B,HF decision
    class O1,O2,O3,O4,O5 output
    class FV,VAL,GEN,MQ,ANN,OG viral
    class SC1,SC2,SCM,SCC,SCB singlecell
```

## Workflow Description

### 1. **Input Processing**
- Supports multiple read types: paired-end, single-end, long reads, and single-cell data
- Automatic read type detection and routing to appropriate preprocessing pipelines

### 2. **Preprocessing**
- **Short Reads**: fastp (QC, automatic adapter detection, trimming)
- **Long Reads**: fastplong (QC, automatic adapter detection, trimming for ONT/PacBio)
- **Single-Cell Pipeline (Pooled Mode)**: Extracts barcodes/UMIs from R1, pools cDNA (R2) as single-end reads for assembly
- Quality filtering with configurable thresholds (default: Phred 20 for short, Phred 7 for long)
- Optional host read filtering using minimap2 against a reference host genome

### 3. **Assembly**
- **Hybrid**: Combines short and long reads using SPAdes + Flye
- **Short-read only**: SPAdes for Illumina data (standard, metaviral, or rnaviral modes)
- **Long-read only**: Flye for ONT/PacBio data
- **RNA mode**: SPAdes `--rnaviral` for transcriptome assembly (RNA-seq and single-cell RNA-seq)

### 4. **Viral Identification**
- Kaiju classification against viral databases (direct contig analysis)
- Filtering based on confidence threshold (default: 0.8)

### 5. **Quality Assessment**
- **CheckV**: Optimal for bacteriophage genomes (completeness, contamination)
- **geNomad**: Optimal for RNA viruses and eukaryotic DNA viruses
- Results automatically merged based on Kaiju taxonomy classification

### 6. **Annotation**
- Gene prediction using Prodigal
- VOG (Viral Orthologous Groups) annotation
- Genes organized by viral taxonomy (family/genus/species)

### 7. **Quantification and Visualization**
- Read mapping to viral contigs (BWA-MEM for short, minimap2 for long)
- Coverage depth analysis
- Abundance plots, quality distribution, taxonomy sunburst, coverage plots

### 8. **Single-Cell Trace-Back**
- Barcoded reads mapped back to discovered viral contigs
- UMI-aware counting per cell per viral contig
- Output: Cell-by-virus MTX matrix compatible with Scanpy/Seurat
- Supports 10x Genomics v1, v2, and v3 chemistries (3' and 5')

## Key Features

- **Modular Design**: Each step can be run independently
- **Configurable**: YAML configuration file for parameter tuning
- **Multi-format Support**: Handles various sequencing technologies
- **Single-cell Ready**: 10x Genomics barcode extraction, pooled assembly, per-cell viral quantification
- **Dual Quality Assessment**: CheckV for phages, geNomad for RNA/eukaryotic viruses
- **Quality Focused**: Multiple validation and QC steps
- **Scalable**: Nextflow DSL2 with Docker/Singularity support
