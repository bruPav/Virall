graph LR
    %% Class Definitions
    classDef input fill:#f9f,stroke:#333,stroke-width:2px,color:#000
    classDef core fill:#e1f5fe,stroke:#01579b,stroke-width:2px,color:#000
    classDef quality fill:#fff3e0,stroke:#e65100,stroke-width:2px,color:#000
    classDef sc fill:#f3e5f5,stroke:#4a148c,stroke-width:2px,color:#000
    classDef output fill:#e8f5e9,stroke:#1b5e20,stroke-width:2px,color:#000

    subgraph "<b>Input Data</b>"
        FASTQ["<b>FASTQ Reads</b><br/>(Bulk or SC)"]:::input
    end

    subgraph "<b>Single-Cell Prep</b>"
        SC_EXTRACT["<b>SC_EXTRACT</b><br/>umi_tools"]:::sc
        SC_POOL["<b>SC_POOL</b><br/>cDNA Pooling"]:::sc
    end

    subgraph "<b>Core Processing</b>"
        PREPROCESS["<b>PREPROCESS</b><br/>fastp QC & Trim"]:::core
        HOST_FILTER["<b>HOST_FILTER</b><br/>Depletion"]:::core
        ASSEMBLE["<b>ASSEMBLE</b><br/>SPAdes / Flye"]:::core
        KAIJU["<b>KAIJU</b><br/>Taxonomy"]:::core
        FILTER_VIRAL["<b>FILTER_VIRAL</b><br/>Extract Contigs"]:::core
    end

    subgraph "<b>Quality & Validation</b>"
        VALIDATE["<b>CheckV</b><br/>(Phages)"]:::quality
        GENOMAD["<b>geNomad</b><br/>(RNA/Euk)"]:::quality
        MERGE_QUALITY["<b>MERGE_QUALITY</b>"]:::quality
    end

    subgraph "<b>Annotation & Stats</b>"
        ANNOTATE["<b>ANNOTATE</b><br/>VOGs"]:::output
        QUANTIFY["<b>QUANTIFY</b><br/>Read Mapping"]:::output
        PLOT["<b>VISUALIZATION</b><br/>Final Reports"]:::output
    end

    subgraph "<b>Single-Cell Trace-Back</b>"
        SC_MAP["<b>SC_MAP</b><br/>BWA-mem"]:::sc
        SC_MATRIX["<b>SC_MATRIX</b><br/>Count Matrix"]:::sc
    end

    %% Connections
    FASTQ -->|Single-Cell| SC_EXTRACT
    SC_EXTRACT --> SC_POOL
    SC_POOL --> PREPROCESS
    
    FASTQ -->|Bulk| PREPROCESS
    
    PREPROCESS --> HOST_FILTER
    HOST_FILTER --> ASSEMBLE
    ASSEMBLE --> KAIJU
    KAIJU --> FILTER_VIRAL

    FILTER_VIRAL --> VALIDATE
    FILTER_VIRAL --> GENOMAD
    VALIDATE & GENOMAD --> MERGE_QUALITY
    
    FILTER_VIRAL --> ANNOTATE
    FILTER_VIRAL --> QUANTIFY
    
    ANNOTATE & QUANTIFY & MERGE_QUALITY --> PLOT

    %% SC Traceback Connections
    SC_EXTRACT -.->|Tagged R2| SC_MAP
    FILTER_VIRAL --> SC_MAP
    SC_MAP --> SC_MATRIX