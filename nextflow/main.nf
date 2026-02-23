/*
 * Virall Nextflow Pipeline – tools directly (no Virall Python wrappers).
 * Each process calls underlying tools (fastp, SPAdes, Flye, Kaiju, CheckV, Prodigal, BWA, etc.).
 * Rules (strategy, params) are driven by user params.
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameters (user requirements)
// ---------------------------------------------------------------------------
params.samples         = "${projectDir}/samples.csv"
params.outdir          = "results"
params.threads         = 8
params.memory          = "16G"
params.reference       = null
params.host_genome     = null   // optional: FASTA of host genome to filter out host reads (like Virall --filter)
params.rna_mode        = false
params.assembly_strategy = "auto"     // auto (detect from input), hybrid, short_only, long_only
params.min_contig_len  = 1000
params.quality_phred   = 20        // short reads (fastp)
params.quality_phred_long = 7     // long reads (fastplong); Virall default (lower for ONT/PacBio)
params.min_read_len    = 50       // short reads (fastp)
params.min_read_len_long = 1000   // long reads (fastplong); Virall default
params.long_read_tech  = "nanopore"  // "nanopore" or "pacbio" – affects SPAdes, Flye, minimap2 presets
params.flye_min_overlap = 1000       // Flye --min-overlap; lower values recover shorter viral genomes
params.flye_genome_size = null       // Flye --genome-size (e.g. "30k"); null = let Flye auto-estimate
params.iontorrent      = false       // set true for Ion Torrent short reads (adds --iontorrent to SPAdes)
params.metaviral_mode  = false       // set true to use metaviralSPAdes (--metaviral) instead of regular SPAdes; note: disables --trusted-contigs
// Database paths (set via -params-file or env VIRALL_DATABASE_DIR)
params.kaiju_db        = null
params.checkv_db       = null
params.genomad_db      = null        // geNomad database (for RNA viruses & eukaryotic viruses)
params.vog_db          = null

// Single-cell sequencing parameters
params.single_cell_mode = false       // Enable single-cell processing
params.sc_chemistry     = "10x_v3"  // 10x chemistry: "10x_v3" (default), "10x_v2", or "10x_v1"
params.sc_barcode_whitelist = null    // Optional: custom barcode whitelist file
params.sc_min_reads_per_cell = 100    // Minimum reads to keep a cell
params.sc_min_viral_umis = 1          // Minimum viral UMIs to call cell "infected"

// Resolve DB paths from env if not set (params.xxx are read-only; resolve in workflow)
def db_dir = System.getenv("VIRALL_DATABASE_DIR") ?: ""

// Expand ~ in outdir so "~/analysis/..." becomes /home/user/analysis/...
if (params.outdir?.toString()?.trim()?.startsWith('~')) {
    def home = System.getenv('HOME') ?: System.getProperty('user.home') ?: ''
    params.outdir = (home ?: '') + params.outdir.toString().trim().substring(1)
}

// ---------------------------------------------------------------------------
// Sample sheet
// ---------------------------------------------------------------------------
def expand_path = { String s ->
    if (!s?.trim()) return ''
    def p = s.trim()
    if (p.startsWith('~')) {
        def home = System.getenv('HOME') ?: System.getProperty('user.home') ?: ''
        p = (home ?: '') + p.substring(1)
    }
    return p
}
def to_file = { String s, Path stub -> (s?.trim()) ? file(expand_path(s)) : stub }
def stub_r1     = file("${projectDir}/.placeholder_r1")
def stub_r2     = file("${projectDir}/.placeholder_r2")
def stub_single = file("${projectDir}/.placeholder_single")
def stub_long   = file("${projectDir}/.placeholder_long")
def stub_sc_r1  = file("${projectDir}/.placeholder_sc_r1")
def stub_sc_r2  = file("${projectDir}/.placeholder_sc_r2")

Channel
    .fromPath(params.samples, checkIfExists: true)
    .splitCsv(header: true, strip: true)
    | map { row -> tuple(
        row.sample_id,
        to_file(row.read1, stub_r1),
        to_file(row.read2, stub_r2),
        to_file(row.single, stub_single),
        to_file(row.long, stub_long),
        to_file(row.sc_read1, stub_sc_r1),
        to_file(row.sc_read2, stub_sc_r2)
      )
    }
    | set { ch_samples }

// ---------------------------------------------------------------------------
// SINGLE-CELL: SC_EXTRACT_BARCODES – Extract cell barcodes and UMIs from 10x data
// ---------------------------------------------------------------------------
process SC_EXTRACT_BARCODES {
    tag "${sample_id}"
    label "single_cell"
    publishDir "${params.outdir}/${sample_id}/08_single_cell/barcoded_reads", mode: "copy"

    input:
    tuple val(sample_id), path(sc_read1), path(sc_read2)

    output:
    tuple val(sample_id), path("tagged_R1.fastq.gz"), path("tagged_R2.fastq.gz"), emit: tagged
    tuple val(sample_id), path("barcode_counts.tsv"), emit: barcode_counts

    script:
    // Select barcode pattern based on chemistry
    // v3 (3'/5'): 16bp barcode + 12bp UMI
    // v2 (3'/5'): 16bp barcode + 10bp UMI  
    // v1 (3'): 14bp barcode + 10bp UMI
    def bc_pattern
    switch(params.sc_chemistry) {
        case "10x_v3":
        case "10x_3prime_v3":
        case "10x_5prime_v3":
            bc_pattern = "CCCCCCCCCCCCCCCCNNNNNNNNNNNN"  // 16C + 12N
            break
        case "10x_v2":
        case "10x_3prime_v2":
        case "10x_5prime_v2":
            bc_pattern = "CCCCCCCCCCCCCCCCNNNNNNNNNN"    // 16C + 10N
            break
        case "10x_v1":
        case "10x_3prime_v1":
            bc_pattern = "CCCCCCCCCCCCCCNNNNNNNNNN"      // 14C + 10N
            break
        default:
            bc_pattern = "CCCCCCCCCCCCCCCCNNNNNNNNNNNN"  // Default to v3
    }
    def whitelist_opt = params.sc_barcode_whitelist ? "--whitelist ${params.sc_barcode_whitelist}" : ""
    """
    # Extract cell barcodes and UMIs from R1, tag R2 reads
    # umi_tools extract will add CB and UMI to read names
    # Chemistry: ${params.sc_chemistry} -> pattern: ${bc_pattern}
    umi_tools extract \
        --bc-pattern=${bc_pattern} \
        --stdin=${sc_read1} \
        --read2-in=${sc_read2} \
        --stdout=tagged_R1.fastq.gz \
        --read2-out=tagged_R2.fastq.gz \
        --log=extract.log \
        ${whitelist_opt} \
        --extract-method=string

    # Count barcodes for QC
    zcat tagged_R1.fastq.gz | awk 'NR%4==1 {split(\$1,a,"_"); print a[2]}' | sort | uniq -c | sort -rn > barcode_counts.tsv
    """
}

// ---------------------------------------------------------------------------
// SINGLE-CELL: SC_POOL_READS – Pool barcoded reads for assembly
// ---------------------------------------------------------------------------
process SC_POOL_READS {
    tag "${sample_id}"
    label "single_cell"

    input:
    tuple val(sample_id), path(tagged_r1), path(tagged_r2)

    output:
    tuple val(sample_id), path("pooled_single.fastq.gz"), emit: pooled

    script:
    """
    # For 10x single-cell, R1 only contains barcode+UMI (no biological sequence)
    # After umi_tools extract, R1 is essentially empty
    # We use R2 (cDNA) as single-end reads for assembly
    # The barcode info is preserved in the read names for later tracing
    cp ${tagged_r2} pooled_single.fastq.gz
    """
}

// ---------------------------------------------------------------------------
// PREPROCESS – fastp (short) / fastplong (long)
// ---------------------------------------------------------------------------
process PREPROCESS {
    tag "${sample_id}"
    label "preprocess"
    publishDir "${params.outdir}/${sample_id}/00_preprocess", mode: "copy"

    input:
    tuple val(sample_id), path(read1), path(read2), path(single), path(long_reads), path(sc_read1), path(sc_read2)

    output:
    tuple val(sample_id),
          path("preprocess_dir/trimmed_R1.fastq.gz"),
          path("preprocess_dir/trimmed_R2.fastq.gz"),
          path("preprocess_dir/trimmed_single.fastq.gz"),
          path("preprocess_dir/trimmed_long.fastq.gz"),
          path("preprocess_dir/fastp_pe.html"),
          path("preprocess_dir/fastp_pe.json"),
          path("preprocess_dir/fastp_single.html"),
          path("preprocess_dir/fastp_single.json"),
          path("preprocess_dir/fastplong.html"),
          path("preprocess_dir/fastplong.json"),
          emit: preprocessed

    script:
    def hasShort = (read1.name != ".placeholder_r1" && read1.size() > 0) || (single.name != ".placeholder_single" && single.size() > 0)
    def hasLong  = long_reads.name != ".placeholder_long" && long_reads.size() > 0
    def q = params.quality_phred
    def ml = params.min_read_len
    // Ion Torrent: skip poly-G trimming (Illumina artifact) and adapter detection (already stripped by Torrent Suite)
    def ion_opts = params.iontorrent ? "--disable_adapter_trimming --disable_trim_poly_g" : "--trim_poly_g"
    def pe_adapter = params.iontorrent ? "" : "--detect_adapter_for_pe"
    """
    mkdir -p preprocess_dir
    if [ -s "${read1}" ] && [ -s "${read2}" ]; then
      fastp -i ${read1} -I ${read2} -o preprocess_dir/trimmed_R1.fastq.gz -O preprocess_dir/trimmed_R2.fastq.gz \\
        --html preprocess_dir/fastp_pe.html --json preprocess_dir/fastp_pe.json \\
        --thread ${task.cpus} --qualified_quality_phred ${q} --length_required ${ml} \\
        ${pe_adapter} --cut_front --cut_tail --cut_mean_quality ${q} ${ion_opts}
    fi
    if [ -s "${single}" ] && [ "${single.name}" != ".placeholder_single" ]; then
      fastp -i ${single} -o preprocess_dir/trimmed_single.fastq.gz \\
        --html preprocess_dir/fastp_single.html --json preprocess_dir/fastp_single.json \\
        --thread ${task.cpus} --qualified_quality_phred ${q} --length_required ${ml} ${ion_opts}
    fi
    if [ -s "${long_reads}" ] && [ "${long_reads.name}" != ".placeholder_long" ]; then
      fastplong -i ${long_reads} -o preprocess_dir/trimmed_long.fastq.gz \\
        --html preprocess_dir/fastplong.html --json preprocess_dir/fastplong.json \\
        --thread ${task.cpus} --qualified_quality_phred ${params.quality_phred_long} \\
        --length_required ${params.min_read_len_long} \\
        --cut_front --cut_tail --cut_mean_quality ${params.quality_phred_long} --cut_window_size 10
    fi
    touch preprocess_dir/trimmed_R1.fastq.gz preprocess_dir/trimmed_R2.fastq.gz preprocess_dir/trimmed_single.fastq.gz preprocess_dir/trimmed_long.fastq.gz
    touch preprocess_dir/fastp_pe.html preprocess_dir/fastp_pe.json preprocess_dir/fastp_single.html preprocess_dir/fastp_single.json preprocess_dir/fastplong.html preprocess_dir/fastplong.json
    """
}

// ---------------------------------------------------------------------------
// HOST_FILTER – optional: remove reads mapping to host genome (minimap2 + samtools, like Virall --filter)
// ---------------------------------------------------------------------------
process HOST_FILTER {
    tag "${sample_id}"
    label "host_filter"
    publishDir "${params.outdir}/${sample_id}/00_preprocess/host_filtered", mode: "copy"

    input:
    tuple val(sample_id),
          path(trimmed_r1),
          path(trimmed_r2),
          path(trimmed_single),
          path(trimmed_long),
          path(qc_pe_html),
          path(qc_pe_json),
          path(qc_single_html),
          path(qc_single_json),
          path(qc_long_html),
          path(qc_long_json),
          path(host_ref)

    output:
    tuple val(sample_id),
          path("host_filter_dir/trimmed_R1.fastq.gz"),
          path("host_filter_dir/trimmed_R2.fastq.gz"),
          path("host_filter_dir/trimmed_single.fastq.gz"),
          path("host_filter_dir/trimmed_long.fastq.gz"),
          path("host_filter_dir/fastp_pe.html"),
          path("host_filter_dir/fastp_pe.json"),
          path("host_filter_dir/fastp_single.html"),
          path("host_filter_dir/fastp_single.json"),
          path("host_filter_dir/fastplong.html"),
          path("host_filter_dir/fastplong.json"),
          emit: host_filtered

    script:
    """
    mkdir -p host_filter_dir
    # Copy QC reports through
    cp "${qc_pe_html}" host_filter_dir/fastp_pe.html 2>/dev/null || touch host_filter_dir/fastp_pe.html
    cp "${qc_pe_json}" host_filter_dir/fastp_pe.json 2>/dev/null || touch host_filter_dir/fastp_pe.json
    cp "${qc_single_html}" host_filter_dir/fastp_single.html 2>/dev/null || touch host_filter_dir/fastp_single.html
    cp "${qc_single_json}" host_filter_dir/fastp_single.json 2>/dev/null || touch host_filter_dir/fastp_single.json
    cp "${qc_long_html}" host_filter_dir/fastplong.html 2>/dev/null || touch host_filter_dir/fastplong.html
    cp "${qc_long_json}" host_filter_dir/fastplong.json 2>/dev/null || touch host_filter_dir/fastplong.json

    # Paired-end: keep pairs where BOTH reads are unmapped (-f 12)
    if [ -s "${trimmed_r1}" ] && [ -s "${trimmed_r2}" ] && [ "${trimmed_r1.name}" != ".placeholder_r1" ]; then
      minimap2 -ax sr -t ${task.cpus} ${host_ref} ${trimmed_r1} ${trimmed_r2} 2>host_filter_dir/host_filter_pe.log | \\
        samtools sort -n - 2>/dev/null | samtools fastq -f 12 -1 host_filter_dir/R1.fq -2 host_filter_dir/R2.fq - 2>/dev/null || true
      if [ -s host_filter_dir/R1.fq ]; then
        gzip -c host_filter_dir/R1.fq > host_filter_dir/trimmed_R1.fastq.gz
        gzip -c host_filter_dir/R2.fq > host_filter_dir/trimmed_R2.fastq.gz
      else
        cp ${trimmed_r1} host_filter_dir/trimmed_R1.fastq.gz
        cp ${trimmed_r2} host_filter_dir/trimmed_R2.fastq.gz
      fi
    else
      cp ${trimmed_r1} host_filter_dir/trimmed_R1.fastq.gz 2>/dev/null || touch host_filter_dir/trimmed_R1.fastq.gz
      cp ${trimmed_r2} host_filter_dir/trimmed_R2.fastq.gz 2>/dev/null || touch host_filter_dir/trimmed_R2.fastq.gz
    fi

    # Single-end short
    if [ -s "${trimmed_single}" ] && [ "${trimmed_single.name}" != ".placeholder_single" ]; then
      minimap2 -ax sr -t ${task.cpus} ${host_ref} ${trimmed_single} 2>host_filter_dir/host_filter_single.log | \\
        samtools fastq -f 4 - > host_filter_dir/single.fq 2>/dev/null || true
      if [ -s host_filter_dir/single.fq ]; then
        gzip -c host_filter_dir/single.fq > host_filter_dir/trimmed_single.fastq.gz
      else
        cp ${trimmed_single} host_filter_dir/trimmed_single.fastq.gz
      fi
    else
      cp ${trimmed_single} host_filter_dir/trimmed_single.fastq.gz 2>/dev/null || touch host_filter_dir/trimmed_single.fastq.gz
    fi

    # Long reads (ONT or PacBio)
    MINIMAP2_LONG_PRESET=\$( [ "${params.long_read_tech}" = "pacbio" ] && echo "map-pb" || echo "map-ont" )
    if [ -s "${trimmed_long}" ] && [ "${trimmed_long.name}" != ".placeholder_long" ]; then
      minimap2 -ax \$MINIMAP2_LONG_PRESET -t ${task.cpus} ${host_ref} ${trimmed_long} 2>host_filter_dir/host_filter_long.log | \\
        samtools fastq -f 4 - > host_filter_dir/long.fq 2>/dev/null || true
      if [ -s host_filter_dir/long.fq ]; then
        gzip -c host_filter_dir/long.fq > host_filter_dir/trimmed_long.fastq.gz
      else
        cp ${trimmed_long} host_filter_dir/trimmed_long.fastq.gz
      fi
    else
      cp ${trimmed_long} host_filter_dir/trimmed_long.fastq.gz 2>/dev/null || touch host_filter_dir/trimmed_long.fastq.gz
    fi
    touch host_filter_dir/trimmed_R1.fastq.gz host_filter_dir/trimmed_R2.fastq.gz host_filter_dir/trimmed_single.fastq.gz host_filter_dir/trimmed_long.fastq.gz
    """
}

// ---------------------------------------------------------------------------
// ASSEMBLE – SPAdes (short/hybrid) and/or Flye (long)
// ---------------------------------------------------------------------------
process ASSEMBLE {
    tag "${sample_id}"
    label "assemble"
    publishDir "${params.outdir}/${sample_id}/01_assembly", mode: "copy", pattern: "contigs*"
    publishDir "${params.outdir}/${sample_id}/01_assembly", mode: "copy", pattern: "scaffolds"

    input:
    tuple val(sample_id),
          path(trimmed_r1),
          path(trimmed_r2),
          path(trimmed_single),
          path(trimmed_long),
          path(qc_pe_html),
          path(qc_pe_json),
          path(qc_single_html),
          path(qc_single_json),
          path(qc_long_html),
          path(qc_long_json)
    path(reference_file)

    output:
    tuple val(sample_id), path("contigs"), path("scaffolds"), path("preprocess_dir"), emit: assembled
    path("contigs_ref_guided.fasta"), optional: true

    script:
    def mem = task.memory ? task.memory.toGiga().toString() : params.memory.replaceAll(/[Gg]/, '')
    def ion = params.iontorrent ? '--iontorrent' : ''
    // --trusted-contigs is incompatible with --metaviral and --rnaviral; reference takes priority
    def has_ref = reference_file.name != '.placeholder_ref'
    // SPAdes only accepts .fa/.fasta extensions for --trusted-contigs; symlink if needed
    def ref_fasta = has_ref ? "reference.fasta" : ''
    def ref = has_ref ? "--trusted-contigs ${ref_fasta}" : ''
    def rna = (params.rna_mode && !has_ref) ? '--rnaviral' : ''
    def metaviral = (params.metaviral_mode && !has_ref) ? '--metaviral' : ''
    def minimap2_lr_preset = params.long_read_tech == "pacbio" ? "map-pb" : "map-ont"
    """
    mkdir -p assembly_dir preprocess_dir
    if [ -n "${ref_fasta}" ] && [ -f "${reference_file}" ]; then
      ln -sf "${reference_file}" "${ref_fasta}"
    fi

    if [ -n "${ref}" ]; then
      echo "NOTE: Reference genome provided — using --trusted-contigs for reference-guided assembly."
      if [ "${params.rna_mode}" = "true" ] || [ "${params.metaviral_mode}" = "true" ]; then
        echo "      Skipping --rnaviral/--metaviral (incompatible with --trusted-contigs)."
      fi
    fi
    
    cp "${trimmed_r1}" preprocess_dir/trimmed_R1.fastq.gz 2>/dev/null || true
    cp "${trimmed_r2}" preprocess_dir/trimmed_R2.fastq.gz 2>/dev/null || true
    cp "${trimmed_single}" preprocess_dir/trimmed_single.fastq.gz 2>/dev/null || true
    cp "${trimmed_long}" preprocess_dir/trimmed_long.fastq.gz 2>/dev/null || true

    # Detect what inputs are present
    HAS_SHORT=\$( [ -s preprocess_dir/trimmed_R1.fastq.gz ] && [ -s preprocess_dir/trimmed_R2.fastq.gz ] && echo 1 || echo 0 )
    HAS_SINGLE=\$( [ -s preprocess_dir/trimmed_single.fastq.gz ] && echo 1 || echo 0 )
    HAS_LONG=\$( [ -s preprocess_dir/trimmed_long.fastq.gz ] && echo 1 || echo 0 )

    if [ "\$HAS_SHORT" = "0" ] && [ "\$HAS_SINGLE" = "0" ] && [ "\$HAS_LONG" = "0" ]; then
      echo "ASSEMBLE: no reads remaining after preprocessing/host filtering for ${sample_id} – all reads may be host. Producing empty assembly."
      touch contigs scaffolds
    else
      # Determine strategy: auto-detect from inputs or use explicit setting
      USER_STRATEGY="${params.assembly_strategy}"
      if [ "\$USER_STRATEGY" = "auto" ] || [ -z "\$USER_STRATEGY" ]; then
        # Auto-detect: short/single only -> short_only; long only -> long_only; both -> hybrid
        if [ "\$HAS_LONG" = "1" ] && [ "\$HAS_SHORT" = "0" ] && [ "\$HAS_SINGLE" = "0" ]; then
          STRATEGY="long_only"
        elif [ "\$HAS_LONG" = "0" ]; then
          STRATEGY="short_only"
        else
          STRATEGY="hybrid"
        fi
        echo "Auto-detected assembly strategy: \$STRATEGY (short=\$HAS_SHORT, single=\$HAS_SINGLE, long=\$HAS_LONG)"
      else
        STRATEGY="\$USER_STRATEGY"
        echo "Using user-specified assembly strategy: \$STRATEGY"
      fi

      # Determine long-read flags based on technology (nanopore vs pacbio)
      LONG_READ_TECH="${params.long_read_tech}"
      if [ "\$LONG_READ_TECH" = "pacbio" ]; then
        SPADES_LONG_FLAG="--pacbio"
        FLYE_INPUT_FLAG="--pacbio-raw"
      else
        SPADES_LONG_FLAG="--nanopore"
        FLYE_INPUT_FLAG="--nano-raw"
      fi

      # Run assembler(s) based on strategy
      if [ "\$STRATEGY" = "short_only" ]; then
        if [ "\$HAS_SHORT" = "1" ] || [ "\$HAS_SINGLE" = "1" ]; then
          SPADES_OPTS="-o assembly_dir/spades -t ${task.cpus} -m ${mem} --only-assembler ${rna} ${ion} ${metaviral} ${ref}"
          [ -s preprocess_dir/trimmed_R1.fastq.gz ] && SPADES_OPTS="\$SPADES_OPTS -1 preprocess_dir/trimmed_R1.fastq.gz -2 preprocess_dir/trimmed_R2.fastq.gz"
          [ -s preprocess_dir/trimmed_single.fastq.gz ] && SPADES_OPTS="\$SPADES_OPTS -s preprocess_dir/trimmed_single.fastq.gz"
          spades.py \$SPADES_OPTS
          cp assembly_dir/spades/contigs.fasta contigs 2>/dev/null || cp assembly_dir/spades/transcripts.fasta contigs 2>/dev/null || touch contigs
          cp assembly_dir/spades/scaffolds.fasta scaffolds 2>/dev/null || cp assembly_dir/spades/hard_filtered_transcripts.fasta scaffolds 2>/dev/null || cp contigs scaffolds
        else
          echo "Warning: strategy=short_only but no short/single reads available"
        fi
      fi

      if [ "\$STRATEGY" = "hybrid" ]; then
        if [ "\$HAS_LONG" = "1" ]; then
          echo "Hybrid assembly: Flye --meta → Medaka (Nanopore only) → Polypolish → Pypolca"
          FLYE_EXTRA="--meta --min-overlap ${params.flye_min_overlap}"
          [ -n "${params.flye_genome_size ?: ''}" ] && FLYE_EXTRA="\$FLYE_EXTRA --genome-size ${params.flye_genome_size}"
          flye \$FLYE_INPUT_FLAG preprocess_dir/trimmed_long.fastq.gz --out-dir assembly_dir/flye \$FLYE_EXTRA -t ${task.cpus}
          cp assembly_dir/flye/assembly.fasta contigs 2>/dev/null || true

          # LR polishing: Medaka fixes systematic Nanopore errors (skip for PacBio)
          if [ "\$LONG_READ_TECH" = "nanopore" ] && [ -s contigs ]; then
            echo "Polishing Nanopore assembly with Medaka..."
            medaka_polish -i preprocess_dir/trimmed_long.fastq.gz -d contigs -o assembly_dir/medaka -t ${task.cpus} \
              && cp assembly_dir/medaka/consensus.fasta contigs \
              || echo "WARNING: Medaka polishing failed, continuing with unpolished assembly"
          fi

          # SR polishing: Polypolish (repeat-aware, uses all BWA alignments)
          if [ -s contigs ] && ([ "\$HAS_SHORT" = "1" ] || [ "\$HAS_SINGLE" = "1" ]); then
            if command -v polypolish &>/dev/null; then
              echo "Polishing assembly with short reads (Polypolish)..."
              bwa index contigs
              if [ "\$HAS_SHORT" = "1" ]; then
                bwa mem -t ${task.cpus} -a contigs preprocess_dir/trimmed_R1.fastq.gz > assembly_dir/alignments_1.sam
                bwa mem -t ${task.cpus} -a contigs preprocess_dir/trimmed_R2.fastq.gz > assembly_dir/alignments_2.sam
                polypolish polish contigs assembly_dir/alignments_1.sam assembly_dir/alignments_2.sam > assembly_dir/polished.fasta \
                  && cp assembly_dir/polished.fasta contigs \
                  || echo "WARNING: Polypolish failed, using previous assembly"
                rm -f assembly_dir/alignments_1.sam assembly_dir/alignments_2.sam
              elif [ "\$HAS_SINGLE" = "1" ]; then
                bwa mem -t ${task.cpus} -a contigs preprocess_dir/trimmed_single.fastq.gz > assembly_dir/alignments.sam
                polypolish polish contigs assembly_dir/alignments.sam > assembly_dir/polished.fasta \
                  && cp assembly_dir/polished.fasta contigs \
                  || echo "WARNING: Polypolish failed, using previous assembly"
                rm -f assembly_dir/alignments.sam
              fi
            else
              echo "WARNING: Polypolish not found – skipping. Install with: conda install -c bioconda polypolish"
            fi
          fi

          # Final polish: Pypolca (k-mer-based cleanup of residual SNPs/indels)
          if [ -s contigs ] && [ "\$HAS_SHORT" = "1" ]; then
            if command -v pypolca &>/dev/null; then
              echo "Final polishing with Pypolca..."
              pypolca run -a contigs \
                -1 preprocess_dir/trimmed_R1.fastq.gz \
                -2 preprocess_dir/trimmed_R2.fastq.gz \
                -t ${task.cpus} --careful \
                -o assembly_dir/pypolca 2>&1 \
                && { PYPOLCA_OUT=\$(ls assembly_dir/pypolca/*corrected*.fasta 2>/dev/null | head -1); \
                     [ -s "\$PYPOLCA_OUT" ] && cp "\$PYPOLCA_OUT" contigs; } \
                || echo "WARNING: Pypolca failed, using Polypolish-polished assembly"
            else
              echo "WARNING: Pypolca not found – skipping final polish. Install with: pip install pypolca"
            fi
          fi

          # Reference-guided consensus from long reads
          if [ -n "${ref_fasta}" ] && [ -f "${ref_fasta}" ]; then
            echo "Reference-guided long-read consensus assembly..."
            minimap2 -t ${task.cpus} -ax ${minimap2_lr_preset} ${ref_fasta} \
              preprocess_dir/trimmed_long.fastq.gz 2>/dev/null | \
              samtools sort -@ ${task.cpus} -o assembly_dir/ref_aligned.bam -
            samtools index assembly_dir/ref_aligned.bam

            if [ "\$LONG_READ_TECH" = "nanopore" ]; then
              medaka_polish -i preprocess_dir/trimmed_long.fastq.gz \
                -d ${ref_fasta} -o assembly_dir/ref_medaka -t ${task.cpus} \
                && cp assembly_dir/ref_medaka/consensus.fasta assembly_dir/ref_consensus.fasta \
                || samtools consensus assembly_dir/ref_aligned.bam \
                     -o assembly_dir/ref_consensus.fasta --show-ins no -a
            else
              samtools consensus assembly_dir/ref_aligned.bam \
                -o assembly_dir/ref_consensus.fasta --show-ins no -a
            fi

            if [ -s assembly_dir/ref_consensus.fasta ] && [ -s contigs ]; then
              cat contigs assembly_dir/ref_consensus.fasta > assembly_dir/combined.fasta
              cp assembly_dir/combined.fasta contigs
            elif [ -s assembly_dir/ref_consensus.fasta ]; then
              cp assembly_dir/ref_consensus.fasta contigs
            fi

            cp assembly_dir/ref_consensus.fasta contigs_ref_guided.fasta 2>/dev/null || true
          fi

          cp contigs scaffolds 2>/dev/null || true
        else
          echo "Warning: strategy=hybrid but no long reads available, falling back to short-read assembly"
          if [ "\$HAS_SHORT" = "1" ] || [ "\$HAS_SINGLE" = "1" ]; then
            SPADES_OPTS="-o assembly_dir/spades -t ${task.cpus} -m ${mem} --only-assembler ${rna} ${ion} ${metaviral} ${ref}"
            [ -s preprocess_dir/trimmed_R1.fastq.gz ] && SPADES_OPTS="\$SPADES_OPTS -1 preprocess_dir/trimmed_R1.fastq.gz -2 preprocess_dir/trimmed_R2.fastq.gz"
            [ -s preprocess_dir/trimmed_single.fastq.gz ] && SPADES_OPTS="\$SPADES_OPTS -s preprocess_dir/trimmed_single.fastq.gz"
            spades.py \$SPADES_OPTS
            cp assembly_dir/spades/contigs.fasta contigs 2>/dev/null || cp assembly_dir/spades/transcripts.fasta contigs 2>/dev/null || touch contigs
            cp assembly_dir/spades/scaffolds.fasta scaffolds 2>/dev/null || cp assembly_dir/spades/hard_filtered_transcripts.fasta scaffolds 2>/dev/null || cp contigs scaffolds
          fi
        fi
      fi

      if [ "\$STRATEGY" = "long_only" ]; then
        if [ "\$HAS_LONG" = "1" ]; then
          FLYE_EXTRA="--meta --min-overlap ${params.flye_min_overlap}"
          [ -n "${params.flye_genome_size ?: ''}" ] && FLYE_EXTRA="\$FLYE_EXTRA --genome-size ${params.flye_genome_size}"
          flye \$FLYE_INPUT_FLAG preprocess_dir/trimmed_long.fastq.gz --out-dir assembly_dir/flye \$FLYE_EXTRA -t ${task.cpus}
          cp assembly_dir/flye/assembly.fasta contigs 2>/dev/null || true

          # Polish Nanopore assemblies with Medaka (PacBio HiFi is already high-accuracy)
          if [ "\$LONG_READ_TECH" = "nanopore" ] && [ -s contigs ]; then
            echo "Polishing Nanopore assembly with Medaka..."
            medaka_polish -i preprocess_dir/trimmed_long.fastq.gz -d contigs -o assembly_dir/medaka -t ${task.cpus} \
              && cp assembly_dir/medaka/consensus.fasta contigs \
              || echo "WARNING: Medaka polishing failed, using unpolished Flye assembly"
          fi

          # Reference-guided consensus from long reads
          if [ -n "${ref_fasta}" ] && [ -f "${ref_fasta}" ]; then
            echo "Reference-guided long-read consensus assembly..."
            minimap2 -t ${task.cpus} -ax ${minimap2_lr_preset} ${ref_fasta} \
              preprocess_dir/trimmed_long.fastq.gz 2>/dev/null | \
              samtools sort -@ ${task.cpus} -o assembly_dir/ref_aligned.bam -
            samtools index assembly_dir/ref_aligned.bam

            if [ "\$LONG_READ_TECH" = "nanopore" ]; then
              medaka_polish -i preprocess_dir/trimmed_long.fastq.gz \
                -d ${ref_fasta} -o assembly_dir/ref_medaka -t ${task.cpus} \
                && cp assembly_dir/ref_medaka/consensus.fasta assembly_dir/ref_consensus.fasta \
                || samtools consensus assembly_dir/ref_aligned.bam \
                     -o assembly_dir/ref_consensus.fasta --show-ins no -a
            else
              samtools consensus assembly_dir/ref_aligned.bam \
                -o assembly_dir/ref_consensus.fasta --show-ins no -a
            fi

            if [ -s assembly_dir/ref_consensus.fasta ] && [ -s contigs ]; then
              cat contigs assembly_dir/ref_consensus.fasta > assembly_dir/combined.fasta
              cp assembly_dir/combined.fasta contigs
            elif [ -s assembly_dir/ref_consensus.fasta ]; then
              cp assembly_dir/ref_consensus.fasta contigs
            fi

            cp assembly_dir/ref_consensus.fasta contigs_ref_guided.fasta 2>/dev/null || true
          fi

          cp contigs scaffolds 2>/dev/null || true
        else
          echo "Warning: strategy=long_only but no long reads available"
        fi
      fi

      [ -f contigs ] || touch contigs scaffolds
    fi
    """
}

// ---------------------------------------------------------------------------
// KAIJU – classify all contigs (contigs.fasta used for viral filtering)
// ---------------------------------------------------------------------------
process KAIJU {
    tag "${sample_id}"
    label "kaiju"
    publishDir "${params.outdir}/${sample_id}/03_classifications", mode: "copy"

    input:
    tuple val(sample_id), path(contigs), path(scaffolds), path(preprocess_dir)
    val(kaiju_db)

    output:
    tuple val(sample_id), path("kaiju_dir"), path(contigs), path(scaffolds), path(preprocess_dir), emit: kaiju_done

    script:
    """
  mkdir -p kaiju_dir
  CONTIG_COUNT=\$(grep -c "^>" contigs 2>/dev/null || true)
  CONTIG_COUNT=\${CONTIG_COUNT:-0}
  if [ "\$CONTIG_COUNT" -eq 0 ]; then
    echo "KAIJU: skipping – no contigs to classify for ${sample_id}"
    touch kaiju_dir/kaiju_results.tsv kaiju_dir/kaiju_results_with_names.tsv
  else
    FMI=\$(find ${kaiju_db} -name "*.fmi" 2>/dev/null | head -1)
    [ -z "\$FMI" ] && { echo "Kaiju .fmi not found under ${kaiju_db}"; exit 1; }
    NODES=\$(find ${kaiju_db} -name "*nodes.dmp" 2>/dev/null | head -1)
    NAMES=\$(find ${kaiju_db} -name "*names.dmp" 2>/dev/null | head -1)
    [ -z "\$NODES" ] && { echo "Kaiju nodes.dmp not found under ${kaiju_db}"; exit 1; }
    [ -z "\$NAMES" ] && { echo "Kaiju names.dmp not found under ${kaiju_db}"; exit 1; }
    kaiju -t "\$NODES" -f "\$FMI" -i contigs -o kaiju_dir/kaiju_results.tsv -z ${task.cpus}
    kaiju-addTaxonNames -t "\$NODES" -n "\$NAMES" \\
        -i kaiju_dir/kaiju_results.tsv -o kaiju_dir/kaiju_results_with_names.tsv \\
        -r superkingdom,phylum,class,order,family,genus,species
  fi
    """
}

// ---------------------------------------------------------------------------
// FILTER_VIRAL – keep classified (non-Unclassified) contigs as viral
// ---------------------------------------------------------------------------
process FILTER_VIRAL {
    tag "${sample_id}"
    label "filter_viral"
    publishDir "${params.outdir}/${sample_id}/02_viral_contigs", mode: "copy"

    input:
    tuple val(sample_id), path(kaiju_dir), path(contigs), path(scaffolds), path(preprocess_dir)
    path(extract_script)

    output:
    tuple val(sample_id), path("viral_contigs.fasta"), path(kaiju_dir), path(preprocess_dir), emit: viral

    script:
    """
    python ${extract_script} \\
      --kaiju kaiju_dir/kaiju_results_with_names.tsv \\
      --fasta contigs \\
      --min-len ${params.min_contig_len} \\
      --out viral_contigs.fasta
    if [ ! -s viral_contigs.fasta ]; then
      echo "NOTE: No viral contigs >= ${params.min_contig_len} bp found for ${sample_id}. Downstream steps will report empty results."
      touch viral_contigs.fasta
    fi
    """
}

// ---------------------------------------------------------------------------
// RENAME_CONTIGS – Rename viral contigs by lowest Kaiju taxonomy level
// ---------------------------------------------------------------------------
process RENAME_CONTIGS {
    tag "${sample_id}"
    label "rename_contigs"
    publishDir "${params.outdir}/${sample_id}/02_viral_contigs", mode: "copy", pattern: "{viral_contigs.fasta,name_mapping.tsv}"

    input:
    tuple val(sample_id), path("input_viral_contigs.fasta"), path("input_kaiju_dir"), path(preprocess_dir)

    output:
    tuple val(sample_id), path("viral_contigs.fasta"), path("kaiju_dir"), path(preprocess_dir), emit: renamed

    script:
    """
    mkdir -p kaiju_dir
    if [ -s input_viral_contigs.fasta ]; then
      rename_contigs.py \\
        --fasta input_viral_contigs.fasta \\
        --kaiju-dir input_kaiju_dir \\
        --out-fasta viral_contigs.fasta \\
        --out-kaiju-dir kaiju_dir \\
        --out-mapping name_mapping.tsv
    else
      echo "RENAME_CONTIGS: no viral contigs to rename for ${sample_id}"
      touch viral_contigs.fasta name_mapping.tsv
      cp input_kaiju_dir/kaiju_results_with_names.tsv kaiju_dir/ 2>/dev/null || touch kaiju_dir/kaiju_results_with_names.tsv
      cp input_kaiju_dir/kaiju_results.tsv kaiju_dir/ 2>/dev/null || touch kaiju_dir/kaiju_results.tsv
    fi
    """
}

// ---------------------------------------------------------------------------
// VALIDATE – CheckV on viral contigs
// ---------------------------------------------------------------------------
process VALIDATE {
    tag "${sample_id}"
    label "validate"
    publishDir "${params.outdir}/${sample_id}/04_quality_assessment", mode: "copy"

    input:
    tuple val(sample_id), path(viral_contigs), path(kaiju_dir), path(preprocess_dir)
    val(checkv_db)

    output:
    tuple val(sample_id), path("checkv_dir"), path(viral_contigs), path(kaiju_dir), path(preprocess_dir), emit: validated

    script:
    """
    mkdir -p checkv_dir
    CHECKV_D=""
    for cand in "${checkv_db}/checkv-db-v1.5" "${checkv_db}/checkv-db-v1.0" "${checkv_db}"; do
      if [ -d "\$cand/genome_db" ] && [ -f "\$cand/genome_db/checkv_reps.dmnd" ]; then
        CHECKV_D="\$cand"
        break
      fi
    done
    if [ -z "\$CHECKV_D" ]; then
      echo "CheckV: no valid database found under ${checkv_db}. Expected genome_db/checkv_reps.dmnd (e.g. from checkv-db-v1.5)."
      echo "Rebuild the container so the CheckV DB download completes, or set checkv_db in run_params to a path with the DB."
      exit 1
    fi
    # Guard: if contigs are empty or too short for gene prediction, CheckV's
    # hmmsearch will fail on the empty proteins file.  Produce a valid but
    # empty quality_summary.tsv so downstream steps continue gracefully.
    CONTIG_COUNT=\$(grep -c "^>" viral_contigs.fasta 2>/dev/null || true)
    CONTIG_COUNT=\${CONTIG_COUNT:-0}
    if [ "\$CONTIG_COUNT" -eq 0 ]; then
      echo "CheckV: skipping – no viral contigs found for ${sample_id}"
      printf "contig_id\\tcontig_length\\tprovirus\\tproviral_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tmiuvig_quality\\tcompleteness\\tcompleteness_method\\tcontamination\\tkmer_freq\\twarnings\\n" > checkv_dir/quality_summary.tsv
    else
      TOTAL_BP=\$(grep -v "^>" viral_contigs.fasta | tr -d '\\n' | wc -c)
      if [ "\$TOTAL_BP" -lt 200 ]; then
        echo "CheckV: skipping – viral contigs too short (\${TOTAL_BP} bp) for meaningful quality assessment"
        printf "contig_id\\tcontig_length\\tprovirus\\tproviral_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tmiuvig_quality\\tcompleteness\\tcompleteness_method\\tcontamination\\tkmer_freq\\twarnings\\n" > checkv_dir/quality_summary.tsv
      else
        checkv end_to_end viral_contigs.fasta checkv_dir -d "\$CHECKV_D" -t ${task.cpus}
      fi
    fi
    """
}

// ---------------------------------------------------------------------------
// GENOMAD – geNomad for RNA viruses & eukaryotic viruses
// ---------------------------------------------------------------------------
process GENOMAD {
    tag "${sample_id}"
    label "genomad"
    publishDir "${params.outdir}/${sample_id}/04_quality_assessment/genomad", mode: "copy"

    input:
    tuple val(sample_id), path(viral_contigs), path(kaiju_dir), path(preprocess_dir)
    val(genomad_db)

    output:
    tuple val(sample_id), path("genomad_out"), path(viral_contigs), path(kaiju_dir), path(preprocess_dir), emit: genomad_done

    script:
    """
    mkdir -p genomad_out

    # Find geNomad database
    GENOMAD_D=""
    for cand in "${genomad_db}/genomad_db" "${genomad_db}"; do
      if [ -d "\$cand" ] && [ -f "\$cand/genomad_db" ] || [ -f "\$cand/virus_hallmark_annotation.txt" ] || [ -d "\$cand/mmseqs2" ]; then
        GENOMAD_D="\$cand"
        break
      fi
    done

    # Skip if no viral contigs
    CONTIG_COUNT=\$(grep -c "^>" ${viral_contigs} 2>/dev/null || true)
    CONTIG_COUNT=\${CONTIG_COUNT:-0}
    if [ "\$CONTIG_COUNT" -eq 0 ]; then
      echo "geNomad: skipping – no viral contigs found for ${sample_id}"
      touch genomad_out/virus_summary.tsv
    elif [ -z "\$GENOMAD_D" ]; then
      echo "WARNING: geNomad database not found at ${genomad_db}. Skipping geNomad analysis."
      echo "To use geNomad, download the database: genomad download-database genomad_db"
      touch genomad_out/virus_summary.tsv
      touch genomad_out/skipped.txt
    else
      # Run geNomad end-to-end (with neural network classification enabled)
      genomad end-to-end \\
        --cleanup \\
        --splits ${task.cpus} \\
        ${viral_contigs} \\
        genomad_out \\
        "\$GENOMAD_D"

      # Move results to expected location (geNomad creates subdir with input filename)
      if [ -d genomad_out/${viral_contigs.baseName}_summary ]; then
        mv genomad_out/${viral_contigs.baseName}_summary/* genomad_out/ 2>/dev/null || true
      fi
      if [ -d genomad_out/viral_contigs_summary ]; then
        mv genomad_out/viral_contigs_summary/* genomad_out/ 2>/dev/null || true
      fi

      # Rename to standard name if needed
      for f in genomad_out/*_virus_summary.tsv; do
        [ -f "\$f" ] && mv "\$f" genomad_out/virus_summary.tsv 2>/dev/null || true
        break
      done
    fi

    # Ensure output file exists
    touch genomad_out/virus_summary.tsv
    """
}

// ---------------------------------------------------------------------------
// MERGE_QUALITY – Combine CheckV (phages) + geNomad (other viruses) results
// ---------------------------------------------------------------------------
process MERGE_QUALITY {
    tag "${sample_id}"
    label "merge_quality"
    publishDir "${params.outdir}/${sample_id}/04_quality_assessment", mode: "copy"

    input:
    tuple val(sample_id), path(checkv_dir), path(genomad_dir), path(viral_contigs), path(kaiju_dir), path(preprocess_dir)
    path(merge_script)

    output:
    tuple val(sample_id), path("merged_quality"), path(viral_contigs), path(kaiju_dir), path(preprocess_dir), emit: merged

    script:
    """
    python ${merge_script} \\
      --checkv-dir ${checkv_dir} \\
      --genomad-dir ${genomad_dir} \\
      --kaiju-dir ${kaiju_dir} \\
      --out-dir merged_quality
    """
}

// ---------------------------------------------------------------------------
// ANNOTATE – prodigal-gv + hmmscan (VOG)
// ---------------------------------------------------------------------------
process ANNOTATE {
    tag "${sample_id}"
    label "annotate"
    publishDir "${params.outdir}/${sample_id}/05_gene_predictions", mode: "copy"

    input:
    tuple val(sample_id), path(viral_contigs), path(kaiju_dir), path(preprocess_dir)
    val(vog_db)

    output:
    tuple val(sample_id), path("annotation_dir"), emit: annotated

    script:
    """
    mkdir -p annotation_dir
    CONTIG_COUNT=\$(grep -c "^>" viral_contigs.fasta 2>/dev/null || true)
    CONTIG_COUNT=\${CONTIG_COUNT:-0}
    if [ "\$CONTIG_COUNT" -eq 0 ]; then
      echo "ANNOTATE: skipping – no viral contigs found for ${sample_id}"
      touch annotation_dir/proteins.faa
    else
      prodigal-gv -i viral_contigs.fasta -a annotation_dir/proteins.faa -p meta -q
      VOG_HMM=\$(find ${vog_db} -maxdepth 1 \\( -name 'vog_all.hmm' -o -name 'vog.hmm' \\) 2>/dev/null | head -1)
      if [ -n "\$VOG_HMM" ] && [ -f "\$VOG_HMM" ]; then
        hmmscan --cpu ${task.cpus} -o annotation_dir/vog_out.txt --domtblout annotation_dir/vog_domains.txt "\$VOG_HMM" annotation_dir/proteins.faa
      fi
    fi
    """
}

// ---------------------------------------------------------------------------
// ORGANIZE_GENES – Organize genes by Kaiju taxonomy + VOG functional annotation
// ---------------------------------------------------------------------------
process ORGANIZE_GENES {
    tag "${sample_id}"
    label "annotate"
    publishDir "${params.outdir}/${sample_id}/05_gene_predictions", mode: "copy"

    input:
    tuple val(sample_id), path(annotation_dir), path(viral_contigs), path(kaiju_dir)
    val(vog_db)

    output:
    tuple val(sample_id), path("by_taxonomy"), emit: organized_genes

    script:
    """
    mkdir -p by_taxonomy
    if [ -s ${annotation_dir}/proteins.faa ] && [ -s ${viral_contigs} ]; then
      # Build VOG annotation flags if metadata files exist
      VOG_FLAGS=""
      if [ -f "${annotation_dir}/vog_domains.txt" ] && [ -s "${annotation_dir}/vog_domains.txt" ]; then
        VOG_FLAGS="--vog-domains ${annotation_dir}/vog_domains.txt"
      fi
      if [ -f "${vog_db}/vog.annotations.tsv" ]; then
        VOG_FLAGS="\$VOG_FLAGS --vog-annotations ${vog_db}/vog.annotations.tsv"
      fi
      if [ -f "${vog_db}/vog.virusonly.tsv" ]; then
        VOG_FLAGS="\$VOG_FLAGS --vog-virusonly ${vog_db}/vog.virusonly.tsv"
      fi
      if [ -f "${vog_db}/vogdb.functional_categories.txt" ]; then
        VOG_FLAGS="\$VOG_FLAGS --vog-categories ${vog_db}/vogdb.functional_categories.txt"
      fi

      organize_genes_by_taxonomy.py \
          --proteins ${annotation_dir}/proteins.faa \
          --contigs ${viral_contigs} \
          --kaiju ${kaiju_dir}/kaiju_results_with_names.tsv \
          \$VOG_FLAGS \
          --outdir by_taxonomy
    else
      echo "ORGANIZE_GENES: skipping – no viral contigs or proteins found for ${sample_id}"
      touch by_taxonomy/gene_annotation.tsv
      ln -sf gene_annotation.tsv by_taxonomy/gene_taxonomy_mapping.tsv
    fi
    """
}

// ---------------------------------------------------------------------------
// QUANTIFY – BWA/minimap2 + samtools depth
// ---------------------------------------------------------------------------
process QUANTIFY {
    tag "${sample_id}"
    label "quantify"
    publishDir "${params.outdir}/${sample_id}/06_quantification", mode: "copy"

    input:
    tuple val(sample_id), path(viral_contigs), path(kaiju_dir), path(preprocess_dir)

    output:
    tuple val(sample_id), path("quant_dir"), emit: quantified

    script:
    """
    mkdir -p quant_dir
    CONTIG_COUNT=\$(grep -c "^>" viral_contigs.fasta 2>/dev/null || true)
    CONTIG_COUNT=\${CONTIG_COUNT:-0}
    if [ "\$CONTIG_COUNT" -eq 0 ]; then
      echo "QUANTIFY: skipping – no viral contigs found for ${sample_id}"
      touch quant_dir/mapped.bam quant_dir/mapped.bam.bai quant_dir/depth.txt
    else
      bwa index -p quant_dir/idx viral_contigs.fasta
      BAMS_TO_MERGE=""
      if [ -s preprocess_dir/trimmed_R1.fastq.gz ] && [ -s preprocess_dir/trimmed_R2.fastq.gz ]; then
        bwa mem -t ${task.cpus} quant_dir/idx preprocess_dir/trimmed_R1.fastq.gz preprocess_dir/trimmed_R2.fastq.gz | samtools sort -o quant_dir/mapped_pe.bam -
        BAMS_TO_MERGE="\$BAMS_TO_MERGE quant_dir/mapped_pe.bam"
      fi
      if [ -s preprocess_dir/trimmed_single.fastq.gz ]; then
        bwa mem -t ${task.cpus} quant_dir/idx preprocess_dir/trimmed_single.fastq.gz | samtools sort -o quant_dir/mapped_single.bam -
        BAMS_TO_MERGE="\$BAMS_TO_MERGE quant_dir/mapped_single.bam"
      fi
      HAS_SHORT=\$([ -n "\$BAMS_TO_MERGE" ] && echo "1" || echo "0")
      if [ "\$HAS_SHORT" = "0" ] && [ -s preprocess_dir/trimmed_long.fastq.gz ]; then
        MINIMAP2_LONG_PRESET=\$( [ "${params.long_read_tech}" = "pacbio" ] && echo "map-pb" || echo "map-ont" )
        minimap2 -t ${task.cpus} -ax \$MINIMAP2_LONG_PRESET viral_contigs.fasta preprocess_dir/trimmed_long.fastq.gz | samtools sort -o quant_dir/mapped_long.bam -
        BAMS_TO_MERGE="\$BAMS_TO_MERGE quant_dir/mapped_long.bam"
      fi
      BAM_COUNT=\$(echo \$BAMS_TO_MERGE | wc -w)
      if [ "\$BAM_COUNT" -gt 1 ]; then
        samtools merge -f quant_dir/mapped.bam \$BAMS_TO_MERGE
      elif [ "\$BAM_COUNT" -eq 1 ]; then
        mv \$BAMS_TO_MERGE quant_dir/mapped.bam
      else
        touch quant_dir/mapped.bam
      fi
      samtools index quant_dir/mapped.bam
      samtools depth -a quant_dir/mapped.bam > quant_dir/depth.txt 2>/dev/null || true
    fi
    """
}

// ---------------------------------------------------------------------------
// REFERENCE_CHECK – optional: map reads to reference, report detection, plot coverage
// ---------------------------------------------------------------------------
process REFERENCE_CHECK {
    tag "${sample_id}"
    label "reference_check"
    publishDir "${params.outdir}/${sample_id}/08_reference_check", mode: "copy"

    input:
    tuple val(sample_id),
          path(trimmed_r1),
          path(trimmed_r2),
          path(trimmed_single),
          path(trimmed_long),
          path(qc_pe_html),
          path(qc_pe_json),
          path(qc_single_html),
          path(qc_single_json),
          path(qc_long_html),
          path(qc_long_json)
    path(reference)

    output:
    tuple val(sample_id), path("ref_check_dir"), emit: ref_checked

    script:
    def minimap2_preset = params.long_read_tech == "pacbio" ? "map-pb" : "map-ont"
    """
    mkdir -p ref_check_dir

    # Index reference
    bwa index -p ref_check_dir/ref_idx ${reference} 2>/dev/null || true

    # Map reads to reference and collect stats
    MAPPED_READS=0
    TOTAL_READS=0

    if [ -s "${trimmed_r1}" ] && [ -s "${trimmed_r2}" ] && [ "${trimmed_r1.name}" != ".placeholder_r1" ]; then
      bwa mem -t ${task.cpus} ref_check_dir/ref_idx ${trimmed_r1} ${trimmed_r2} 2>/dev/null | \\
        samtools sort -o ref_check_dir/mapped_pe.bam - 2>/dev/null
      samtools index ref_check_dir/mapped_pe.bam 2>/dev/null || true
      PE_MAPPED=\$(samtools view -c -F 4 ref_check_dir/mapped_pe.bam 2>/dev/null || echo 0)
      PE_TOTAL=\$(samtools view -c ref_check_dir/mapped_pe.bam 2>/dev/null || echo 0)
      MAPPED_READS=\$((MAPPED_READS + PE_MAPPED))
      TOTAL_READS=\$((TOTAL_READS + PE_TOTAL))
    fi

    if [ -s "${trimmed_single}" ] && [ "${trimmed_single.name}" != ".placeholder_single" ]; then
      bwa mem -t ${task.cpus} ref_check_dir/ref_idx ${trimmed_single} 2>/dev/null | \\
        samtools sort -o ref_check_dir/mapped_single.bam - 2>/dev/null
      samtools index ref_check_dir/mapped_single.bam 2>/dev/null || true
      SINGLE_MAPPED=\$(samtools view -c -F 4 ref_check_dir/mapped_single.bam 2>/dev/null || echo 0)
      SINGLE_TOTAL=\$(samtools view -c ref_check_dir/mapped_single.bam 2>/dev/null || echo 0)
      MAPPED_READS=\$((MAPPED_READS + SINGLE_MAPPED))
      TOTAL_READS=\$((TOTAL_READS + SINGLE_TOTAL))
    fi

    if [ -s "${trimmed_long}" ] && [ "${trimmed_long.name}" != ".placeholder_long" ]; then
      minimap2 -t ${task.cpus} -ax ${minimap2_preset} ${reference} ${trimmed_long} 2>/dev/null | \\
        samtools sort -o ref_check_dir/mapped_long.bam - 2>/dev/null
      samtools index ref_check_dir/mapped_long.bam 2>/dev/null || true
      LONG_MAPPED=\$(samtools view -c -F 4 ref_check_dir/mapped_long.bam 2>/dev/null || echo 0)
      LONG_TOTAL=\$(samtools view -c ref_check_dir/mapped_long.bam 2>/dev/null || echo 0)
      MAPPED_READS=\$((MAPPED_READS + LONG_MAPPED))
      TOTAL_READS=\$((TOTAL_READS + LONG_TOTAL))
    fi

    # Merge BAMs if multiple exist
    BAM_FILES=""
    [ -f ref_check_dir/mapped_pe.bam ] && BAM_FILES="\$BAM_FILES ref_check_dir/mapped_pe.bam"
    [ -f ref_check_dir/mapped_single.bam ] && BAM_FILES="\$BAM_FILES ref_check_dir/mapped_single.bam"
    [ -f ref_check_dir/mapped_long.bam ] && BAM_FILES="\$BAM_FILES ref_check_dir/mapped_long.bam"

    if [ -n "\$BAM_FILES" ]; then
      samtools merge -f ref_check_dir/merged.bam \$BAM_FILES 2>/dev/null || cp \$(echo \$BAM_FILES | awk '{print \$1}') ref_check_dir/merged.bam
      samtools index ref_check_dir/merged.bam 2>/dev/null || true
      samtools depth -a ref_check_dir/merged.bam > ref_check_dir/reference_depth.txt 2>/dev/null || true
    fi

    # Calculate mapping rate and coverage (awk instead of bc for container compatibility)
    if [ "\$TOTAL_READS" -gt 0 ]; then
      MAPPING_RATE=\$(awk -v m="\$MAPPED_READS" -v t="\$TOTAL_READS" 'BEGIN {printf "%.4f", m / t * 100}')
    else
      MAPPING_RATE="0.0000"
    fi

    # Calculate reference coverage breadth
    REF_LENGTH=\$(grep -v "^>" ${reference} | tr -d '\\n' | wc -c)
    if [ -f ref_check_dir/reference_depth.txt ]; then
      COVERED_BASES=\$(awk '\$3 > 0 {count++} END {print count+0}' ref_check_dir/reference_depth.txt)
      if [ "\$REF_LENGTH" -gt 0 ]; then
        COVERAGE_BREADTH=\$(awk -v c="\$COVERED_BASES" -v r="\$REF_LENGTH" 'BEGIN {printf "%.4f", c / r * 100}')
      else
        COVERAGE_BREADTH="0.0000"
      fi
      MEAN_DEPTH=\$(awk '{sum+=\$3} END {if(NR>0) printf "%.2f", sum/NR; else print "0.00"}' ref_check_dir/reference_depth.txt)
    else
      COVERED_BASES=0
      COVERAGE_BREADTH="0.0000"
      MEAN_DEPTH="0.00"
    fi

    # Determine detection status
    if [ "\$(awk -v mr="\$MAPPING_RATE" -v cb="\$COVERAGE_BREADTH" 'BEGIN {print (mr > 1 && cb > 10) ? 1 : 0}')" -eq 1 ]; then
      DETECTION_STATUS="DETECTED"
    elif [ "\$(awk -v mr="\$MAPPING_RATE" -v cb="\$COVERAGE_BREADTH" 'BEGIN {print (mr > 0.1 || cb > 1) ? 1 : 0}')" -eq 1 ]; then
      DETECTION_STATUS="LOW_SIGNAL"
    else
      DETECTION_STATUS="NOT_DETECTED"
    fi

    # Write summary report
    cat > ref_check_dir/reference_report.txt << EOF
Reference Genome Detection Report
==================================
Sample: ${sample_id}
Reference: ${reference}

Detection Status: \$DETECTION_STATUS

Mapping Statistics:
  Total reads processed: \$TOTAL_READS
  Reads mapped to reference: \$MAPPED_READS
  Mapping rate: \${MAPPING_RATE}%

Coverage Statistics:
  Reference length: \$REF_LENGTH bp
  Bases covered (depth > 0): \$COVERED_BASES bp
  Coverage breadth: \${COVERAGE_BREADTH}%
  Mean depth: \${MEAN_DEPTH}x

Interpretation:
EOF

    if [ "\$DETECTION_STATUS" = "DETECTED" ]; then
      echo "  The reference genome was DETECTED in this sample." >> ref_check_dir/reference_report.txt
      echo "  Mapping rate >1% and coverage breadth >10% indicate presence." >> ref_check_dir/reference_report.txt
    elif [ "\$DETECTION_STATUS" = "LOW_SIGNAL" ]; then
      echo "  LOW SIGNAL: Some reads map to reference but coverage is limited." >> ref_check_dir/reference_report.txt
      echo "  This may indicate low viral load or partial genome presence." >> ref_check_dir/reference_report.txt
    else
      echo "  The reference genome was NOT DETECTED in this sample." >> ref_check_dir/reference_report.txt
      echo "  Very few or no reads mapped to the reference genome." >> ref_check_dir/reference_report.txt
    fi

    echo "" >> ref_check_dir/reference_report.txt

    # Generate coverage plot and per-segment stats
    if command -v python3 &>/dev/null && [ -f ref_check_dir/reference_depth.txt ] && [ -s ref_check_dir/reference_depth.txt ]; then
      python3 << 'PYEOF'
import sys
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd

    depth_file = "ref_check_dir/reference_depth.txt"
    df = pd.read_csv(depth_file, sep='\t', header=None, names=['contig', 'pos', 'depth'])

    if not df.empty:
        contigs = df['contig'].unique()
        n = len(contigs)

        if n > 1:
            fig, axes = plt.subplots(n, 1, figsize=(12, 2.5 * n), sharex=False)
            if n == 1:
                axes = [axes]
            colors = plt.cm.tab10.colors
            for i, (ctg, ax) in enumerate(zip(contigs, axes)):
                sub = df[df['contig'] == ctg]
                c = colors[i % len(colors)]
                ax.fill_between(sub['pos'], sub['depth'], alpha=0.7, color=c)
                ax.set_xlim(0, sub['pos'].max())
                ax.set_ylim(0, None)
                ax.set_ylabel('Depth')
                label = ctg if len(ctg) <= 60 else ctg[:57] + '...'
                ax.set_title(label, fontsize=9, loc='left')
            axes[-1].set_xlabel('Position (bp)')
            fig.suptitle('Reference Genome Coverage (per segment)', fontsize=11, y=1.0)
        else:
            fig, ax = plt.subplots(figsize=(12, 4))
            ax.fill_between(df['pos'], df['depth'], alpha=0.7, color='steelblue')
            ax.set_xlabel('Position (bp)')
            ax.set_ylabel('Coverage Depth')
            ax.set_title('Reference Genome Coverage')
            ax.set_xlim(0, df['pos'].max())
            ax.set_ylim(0, None)

        plt.tight_layout()
        plt.savefig('ref_check_dir/reference_coverage.png', dpi=150, bbox_inches='tight')
        plt.close()

        # Write per-segment stats
        with open('ref_check_dir/per_segment_stats.tsv', 'w') as f:
            f.write('segment\\tlength\\tcovered_bases\\tbreadth_%\\tmean_depth\\n')
            for ctg in contigs:
                sub = df[df['contig'] == ctg]
                seg_len = int(sub['pos'].max())
                covered = int((sub['depth'] > 0).sum())
                breadth = covered / seg_len * 100 if seg_len > 0 else 0
                mean_d = sub['depth'].mean()
                f.write(f'{ctg}\\t{seg_len}\\t{covered}\\t{breadth:.2f}\\t{mean_d:.2f}\\n')
        print(f"Coverage plot generated ({n} segment{'s' if n > 1 else ''})", file=sys.stderr)
except Exception as e:
    print(f"Could not generate coverage plot: {e}", file=sys.stderr)
PYEOF
    fi

    # Append per-segment stats to report if available
    if [ -f ref_check_dir/per_segment_stats.tsv ]; then
      SEG_COUNT=\$(tail -n +2 ref_check_dir/per_segment_stats.tsv | wc -l)
      if [ "\$SEG_COUNT" -gt 1 ]; then
        echo "Per-Segment Statistics:" >> ref_check_dir/reference_report.txt
        echo "  Segment                Length   Covered  Breadth%  MeanDepth" >> ref_check_dir/reference_report.txt
        tail -n +2 ref_check_dir/per_segment_stats.tsv | while IFS=\$'\\t' read -r seg slen cov br md; do
          printf "  %-22s %6s   %6s   %7s   %8s\\n" "\$seg" "\$slen" "\$cov" "\$br" "\$md" >> ref_check_dir/reference_report.txt
        done
        echo "" >> ref_check_dir/reference_report.txt
      fi
    fi
    echo "Generated by Virall Nextflow pipeline" >> ref_check_dir/reference_report.txt

    # Print summary to stdout for Nextflow log
    echo "=========================================="
    echo "REFERENCE CHECK: \$DETECTION_STATUS"
    echo "  Mapping rate: \${MAPPING_RATE}%"
    echo "  Coverage breadth: \${COVERAGE_BREADTH}%"
    echo "  Mean depth: \${MEAN_DEPTH}x"
    echo "=========================================="
    """
}

// ---------------------------------------------------------------------------
// PLOT – Virall-style plots (abundance, contig lengths, taxonomy, quality, coverage)
// ---------------------------------------------------------------------------
process PLOT {
    tag "${sample_id}"
    label "plot"
    publishDir "${params.outdir}/${sample_id}/07_plots", mode: "copy"

    input:
    tuple val(sample_id), path(quant_dir), path(quality_dir), path(viral_contigs), path(kaiju_dir), path(preprocess_dir)
    path(plot_script)

    output:
    tuple val(sample_id), path("plots_dir"), emit: plotted

    script:
    def viral_fasta = viral_contigs.name
    """
    mkdir -p plots_dir
    CONTIG_COUNT=\$(grep -c "^>" ${viral_fasta} 2>/dev/null || true)
    CONTIG_COUNT=\${CONTIG_COUNT:-0}
    if [ "\$CONTIG_COUNT" -eq 0 ]; then
      echo "PLOT: skipping – no viral contigs found for ${sample_id}"
      echo "No viral contigs >= ${params.min_contig_len} bp were found in this sample." > plots_dir/NO_VIRAL_SEQUENCES.txt
    else
      python ${plot_script} \\
        --quant-dir quant_dir \\
        --viral-contigs ${viral_fasta} \\
        --kaiju-dir kaiju_dir \\
        --checkv-dir ${quality_dir} \\
        --out-dir plots_dir
    fi
    """
}

// ---------------------------------------------------------------------------
// SINGLE-CELL: SC_MAP_VIRAL – Map barcoded reads to viral contigs
// ---------------------------------------------------------------------------
process SC_MAP_VIRAL {
    tag "${sample_id}"
    label "single_cell"
    publishDir "${params.outdir}/${sample_id}/08_single_cell/viral_mapping", mode: "copy"

    input:
    tuple val(sample_id), path(tagged_r2), path(viral_contigs)

    output:
    tuple val(sample_id), path("viral_aligned.bam"), path("viral_aligned.bam.bai"), emit: bam
    tuple val(sample_id), path("mapping_stats.txt"), emit: stats

    script:
    """
    # Index viral contigs
    bwa index ${viral_contigs}

    # Map barcoded reads (R2 = cDNA with barcode in read name) to viral contigs
    # Read names contain cell barcode and UMI from umi_tools extract
    bwa mem -t ${task.cpus} ${viral_contigs} ${tagged_r2} | \\
        samtools view -bS -F 4 - | \\
        samtools sort -@ ${task.cpus} -o viral_aligned.bam -

    samtools index viral_aligned.bam

    # Generate mapping stats
    samtools flagstat viral_aligned.bam > mapping_stats.txt
    echo "Unique cell barcodes with viral reads:" >> mapping_stats.txt
    # Barcode is the second-to-last underscore-separated field in read name
    samtools view viral_aligned.bam | awk -F'\\t' '{n=split(\$1,a,"_"); print a[n-1]}' | sort -u | wc -l >> mapping_stats.txt
    """
}

// ---------------------------------------------------------------------------
// SINGLE-CELL: SC_COUNT_CELLS – UMI-aware counting per cell
// ---------------------------------------------------------------------------
process SC_COUNT_CELLS {
    tag "${sample_id}"
    label "single_cell"
    publishDir "${params.outdir}/${sample_id}/08_single_cell", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("umi_counts.tsv"), emit: counts
    tuple val(sample_id), path("count_stats.log"), emit: stats

    script:
    """
    # umi_tools count requires CB, UB, and XT tags in BAM
    # umi_tools extract appends barcode and UMI to read name as:
    #   ORIGINALREADNAME_BARCODE_UMI
    # So barcode = second-to-last and UMI = last underscore-separated field
    # We also need XT tag = reference contig name (column 3) for --per-gene

    samtools view -h ${bam} | awk 'BEGIN{OFS="\\t"} 
        /^@/ {print; next}
        {
            n = split(\$1, parts, "_");
            cb = parts[n-1];
            ub = parts[n];
            print \$0, "CB:Z:"cb, "UB:Z:"ub, "XT:Z:"\$3
        }' | samtools view -bS - > tagged.bam
    
    samtools index tagged.bam

    # Count UMIs per cell per gene (viral contig)
    umi_tools count \
        --per-gene \
        --gene-tag=XT \
        --per-cell \
        --cell-tag=CB \
        --extract-umi-method=tag \
        --umi-tag=UB \
        -I tagged.bam \
        -S umi_counts.tsv \
        -L count_stats.log || {
            # Fallback: simple counting without umi_tools if it fails
            echo "umi_tools count failed, using simple counting" >> count_stats.log
            samtools view ${bam} | awk '{
                n = split(\$1, parts, "_");
                cb = parts[n-1];
                umi = parts[n];
                contig = \$3;
                key = cb"\\t"contig"\\t"umi;
                if (!(key in seen)) {
                    seen[key] = 1;
                    counts[cb"\\t"contig]++;
                }
            } END {
                print "cell\\tgene\\tcount";
                for (k in counts) print k"\\t"counts[k];
            }' > umi_counts.tsv
        }
    """
}

// ---------------------------------------------------------------------------
// SINGLE-CELL: SC_BUILD_MATRIX – Generate 10x-compatible MTX matrix
// ---------------------------------------------------------------------------
process SC_BUILD_MATRIX {
    tag "${sample_id}"
    label "single_cell"
    publishDir "${params.outdir}/${sample_id}/08_single_cell/matrix", mode: "copy"

    input:
    tuple val(sample_id), path(counts), path(viral_contigs), path(kaiju_dir)

    output:
    tuple val(sample_id), path("matrix.mtx"), path("barcodes.tsv"), path("features.tsv"), emit: matrix
    tuple val(sample_id), path("sc_viral_summary.tsv"), emit: summary

    script:
    """
    sc_build_matrix.py \
        --counts ${counts} \
        --contigs ${viral_contigs} \
        --kaiju-dir ${kaiju_dir} \
        --min-reads ${params.sc_min_reads_per_cell} \
        --min-viral-umis ${params.sc_min_viral_umis} \
        --output-dir .
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {
    // Default DB base dir: in-container path so -profile docker/singularity works without extra config
    def container_db = "/opt/virall/databases"
    def k_db = params.kaiju_db ?: (db_dir ? "${db_dir}/kaiju_db" : "${container_db}/kaiju_db")
    def c_db = params.checkv_db ?: (db_dir ? "${db_dir}/checkv_db" : "${container_db}/checkv_db")
    def v_db = params.vog_db ?: (db_dir ? "${db_dir}/vog_db" : "${container_db}/vog_db")
    def extract_script = file("${projectDir}/bin/extract_fasta_by_ids.py")
    def merge_script  = file("${projectDir}/bin/merge_quality.py")
    def plot_script   = file("${projectDir}/bin/run_plots.py")

    def g_db = params.genomad_db ?: (db_dir ? "${db_dir}/genomad_db" : "${container_db}/genomad_db")

    def ch_kaiju_db   = Channel.value(k_db)
    def ch_checkv_db  = Channel.value(c_db)
    def ch_genomad_db = Channel.value(g_db)
    def ch_vog_db     = v_db ? Channel.value(v_db) : Channel.value(null)

    // =========================================================================
    // SINGLE-CELL MODE: Extract barcodes, pool, run pipeline, then trace back
    // =========================================================================
    if (params.single_cell_mode) {
        // Extract single-cell reads from sample channel
        // ch_samples: (sample_id, read1, read2, single, long, sc_read1, sc_read2)
        ch_sc_input = ch_samples.map { t -> tuple(t[0], t[5], t[6]) }  // (sample_id, sc_read1, sc_read2)
        
        // Extract cell barcodes and UMIs
        SC_EXTRACT_BARCODES(ch_sc_input)
        
        // Pool reads for assembly (barcode info preserved in read names)
        SC_POOL_READS(SC_EXTRACT_BARCODES.out.tagged)
        
        // Create channel in PREPROCESS format: (sample_id, read1, read2, single, long, ...)
        // For 10x single-cell, R2 (cDNA) is used as single-end reads for assembly
        ch_sc_for_preprocess = SC_POOL_READS.out.pooled.map { t -> 
            tuple(
                t[0],              // sample_id
                file("${projectDir}/.placeholder_r1"),      // no paired R1
                file("${projectDir}/.placeholder_r2"),      // no paired R2
                t[1],              // pooled_single.fastq.gz as single-end
                file("${projectDir}/.placeholder_long"),    // no long
                file("${projectDir}/.placeholder_sc_r1"),   // placeholder
                file("${projectDir}/.placeholder_sc_r2")    // placeholder
            )
        }
        
        PREPROCESS(ch_sc_for_preprocess)
    } else {
        // Standard bulk mode
        PREPROCESS(ch_samples)
    }

    // Optional host filtering: when host_genome is set, filter out host reads
    if (params.host_genome) {
        def host_file = file(expand_path(params.host_genome))
        HOST_FILTER(
            PREPROCESS.out.preprocessed.map { t -> tuple(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10], host_file) }
        )
        ch_for_assemble = HOST_FILTER.out.host_filtered
    } else {
        ch_for_assemble = PREPROCESS.out.preprocessed
    }

    def ref_file = params.reference ? file(expand_path(params.reference)) : file("${projectDir}/.placeholder_ref")
    ASSEMBLE(ch_for_assemble, ref_file)

    // Optional reference check: when reference is set, map reads to reference and report detection
    // Uses host-filtered reads if host_genome was provided, otherwise uses preprocessed reads
    if (params.reference) {
        REFERENCE_CHECK(ch_for_assemble, ref_file)
    }

    KAIJU(ASSEMBLE.out.assembled, ch_kaiju_db)
    FILTER_VIRAL(KAIJU.out.kaiju_done, extract_script)
    RENAME_CONTIGS(FILTER_VIRAL.out.viral)

    // Run CheckV and geNomad in parallel on viral contigs
    VALIDATE(RENAME_CONTIGS.out.renamed, ch_checkv_db)
    GENOMAD(RENAME_CONTIGS.out.renamed, ch_genomad_db)

    // Merge quality assessments: CheckV for phages, geNomad for RNA/eukaryotic viruses
    // Join by sample_id: VALIDATE output is (sample_id, checkv_dir, viral_contigs, kaiju_dir, preprocess_dir)
    //                    GENOMAD output is (sample_id, genomad_out, viral_contigs, kaiju_dir, preprocess_dir)
    ch_checkv = VALIDATE.out.validated.map { t -> tuple(t[0], t[1]) }  // (sample_id, checkv_dir)
    ch_genomad = GENOMAD.out.genomad_done.map { t -> tuple(t[0], t[1]) }  // (sample_id, genomad_out)
    ch_viral_meta = RENAME_CONTIGS.out.renamed.map { t -> tuple(t[0], t[1], t[2], t[3]) }  // (sample_id, viral_contigs, kaiju_dir, preprocess_dir)

    ch_merge_input = ch_checkv
        .join(ch_genomad)
        .join(ch_viral_meta)
        .map { t -> tuple(t[0], t[1], t[2], t[3], t[4], t[5]) }
        // Result: (sample_id, checkv_dir, genomad_out, viral_contigs, kaiju_dir, preprocess_dir)

    MERGE_QUALITY(ch_merge_input, merge_script)

    if (v_db) {
        ANNOTATE(RENAME_CONTIGS.out.renamed, ch_vog_db)
        
        // Organize genes by taxonomy
        // Join ANNOTATE output with RENAME_CONTIGS output to get viral_contigs and kaiju_dir
        ch_organize_input = ANNOTATE.out.annotated
            .join(RENAME_CONTIGS.out.renamed)
            .map { sample_id, annotation_dir, viral_contigs, kaiju_dir, preprocess_dir ->
                tuple(sample_id, annotation_dir, viral_contigs, kaiju_dir)
            }
        ORGANIZE_GENES(ch_organize_input, ch_vog_db)
    }
    QUANTIFY(RENAME_CONTIGS.out.renamed)

    // Use merged quality for plotting
    ch_plot_input = QUANTIFY.out.quantified.join(MERGE_QUALITY.out.merged)
    PLOT(ch_plot_input, plot_script)

    // =========================================================================
    // SINGLE-CELL: Map barcoded reads back to viral contigs, count, build matrix
    // =========================================================================
    if (params.single_cell_mode) {
        // Join barcoded reads with viral contigs
        // For 10x, we only use R2 (cDNA) for mapping - R1 is just barcode/UMI
        ch_sc_map_input = SC_EXTRACT_BARCODES.out.tagged
            .join(RENAME_CONTIGS.out.renamed)
            .map { t -> tuple(t[0], t[2], t[3]) }
            // (sample_id, tagged_R2, viral_contigs)
        
        SC_MAP_VIRAL(ch_sc_map_input)
        
        SC_COUNT_CELLS(SC_MAP_VIRAL.out.bam)
        
        // Join counts with viral contigs and kaiju for matrix building
        ch_matrix_input = SC_COUNT_CELLS.out.counts
            .join(RENAME_CONTIGS.out.renamed)
            .map { t -> tuple(t[0], t[1], t[2], t[3]) }
            // (sample_id, umi_counts, viral_contigs, kaiju_dir)
        
        SC_BUILD_MATRIX(ch_matrix_input)
    }
}
