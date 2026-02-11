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
    """
    mkdir -p preprocess_dir
    if [ -s "${read1}" ] && [ -s "${read2}" ]; then
      fastp -i ${read1} -I ${read2} -o preprocess_dir/trimmed_R1.fastq.gz -O preprocess_dir/trimmed_R2.fastq.gz \\
        --html preprocess_dir/fastp_pe.html --json preprocess_dir/fastp_pe.json \\
        --thread ${params.threads} --qualified_quality_phred ${q} --length_required ${ml} \\
        --detect_adapter_for_pe --cut_front --cut_tail --cut_mean_quality ${q} --trim_poly_g
    fi
    if [ -s "${single}" ] && [ "${single.name}" != ".placeholder_single" ]; then
      fastp -i ${single} -o preprocess_dir/trimmed_single.fastq.gz \\
        --html preprocess_dir/fastp_single.html --json preprocess_dir/fastp_single.json \\
        --thread ${params.threads} --qualified_quality_phred ${q} --length_required ${ml} --trim_poly_g
    fi
    if [ -s "${long_reads}" ] && [ "${long_reads.name}" != ".placeholder_long" ]; then
      fastplong -i ${long_reads} -o preprocess_dir/trimmed_long.fastq.gz \\
        --html preprocess_dir/fastplong.html --json preprocess_dir/fastplong.json \\
        --thread ${params.threads} --qualified_quality_phred ${params.quality_phred_long} \\
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
      minimap2 -ax sr -t ${params.threads} ${host_ref} ${trimmed_r1} ${trimmed_r2} 2>host_filter_dir/host_filter_pe.log | \\
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
      minimap2 -ax sr -t ${params.threads} ${host_ref} ${trimmed_single} 2>host_filter_dir/host_filter_single.log | \\
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
      minimap2 -ax \$MINIMAP2_LONG_PRESET -t ${params.threads} ${host_ref} ${trimmed_long} 2>host_filter_dir/host_filter_long.log | \\
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

    output:
    tuple val(sample_id), path("contigs"), path("scaffolds"), path("preprocess_dir"), emit: assembled

    script:
    def mem = task.memory ? task.memory.toGiga().toString() : params.memory.replaceAll(/[Gg]/, '')
    def rna = params.rna_mode ? '--rnaviral' : ''
    def ion = params.iontorrent ? '--iontorrent' : ''
    def metaviral = params.metaviral_mode ? '--metaviral' : ''
    // Note: --metaviral and --trusted-contigs may conflict; when metaviral_mode is on, skip trusted-contigs
    def ref = (params.reference && !params.metaviral_mode) ? "--trusted-contigs ${file(params.reference)}" : ''
    """
    mkdir -p assembly_dir preprocess_dir
    
    # Warn if metaviral mode is enabled but reference is also provided
    if [ -n "${metaviral}" ] && [ -n "${params.reference ?: ''}" ]; then
      echo "WARNING: metaviral_mode is enabled, but a reference genome is also provided."
      echo "         --trusted-contigs is disabled when using --metaviral to avoid potential conflicts."
      echo "         The reference genome will still be used for REFERENCE_CHECK if enabled."
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
      if [ "\$STRATEGY" = "short_only" ] || [ "\$STRATEGY" = "hybrid" ]; then
        if [ "\$HAS_SHORT" = "1" ] || [ "\$HAS_SINGLE" = "1" ]; then
          SPADES_OPTS="-o assembly_dir/spades -t ${params.threads} -m ${mem} --only-assembler ${rna} ${ion} ${metaviral} ${ref}"
          [ -s preprocess_dir/trimmed_R1.fastq.gz ] && SPADES_OPTS="\$SPADES_OPTS -1 preprocess_dir/trimmed_R1.fastq.gz -2 preprocess_dir/trimmed_R2.fastq.gz"
          [ -s preprocess_dir/trimmed_single.fastq.gz ] && SPADES_OPTS="\$SPADES_OPTS -s preprocess_dir/trimmed_single.fastq.gz"
          # Add long reads for hybrid mode (--nanopore or --pacbio)
          [ "\$STRATEGY" = "hybrid" ] && [ "\$HAS_LONG" = "1" ] && SPADES_OPTS="\$SPADES_OPTS \$SPADES_LONG_FLAG preprocess_dir/trimmed_long.fastq.gz"
          spades.py \$SPADES_OPTS
          cp assembly_dir/spades/contigs.fasta contigs 2>/dev/null || cp assembly_dir/spades/transcripts.fasta contigs 2>/dev/null || touch contigs
          cp assembly_dir/spades/scaffolds.fasta scaffolds 2>/dev/null || cp assembly_dir/spades/hard_filtered_transcripts.fasta scaffolds 2>/dev/null || cp contigs scaffolds
        else
          echo "Warning: strategy=\$STRATEGY but no short/single reads available"
        fi
      fi

      if [ "\$STRATEGY" = "long_only" ]; then
        if [ "\$HAS_LONG" = "1" ]; then
          flye \$FLYE_INPUT_FLAG preprocess_dir/trimmed_long.fastq.gz --out-dir assembly_dir/flye -t ${params.threads}
          cp assembly_dir/flye/assembly.fasta contigs 2>/dev/null || true
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
  CONTIG_COUNT=\$(grep -c "^>" contigs 2>/dev/null || echo 0)
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
    kaiju -t "\$NODES" -f "\$FMI" -i contigs -o kaiju_dir/kaiju_results.tsv -z ${params.threads}
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
    CONTIG_COUNT=\$(grep -c "^>" viral_contigs.fasta 2>/dev/null || echo 0)
    if [ "\$CONTIG_COUNT" -eq 0 ]; then
      echo "CheckV: skipping – no viral contigs found for ${sample_id}"
      printf "contig_id\\tcontig_length\\tprovirus\\tproviral_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tmiuvig_quality\\tcompleteness\\tcompleteness_method\\tcontamination\\tkmer_freq\\twarnings\\n" > checkv_dir/quality_summary.tsv
    else
      TOTAL_BP=\$(grep -v "^>" viral_contigs.fasta | tr -d '\\n' | wc -c)
      if [ "\$TOTAL_BP" -lt 200 ]; then
        echo "CheckV: skipping – viral contigs too short (\${TOTAL_BP} bp) for meaningful quality assessment"
        printf "contig_id\\tcontig_length\\tprovirus\\tproviral_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tmiuvig_quality\\tcompleteness\\tcompleteness_method\\tcontamination\\tkmer_freq\\twarnings\\n" > checkv_dir/quality_summary.tsv
      else
        checkv end_to_end viral_contigs.fasta checkv_dir -d "\$CHECKV_D" -t 1
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
    CONTIG_COUNT=\$(grep -c "^>" ${viral_contigs} 2>/dev/null || echo 0)
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
        --splits ${params.threads} \\
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

    output:
    tuple val(sample_id), path("merged_quality"), path(viral_contigs), path(kaiju_dir), path(preprocess_dir), emit: merged

    script:
    """
    mkdir -p merged_quality

    python3 << 'PYMERGE'
import pandas as pd
import sys
from pathlib import Path

# Input files
checkv_file = Path("${checkv_dir}/quality_summary.tsv")
genomad_file = Path("${genomad_dir}/virus_summary.tsv")
kaiju_file = Path("${kaiju_dir}/kaiju_results_with_names.tsv")

# Output file
output_file = Path("merged_quality/quality_summary.tsv")

# Read CheckV results
checkv_df = pd.DataFrame()
if checkv_file.exists() and checkv_file.stat().st_size > 0:
    try:
        checkv_df = pd.read_csv(checkv_file, sep='\\t')
        checkv_df['quality_source'] = 'checkv'
        print(f"Loaded {len(checkv_df)} CheckV results", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not read CheckV file: {e}", file=sys.stderr)

# Read geNomad results
genomad_df = pd.DataFrame()
if genomad_file.exists() and genomad_file.stat().st_size > 0:
    try:
        genomad_df = pd.read_csv(genomad_file, sep='\\t')
        print(f"Loaded {len(genomad_df)} geNomad results", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not read geNomad file: {e}", file=sys.stderr)

# Read Kaiju taxonomy to determine phage vs other viruses
kaiju_taxonomy = {}
if kaiju_file.exists():
    try:
        with open(kaiju_file) as f:
            for line in f:
                if line.startswith('C'):
                    parts = line.strip().split('\\t')
                    if len(parts) >= 4:
                        contig_id = parts[1].strip()
                        lineage = parts[3] if len(parts) == 4 else ';'.join(parts[3:10])
                        kaiju_taxonomy[contig_id] = lineage.lower()
        print(f"Loaded taxonomy for {len(kaiju_taxonomy)} contigs", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not read Kaiju file: {e}", file=sys.stderr)

def is_phage(contig_id, lineage_str):
    \"\"\"Determine if a contig is likely a bacteriophage based on Kaiju taxonomy.\"\"\"
    if not lineage_str:
        return False
    lineage_lower = lineage_str.lower()
    
    # First, check for known NON-phage viral groups that might contain phage-like substrings
    # Kitrinoviricota contains 'inovir' but is actually an RNA virus realm (Flaviviridae, etc.)
    non_phage_indicators = [
        'kitrinoviricota',  # RNA viruses (Flaviviridae, etc.)
        'pisuviricota',     # RNA viruses (Picornavirales, etc.)
        'negarnaviricota',  # RNA viruses (negative-sense)
        'nucleocytoviricota',  # Giant DNA viruses (not phages)
        'artverviricota',   # Retroviruses
        'flaviviridae', 'togaviridae', 'coronaviridae', 'picornaviridae',
        'rhabdoviridae', 'paramyxoviridae', 'orthomyxoviridae', 'bunyavirales',
        'herpesviridae', 'poxviridae', 'adenoviridae', 'papillomaviridae',
        'polyomaviridae', 'retroviridae', 'hepadnaviridae', 'parvoviridae',
        'mimiviridae', 'phycodnaviridae', 'iridoviridae', 'ascoviridae',
    ]
    
    if any(indicator in lineage_lower for indicator in non_phage_indicators):
        return False
    
    # Now check for phage-specific indicators
    # Use more specific patterns to avoid false positives
    phage_indicators = [
        'caudovir',        # Tailed phages (Caudovirales, Caudoviricetes)
        'phage',           # Explicit "phage" in name
        'bacteriophage',   # Explicit bacteriophage
        'siphovir',        # Siphoviridae
        'myovir',          # Myoviridae
        'podovir',         # Podoviridae
        'autographivir',   # Autographiviridae
        'demerecvir',      # Demerecviridae
        'herellevir',      # Herelleviridae
        'inoviridae',      # Inoviridae (filamentous phages) - use full family name
        'microviridae',    # Microviridae - use full family name
        'leviviricetes',   # RNA phages
        'cystoviridae',    # dsRNA phages
        'tectiviridae',    # Tectiviridae
        'corticoviridae',  # Corticoviridae
        'plasmaviridae',   # Plasmaviridae
        'sphaerolipoviridae',  # Sphaerolipoviridae
        'uroviricota',     # Phage realm
        'peduoviridae',    # Peduoviridae
        'drexlerviridae',  # Drexlerviridae
        'ackermannviridae',# Ackermannviridae
    ]
    return any(indicator in lineage_lower for indicator in phage_indicators)

def genomad_to_quality_tier(row):
    \"\"\"Convert geNomad metrics to CheckV-style quality tier.\"\"\"
    # geNomad columns: seq_name, length, topology, coordinates, n_genes, genetic_code,
    #                  virus_score, fdr, n_hallmarks, marker_enrichment, taxonomy
    topology = str(row.get('topology', '')).lower()
    n_hallmarks = int(row.get('n_hallmarks', 0)) if pd.notna(row.get('n_hallmarks')) else 0
    virus_score = float(row.get('virus_score', 0)) if pd.notna(row.get('virus_score')) else 0

    # Circular genomes with hallmarks are likely complete
    if topology == 'circular' and n_hallmarks >= 1:
        return 'Complete', 100.0
    elif topology == 'circular':
        return 'High-quality', 90.0
    elif n_hallmarks >= 3:
        return 'High-quality', 80.0
    elif n_hallmarks >= 1:
        return 'Medium-quality', 50.0
    elif virus_score >= 0.9:
        return 'Low-quality', 30.0
    elif virus_score >= 0.7:
        return 'Low-quality', 20.0
    else:
        return 'Not-determined', None

# Create standardized output
merged_rows = []

# Get all contig IDs from CheckV (primary source)
if not checkv_df.empty and 'contig_id' in checkv_df.columns:
    for _, row in checkv_df.iterrows():
        contig_id = row['contig_id']
        lineage = kaiju_taxonomy.get(contig_id, '')

        # Determine which quality source to use
        use_checkv = is_phage(contig_id, lineage)

        if use_checkv:
            # Use CheckV results for phages
            merged_rows.append({
                'contig_id': contig_id,
                'contig_length': row.get('contig_length', 0),
                'quality_source': 'checkv',
                'checkv_quality': row.get('checkv_quality', 'Not-determined'),
                'completeness': row.get('completeness', None),
                'completeness_method': row.get('completeness_method', ''),
                'contamination': row.get('contamination', None),
                'viral_genes': row.get('viral_genes', 0),
                'host_genes': row.get('host_genes', 0),
                'taxonomy_hint': 'phage' if is_phage(contig_id, lineage) else 'other',
                'provirus': row.get('provirus', 'No'),
            })
        else:
            # For non-phages: combine best data from CheckV and geNomad
            genomad_row = None
            if not genomad_df.empty and 'seq_name' in genomad_df.columns:
                match = genomad_df[genomad_df['seq_name'] == contig_id]
                if not match.empty:
                    genomad_row = match.iloc[0]

            # Get CheckV completeness if available
            checkv_completeness = row.get('completeness', None)
            checkv_quality = row.get('checkv_quality', 'Not-determined')
            checkv_method = row.get('completeness_method', '')
            
            # Determine best completeness source:
            # - If CheckV has AAI-based completeness, use it (more reliable for viruses in DB)
            # - Otherwise, use geNomad's estimate
            use_checkv_completeness = pd.notna(checkv_completeness) and checkv_completeness != ''
            
            if genomad_row is not None:
                genomad_quality, genomad_completeness = genomad_to_quality_tier(genomad_row)
                
                # Choose best completeness
                if use_checkv_completeness:
                    final_completeness = checkv_completeness
                    final_quality = checkv_quality
                    final_method = checkv_method
                    source = 'checkv+genomad'  # Combined: CheckV quality, geNomad identification
                else:
                    final_completeness = genomad_completeness
                    final_quality = genomad_quality
                    final_method = 'genomad_hallmarks'
                    source = 'genomad'
                
                merged_rows.append({
                    'contig_id': contig_id,
                    'contig_length': row.get('contig_length', 0),
                    'quality_source': source,
                    'checkv_quality': final_quality,
                    'completeness': final_completeness,
                    'completeness_method': final_method,
                    'contamination': row.get('contamination', None),
                    'viral_genes': max(row.get('viral_genes', 0), int(genomad_row.get('n_genes', 0) or 0)),
                    'host_genes': row.get('host_genes', 0),
                    'taxonomy_hint': 'rna_or_eukaryotic',
                    'provirus': row.get('provirus', 'No'),
                    'genomad_score': genomad_row.get('virus_score', None),
                    'genomad_topology': genomad_row.get('topology', ''),
                    'genomad_hallmarks': genomad_row.get('n_hallmarks', 0),
                    'genomad_taxonomy': genomad_row.get('taxonomy', ''),
                })
            else:
                # Fallback to CheckV if geNomad didn't process this contig
                merged_rows.append({
                    'contig_id': contig_id,
                    'contig_length': row.get('contig_length', 0),
                    'quality_source': 'checkv_fallback',
                    'checkv_quality': row.get('checkv_quality', 'Not-determined'),
                    'completeness': row.get('completeness', None),
                    'completeness_method': row.get('completeness_method', ''),
                    'contamination': row.get('contamination', None),
                    'viral_genes': row.get('viral_genes', 0),
                    'host_genes': row.get('host_genes', 0),
                    'taxonomy_hint': 'other',
                    'provirus': row.get('provirus', 'No'),
                })

# If CheckV failed but geNomad worked, use geNomad only
if not merged_rows and not genomad_df.empty:
    print("Using geNomad results only (CheckV results empty)", file=sys.stderr)
    for _, row in genomad_df.iterrows():
        contig_id = row.get('seq_name', '')
        quality_tier, completeness = genomad_to_quality_tier(row)
        merged_rows.append({
            'contig_id': contig_id,
            'contig_length': row.get('length', 0),
            'quality_source': 'genomad',
            'checkv_quality': quality_tier,
            'completeness': completeness,
            'completeness_method': 'genomad_hallmarks',
            'contamination': None,
            'viral_genes': row.get('n_genes', 0),
            'host_genes': 0,
            'taxonomy_hint': 'unknown',
            'provirus': 'No',
            'genomad_score': row.get('virus_score', None),
            'genomad_topology': row.get('topology', ''),
            'genomad_hallmarks': row.get('n_hallmarks', 0),
            'genomad_taxonomy': row.get('taxonomy', ''),
        })

# Create output DataFrame
if merged_rows:
    merged_df = pd.DataFrame(merged_rows)
    merged_df.to_csv(output_file, sep='\\t', index=False)
    print(f"Wrote {len(merged_df)} rows to merged quality summary", file=sys.stderr)

    # Summary stats
    source_counts = merged_df['quality_source'].value_counts()
    quality_counts = merged_df['checkv_quality'].value_counts()
    print(f"Quality sources: {source_counts.to_dict()}", file=sys.stderr)
    print(f"Quality tiers: {quality_counts.to_dict()}", file=sys.stderr)
else:
    # Create empty file with header
    pd.DataFrame(columns=['contig_id', 'contig_length', 'quality_source', 'checkv_quality',
                          'completeness', 'completeness_method', 'contamination']).to_csv(output_file, sep='\\t', index=False)
    print("No quality data to merge", file=sys.stderr)

# Also copy original files to merged_quality for reference
import shutil
if checkv_file.exists():
    shutil.copy(checkv_file, "merged_quality/checkv_quality_summary.tsv")
if genomad_file.exists() and genomad_file.stat().st_size > 0:
    shutil.copy(genomad_file, "merged_quality/genomad_virus_summary.tsv")

PYMERGE
    """
}

// ---------------------------------------------------------------------------
// ANNOTATE – Prodigal + hmmscan (VOG)
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
    CONTIG_COUNT=\$(grep -c "^>" viral_contigs.fasta 2>/dev/null || echo 0)
    if [ "\$CONTIG_COUNT" -eq 0 ]; then
      echo "ANNOTATE: skipping – no viral contigs found for ${sample_id}"
      touch annotation_dir/proteins.faa
    else
      prodigal -i viral_contigs.fasta -a annotation_dir/proteins.faa -p meta -q
      VOG_HMM=\$(find ${vog_db} -maxdepth 1 \\( -name 'vog_all.hmm' -o -name 'vog.hmm' \\) 2>/dev/null | head -1)
      if [ -n "\$VOG_HMM" ] && [ -f "\$VOG_HMM" ]; then
        hmmscan --cpu ${params.threads} -o annotation_dir/vog_out.txt --domtblout annotation_dir/vog_domains.txt "\$VOG_HMM" annotation_dir/proteins.faa
      fi
    fi
    """
}

// ---------------------------------------------------------------------------
// ORGANIZE_GENES – Organize genes by Kaiju taxonomy
// ---------------------------------------------------------------------------
process ORGANIZE_GENES {
    tag "${sample_id}"
    label "annotate"
    publishDir "${params.outdir}/${sample_id}/05_gene_predictions", mode: "copy"

    input:
    tuple val(sample_id), path(annotation_dir), path(viral_contigs), path(kaiju_dir)

    output:
    tuple val(sample_id), path("by_taxonomy"), emit: organized_genes

    script:
    """
    mkdir -p by_taxonomy
    if [ -s ${annotation_dir}/proteins.faa ] && [ -s ${viral_contigs} ]; then
      organize_genes_by_taxonomy.py \
          --proteins ${annotation_dir}/proteins.faa \
          --contigs ${viral_contigs} \
          --kaiju ${kaiju_dir}/kaiju_results_with_names.tsv \
          --outdir by_taxonomy
    else
      echo "ORGANIZE_GENES: skipping – no viral contigs or proteins found for ${sample_id}"
      touch by_taxonomy/gene_taxonomy_mapping.tsv
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
    CONTIG_COUNT=\$(grep -c "^>" viral_contigs.fasta 2>/dev/null || echo 0)
    if [ "\$CONTIG_COUNT" -eq 0 ]; then
      echo "QUANTIFY: skipping – no viral contigs found for ${sample_id}"
      touch quant_dir/mapped.bam quant_dir/mapped.bam.bai quant_dir/depth.txt
    else
      bwa index -p quant_dir/idx viral_contigs.fasta
      # Use -s (non-empty file) instead of -f to avoid matching 0-byte placeholders
      # In single-cell mode, R1/R2 are empty placeholders; actual data is in trimmed_single
      if [ -s preprocess_dir/trimmed_R1.fastq.gz ] && [ -s preprocess_dir/trimmed_R2.fastq.gz ]; then
        bwa mem -t ${params.threads} quant_dir/idx preprocess_dir/trimmed_R1.fastq.gz preprocess_dir/trimmed_R2.fastq.gz | samtools sort -o quant_dir/mapped.bam -
      elif [ -s preprocess_dir/trimmed_single.fastq.gz ]; then
        bwa mem -t ${params.threads} quant_dir/idx preprocess_dir/trimmed_single.fastq.gz | samtools sort -o quant_dir/mapped.bam -
      elif [ -s preprocess_dir/trimmed_long.fastq.gz ]; then
        MINIMAP2_LONG_PRESET=\$( [ "${params.long_read_tech}" = "pacbio" ] && echo "map-pb" || echo "map-ont" )
        minimap2 -t ${params.threads} -ax \$MINIMAP2_LONG_PRESET viral_contigs.fasta preprocess_dir/trimmed_long.fastq.gz | samtools sort -o quant_dir/mapped.bam -
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
      bwa mem -t ${params.threads} ref_check_dir/ref_idx ${trimmed_r1} ${trimmed_r2} 2>/dev/null | \\
        samtools sort -o ref_check_dir/mapped_pe.bam - 2>/dev/null
      samtools index ref_check_dir/mapped_pe.bam 2>/dev/null || true
      PE_MAPPED=\$(samtools view -c -F 4 ref_check_dir/mapped_pe.bam 2>/dev/null || echo 0)
      PE_TOTAL=\$(samtools view -c ref_check_dir/mapped_pe.bam 2>/dev/null || echo 0)
      MAPPED_READS=\$((MAPPED_READS + PE_MAPPED))
      TOTAL_READS=\$((TOTAL_READS + PE_TOTAL))
    fi

    if [ -s "${trimmed_single}" ] && [ "${trimmed_single.name}" != ".placeholder_single" ]; then
      bwa mem -t ${params.threads} ref_check_dir/ref_idx ${trimmed_single} 2>/dev/null | \\
        samtools sort -o ref_check_dir/mapped_single.bam - 2>/dev/null
      samtools index ref_check_dir/mapped_single.bam 2>/dev/null || true
      SINGLE_MAPPED=\$(samtools view -c -F 4 ref_check_dir/mapped_single.bam 2>/dev/null || echo 0)
      SINGLE_TOTAL=\$(samtools view -c ref_check_dir/mapped_single.bam 2>/dev/null || echo 0)
      MAPPED_READS=\$((MAPPED_READS + SINGLE_MAPPED))
      TOTAL_READS=\$((TOTAL_READS + SINGLE_TOTAL))
    fi

    if [ -s "${trimmed_long}" ] && [ "${trimmed_long.name}" != ".placeholder_long" ]; then
      minimap2 -t ${params.threads} -ax ${minimap2_preset} ${reference} ${trimmed_long} 2>/dev/null | \\
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

    # Calculate mapping rate and coverage
    if [ "\$TOTAL_READS" -gt 0 ]; then
      MAPPING_RATE=\$(echo "scale=4; \$MAPPED_READS / \$TOTAL_READS * 100" | bc)
    else
      MAPPING_RATE="0.0000"
    fi

    # Calculate reference coverage breadth
    REF_LENGTH=\$(grep -v "^>" ${reference} | tr -d '\\n' | wc -c)
    if [ -f ref_check_dir/reference_depth.txt ]; then
      COVERED_BASES=\$(awk '\$3 > 0 {count++} END {print count+0}' ref_check_dir/reference_depth.txt)
      if [ "\$REF_LENGTH" -gt 0 ]; then
        COVERAGE_BREADTH=\$(echo "scale=4; \$COVERED_BASES / \$REF_LENGTH * 100" | bc)
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
    if [ "\$(echo "\$MAPPING_RATE > 1" | bc)" -eq 1 ] && [ "\$(echo "\$COVERAGE_BREADTH > 10" | bc)" -eq 1 ]; then
      DETECTION_STATUS="DETECTED"
    elif [ "\$(echo "\$MAPPING_RATE > 0.1" | bc)" -eq 1 ] || [ "\$(echo "\$COVERAGE_BREADTH > 1" | bc)" -eq 1 ]; then
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
      echo "  Mapping rate >\$1% and coverage breadth >10% indicate presence." >> ref_check_dir/reference_report.txt
    elif [ "\$DETECTION_STATUS" = "LOW_SIGNAL" ]; then
      echo "  LOW SIGNAL: Some reads map to reference but coverage is limited." >> ref_check_dir/reference_report.txt
      echo "  This may indicate low viral load or partial genome presence." >> ref_check_dir/reference_report.txt
    else
      echo "  The reference genome was NOT DETECTED in this sample." >> ref_check_dir/reference_report.txt
      echo "  Very few or no reads mapped to the reference genome." >> ref_check_dir/reference_report.txt
    fi

    echo "" >> ref_check_dir/reference_report.txt
    echo "Generated by Virall Nextflow pipeline" >> ref_check_dir/reference_report.txt

    # Generate coverage plot (simple text-based if no python, or use python if available)
    if command -v python3 &>/dev/null && [ -f ref_check_dir/reference_depth.txt ] && [ -s ref_check_dir/reference_depth.txt ]; then
      python3 << 'PYEOF'
import sys
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd

    depth_file = "ref_check_dir/reference_depth.txt"
    df = pd.read_csv(depth_file, sep='\\t', header=None, names=['contig', 'pos', 'depth'])

    if not df.empty:
        fig, ax = plt.subplots(figsize=(12, 4))
        ax.fill_between(df['pos'], df['depth'], alpha=0.7, color='steelblue')
        ax.set_xlabel('Position (bp)')
        ax.set_ylabel('Coverage Depth')
        ax.set_title('Reference Genome Coverage')
        ax.set_xlim(0, df['pos'].max())
        ax.set_ylim(0, None)
        plt.tight_layout()
        plt.savefig('ref_check_dir/reference_coverage.png', dpi=150)
        plt.close()
        print("Coverage plot generated", file=sys.stderr)
except Exception as e:
    print(f"Could not generate coverage plot: {e}", file=sys.stderr)
PYEOF
    fi

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
    CONTIG_COUNT=\$(grep -c "^>" ${viral_fasta} 2>/dev/null || echo 0)
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
    bwa mem -t ${params.threads} ${viral_contigs} ${tagged_r2} | \\
        samtools view -bS -F 4 - | \\
        samtools sort -@ ${params.threads} -o viral_aligned.bam -

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

    ASSEMBLE(ch_for_assemble)

    // Optional reference check: when reference is set, map reads to reference and report detection
    // Uses host-filtered reads if host_genome was provided, otherwise uses preprocessed reads
    if (params.reference) {
        def ref_file = file(expand_path(params.reference))
        REFERENCE_CHECK(ch_for_assemble, ref_file)
    }

    KAIJU(ASSEMBLE.out.assembled, ch_kaiju_db)
    FILTER_VIRAL(KAIJU.out.kaiju_done, extract_script)

    // Run CheckV and geNomad in parallel on viral contigs
    VALIDATE(FILTER_VIRAL.out.viral, ch_checkv_db)
    GENOMAD(FILTER_VIRAL.out.viral, ch_genomad_db)

    // Merge quality assessments: CheckV for phages, geNomad for RNA/eukaryotic viruses
    // Join by sample_id: VALIDATE output is (sample_id, checkv_dir, viral_contigs, kaiju_dir, preprocess_dir)
    //                    GENOMAD output is (sample_id, genomad_out, viral_contigs, kaiju_dir, preprocess_dir)
    ch_checkv = VALIDATE.out.validated.map { t -> tuple(t[0], t[1]) }  // (sample_id, checkv_dir)
    ch_genomad = GENOMAD.out.genomad_done.map { t -> tuple(t[0], t[1]) }  // (sample_id, genomad_out)
    ch_viral_meta = FILTER_VIRAL.out.viral.map { t -> tuple(t[0], t[1], t[2], t[3]) }  // (sample_id, viral_contigs, kaiju_dir, preprocess_dir)

    ch_merge_input = ch_checkv
        .join(ch_genomad)
        .join(ch_viral_meta)
        .map { t -> tuple(t[0], t[1], t[2], t[3], t[4], t[5]) }
        // Result: (sample_id, checkv_dir, genomad_out, viral_contigs, kaiju_dir, preprocess_dir)

    MERGE_QUALITY(ch_merge_input)

    if (v_db) {
        ANNOTATE(FILTER_VIRAL.out.viral, ch_vog_db)
        
        // Organize genes by taxonomy
        // Join ANNOTATE output with FILTER_VIRAL output to get viral_contigs and kaiju_dir
        ch_organize_input = ANNOTATE.out.annotated
            .join(FILTER_VIRAL.out.viral)
            .map { sample_id, annotation_dir, viral_contigs, kaiju_dir, preprocess_dir ->
                tuple(sample_id, annotation_dir, viral_contigs, kaiju_dir)
            }
        ORGANIZE_GENES(ch_organize_input)
    }
    QUANTIFY(FILTER_VIRAL.out.viral)

    // Use merged quality for plotting
    ch_plot_input = QUANTIFY.out.quantified.join(MERGE_QUALITY.out.merged)
    PLOT(ch_plot_input, plot_script)

    // =========================================================================
    // SINGLE-CELL: Map barcoded reads back to viral contigs, count, build matrix
    // =========================================================================
    if (params.single_cell_mode) {
        // Join barcoded reads with viral contigs
        // SC_EXTRACT_BARCODES.out.tagged: (sample_id, tagged_R1, tagged_R2)
        // FILTER_VIRAL.out.viral: (sample_id, viral_contigs, kaiju_dir, preprocess_dir)
        // SC_EXTRACT_BARCODES.out.tagged: (sample_id, tagged_R1, tagged_R2)
        // FILTER_VIRAL.out.viral: (sample_id, viral_contigs, kaiju_dir, preprocess_dir)
        // For 10x, we only use R2 (cDNA) for mapping - R1 is just barcode/UMI
        ch_sc_map_input = SC_EXTRACT_BARCODES.out.tagged
            .join(FILTER_VIRAL.out.viral)
            .map { t -> tuple(t[0], t[2], t[3]) }
            // (sample_id, tagged_R2, viral_contigs)
        
        SC_MAP_VIRAL(ch_sc_map_input)
        
        SC_COUNT_CELLS(SC_MAP_VIRAL.out.bam)
        
        // Join counts with viral contigs and kaiju for matrix building
        ch_matrix_input = SC_COUNT_CELLS.out.counts
            .join(FILTER_VIRAL.out.viral)
            .map { t -> tuple(t[0], t[1], t[2], t[3]) }
            // (sample_id, umi_counts, viral_contigs, kaiju_dir)
        
        SC_BUILD_MATRIX(ch_matrix_input)
    }
}
