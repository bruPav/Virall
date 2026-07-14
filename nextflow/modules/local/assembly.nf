/*
 * Virall Nextflow pipeline — ASSEMBLY module
 * ASSEMBLE_SHORT (SPAdes), ASSEMBLE_LONG (Flye), ASSEMBLE_HYBRID, REF_ASSEMBLE
 */

// ---------------------------------------------------------------------------
// ASSEMBLE_SHORT – SPAdes for short/single-end reads only
// ---------------------------------------------------------------------------
process ASSEMBLE_SHORT {
    tag { sample_id }
    label "assemble"
    publishDir { "${params.outdir}/${sample_id}/01_assembly" }, mode: "copy", pattern: "contigs*"
    publishDir { "${params.outdir}/${sample_id}/01_assembly" }, mode: "copy", pattern: "scaffolds"

    input:
    tuple val(sample_id), path(r1), path(r2), path(single), val(short_tech)
    path(reference_file)
    path(assembly_utils_script)

    output:
    tuple val(sample_id), path("contigs"), path("scaffolds"), emit: assembled

    script:
    def mem = task.memory.toGiga().toString()
    def ion = (short_tech == "iontorrent") ? '--iontorrent' : ''
    def has_ref = reference_file.name != '.placeholder_ref'
    def ref_fasta = has_ref ? "reference.fasta" : ''
    def ref = has_ref ? "--trusted-contigs ${ref_fasta}" : ''
    def rna = (params.rna_mode && !has_ref) ? '--rnaviral' : ''
    def metaviral = (params.metaviral_mode && !has_ref) ? '--metaviral' : ''
    """
    source ${assembly_utils_script}

    if [ -n "${ref_fasta}" ] && [ -f "${reference_file}" ]; then
      ln -sf "${reference_file}" "${ref_fasta}"
    fi
    if [ -n "${ref}" ]; then
      echo "NOTE: Reference genome provided — using --trusted-contigs for reference-guided assembly."
      if [ "${params.rna_mode}" = "true" ] || [ "${params.metaviral_mode}" = "true" ]; then
        echo "      Skipping --rnaviral/--metaviral (incompatible with --trusted-contigs)."
      fi
    fi

    HAS_SHORT=\$( [ -s "${r1}" ] && [ -s "${r2}" ] && echo 1 || echo 0 )
    HAS_SINGLE=\$( [ -s "${single}" ] && [ "${single.name}" != ".placeholder_single" ] && echo 1 || echo 0 )

    if [ "\$HAS_SHORT" = "0" ] && [ "\$HAS_SINGLE" = "0" ]; then
      echo "ASSEMBLE_SHORT: no short/single reads for ${sample_id}"
      touch contigs scaffolds
    else
      SPADES_OPTS="-o assembly_dir/spades -t ${task.cpus} -m ${mem} --only-assembler ${rna} ${ion} ${metaviral} ${ref}"
      [ -s "${r1}" ] && SPADES_OPTS="\$SPADES_OPTS -1 ${r1} -2 ${r2}"
      [ -s "${single}" ] && [ "${single.name}" != ".placeholder_single" ] && SPADES_OPTS="\$SPADES_OPTS -s ${single}"
      spades.py \$SPADES_OPTS
      check_spades_coverage "assembly_dir/spades/spades.log" "${sample_id}"
      copy_spades_outputs "assembly_dir/spades"
    fi
    """
}

// ---------------------------------------------------------------------------
// ASSEMBLE_LONG – Flye for long reads only
// ---------------------------------------------------------------------------
process ASSEMBLE_LONG {
    tag { sample_id }
    label "assemble"
    publishDir { "${params.outdir}/${sample_id}/01_assembly" }, mode: "copy", pattern: "contigs*"
    publishDir { "${params.outdir}/${sample_id}/01_assembly" }, mode: "copy", pattern: "scaffolds"

    input:
    tuple val(sample_id), path(long_reads), val(long_tech)
    path(reference_file)
    path(assembly_utils_script)

    output:
    tuple val(sample_id), path("contigs"), path("scaffolds"), emit: assembled
    path("contigs_ref_guided.fasta"), optional: true

    script:
    def has_ref = reference_file.name != '.placeholder_ref'
    def ref_fasta = has_ref ? "reference.fasta" : ''
    """
    source ${assembly_utils_script}
    mkdir -p assembly_dir

    if [ -n "${ref_fasta}" ] && [ -f "${reference_file}" ]; then
      ln -sf "${reference_file}" "${ref_fasta}"
    fi

    if [ -s "${long_reads}" ] && [ "${long_reads.name}" != ".placeholder_long" ]; then
      FLYE_EXTRA="--meta --min-overlap ${params.flye_min_overlap}"
      [ -n "${params.flye_genome_size ?: ''}" ] && FLYE_EXTRA="\$FLYE_EXTRA --genome-size ${params.flye_genome_size}"

      if [ "${long_tech}" = "pacbio" ]; then
        FLYE_INPUT_FLAG="--pacbio-raw"
      else
        FLYE_INPUT_FLAG="--nano-raw"
      fi

      flye \$FLYE_INPUT_FLAG ${long_reads} --out-dir assembly_dir/flye \$FLYE_EXTRA -t ${task.cpus}
      check_flye_coverage "assembly_dir/flye/flye.log" "${sample_id}"
      cp assembly_dir/flye/assembly.fasta contigs 2>/dev/null || true

      # Polish Nanopore assemblies with Medaka
      if [ "${long_tech}" = "nanopore" ] && [ -s contigs ]; then
        echo "Polishing Nanopore assembly with Medaka..."
        medaka_consensus -i ${long_reads} -d contigs -o assembly_dir/medaka -t ${task.cpus} \
          && cp assembly_dir/medaka/consensus.fasta contigs \
          || echo "WARNING: Medaka polishing failed, using unpolished Flye assembly"
      fi

      # Reference-guided consensus
      if [ -n "${ref_fasta}" ] && [ -f "${ref_fasta}" ]; then
        run_ref_guided_consensus "${long_reads}" "${ref_fasta}" "${long_tech}" "${task.cpus}" "contigs_ref_guided.fasta" "assembly_dir"
      fi

      cp contigs scaffolds 2>/dev/null || true
    else
      echo "ASSEMBLE_LONG: no long reads for ${sample_id}"
      touch contigs scaffolds
    fi
    """
}

// ---------------------------------------------------------------------------
// ASSEMBLE_HYBRID – Flye + short-read polishers
// ---------------------------------------------------------------------------
process ASSEMBLE_HYBRID {
    tag { sample_id }
    label "assemble"
    publishDir { "${params.outdir}/${sample_id}/01_assembly" }, mode: "copy", pattern: "contigs*"
    publishDir { "${params.outdir}/${sample_id}/01_assembly" }, mode: "copy", pattern: "scaffolds"

    input:
    tuple val(sample_id), path(r1), path(r2), path(single), path(long_reads), val(short_tech), val(long_tech)
    path(reference_file)
    path(assembly_utils_script)

    output:
    tuple val(sample_id), path("contigs"), path("scaffolds"), emit: assembled
    path("contigs_ref_guided.fasta"), optional: true

    script:
    def mem = task.memory.toGiga().toString()
    def ion = (short_tech == "iontorrent") ? '--iontorrent' : ''
    def has_ref = reference_file.name != '.placeholder_ref'
    def ref_fasta = has_ref ? "reference.fasta" : ''
    def ref = has_ref ? "--trusted-contigs ${ref_fasta}" : ''
    def rna = (params.rna_mode && !has_ref) ? '--rnaviral' : ''
    def metaviral = (params.metaviral_mode && !has_ref) ? '--metaviral' : ''
    """
    source ${assembly_utils_script}
    mkdir -p assembly_dir

    if [ -n "${ref_fasta}" ] && [ -f "${reference_file}" ]; then
      ln -sf "${reference_file}" "${ref_fasta}"
    fi
    if [ -n "${ref}" ]; then
      echo "NOTE: Reference genome provided — using --trusted-contigs for reference-guided assembly."
      if [ "${params.rna_mode}" = "true" ] || [ "${params.metaviral_mode}" = "true" ]; then
        echo "      Skipping --rnaviral/--metaviral (incompatible with --trusted-contigs)."
      fi
    fi

    if [ -s "${long_reads}" ] && [ "${long_reads.name}" != ".placeholder_long" ]; then
      echo "Hybrid assembly: Flye --meta → Medaka (Nanopore only) → Polypolish → Pypolca"
      FLYE_EXTRA="--meta --min-overlap ${params.flye_min_overlap}"
      [ -n "${params.flye_genome_size ?: ''}" ] && FLYE_EXTRA="\$FLYE_EXTRA --genome-size ${params.flye_genome_size}"

      if [ "${long_tech}" = "pacbio" ]; then
        FLYE_INPUT_FLAG="--pacbio-raw"
      else
        FLYE_INPUT_FLAG="--nano-raw"
      fi

      flye \$FLYE_INPUT_FLAG ${long_reads} --out-dir assembly_dir/flye \$FLYE_EXTRA -t ${task.cpus}
      check_flye_coverage "assembly_dir/flye/flye.log" "${sample_id}"
      cp assembly_dir/flye/assembly.fasta contigs 2>/dev/null || true

      # LR polishing: Medaka fixes systematic Nanopore errors (skip for PacBio)
      if [ "${long_tech}" = "nanopore" ] && [ -s contigs ]; then
        echo "Polishing Nanopore assembly with Medaka..."
        medaka_consensus -i ${long_reads} -d contigs -o assembly_dir/medaka -t ${task.cpus} \
          && cp assembly_dir/medaka/consensus.fasta contigs \
          || echo "WARNING: Medaka polishing failed, continuing with unpolished assembly"
      fi

      HAS_SHORT=\$( [ -s "${r1}" ] && [ -s "${r2}" ] && echo 1 || echo 0 )
      HAS_SINGLE=\$( [ -s "${single}" ] && [ "${single.name}" != ".placeholder_single" ] && echo 1 || echo 0 )

      # SR polishing: Polypolish
      if [ -s contigs ] && ([ "\$HAS_SHORT" = "1" ] || [ "\$HAS_SINGLE" = "1" ]); then
        if command -v polypolish &>/dev/null; then
          echo "Polishing assembly with short reads (Polypolish)..."
          bwa index contigs
          if [ "\$HAS_SHORT" = "1" ]; then
            bwa mem -t ${task.cpus} -a contigs ${r1} > assembly_dir/alignments_1.sam
            bwa mem -t ${task.cpus} -a contigs ${r2} > assembly_dir/alignments_2.sam
            polypolish polish contigs assembly_dir/alignments_1.sam assembly_dir/alignments_2.sam > assembly_dir/polished.fasta \
              && cp assembly_dir/polished.fasta contigs \
              || echo "WARNING: Polypolish failed, using previous assembly"
            rm -f assembly_dir/alignments_1.sam assembly_dir/alignments_2.sam
          elif [ "\$HAS_SINGLE" = "1" ]; then
            bwa mem -t ${task.cpus} -a contigs ${single} > assembly_dir/alignments.sam
            polypolish polish contigs assembly_dir/alignments.sam > assembly_dir/polished.fasta \
              && cp assembly_dir/polished.fasta contigs \
              || echo "WARNING: Polypolish failed, using previous assembly"
            rm -f assembly_dir/alignments.sam
          fi
        else
          echo "WARNING: Polypolish not found – skipping. Install with: conda install -c bioconda polypolish"
        fi
      fi

      # Final polish: Pypolca
      if [ -s contigs ] && [ "\$HAS_SHORT" = "1" ]; then
        if command -v pypolca &>/dev/null; then
          echo "Final polishing with Pypolca..."
          pypolca run -a contigs \
            -1 ${r1} \
            -2 ${r2} \
            -t ${task.cpus} --careful \
            -o assembly_dir/pypolca 2>&1 \
            && { PYPOLCA_OUT=\$(ls assembly_dir/pypolca/*corrected*.fasta 2>/dev/null | head -1); \
                 [ -s "\$PYPOLCA_OUT" ] && cp "\$PYPOLCA_OUT" contigs; } \
            || echo "WARNING: Pypolca failed, using Polypolish-polished assembly"
        else
          echo "WARNING: Pypolca not found – skipping final polish. Install with: pip install pypolca"
        fi
      fi

      # Reference-guided consensus
      if [ -n "${ref_fasta}" ] && [ -f "${ref_fasta}" ]; then
        run_ref_guided_consensus "${long_reads}" "${ref_fasta}" "${long_tech}" "${task.cpus}" "contigs_ref_guided.fasta" "assembly_dir"
      fi

      cp contigs scaffolds 2>/dev/null || true
    else
      echo "Warning: strategy=hybrid but no long reads available, falling back to short-read assembly"
      if [ -s "${r1}" ] || [ -s "${single}" ]; then
        SPADES_OPTS="-o assembly_dir/spades -t ${task.cpus} -m ${mem} --only-assembler ${rna} ${ion} ${metaviral} ${ref}"
        [ -s "${r1}" ] && SPADES_OPTS="\$SPADES_OPTS -1 ${r1} -2 ${r2}"
        [ -s "${single}" ] && [ "${single.name}" != ".placeholder_single" ] && SPADES_OPTS="\$SPADES_OPTS -s ${single}"
        spades.py \$SPADES_OPTS
        check_spades_coverage "assembly_dir/spades/spades.log" "${sample_id}"
        copy_spades_outputs "assembly_dir/spades"
      else
        touch contigs scaffolds
      fi
    fi
    """
}

// ---------------------------------------------------------------------------
// REF_ASSEMBLE – reference-guided consensus only (no de novo assembly)
// ---------------------------------------------------------------------------
process REF_ASSEMBLE {
    tag { sample_id }
    label "assemble"
    publishDir { "${params.outdir}/${sample_id}/01_assembly" }, mode: "copy", pattern: "contigs*"
    publishDir { "${params.outdir}/${sample_id}/01_assembly" }, mode: "copy", pattern: "scaffolds"

    input:
    tuple val(sample_id),
          path(r1),
          path(r2),
          path(single),
          path(long_reads),
          val(short_tech),
          val(long_tech)
    path(reference_file)

    output:
    tuple val(sample_id), path("contigs"), path("scaffolds"), emit: assembled

    script:
    def minimap2_lr_preset = (long_tech == "pacbio") ? "map-pb" : "map-ont"
    """
    mkdir -p assembly_dir preprocess_dir

    ln -sf "${reference_file}" reference.fasta

    cp "${r1}" preprocess_dir/trimmed_R1.fastq.gz 2>/dev/null || true
    cp "${r2}" preprocess_dir/trimmed_R2.fastq.gz 2>/dev/null || true
    cp "${single}" preprocess_dir/trimmed_single.fastq.gz 2>/dev/null || true
    cp "${long_reads}" preprocess_dir/trimmed_long.fastq.gz 2>/dev/null || true

    HAS_SHORT=\$( [ -s preprocess_dir/trimmed_R1.fastq.gz ] && [ -s preprocess_dir/trimmed_R2.fastq.gz ] && echo 1 || echo 0 )
    HAS_SINGLE=\$( [ -s preprocess_dir/trimmed_single.fastq.gz ] && echo 1 || echo 0 )
    HAS_LONG=\$( [ -s preprocess_dir/trimmed_long.fastq.gz ] && echo 1 || echo 0 )

    echo "Reference-only mode: generating consensus from reference alignment..."

    # Align reads to reference
    if [ "\$HAS_SHORT" = "1" ] || [ "\$HAS_SINGLE" = "1" ]; then
      bwa index reference.fasta
      if [ "\$HAS_SHORT" = "1" ]; then
        bwa mem -t ${task.cpus} reference.fasta \
          preprocess_dir/trimmed_R1.fastq.gz preprocess_dir/trimmed_R2.fastq.gz 2>/dev/null | \
          samtools view -q 30 -F 256 -b - 2>/dev/null | \
          samtools sort -@ ${task.cpus} -o assembly_dir/ref_short_pe.bam -
      fi
      if [ "\$HAS_SINGLE" = "1" ]; then
        bwa mem -t ${task.cpus} reference.fasta \
          preprocess_dir/trimmed_single.fastq.gz 2>/dev/null | \
          samtools view -q 30 -F 256 -b - 2>/dev/null | \
          samtools sort -@ ${task.cpus} -o assembly_dir/ref_short_single.bam -
      fi
    fi

    if [ "\$HAS_LONG" = "1" ]; then
      minimap2 -t ${task.cpus} -ax ${minimap2_lr_preset} reference.fasta \
        preprocess_dir/trimmed_long.fastq.gz 2>/dev/null | \
        samtools sort -@ ${task.cpus} -o assembly_dir/ref_long.bam -
    fi

    # Merge all BAMs
    BAMS=""
    [ -f assembly_dir/ref_short_pe.bam ] && BAMS="\$BAMS assembly_dir/ref_short_pe.bam"
    [ -f assembly_dir/ref_short_single.bam ] && BAMS="\$BAMS assembly_dir/ref_short_single.bam"
    [ -f assembly_dir/ref_long.bam ] && BAMS="\$BAMS assembly_dir/ref_long.bam"

    BAM_COUNT=\$(echo \$BAMS | wc -w)
    if [ "\$BAM_COUNT" -gt 1 ]; then
      samtools merge -f assembly_dir/ref_merged.bam \$BAMS
    elif [ "\$BAM_COUNT" -eq 1 ]; then
      mv \$BAMS assembly_dir/ref_merged.bam
    else
      echo "REF_ASSEMBLE: no reads aligned to reference for ${sample_id}"
      touch contigs scaffolds
      exit 0
    fi
    samtools index assembly_dir/ref_merged.bam

    # Generate consensus
    LONG_READ_TECH="${long_tech}"
    if [ "\$HAS_LONG" = "1" ] && [ "\$LONG_READ_TECH" = "nanopore" ]; then
      medaka_consensus -i preprocess_dir/trimmed_long.fastq.gz \
        -d reference.fasta -o assembly_dir/ref_medaka -t ${task.cpus} \
        && cp assembly_dir/ref_medaka/consensus.fasta contigs \
        || { echo "WARNING: Medaka failed, falling back to samtools consensus"; \
             samtools consensus assembly_dir/ref_merged.bam -o contigs --show-ins no -a; }
    else
      samtools consensus assembly_dir/ref_merged.bam -o contigs --show-ins no -a
    fi

    cp contigs scaffolds 2>/dev/null || touch contigs scaffolds
    touch preprocess_dir/trimmed_R1.fastq.gz preprocess_dir/trimmed_R2.fastq.gz preprocess_dir/trimmed_single.fastq.gz preprocess_dir/trimmed_long.fastq.gz

    echo "Reference-only assembly complete for ${sample_id}"
    """
}


