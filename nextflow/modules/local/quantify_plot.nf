/*
 * Virall Nextflow pipeline -- QUANTIFY + PLOT module
 * QUANTIFY (BWA/minimap2 + samtools depth), REFERENCE_CHECK, PLOT
 */

// QUANTIFY – BWA/minimap2 + samtools depth
// ---------------------------------------------------------------------------
process QUANTIFY {
    tag { sample_id }
    label "quantify"
    publishDir { "${params.outdir}/${sample_id}/06_quantification" }, mode: "copy"

    input:
    tuple val(sample_id), path(viral_contigs), path(kaiju_dir), path(r1), path(r2), path(single), path(long_reads), val(long_tech)

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
      if [ -s "${r1}" ] && [ -s "${r2}" ] && [ "${r1.name}" != ".placeholder_r1" ]; then
        bwa mem -t ${task.cpus} quant_dir/idx ${r1} ${r2} | samtools sort -o quant_dir/mapped_pe.bam -
        BAMS_TO_MERGE="\$BAMS_TO_MERGE quant_dir/mapped_pe.bam"
      fi
      if [ -s "${single}" ] && [ "${single.name}" != ".placeholder_single" ]; then
        bwa mem -t ${task.cpus} quant_dir/idx ${single} | samtools sort -o quant_dir/mapped_single.bam -
        BAMS_TO_MERGE="\$BAMS_TO_MERGE quant_dir/mapped_single.bam"
      fi
      HAS_SHORT=\$([ -n "\$BAMS_TO_MERGE" ] && echo "1" || echo "0")
      if [ "\$HAS_SHORT" = "0" ] && [ -s "${long_reads}" ] && [ "${long_reads.name}" != ".placeholder_long" ]; then
        MINIMAP2_LONG_PRESET=\$( [ "${long_tech}" = "pacbio" ] && echo "map-pb" || echo "map-ont" )
        minimap2 -t ${task.cpus} -ax \$MINIMAP2_LONG_PRESET viral_contigs.fasta ${long_reads} | samtools sort -o quant_dir/mapped_long.bam -
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
      samtools depth -a quant_dir/mapped.bam > quant_dir/depth.txt 2>/dev/null || { echo "WARNING: samtools depth failed for ${sample_id}" >&2; }

      # High-confidence subset for abundance estimation:
      # -q: keep alignments with MAPQ >= params.quant_mapq
      # -F 2308: remove unmapped(4), secondary(256), supplementary(2048)
      samtools view -b -q ${params.quant_mapq} -F 2308 quant_dir/mapped.bam > quant_dir/mapped_primary_mapq.bam 2>/dev/null || { echo "WARNING: samtools view MAPQ filter failed for ${sample_id}" >&2; }
      if [ -s quant_dir/mapped_primary_mapq.bam ]; then
        samtools index quant_dir/mapped_primary_mapq.bam 2>/dev/null || true
        samtools depth -a quant_dir/mapped_primary_mapq.bam > quant_dir/depth_primary_mapq${params.quant_mapq}.txt 2>/dev/null || { echo "WARNING: samtools depth (MAPQ) failed for ${sample_id}" >&2; }
      else
        touch quant_dir/mapped_primary_mapq.bam quant_dir/depth_primary_mapq${params.quant_mapq}.txt
      fi
    fi
    """
}

// ---------------------------------------------------------------------------
// REFERENCE_CHECK – optional: map reads to reference, report detection, plot coverage
// ---------------------------------------------------------------------------
process REFERENCE_CHECK {
    tag { sample_id }
    label "reference_check"
    publishDir { "${params.outdir}/${sample_id}/08_reference_check" }, mode: "symlink"

    input:
    tuple val(sample_id),
          path(r1),
          path(r2),
          path(single),
          path(long_reads),
          val(short_tech),
          val(long_tech)
    path(reference)
    path(ref_check_plot_script)

    output:
    tuple val(sample_id), path("ref_check_dir"), emit: ref_checked

    script:
    def minimap2_preset = (long_tech == "pacbio") ? "map-pb" : "map-ont"
    """
    mkdir -p ref_check_dir

    # Index reference
    bwa index -p ref_check_dir/ref_idx ${reference} 2>/dev/null || true

    # Map reads to reference and collect stats
    MAPPED_READS=0
    TOTAL_READS=0

    if [ -s "${r1}" ] && [ -s "${r2}" ] && [ "${r1.name}" != ".placeholder_r1" ]; then
      bwa mem -t ${task.cpus} ref_check_dir/ref_idx ${r1} ${r2} 2>/dev/null | \
        samtools sort -o ref_check_dir/mapped_pe.bam - 2>/dev/null
      samtools index ref_check_dir/mapped_pe.bam 2>/dev/null || true
      PE_MAPPED=\$(samtools view -c -F 4 ref_check_dir/mapped_pe.bam 2>/dev/null || echo 0)
      PE_TOTAL=\$(samtools view -c ref_check_dir/mapped_pe.bam 2>/dev/null || echo 0)
      MAPPED_READS=\$((MAPPED_READS + PE_MAPPED))
      TOTAL_READS=\$((TOTAL_READS + PE_TOTAL))
    fi

    if [ -s "${single}" ] && [ "${single.name}" != ".placeholder_single" ]; then
      bwa mem -t ${task.cpus} ref_check_dir/ref_idx ${single} 2>/dev/null | \
        samtools sort -o ref_check_dir/mapped_single.bam - 2>/dev/null
      samtools index ref_check_dir/mapped_single.bam 2>/dev/null || true
      SINGLE_MAPPED=\$(samtools view -c -F 4 ref_check_dir/mapped_single.bam 2>/dev/null || echo 0)
      SINGLE_TOTAL=\$(samtools view -c ref_check_dir/mapped_single.bam 2>/dev/null || echo 0)
      MAPPED_READS=\$((MAPPED_READS + SINGLE_MAPPED))
      TOTAL_READS=\$((TOTAL_READS + SINGLE_TOTAL))
    fi

    if [ -s "${long_reads}" ] && [ "${long_reads.name}" != ".placeholder_long" ]; then
      minimap2 -t ${task.cpus} -ax ${minimap2_preset} ${reference} ${long_reads} 2>/dev/null | \
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
      samtools depth -a ref_check_dir/merged.bam > ref_check_dir/reference_depth.txt 2>/dev/null || { echo "WARNING: samtools depth (ref) failed for ${sample_id}" >&2; }
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
      python3 ${ref_check_plot_script} \
        --depth ref_check_dir/reference_depth.txt \
        --out-prefix ref_check_dir/reference || true
    fi

    # Append per-segment stats to report if available
    if [ -f ref_check_dir/reference_stats.tsv ]; then
      SEG_COUNT=\$(tail -n +2 ref_check_dir/reference_stats.tsv | wc -l)
      if [ "\$SEG_COUNT" -gt 1 ]; then
        echo "Per-Segment Statistics:" >> ref_check_dir/reference_report.txt
        echo "  Segment                Length   Covered  Breadth%  MeanDepth" >> ref_check_dir/reference_report.txt
        tail -n +2 ref_check_dir/reference_stats.tsv | while IFS=\$'\\t' read -r seg slen cov br md; do
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
    tag { sample_id }
    label "plot"
    publishDir { "${params.outdir}/${sample_id}/07_plots" }, mode: "copy"

    input:
    tuple val(sample_id), path(quant_dir), path(merged_quality), path(viral_contigs), path(kaiju_dir)
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
        --checkv-dir ${merged_quality} \\
        --min-breadth ${params.quant_min_breadth} \\
        --out-dir plots_dir
    fi
    """
}

// ---------------------------------------------------------------------------
// SINGLE-CELL: SC_MAP_VIRAL – Map barcoded reads to viral contigs
// ---------------------------------------------------------------------------

