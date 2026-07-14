/*
 * Virall Nextflow pipeline — HOST_FILTER module
 * Optional host read filtering via minimap2 + samtools
 */

process HOST_FILTER {
    tag { sample_id }
    label "host_filter"
    publishDir { "${params.outdir}/${sample_id}/00_preprocess/host_filtered" }, mode: "symlink"

    input:
    tuple val(sample_id),
          path(r1),
          path(r2),
          path(single),
          path(long_reads),
          val(short_tech),
          val(long_tech)
    path(host_ref)

    output:
    tuple val(sample_id),
          path("host_filter_dir/trimmed_R1.fastq.gz"),
          path("host_filter_dir/trimmed_R2.fastq.gz"),
          path("host_filter_dir/trimmed_single.fastq.gz"),
          path("host_filter_dir/trimmed_long.fastq.gz"),
          emit: reads

    script:
    """
    mkdir -p host_filter_dir

    # Paired-end: keep pairs where BOTH reads are unmapped (-f 12)
    if [ -s "${r1}" ] && [ -s "${r2}" ] && [ "${r1.name}" != ".placeholder_r1" ]; then
      minimap2 -ax sr -t ${task.cpus} ${host_ref} ${r1} ${r2} 2>host_filter_dir/host_filter_pe.log | \
        samtools sort -n - 2>/dev/null | samtools fastq -f 12 -1 host_filter_dir/R1.fq -2 host_filter_dir/R2.fq - 2>/dev/null || { echo "WARNING: samtools fastq (PE) failed for ${sample_id}" >&2; }
      if [ -s host_filter_dir/R1.fq ]; then
        gzip -c host_filter_dir/R1.fq > host_filter_dir/trimmed_R1.fastq.gz
        gzip -c host_filter_dir/R2.fq > host_filter_dir/trimmed_R2.fastq.gz
      else
        cp ${r1} host_filter_dir/trimmed_R1.fastq.gz
        cp ${r2} host_filter_dir/trimmed_R2.fastq.gz
      fi
    else
      cp ${r1} host_filter_dir/trimmed_R1.fastq.gz 2>/dev/null || touch host_filter_dir/trimmed_R1.fastq.gz
      cp ${r2} host_filter_dir/trimmed_R2.fastq.gz 2>/dev/null || touch host_filter_dir/trimmed_R2.fastq.gz
    fi

    # Single-end short
    if [ -s "${single}" ] && [ "${single.name}" != ".placeholder_single" ]; then
      minimap2 -ax sr -t ${task.cpus} ${host_ref} ${single} 2>host_filter_dir/host_filter_single.log | \
        samtools fastq -f 4 - > host_filter_dir/single.fq 2>/dev/null || { echo "WARNING: samtools fastq (single) failed for ${sample_id}" >&2; }
      if [ -s host_filter_dir/single.fq ]; then
        gzip -c host_filter_dir/single.fq > host_filter_dir/trimmed_single.fastq.gz
      else
        cp ${single} host_filter_dir/trimmed_single.fastq.gz
      fi
    else
      cp ${single} host_filter_dir/trimmed_single.fastq.gz 2>/dev/null || touch host_filter_dir/trimmed_single.fastq.gz
    fi

    # Long reads (ONT or PacBio)
    MINIMAP2_LONG_PRESET=\$( [ "${long_tech}" = "pacbio" ] && echo "map-pb" || echo "map-ont" )
    if [ -s "${long_reads}" ] && [ "${long_reads.name}" != ".placeholder_long" ]; then
      minimap2 -ax \$MINIMAP2_LONG_PRESET -t ${task.cpus} ${host_ref} ${long_reads} 2>host_filter_dir/host_filter_long.log | \
        samtools fastq -f 4 - > host_filter_dir/long.fq 2>/dev/null || { echo "WARNING: samtools fastq (long) failed for ${sample_id}" >&2; }
      if [ -s host_filter_dir/long.fq ]; then
        gzip -c host_filter_dir/long.fq > host_filter_dir/trimmed_long.fastq.gz
      else
        cp ${long_reads} host_filter_dir/trimmed_long.fastq.gz
      fi
    else
      cp ${long_reads} host_filter_dir/trimmed_long.fastq.gz 2>/dev/null || touch host_filter_dir/trimmed_long.fastq.gz
    fi
    touch host_filter_dir/trimmed_R1.fastq.gz host_filter_dir/trimmed_R2.fastq.gz host_filter_dir/trimmed_single.fastq.gz host_filter_dir/trimmed_long.fastq.gz
    """

}


