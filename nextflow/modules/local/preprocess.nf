/*
 * Virall Nextflow pipeline — PREPROCESS module
 * fastp (short) / fastplong (long) read preprocessing
 */

process PREPROCESS {
    tag { sample_id }
    label "preprocess"
    publishDir { "${params.outdir}/${sample_id}/00_preprocess" }, mode: "symlink"

    input:
    tuple val(sample_id),
          path(read1),
          path(read2),
          path(single),
          path(long_reads),
          val(short_tech),
          val(long_tech)

    output:
    tuple val(sample_id),
          path("preprocess_dir/trimmed_R1.fastq.gz"),
          path("preprocess_dir/trimmed_R2.fastq.gz"),
          path("preprocess_dir/trimmed_single.fastq.gz"),
          path("preprocess_dir/trimmed_long.fastq.gz"),
          emit: reads
    tuple val(sample_id),
          path("preprocess_dir/fastp_pe.html"),
          path("preprocess_dir/fastp_pe.json"),
          path("preprocess_dir/fastp_single.html"),
          path("preprocess_dir/fastp_single.json"),
          path("preprocess_dir/fastplong.html"),
          path("preprocess_dir/fastplong.json"),
          emit: qc

    script:
    def hasShort = (read1.name != ".placeholder_r1" && read1.size() > 0) || (single.name != ".placeholder_single" && single.size() > 0)
    def hasLong  = long_reads.name != ".placeholder_long" && long_reads.size() > 0
    def q = params.quality_phred
    def ml = params.min_read_len
    def ion_opts = (short_tech == "iontorrent") ? "--disable_adapter_trimming --disable_trim_poly_g" : "--trim_poly_g"
    def pe_adapter = (short_tech == "iontorrent") ? "" : "--detect_adapter_for_pe"
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
