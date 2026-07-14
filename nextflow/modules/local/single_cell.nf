/*
 * Virall Nextflow pipeline -- SINGLE-CELL module
 * SC_EXTRACT_BARCODES, SC_POOL_READS,
 * SC_MAP_VIRAL, SC_COUNT_CELLS, SC_BUILD_MATRIX
 */

params.sc_min_reads_per_cell = 100    // Minimum reads to keep a cell
params.sc_min_viral_umis = 1          // Minimum viral UMIs to call cell "infected"

// ---------------------------------------------------------------------------
// SINGLE-CELL: SC_EXTRACT_BARCODES – Extract cell barcodes and UMIs from 10x data
// ---------------------------------------------------------------------------
process SC_EXTRACT_BARCODES {
    tag {"SC_BC: $sample_id"}
    label "single_cell"
    publishDir { "${params.outdir}/${sample_id}/08_single_cell/barcoded_reads" }, mode: "symlink"

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
    def chem = params.sc_chemistry
    if (chem == "10x_v3" || chem == "10x_3prime_v3" || chem == "10x_5prime_v3") {
        bc_pattern = "CCCCCCCCCCCCCCCCNNNNNNNNNNNN"  // 16C + 12N
    } else if (chem == "10x_v2" || chem == "10x_3prime_v2" || chem == "10x_5prime_v2") {
        bc_pattern = "CCCCCCCCCCCCCCCCNNNNNNNNNN"    // 16C + 10N
    } else if (chem == "10x_v1" || chem == "10x_3prime_v1") {
        bc_pattern = "CCCCCCCCCCCCCCNNNNNNNNNN"      // 14C + 10N
    } else {
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
    tag { sample_id }
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

process SC_MAP_VIRAL {
    tag { sample_id }
    label "single_cell"
    publishDir { "${params.outdir}/${sample_id}/08_single_cell/viral_mapping" }, mode: "symlink"

    input:
    tuple val(sample_id), path(tagged_r2), path(viral_contigs)

    output:
    tuple val(sample_id), path("viral_aligned.bam"), path("viral_aligned.bam.bai"), emit: bam
    tuple val(sample_id), path("mapping_stats.txt"), emit: stats

    script:
    """
    if [ ! -s ${viral_contigs} ]; then
        echo "WARNING: No viral contigs found, skipping SC_MAP_VIRAL for ${sample_id}" >&2
        printf '@HD\tVN:1.6\tSO:coordinate\n' | samtools view -bS -o viral_aligned.bam -
        samtools index viral_aligned.bam
        echo "## No viral contigs in input" > mapping_stats.txt
        exit 0
    fi

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
    tag { sample_id }
    label "single_cell"
    publishDir { "${params.outdir}/${sample_id}/08_single_cell" }, mode: "symlink"

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
    tag { sample_id }
    label "single_cell"
    publishDir { "${params.outdir}/${sample_id}/08_single_cell/matrix" }, mode: "copy"

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

