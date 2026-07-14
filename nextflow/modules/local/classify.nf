/*
 * Virall Nextflow pipeline -- CLASSIFY module
 * KAIJU (taxonomic classification), FILTER_VIRAL, RENAME_CONTIGS
 */

// ---------------------------------------------------------------------------
// KAIJU – classify all contigs (contigs.fasta used for viral filtering)
// ---------------------------------------------------------------------------
process KAIJU {
    tag { sample_id }
    label "kaiju"
    publishDir { "${params.outdir}/${sample_id}/03_classifications" }, mode: "symlink"

    input:
    tuple val(sample_id), path(contigs)
    val(kaiju_db)

    output:
    tuple val(sample_id), path("kaiju_dir"), emit: kaiju_done

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
    tag { sample_id }
    label "filter_viral"
    publishDir { "${params.outdir}/${sample_id}/02_viral_contigs" }, mode: "symlink"

    input:
    tuple val(sample_id), path(kaiju_dir), path(contigs)
    path(extract_script)

    output:
    tuple val(sample_id), path("viral_contigs.fasta"), path(kaiju_dir), emit: viral

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
    tag { sample_id }
    label "rename_contigs"
    publishDir { "${params.outdir}/${sample_id}/02_viral_contigs" }, mode: "copy", pattern: "{viral_contigs.fasta,name_mapping.tsv}"

    input:
    tuple val(sample_id), path("input_viral_contigs.fasta"), path("input_kaiju_dir")

    output:
    tuple val(sample_id), path("viral_contigs.fasta"), path("kaiju_dir"), emit: renamed

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

