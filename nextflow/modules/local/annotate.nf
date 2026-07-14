/*
 * Virall Nextflow pipeline -- ANNOTATE module
 * ANNOTATE (prodigal-gv + hmmscan VOG), ORGANIZE_GENES
 */

// ANNOTATE – prodigal-gv + hmmscan (VOG)
// ---------------------------------------------------------------------------
process ANNOTATE {
    tag { sample_id }
    label "annotate"
    publishDir { "${params.outdir}/${sample_id}/05_gene_predictions" }, mode: "symlink"

    input:
    tuple val(sample_id), path(viral_contigs), path(kaiju_dir)
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
        hmmsearch --cpu ${task.cpus} -o annotation_dir/vog_out.txt --domtblout annotation_dir/vog_domains.txt "\$VOG_HMM" annotation_dir/proteins.faa
      fi
    fi
    """
}

// ---------------------------------------------------------------------------
// ORGANIZE_GENES – Organize genes by Kaiju taxonomy + VOG functional annotation
// ---------------------------------------------------------------------------
process ORGANIZE_GENES {
    tag { sample_id }
    label "annotate"
    publishDir { "${params.outdir}/${sample_id}/05_gene_predictions" }, mode: "symlink"

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

