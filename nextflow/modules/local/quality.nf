/*
 * Virall Nextflow pipeline -- QUALITY module
 * VALIDATE (CheckV), GENOMAD, MERGE_QUALITY
 */

// VALIDATE – CheckV on viral contigs
// ---------------------------------------------------------------------------
process VALIDATE {
    tag { sample_id }
    label "validate"
    publishDir { "${params.outdir}/${sample_id}/04_quality_assessment" }, mode: "symlink"

    input:
    tuple val(sample_id), path(viral_contigs), path(kaiju_dir)
    val(checkv_db)

    output:
    tuple val(sample_id), path("checkv_dir"), emit: validated

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
    tag { sample_id }
    label "genomad"
    publishDir { "${params.outdir}/${sample_id}/04_quality_assessment/genomad" }, mode: "symlink"

    input:
    tuple val(sample_id), path(viral_contigs), path(kaiju_dir)
    val(genomad_db)

    output:
    tuple val(sample_id), path("genomad_out"), emit: genomad_done

    script:
    """
    mkdir -p genomad_out

    # Find geNomad database
    GENOMAD_D=""
    for cand in "${genomad_db}/genomad_db" "${genomad_db}"; do
      if [ -d "\$cand" ] && { [ -f "\$cand/genomad_db" ] || [ -f "\$cand/virus_hallmark_annotation.txt" ] || [ -d "\$cand/mmseqs2" ]; }; then
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
    tag { sample_id }
    label "merge_quality"
    publishDir { "${params.outdir}/${sample_id}/04_quality_assessment" }, mode: "copy"

    input:
    tuple val(sample_id), path(checkv_dir), path(genomad_dir), path(viral_contigs), path(kaiju_dir)
    path(merge_script)

    output:
    tuple val(sample_id), path("merged_quality"), emit: merged

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

