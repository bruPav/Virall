/*
 * Virall Nextflow Pipeline – tools directly (no Virall Python wrappers).
 * Each process calls underlying tools (fastp, SPAdes, Flye, Kaiju, CheckV, Prodigal, BWA, etc.).
 * Rules (strategy, params) are driven by user params.
 */

nextflow.enable.dsl = 2

include { PREPROCESS } from './modules/local/preprocess'
include { HOST_FILTER } from './modules/local/host_filter'
include { ASSEMBLE_SHORT; ASSEMBLE_LONG; ASSEMBLE_HYBRID; REF_ASSEMBLE } from './modules/local/assembly'
include { KAIJU; FILTER_VIRAL; RENAME_CONTIGS } from './modules/local/classify'
include { VALIDATE; GENOMAD; MERGE_QUALITY } from './modules/local/quality'
include { ANNOTATE; ORGANIZE_GENES } from './modules/local/annotate'
include { QUANTIFY; REFERENCE_CHECK; PLOT } from './modules/local/quantify_plot'
include { SC_EXTRACT_BARCODES; SC_POOL_READS; SC_MAP_VIRAL; SC_COUNT_CELLS; SC_BUILD_MATRIX } from './modules/local/single_cell'

// ---------------------------------------------------------------------------
// Helper: resolve ~ in paths without mutating params
// ---------------------------------------------------------------------------
def resolvePath(String path) {
    if (path?.trim()?.startsWith('~')) {
        def home = System.getenv('HOME') ?: System.getProperty('user.home') ?: ''
        return (home ?: '') + path.trim().substring(1)
    }
    return path
}

// ---------------------------------------------------------------------------
// Parameters (user requirements)
// ---------------------------------------------------------------------------
params.samples         = "${projectDir}/samples.csv"
params.outdir          = "results"
params.threads         = 8
params.memory          = "16G"
params.reference       = null
params.reference_only  = false  // skip de novo assembly; use reference-guided consensus only
params.host_genome     = null   // optional: FASTA of host genome to filter out host reads (like Virall --filter)
params.rna_mode        = false
params.assembly_strategy = "auto"     // auto (detect from input), hybrid, short_only, long_only
params.min_contig_len  = 1000
params.quality_phred   = 20        // short reads (fastp)
params.quality_phred_long = 7     // long reads (fastplong); Virall default (lower for ONT/PacBio)
params.min_read_len    = 50       // short reads (fastp)
params.min_read_len_long = 1000   // long reads (fastplong); Virall default
params.long_read_tech  = "nanopore"  // "nanopore" or "pacbio" – affects SPAdes, Flye, minimap2 presets
params.quant_mapq      = 20       // MAPQ threshold for high-confidence abundance estimation
params.quant_min_breadth = 0.10   // Minimum contig coverage breadth (fraction 0-1) for abundance estimation
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

workflow {
    if (params.reference_only && !params.reference) {
        error "ERROR: reference_only=true requires a reference genome. Set 'reference' in your params file."
    }

    // Resolve DB paths from env if not set (params.xxx are read-only; resolve in workflow)
    def db_dir = System.getenv("VIRALL_DATABASE_DIR") ?: ""

    // ---------------------------------------------------------------------------
    // Helper: resolve ~ in paths
    // ---------------------------------------------------------------------------
    def home_dir = System.getenv('HOME') ?: System.getProperty('user.home') ?: ''
    if (params.outdir?.startsWith('~'))
        params.outdir = home_dir + params.outdir.substring(1)

    // ---------------------------------------------------------------------------
    // Placeholders and sample sheet channel
    // ---------------------------------------------------------------------------
    def stub_r1     = file("${projectDir}/.placeholder_r1")
    def stub_r2     = file("${projectDir}/.placeholder_r2")
    def stub_single = file("${projectDir}/.placeholder_single")
    def stub_long   = file("${projectDir}/.placeholder_long")
    def stub_sc_r1  = file("${projectDir}/.placeholder_sc_r1")
    def stub_sc_r2  = file("${projectDir}/.placeholder_sc_r2")

    Channel
        .fromPath(params.samples, checkIfExists: true)
        .map { file ->
            def raw = file.getText('UTF-8')
            if (raw.startsWith('\uFEFF'))
                raw = raw.substring(1)
            raw
        }
        .splitCsv(header: true, strip: true)
        .map { row ->
            // Map short-read columns to unified format with tech detection
            def read1
            def read2
            def single
            def short_tech
            if (row.illumina_r1?.trim()) {
                read1 = row.illumina_r1
                read2 = row.illumina_r2 ?: ""
                single = row.illumina_single ?: ""
                short_tech = "illumina"
            } else if (row.iontorrent_r1?.trim()) {
                read1 = row.iontorrent_r1
                read2 = row.iontorrent_r2 ?: ""
                single = row.iontorrent_single ?: ""
                short_tech = "iontorrent"
            } else if (row.read1?.trim()) {
                // Backward compatibility: fall back to generic columns
                log.warn "Sample ${row.sample_id}: Using legacy read1/read2/single columns. Consider migrating to tech-specific columns (illumina_r1, iontorrent_r1, etc.)"
                read1 = row.read1
                read2 = row.read2 ?: ""
                single = row.single ?: ""
                short_tech = params.iontorrent ? "iontorrent" : "illumina"
            } else if (row.illumina_single?.trim()) {
                read1 = ""
                read2 = ""
                single = row.illumina_single
                short_tech = "illumina"
            } else if (row.iontorrent_single?.trim()) {
                read1 = ""
                read2 = ""
                single = row.iontorrent_single
                short_tech = "iontorrent"
            } else if (row.single?.trim()) {
                log.warn "Sample ${row.sample_id}: Using legacy single column. Consider migrating to tech-specific columns."
                read1 = ""
                read2 = ""
                single = row.single
                short_tech = params.iontorrent ? "iontorrent" : "illumina"
            } else {
                read1 = ""
                read2 = ""
                single = ""
                short_tech = "illumina"
            }

            // Map long-read columns to unified format with tech detection
            def long_reads
            def long_tech
            if (row.nanopore_reads?.trim()) {
                long_reads = row.nanopore_reads
                long_tech = "nanopore"
            } else if (row.pacbio_reads?.trim()) {
                long_reads = row.pacbio_reads
                long_tech = "pacbio"
            } else if (row.long?.trim()) {
                // Backward compatibility: fall back to generic column
                log.warn "Sample ${row.sample_id}: Using legacy long column. Consider migrating to tech-specific columns (nanopore_reads, pacbio_reads)."
                long_reads = row.long
                long_tech = params.long_read_tech ?: "nanopore"
            } else {
                long_reads = ""
                long_tech = ""
            }

            tuple(
                row.sample_id,
                (read1?.trim()) ? file(read1.trim()) : stub_r1,
                (read2?.trim()) ? file(read2.trim()) : stub_r2,
                (single?.trim()) ? file(single.trim()) : stub_single,
                (long_reads?.trim()) ? file(long_reads.trim()) : stub_long,
                (row.sc_read1?.trim()) ? file(row.sc_read1.trim()) : stub_sc_r1,
                (row.sc_read2?.trim()) ? file(row.sc_read2.trim()) : stub_sc_r2,
                short_tech,
                long_tech
            )
        }
        .set { ch_samples }

    // Default DB base dir: in-container path so -profile docker/singularity works without extra config
    def container_db = "/opt/virall/databases"
    def k_db = params.kaiju_db ?: (db_dir ? "${db_dir}/kaiju_db" : "${container_db}/kaiju_db")
    def c_db = params.checkv_db ?: (db_dir ? "${db_dir}/checkv_db" : "${container_db}/checkv_db")
    def v_db = params.vog_db ?: (db_dir ? "${db_dir}/vog_db" : "${container_db}/vog_db")
    def extract_script = file("${projectDir}/bin/extract_fasta_by_ids.py")
    def merge_script  = file("${projectDir}/bin/merge_quality.py")
    def plot_script   = file("${projectDir}/bin/run_plots_genome_corrected.py")
    def assembly_utils_script = file("${projectDir}/bin/assembly_utils.sh")
    def ref_check_plot_script = file("${projectDir}/bin/reference_check_plot.py")

    def g_db = params.genomad_db ?: (db_dir ? "${db_dir}/genomad_db" : "${container_db}/genomad_db")

    def ch_kaiju_db   = Channel.value(k_db)
    def ch_checkv_db  = Channel.value(c_db)
    def ch_genomad_db = Channel.value(g_db)
    def ch_vog_db     = v_db ? Channel.value(v_db) : Channel.value(null)

    // ---------------------------------------------------------------------------
    // Split sample channel into logical groups
    // ---------------------------------------------------------------------------
    ch_meta      = ch_samples.map { [it[0], it[7], it[8]] }       // sample_id, short_tech, long_tech
    ch_raw_reads = ch_samples.map { [it[0], it[1], it[2], it[3], it[4]] }  // sample_id, r1, r2, single, long
    ch_sc_input  = ch_samples.map { [it[0], it[5], it[6]] }       // sample_id, sc_r1, sc_r2

    // =========================================================================
    // SINGLE-CELL MODE: Extract barcodes, pool, run pipeline, then trace back
    // =========================================================================
    if (params.single_cell_mode) {
        // Extract cell barcodes and UMIs
        SC_EXTRACT_BARCODES(ch_sc_input)
        
        // Pool reads for assembly (barcode info preserved in read names)
        SC_POOL_READS(SC_EXTRACT_BARCODES.out.tagged)
        
        // For 10x single-cell, R2 (cDNA) is used as single-end reads for assembly
        ch_sc_for_preprocess = SC_POOL_READS.out.pooled.map { t -> 
            tuple(
                t[0],                                              // sample_id
                file("${projectDir}/.placeholder_r1"),             // no paired R1
                file("${projectDir}/.placeholder_r2"),             // no paired R2
                t[1],                                              // pooled_single.fastq.gz as single-end
                file("${projectDir}/.placeholder_long")            // no long
            )
        }
        PREPROCESS(ch_sc_for_preprocess.join(ch_meta, remainder: true))
    } else {
        // Standard bulk mode
        PREPROCESS(ch_raw_reads.join(ch_meta, remainder: true))
    }
    ch_trimmed_reads = PREPROCESS.out.reads

    // Optional host filtering: when host_genome is set, filter out host reads
    if (params.host_genome) {
        def host_file = file(resolvePath(params.host_genome))
        HOST_FILTER(ch_trimmed_reads.join(ch_meta, remainder: true), host_file)
        ch_filtered_reads = HOST_FILTER.out.reads
    } else {
        ch_filtered_reads = ch_trimmed_reads
    }

    def ref_file = params.reference ? file(resolvePath(params.reference)) : file("${projectDir}/.placeholder_ref")

    if (params.reference_only && params.reference) {
        REF_ASSEMBLE(ch_filtered_reads.join(ch_meta, remainder: true), ref_file)
        ch_assembly = REF_ASSEMBLE.out.assembled
    } else if (params.assembly_strategy == "short_only") {
        ch_short_input = ch_filtered_reads.join(ch_meta, remainder: true)
            .map { [it[0], it[1], it[2], it[3], it[5]] }
        ASSEMBLE_SHORT(ch_short_input, ref_file, assembly_utils_script)
        ch_assembly = ASSEMBLE_SHORT.out.assembled
    } else if (params.assembly_strategy == "long_only") {
        ch_long_input = ch_filtered_reads.join(ch_meta, remainder: true)
            .map { [it[0], it[4], it[6]] }
        ASSEMBLE_LONG(ch_long_input, ref_file, assembly_utils_script)
        ch_assembly = ASSEMBLE_LONG.out.assembled
    } else if (params.assembly_strategy == "hybrid") {
        ASSEMBLE_HYBRID(ch_filtered_reads.join(ch_meta, remainder: true), ref_file, assembly_utils_script)
        ch_assembly = ASSEMBLE_HYBRID.out.assembled
    } else {
        // auto-detect per sample based on available reads
        ch_reads_meta = ch_filtered_reads.join(ch_meta, remainder: true)

        // classify each sample by available reads → (mode, sample_tuple)
        ch_reads_classified = ch_reads_meta.map { t ->
            def has_pe   = t[1].size() > 0 && t[2].size() > 0
            def has_se   = t[3].size() > 0
            def has_long = t[4].size() > 0
            if (!has_pe && !has_se && !has_long) return tuple("empty", t)
            if ((has_pe || has_se) && !has_long) return tuple("short", t)
            if (!has_pe && !has_se && has_long)  return tuple("long", t)
            return tuple("hybrid", t)
        }

        ch_auto_short = ch_reads_classified
            .filter { it[0] == "short" || it[0] == "empty" }
            .map { [it[1][0], it[1][1], it[1][2], it[1][3], it[1][5]] }
        ch_auto_long = ch_reads_classified
            .filter { it[0] == "long" }
            .map { [it[1][0], it[1][4], it[1][6]] }
        ch_auto_hybrid = ch_reads_classified
            .filter { it[0] == "hybrid" }
            .map { it[1] }

        ASSEMBLE_SHORT(ch_auto_short, ref_file, assembly_utils_script)
        ASSEMBLE_LONG(ch_auto_long, ref_file, assembly_utils_script)
        ASSEMBLE_HYBRID(ch_auto_hybrid, ref_file, assembly_utils_script)

        ch_assembly = ASSEMBLE_SHORT.out.assembled
            .mix(ASSEMBLE_LONG.out.assembled)
            .mix(ASSEMBLE_HYBRID.out.assembled)
    }

    // Optional reference check: when reference is set, map reads to reference and report detection
    // Uses host-filtered reads if host_genome was provided, otherwise uses preprocessed reads
    if (params.reference) {
        REFERENCE_CHECK(ch_filtered_reads.join(ch_meta, remainder: true), ref_file, ref_check_plot_script)
    }

    KAIJU(ch_assembly.map { [it[0], it[1]] }, ch_kaiju_db)
    ch_kaiju_with_contigs = KAIJU.out.kaiju_done
        .join(ch_assembly.map { [it[0], it[1]] }, remainder: true)
    FILTER_VIRAL(ch_kaiju_with_contigs, extract_script)
    RENAME_CONTIGS(FILTER_VIRAL.out.viral)

    // Run CheckV and geNomad in parallel on viral contigs
    VALIDATE(RENAME_CONTIGS.out.renamed, ch_checkv_db)
    GENOMAD(RENAME_CONTIGS.out.renamed, ch_genomad_db)

    // Merge quality assessments: CheckV for phages, geNomad for RNA/eukaryotic viruses
    ch_merge_input = VALIDATE.out.validated
        .join(GENOMAD.out.genomad_done, remainder: true)
        .join(RENAME_CONTIGS.out.renamed, remainder: true)
        .map { t -> tuple(t[0], t[1], t[2], t[3], t[4]) }
        // Result: (sample_id, checkv_dir, genomad_out, viral_contigs, kaiju_dir)

    MERGE_QUALITY(ch_merge_input, merge_script)

    if (v_db) {
        ANNOTATE(RENAME_CONTIGS.out.renamed, ch_vog_db)
        
        // Organize genes by taxonomy
        ch_organize_input = ANNOTATE.out.annotated
            .join(RENAME_CONTIGS.out.renamed, remainder: true)
            .map { sample_id, annotation_dir, viral_contigs, kaiju_dir ->
                tuple(sample_id, annotation_dir, viral_contigs, kaiju_dir)
            }
        ORGANIZE_GENES(ch_organize_input, ch_vog_db)
    }

    // Quantification: join viral contigs with original filtered reads
    ch_quantify_input = RENAME_CONTIGS.out.renamed
        .join(ch_filtered_reads, remainder: true)
        .join(ch_meta, remainder: true)
        .map { t -> tuple(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[8]) }
        // (sample_id, viral_contigs, kaiju_dir, r1, r2, single, long, long_tech)
    QUANTIFY(ch_quantify_input)

    // Plotting: join quant, merged quality, and viral metadata
    ch_plot_input = QUANTIFY.out.quantified
        .join(MERGE_QUALITY.out.merged, remainder: true)
        .join(RENAME_CONTIGS.out.renamed, remainder: true)
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
