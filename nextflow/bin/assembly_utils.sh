#!/usr/bin/env bash
# Shared helper functions for Virall assembly processes.
# Sourced by ASSEMBLE_SHORT, ASSEMBLE_LONG, ASSEMBLE_HYBRID, and REF_ASSEMBLE.

# ── Coverage gates ──────────────────────────────────────────────────────────

check_spades_coverage() {
    local LOGFILE="$1"
    local SAMPLE_ID="$2"
    local SPADES_COV=""
    if [ -f "$LOGFILE" ]; then
        SPADES_COV=$(grep -a "Average coverage =" "$LOGFILE" | tail -1 | sed 's/.*Average coverage = //' | awk '{print $1}' | tr -d ',')
    fi
    # Fallback: metaviralSPAdes encodes coverage in contig headers
    if [ -z "$SPADES_COV" ] || [ "$SPADES_COV" = "0" ]; then
        if [ -f contigs ] && [ -s contigs ]; then
            SPADES_COV=$(grep "^>" contigs | grep -Eo 'cov_[0-9.]+' | sed 's/cov_//' | awk '{sum+=$1; n++} END {if(n>0) printf "%.1f", sum/n}')
        fi
    fi
    if [ -n "$SPADES_COV" ]; then
        # Use awk instead of bc for container portability
        local COV_OK=$(awk -v cov="$SPADES_COV" 'BEGIN {print (cov >= 10) ? 1 : 0}')
        if [ "$COV_OK" = "0" ]; then
            echo "========================================"
            echo "  WARNING: Low SPAdes coverage (${SPADES_COV}x)"
            echo "  Coverage <10x may yield fragmented assemblies"
            echo "  and unreliable abundance estimates."
            echo "========================================"
        fi
    fi
}

check_flye_coverage() {
    local LOGFILE="$1"
    local SAMPLE_ID="$2"
    if [ -f "$LOGFILE" ]; then
        local FLYE_COV=$(grep -a "Overlap-based coverage:" "$LOGFILE" | grep -Eo '[0-9]+' | tail -1)
        if [ -n "$FLYE_COV" ] && [ "$FLYE_COV" -gt 0 ] 2>/dev/null && [ "$FLYE_COV" -lt 10 ] 2>/dev/null; then
            echo "========================================"
            echo "  WARNING: Low Flye coverage (${FLYE_COV}x)"
            echo "  Coverage <10x may yield fragmented assemblies"
            echo "  and unreliable abundance estimates."
            echo "========================================"
        fi
    fi
}

# ── Output helpers ──────────────────────────────────────────────────────────

copy_spades_outputs() {
    local SPADES_DIR="$1"
    cp "${SPADES_DIR}/contigs.fasta" contigs 2>/dev/null || \
        cp "${SPADES_DIR}/transcripts.fasta" contigs 2>/dev/null || \
        touch contigs
    cp "${SPADES_DIR}/scaffolds.fasta" scaffolds 2>/dev/null || \
        cp "${SPADES_DIR}/hard_filtered_transcripts.fasta" scaffolds 2>/dev/null || \
        cp contigs scaffolds
}

# ── Reference-guided consensus ──────────────────────────────────────────────

run_ref_guided_consensus() {
    local LONG_READS="$1"
    local REF_FASTA="$2"
    local LONG_TECH="$3"
    local CPUS="$4"
    local OUT_FASTA="$5"
    local WORK_DIR="$6"

    mkdir -p "$WORK_DIR"

    local MINIMAP2_PRESET
    if [ "$LONG_TECH" = "pacbio" ]; then
        MINIMAP2_PRESET="map-pb"
    else
        MINIMAP2_PRESET="map-ont"
    fi

    minimap2 -t "$CPUS" -ax "$MINIMAP2_PRESET" "$REF_FASTA" \
        "$LONG_READS" 2>/dev/null | \
        samtools sort -@ "$CPUS" -o "${WORK_DIR}/ref_aligned.bam" -
    samtools index "${WORK_DIR}/ref_aligned.bam"

    if [ "$LONG_TECH" = "nanopore" ]; then
        medaka_consensus -i "$LONG_READS" \
            -d "$REF_FASTA" -o "${WORK_DIR}/ref_medaka" -t "$CPUS" \
            && cp "${WORK_DIR}/ref_medaka/consensus.fasta" "${WORK_DIR}/ref_consensus.fasta" \
            || samtools consensus "${WORK_DIR}/ref_aligned.bam" \
                 -o "${WORK_DIR}/ref_consensus.fasta" --show-ins no -a
    else
        samtools consensus "${WORK_DIR}/ref_aligned.bam" \
            -o "${WORK_DIR}/ref_consensus.fasta" --show-ins no -a
    fi

    cp "${WORK_DIR}/ref_consensus.fasta" "$OUT_FASTA" 2>/dev/null || true
}
