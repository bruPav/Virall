#!/usr/bin/env python3
"""
Merge CheckV and geNomad quality assessment results.

CheckV is optimal for bacteriophages; geNomad is better for RNA viruses
and eukaryotic DNA viruses.  This script merges both quality assessments
based on Kaiju taxonomy: phage contigs use CheckV results, others prefer
geNomad when available.

Used by the Nextflow MERGE_QUALITY process.
"""
import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd


# ── Phage detection ──────────────────────────────────────────────────────────

# Known NON-phage viral groups (checked first to avoid false positives from
# substring matches, e.g. Kitrinoviricota contains 'inovir').
NON_PHAGE_INDICATORS = [
    "kitrinoviricota",       # RNA viruses (Flaviviridae, etc.)
    "pisuviricota",          # RNA viruses (Picornavirales, etc.)
    "negarnaviricota",       # RNA viruses (negative-sense)
    "nucleocytoviricota",    # Giant DNA viruses (not phages)
    "artverviricota",        # Retroviruses
    "flaviviridae", "togaviridae", "coronaviridae", "picornaviridae",
    "rhabdoviridae", "paramyxoviridae", "orthomyxoviridae", "bunyavirales",
    "herpesviridae", "poxviridae", "adenoviridae", "papillomaviridae",
    "polyomaviridae", "retroviridae", "hepadnaviridae", "parvoviridae",
    "mimiviridae", "phycodnaviridae", "iridoviridae", "ascoviridae",
]

PHAGE_INDICATORS = [
    "caudovir",          # Tailed phages (Caudovirales, Caudoviricetes)
    "phage",             # Explicit "phage" in name
    "bacteriophage",     # Explicit bacteriophage
    "siphovir",          # Siphoviridae
    "myovir",            # Myoviridae
    "podovir",           # Podoviridae
    "autographivir",     # Autographiviridae
    "demerecvir",        # Demerecviridae
    "herellevir",        # Herelleviridae
    "inoviridae",        # Inoviridae (filamentous phages)
    "microviridae",      # Microviridae
    "leviviricetes",     # RNA phages
    "cystoviridae",      # dsRNA phages
    "tectiviridae",      # Tectiviridae
    "corticoviridae",    # Corticoviridae
    "plasmaviridae",     # Plasmaviridae
    "sphaerolipoviridae",  # Sphaerolipoviridae
    "uroviricota",       # Phage realm
    "peduoviridae",      # Peduoviridae
    "drexlerviridae",    # Drexlerviridae
    "ackermannviridae",  # Ackermannviridae
]


def is_phage(lineage_str: str) -> bool:
    """Determine if a contig is likely a bacteriophage based on Kaiju taxonomy."""
    if not lineage_str:
        return False
    lineage_lower = lineage_str.lower()
    if any(ind in lineage_lower for ind in NON_PHAGE_INDICATORS):
        return False
    return any(ind in lineage_lower for ind in PHAGE_INDICATORS)


# ── geNomad → quality tier mapping ──────────────────────────────────────────

# Expected hallmark density (hallmarks per kb) for a "complete" genome.
# Typical dsDNA viruses encode ~1-1.5 hallmark genes per kb.  We use 1.0 as
# a conservative baseline so the heuristic never overestimates completeness.
_HALLMARK_DENSITY_PER_KB = 1.0


def genomad_to_quality_tier(row: pd.Series, *,
                            contig_length: int = 0) -> tuple:
    """Convert geNomad metrics to a (quality_tier, completeness%) tuple.

    When *contig_length* is provided (in bp), the hallmark count is normalised
    by genome size so that a 200 kb virus with 3 hallmarks is not rated the
    same as a 5 kb phage with 3 hallmarks.
    """
    topology = str(row.get("topology", "")).lower()
    n_hallmarks = int(row.get("n_hallmarks", 0)) if pd.notna(row.get("n_hallmarks")) else 0
    virus_score = float(row.get("virus_score", 0)) if pd.notna(row.get("virus_score")) else 0

    if topology == "circular" and n_hallmarks >= 1:
        return "Complete", 100.0
    if topology == "circular":
        return "High-quality", 90.0

    # Length-normalised hallmark density (hallmarks per kb)
    length_kb = contig_length / 1000.0 if contig_length > 0 else 0.0
    if length_kb > 0 and n_hallmarks > 0:
        density = n_hallmarks / length_kb
        estimated_pct = min((density / _HALLMARK_DENSITY_PER_KB) * 100.0, 100.0)

        if estimated_pct >= 80.0:
            return "High-quality", round(estimated_pct, 1)
        if estimated_pct >= 40.0:
            return "Medium-quality", round(estimated_pct, 1)
        if estimated_pct >= 10.0:
            return "Low-quality", round(estimated_pct, 1)
        return "Low-quality", round(estimated_pct, 1)

    # Fallback when length is unavailable: use original fixed thresholds
    if n_hallmarks >= 3:
        return "High-quality", 80.0
    if n_hallmarks >= 1:
        return "Medium-quality", 50.0
    if virus_score >= 0.9:
        return "Low-quality", 30.0
    if virus_score >= 0.7:
        return "Low-quality", 20.0
    return "Not-determined", None


# ── I/O helpers ──────────────────────────────────────────────────────────────

def load_kaiju_taxonomy(kaiju_file: Path) -> dict:
    """Parse Kaiju results-with-names into {contig_id: lineage_str}."""
    taxonomy = {}
    if not kaiju_file.exists():
        return taxonomy
    try:
        with open(kaiju_file) as fh:
            for line in fh:
                if not line.startswith("C"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    contig_id = parts[1].strip()
                    lineage = parts[3] if len(parts) == 4 else ";".join(parts[3:10])
                    taxonomy[contig_id] = lineage.lower()
        print(f"Loaded taxonomy for {len(taxonomy)} contigs", file=sys.stderr)
    except Exception as exc:
        print(f"Warning: Could not read Kaiju file: {exc}", file=sys.stderr)
    return taxonomy


def load_tsv(path: Path, label: str) -> pd.DataFrame:
    """Read a TSV, returning an empty DataFrame on failure."""
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t")
        print(f"Loaded {len(df)} {label} results", file=sys.stderr)
        return df
    except Exception as exc:
        print(f"Warning: Could not read {label} file: {exc}", file=sys.stderr)
        return pd.DataFrame()


# ── Merge logic ──────────────────────────────────────────────────────────────

def _taxonomy_confidence(lineage: str) -> str:
    """Classify how trustworthy the Kaiju taxonomy call is.

    Returns one of:
      - "high"          – Kaiju returned a lineage that matched a known
                          phage or non-phage indicator.
      - "low_unclassified" – Kaiju returned no lineage at all; the routing
                          decision (phage vs non-phage) is essentially a guess.
      - "low_ambiguous"  – Kaiju returned a lineage but it matched neither
                          the phage nor the non-phage indicator lists.
    """
    if not lineage:
        return "low_unclassified"
    lineage_lower = lineage.lower()
    has_non_phage = any(ind in lineage_lower for ind in NON_PHAGE_INDICATORS)
    has_phage = any(ind in lineage_lower for ind in PHAGE_INDICATORS)
    if has_non_phage or has_phage:
        return "high"
    return "low_ambiguous"


def merge(checkv_df: pd.DataFrame,
          genomad_df: pd.DataFrame,
          kaiju_taxonomy: dict) -> list[dict]:
    """Return a list of merged-quality row dicts."""
    rows: list[dict] = []

    if not checkv_df.empty and "contig_id" in checkv_df.columns:
        for _, row in checkv_df.iterrows():
            cid = row["contig_id"]
            lineage = kaiju_taxonomy.get(cid, "")
            use_checkv = is_phage(lineage)
            confidence = _taxonomy_confidence(lineage)

            if use_checkv:
                rows.append({
                    "contig_id": cid,
                    "contig_length": row.get("contig_length", 0),
                    "quality_source": "checkv",
                    "taxonomy_confidence": confidence,
                    "checkv_quality": row.get("checkv_quality", "Not-determined"),
                    "completeness": row.get("completeness", None),
                    "completeness_method": row.get("completeness_method", ""),
                    "contamination": row.get("contamination", None),
                    "viral_genes": row.get("viral_genes", 0),
                    "host_genes": row.get("host_genes", 0),
                    "taxonomy_hint": "phage",
                    "provirus": row.get("provirus", "No"),
                })
            else:
                genomad_row = _match_genomad(genomad_df, cid)
                if genomad_row is not None:
                    entry = _combined_row(row, genomad_row)
                    entry["taxonomy_confidence"] = confidence
                    rows.append(entry)
                else:
                    rows.append({
                        "contig_id": cid,
                        "contig_length": row.get("contig_length", 0),
                        "quality_source": "checkv_fallback",
                        "taxonomy_confidence": confidence,
                        "checkv_quality": row.get("checkv_quality", "Not-determined"),
                        "completeness": row.get("completeness", None),
                        "completeness_method": row.get("completeness_method", ""),
                        "contamination": row.get("contamination", None),
                        "viral_genes": row.get("viral_genes", 0),
                        "host_genes": row.get("host_genes", 0),
                        "taxonomy_hint": "other",
                        "provirus": row.get("provirus", "No"),
                    })

    # Fallback: if CheckV was empty but geNomad produced results
    if not rows and not genomad_df.empty:
        print("Using geNomad results only (CheckV results empty)", file=sys.stderr)
        for _, grow in genomad_df.iterrows():
            contig_length = int(grow.get("length", 0) or 0)
            quality_tier, completeness = genomad_to_quality_tier(
                grow, contig_length=contig_length,
            )
            rows.append({
                "contig_id": grow.get("seq_name", ""),
                "contig_length": contig_length,
                "quality_source": "genomad",
                "taxonomy_confidence": "low_unclassified",
                "checkv_quality": quality_tier,
                "completeness": completeness,
                "completeness_method": "genomad_hallmarks",
                "contamination": None,
                "viral_genes": grow.get("n_genes", 0),
                "host_genes": 0,
                "taxonomy_hint": "unknown",
                "provirus": "No",
                "genomad_score": grow.get("virus_score", None),
                "genomad_topology": grow.get("topology", ""),
                "genomad_hallmarks": grow.get("n_hallmarks", 0),
                "genomad_taxonomy": grow.get("taxonomy", ""),
            })

    return rows


def _match_genomad(genomad_df: pd.DataFrame, contig_id: str):
    """Return the first matching geNomad row for *contig_id*, or None."""
    if genomad_df.empty or "seq_name" not in genomad_df.columns:
        return None
    match = genomad_df[genomad_df["seq_name"] == contig_id]
    return match.iloc[0] if not match.empty else None


def _combined_row(checkv_row, genomad_row) -> dict:
    """Merge a CheckV row with a matching geNomad row for non-phage contigs.

    For non-phage viruses, geNomad's length-normalised hallmark assessment is
    preferred over CheckV (whose reference DB is phage-centric).  CheckV's
    completeness is kept as a secondary column for reference, and its
    contamination estimate is always retained (geNomad does not provide one).
    """
    checkv_completeness = checkv_row.get("completeness", None)
    has_checkv_compl = pd.notna(checkv_completeness) and checkv_completeness != ""

    contig_length = int(genomad_row.get("length", 0) or 0)
    genomad_quality, genomad_completeness = genomad_to_quality_tier(
        genomad_row, contig_length=contig_length,
    )

    if genomad_completeness is not None:
        final_completeness = genomad_completeness
        final_quality = genomad_quality
        final_method = "genomad_hallmarks"
        source = "checkv+genomad"
    elif has_checkv_compl:
        final_completeness = checkv_completeness
        final_quality = checkv_row.get("checkv_quality", "Not-determined")
        final_method = checkv_row.get("completeness_method", "")
        source = "checkv+genomad"
    else:
        final_completeness = None
        final_quality = "Not-determined"
        final_method = ""
        source = "genomad"

    return {
        "contig_id": checkv_row["contig_id"],
        "contig_length": checkv_row.get("contig_length", 0),
        "quality_source": source,
        "checkv_quality": final_quality,
        "completeness": final_completeness,
        "completeness_method": final_method,
        "checkv_completeness": checkv_completeness if has_checkv_compl else None,
        "contamination": checkv_row.get("contamination", None),
        "viral_genes": max(
            checkv_row.get("viral_genes", 0),
            int(genomad_row.get("n_genes", 0) or 0),
        ),
        "host_genes": checkv_row.get("host_genes", 0),
        "taxonomy_hint": "rna_or_eukaryotic",
        "provirus": checkv_row.get("provirus", "No"),
        "genomad_score": genomad_row.get("virus_score", None),
        "genomad_topology": genomad_row.get("topology", ""),
        "genomad_hallmarks": genomad_row.get("n_hallmarks", 0),
        "genomad_taxonomy": genomad_row.get("taxonomy", ""),
    }


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Merge CheckV + geNomad quality assessments"
    )
    parser.add_argument("--checkv-dir", type=Path, required=True,
                        help="Directory containing CheckV quality_summary.tsv")
    parser.add_argument("--genomad-dir", type=Path, required=True,
                        help="Directory containing geNomad virus_summary.tsv")
    parser.add_argument("--kaiju-dir", type=Path, required=True,
                        help="Directory containing kaiju_results_with_names.tsv")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="Output directory for merged results")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    checkv_file = args.checkv_dir / "quality_summary.tsv"
    genomad_file = args.genomad_dir / "virus_summary.tsv"
    kaiju_file = args.kaiju_dir / "kaiju_results_with_names.tsv"
    output_file = args.out_dir / "quality_summary.tsv"

    checkv_df = load_tsv(checkv_file, "CheckV")
    if not checkv_df.empty:
        checkv_df["quality_source"] = "checkv"
    genomad_df = load_tsv(genomad_file, "geNomad")
    kaiju_taxonomy = load_kaiju_taxonomy(kaiju_file)

    merged_rows = merge(checkv_df, genomad_df, kaiju_taxonomy)

    if merged_rows:
        merged_df = pd.DataFrame(merged_rows)
        merged_df.to_csv(output_file, sep="\t", index=False)
        print(f"Wrote {len(merged_df)} rows to merged quality summary",
              file=sys.stderr)
        print(f"Quality sources: {merged_df['quality_source'].value_counts().to_dict()}",
              file=sys.stderr)
        print(f"Quality tiers: {merged_df['checkv_quality'].value_counts().to_dict()}",
              file=sys.stderr)
        if "taxonomy_confidence" in merged_df.columns:
            print(f"Taxonomy confidence: {merged_df['taxonomy_confidence'].value_counts().to_dict()}",
                  file=sys.stderr)
    else:
        empty_cols = [
            "contig_id", "contig_length", "quality_source",
            "taxonomy_confidence", "checkv_quality", "completeness",
            "completeness_method", "contamination",
        ]
        pd.DataFrame(columns=empty_cols).to_csv(output_file, sep="\t",
                                                  index=False)
        print("No quality data to merge", file=sys.stderr)

    # Copy originals for reference
    if checkv_file.exists():
        shutil.copy(checkv_file, args.out_dir / "checkv_quality_summary.tsv")
    if genomad_file.exists() and genomad_file.stat().st_size > 0:
        shutil.copy(genomad_file, args.out_dir / "genomad_virus_summary.tsv")


if __name__ == "__main__":
    main()
