#!/usr/bin/env python3
"""Shared Kaiju taxonomy parser used by all Virall pipeline Python helpers."""

from pathlib import Path


def parse_kaiju(tsv_path: Path) -> dict:
    """Parse kaiju_results_with_names.tsv into {contig_id: {...}}.

    Handles both 4-column (standard) and 10+column (addTaxonNames expanded)
    formats.  Returns an empty dict on missing/empty file or zero classified
    lines — callers already handle this gracefully.
    """
    if not tsv_path.exists():
        return {}

    taxonomy: dict = {}

    with open(tsv_path) as f:
        for line in f:
            line = line.strip()
            if not line or not line.startswith("C"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue

            contig_id = parts[1].strip()
            taxon_id = parts[2].strip() if len(parts) > 2 else None

            if len(parts) >= 10:
                lineage = ";".join(parts[3:10])
                parsed = [p.strip() for p in parts[3:10] if p.strip()]
            else:
                lineage = parts[3].strip()
                parsed = [p.strip() for p in lineage.split(";") if p.strip()]

            filtered = [p for p in parsed if p and p != "NA" and p != "Viruses"]
            lowest = filtered[-1] if filtered else contig_id

            taxonomy[contig_id] = {
                "taxon_id": taxon_id,
                "lineage": lineage,
                "parsed_lineage": filtered,
                "lowest_taxon": lowest,
            }

    return taxonomy
