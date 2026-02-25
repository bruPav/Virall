#!/usr/bin/env python3
"""
Generate plots from Nextflow pipeline outputs.
Self-contained plotting - no dependency on virall Python package.
Builds contig_abundance.tsv and kaiju_summary.tsv from pipeline files,
then generates plots using matplotlib/plotly.
"""
import argparse
import shutil
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("run_plots.py requires pandas: pip install pandas")

def main():
    parser = argparse.ArgumentParser(description="Generate plots from Virall Nextflow outputs")
    parser.add_argument("--quant-dir", type=Path, required=True, help="06_quantification directory (depth.txt)")
    parser.add_argument("--viral-contigs", type=Path, required=True, help="Viral contigs FASTA")
    parser.add_argument("--kaiju-dir", type=Path, required=True, help="03_classifications directory (kaiju_results_with_names.tsv)")
    parser.add_argument("--checkv-dir", type=Path, required=True, help="04_quality_assessment directory")
    parser.add_argument("--min-breadth", type=float, default=0.10, help="Minimum contig coverage breadth (fraction 0-1) used for abundance estimation")
    parser.add_argument("--out-dir", type=Path, required=True, help="07_plots output directory")
    args = parser.parse_args()

    work_dir = args.out_dir / ".work"
    work_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    depth_primary_candidates = sorted(args.quant_dir.glob("depth_primary_mapq*.txt"))
    depth_file = depth_primary_candidates[0] if depth_primary_candidates else (args.quant_dir / "depth.txt")
    print(f"[run_plots] using depth file: {depth_file}", file=sys.stderr)
    kaiju_file = args.kaiju_dir / "kaiju_results_with_names.tsv"

    # 1. Build contig_abundance.tsv from depth.txt + viral_contigs + kaiju
    contig_lengths = {}
    try:
        from Bio import SeqIO
        for rec in SeqIO.parse(args.viral_contigs, "fasta"):
            contig_lengths[rec.id.split()[0]] = len(rec.seq)
    except ImportError:
        pass
    if not contig_lengths and depth_file.exists():
        # Infer lengths from max position in depth file
        max_pos = {}
        with open(depth_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    cid, pos, depth = parts[0], int(parts[1]), float(parts[2])
                    max_pos[cid] = max(max_pos.get(cid, 0), pos)
        contig_lengths = {c: p + 1 for c, p in max_pos.items()}

    # Parse Kaiju addTaxonNames output. Format: C, contig_id, taxon_id, [lineage or rank columns]
    # - 4 cols: 4th = full lineage (semicolon-separated) -> species = last part
    # - 10 cols: 4th..10th = superkingdom,phylum,...,species -> species = last column
    species_by_contig = {}
    kaiju_lineage_by_contig = {}
    if kaiju_file.exists():
        with open(kaiju_file) as f:
            for line in f:
                if line.startswith("C"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        cid = parts[1].strip()
                        if len(parts) >= 10:
                            lineage = ";".join(p.strip() for p in parts[3:10])
                            name = (parts[9].strip() if len(parts) > 9 else "") or "Unknown"
                        else:
                            lineage = parts[3].strip()
                            lineage_parts = [p.strip() for p in lineage.split(";") if p.strip()]
                            name = lineage_parts[-1] if lineage_parts else "Unknown"
                        species_by_contig[cid] = name or "Unknown"
                        kaiju_lineage_by_contig[cid] = lineage or name or "Unknown"

    def species_for_contig(depth_cid):
        cid = depth_cid.strip() if isinstance(depth_cid, str) else str(depth_cid)
        if cid in species_by_contig:
            return species_by_contig[cid]
        cid_first = cid.split()[0]
        if cid_first in species_by_contig:
            return species_by_contig[cid_first]
        for kaiju_cid, sp in species_by_contig.items():
            if kaiju_cid.startswith(cid) or cid.startswith(kaiju_cid):
                return sp
        return "Unknown"

    # Debug: report Kaiju parsing and ID matching (stderr)
    print(f"[run_plots] kaiju_file={kaiju_file} exists={kaiju_file.exists()}", file=sys.stderr)
    if kaiju_file.exists():
        with open(kaiju_file) as f:
            for L in f:
                if L.startswith("C"):
                    parts = L.strip().split("\t")
                    if len(parts) >= 4:
                        print(f"[run_plots] first C line: ncols={len(parts)} cid={parts[1][:50]!r} col3={parts[3][:80]!r}", file=sys.stderr)
                    break
        print(f"[run_plots] species_by_contig entries: {len(species_by_contig)}", file=sys.stderr)
        if species_by_contig:
            sample = list(species_by_contig.items())[:3]
            print(f"[run_plots] sample Kaiju IDs->species: {sample}", file=sys.stderr)

    abundance_file = work_dir / "contig_abundance.tsv"
    species_abundance_file = work_dir / "species_abundance.tsv"
    if depth_file.exists() and contig_lengths:
        depth_df = pd.read_csv(depth_file, sep="\t", header=None, names=["contig_id", "pos", "depth"])
        depth_ids = depth_df["contig_id"].unique().tolist()
        print(f"[run_plots] depth contig_ids (first 5): {depth_ids[:5]}", file=sys.stderr)
        total_cov = depth_df.groupby("contig_id").agg(
            total_coverage=("depth", "sum"),
            mapped_positions=("depth", lambda x: (x > 0).sum()),
        ).reset_index()
        total_cov["contig_length"] = total_cov["contig_id"].map(contig_lengths).fillna(0).astype(int)
        total_cov["contig_length_kb"] = total_cov["contig_length"].replace(0, 1) / 1000.0
        total_cov["mean_coverage"] = (total_cov["total_coverage"] / total_cov["contig_length"].replace(0, 1)).round(2)
        total_cov["max_coverage"] = depth_df.groupby("contig_id")["depth"].max().reindex(total_cov["contig_id"]).values
        total_cov["coverage_breadth"] = (total_cov["mapped_positions"] / total_cov["contig_length"].replace(0, 1)).round(4)
        min_breadth = min(max(float(args.min_breadth), 0.0), 1.0)
        before_filter = len(total_cov)
        total_cov = total_cov[total_cov["coverage_breadth"] >= min_breadth].copy()
        print(f"[run_plots] breadth filter >= {min_breadth:.4f}: kept {len(total_cov)}/{before_filter} contigs", file=sys.stderr)

        sum_all = total_cov["total_coverage"].sum()
        total_cov["relative_abundance_coverage"] = (total_cov["total_coverage"] / sum_all).round(6) if sum_all > 0 else 0.0
        total_cov["rpk"] = total_cov["total_coverage"] / total_cov["contig_length_kb"].replace(0, 1)
        rpk_sum = total_cov["rpk"].sum()
        total_cov["tpm"] = (total_cov["rpk"] / rpk_sum * 1_000_000).round(3) if rpk_sum > 0 else 0.0
        total_cov["relative_abundance"] = (total_cov["tpm"] / 1_000_000).round(6)
        total_cov["species"] = total_cov["contig_id"].map(species_for_contig)
        matched = (total_cov["species"] != "Unknown").sum()
        print(f"[run_plots] abundance rows: {len(total_cov)} species matched: {matched} unknown: {(total_cov['species'] == 'Unknown').sum()}", file=sys.stderr)
        total_cov = total_cov.sort_values("relative_abundance", ascending=False)
        total_cov.to_csv(abundance_file, sep="\t", index=False)
        shutil.copy2(abundance_file, args.out_dir / "contig_abundance.tsv")

        species_cov = total_cov.groupby("species", dropna=False).agg(
            contig_count=("contig_id", "count"),
            total_coverage=("total_coverage", "sum"),
            total_contig_length=("contig_length", "sum"),
            rpk=("rpk", "sum"),
        ).reset_index()
        species_rpk_sum = species_cov["rpk"].sum()
        species_cov["tpm"] = (species_cov["rpk"] / species_rpk_sum * 1_000_000).round(3) if species_rpk_sum > 0 else 0.0
        species_cov["relative_abundance"] = (species_cov["tpm"] / 1_000_000).round(6)
        species_cov = species_cov.sort_values("tpm", ascending=False)
        species_cov.to_csv(species_abundance_file, sep="\t", index=False)
        shutil.copy2(species_abundance_file, args.out_dir / "species_abundance.tsv")

    # 2. Build kaiju_summary.tsv for sunburst (contig_id, taxon_id, taxon_name, classification, lineage)
    #    Only include contigs present in viral_contigs.fasta (i.e. those that passed FILTER_VIRAL)
    kaiju_summary_file = work_dir / "kaiju_summary.tsv"
    print(f"[run_plots] Building kaiju_summary: kaiju_file={kaiju_file}, exists={kaiju_file.exists()}, lineages={len(kaiju_lineage_by_contig)}", file=sys.stderr)

    # Build set of viral contig IDs from the filtered FASTA
    viral_contig_ids = set()
    if args.viral_contigs.exists():
        with open(args.viral_contigs) as fasta_f:
            for line in fasta_f:
                if line.startswith(">"):
                    viral_contig_ids.add(line[1:].strip().split()[0])
    print(f"[run_plots] viral_contigs FASTA IDs: {len(viral_contig_ids)}", file=sys.stderr)

    if kaiju_file.exists() and kaiju_lineage_by_contig:
        rows = []
        skipped = 0
        with open(kaiju_file) as f:
            for line in f:
                if line.startswith("C"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        cid, taxon_id = parts[1], parts[2]
                        # Only include contigs that passed FILTER_VIRAL
                        if viral_contig_ids and cid.strip() not in viral_contig_ids:
                            skipped += 1
                            continue
                        if len(parts) >= 10:
                            lineage = ";".join(p.strip() for p in parts[3:10])
                            name = (parts[9].strip() if len(parts) > 9 else "") or "Unknown"
                        else:
                            lineage = parts[3].strip()
                            lineage_parts = [p.strip() for p in lineage.split(";") if p.strip()]
                            name = lineage_parts[-1] if lineage_parts else "Unknown"
                        rows.append({"contig_id": cid, "taxon_id": taxon_id, "taxon_name": name or "Unknown", "classification": name or "Unknown", "lineage": lineage or name or "Unknown"})
        print(f"[run_plots] kaiju_summary: found {len(rows)} classified contigs ({skipped} skipped – not in viral_contigs.fasta)", file=sys.stderr)
        if rows:
            pd.DataFrame(rows).to_csv(kaiju_summary_file, sep="\t", index=False)
            print(f"[run_plots] Wrote {kaiju_summary_file}", file=sys.stderr)
        else:
            print("[run_plots] No classified contigs found in kaiju results", file=sys.stderr)
    else:
        if not kaiju_file.exists():
            print(f"[run_plots] kaiju_file not found: {kaiju_file}", file=sys.stderr)
        if not kaiju_lineage_by_contig:
            print("[run_plots] No lineages parsed from kaiju results (all contigs unclassified?)", file=sys.stderr)

    # 3. Find quality_summary.tsv (CheckV writes in output dir; may be quality_summary.tsv or in subdir)
    quality_file = args.checkv_dir / "quality_summary.tsv"
    if not quality_file.exists():
        for sub in ["contigs", "viral_contigs", ""]:
            cand = args.checkv_dir / sub / "quality_summary.tsv" if sub else quality_file
            if cand.exists():
                quality_file = cand
                break

    # 4. Generate plots using standalone matplotlib/seaborn implementation
    # (No dependency on virall package - all logic is self-contained)
    # Matches the style of virall/core/plotter.py
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
        HAS_SEABORN = True
    except ImportError:
        HAS_SEABORN = False
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not available; skipping plots.", file=sys.stderr)
            print(f"Plots written to {args.out_dir}", file=sys.stderr)
            return
    
    try:
        import plotly.express as px
        HAS_PLOTLY = True
    except ImportError:
        HAS_PLOTLY = False
    
    # Set plot style to match plotter.py
    if HAS_SEABORN:
        sns.set_theme(style="whitegrid")
    plt.rcParams.update({'figure.autolayout': True})
    
    # =========================================================================
    # Abundance Plot (matches plotter.py style)
    # =========================================================================
    if abundance_file.exists():
        df = pd.read_csv(abundance_file, sep="\t")
        if species_abundance_file.exists():
            sdf = pd.read_csv(species_abundance_file, sep="\t")
        else:
            sdf = pd.DataFrame()
        if not sdf.empty and "tpm" in sdf.columns and "species" in sdf.columns:
            top_species = sdf.sort_values("tpm", ascending=False).head(20).set_index("species")["tpm"]
            abundance_title = "Top 20 Most Abundant Viral Species (TPM-like)"
            abundance_xlabel = "TPM (length-normalized)"
        elif not df.empty and "relative_abundance" in df.columns:
            # Fallback: derive from contig-level table
            if "species" in df.columns:
                species_abundance = df.groupby("species")["relative_abundance"].sum().sort_values(ascending=False)
                top_species = species_abundance.head(20)
            else:
                top = df.head(20)
                top_species = pd.Series(top["relative_abundance"].values, index=top.get("contig_id", range(len(top))))
            abundance_title = "Top 20 Most Abundant Viral Species"
            abundance_xlabel = "Relative Abundance"
        else:
            top_species = pd.Series(dtype=float)
            
        if not top_species.empty:
            plt.figure(figsize=(10, 8))
            if HAS_SEABORN:
                sns.barplot(x=top_species.values, y=top_species.index, 
                           hue=top_species.index, legend=False, palette="viridis")
            else:
                plt.barh(range(len(top_species)), top_species.values, color="steelblue")
                plt.yticks(range(len(top_species)), top_species.index.tolist(), fontsize=8)
            
            plt.title(abundance_title, fontsize=16)
            plt.xlabel(abundance_xlabel, fontsize=12)
            plt.ylabel("Species", fontsize=12)
            plt.tight_layout()
            plt.savefig(args.out_dir / "viral_abundance.png", dpi=300, bbox_inches="tight")
            plt.close()
            print(f"Generated viral_abundance.png", file=sys.stderr)
        
        # Contig Length Distribution (matches plotter.py style with log scale, KDE, N50)
        if not df.empty and "contig_length" in df.columns:
            valid = df[df["contig_length"] > 0].copy()
            if not valid.empty:
                plt.figure(figsize=(10, 6))
                
                # Determine if log scale and KDE are appropriate
                use_log = True
                use_kde = len(valid) >= 5
                
                # If range is small, log scale might not be useful
                if valid["contig_length"].max() / valid["contig_length"].min() < 10:
                    use_log = False
                
                if HAS_SEABORN:
                    try:
                        sns.histplot(data=valid, x="contig_length", log_scale=use_log, 
                                    kde=use_kde, color="teal")
                    except Exception:
                        # Fallback if seaborn fails
                        plt.hist(valid["contig_length"], bins=30, color="teal", alpha=0.7)
                        if use_log:
                            plt.xscale("log")
                else:
                    plt.hist(valid["contig_length"], bins=30, color="teal", alpha=0.7)
                    if use_log:
                        plt.xscale("log")
                
                plt.title("Distribution of Viral Contig Lengths", fontsize=16)
                plt.xlabel(f"Contig Length (bp{', log scale' if use_log else ''})", fontsize=12)
                plt.ylabel("Count", fontsize=12)
                
                # Add N50 line (like plotter.py)
                lengths = sorted(valid["contig_length"].tolist(), reverse=True)
                total_len = sum(lengths)
                cum_len = 0
                n50 = 0
                for l in lengths:
                    cum_len += l
                    if cum_len >= total_len / 2:
                        n50 = l
                        break
                
                if n50 > 0:
                    plt.axvline(n50, color="red", linestyle="--", label=f"N50: {n50} bp")
                    plt.legend()
                
                plt.savefig(args.out_dir / "contig_length_distribution.png", dpi=300, bbox_inches="tight")
                plt.close()
                print(f"Generated contig_length_distribution.png", file=sys.stderr)
    
    # =========================================================================
    # Taxonomy Sunburst (requires plotly) - matches plotter.py style
    # =========================================================================
    print(f"[run_plots] Sunburst: HAS_PLOTLY={HAS_PLOTLY}, kaiju_summary_file={kaiju_summary_file}, exists={kaiju_summary_file.exists()}", file=sys.stderr)
    
    if not HAS_PLOTLY:
        print("[run_plots] Skipping sunburst: plotly not available", file=sys.stderr)
    elif not kaiju_summary_file.exists():
        print(f"[run_plots] Skipping sunburst: {kaiju_summary_file} not found", file=sys.stderr)
        print(f"[run_plots] kaiju_lineage_by_contig has {len(kaiju_lineage_by_contig)} entries", file=sys.stderr)
    
    if HAS_PLOTLY and kaiju_summary_file.exists():
        try:
            kdf = pd.read_csv(kaiju_summary_file, sep="\t")
            print(f"[run_plots] Sunburst: loaded {len(kdf)} rows from kaiju_summary", file=sys.stderr)
            if not kdf.empty and "lineage" in kdf.columns:
                data = []
                for lineage in kdf["lineage"].dropna():
                    # Split by semicolon and filter out ONLY empty parts (keep "NA" as valid taxonomy)
                    # Lineages often end with ";" which creates empty strings
                    parts = [p.strip() for p in lineage.split(";") if p.strip()]
                    if parts:
                        # Add 'Viruses' as root if not present (like plotter.py)
                        if parts[0] != "Viruses":
                            parts.insert(0, "Viruses")
                        data.append(parts)
                
                print(f"[run_plots] Sunburst: parsed {len(data)} lineages with data", file=sys.stderr)
                
                if not data:
                    print("[run_plots] Skipping sunburst: no valid lineage data found", file=sys.stderr)
                
                if data:
                    max_depth = max(len(d) for d in data)
                    cols = [f"Level_{i}" for i in range(max_depth)]
                    
                    # Pad shorter lineages with None (plotly handles None/NaN fine)
                    padded_data = [d + [None] * (max_depth - len(d)) for d in data]
                    df_lineage = pd.DataFrame(padded_data, columns=cols)
                    df_lineage["count"] = 1
                    
                    # Only replace empty strings with None - keep "NA" as valid taxonomy level
                    # Empty strings cause plotly's "non-leaf" error
                    for col in cols:
                        df_lineage[col] = df_lineage[col].replace("", None)
                    
                    # Find valid columns: plotly sunburst can't handle None in paths,
                    # so we only include columns where ALL rows have data (no nulls).
                    # We stop at the first column that has any nulls.
                    valid_cols = []
                    for col in cols:
                        if df_lineage[col].isna().sum() == 0:
                            valid_cols.append(col)
                        else:
                            break  # Stop at first column with nulls
                    
                    print(f"[run_plots] Sunburst: max_depth={max_depth}, valid_cols={len(valid_cols)}", file=sys.stderr)
                    
                    if valid_cols and len(df_lineage) > 0:
                        # Find best column for color
                        color_col = valid_cols[0]
                        start_idx = 1 if len(valid_cols) > 1 else 0
                        
                        for col in valid_cols[start_idx:]:
                            n_unique = df_lineage[col].nunique()
                            if n_unique > 1:
                                color_col = col
                                break
                        
                        fig = px.sunburst(
                            df_lineage, 
                            path=valid_cols, 
                            values="count",
                            color=color_col,
                            title="Viral Taxonomy Sunburst", 
                            width=800, 
                            height=800
                        )
                        fig.update_traces(textinfo="label+percent entry")
                        fig.write_html(str(args.out_dir / "viral_taxonomy_sunburst.html"))
                        print(f"Generated viral_taxonomy_sunburst.html", file=sys.stderr)
                        
                        # Try to save static image (like plotter.py)
                        try:
                            fig.write_image(str(args.out_dir / "viral_taxonomy_sunburst.png"))
                            print(f"Generated viral_taxonomy_sunburst.png", file=sys.stderr)
                        except Exception:
                            pass  # kaleido might not be available
        except Exception as e:
            print(f"Could not generate sunburst plot: {e}", file=sys.stderr)
    
    # =========================================================================
    # Quality Plots (merged CheckV + geNomad data) - matches plotter.py style
    # =========================================================================
    if quality_file.exists():
        qdf = pd.read_csv(quality_file, sep="\t")
        if not qdf.empty and "checkv_quality" in qdf.columns:
            # Determine data sources available
            has_quality_source = "quality_source" in qdf.columns
            has_checkv = False
            has_genomad = False
            
            if has_quality_source:
                sources = qdf["quality_source"].unique()
                has_checkv = any(s in ["checkv", "checkv_fallback"] for s in sources if pd.notna(s))
                has_genomad = "genomad" in [s for s in sources if pd.notna(s)]
            
            # =================================================================
            # 1. Quality Distribution (both tools, or one if only one ran)
            # =================================================================
            plt.figure(figsize=(10, 6))
            if has_quality_source and (has_checkv or has_genomad):
                grouped = qdf.groupby(["checkv_quality", "quality_source"]).size().unstack(fill_value=0)
                if not grouped.empty:
                    ax = grouped.plot(kind="bar", stacked=True, figsize=(10, 6), 
                                     colormap="viridis", alpha=0.8)
                    if has_checkv and has_genomad:
                        title = "Distribution of Viral Genome Quality\n(CheckV for phages, geNomad for RNA/eukaryotic viruses)"
                    elif has_genomad:
                        title = "Distribution of Viral Genome Quality\n(geNomad assessment)"
                    else:
                        title = "Distribution of Viral Genome Quality\n(CheckV assessment)"
                    plt.title(title, fontsize=14)
                    plt.xlabel("Quality Tier", fontsize=12)
                    plt.ylabel("Number of Genomes", fontsize=12)
                    plt.legend(title="Quality Source", loc="upper right")
                    plt.xticks(rotation=45, ha="right")
                    plt.tight_layout()
            else:
                # Simple bar chart (no quality_source column) - use seaborn like plotter.py
                quality_counts = qdf["checkv_quality"].value_counts()
                if HAS_SEABORN:
                    sns.barplot(x=quality_counts.index, y=quality_counts.values,
                               hue=quality_counts.index, legend=False, palette="viridis")
                else:
                    plt.bar(quality_counts.index, quality_counts.values, color="teal", alpha=0.7)
                plt.title("Distribution of Viral Genome Quality", fontsize=16)
                plt.xlabel("Quality Tier", fontsize=12)
                plt.ylabel("Number of Genomes", fontsize=12)
                plt.xticks(rotation=45)
                plt.tight_layout()
            plt.savefig(args.out_dir / "genome_quality_distribution.png", dpi=300, bbox_inches="tight")
            plt.close()
            print(f"Generated genome_quality_distribution.png", file=sys.stderr)
            
            # =================================================================
            # 2. Scatter Plot: Completeness vs Contamination (CheckV only)
            # =================================================================
            if "completeness" in qdf.columns and "contamination" in qdf.columns:
                # Filter to CheckV results only (they have contamination data)
                if has_quality_source:
                    checkv_df = qdf[qdf["quality_source"].isin(["checkv", "checkv_fallback"])].copy()
                else:
                    checkv_df = qdf.copy()
                
                scatter_df = checkv_df[checkv_df["completeness"].notna() & checkv_df["contamination"].notna()].copy()
                if not scatter_df.empty:
                    plt.figure(figsize=(10, 8))
                    
                    if HAS_SEABORN:
                        sns.scatterplot(
                            data=scatter_df,
                            x="completeness",
                            y="contamination",
                            hue="checkv_quality",
                            style="checkv_quality",
                            s=100,
                            alpha=0.7,
                            palette="viridis"
                        )
                    else:
                        quality_tiers = scatter_df["checkv_quality"].unique()
                        colors = plt.cm.viridis([i / max(len(quality_tiers), 1) for i in range(len(quality_tiers))])
                        for i, tier in enumerate(quality_tiers):
                            tier_df = scatter_df[scatter_df["checkv_quality"] == tier]
                            plt.scatter(tier_df["completeness"], tier_df["contamination"], 
                                       label=tier, alpha=0.7, s=100, c=[colors[i]])
                    
                    plt.title("Viral Genome Quality: Completeness vs. Contamination\n(CheckV assessment - phage genomes)", fontsize=14)
                    plt.xlabel("Completeness (%)", fontsize=12)
                    plt.ylabel("Contamination (%)", fontsize=12)
                    plt.grid(True, linestyle="--", alpha=0.3)
                    plt.xlim(0, 105)
                    plt.ylim(-5, 105)
                    
                    # Add quality zone lines (like plotter.py)
                    plt.axvline(x=90, color="green", linestyle=":", alpha=0.5, label="High Quality (>90%)")
                    plt.axvline(x=50, color="orange", linestyle=":", alpha=0.5, label="Medium Quality (>50%)")
                    
                    plt.tight_layout()
                    plt.savefig(args.out_dir / "genome_quality_scatter.png", dpi=300, bbox_inches="tight")
                    plt.close()
                    print(f"Generated genome_quality_scatter.png ({len(scatter_df)} CheckV contigs)", file=sys.stderr)
            
            # =================================================================
            # 3. Completeness Distribution (both tools, or one if only one has data)
            # =================================================================
            if "completeness" in qdf.columns:
                comp_df = qdf[qdf["completeness"].notna()].copy()
                if not comp_df.empty:
                    plt.figure(figsize=(10, 6))
                    
                    if has_quality_source:
                        sources_with_data = [s for s in comp_df["quality_source"].unique() if pd.notna(s)]
                        n_sources = len(sources_with_data)
                        
                        if n_sources > 1:
                            # Multiple sources - box plot comparison (like plotter.py)
                            if HAS_SEABORN:
                                sns.boxplot(
                                    data=comp_df,
                                    x="quality_source",
                                    y="completeness",
                                    hue="quality_source",
                                    palette="viridis",
                                    legend=False
                                )
                            else:
                                data_by_source = [comp_df[comp_df["quality_source"] == s]["completeness"].values 
                                                 for s in sources_with_data]
                                plt.boxplot(data_by_source, labels=sources_with_data)
                            plt.title("Completeness Distribution by Assessment Tool", fontsize=16)
                            plt.xlabel("Quality Assessment Tool", fontsize=12)
                        elif n_sources == 1:
                            # Single source - simple box plot
                            source_name = sources_with_data[0]
                            if HAS_SEABORN:
                                sns.boxplot(y=comp_df["completeness"], color="steelblue")
                            else:
                                plt.boxplot([comp_df["completeness"].values])
                            plt.title(f"Completeness Distribution ({source_name})", fontsize=16)
                            plt.xlabel("")
                        else:
                            if HAS_SEABORN:
                                sns.boxplot(y=comp_df["completeness"], color="steelblue")
                            else:
                                plt.boxplot([comp_df["completeness"].values])
                            plt.title("Completeness Distribution", fontsize=16)
                            plt.xlabel("")
                    else:
                        # No quality_source column - simple box plot
                        if HAS_SEABORN:
                            sns.boxplot(y=comp_df["completeness"], color="steelblue")
                        else:
                            plt.boxplot([comp_df["completeness"].values])
                        plt.title("Completeness Distribution", fontsize=16)
                        plt.xlabel("")
                    
                    plt.ylabel("Completeness (%)", fontsize=12)
                    plt.ylim(0, 105)
                    plt.tight_layout()
                    plt.savefig(args.out_dir / "completeness_by_source.png", dpi=300, bbox_inches="tight")
                    plt.close()
                    print(f"Generated completeness_by_source.png ({len(comp_df)} contigs)", file=sys.stderr)
    
    # =========================================================================
    # Coverage Plots for High-Quality/Complete Contigs (matches plotter.py style)
    # =========================================================================
    if quality_file.exists() and depth_file.exists():
        try:
            qdf = pd.read_csv(quality_file, sep="\t")
            if not qdf.empty and "checkv_quality" in qdf.columns:
                # Find Complete or High-quality contigs
                hq_contigs = qdf[qdf["checkv_quality"].isin(["Complete", "High-quality"])]
                if not hq_contigs.empty:
                    # Sort by length (longest first) if available
                    if "contig_length" in hq_contigs.columns:
                        hq_contigs = hq_contigs.sort_values("contig_length", ascending=False)
                    
                    target_contigs = hq_contigs["contig_id"].head(10).tolist()
                    
                    # Read depth data for these contigs
                    coverage_dir = args.out_dir / "coverage_plots"
                    coverage_dir.mkdir(exist_ok=True)
                    
                    depth_df = pd.read_csv(depth_file, sep="\t", header=None, names=["contig_id", "pos", "depth"])
                    
                    plots_generated = 0
                    for contig_id in target_contigs:
                        contig_depth = depth_df[depth_df["contig_id"] == contig_id].sort_values("pos")
                        if not contig_depth.empty:
                            plt.figure(figsize=(12, 4))
                            
                            # Area plot (like plotter.py)
                            plt.fill_between(contig_depth["pos"], contig_depth["depth"], 
                                           color="#2ecc71", alpha=0.6)
                            plt.plot(contig_depth["pos"], contig_depth["depth"], 
                                    color="#27ae60", linewidth=1)
                            
                            plt.title(f"Read Coverage: {contig_id}", fontsize=14)
                            plt.xlabel("Position (bp)", fontsize=12)
                            plt.ylabel("Depth", fontsize=12)
                            plt.grid(True, linestyle="--", alpha=0.3)
                            
                            # Add mean coverage line (like plotter.py)
                            mean_cov = contig_depth["depth"].mean()
                            plt.axhline(y=mean_cov, color="#e74c3c", linestyle="--", 
                                       alpha=0.8, label=f"Mean: {mean_cov:.1f}x")
                            plt.legend()
                            
                            # Sanitize filename
                            safe_name = contig_id.replace("|", "_").replace("/", "_")
                            plt.savefig(coverage_dir / f"coverage_{safe_name}.png", 
                                       dpi=300, bbox_inches="tight")
                            plt.close()
                            plots_generated += 1
                    
                    if plots_generated > 0:
                        print(f"Generated {plots_generated} coverage plots in {coverage_dir}", file=sys.stderr)
        except Exception as e:
            print(f"Could not generate coverage plots: {e}", file=sys.stderr)

    print(f"Plots written to {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
