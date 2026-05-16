#!/usr/bin/env python3
"""
Generate plots from Nextflow pipeline outputs WITH GENOME FRACTION CORRECTION.
Self-contained plotting - no dependency on virall Python package.
Builds contig_abundance.tsv and kaiju_summary.tsv from pipeline files,
then generates plots using matplotlib/plotly.

KEY IMPROVEMENT: Applies genome fraction correction to abundance estimation
to account for incomplete assembly.
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
    parser = argparse.ArgumentParser(
        description="Generate plots from Virall Nextflow outputs with genome fraction correction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
ABUNDANCE CORRECTION METHOD:
  1. For each contig, calculate genome fraction = max(completeness, checkv_completeness) / 100
  2. Adjusted length = contig_length / genome_fraction (if genome_fraction > 0)
  3. Calculate rpk = total_coverage / (adjusted_length_kb)
  4. TPM = (rpk / Σ(rpk_all)) × 1,000,000
  5. Relative abundance = TPM / 1,000,000
  
This corrects for incomplete assembly by estimating what coverage would be
if the entire viral genome were assembled.
"""
    )
    parser.add_argument("--quant-dir", type=Path, required=True, help="06_quantification directory (depth.txt)")
    parser.add_argument("--viral-contigs", type=Path, required=True, help="Viral contigs FASTA")
    parser.add_argument("--kaiju-dir", type=Path, required=True, help="03_classifications directory (kaiju_results_with_names.tsv)")
    parser.add_argument("--checkv-dir", type=Path, required=True, help="04_quality_assessment directory")
    parser.add_argument("--min-breadth", type=float, default=0.10, help="Minimum contig coverage breadth (fraction 0-1) used for abundance estimation")
    parser.add_argument("--out-dir", type=Path, required=True, help="07_plots output directory")
    parser.add_argument("--apply-genome-correction", action="store_true", default=True, help="Apply genome fraction correction (default: True)")
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

    # Parse Kaiju addTaxonNames output
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

    # Load quality data for genome fraction correction
    quality_file = args.checkv_dir / "quality_summary.tsv"
    if not quality_file.exists():
        for sub in ["contigs", "viral_contigs", ""]:
            cand = args.checkv_dir / sub / "quality_summary.tsv" if sub else quality_file
            if cand.exists():
                quality_file = cand
                break
    
    genome_fraction_by_contig = {}
    if quality_file.exists() and args.apply_genome_correction:
        try:
            qdf = pd.read_csv(quality_file, sep="\t")
            print(f"[run_plots] Loaded quality data from {quality_file}", file=sys.stderr)
            for _, row in qdf.iterrows():
                cid = str(row["contig_id"]).strip()
                # Prefer CheckV's checkv_completeness (AAI-based genome fraction).
                # geNomad's completeness measures hallmark gene presence, not genome
                # fraction — a 2 kb fragment with a marker gene is not 100% complete.
                checkv_val = float(row.get("checkv_completeness", 0) or 0)
                genomad_val = float(row.get("completeness", 0) or 0)
                if pd.notna(checkv_val) and checkv_val > 0:
                    completeness = checkv_val
                elif genomad_val > 0:
                    completeness = genomad_val
                else:
                    completeness = 0
                genome_fraction = completeness / 100.0 if completeness > 0 else 0.01
                genome_fraction_by_contig[cid] = genome_fraction
                print(f"[run_plots] {cid}: completeness={completeness}%, genome_fraction={genome_fraction:.3f}", file=sys.stderr)
        except Exception as e:
            print(f"[run_plots] Warning: Could not load genome fraction data: {e}", file=sys.stderr)
            args.apply_genome_correction = False

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
        
        # Apply genome fraction correction if available
        if args.apply_genome_correction and genome_fraction_by_contig:
            total_cov["genome_fraction"] = total_cov["contig_id"].map(
                lambda x: genome_fraction_by_contig.get(str(x).strip(), 1.0)
            )
            total_cov["genome_fraction"] = total_cov["genome_fraction"].clip(lower=0.01, upper=1.0)
            total_cov["adjusted_length"] = total_cov["contig_length"] / total_cov["genome_fraction"]
            print(f"[run_plots] Applied genome fraction correction to {len(total_cov)} contigs", file=sys.stderr)
            print(f"[run_plots] Genome fractions: {dict(zip(total_cov['contig_id'], total_cov['genome_fraction'].round(3)))}", file=sys.stderr)
        else:
            total_cov["genome_fraction"] = 1.0
            total_cov["adjusted_length"] = total_cov["contig_length"]
            print(f"[run_plots] No genome fraction correction applied", file=sys.stderr)
        
        total_cov["contig_length_kb"] = total_cov["contig_length"].replace(0, 1) / 1000.0
        total_cov["adjusted_length_kb"] = total_cov["adjusted_length"].replace(0, 1) / 1000.0
        
        total_cov["mean_coverage"] = (total_cov["total_coverage"] / total_cov["contig_length"].replace(0, 1)).round(2)
        total_cov["adjusted_mean_coverage"] = (total_cov["total_coverage"] / total_cov["adjusted_length"].replace(0, 1)).round(2)
        
        total_cov["max_coverage"] = depth_df.groupby("contig_id")["depth"].max().reindex(total_cov["contig_id"]).values
        total_cov["coverage_breadth"] = (total_cov["mapped_positions"] / total_cov["contig_length"].replace(0, 1)).round(4)
        
        min_breadth = min(max(float(args.min_breadth), 0.0), 1.0)
        before_filter = len(total_cov)
        total_cov = total_cov[total_cov["coverage_breadth"] >= min_breadth].copy()
        print(f"[run_plots] breadth filter >= {min_breadth:.4f}: kept {len(total_cov)}/{before_filter} contigs", file=sys.stderr)

        # Calculate abundances with genome fraction correction
        sum_all = total_cov["total_coverage"].sum()
        total_cov["relative_abundance_coverage"] = (total_cov["total_coverage"] / sum_all).round(6) if sum_all > 0 else 0.0
        
        # Use adjusted length for rpk calculation
        total_cov["rpk"] = total_cov["total_coverage"] / total_cov["adjusted_length_kb"].replace(0, 1)
        rpk_sum = total_cov["rpk"].sum()
        total_cov["tpm"] = (total_cov["rpk"] / rpk_sum * 1_000_000).round(3) if rpk_sum > 0 else 0.0
        total_cov["relative_abundance"] = (total_cov["tpm"] / 1_000_000).round(6)
        
        # Also calculate uncorrected for comparison
        total_cov["rpk_uncorrected"] = total_cov["total_coverage"] / total_cov["contig_length_kb"].replace(0, 1)
        rpk_sum_uncorrected = total_cov["rpk_uncorrected"].sum()
        total_cov["tpm_uncorrected"] = (total_cov["rpk_uncorrected"] / rpk_sum_uncorrected * 1_000_000).round(3) if rpk_sum_uncorrected > 0 else 0.0
        total_cov["relative_abundance_uncorrected"] = (total_cov["tpm_uncorrected"] / 1_000_000).round(6)
        
        total_cov["species"] = total_cov["contig_id"].map(species_for_contig)
        matched = (total_cov["species"] != "Unknown").sum()
        print(f"[run_plots] abundance rows: {len(total_cov)} species matched: {matched} unknown: {(total_cov['species'] == 'Unknown').sum()}", file=sys.stderr)
        
        # Show correction impact
        print(f"[run_plots] Genome fraction correction impact:", file=sys.stderr)
        for _, row in total_cov.sort_values("relative_abundance", ascending=False).iterrows():
            gf = row.get("genome_fraction", 1.0)
            uncorrected = row["relative_abundance_uncorrected"]
            corrected = row["relative_abundance"]
            change = corrected - uncorrected
            pct_change = (change / uncorrected * 100) if uncorrected > 0 else 0
            print(f"[run_plots]   {row['contig_id']}: gf={gf:.3f}, uncorrected={uncorrected:.6f}, corrected={corrected:.6f}, change={change:+.6f} ({pct_change:+.1f}%)", file=sys.stderr)
        
        total_cov = total_cov.sort_values("relative_abundance", ascending=False)
        total_cov.to_csv(abundance_file, sep="\t", index=False)
        shutil.copy2(abundance_file, args.out_dir / "contig_abundance.tsv")

        # Species-level abundance (using corrected values)
        species_cov = total_cov.groupby("species", dropna=False).agg(
            contig_count=("contig_id", "count"),
            total_coverage=("total_coverage", "sum"),
            total_contig_length=("contig_length", "sum"),
            total_adjusted_length=("adjusted_length", "sum"),
            genome_fraction_weighted=("genome_fraction", "mean"),
            rpk=("rpk", "sum"),
            rpk_uncorrected=("rpk_uncorrected", "sum"),
        ).reset_index()
        
        # Corrected abundance
        species_rpk_sum = species_cov["rpk"].sum()
        species_cov["tpm"] = (species_cov["rpk"] / species_rpk_sum * 1_000_000).round(3) if species_rpk_sum > 0 else 0.0
        species_cov["relative_abundance"] = (species_cov["tpm"] / 1_000_000).round(6)
        
        # Uncorrected for comparison
        species_rpk_sum_uncorrected = species_cov["rpk_uncorrected"].sum()
        species_cov["tpm_uncorrected"] = (species_cov["rpk_uncorrected"] / species_rpk_sum_uncorrected * 1_000_000).round(3) if species_rpk_sum_uncorrected > 0 else 0.0
        species_cov["relative_abundance_uncorrected"] = (species_cov["tpm_uncorrected"] / 1_000_000).round(6)
        
        species_cov = species_cov.sort_values("tpm", ascending=False)
        species_cov.to_csv(species_abundance_file, sep="\t", index=False)
        shutil.copy2(species_abundance_file, args.out_dir / "species_abundance.tsv")
        
        # Print summary
        print(f"[run_plots] Corrected species abundances:", file=sys.stderr)
        for _, row in species_cov.iterrows():
            print(f"[run_plots]   {row['species']}: corrected={row['relative_abundance']:.6f}, uncorrected={row['relative_abundance_uncorrected']:.6f}, "
                  f"genome_fraction={row['genome_fraction_weighted']:.3f}, contigs={row['contig_count']}", file=sys.stderr)

    # 2. Build kaiju_summary.tsv for sunburst (contig_id, taxon_id, taxon_name, classification, lineage)
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
    # (Already loaded earlier for genome fraction correction)

    # 4. Generate plots using standalone matplotlib/seaborn implementation
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
    # Abundance Plot (matches plotter.py style) - show both corrected and uncorrected
    # =========================================================================
    if abundance_file.exists():
        df = pd.read_csv(abundance_file, sep="\t")
        if species_abundance_file.exists():
            sdf = pd.read_csv(species_abundance_file, sep="\t")
        else:
            sdf = pd.DataFrame()
        
        if not sdf.empty and "tpm" in sdf.columns and "species" in sdf.columns:
            # Plot corrected abundances
            top_species = sdf.sort_values("tpm", ascending=False).head(20).set_index("species")["tpm"]
            abundance_title = "Top 20 Most Abundant Viral Species (Genome-Fraction Corrected TPM)"
            abundance_xlabel = "TPM (genome-fraction corrected)"
            
            plt.figure(figsize=(12, 8))
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
            plt.savefig(args.out_dir / "viral_abundance_corrected.png", dpi=300, bbox_inches="tight")
            plt.close()
            print(f"Generated viral_abundance_corrected.png", file=sys.stderr)
            
            # Also plot comparison of corrected vs uncorrected
            if "tpm_uncorrected" in sdf.columns:
                comp_df = sdf[sdf["tpm"] > 0].copy()
                if len(comp_df) > 1:
                    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
                    
                    # Corrected
                    axes[0].barh(range(len(comp_df)), comp_df["tpm"].values, color="green", alpha=0.7)
                    axes[0].set_yticks(range(len(comp_df)))
                    axes[0].set_yticklabels(comp_df["species"].tolist(), fontsize=8)
                    axes[0].set_xlabel("TPM (corrected)", fontsize=12)
                    axes[0].set_title("Genome-Fraction Corrected", fontsize=14)
                    
                    # Uncorrected
                    axes[1].barh(range(len(comp_df)), comp_df["tpm_uncorrected"].values, color="red", alpha=0.7)
                    axes[1].set_yticks(range(len(comp_df)))
                    axes[1].set_yticklabels(comp_df["species"].tolist(), fontsize=8)
                    axes[1].set_xlabel("TPM (uncorrected)", fontsize=12)
                    axes[1].set_title("Uncorrected (standard)", fontsize=14)
                    
                    plt.suptitle("Abundance Correction: Genome Fraction Impact", fontsize=16, y=0.98)
                    plt.tight_layout()
                    plt.savefig(args.out_dir / "abundance_correction_comparison.png", dpi=300, bbox_inches="tight")
                    plt.close()
                    print(f"Generated abundance_correction_comparison.png", file=sys.stderr)
        
        # Contig Length Distribution
        if not df.empty and "contig_length" in df.columns:
            valid = df[df["contig_length"] > 0].copy()
            if not valid.empty:
                plt.figure(figsize=(10, 6))
                
                use_log = True
                use_kde = len(valid) >= 5
                if valid["contig_length"].max() / valid["contig_length"].min() < 10:
                    use_log = False
                
                if HAS_SEABORN:
                    try:
                        sns.histplot(data=valid, x="contig_length", log_scale=use_log, 
                                    kde=use_kde, color="teal")
                    except Exception:
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
    if HAS_PLOTLY and kaiju_summary_file.exists():
        try:
            kdf = pd.read_csv(kaiju_summary_file, sep="\t")
            if not kdf.empty and "lineage" in kdf.columns:
                data = []
                for lineage in kdf["lineage"].dropna():
                    parts = [p.strip() for p in lineage.split(";") if p.strip()]
                    if parts:
                        if parts[0] != "Viruses":
                            parts.insert(0, "Viruses")
                        data.append(parts)
                
                if data:
                    max_depth = max(len(d) for d in data)
                    cols = [f"Level_{i}" for i in range(max_depth)]
                    padded_data = [d + [None] * (max_depth - len(d)) for d in data]
                    df_lineage = pd.DataFrame(padded_data, columns=cols)
                    df_lineage["count"] = 1
                    
                    for col in cols:
                        df_lineage[col] = df_lineage[col].replace("", None)
                    
                    valid_cols = []
                    for col in cols:
                        if df_lineage[col].isna().sum() == 0:
                            valid_cols.append(col)
                        else:
                            break
                    
                    if valid_cols and len(df_lineage) > 0:
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
                        
                        try:
                            fig.write_image(str(args.out_dir / "viral_taxonomy_sunburst.png"))
                            print(f"Generated viral_taxonomy_sunburst.png", file=sys.stderr)
                        except Exception:
                            pass
        except Exception as e:
            print(f"Could not generate sunburst plot: {e}", file=sys.stderr)
    
    # =========================================================================
    # Quality Plots (merged CheckV + geNomad data) - matches plotter.py style
    # =========================================================================
    if quality_file.exists():
        qdf = pd.read_csv(quality_file, sep="\t")
        if not qdf.empty and "checkv_quality" in qdf.columns:
            # Quality Distribution
            plt.figure(figsize=(10, 6))
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
            
            # Completeness Distribution
            if "completeness" in qdf.columns:
                comp_df = qdf[qdf["completeness"].notna()].copy()
                if not comp_df.empty:
                    plt.figure(figsize=(10, 6))
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
    # Genome Fraction vs Abundance Plot (NEW: shows correction relationship)
    # =========================================================================
    if abundance_file.exists() and quality_file.exists():
        try:
            ab_df = pd.read_csv(abundance_file, sep="\t")
            q_df = pd.read_csv(quality_file, sep="\t")
            
            merged = pd.merge(ab_df, q_df, on="contig_id", how="inner", suffixes=('', '_quality'))
            
            if not merged.empty and "genome_fraction" in merged.columns and "relative_abundance" in merged.columns:
                plt.figure(figsize=(10, 8))
                
                if HAS_SEABORN:
                    sns.scatterplot(
                        data=merged,
                        x="genome_fraction",
                        y="relative_abundance",
                        hue="species",
                        style="species",
                        s=100,
                        alpha=0.8
                    )
                else:
                    for species in merged["species"].unique():
                        subset = merged[merged["species"] == species]
                        plt.scatter(subset["genome_fraction"], subset["relative_abundance"], 
                                   label=species, s=100, alpha=0.8)
                
                plt.xlabel("Genome Fraction (completeness)", fontsize=12)
                plt.ylabel("Relative Abundance (corrected)", fontsize=12)
                plt.title("Genome Fraction vs Abundance Correction", fontsize=16)
                plt.grid(True, alpha=0.3)
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.tight_layout()
                plt.savefig(args.out_dir / "genome_fraction_vs_abundance.png", dpi=300, bbox_inches="tight")
                plt.close()
                print(f"Generated genome_fraction_vs_abundance.png", file=sys.stderr)
        except Exception as e:
            print(f"Could not generate genome fraction vs abundance plot: {e}", file=sys.stderr)

    print(f"Plots written to {args.out_dir}", file=sys.stderr)
    print(f"Genome fraction correction {'APPLIED' if args.apply_genome_correction else 'NOT APPLIED'}", file=sys.stderr)


if __name__ == "__main__":
    main()