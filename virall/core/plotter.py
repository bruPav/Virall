"""
Module for generating plots from viral analysis results.
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Union
from loguru import logger

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Set backend to non-interactive to avoid Qt warnings
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.io as pio

class ViralPlotter:
    """
    Generates visualization plots for viral analysis results.
    """
    
    def __init__(self, output_dir: Union[str, Path]):
        """
        Initialize the viral plotter.
        
        Args:
            output_dir: Base directory for output files
        """
        self.output_dir = Path(output_dir)
        self.plots_dir = self.output_dir / "07_plots"
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        
        # Set plot style
        sns.set_theme(style="whitegrid")
        plt.rcParams.update({'figure.autolayout': True})
        
    def plot_abundance(self, abundance_file: Union[str, Path]) -> Optional[str]:
        """
        Generate a bar chart of viral abundance.
        
        Args:
            abundance_file: Path to contig_abundance.tsv file
            
        Returns:
            Path to the generated plot file, or None if failed
        """
        if not os.path.exists(abundance_file):
            logger.warning(f"Abundance file not found: {abundance_file}")
            return None
            
        try:
            # Read data
            df = pd.read_csv(abundance_file, sep='\t')
            
            if df.empty:
                logger.warning("Abundance file is empty")
                return None
                
            # Check required columns
            required_cols = ['species', 'relative_abundance']
            if not all(col in df.columns for col in required_cols):
                logger.warning(f"Abundance file missing required columns: {required_cols}")
                return None
                
            # Group by species and sum abundance (in case of multiple contigs per species)
            species_abundance = df.groupby('species')['relative_abundance'].sum().sort_values(ascending=False)
            
            # Take top 20
            top_species = species_abundance.head(20)
            
            # Create plot
            plt.figure(figsize=(10, 8))
            sns.barplot(x=top_species.values, y=top_species.index, hue=top_species.index, legend=False, palette="viridis")
            
            plt.title("Top 20 Most Abundant Viral Species", fontsize=16)
            plt.xlabel("Relative Abundance", fontsize=12)
            plt.ylabel("Species", fontsize=12)
            
            # Save plot
            output_file = self.plots_dir / "viral_abundance.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Generated abundance plot: {output_file}")
            return str(output_file)
            
        except Exception as e:
            logger.error(f"Failed to generate abundance plot: {e}")
            return None
            
    def plot_contig_lengths(self, abundance_file: Union[str, Path]) -> Optional[str]:
        """
        Generate a histogram of contig lengths.
        
        Args:
            abundance_file: Path to contig_abundance.tsv file (contains length info)
            
        Returns:
            Path to the generated plot file, or None if failed
        """
        if not os.path.exists(abundance_file):
            logger.warning(f"Abundance file not found: {abundance_file}")
            return None
            
        try:
            # Read data
            df = pd.read_csv(abundance_file, sep='\t')
            
            if df.empty:
                logger.warning("Abundance file is empty")
                return None
                
            # Check required columns
            if 'contig_length' not in df.columns:
                logger.warning("Abundance file missing 'contig_length' column")
                return None
                
            # Filter valid lengths (positive values only for log scale)
            valid_df = df[df['contig_length'] > 0].copy()
            
            if valid_df.empty:
                logger.warning("No valid contig lengths found (all <= 0)")
                return None
                
            # Create plot
            plt.figure(figsize=(10, 6))
            
            # Check if we have enough data for KDE and log scale
            use_log = True
            use_kde = True
            
            if len(valid_df) < 5:
                use_kde = False
            
            # If range is small, log scale might not be useful or could cause issues
            if valid_df['contig_length'].max() / valid_df['contig_length'].min() < 10:
                use_log = False
            
            try:
                sns.histplot(data=valid_df, x="contig_length", log_scale=use_log, kde=use_kde, color="teal")
            except Exception as e:
                logger.warning(f"Seaborn plot failed with log_scale={use_log}, kde={use_kde}: {e}. Retrying with simple hist.")
                # Fallback to simple histogram
                plt.clf()
                plt.hist(valid_df['contig_length'], bins=30, color="teal", alpha=0.7)
                if use_log:
                    plt.xscale('log')
            
            plt.title("Distribution of Viral Contig Lengths", fontsize=16)
            plt.xlabel(f"Contig Length (bp{' , log scale' if use_log else ''})", fontsize=12)
            plt.ylabel("Count", fontsize=12)
            
            # Add N50 line if possible
            # (Simple calculation for visualization purposes)
            lengths = sorted(valid_df['contig_length'].tolist(), reverse=True)
            total_len = sum(lengths)
            cum_len = 0
            n50 = 0
            for l in lengths:
                cum_len += l
                if cum_len >= total_len / 2:
                    n50 = l
                    break
            
            if n50 > 0:
                plt.axvline(n50, color='red', linestyle='--', label=f'N50: {n50} bp')
                plt.legend()
            
            # Save plot
            output_file = self.plots_dir / "contig_length_distribution.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Generated contig length plot: {output_file}")
            return str(output_file)
            
        except Exception as e:
            logger.error(f"Failed to generate contig length plot: {e}")
            return None

    def plot_taxonomy_sunburst(self, summary_file: Union[str, Path]) -> Optional[str]:
        """
        Generate an interactive sunburst chart of viral taxonomy.
        
        Args:
            summary_file: Path to kaiju_summary.tsv file
            
        Returns:
            Path to the generated HTML plot file, or None if failed
        """
        if not os.path.exists(summary_file):
            logger.warning(f"Summary file not found: {summary_file}")
            return None
            
        try:
            # Read data
            df = pd.read_csv(summary_file, sep='\t')
            
            if df.empty:
                logger.warning("Summary file is empty")
                return None
                
            # Check for lineage column
            if 'lineage' not in df.columns:
                logger.warning("Summary file missing 'lineage' column. Cannot generate sunburst.")
                return None
                
            # Process lineage data
            # Lineage format: superkingdom,phylum,class,order,family,genus,species
            # But it's a semicolon separated string, and some parts might be empty or missing
            
            data = []
            
            for lineage in df['lineage'].dropna():
                # Split by semicolon and filter out ONLY empty parts (keep "NA" as valid taxonomy)
                # Lineages often end with ";" which creates empty strings
                parts = [p.strip() for p in lineage.split(';') if p.strip()]
                if not parts:
                    continue
                    
                # We need at least a few levels for a meaningful sunburst
                # Let's normalize to a fixed depth or just use what we have
                # For simplicity, we'll take the available parts up to species
                
                # Kaiju with -r outputs specific ranks. 
                # If a rank is missing, it might be empty string or just skipped.
                # We'll treat the parts as a hierarchy regardless of exact rank name for visualization
                
                # Add 'Viruses' as root if not present (usually first part is 'Viruses' or 'cellular organisms' etc)
                if parts[0] != 'Viruses':
                    parts.insert(0, 'Viruses')
                
                data.append(parts)
            
            if not data:
                logger.warning("No valid lineage data found")
                return None
                
            # Convert to DataFrame for plotly
            # We need columns for the hierarchy
            max_depth = max(len(d) for d in data)
            cols = [f"Level_{i}" for i in range(max_depth)]
            
            # Pad shorter lineages with None (plotly handles None/NaN fine)
            padded_data = [d + [None] * (max_depth - len(d)) for d in data]
            df_lineage = pd.DataFrame(padded_data, columns=cols)
            df_lineage['count'] = 1
            
            # Only replace empty strings with None - keep "NA" as valid taxonomy level
            for col in cols:
                df_lineage[col] = df_lineage[col].replace('', None)
            
            # Find valid columns: plotly sunburst can't handle None in paths,
            # so we only include columns where ALL rows have data (no nulls).
            # We stop at the first column that has any nulls.
            valid_cols = []
            for col in cols:
                if df_lineage[col].isna().sum() == 0:
                    valid_cols.append(col)
                else:
                    break  # Stop at first column with nulls
            
            if not valid_cols:
                return None
                
            # Find the best column for color (first one with >1 unique values)
            color_col = valid_cols[0]  # Default to root
            
            # Skip the first column (usually 'Viruses') if possible
            start_idx = 1 if len(valid_cols) > 1 else 0
            
            for col in valid_cols[start_idx:]:
                # Check number of unique values (excluding None/NaN)
                n_unique = df_lineage[col].nunique()
                if n_unique > 1:
                    color_col = col
                    break
                
            fig = px.sunburst(
                df_lineage, 
                path=valid_cols, 
                values='count',
                color=color_col,
                title="Viral Taxonomy Sunburst",
                width=800,
                height=800
            )
            
            fig.update_traces(textinfo="label+percent entry")
            
            # Save interactive HTML
            output_html = self.plots_dir / "viral_taxonomy_sunburst.html"
            fig.write_html(str(output_html))
            
            # Save static image (optional, requires kaleido)
            try:
                output_png = self.plots_dir / "viral_taxonomy_sunburst.png"
                fig.write_image(str(output_png))
            except Exception:
                logger.debug("Could not save static image (kaleido might be missing)")
            
            logger.info(f"Generated sunburst plot: {output_html}")
            return str(output_html)
            
        except Exception as e:
            logger.error(f"Failed to generate sunburst plot: {e}")
            return None

    def plot_genome_quality(self, quality_file: Union[str, Path]) -> List[str]:
        """
        Generate genome quality plots from merged CheckV + geNomad results.
        
        Plots generated:
        1. genome_quality_distribution.png - Quality tiers from both tools (or one if only one ran)
        2. genome_quality_scatter.png - Completeness vs Contamination (CheckV only, if available)
        3. completeness_by_source.png - Completeness distribution (both tools, or one if only one has data)
        
        Args:
            quality_file: Path to quality_summary.tsv (merged or original)
            
        Returns:
            List of paths to generated plot files
        """
        generated_plots = []
        
        if not os.path.exists(quality_file):
            logger.warning(f"Quality summary file not found: {quality_file}")
            return generated_plots
            
        try:
            df = pd.read_csv(quality_file, sep='\t')
            
            if df.empty:
                logger.warning("Quality summary file is empty")
                return generated_plots
                
            # Check for minimum required columns (checkv_quality is essential)
            if 'checkv_quality' not in df.columns:
                logger.warning("Quality file missing 'checkv_quality' column")
                return generated_plots
            
            # Determine data sources available
            has_quality_source = 'quality_source' in df.columns
            has_checkv = False
            has_genomad = False
            
            if has_quality_source:
                sources = df['quality_source'].unique()
                has_checkv = any(s in ['checkv', 'checkv_fallback'] for s in sources if pd.notna(s))
                has_genomad = 'genomad' in sources
            
            # =====================================================================
            # 1. Quality Distribution Plot (both tools, or one if only one ran)
            # =====================================================================
            plt.figure(figsize=(10, 6))
            
            if has_quality_source and (has_checkv or has_genomad):
                # Create grouped counts by source
                grouped = df.groupby(['checkv_quality', 'quality_source']).size().unstack(fill_value=0)
                
                if not grouped.empty:
                    ax = grouped.plot(kind='bar', stacked=True, figsize=(10, 6), 
                                     colormap='viridis', alpha=0.8)
                    
                    # Dynamic title based on what's available
                    if has_checkv and has_genomad:
                        title = "Distribution of Viral Genome Quality\n(CheckV for phages, geNomad for RNA/eukaryotic viruses)"
                    elif has_genomad:
                        title = "Distribution of Viral Genome Quality\n(geNomad assessment)"
                    else:
                        title = "Distribution of Viral Genome Quality\n(CheckV assessment)"
                    
                    plt.title(title, fontsize=14)
                    plt.xlabel("Quality Tier", fontsize=12)
                    plt.ylabel("Number of Genomes", fontsize=12)
                    plt.xticks(rotation=45, ha='right')
                    plt.legend(title='Quality Source', loc='upper right')
                    plt.tight_layout()
            else:
                # Simple bar chart (no quality_source column)
                quality_counts = df['checkv_quality'].value_counts()
                sns.barplot(
                    x=quality_counts.index, 
                    y=quality_counts.values,
                    hue=quality_counts.index,
                    legend=False,
                    palette='viridis'
                )
                plt.title("Distribution of Viral Genome Quality", fontsize=16)
                plt.xlabel("Quality Tier", fontsize=12)
                plt.ylabel("Number of Genomes", fontsize=12)
                plt.xticks(rotation=45)
            
            output_bar = self.plots_dir / "genome_quality_distribution.png"
            plt.savefig(output_bar, dpi=300, bbox_inches='tight')
            plt.close()
            
            generated_plots.append(str(output_bar))
            logger.info(f"Generated quality distribution plot: {output_bar}")
            
            # =====================================================================
            # 2. Scatter Plot: Completeness vs Contamination (CheckV only)
            # =====================================================================
            if 'completeness' in df.columns and 'contamination' in df.columns:
                # Filter to CheckV results only (they have contamination data)
                if has_quality_source:
                    checkv_df = df[df['quality_source'].isin(['checkv', 'checkv_fallback'])].copy()
                else:
                    checkv_df = df.copy()
                
                # Further filter to rows with valid completeness and contamination
                scatter_df = checkv_df[checkv_df['completeness'].notna() & checkv_df['contamination'].notna()].copy()
                
                if not scatter_df.empty:
                    plt.figure(figsize=(10, 8))
                    
                    sns.scatterplot(
                        data=scatter_df,
                        x='completeness',
                        y='contamination',
                        hue='checkv_quality',
                        style='checkv_quality',
                        s=100,
                        alpha=0.7,
                        palette='viridis'
                    )
                    
                    plt.title("Viral Genome Quality: Completeness vs. Contamination\n(CheckV assessment - phage genomes)", fontsize=14)
                    plt.xlabel("Completeness (%)", fontsize=12)
                    plt.ylabel("Contamination (%)", fontsize=12)
                    plt.grid(True, linestyle='--', alpha=0.3)
                    plt.xlim(0, 105)
                    plt.ylim(-5, 105)
                    
                    # Add quality zone lines
                    plt.axvline(x=90, color='green', linestyle=':', alpha=0.5, label='High Quality (>90%)')
                    plt.axvline(x=50, color='orange', linestyle=':', alpha=0.5, label='Medium Quality (>50%)')
                    
                    output_scatter = self.plots_dir / "genome_quality_scatter.png"
                    plt.savefig(output_scatter, dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    generated_plots.append(str(output_scatter))
                    logger.info(f"Generated quality scatter plot: {output_scatter} ({len(scatter_df)} CheckV contigs)")
                else:
                    logger.info("No CheckV contigs with completeness/contamination data for scatter plot")
            
            # =====================================================================
            # 3. Completeness Distribution by Source (both tools, or one if only one has data)
            # =====================================================================
            if 'completeness' in df.columns:
                completeness_df = df[df['completeness'].notna()].copy()
                
                if not completeness_df.empty:
                    plt.figure(figsize=(10, 6))
                    
                    if has_quality_source:
                        sources_with_data = completeness_df['quality_source'].unique()
                        n_sources = len(sources_with_data)
                        
                        if n_sources > 1:
                            # Multiple sources - box plot comparison
                            sns.boxplot(
                                data=completeness_df,
                                x='quality_source',
                                y='completeness',
                                hue='quality_source',
                                palette='viridis',
                                legend=False
                            )
                            plt.title("Completeness Distribution by Assessment Tool", fontsize=16)
                            plt.xlabel("Quality Assessment Tool", fontsize=12)
                        else:
                            # Single source - simple box plot
                            source_name = sources_with_data[0] if len(sources_with_data) > 0 else "Unknown"
                            sns.boxplot(
                                y=completeness_df['completeness'],
                                color='steelblue'
                            )
                            plt.title(f"Completeness Distribution ({source_name})", fontsize=16)
                            plt.xlabel("")
                    else:
                        # No quality_source column - simple box plot
                        sns.boxplot(
                            y=completeness_df['completeness'],
                            color='steelblue'
                        )
                        plt.title("Completeness Distribution", fontsize=16)
                        plt.xlabel("")
                    
                    plt.ylabel("Completeness (%)", fontsize=12)
                    plt.ylim(0, 105)
                    
                    output_completeness = self.plots_dir / "completeness_by_source.png"
                    plt.savefig(output_completeness, dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    generated_plots.append(str(output_completeness))
                    logger.info(f"Generated completeness plot: {output_completeness}")
            
            return generated_plots
            
        except Exception as e:
            logger.error(f"Failed to generate quality plots: {e}")
            return generated_plots

    def plot_gene_predictions(self, annotations_file: Union[str, Path], summary_file: Union[str, Path] = None, kaiju_summary_file: Union[str, Path] = None) -> List[str]:
        """
        Generate plots from gene prediction results.
        
        Args:
            annotations_file: Path to protein_annotations.tsv
            summary_file: Path to gene_prediction_summary.tsv (optional)
            kaiju_summary_file: Path to kaiju_summary.tsv (optional, for family coloring)
            
        Returns:
            List of paths to generated plot files
        """
        generated_plots = []
        
        # 1. Proteome Scatter Plot (MW vs pI) & Gene Lengths
        if os.path.exists(annotations_file):
            try:
                df_ann = pd.read_csv(annotations_file, sep='\t')
                
                if not df_ann.empty and 'molecular_weight' in df_ann.columns and 'isoelectric_point' in df_ann.columns:
                    plt.figure(figsize=(12, 8))
                    
                    # Try to add family info if available
                    hue_col = None
                    palette = None
                    
                    if kaiju_summary_file and os.path.exists(kaiju_summary_file):
                        try:
                            df_kaiju = pd.read_csv(kaiju_summary_file, sep='\t')
                            if 'contig_id' in df_kaiju.columns and 'lineage' in df_kaiju.columns:
                                # Extract family from lineage (5th element, index 4)
                                # Lineage: root; kingdom; phylum; class; order; family; ...
                                # But sometimes it's just a list of names.
                                # Let's try to parse it robustly.
                                
                                contig_to_family = {}
                                for _, row in df_kaiju.iterrows():
                                    lineage = str(row['lineage'])
                                    parts = [p.strip() for p in lineage.split(';') if p.strip()]
                                    
                                    # Heuristic: Family is usually the one ending in 'viridae'
                                    # Or we can just take a specific index if we trust the format.
                                    # The user mentioned "dots colored by viral Family".
                                    
                                    family = "Unknown"
                                    for part in parts:
                                        if part.endswith('viridae'):
                                            family = part
                                            break
                                    
                                    # If no viridae found, maybe use the last classified level if it's not NA?
                                    # But let's stick to explicit families for now to avoid clutter.
                                    
                                    contig_to_family[row['contig_id']] = family
                                
                                # Map to annotations
                                # Annotations have 'contig_id' usually? 
                                # Let's check the file content I viewed earlier.
                                # Yes: contig_id	gene_id	length ...
                                
                                df_ann['Family'] = df_ann['contig_id'].map(contig_to_family).fillna('Unknown')
                                
                                # Filter out Unknowns for the legend if too many?
                                # Or just plot them.
                                hue_col = 'Family'
                                # Use a categorical palette
                                palette = 'tab20'
                                
                        except Exception as e:
                            logger.warning(f"Failed to process Kaiju summary for plotting: {e}")
                    
                    # Scatter plot
                    sns.scatterplot(
                        data=df_ann,
                        x='isoelectric_point',
                        y='molecular_weight',
                        hue=hue_col,
                        palette=palette,
                        alpha=0.7,
                        edgecolor=None,
                        s=40
                    )
                    
                    if hue_col:
                        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
                    
                    plt.title("Viral Proteome Properties: MW vs. pI", fontsize=16)
                    plt.xlabel("Isoelectric Point (pI)", fontsize=12)
                    plt.ylabel("Molecular Weight (Da)", fontsize=12)
                    plt.grid(True, linestyle='--', alpha=0.3)
                    
                    # Add density contours if enough data
                    if len(df_ann) > 50:
                        sns.kdeplot(
                            data=df_ann,
                            x='isoelectric_point',
                            y='molecular_weight',
                            levels=5,
                            color="#e74c3c",
                            linewidths=1,
                            alpha=0.5
                        )
                    
                    output_proteome = self.plots_dir / "viral_proteome_scatter.png"
                    plt.savefig(output_proteome, dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    generated_plots.append(str(output_proteome))
                    logger.info(f"Generated proteome scatter plot: {output_proteome}")
                
                if not df_ann.empty and 'length' in df_ann.columns:
                    plt.figure(figsize=(10, 6))
                    sns.histplot(data=df_ann, x='length', bins=30, color='#3498db', kde=True)
                    plt.title("Distribution of Viral Gene Lengths", fontsize=16)
                    plt.xlabel("Gene Length (aa)", fontsize=12)
                    plt.ylabel("Count", fontsize=12)
                    
                    output_lengths = self.plots_dir / "gene_length_distribution.png"
                    plt.savefig(output_lengths, dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    generated_plots.append(str(output_lengths))
                    logger.info(f"Generated gene length plot: {output_lengths}")
                    
            except Exception as e:
                logger.error(f"Failed to process annotations file: {e}")
        
        # 2. Genes per Contig
        if summary_file and os.path.exists(summary_file):
            try:
                df_sum = pd.read_csv(summary_file, sep='\t')
                
                if not df_sum.empty and 'total_genes' in df_sum.columns:
                    plt.figure(figsize=(10, 6))
                    sns.histplot(data=df_sum, x='total_genes', bins=20, color='#9b59b6', discrete=True)
                    plt.title("Distribution of Genes per Viral Contig", fontsize=16)
                    plt.xlabel("Number of Genes", fontsize=12)
                    plt.ylabel("Count of Contigs", fontsize=12)
                    
                    output_genes = self.plots_dir / "genes_per_contig_distribution.png"
                    plt.savefig(output_genes, dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    generated_plots.append(str(output_genes))
                    logger.info(f"Generated genes per contig plot: {output_genes}")
                    
            except Exception as e:
                logger.error(f"Failed to process summary file: {e}")
                
        return generated_plots

    def plot_high_quality_coverage(self, quality_file: Union[str, Path], depth_file: Union[str, Path, List[Union[str, Path]]], max_plots: int = 10, plot_all: bool = False) -> List[str]:
        """
        Generate coverage plots for high-quality contigs (or all contigs if plot_all=True).
        
        For hybrid assemblies, combines coverage from multiple read types (long_reads, paired_reads, single_reads).
        
        Args:
            quality_file: Path to quality_summary.tsv
            depth_file: Path to contig_depth.txt, or list of paths to combine coverage from multiple sources
            max_plots: Maximum number of plots to generate (default: 10)
            plot_all: If True, plot all contigs regardless of quality (default: False)
            
        Returns:
            List of paths to generated plot files
        """
        generated_plots = []
        
        if not os.path.exists(quality_file):
            logger.warning("Quality summary file not found, skipping coverage plots")
            return generated_plots
        
        # Handle single file or list of files
        if isinstance(depth_file, (str, Path)):
            depth_files = [depth_file]
        else:
            depth_files = list(depth_file)
        
        # Filter to existing, non-empty files
        valid_depth_files = []
        for df in depth_files:
            if os.path.exists(df) and os.path.getsize(df) > 0:
                valid_depth_files.append(df)
        
        if not valid_depth_files:
            logger.warning("No valid depth files found, skipping coverage plots")
            return generated_plots
        
        logger.info(f"Combining coverage from {len(valid_depth_files)} source(s): {[Path(f).parent.name for f in valid_depth_files]}")
            
        try:
            # 1. Identify High-Quality Contigs
            df_quality = pd.read_csv(quality_file, sep='\t')
            
            if df_quality.empty or 'checkv_quality' not in df_quality.columns:
                logger.warning("Quality summary empty or missing checkv_quality column")
                return generated_plots
                
            # Filter for High-quality and Complete genomes unless plot_all is True
            if plot_all:
                hq_contigs = df_quality
                logger.info("Plotting coverage for all contigs (reference mode enabled)")
            else:
                hq_contigs = df_quality[df_quality['checkv_quality'].isin(['High-quality', 'Complete'])]
            
            if hq_contigs.empty:
                logger.info("No High-quality or Complete contigs found to plot coverage for")
                return generated_plots
            
            # Sort by length (longest first) and take top N
            if 'contig_length' in hq_contigs.columns:
                hq_contigs = hq_contigs.sort_values('contig_length', ascending=False)
            
            target_contigs = hq_contigs['contig_id'].head(max_plots).tolist()
            logger.info(f"Generating coverage plots for top {len(target_contigs)} high-quality/complete contigs")
            
            # 2. Read Depth Data
            # contig_depth.txt usually has no header: contig_id, pos, depth
            # It can be large, so we'll read it carefully or use pandas with chunking if needed
            # For now assuming it fits in memory or is reasonable for the filtered set
            
            # Create subdirectory for coverage plots
            coverage_dir = self.plots_dir / "coverage_plots"
            coverage_dir.mkdir(exist_ok=True)
            
            # 2. Read and combine depth data from all sources
            # For hybrid assemblies, we combine coverage from multiple read types
            # We'll aggregate depth by (contig_id, pos) and sum the depths
            
            chunk_size = 100000
            # Dictionary to store combined coverage: {(contig_id, pos): total_depth}
            combined_coverage = {}
            
            for depth_file in valid_depth_files:
                logger.debug(f"Reading coverage from {Path(depth_file).parent.name}")
                chunks = pd.read_csv(depth_file, sep='\t', header=None, names=['contig_id', 'pos', 'depth'], chunksize=chunk_size)
                
                for chunk in chunks:
                    # Filter for target contigs only
                    relevant_data = chunk[chunk['contig_id'].isin(target_contigs)]
                    
                    if not relevant_data.empty:
                        # Group by contig and position, sum depths
                        for (contig_id, pos), group in relevant_data.groupby(['contig_id', 'pos']):
                            key = (contig_id, pos)
                            total_depth = group['depth'].sum()
                            combined_coverage[key] = combined_coverage.get(key, 0) + total_depth
            
            # Convert combined coverage to per-contig data structure
            current_contig_data = {contig: {'pos': [], 'depth': []} for contig in target_contigs}
            
            for (contig_id, pos), depth in combined_coverage.items():
                current_contig_data[contig_id]['pos'].append(pos)
                current_contig_data[contig_id]['depth'].append(depth)
            
            # 3. Generate Plots
            for contig_id in target_contigs:
                data = current_contig_data[contig_id]
                if not data['pos']:
                    continue
                    
                # Sort by position just in case
                df_plot = pd.DataFrame(data).sort_values('pos')
                
                plt.figure(figsize=(12, 4))
                
                # Area plot
                plt.fill_between(df_plot['pos'], df_plot['depth'], color='#2ecc71', alpha=0.6)
                plt.plot(df_plot['pos'], df_plot['depth'], color='#27ae60', linewidth=1)
                
                plt.title(f"Read Coverage: {contig_id}", fontsize=14)
                plt.xlabel("Position (bp)", fontsize=12)
                plt.ylabel("Depth", fontsize=12)
                plt.grid(True, linestyle='--', alpha=0.3)
                
                # Add mean coverage line
                mean_cov = df_plot['depth'].mean()
                plt.axhline(y=mean_cov, color='#e74c3c', linestyle='--', alpha=0.8, label=f"Mean: {mean_cov:.1f}x")
                plt.legend()
                
                # Sanitize filename
                safe_name = contig_id.replace('|', '_').replace('/', '_')
                output_path = coverage_dir / f"coverage_{safe_name}.png"
                
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                generated_plots.append(str(output_path))
            
            if generated_plots:
                logger.info(f"Generated {len(generated_plots)} coverage plots in {coverage_dir}")
                
        except Exception as e:
            logger.error(f"Failed to generate coverage plots: {e}")
            
        return generated_plots
