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
                
            # Create plot
            plt.figure(figsize=(10, 6))
            sns.histplot(data=df, x="contig_length", log_scale=True, kde=True, color="teal")
            
            plt.title("Distribution of Viral Contig Lengths", fontsize=16)
            plt.xlabel("Contig Length (bp, log scale)", fontsize=12)
            plt.ylabel("Count", fontsize=12)
            
            # Add N50 line if possible
            # (Simple calculation for visualization purposes)
            lengths = sorted(df['contig_length'].tolist(), reverse=True)
            total_len = sum(lengths)
            cum_len = 0
            n50 = 0
            for l in lengths:
                cum_len += l
                if cum_len >= total_len / 2:
                    n50 = l
                    break
            
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
            
            # Pad shorter lineages with None
            padded_data = [d + [None] * (max_depth - len(d)) for d in data]
            
            df_lineage = pd.DataFrame(padded_data, columns=cols)
            
            # Add a count column
            df_lineage['count'] = 1
            
            # Create sunburst
            # We use the columns that have enough data
            valid_cols = []
            for col in cols:
                if df_lineage[col].notna().sum() > 0:
                    valid_cols.append(col)
            
            if not valid_cols:
                return None
                
            # Find the best column for color (first one with >1 unique values)
            color_col = valid_cols[0]  # Default to root
            
            # Skip the first column (usually 'Viruses' or 'NA') if possible
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
        Generate genome quality plots from CheckV results.
        
        Args:
            quality_file: Path to quality_summary.tsv
            
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
                
            # Check for required columns
            required_cols = ['completeness', 'contamination', 'checkv_quality']
            if not all(col in df.columns for col in required_cols):
                logger.warning(f"Quality file missing required columns: {required_cols}")
                return generated_plots
                
            # 1. Scatter Plot: Completeness vs Contamination
            plt.figure(figsize=(10, 8))
            
            # Create scatter plot
            sns.scatterplot(
                data=df,
                x='completeness',
                y='contamination',
                hue='checkv_quality',
                style='checkv_quality',
                s=100,
                alpha=0.7,
                palette='viridis'
            )
            
            plt.title("Viral Genome Quality: Completeness vs. Contamination", fontsize=16)
            plt.xlabel("Completeness (%)", fontsize=12)
            plt.ylabel("Contamination (%)", fontsize=12)
            plt.grid(True, linestyle='--', alpha=0.3)
            plt.xlim(0, 105)  # Completeness is 0-100
            plt.ylim(-5, 105) # Contamination is 0-100
            
            # Add quality zones
            plt.axvline(x=90, color='green', linestyle=':', alpha=0.5, label='High Quality (>90%)')
            plt.axvline(x=50, color='orange', linestyle=':', alpha=0.5, label='Medium Quality (>50%)')
            
            output_scatter = self.plots_dir / "genome_quality_scatter.png"
            plt.savefig(output_scatter, dpi=300, bbox_inches='tight')
            plt.close()
            
            generated_plots.append(str(output_scatter))
            logger.info(f"Generated quality scatter plot: {output_scatter}")
            
            # 2. Bar Chart: Quality Tier Distribution
            plt.figure(figsize=(10, 6))
            
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
            
            return generated_plots
            
        except Exception as e:
            logger.error(f"Failed to generate quality plots: {e}")
            return generated_plots
