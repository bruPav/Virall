#!/usr/bin/env python3
"""
Generate reference coverage plot and per-segment stats from samtools depth output.
Used by the Nextflow REFERENCE_CHECK process.
"""
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Generate reference coverage plot and per-segment stats"
    )
    parser.add_argument("--depth", required=True, help="samtools depth TSV file")
    parser.add_argument("--out-prefix", required=True,
                        help="Output prefix (e.g. ref_check_dir/reference)")
    args = parser.parse_args()

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import pandas as pd
    except ImportError as exc:
        print(f"Could not generate coverage plot: {exc}", file=sys.stderr)
        sys.exit(1)

    try:
        df = pd.read_csv(args.depth, sep='\t', header=None,
                         names=['contig', 'pos', 'depth'])

        if df.empty:
            print("Depth file is empty – no plot generated", file=sys.stderr)
            sys.exit(1)

        contigs = df['contig'].unique()
        n = len(contigs)

        if n > 1:
            fig, axes = plt.subplots(n, 1, figsize=(12, 2.5 * n), sharex=False)
            colors = plt.cm.tab10.colors
            for i, (ctg, ax) in enumerate(zip(contigs, axes)):
                sub = df[df['contig'] == ctg]
                c = colors[i % len(colors)]
                ax.fill_between(sub['pos'], sub['depth'], alpha=0.7, color=c)
                ax.set_xlim(0, sub['pos'].max())
                ax.set_ylim(0, None)
                ax.set_ylabel('Depth')
                label = ctg if len(ctg) <= 60 else ctg[:57] + '...'
                ax.set_title(label, fontsize=9, loc='left')
            axes[-1].set_xlabel('Position (bp)')
            fig.suptitle('Reference Genome Coverage (per segment)',
                         fontsize=11, y=1.0)
        else:
            fig, ax = plt.subplots(figsize=(12, 4))
            ax.fill_between(df['pos'], df['depth'], alpha=0.7,
                           color='steelblue')
            ax.set_xlabel('Position (bp)')
            ax.set_ylabel('Coverage Depth')
            ax.set_title('Reference Genome Coverage')
            ax.set_xlim(0, df['pos'].max())
            ax.set_ylim(0, None)

        plt.tight_layout()
        plt.savefig(f'{args.out_prefix}_coverage.png', dpi=150,
                   bbox_inches='tight')
        plt.close()

        # Write per-segment stats
        stats_file = f'{args.out_prefix}_stats.tsv'
        with open(stats_file, 'w') as f:
            f.write('segment\tlength\tcovered_bases\tbreadth_%\tmean_depth\n')
            for ctg in contigs:
                sub = df[df['contig'] == ctg]
                seg_len = int(sub['pos'].max())
                covered = int((sub['depth'] > 0).sum())
                breadth = covered / seg_len * 100 if seg_len > 0 else 0
                mean_d = sub['depth'].mean()
                f.write(f'{ctg}\t{seg_len}\t{covered}\t{breadth:.2f}\t{mean_d:.2f}\n')

        print(f"Coverage plot generated ({n} segment{'s' if n > 1 else ''})",
              file=sys.stderr)
    except Exception as e:
        print(f"Could not generate coverage plot: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
