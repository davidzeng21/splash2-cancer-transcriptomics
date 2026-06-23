#!/usr/bin/env python3
"""
Extendor Analysis Script for CPTAC Lung Adenocarcinoma
======================================================
Extracts anchor+target (extendor) counts from SATC files,
aggregates by tumor vs normal groups using proportions,
and generates a scatter plot comparing extendor representation.

Author: Generated for CPTAC Lung Adeno matched tumor/normal analysis
"""

import numpy as np
import os
import glob
import subprocess
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
from tqdm import tqdm


def get_args():
    parser = argparse.ArgumentParser(
        description='Analyze extendor (anchor+target) proportions between tumor and normal samples'
    )
    parser.add_argument(
        "--inputFile",
        type=str,
        default="input_CPTAC.txt",
        help='Input file containing sample names and fastq paths (default: input_CPTAC.txt)'
    )
    parser.add_argument(
        "--scoresFile",
        type=str,
        default="2025-12-04_test.after_correction.scores.tsv",
        help='Scores TSV file from SPLASH output'
    )
    parser.add_argument(
        "--satcFolder",
        type=str,
        default="2025-12-04_test_satc",
        help='Folder containing .satc files'
    )
    parser.add_argument(
        "--sampleMapping",
        type=str,
        default="sample_name_to_id.mapping.txt",
        help='Path to sample_name_to_id.mapping.txt file'
    )
    parser.add_argument(
        "--satcDumpWrapper",
        type=str,
        default="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/tcga_endometrial/satc_dump_wrapper.sh",
        help='Path to satc_dump wrapper script'
    )
    parser.add_argument(
        "--outFolder",
        type=str,
        default="extendor_analysis_results",
        help='Output folder for results'
    )
    parser.add_argument(
        "--topN",
        type=int,
        default=10000,
        help='Number of top anchors to analyze (by effect size, default: 10000)'
    )
    parser.add_argument(
        "--skipSATCDump",
        action='store_true',
        default=False,
        help='Skip SATC dump step if counts already extracted'
    )
    parser.add_argument(
        "--anchorFile",
        type=str,
        default=None,
        help='Optional: provide a pre-made anchor file instead of extracting top anchors'
    )
    parser.add_argument(
        "--skipPropCalc",
        action='store_true',
        default=False,
        help='Skip SATC dump and proportion calculation; load existing extendor_proportions.tsv and re-run plotting and summary stats only'
    )

    return parser.parse_args()


def parse_sample_metadata(input_file):
    """
    Parse the input file to classify samples as tumor or normal.
    Returns two sets: tumor_samples, normal_samples
    """
    tumor_samples = set()
    normal_samples = set()
    
    with open(input_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 1:
                sample_name = parts[0]
                if sample_name.endswith('-tumor'):
                    tumor_samples.add(sample_name)
                elif sample_name.endswith('-normal'):
                    normal_samples.add(sample_name)
    
    print(f"Found {len(tumor_samples)} tumor samples and {len(normal_samples)} normal samples")
    return tumor_samples, normal_samples


def extract_top_anchors(scores_file, top_n, out_folder):
    """
    Extract top N anchors by effect size from the scores file.
    Creates an anchor file for satc_dump.
    """
    print(f"Reading scores file to extract top {top_n} anchors...")
    
    # Read only necessary columns to save memory
    # Use iterator to handle large file
    chunks = pd.read_csv(
        scores_file, 
        sep='\t', 
        usecols=['anchor', 'effect_size_bin', 'pval_opt_corrected'],
        chunksize=100000
    )
    
    # Collect anchor data
    anchor_data = []
    for chunk in tqdm(chunks, desc="Reading scores file"):
        anchor_data.append(chunk)
    
    df = pd.concat(anchor_data, ignore_index=True)
    print(f"Total anchors in scores file: {len(df)}")
    
    # Sort by effect size (descending) and get top N
    df_sorted = df.nlargest(top_n, 'effect_size_bin')
    
    # Save anchor list
    anchor_file = os.path.join(out_folder, 'anchors_to_analyze.txt')
    df_sorted['anchor'].to_csv(anchor_file, index=False, header=False)
    print(f"Saved {len(df_sorted)} anchors to {anchor_file}")
    
    return anchor_file, df_sorted


def run_satc_dump(satc_folder, sample_mapping, anchor_file, satc_dump_wrapper, out_folder):
    """
    Run satc_dump to extract per-sample counts for specified anchors.
    """
    dump_folder = os.path.join(out_folder, 'satc_dumps')
    Path(dump_folder).mkdir(parents=True, exist_ok=True)
    
    satc_files = glob.glob(os.path.join(satc_folder, 'bin*.satc'))
    print(f"Found {len(satc_files)} SATC files to process")
    
    for satc_file in tqdm(satc_files, desc="Running satc_dump"):
        bin_name = os.path.basename(satc_file).replace('.satc', '')
        out_file = os.path.join(dump_folder, f'{bin_name}.dump')
        
        # Skip if already exists
        if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
            continue
        
        cmd = [
            satc_dump_wrapper,
            '--sample_names', sample_mapping,
            '--anchor_list', anchor_file,
            '--format', 'splash',
            satc_file,
            out_file
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Warning: satc_dump failed for {satc_file}: {e.stderr}")
    
    return dump_folder


def load_dump_files(dump_folder):
    """
    Load all dump files and combine into a single DataFrame.
    """
    dump_files = glob.glob(os.path.join(dump_folder, '*.dump'))
    print(f"Loading {len(dump_files)} dump files...")
    
    df_list = []
    for dump_file in tqdm(dump_files, desc="Loading dump files"):
        if os.path.getsize(dump_file) > 0:
            try:
                df = pd.read_csv(
                    dump_file, 
                    sep='\t', 
                    names=['sample', 'anchor', 'target', 'count'],
                    dtype={'sample': str, 'anchor': str, 'target': str, 'count': int}
                )
                df_list.append(df)
            except Exception as e:
                print(f"Warning: Could not read {dump_file}: {e}")
    
    if not df_list:
        raise ValueError("No valid dump files found!")
    
    counts_df = pd.concat(df_list, ignore_index=True)
    print(f"Loaded {len(counts_df)} sample-anchor-target rows")
    
    return counts_df


def calculate_extendor_proportions(counts_df, tumor_samples, normal_samples):
    """
    Calculate proportions for each extendor (anchor+target) in tumor vs normal.
    """
    # Create extendor identifier
    counts_df['extendor'] = counts_df['anchor'] + '_' + counts_df['target']
    
    # Classify samples
    counts_df['group'] = counts_df['sample'].apply(
        lambda x: 'tumor' if x in tumor_samples else ('normal' if x in normal_samples else 'unknown')
    )
    
    # Filter out unknown samples
    counts_df = counts_df[counts_df['group'] != 'unknown']
    
    # Aggregate counts by extendor and group
    extendor_counts = counts_df.groupby(['extendor', 'anchor', 'target', 'group'])['count'].sum().reset_index()
    
    # Pivot to get tumor and normal counts side by side
    pivot_df = extendor_counts.pivot_table(
        index=['extendor', 'anchor', 'target'], 
        columns='group', 
        values='count', 
        fill_value=0
    ).reset_index()
    
    # Flatten column names - handle case where only one group has data
    pivot_df.columns.name = None
    cols = list(pivot_df.columns)
    # Ensure both tumor and normal columns exist
    if 'tumor' not in cols:
        pivot_df['tumor'] = 0
    if 'normal' not in cols:
        pivot_df['normal'] = 0
    pivot_df = pivot_df.rename(columns={'normal': 'count_normal', 'tumor': 'count_tumor'})
    
    # Calculate total counts per group
    total_tumor = pivot_df['count_tumor'].sum()
    total_normal = pivot_df['count_normal'].sum()
    
    print(f"Total counts - Tumor: {total_tumor:,}, Normal: {total_normal:,}")
    
    # Calculate proportions
    pivot_df['prop_tumor'] = pivot_df['count_tumor'] / total_tumor
    pivot_df['prop_normal'] = pivot_df['count_normal'] / total_normal
    
    # Calculate log10 fold change (with pseudocount to avoid division by zero)
    pseudocount = 1e-10
    pivot_df['log10FC'] = np.log10(
        (pivot_df['prop_tumor'] + pseudocount) / (pivot_df['prop_normal'] + pseudocount)
    )
    
    return pivot_df


def create_scatter_plot(prop_df, out_folder):
    """
    Create 2x2 scatter plot comparing extendor proportions and raw counts
    in tumor vs normal.
      Row 1: proportion linear | proportion log
      Row 2: raw counts linear | raw counts log
    All panels colored by log10(Tumor/Normal) with a dynamic symmetric colorbar.
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # Filter to extendors with at least some counts
    df_plot = prop_df[(prop_df['count_tumor'] > 0) | (prop_df['count_normal'] > 0)].copy()

    # Downsample for plotting to keep rendering fast
    MAX_PLOT_POINTS = 2_000_000
    if len(df_plot) > MAX_PLOT_POINTS:
        df_plot = df_plot.sample(n=MAX_PLOT_POINTS, random_state=42)
        print(f"Downsampled to {MAX_PLOT_POINTS:,} points for plotting "
              f"(original: {len(prop_df):,} extendors)")

    # Dynamic symmetric colorbar range covers the full data range
    colors = df_plot['log10FC'].values
    vlim = float(np.nanmax(np.abs(colors)))

    scatter_kwargs = dict(c=colors, cmap='RdBu_r', alpha=0.5, s=10,
                          vmin=-vlim, vmax=vlim)

    # ------------------------------------------------------------------ #
    # Row 1, Col 1 — Proportions, linear scale
    # ------------------------------------------------------------------ #
    ax = axes[0, 0]
    sc = ax.scatter(df_plot['prop_normal'], df_plot['prop_tumor'], **scatter_kwargs)
    max_val = max(df_plot['prop_normal'].max(), df_plot['prop_tumor'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='y=x')
    ax.set_xlabel('Proportion in Normal Tissue', fontsize=12)
    ax.set_ylabel('Proportion in Tumor Tissue', fontsize=12)
    ax.set_title('Extendor Proportions: Tumor vs Normal\n(Linear Scale)', fontsize=13)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('log10(Tumor/Normal)', fontsize=10)
    ax.legend(loc='upper left')

    # ------------------------------------------------------------------ #
    # Row 1, Col 2 — Proportions, log scale
    # ------------------------------------------------------------------ #
    ax = axes[0, 1]
    pseudo_prop = 1e-8
    x_log = df_plot['prop_normal'] + pseudo_prop
    y_log = df_plot['prop_tumor'] + pseudo_prop
    sc = ax.scatter(x_log, y_log, **scatter_kwargs)
    min_val = min(x_log.min(), y_log.min())
    max_val = max(x_log.max(), y_log.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='y=x')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Proportion in Normal Tissue (log scale)', fontsize=12)
    ax.set_ylabel('Proportion in Tumor Tissue (log scale)', fontsize=12)
    ax.set_title('Extendor Proportions: Tumor vs Normal\n(Log Scale)', fontsize=13)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('log10(Tumor/Normal)', fontsize=10)
    ax.legend(loc='upper left')

    # ------------------------------------------------------------------ #
    # Row 2, Col 1 — Raw counts, linear scale
    # ------------------------------------------------------------------ #
    ax = axes[1, 0]
    sc = ax.scatter(df_plot['count_normal'], df_plot['count_tumor'], **scatter_kwargs)
    max_val = max(df_plot['count_normal'].max(), df_plot['count_tumor'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='y=x')
    ax.set_xlabel('Count in Normal Tissue', fontsize=12)
    ax.set_ylabel('Count in Tumor Tissue', fontsize=12)
    ax.set_title('Extendor Raw Counts: Tumor vs Normal\n(Linear Scale)', fontsize=13)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('log10(Tumor/Normal)', fontsize=10)
    ax.legend(loc='upper left')

    # ------------------------------------------------------------------ #
    # Row 2, Col 2 — Raw counts, log scale
    # ------------------------------------------------------------------ #
    ax = axes[1, 1]
    pseudo_count = 1
    x_cnt_log = df_plot['count_normal'] + pseudo_count
    y_cnt_log = df_plot['count_tumor'] + pseudo_count
    sc = ax.scatter(x_cnt_log, y_cnt_log, **scatter_kwargs)
    min_val = min(x_cnt_log.min(), y_cnt_log.min())
    max_val = max(x_cnt_log.max(), y_cnt_log.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='y=x')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Count in Normal Tissue (log scale)', fontsize=12)
    ax.set_ylabel('Count in Tumor Tissue (log scale)', fontsize=12)
    ax.set_title('Extendor Raw Counts: Tumor vs Normal\n(Log Scale)', fontsize=13)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('log10(Tumor/Normal)', fontsize=10)
    ax.legend(loc='upper left')

    plt.tight_layout()

    out_path = os.path.join(out_folder, 'tumor_vs_normal_scatter.png')
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Saved scatter plot to {out_path}")

    plt.close()
    return out_path


def create_summary_stats(prop_df, tumor_samples, normal_samples, out_folder):
    """
    Generate summary statistics and additional visualizations.
    """
    # Summary stats
    stats = {
        'total_extendors': len(prop_df),
        'extendors_both_groups': len(prop_df[(prop_df['count_tumor'] > 0) & (prop_df['count_normal'] > 0)]),
        'extendors_tumor_only': len(prop_df[(prop_df['count_tumor'] > 0) & (prop_df['count_normal'] == 0)]),
        'extendors_normal_only': len(prop_df[(prop_df['count_tumor'] == 0) & (prop_df['count_normal'] > 0)]),
        'tumor_samples': len(tumor_samples),
        'normal_samples': len(normal_samples),
        'total_count_tumor': prop_df['count_tumor'].sum(),
        'total_count_normal': prop_df['count_normal'].sum(),
        'mean_log10FC': prop_df['log10FC'].mean(),
        'median_log10FC': prop_df['log10FC'].median(),
    }
    
    # Save summary
    summary_path = os.path.join(out_folder, 'summary_stats.txt')
    with open(summary_path, 'w') as f:
        f.write("Extendor Analysis Summary\n")
        f.write("=" * 50 + "\n\n")
        for key, value in stats.items():
            if isinstance(value, float):
                f.write(f"{key}: {value:.4f}\n")
            else:
                f.write(f"{key}: {value:,}\n")
    
    print(f"\nSummary Statistics:")
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.4f}")
        else:
            print(f"  {key}: {value:,}")
    
    # Log10FC histogram
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Filter reasonable values for histogram
    log10fc_filtered = prop_df['log10FC'][(prop_df['log10FC'] > -5) & (prop_df['log10FC'] < 5)]
    
    ax.hist(log10fc_filtered, bins=100, alpha=0.7, edgecolor='black')
    ax.axvline(x=0, color='red', linestyle='--', label='No change')
    ax.axvline(x=np.log10(2), color='green', linestyle=':', alpha=0.7, label='2x tumor enriched')
    ax.axvline(x=-np.log10(2), color='blue', linestyle=':', alpha=0.7, label='2x normal enriched')
    
    ax.set_xlabel('log10(Tumor/Normal)', fontsize=12)
    ax.set_ylabel('Number of Extendors', fontsize=12)
    ax.set_title('Distribution of Tumor vs Normal Fold Changes', fontsize=14)
    ax.legend()
    
    plt.tight_layout()
    hist_path = os.path.join(out_folder, 'log10FC_histogram.png')
    plt.savefig(hist_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return stats


def main():
    args = get_args()
    
    # Create output folder
    Path(args.outFolder).mkdir(parents=True, exist_ok=True)
    
    # Parse sample metadata (always needed for summary stats)
    tumor_samples, normal_samples = parse_sample_metadata(args.inputFile)

    prop_file = os.path.join(args.outFolder, 'extendor_proportions.tsv')

    if args.skipPropCalc:
        if not os.path.exists(prop_file):
            raise FileNotFoundError(
                f"--skipPropCalc set but proportions file not found: {prop_file}"
            )
        # Only load the 5 numeric columns needed for plotting and summary stats.
        # The string columns (extendor, anchor, target) are not used downstream
        # and would consume ~16 GB of memory for 209M rows.
        NEEDED_COLS = ['count_tumor', 'count_normal', 'prop_tumor', 'prop_normal', 'log10FC']
        import polars as pl
        parquet_path = prop_file.replace('.tsv', '.parquet')
        if os.path.exists(parquet_path):
            print(f"Loading numeric columns from parquet: {parquet_path}")
            prop_df = (
                pl.scan_parquet(parquet_path)
                  .select(NEEDED_COLS)
                  .collect(streaming=True)
                  .to_pandas()
            )
        else:
            print(f"Loading numeric columns from TSV: {prop_file}")
            prop_df = (
                pl.scan_csv(prop_file, separator='\t')
                  .select(NEEDED_COLS)
                  .collect(streaming=True)
                  .to_pandas()
            )
        print(f"Loaded {len(prop_df):,} extendors ({prop_df.memory_usage(deep=True).sum() / 1e9:.1f} GB in memory)")
    else:
        # Extract or use provided anchor file
        if args.anchorFile and os.path.exists(args.anchorFile):
            anchor_file = args.anchorFile
            print(f"Using provided anchor file: {anchor_file}")
        else:
            anchor_file, anchor_df = extract_top_anchors(args.scoresFile, args.topN, args.outFolder)

        # Run satc_dump if needed
        dump_folder = os.path.join(args.outFolder, 'satc_dumps')
        if not args.skipSATCDump:
            dump_folder = run_satc_dump(
                args.satcFolder,
                args.sampleMapping,
                anchor_file,
                args.satcDumpWrapper,
                args.outFolder
            )
        else:
            print(f"Skipping SATC dump, loading from {dump_folder}")

        # Load dump files
        counts_df = load_dump_files(dump_folder)

        # Calculate proportions
        prop_df = calculate_extendor_proportions(counts_df, tumor_samples, normal_samples)

        # Save proportions table (TSV + parquet for fast future reads)
        prop_df.to_csv(prop_file, sep='\t', index=False)
        parquet_path = prop_file.replace('.tsv', '.parquet')
        import polars as pl
        (
            pl.from_pandas(prop_df)
            .with_columns([
                pl.col('count_normal').cast(pl.Int32),
                pl.col('count_tumor').cast(pl.Int32),
                pl.col('prop_tumor').cast(pl.Float32),
                pl.col('prop_normal').cast(pl.Float32),
                pl.col('log10FC').cast(pl.Float32),
            ])
            .write_parquet(parquet_path)
        )
        print(f"Saved extendor proportions to {prop_file} and {parquet_path}")

    # Create scatter plot
    create_scatter_plot(prop_df, args.outFolder)
    
    # Create summary statistics
    create_summary_stats(prop_df, tumor_samples, normal_samples, args.outFolder)
    
    print("\nAnalysis complete!")


if __name__ == '__main__':
    main()

