#!/usr/bin/env python3
"""
Threshold Sweep for Extendor Analysis
======================================
Loads pre-computed extendor proportions and scores, then filters by
different combinations of number_nonzero_samples, effect_size_bin,
and target_entropy to generate tumor-vs-normal scatter plots.

Usage:
    python threshold_sweep.py \
        --proportions extendor_proportions.tsv \
        --scores superset_scores.tsv \
        --inputFile input_CPTAC.txt \
        --outFolder extendor_threshold_sweep
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from itertools import product


def get_args():
    parser = argparse.ArgumentParser(description='Threshold sweep for extendor analysis')
    parser.add_argument('--proportions', required=True, help='Path to extendor_proportions.tsv')
    parser.add_argument('--scores', required=True, help='Path to superset_scores.tsv')
    parser.add_argument('--inputFile', required=True, help='Path to input_CPTAC.txt')
    parser.add_argument('--outFolder', required=True, help='Output folder')
    return parser.parse_args()


def load_data(proportions_path, scores_path):
    """Load proportions and scores data."""
    print("Loading extendor proportions...")
    prop_df = pd.read_csv(proportions_path, sep='\t')
    print(f"  Loaded {len(prop_df)} extendors")

    print("Loading scores...")
    scores_df = pd.read_csv(
        scores_path, sep='\t',
        usecols=['anchor', 'effect_size_bin', 'number_nonzero_samples', 'target_entropy']
    )
    print(f"  Loaded {len(scores_df)} anchors")

    return prop_df, scores_df


def filter_and_compute(prop_df, scores_df, min_samples, min_es, min_te):
    """
    Filter anchors by thresholds, then return filtered proportions.
    """
    # Filter scores
    mask = (scores_df['number_nonzero_samples'] >= min_samples) & \
           (scores_df['effect_size_bin'] >= min_es)
    if min_te is not None:
        mask = mask & (scores_df['target_entropy'] >= min_te)

    filtered_anchors = set(scores_df.loc[mask, 'anchor'])

    # Filter proportions to only these anchors
    filtered_prop = prop_df[prop_df['anchor'].isin(filtered_anchors)].copy()

    return filtered_prop, len(filtered_anchors)


def make_scatter(ax, prop_df, title, log_scale=True):
    """Create a scatter plot on the given axes."""
    df_plot = prop_df[(prop_df['count_tumor'] > 0) | (prop_df['count_normal'] > 0)].copy()

    if len(df_plot) == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title, fontsize=9)
        return

    # Recalculate proportions for this subset
    total_tumor = df_plot['count_tumor'].sum()
    total_normal = df_plot['count_normal'].sum()

    if total_tumor == 0 or total_normal == 0:
        ax.text(0.5, 0.5, 'No data in one group', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title, fontsize=9)
        return

    pt = df_plot['count_tumor'] / total_tumor
    pn = df_plot['count_normal'] / total_normal

    pseudo = 1e-10
    log2fc = np.log2((pt + pseudo) / (pn + pseudo))

    if log_scale:
        pseudo_plot = 1e-8
        x = pn + pseudo_plot
        y = pt + pseudo_plot
    else:
        x = pn
        y = pt

    sc = ax.scatter(
        x, y,
        c=log2fc.clip(-10, 10),
        cmap='RdBu_r',
        alpha=0.3,
        s=3,
        vmin=-10, vmax=10,
        rasterized=True
    )

    # Diagonal
    if log_scale:
        mn = min(x.min(), y.min())
        mx = max(x.max(), y.max())
        ax.plot([mn, mx], [mn, mx], 'k--', alpha=0.4, lw=0.5)
        ax.set_xscale('log')
        ax.set_yscale('log')
    else:
        mx = max(x.max(), y.max())
        ax.plot([0, mx], [0, mx], 'k--', alpha=0.4, lw=0.5)

    n_anchors = df_plot['anchor'].nunique()
    n_ext = len(df_plot)
    ax.set_title(title + f"\n{n_anchors} anchors, {n_ext:,} extendors", fontsize=8)
    ax.set_xlabel('Normal prop', fontsize=7)
    ax.set_ylabel('Tumor prop', fontsize=7)
    ax.tick_params(labelsize=6)

    return sc


def main():
    args = get_args()
    os.makedirs(args.outFolder, exist_ok=True)

    # Load data
    prop_df, scores_df = load_data(args.proportions, args.scores)

    # Define threshold grid
    samples_thresholds = [100, 200, 300, 400]
    es_thresholds = [0.80, 0.85, 0.90, 0.95]
    te_thresholds = [None, 0.5, 1.0, 1.5]  # None = no filter

    te_labels = {None: 'no filter', 0.5: '≥0.5', 1.0: '≥1.0', 1.5: '≥1.5'}

    # ================================================================
    # Summary table
    # ================================================================
    print("\nComputing summary table...")
    summary_rows = []

    for min_samples, min_es, min_te in product(samples_thresholds, es_thresholds, te_thresholds):
        filtered_prop, n_anchors = filter_and_compute(prop_df, scores_df, min_samples, min_es, min_te)

        n_extendors = len(filtered_prop)
        n_both = len(filtered_prop[(filtered_prop['count_tumor'] > 0) & (filtered_prop['count_normal'] > 0)])
        n_tumor_only = len(filtered_prop[(filtered_prop['count_tumor'] > 0) & (filtered_prop['count_normal'] == 0)])
        n_normal_only = len(filtered_prop[(filtered_prop['count_tumor'] == 0) & (filtered_prop['count_normal'] > 0)])

        summary_rows.append({
            'min_samples': min_samples,
            'min_effect_size': min_es,
            'min_target_entropy': min_te if min_te is not None else 'none',
            'n_anchors': n_anchors,
            'n_extendors': n_extendors,
            'n_both_groups': n_both,
            'n_tumor_only': n_tumor_only,
            'n_normal_only': n_normal_only,
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(args.outFolder, 'threshold_summary.tsv')
    summary_df.to_csv(summary_path, sep='\t', index=False)
    print(f"Saved threshold summary to {summary_path}")
    print(summary_df.to_string(index=False))

    # ================================================================
    # Grid of scatter plots: one page per target_entropy threshold
    # Rows = samples thresholds, Cols = effect_size thresholds
    # ================================================================
    print("\nGenerating scatter plot grids...")

    # Log scale plots
    pdf_path = os.path.join(args.outFolder, 'scatter_grid_log.pdf')
    with PdfPages(pdf_path) as pdf:
        for min_te in te_thresholds:
            fig, axes = plt.subplots(
                len(samples_thresholds), len(es_thresholds),
                figsize=(16, 16),
                squeeze=False
            )
            fig.suptitle(
                f'Extendor Proportions: Tumor vs Normal (Log Scale)\n'
                f'target_entropy {te_labels[min_te]}',
                fontsize=14, y=0.98
            )

            for i, min_samples in enumerate(samples_thresholds):
                for j, min_es in enumerate(es_thresholds):
                    ax = axes[i][j]
                    filtered_prop, n_anchors = filter_and_compute(
                        prop_df, scores_df, min_samples, min_es, min_te
                    )
                    title = f'samples≥{min_samples}, es≥{min_es}'
                    make_scatter(ax, filtered_prop, title, log_scale=True)

            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig(fig, dpi=150)
            plt.close(fig)
            print(f"  Page: target_entropy {te_labels[min_te]}")

    print(f"Saved log-scale scatter grid to {pdf_path}")

    # Also save individual PNGs for the "no TE filter" page
    print("\nGenerating individual scatter PNGs (no TE filter)...")
    png_dir = os.path.join(args.outFolder, 'scatter_plots')
    os.makedirs(png_dir, exist_ok=True)

    for min_samples, min_es in product(samples_thresholds, es_thresholds):
        filtered_prop, n_anchors = filter_and_compute(
            prop_df, scores_df, min_samples, min_es, None
        )

        fig, axes = plt.subplots(1, 2, figsize=(14, 7))

        # Log scale
        title = f'samples≥{min_samples}, es≥{min_es}'
        sc = make_scatter(axes[0], filtered_prop, title + ' (log)', log_scale=True)
        sc2 = make_scatter(axes[1], filtered_prop, title + ' (linear)', log_scale=False)

        # Add colorbar below the plots
        if sc is not None:
            fig.subplots_adjust(bottom=0.15)
            cbar_ax = fig.add_axes([0.15, 0.04, 0.7, 0.02])  # [left, bottom, width, height]
            cbar = fig.colorbar(sc, cax=cbar_ax, orientation='horizontal')
            cbar.set_label('log2(Tumor/Normal)', fontsize=10)

        fname = f'scatter_s{min_samples}_es{min_es:.2f}.png'
        fig.savefig(os.path.join(png_dir, fname), dpi=200, bbox_inches='tight')
        plt.close(fig)

    print(f"Saved individual PNGs to {png_dir}/")

    # ================================================================
    # Summary grid image (compact overview)
    # ================================================================
    print("\nGenerating compact summary grid...")
    fig, axes = plt.subplots(
        len(samples_thresholds), len(es_thresholds),
        figsize=(20, 20),
        squeeze=False
    )
    fig.suptitle(
        'Threshold Sweep: Tumor vs Normal Extendor Proportions (Log Scale, no TE filter)',
        fontsize=16, y=0.99
    )

    for i, min_samples in enumerate(samples_thresholds):
        for j, min_es in enumerate(es_thresholds):
            ax = axes[i][j]
            filtered_prop, n_anchors = filter_and_compute(
                prop_df, scores_df, min_samples, min_es, None
            )
            title = f'samples≥{min_samples}, es≥{min_es}'
            make_scatter(ax, filtered_prop, title, log_scale=True)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    grid_path = os.path.join(args.outFolder, 'scatter_grid_overview.png')
    fig.savefig(grid_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved overview grid to {grid_path}")

    print("\n=== Threshold sweep complete ===")


if __name__ == '__main__':
    main()
