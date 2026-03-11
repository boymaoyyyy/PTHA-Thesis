"""
Phase 4 Visualization: Hazard Curves and Return Period Maps

This script produces visualizations of the aggregated hazard analysis:
  1. Hazard curves (AEP vs intensity) at selected sites
  2. Return period heatmaps across the observation domain
  3. Displacement distribution comparisons

Usage:
    python phase4_visualize.py
"""

import numpy as np
import json
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, str(Path(__file__).parent))


def load_hazard_results(hazard_curves_file='output/phase4_hazard_curves.json',
                        site_stats_file='output/phase4_site_statistics.json'):
    """Load hazard results from JSON files."""
    with open(hazard_curves_file) as f:
        hazard_data = json.load(f)

    with open(site_stats_file) as f:
        stats_data = json.load(f)

    return hazard_data, stats_data


def hazard_curves_plot(hazard_data, output_file='output/phase4_hazard_curves.png'):
    """Plot hazard curves at selected sites."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Reshape observation grid to map indices to lat/lon
    # 20x20 grid: index = i*20 + j
    n_obs = len(hazard_data)
    grid_size = int(np.sqrt(n_obs))

    # Select sites
    indices = [0, grid_size // 2, grid_size * grid_size // 2, grid_size * grid_size - 1]

    colors = ['C0', 'C1', 'C2', 'C3']
    labels = ['Southwest corner', 'Top-left mid', 'Center', 'Northeast corner']

    for idx, color, label in zip(indices, colors, labels):
        if str(idx) not in hazard_data:
            continue

        result = hazard_data[str(idx)]
        intensities = np.array(result['intensities_m'])
        aep = np.array(result['aep_per_year'])

        # Filter out zero AEP values for log scale
        valid = aep > 0
        ax.semilogy(intensities[valid], aep[valid], 'o-', color=color,
                   label=label, linewidth=2, markersize=5)

    ax.set_xlabel('Displacement Intensity [m]', fontsize=12)
    ax.set_ylabel('Annual Exceedance Probability [1/year]', fontsize=12)
    ax.set_title('Probabilistic Tsunami Hazard Curves (KL-based RQMC)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Saved: {output_file}")
    plt.close()


def return_period_heatmap(hazard_data, target_return_periods=[475, 2500],
                          output_prefix='output/phase4_'):
    """Create heatmaps of displacement at target return periods."""
    n_obs = len(hazard_data)
    grid_size = int(np.sqrt(n_obs))

    for return_period in target_return_periods:
        displacements = []

        for obs_idx in range(n_obs):
            result = hazard_data[str(obs_idx)]
            intensities = np.array(result['intensities_m'])
            return_periods = np.array(result['return_periods_yr'])

            # Interpolate: find displacement at this return period
            # (return_periods is sorted descending)
            if return_period > return_periods.min() and return_period < return_periods.max():
                disp_at_rp = np.interp(return_period, return_periods[::-1], intensities[::-1])
            else:
                disp_at_rp = np.nan

            displacements.append(disp_at_rp)

        displacements = np.array(displacements)

        # Reshape to grid
        displacement_grid = displacements.reshape(grid_size, grid_size)

        # Plot
        fig, ax = plt.subplots(figsize=(10, 9))

        im = ax.imshow(displacement_grid, origin='lower', cmap='viridis',
                      aspect='auto', interpolation='nearest')

        ax.set_xlabel('Lon index', fontsize=11)
        ax.set_ylabel('Lat index', fontsize=11)
        ax.set_title(f'Displacement at {return_period}-year Return Period\n(KL-based RQMC)', 
                    fontsize=12, fontweight='bold')

        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Displacement [m]', fontsize=11)

        plt.tight_layout()
        output_file = f"{output_prefix}return_period_{return_period}yr.png"
        plt.savefig(output_file, dpi=150)
        print(f"Saved: {output_file}")
        plt.close()


def displacement_statistics_plot(hazard_data, output_file='output/phase4_displacement_stats.png'):
    """Plot median, mean, and percentile displacements across domain."""
    n_obs = len(hazard_data)
    grid_size = int(np.sqrt(n_obs))

    medians = []
    means = []
    p05 = []
    p95 = []

    for obs_idx in range(n_obs):
        result = hazard_data[str(obs_idx)]

        medians.append(result.get('median_displacement_m', np.nan))
        means.append(result.get('mean_displacement_m', np.nan))
        p05.append(result.get('p05_displacement_m', np.nan))
        p95.append(result.get('p95_displacement_m', np.nan))

    medians = np.array(medians).reshape(grid_size, grid_size)
    means = np.array(means).reshape(grid_size, grid_size)
    p05 = np.array(p05).reshape(grid_size, grid_size)
    p95 = np.array(p95).reshape(grid_size, grid_size)

    fig, axes = plt.subplots(2, 2, figsize=(13, 11))

    plots = [
        (medians, 'Median', axes[0, 0]),
        (means, 'Mean', axes[0, 1]),
        (p05, '5th Percentile', axes[1, 0]),
        (p95, '95th Percentile', axes[1, 1]),
    ]

    for data, title, ax in plots:
        im = ax.imshow(data, origin='lower', cmap='viridis', aspect='auto')
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.set_xlabel('Lon index')
        ax.set_ylabel('Lat index')
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('[m]')

    fig.suptitle('Displacement Ensemble Statistics (KL-based RQMC)', 
                fontsize=13, fontweight='bold', y=1.00)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Saved: {output_file}")
    plt.close()



def main():
    print("=" * 70)
    print("Phase 4: Visualization")
    print("=" * 70)

    # Check if hazard curves file exists
    if not Path('output/phase4_hazard_curves.json').exists():
        print("Error: phase4_hazard_curves.json not found. Run phase4_hazard_aggregation.py first.")
        return

    # Load results
    print("\nLoading hazard results...")
    hazard_data, stats_data = load_hazard_results()
    print(f"  Loaded {len(hazard_data)} observation points")

    # Generate plots
    print("\nGenerating visualizations...")

    hazard_curves_plot(hazard_data)
    return_period_heatmap(hazard_data, target_return_periods=[475, 2500])
    displacement_statistics_plot(hazard_data)

    print("\nPhase 4 visualization complete.")


if __name__ == '__main__':
    import sys
    main()
