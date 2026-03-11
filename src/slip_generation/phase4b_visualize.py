"""
Phase 4b Visualization: Coastal Hazard Curves

Produces:
  1. Coastal hazard curves at each tide gauge site
  2. Comparison of hazard between sites
  3. Magnitude contribution breakdown

Usage:
    python phase4b_visualize.py
"""

import numpy as np
import json
from pathlib import Path
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, str(Path(__file__).parent))


def load_coastal_hazard_data(hazard_file='output/coastal_hazard/coastal_hazard_curves.json'):
    """Load coastal hazard curves from JSON."""
    with open(hazard_file) as f:
        return json.load(f)


def main():
    print("=" * 70)
    print("Phase 4b: Coastal Hazard Visualization")
    print("=" * 70)

    if not Path('output/coastal_hazard/coastal_hazard_curves.json').exists():
        print("Error: coastal hazard curves not found. Run phase4b_coastal_hazard_aggregation.py first.")
        return

    # Load data
    print("\nLoading coastal hazard curves...")
    site_hazards = load_coastal_hazard_data()

    # Plot 1: Hazard curves for all sites
    fig, ax = plt.subplots(figsize=(12, 7))

    colors = ['C0', 'C1', 'C2']
    for (site, color) in zip(sorted(site_hazards.keys()), colors):
        data = site_hazards[site]
        thresholds = np.array(data['inundation_thresholds_m'])
        aep = np.array(data['aep_per_year'])

        # Filter out zero AEP
        valid = aep > 0
        ax.semilogy(thresholds[valid], aep[valid], 'o-', color=color,
                   label=f"{site.replace('_', ' ')}", linewidth=2.5, markersize=6)

    ax.set_xlabel('Inundation Depth [m]', fontsize=13)
    ax.set_ylabel('Annual Exceedance Probability [1/year]', fontsize=13)
    ax.set_title('Coastal Tsunami Hazard Curves (KL-RQMC with Propagation)',
                fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(loc='upper right', fontsize=11)
    ax.set_ylim([1e-5, 1e-0])

    plt.tight_layout()
    output_file = 'output/phase4b_coastal_hazard_curves.png'
    plt.savefig(output_file, dpi=150)
    print(f"Saved: {output_file}")
    plt.close()

    # Plot 2: Ensemble distributions for each site
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, (site, ax) in enumerate(zip(sorted(site_hazards.keys()), axes)):
        data = site_hazards[site]
        samples = np.array(data['samples'])

        # Plot histogram
        ax.hist(samples, bins=40, color='steelblue', alpha=0.7, edgecolor='black')
        ax.set_xlabel('Inundation Depth [m]', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title(f"{site.replace('_', ' ')}", fontsize=12, fontweight='bold')
        ax.text(0.98, 0.97,
                f"Mean: {data['mean']:.2e} m\nMedian: {data['median']:.2e} m\nStd: {data['std']:.2e} m",
                transform=ax.transAxes, ha='right', va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                fontsize=10)

    fig.suptitle('Coastal Inundation Ensemble Distributions (Aggregated over Mw 8.0, 8.5, 9.0)',
                fontsize=13, fontweight='bold')
    plt.tight_layout()
    output_file = 'output/phase4b_coastal_distributions.png'
    plt.savefig(output_file, dpi=150)
    print(f"Saved: {output_file}")
    plt.close()

    # Plot 3: Site comparison at fixed AEP levels
    fig, ax = plt.subplots(figsize=(11, 6))

    aep_levels = [1e-1, 1e-2, 1e-3, 1e-4]
    site_names = sorted(site_hazards.keys())
    x_pos = np.arange(len(aep_levels))
    width = 0.25

    for i, site in enumerate(site_names):
        data = site_hazards[site]
        thresholds = np.array(data['inundation_thresholds_m'])
        aep = np.array(data['aep_per_year'])

        inunds_at_aep = []
        for aep_target in aep_levels:
            # Find inundation depth at this AEP
            idx = np.argmin(np.abs(aep - aep_target))
            inunds_at_aep.append(thresholds[idx])

        ax.bar(x_pos + i*width, inunds_at_aep, width, label=site.replace('_', ' '))

    ax.set_xlabel('Annual Exceedance Probability', fontsize=12)
    ax.set_ylabel('Inundation Depth [m]', fontsize=12)
    ax.set_title('Inundation Depth at Fixed Exceedance Probabilities',
                fontsize=13, fontweight='bold')
    ax.set_xticks(x_pos + width)
    ax.set_xticklabels([f"{aep:.0e}" for aep in aep_levels])
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    output_file = 'output/phase4b_coastal_comparison.png'
    plt.savefig(output_file, dpi=150)
    print(f"Saved: {output_file}")
    plt.close()

    print("\nPhase 4b visualization complete.")


if __name__ == '__main__':
    main()
