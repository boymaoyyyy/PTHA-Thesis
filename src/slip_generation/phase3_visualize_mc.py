"""
Phase 3 Results Visualization (MC)

This mirrors `phase3_visualize.py` but reads displacement fields
created from standard Monte Carlo and saves equivalent plots.
"""

import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))

from philippine_trench_config import RUPTURE_SCENARIOS


def load_displacement_fields_mc(mag):
    mag_str = f"Mw{mag:.1f}"
    disp_dir = Path('output') / 'tsunami_sources_mc' / mag_str
    if not disp_dir.exists():
        return None, None
    disp_files = sorted(disp_dir.glob('*_displacement.npy'))
    displacements = [np.load(f) for f in disp_files]
    return np.array(displacements), disp_files


def print_stats(mag, displacements):
    if displacements is None or len(displacements) == 0:
        return None
    max_disps = np.max(np.abs(displacements), axis=1)
    print(f"\n{'='*70}")
    print(f"Magnitude: Mw {mag:.1f} (MC)")
    print(f"{'='*70}")
    print(f"Total samples: {len(displacements)}")
    print(f"Observation points per sample: {displacements.shape[1]}")
    print(f"\nMaximum displacement statistics (across all observation points):")
    print(f"  Min:    {max_disps.min():12.3e} m")
    print(f"  Max:    {max_disps.max():12.3e} m")
    print(f"  Mean:   {max_disps.mean():12.3e} m")
    print(f"  Median: {np.median(max_disps):12.3e} m")
    print(f"  Std:    {max_disps.std():12.3e} m")
    return max_disps


def plot_distributions_mc(stats_dict):
    mags = sorted(stats_dict.keys())
    data = [stats_dict[mag] for mag in mags]
    fig, ax = plt.subplots(figsize=(10, 6))
    bp = ax.boxplot(data, labels=[f"Mw {m:.1f}" for m in mags], patch_artist=True)
    colors_list = ['lightcoral', 'lightsalmon', 'lightgray']
    for patch, color in zip(bp['boxes'], colors_list):
        patch.set_facecolor(color)
    ax.set_ylabel('Maximum Displacement [m]', fontsize=12)
    ax.set_xlabel('Magnitude', fontsize=12)
    ax.set_title('MC: Distribution of Maximum Seafloor Displacements', fontsize=13)
    ax.grid(axis='y', alpha=0.3)
    ax.set_yscale('log')
    plt.tight_layout()
    out_path = Path('output') / 'phase3_displacement_distributions_mc.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {out_path}")
    plt.close()


def plot_example_fields_mc(mag, displacements):
    if displacements is None or len(displacements) == 0:
        return
    indices = [0, min(100, len(displacements)-1), min(500, len(displacements)-1), len(displacements)-1]
    indices = list(sorted(set(indices)))
    n_cols = min(2, len(indices))
    n_rows = (len(indices) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows))
    if len(indices) == 1:
        axes = [axes]
    else:
        axes = axes.ravel()
    for i, idx in enumerate(indices):
        disp = np.abs(displacements[idx])
        axes[i].hist(disp, bins=30, color='teal', edgecolor='black', alpha=0.7)
        axes[i].set_xlabel('Displacement [m]', fontsize=10)
        axes[i].set_ylabel('Frequency', fontsize=10)
        axes[i].set_title(f'Sample {idx}: max = {disp.max():.2e} m', fontsize=11)
        axes[i].set_yscale('log')
    for i in range(len(indices), len(axes)):
        axes[i].set_visible(False)
    fig.suptitle(f'MC: Example Displacement Distributions — Mw {mag:.1f}', fontsize=13, y=1.00)
    plt.tight_layout()
    out_path = Path('output') / f'phase3_examples_Mw{mag:.1f}_mc.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {out_path}")
    plt.close()


def main():
    print("\nPhase 3: Tsunami Source Results (MC)")
    print("=" * 70)
    stats_dict = {}
    for scenario in RUPTURE_SCENARIOS.values():
        if scenario.magnitude not in [8.0, 8.5, 9.0]:
            continue
        mag = scenario.magnitude
        displacements, disp_files = load_displacement_fields_mc(mag)
        if displacements is None:
            print(f"No MC displacement data for Mw {mag:.1f}")
            continue
        max_disps = print_stats(mag, displacements)
        if max_disps is not None:
            stats_dict[mag] = max_disps
        plot_example_fields_mc(mag, displacements)
    if len(stats_dict) > 1:
        plot_distributions_mc(stats_dict)
    print("\nVisualization (MC) complete.")


if __name__ == '__main__':
    main()
