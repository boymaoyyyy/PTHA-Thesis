"""
Phase 3 Results Visualization

This script loads the displacement fields generated in Phase 3 and produces:
  1. Terminal output with detailed statistics per magnitude
  2. PNG plots showing:
     - Distribution of maximum displacements per scenario
     - Example displacement field for each magnitude
     - Comparison across magnitudes
"""

import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors

sys.path.insert(0, str(Path(__file__).parent))

from philippine_trench_config import RUPTURE_SCENARIOS


def load_displacement_fields(mag):
    """Load all displacement .npy files for a given magnitude."""
    mag_str = f"Mw{mag:.1f}"
    disp_dir = Path('output') / 'tsunami_sources' / mag_str
    
    if not disp_dir.exists():
        return None
    
    disp_files = sorted(disp_dir.glob('*_displacement.npy'))
    displacements = [np.load(f) for f in disp_files]
    return np.array(displacements), disp_files


def print_stats(mag, displacements, disp_files):
    """Print statistics for a magnitude scenario."""
    if displacements is None or len(displacements) == 0:
        return
    
    # Compute max displacement for each sample
    max_disps = np.max(np.abs(displacements), axis=1)
    
    print(f"\n{'='*70}")
    print(f"Magnitude: Mw {mag:.1f}")
    print(f"{'='*70}")
    print(f"Total samples: {len(displacements)}")
    print(f"Observation points per sample: {displacements.shape[1]}")
    print(f"\nMaximum displacement statistics (across all observation points):")
    print(f"  Min:    {max_disps.min():12.3e} m")
    print(f"  Max:    {max_disps.max():12.3e} m")
    print(f"  Mean:   {max_disps.mean():12.3e} m")
    print(f"  Median: {np.median(max_disps):12.3e} m")
    print(f"  Std:    {max_disps.std():12.3e} m")
    print(f"\nPercentiles:")
    for p in [5, 25, 50, 75, 95]:
        val = np.percentile(max_disps, p)
        print(f"  {p:3d}%:  {val:12.3e} m")
    
    return max_disps


def plot_distributions(stats_dict):
    """Create boxplot comparing max displacement distributions across magnitudes."""
    mags = sorted(stats_dict.keys())
    data = [stats_dict[mag] for mag in mags]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    bp = ax.boxplot(data, labels=[f"Mw {m:.1f}" for m in mags],
                     patch_artist=True)
    
    # Color boxes
    colors_list = ['lightblue', 'lightgreen', 'lightyellow']
    for patch, color in zip(bp['boxes'], colors_list):
        patch.set_facecolor(color)
    
    ax.set_ylabel('Maximum Displacement [m]', fontsize=12)
    ax.set_xlabel('Magnitude', fontsize=12)
    ax.set_title('Distribution of Maximum Seafloor Displacements', fontsize=13)
    ax.grid(axis='y', alpha=0.3)
    ax.set_yscale('log')
    
    plt.tight_layout()
    out_path = Path('output') / 'phase3_displacement_distributions.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved plot: {out_path}")
    plt.close()


def plot_example_fields(mag, displacements, disp_files):
    """Plot example displacement fields for a magnitude (sample 0, 100, 500, 999)."""
    if displacements is None or len(displacements) == 0:
        return
    
    indices = [0, min(100, len(displacements)-1), 
               min(500, len(displacements)-1), 
               len(displacements)-1]
    indices = list(set(indices))  # Remove duplicates
    indices.sort()
    
    n_cols = min(2, len(indices))
    n_rows = (len(indices) + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows))
    if len(indices) == 1:
        axes = [axes]
    else:
        axes = axes.ravel()
    
    for i, idx in enumerate(indices):
        disp = np.abs(displacements[idx])
        im = axes[i].hist(disp, bins=30, color='steelblue', edgecolor='black', alpha=0.7)
        axes[i].set_xlabel('Displacement [m]', fontsize=10)
        axes[i].set_ylabel('Frequency', fontsize=10)
        axes[i].set_title(f'Sample {idx}: max = {disp.max():.2e} m', fontsize=11)
        axes[i].set_yscale('log')
    
    # Hide unused subplots
    for i in range(len(indices), len(axes)):
        axes[i].set_visible(False)
    
    fig.suptitle(f'Example Displacement Distributions — Mw {mag:.1f}', fontsize=13, y=1.00)
    plt.tight_layout()
    out_path = Path('output') / f'phase3_examples_Mw{mag:.1f}.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {out_path}")
    plt.close()


def main():
    import sys
    
    print("\nPhase 3: Tsunami Source Results")
    print("=" * 70)
    
    stats_dict = {}
    
    for scenario in RUPTURE_SCENARIOS.values():
        if scenario.magnitude not in [8.0, 8.5, 9.0]:
            continue
        
        mag = scenario.magnitude
        displacements, disp_files = load_displacement_fields(mag)
        
        if displacements is None:
            print(f"\nNo displacement data found for Mw {mag:.1f}")
            continue
        
        max_disps = print_stats(mag, displacements, disp_files)
        if max_disps is not None:
            stats_dict[mag] = max_disps
        
        plot_example_fields(mag, displacements, disp_files)
    
    if len(stats_dict) > 1:
        plot_distributions(stats_dict)
    
    print("\n" + "=" * 70)
    print("Visualizations complete.")


if __name__ == '__main__':
    import sys
    main()
