"""
RQMC vs Standard Monte Carlo Comparison

This script loads displacement results from both RQMC and MC Phase 3 runs
and produces comparison visualizations and statistics.
"""

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

def load_displacement_stats(source_dir):
    """Load all displacement .npy files and extract max displacement statistics."""
    results = {}
    
    for mag_dir in sorted(source_dir.glob('Mw*')):
        mag_str = mag_dir.name
        mag = float(mag_str[2:])
        
        disp_files = sorted(mag_dir.glob('*_displacement.npy'))
        if not disp_files:
            continue
        
        max_disps = [np.max(np.abs(np.load(f))) for f in disp_files]
        results[mag] = np.array(max_disps)
    
    return results


def print_comparison(rqmc_stats, mc_stats):
    """Print side-by-side comparison of RQMC vs MC statistics."""
    print("\nRQMC vs Monte Carlo Comparison")
    print("=" * 90)
    print(f"{'Magnitude':<12} | {'Method':<6} | {'Min [m]':<12} | {'Max [m]':<12} | "
          f"{'Mean [m]':<12} | {'Std [m]':<12}")
    print("-" * 90)
    
    all_mags = sorted(set(list(rqmc_stats.keys()) + list(mc_stats.keys())))
    
    for mag in all_mags:
        for method, stats_dict in [('RQMC', rqmc_stats), ('MC', mc_stats)]:
            if mag in stats_dict:
                data = stats_dict[mag]
                print(f"Mw {mag:.1f}      | {method:<6} | {data.min():.3e} | "
                      f"{data.max():.3e} | {data.mean():.3e} | {data.std():.3e}")


def plot_comparison(rqmc_stats, mc_stats):
    """Create comparison plots."""
    mags = sorted(set(list(rqmc_stats.keys()) + list(mc_stats.keys())))
    
    # Plot 1: Boxplot comparison
    fig, ax = plt.subplots(figsize=(12, 6))
    
    positions = []
    data_to_plot = []
    labels = []
    colors_list = []
    
    pos = 1
    for mag in mags:
        if mag in rqmc_stats:
            positions.append(pos)
            data_to_plot.append(rqmc_stats[mag])
            labels.append(f'RQMC\nMw{mag:.1f}')
            colors_list.append('lightblue')
            pos += 1
        
        if mag in mc_stats:
            positions.append(pos)
            data_to_plot.append(mc_stats[mag])
            labels.append(f'MC\nMw{mag:.1f}')
            colors_list.append('lightcoral')
            pos += 1
        
        pos += 0.5  # Gap between magnitude groups
    
    bp = ax.boxplot(data_to_plot, positions=positions, widths=0.6, patch_artist=True)
    
    for patch, color in zip(bp['boxes'], colors_list):
        patch.set_facecolor(color)
    
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('Maximum Displacement [m]', fontsize=11)
    ax.set_title('RQMC vs Monte Carlo: Displacement Distributions', fontsize=12)
    ax.grid(axis='y', alpha=0.3)
    ax.set_yscale('log')
    
    plt.tight_layout()
    out = Path('output') / 'rqmc_vs_mc_comparison.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {out}")
    plt.close()
    
    # Plot 2: Mean and std deviation comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    rqmc_means = [rqmc_stats[m].mean() for m in mags if m in rqmc_stats]
    mc_means = [mc_stats[m].mean() for m in mags if m in mc_stats]
    rqmc_stds = [rqmc_stats[m].std() for m in mags if m in rqmc_stats]
    mc_stds = [mc_stats[m].std() for m in mags if m in mc_stats]
    mags_rqmc = [m for m in mags if m in rqmc_stats]
    mags_mc = [m for m in mags if m in mc_stats]
    
    x = np.arange(len(mags_rqmc))
    width = 0.35
    
    ax1.bar(x - width/2, rqmc_means, width, label='RQMC', color='steelblue')
    ax1.bar(x + width/2, mc_means, width, label='MC', color='coral')
    ax1.set_ylabel('Mean Maximum Displacement [m]', fontsize=11)
    ax1.set_xlabel('Magnitude', fontsize=11)
    ax1.set_title('Mean Maximum Displacement', fontsize=12)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f'Mw {m:.1f}' for m in mags_rqmc])
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    ax1.set_yscale('log')
    
    ax2.bar(x - width/2, rqmc_stds, width, label='RQMC', color='steelblue')
    ax2.bar(x + width/2, mc_stds, width, label='MC', color='coral')
    ax2.set_ylabel('Standard Deviation [m]', fontsize=11)
    ax2.set_xlabel('Magnitude', fontsize=11)
    ax2.set_title('Variability in Maximum Displacement', fontsize=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'Mw {m:.1f}' for m in mags_rqmc])
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)
    ax2.set_yscale('log')
    
    plt.tight_layout()
    out = Path('output') / 'rqmc_vs_mc_statistics.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Saved: {out}")
    plt.close()


def main():
    rqmc_dir = Path('output') / 'tsunami_sources'
    mc_dir = Path('output') / 'tsunami_sources_mc'
    
    if not rqmc_dir.exists() or not mc_dir.exists():
        print("Error: Missing tsunami_sources directories")
        return
    
    rqmc_stats = load_displacement_stats(rqmc_dir)
    mc_stats = load_displacement_stats(mc_dir)
    
    print_comparison(rqmc_stats, mc_stats)
    plot_comparison(rqmc_stats, mc_stats)
    
    print("\nComparison complete. Check output/ for visualizations.")


if __name__ == '__main__':
    main()
