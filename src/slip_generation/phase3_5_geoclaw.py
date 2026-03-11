"""
Phase 3.5 (GeoClaw Edition): Tsunami Propagation to Coastal Inundation

This script takes Okada-computed seafloor displacements and propagates them
to coastal sites using the GeoClaw tsunami solver (with linear fallback).

Workflow:
    1. Load Okada displacement fields from Phase 3
    2. For each displacement:
       - Initialize GeoClaw solver
       - Run nonlinear shallow-water propagation
       - Extract coastal inundation at 3 sites
    3. Save ensemble of coastal results
    4. Compute statistics

Output:
    - output/geoclaw_results/geoclaw_coastal_results.json
    - Coastal inundation ensembles at:
        * Dapa (lon=127.82°E, lat=9.20°N, depth=50 m)
        * General Santos (lon=125.37°E, lat=6.11°N, depth=80 m)
        * Zamboanga (lon=122.08°E, lat=6.93°N, depth=120 m)

Expected Runtime:
    - Linear fallback (no GeoClaw): ~5 minutes for 6000 samples
    - GeoClaw (if available): 30-60 minutes for 6000 samples

Author: Assistant
Date: 2025
"""

import numpy as np
from pathlib import Path
import json
from tqdm import tqdm
import sys

sys.path.insert(0, str(Path(__file__).parent))

from geoclaw_wrapper import GeoClawTsunamiSolver, run_geoclaw_batch_simulation
from philippine_trench_config import PHILIPPINE_TRENCH_PTHA_CONFIG


def phase3_5_geoclaw_propagation():
    """
    Execute Phase 3.5 with GeoClaw tsunami propagation.

    Loads Okada displacements from Phase 3 and propagates to coast.
    """
    cfg = PHILIPPINE_TRENCH_PTHA_CONFIG
    
    # Input/output paths
    tsunami_source_dir = Path('output') / 'tsunami_sources_okada'
    geoclaw_output_dir = Path('output') / 'geoclaw_results'
    geoclaw_output_dir.mkdir(parents=True, exist_ok=True)

    # Find all displacement files
    displacement_files = sorted(tsunami_source_dir.glob('Mw*/slip_Mw*_sample*_displacement_okada.npy'))
    
    if not displacement_files:
        print("⚠️  No Okada displacement files found!")
        print(f"   Expected path: {tsunami_source_dir}/Mw*/slip_Mw*_sample*_displacement_okada.npy")
        print("   Run phase3_okada.py first to generate displacements.")
        return None

    print("=" * 80)
    print("PHASE 3.5: GeoClaw Tsunami Propagation (Coastal Inundation)")
    print("=" * 80)
    print(f"Input: {len(displacement_files)} Okada displacement fields")
    print(f"Output: Coastal inundation at 3 sites via GeoClaw\n")

    # Check if GeoClaw is available
    solver_demo = GeoClawTsunamiSolver(np.zeros((10, 10)))
    if solver_demo.geoclaw_available:
        print("✓ GeoClaw (Clawpack) detected. Will use full nonlinear solver.")
        print("  Expected runtime: 30-60 minutes\n")
    else:
        print("⚠️  GeoClaw (Clawpack) not found. Using linear fallback propagation.")
        print("  (Install via: pip install clawpack)\n")

    # Run batch simulation
    coastal_results = run_geoclaw_batch_simulation(
        displacement_files,
        output_dir=str(geoclaw_output_dir)
    )

    print("\n" + "=" * 80)
    print("Phase 3.5 Complete: Coastal inundation ensembles saved")
    print("=" * 80)

    return coastal_results


def phase3_5_geoclaw_visualization():
    """
    Create visualizations of Phase 3.5 (GeoClaw) results.

    Plots coastal inundation distributions and arrival times.
    """
    geoclaw_results_file = Path('output') / 'geoclaw_results' / 'geoclaw_coastal_results.json'
    
    if not geoclaw_results_file.exists():
        print(f"Results file not found: {geoclaw_results_file}")
        return

    with open(geoclaw_results_file, 'r') as f:
        results = json.load(f)

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available. Skipping visualization.")
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Phase 3.5 (GeoClaw): Coastal Inundation Distributions', fontsize=14, fontweight='bold')

    for ax, (site, data) in zip(axes, results.items()):
        inundations = np.array(data['inundations_m'])
        
        ax.hist(inundations, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
        ax.axvline(data['mean_inundation_m'], color='red', linestyle='--', 
                   linewidth=2, label=f"Mean: {data['mean_inundation_m']:.2f} m")
        ax.axvline(data['median_inundation_m'], color='orange', linestyle='--', 
                   linewidth=2, label=f"Median: {data['median_inundation_m']:.2f} m")
        
        ax.set_title(site.replace('_', ' '), fontweight='bold')
        ax.set_xlabel('Coastal Inundation (m)')
        ax.set_ylabel('Frequency (sample count)')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file = Path('output') / 'phase3_5_geoclaw_inundations.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nVisualization saved: {output_file}")
    plt.close()

    # Create second figure: arrival times and amplification
    arrival_times_file = Path('output') / 'geoclaw_results' / 'geoclaw_coastal_results.json'
    if arrival_times_file.exists():
        with open(arrival_times_file, 'r') as f:
            arrival_data = json.load(f)

        fig, ax = plt.subplots(figsize=(10, 6))

        sites = list(arrival_data.keys())
        arrivals = [np.median(arrival_data[s]['arrival_times_min']) for s in sites]
        colors = ['#e74c3c', '#3498db', '#2ecc71']

        bars = ax.bar(sites, arrivals, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)

        for bar, arrival in zip(bars, arrivals):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.1f} min',
                   ha='center', va='bottom', fontweight='bold', fontsize=11)

        ax.set_ylabel('Arrival Time (minutes)', fontsize=12, fontweight='bold')
        ax.set_title('Phase 3.5 (GeoClaw): Median Wave Arrival Times', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_ylim(0, max(arrivals) * 1.2)

        plt.tight_layout()
        output_file = Path('output') / 'phase3_5_geoclaw_arrival_times.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Visualization saved: {output_file}")
        plt.close()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Phase 3.5: GeoClaw tsunami propagation to coastal sites'
    )
    parser.add_argument('--propagate', action='store_true', 
                       help='Run tsunami propagation simulation')
    parser.add_argument('--visualize', action='store_true',
                       help='Create visualization plots')
    parser.add_argument('--all', action='store_true',
                       help='Run both propagation and visualization')

    args = parser.parse_args()

    if args.all:
        args.propagate = True
        args.visualize = True

    if not (args.propagate or args.visualize):
        args.all = True
        args.propagate = True
        args.visualize = True

    if args.propagate:
        phase3_5_geoclaw_propagation()

    if args.visualize:
        phase3_5_geoclaw_visualization()

    print("\n✓ Phase 3.5 (GeoClaw) execution complete")
