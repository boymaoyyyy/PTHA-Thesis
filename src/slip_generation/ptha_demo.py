"""
Comprehensive PTHA demonstration: Full workflow from slip generation to hazard curves.

This example integrates all modules of the PTHA framework:
    Phase 1: RQMC slip generation
    Phase 2: Gutenberg-Richter magnitude frequency
    Phase 3: Tsunami source characterization
    Phase 4: Hazard aggregation and exceedance probability
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Headless backend
import matplotlib.pyplot as plt
from pathlib import Path

# Import PTHA modules
from slip_sampler import StochasticSlipGenerator
from magnitude_frequency import GutenbergRichterRelation, HazardScenarioWeights
from tsunami_source import SimpleDislocationSource, TsunamiPropagationInterface
from hazard_aggregation import HazardAggregator, compute_empirical_pdf_kde


# =====================================================================
# CONFIGURATION
# =====================================================================

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data" / "fault_geometry"
OUT_DIR = REPO_ROOT / "output"
PLOT_DIR = OUT_DIR / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

# PTHA parameters
MAGNITUDES = [7.5, 8.0, 8.5]
N_SLIP_SAMPLES_PER_MAGNITUDE = 50
SUBFAULT_SIZE_KM = 20.0
CORRELATION_LENGTH_KM = 20.0

# Gutenberg-Richter parameters (typical for subduction zones)
GR_A_VALUE = 5.0
GR_B_VALUE = 1.0
GR_MMIN = 7.5
GR_MMAX = 9.5

# Tsunami propagation
OBSERVATION_RANGE_KM = 200.0  # Max distance from fault
N_OBSERVATION_POINTS = 300


def phase1_rqmc_slip_generation():
    """Phase 1: Generate stochastic slip realizations using RQMC."""
    print("\n" + "=" * 70)
    print("PHASE 1: RQMC Slip Generation")
    print("=" * 70)

    slip_samples = {}

    for mag in MAGNITUDES:
        print(f"\nGenerating slips for Mw {mag}...")

        generator = StochasticSlipGenerator(
            magnitude=mag,
            fault_params_csv=DATA_DIR / "heidarzadeh_2025_table3.csv",
            subfault_size_km=SUBFAULT_SIZE_KM,
            covariance_model='exponential',
            correlation_length_km=CORRELATION_LENGTH_KM,
            rigidity_pa=4.0e10,
        )

        # Generate ensemble
        samples = generator.generate_ensemble(
            n_samples=N_SLIP_SAMPLES_PER_MAGNITUDE,
            seed=42 + mag,  # Reproducible seeding
        )

        slip_samples[mag] = {
            'generator': generator,
            'samples': samples,
        }

        # Print statistics
        all_mean_slips = np.array([np.mean(s) for s in samples])
        print(f"  Generated {len(samples)} samples")
        print(f"  Mean slip: {np.mean(all_mean_slips):.3f} m ± {np.std(all_mean_slips):.3f} m")
        print(f"  Max slip range: [{np.min([np.max(s) for s in samples]):.2f}, "
              f"{np.max([np.max(s) for s in samples]):.2f}] m")

    return slip_samples


def phase2_magnitude_frequency():
    """Phase 2: Set up Gutenberg-Richter distribution and scenario weights."""
    print("\n" + "=" * 70)
    print("PHASE 2: Magnitude-Frequency Distribution")
    print("=" * 70)

    # Create Gutenberg-Richter relation
    gr = GutenbergRichterRelation(
        a_value=GR_A_VALUE,
        b_value=GR_B_VALUE,
        M_min=GR_MMIN,
        M_max=GR_MMAX,
    )

    print(f"\nGutenberg-Richter parameters:")
    print(f"  a = {gr.a}, b = {gr.b}")
    print(f"  M_min = {gr.M_min}, M_max = {gr.M_max}")

    # Compute annual rates
    print(f"\nAnnual occurrence rates:")
    for mag in MAGNITUDES:
        rate = gr.annual_rate(mag)
        prob = gr.probability_given_scenario(mag)
        print(f"  M ≥ {mag}: {rate:.4e} events/year (P|M_min = {prob:.4f})")

    # Compute scenario weights
    weighter = HazardScenarioWeights(gr)
    mag_bins = [7.5, 8.0, 8.5, 9.0, 9.5]
    result = weighter.normalize_weights(mag_bins)

    print(f"\nScenario weights (for aggregation):")
    print(f"  Total rate: {result['total_annual_rate']:.4e} events/year")
    for i, (center, rate, prob) in enumerate(zip(
        result['bin_centers'],
        result['bin_weights'],
        result['bin_probabilities'],
    )):
        print(f"  Bin {i} (M~{center:.2f}): rate={rate:.4e}/yr, prob={prob:.4f}")

    # Extract weights for our scenarios
    scenario_weights = {mag: gr.annual_rate(mag) for mag in MAGNITUDES}

    return gr, scenario_weights


def phase3_tsunami_sources(slip_samples):
    """Phase 3: Convert slip distributions to tsunami sources."""
    print("\n" + "=" * 70)
    print("PHASE 3: Tsunami Source Characterization")
    print("=" * 70)

    tsunami_sources = {}

    for mag, slip_dict in slip_samples.items():
        gen = slip_dict['generator']
        samples = slip_dict['samples']

        print(f"\nCharacterizing tsunami sources for Mw {mag}...")

        # Create fault geometry dict
        fault_geom = {
            'length_km': gen.L_km,
            'width_km': gen.W_km,
            'top_depth_km': gen.top_depth_km,
            'strike_deg': gen.strike_deg,
            'dip_deg': gen.dip_deg,
            'rake_deg': gen.rake_deg,
            'n_subfaults_along': gen.n_along,
            'n_subfaults_down': gen.n_down,
        }

        # Create sources for each sample
        sources = []
        for slip in samples:
            source = SimpleDislocationSource(
                fault_geom,
                slip,
                rigidity=gen.rigidity,
            )
            sources.append(source)

        tsunami_sources[mag] = {
            'fault_geometry': fault_geom,
            'sources': sources,
        }

        # Print source metadata
        metadata = sources[0].get_source_metadata()
        print(f"  Fault: L={metadata['fault_length_km']:.0f} km, "
              f"W={metadata['fault_width_km']:.0f} km, "
              f"D={metadata['fault_top_depth_km']:.1f} km")
        print(f"  Mean/Max slip: {metadata['mean_slip_m']:.2f} / {metadata['max_slip_m']:.2f} m")
        print(f"  Seismic moment: {metadata['seismic_moment_nm']:.3e} N·m")

    return tsunami_sources


def phase4_hazard_aggregation(tsunami_sources, scenario_weights):
    """Phase 4: Compute hazard from tsunami sources."""
    print("\n" + "=" * 70)
    print("PHASE 4: Hazard Aggregation & Exceedance Probability")
    print("=" * 70)

    # Create hazard aggregator
    aggregator = HazardAggregator(scenario_weights)

    # For each magnitude scenario, compute max tsunami heights from sources
    print("\nComputing tsunami hazard metrics...")

    intensity_samples = {}

    for mag, source_dict in tsunami_sources.items():
        sources = source_dict['sources']
        fault_geom = source_dict['fault_geometry']

        print(f"\n  Processing Mw {mag} ({len(sources)} realizations)...")

        # Define observation line perpendicular to fault at various distances
        distances = np.linspace(0, OBSERVATION_RANGE_KM, N_OBSERVATION_POINTS)
        obs_points = np.column_stack([distances, np.zeros_like(distances)])

        # For each realization, compute max displacement
        max_displacements = []
        for source in sources:
            disp = source.get_seafloor_displacement(obs_points)
            max_disp = np.max(np.abs(disp))
            max_displacements.append(max_disp)

        max_displacements = np.array(max_displacements)
        intensity_samples[mag] = max_displacements

        # Add to aggregator
        aggregator.add_scenario_realizations(mag, max_displacements)

        print(f"    Max wave height: {np.mean(max_displacements):.3f} ± "
              f"{np.std(max_displacements):.3f} m")
        print(f"    Range: [{np.min(max_displacements):.2f}, "
              f"{np.max(max_displacements):.2f}] m")

    # Compute aggregated hazard curve
    print("\nAggregating hazard across magnitudes...")
    hazard_curve = aggregator.compute_aggregated_hazard()

    # Print hazard statistics
    print("\nHazard curve at selected wave heights:")
    test_heights = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    for h in test_heights:
        prob = hazard_curve.exceedance_probability_at_level(h)
        ret_period = hazard_curve.return_period_at_level(h)
        if ret_period < 1e6:
            print(f"  h = {h:.1f} m: P(exceed) = {prob:.4e}/year, T ≈ {ret_period:.0f} years")
        else:
            print(f"  h = {h:.1f} m: P(exceed) = {prob:.4e}/year, T > 1e6 years")

    return aggregator, hazard_curve, intensity_samples


def generate_plots(aggregator, hazard_curve, intensity_samples, slip_samples):
    """Generate diagnostic plots."""
    print("\n" + "=" * 70)
    print("GENERATING PLOTS")
    print("=" * 70)

    # Plot 1: Hazard curve
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    heights = np.linspace(0.1, 3.0, 100)
    probs = np.array([hazard_curve.exceedance_probability_at_level(h) for h in heights])
    return_periods = 1.0 / np.maximum(probs, 1e-15)

    ax1.loglog(return_periods, heights, 'b-', linewidth=2, label='Aggregated hazard')
    ax1.set_xlabel('Return Period [years]', fontsize=12)
    ax1.set_ylabel('Wave Height [m]', fontsize=12)
    ax1.set_title('Probabilistic Tsunami Hazard Curve', fontsize=13, fontweight='bold')
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(fontsize=11)

    # Plot 2: Annual exceedance probability
    ax2.loglog(heights, probs, 'r-', linewidth=2, label='Aggregated hazard')
    ax2.set_xlabel('Wave Height [m]', fontsize=12)
    ax2.set_ylabel('Annual Exceedance Probability', fontsize=12)
    ax2.set_title('Exceedance Probability vs. Wave Height', fontsize=13, fontweight='bold')
    ax2.grid(True, which='both', alpha=0.3)
    ax2.legend(fontsize=11)

    plt.tight_layout()
    plot_path = PLOT_DIR / "hazard_curve.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {plot_path}")
    plt.close()

    # Plot 3: Tsunami intensity distributions by magnitude
    fig, axes = plt.subplots(1, len(intensity_samples), figsize=(4*len(intensity_samples), 4))
    if not isinstance(axes, np.ndarray):
        axes = [axes]

    for ax, (mag, intensities) in zip(axes, sorted(intensity_samples.items())):
        ax.hist(intensities, bins=20, alpha=0.7, edgecolor='k', color='steelblue')
        ax.set_xlabel('Max Wave Height [m]', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title(f'Mw {mag}', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        stats_text = f"μ={np.mean(intensities):.2f} m\nσ={np.std(intensities):.2f} m\nn={len(intensities)}"
        ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
               verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
               fontsize=10)

    plt.suptitle('Tsunami Intensity Distributions by Magnitude', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plot_path = PLOT_DIR / "intensity_distributions.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {plot_path}")
    plt.close()

    # Plot 4: Slip sample visualization
    fig, axes = plt.subplots(1, 3, figsize=(12, 3))
    
    for ax, mag in enumerate(sorted(slip_samples.keys())):
        slip = slip_samples[mag]['samples'][0]  # Show first sample
        im = axes[ax].imshow(slip, cmap='RdYlBu_r', aspect='auto')
        axes[ax].set_title(f'Mw {mag} (Sample 0)', fontsize=11, fontweight='bold')
        axes[ax].set_xlabel('Along-strike index', fontsize=10)
        axes[ax].set_ylabel('Down-dip index', fontsize=10)
        plt.colorbar(im, ax=axes[ax], label='Slip [m]')

    plt.suptitle('Stochastic Slip Distributions', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plot_path = PLOT_DIR / "slip_samples.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {plot_path}")
    plt.close()


def print_summary():
    """Print summary of PTHA framework."""
    print("\n" + "=" * 70)
    print("PTHA FRAMEWORK SUMMARY")
    print("=" * 70)
    print("""
This demonstration integrates a complete PTHA workflow:

  [Phase 1: Slip Generation]
  - RQMC sampling with Sobol sequences
  - Owen scrambling for independent realizations
  - Exponential covariance model
  - Moment-scaling via Hanks–Kanamori relation

  [Phase 2: Magnitude-Frequency]
  - Gutenberg–Richter log-linear distribution
  - Annual occurrence rates
  - Scenario weighting

  [Phase 3: Tsunami Sources]
  - Conversion of slip to seafloor displacement
  - Simplified elastic dislocation model
  - Interface for coupling to tsunami solvers

  [Phase 4: Hazard Aggregation]
  - Combination of multiple magnitude scenarios
  - Empirical exceedance probability
  - Return period computation

Key outputs:
  - Hazard curves (exceedance probability vs. wave height)
  - Return periods for engineering design
  - Scenario-specific results for detailed analysis

All modules emphasize:
  ✓ Clarity and reproducibility
  ✓ Defensible methodology
  ✓ Research-grade (not production) code quality
    """)


def main():
    """Execute complete PTHA demonstration."""
    print("\n" + "=" * 70)
    print("PHYSICAL FIDELITY-BASED PTHA FRAMEWORK DEMONSTRATION")
    print("=" * 70)

    # Phase 1: Slip generation
    slip_samples = phase1_rqmc_slip_generation()

    # Phase 2: Magnitude-frequency
    gr, scenario_weights = phase2_magnitude_frequency()

    # Phase 3: Tsunami sources
    tsunami_sources = phase3_tsunami_sources(slip_samples)

    # Phase 4: Hazard aggregation
    aggregator, hazard_curve, intensity_samples = phase4_hazard_aggregation(
        tsunami_sources, scenario_weights
    )

    # Generate plots
    generate_plots(aggregator, hazard_curve, intensity_samples, slip_samples)

    # Print summary
    print_summary()

    print("\n" + "=" * 70)
    print("DEMONSTRATION COMPLETE")
    print("=" * 70)
    print(f"Outputs saved to: {PLOT_DIR}")
    print()


if __name__ == "__main__":
    main()
