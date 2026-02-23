"""
QUICK START GUIDE: PTHA Framework
==================================

This guide shows how to use each component of the PTHA framework.

====================================================================
INSTALLATION
====================================================================

Requirements:
  - Python 3.8+
  - numpy, scipy, pandas, matplotlib

Install:
  pip install numpy scipy pandas matplotlib

====================================================================
COMPONENT 1: SLIP GENERATION (RQMC)
====================================================================

Generate stochastic earthquake slip distributions:

    from slip_sampler import StochasticSlipGenerator
    from pathlib import Path
    
    # Set up generator for Mw 8.5
    generator = StochasticSlipGenerator(
        magnitude=8.5,
        fault_params_csv='data/fault_geometry/heidarzadeh_2025_table3.csv',
        subfault_size_km=20.0,
        covariance_model='exponential',
        correlation_length_km=20.0,
    )
    
    # Generate a single RQMC sample
    slip_sample_0 = generator.generate_rqmc_sample(sample_index=0, seed=42)
    
    # Or generate an ensemble
    samples = generator.generate_ensemble(n_samples=100, seed=42)
    
    # Access metadata
    metadata = generator.get_metadata()
    print(metadata['target_moment_nm'])

Output:
  slip_sample_0 : ndarray of shape (n_down, n_along) in meters
  samples : list of 100 such arrays

Key parameters:
  - magnitude: Target Mw (e.g., 8.5)
  - covariance_model: 'exponential' or 'matern'
  - correlation_length_km: Spatial correlation length (typically 15-30 km)
  - subfault_size_km: Grid discretization (typically 10-30 km)

====================================================================
COMPONENT 2: MAGNITUDE-FREQUENCY (GUTENBERG-RICHTER)
====================================================================

Set up earthquake magnitude-frequency relationships:

    from magnitude_frequency import (
        GutenbergRichterRelation, 
        HazardScenarioWeights,
        hanks_kanamori_moment,
    )
    
    # Create GR relation (a=5.0, b=1.0 typical for subduction zones)
    gr = GutenbergRichterRelation(
        a_value=5.0, 
        b_value=1.0, 
        M_min=7.5, 
        M_max=9.5
    )
    
    # Query annual occurrence rate
    rate = gr.annual_rate(8.5)  # events/year
    
    # Get conditional probability given M_min
    prob = gr.probability_given_scenario(8.5)
    
    # Sample magnitudes
    mag_samples = gr.sample_magnitude(n_samples=1000, rng=rng)
    
    # Compute scenario weights for bin [7.5, 8.0)
    weighter = HazardScenarioWeights(gr)
    result = weighter.normalize_weights([7.5, 8.0, 8.5, 9.0, 9.5])
    print(result['bin_weights'])  # Annual rates
    print(result['bin_probabilities'])  # Normalized probabilities
    
    # Moment-magnitude conversion
    M0 = hanks_kanamori_moment(8.5)  # → seismic moment in N·m
    Mw = hanks_kanamori_magnitude(M0)  # → back to magnitude

Output:
  rate : float (events/year)
  prob : float in [0, 1]
  mag_samples : ndarray of magnitudes
  result : dict with 'bin_weights', 'bin_probabilities', etc.

====================================================================
COMPONENT 3: TSUNAMI SOURCES (DISLOCATION)
====================================================================

Convert slip distributions to seafloor displacement:

    from tsunami_source import SimpleDislocationSource
    import numpy as np
    
    # Define fault geometry
    fault_geom = {
        'length_km': 300.0,
        'width_km': 100.0,
        'top_depth_km': 7.6,
        'strike_deg': 164.0,
        'dip_deg': 39.0,
        'rake_deg': 90.0,
        'n_subfaults_along': 5,
        'n_subfaults_down': 3,
    }
    
    # Use slip from Phase 1
    slip = slip_sample_0  # or any realization
    
    # Create tsunami source
    source = SimpleDislocationSource(
        fault_geometry=fault_geom,
        slip_distribution=slip,
        rigidity=4.0e10,  # Pa
    )
    
    # Compute seafloor displacement at observation points
    obs_points = np.array([
        [100, 0],    # 100 km along-strike
        [0, 50],     # 50 km across-strike
        [150, 100],  # Diagonal
    ])
    
    displacement = source.get_seafloor_displacement(obs_points)
    
    # Get source metadata
    metadata = source.get_source_metadata()
    print(metadata['max_slip_m'])
    print(metadata['seismic_moment_nm'])

Output:
  displacement : ndarray of vertical displacements [m] at each observation point
  metadata : dict with source parameters

====================================================================
COMPONENT 4: HAZARD AGGREGATION
====================================================================

Combine slip realizations and magnitude scenarios into hazard curves:

    from hazard_aggregation import HazardAggregator
    
    # Define scenario weights (from Phase 2)
    scenario_weights = {
        8.0: 0.01,    # 0.01 events/year
        8.5: 0.001,   # 0.001 events/year
        9.0: 0.0001,  # 0.0001 events/year
    }
    
    # Create aggregator
    aggregator = HazardAggregator(scenario_weights)
    
    # For each magnitude scenario, add tsunami intensity samples
    # (e.g., from Phase 3 computations)
    aggregator.add_scenario_realizations(8.0, wave_heights_m_8p0)
    aggregator.add_scenario_realizations(8.5, wave_heights_m_8p5)
    aggregator.add_scenario_realizations(9.0, wave_heights_m_9p0)
    
    # Compute aggregated hazard curve
    hazard_curve = aggregator.compute_aggregated_hazard()
    
    # Query exceedance probability at specific heights
    prob_1m = hazard_curve.exceedance_probability_at_level(1.0)  # P(h >= 1 m)
    
    # Query return periods
    return_period_1m = hazard_curve.return_period_at_level(1.0)  # years
    
    # Reverse query: intensity at return period
    height_500yr = hazard_curve.intensity_at_return_period(500)  # m for T=500 years
    
    # Get scenario-specific statistics
    stats = aggregator.get_summary_statistics()

Output:
  prob_1m : float (annual probability, e.g., 1e-3)
  return_period_1m : float (years)
  height_500yr : float (meters)
  stats : dict with per-magnitude statistics

====================================================================
COMPLETE WORKFLOW EXAMPLE
====================================================================

Here's a minimal example connecting all phases:

    import numpy as np
    from pathlib import Path
    from slip_sampler import StochasticSlipGenerator
    from magnitude_frequency import GutenbergRichterRelation
    from tsunami_source import SimpleDislocationSource
    from hazard_aggregation import HazardAggregator
    
    DATA_DIR = Path('data/fault_geometry')
    
    # Phase 1: Generate slips
    magnitudes = [8.0, 8.5]
    slip_samples = {}
    for mag in magnitudes:
        gen = StochasticSlipGenerator(
            magnitude=mag,
            fault_params_csv=DATA_DIR / 'heidarzadeh_2025_table3.csv',
        )
        slip_samples[mag] = gen.generate_ensemble(n_samples=50, seed=42)
    
    # Phase 2: Magnitude-frequency
    gr = GutenbergRichterRelation(a_value=5.0, b_value=1.0)
    scenario_weights = {mag: gr.annual_rate(mag) for mag in magnitudes}
    
    # Phase 3 & 4: Tsunami sources and aggregation
    aggregator = HazardAggregator(scenario_weights)
    
    fault_geom = {  # Define based on your fault model
        'length_km': 300.0,
        'width_km': 100.0,
        'top_depth_km': 7.6,
        'strike_deg': 164.0,
        'dip_deg': 39.0,
        'rake_deg': 90.0,
        'n_subfaults_along': 5,
        'n_subfaults_down': 3,
    }
    
    for mag in magnitudes:
        max_heights = []
        for slip in slip_samples[mag]:
            source = SimpleDislocationSource(fault_geom, slip)
            # Compute at observation points
            obs_pts = np.array([[100, 0], [150, 100]])
            disp = source.get_seafloor_displacement(obs_pts)
            max_heights.append(np.max(np.abs(disp)))
        
        aggregator.add_scenario_realizations(mag, np.array(max_heights))
    
    # Get hazard curve
    hazard = aggregator.compute_aggregated_hazard()
    
    # Query results
    print(hazard.exceedance_probability_at_level(1.5))  # P(h >= 1.5m)
    print(hazard.return_period_at_level(1.5))  # Return period for 1.5m

====================================================================
TESTING & VALIDATION
====================================================================

Run the validation test suite:

    python test_framework.py

This verifies:
  ✓ Moment constraints satisfied
  ✓ Covariance matrices positive-definite
  ✓ Magnitude-frequency distribution correct
  ✓ Hazard curves monotonic
  ✓ Tsunami source physics valid
  ✓ Hazard aggregation consistent

====================================================================
RUNNING THE FULL DEMO
====================================================================

Execute the complete 4-phase demonstration:

    python ptha_demo.py

This generates:
  - plots/hazard_curve.png (main PTHA result)
  - plots/intensity_distributions.png (by magnitude)
  - plots/slip_samples.png (example slip patterns)
  - Console output with detailed statistics

====================================================================
TIPS & BEST PRACTICES
====================================================================

1. Random seeds:
   Always use fixed seeds for reproducibility:
   
   samples = generator.generate_ensemble(n_samples=100, seed=42)

2. Grid resolution:
   Balance accuracy vs. computation:
   - Fine (10 km): ~900 subfaults, slower
   - Medium (20 km): ~225 subfaults, good compromise
   - Coarse (30 km): ~50 subfaults, fast

3. Ensemble size:
   - 50–100 samples: Quick prototyping
   - 500–1000 samples: Research quality
   - 10,000+ samples: Publication quality

4. Correlation length:
   Physically reasonable values: 15–40 km
   Smaller → rougher slip patterns
   Larger → smoother slip patterns

5. Covariance model:
   - exponential: Default, fast, realistic
   - matern: More flexible, can tune smoothness

6. Magnitude range:
   Include scenario magnitudes that contribute significantly
   to hazard (typically M_min to M_max - 0.5)

====================================================================
COMMON ISSUES & SOLUTIONS
====================================================================

Q: Moment constraint error > 1e-5?
A: Check subfault size consistency, ensure fault geometry matches
   CSV parameters exactly.

Q: Covariance matrix not positive-definite?
A: Increase nugget parameter in covariance functions or use
   finer grid discretization.

Q: Hazard curve not monotonic?
A: Increase number of slip samples per scenario (currently trying
   to estimate distribution from too few samples).

Q: Import errors?
A: Ensure all modules are in slip_generation/ directory and
   you're running from correct working directory.

====================================================================
REFERENCES & FURTHER READING
====================================================================

See README_PTHA.md for:
  - Complete methodology references
  - Mathematical formulations
  - Extension guidelines
  - Full module documentation

Key papers:
  1. Mai & Beroza (2002): Spatial slip heterogeneity
  2. Sobol (1967): Low-discrepancy sequences
  3. Goda & Atkinson (2010): Spatially-correlated hazard

====================================================================
GETTING HELP
====================================================================

1. Check test_framework.py for working examples
2. Review docstrings in each module
3. Run ptha_demo.py for complete workflow
4. Consult README_PTHA.md for theory and methods
5. Check fault_geometry CSV for available scenarios
"""


def print_quick_reference():
    """Print quick reference table."""
    print("""
╔════════════════════════════════════════════════════════════════════╗
║                    PTHA Framework Quick Reference                  ║
╠════════════════════════════════════════════════════════════════════╣
║ COMPONENT        │ CLASS/FUNCTION        │ KEY OUTPUT              ║
╠──────────────────┼──────────────────────┼─────────────────────────╣
║ Phase 1: Slip    │ StochasticSlipGen    │ Slip arrays [m]         ║
║ Generation       │ .generate_ensemble() │ (n_down, n_along)       ║
├──────────────────┼──────────────────────┼─────────────────────────┤
║ Phase 2: Mag-    │ GutenbergRichter     │ Annual rates [evt/yr]   ║
║ Frequency        │ .annual_rate(M)      │ Scenario weights        ║
├──────────────────┼──────────────────────┼─────────────────────────┤
║ Phase 3: Tsunami │ SimpleDislocation    │ Displacement [m]        ║
║ Source           │ .get_seafloor_...()  │ At observation points   ║
├──────────────────┼──────────────────────┼─────────────────────────┤
║ Phase 4: Hazard  │ HazardAggregator     │ Hazard curves           ║
║ Aggregation      │ .compute_aggregated_ │ Return periods [years]  ║
║                  │ hazard()             │ Exceedance probabilities║
╚════════════════════════════════════════════════════════════════════╝
    """)


if __name__ == "__main__":
    print_quick_reference()
