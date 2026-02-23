"""
INDEX: PTHA Framework Files & Documentation
==========================================

Navigate this documentation to understand and use the PTHA framework.

📚 DOCUMENTATION FILES (Start Here)
===================================

1. QUICK_START.md ⭐ START HERE
   Quick reference and practical examples for each component
   • Installation instructions
   • 4-phase overview with code examples
   • Common issues & solutions
   • Tips & best practices

2. README_PTHA.md (Comprehensive Reference)
   Complete methodology documentation
   • Full mathematical formulations
   • Literature references for each method
   • Module-by-module deep dive
   • Extension guidelines

3. IMPLEMENTATION_SUMMARY.md (What Was Built)
   Overview of the entire implementation
   • Module descriptions and key classes
   • Feature list and architecture
   • Usage workflow examples
   • Next steps for your thesis

4. This file (INDEX.md)
   Navigation guide and file reference

⚙️ SOURCE CODE MODULES (src/slip_generation/)
==============================================

Core Implementation:

1. __init__.py
   Package initialization and convenience imports
   • Exposes all public APIs
   • Enables: from slip_generation import *

2. covariance.py (PHASE 1A: Spatial Correlation)
   Implements spatial covariance models
   
   Key Classes:
   • (None - functional design)
   
   Key Functions:
   • exponential_covariance(distance_matrix, correlation_length)
   • matern_covariance(distance_matrix, correlation_length, smoothness)
   • build_covariance_matrix(coordinates, model, **kwargs)
   • build_distance_matrix(coordinates)
   
   Usage:
     from covariance import build_covariance_matrix
     C = build_covariance_matrix(coords, model='exponential', 
                                 correlation_length=20.0)

3. slip_sampler.py (PHASE 1B: RQMC Slip Generation) ⭐ MAIN
   Generates stochastic slip via Sobol sequences + Owen scrambling
   
   Key Classes:
   • StochasticSlipGenerator
     - generate_rqmc_sample(sample_index, seed)
     - generate_ensemble(n_samples, seed)
     - get_metadata()
   
   Key Functions:
   • generate_slip_sample(magnitude, sample_index, ...)
   • load_fault_geometry(magnitude, fault_params_csv)
   
   Usage:
     from slip_sampler import StochasticSlipGenerator
     gen = StochasticSlipGenerator(magnitude=8.5, ...)
     slip = gen.generate_rqmc_sample(0)
     samples = gen.generate_ensemble(100)

4. magnitude_frequency.py (PHASE 2: Seismic Hazard Model)
   Gutenberg-Richter distribution and scenario weighting
   
   Key Classes:
   • GutenbergRichterRelation
     - annual_rate(magnitude)
     - probability_given_scenario(magnitude)
     - sample_magnitude(n_samples, rng)
   
   • HazardScenarioWeights
     - weight_magnitude_bin(M_lower, M_upper)
     - normalize_weights(magnitude_bins)
   
   Key Functions:
   • hanks_kanamori_moment(magnitude)
   • hanks_kanamori_magnitude(moment)
   
   Usage:
     from magnitude_frequency import GutenbergRichterRelation
     gr = GutenbergRichterRelation(a_value=5.0, b_value=1.0)
     rate = gr.annual_rate(8.5)

5. tsunami_source.py (PHASE 3: Earthquake-to-Tsunami Conversion)
   Elastic dislocation model for seafloor displacement
   
   Key Classes:
   • TsunamiSource (ABC - abstract base)
   • SimpleDislocationSource(TsunamiSource)
     - get_seafloor_displacement(observation_points)
     - get_source_metadata()
   
   • TsunamiPropagationInterface
     - export_source_to_netcdf(output_path, ...)
     - get_initial_condition_on_grid(lon_range, lat_range, ...)
   
   Usage:
     from tsunami_source import SimpleDislocationSource
     source = SimpleDislocationSource(fault_geom, slip)
     disp = source.get_seafloor_displacement(obs_points)

6. hazard_aggregation.py (PHASE 4: Hazard Curve Generation)
   Exceedance probability and return period computation
   
   Key Classes:
   • HazardCurve
     - exceedance_probability_at_level(intensity_level)
     - return_period_at_level(intensity_level)
     - intensity_at_return_period(return_period)
   
   • HazardAggregator
     - add_scenario_realizations(magnitude, intensities)
     - compute_aggregated_hazard(intensity_levels)
     - get_summary_statistics()
   
   Key Functions:
   • compute_empirical_pdf_kde(samples, support_range, n_points)
   
   Usage:
     from hazard_aggregation import HazardAggregator
     agg = HazardAggregator(scenario_weights)
     agg.add_scenario_realizations(8.5, heights)
     hazard = agg.compute_aggregated_hazard()

Demonstration & Testing:

7. ptha_demo.py
   Complete end-to-end demonstration of all 4 phases
   
   Key Functions:
   • phase1_rqmc_slip_generation()
   • phase2_magnitude_frequency()
   • phase3_tsunami_sources(slip_samples)
   • phase4_hazard_aggregation(tsunami_sources, scenario_weights)
   • generate_plots(aggregator, hazard_curve, ...)
   • main()
   
   Run: python ptha_demo.py
   
   Outputs:
     • output/plots/hazard_curve.png
     • output/plots/intensity_distributions.png
     • output/plots/slip_samples.png
     • Console statistics

8. test_framework.py
   Validation test suite (6 test categories)
   
   Tests:
   1. Moment Constraint Verification
   2. Covariance Matrix Properties
   3. Magnitude-Frequency Distribution
   4. Hazard Curve Properties
   5. Tsunami Source Physics
   6. Hazard Aggregation Logic
   
   Run: python test_framework.py
   
   Output: Pass/Fail report for each test

Configuration:

9. ptha_config.py
   Configuration template and presets
   
   Classes:
   • PTHAConfig (default)
   • QuickTestConfig (fast testing)
   • ResearchConfig (balanced accuracy)
   • ProductionConfig (high accuracy)
   
   Usage:
     from ptha_config import ResearchConfig
     config = ResearchConfig()
     config.print_summary()

Legacy/Utilities:

10. plot_slip.py
    Visualization utilities for slip distributions
    (Pre-existing, integrated with framework)

📊 DATA FILES (data/fault_geometry/)
===================================

1. heidarzadeh_2025_table2.csv
   Alternative fault scenario data
   Columns: Event date, L, W, Top depth, Strike, Dip, Rake, Max slip, Mean slip

2. heidarzadeh_2025_table3.csv ⭐ PRIMARY
   Main fault geometry for this thesis
   Scenarios: Mw 8.0, Mw 8.5
   Used by: slip_sampler.py, ptha_demo.py

📁 OUTPUT FILES (output/)
========================

Generated by running ptha_demo.py:

slip_samples/
  └── slip_Mw*.npy          Saved slip realizations (NumPy format)

plots/
  ├── hazard_curve.png      ⭐ Main PTHA result
  ├── intensity_distributions.png
  └── slip_samples.png

🚀 QUICK START WORKFLOW
======================

Step 1: Understand the framework
  → Read QUICK_START.md (5 minutes)

Step 2: Run the demonstration
  → cd src/slip_generation
  → python ptha_demo.py
  → Check output/plots/

Step 3: Explore individual components
  → Read relevant module docstrings
  → Try code examples from QUICK_START.md

Step 4: Run validation tests
  → python test_framework.py
  → Verify all tests pass

Step 5: Customize for your analysis
  → Edit ptha_config.py or parameters in ptha_demo.py
  → Regenerate with custom settings

Step 6: Integrate into your thesis
  → Extract key results from hazard curves
  → Document methodology with references
  → Cite PTHA literature appropriately

📖 EXAMPLE USE CASES
===================

Use Case 1: Generate slip samples for Mw 8.5
  Code:
    from slip_sampler import StochasticSlipGenerator
    gen = StochasticSlipGenerator(magnitude=8.5, ...)
    samples = gen.generate_ensemble(n_samples=100, seed=42)
  Reference: QUICK_START.md section "Component 1"

Use Case 2: Compute hazard curve from scenarios
  Code:
    from hazard_aggregation import HazardAggregator
    aggregator = HazardAggregator(scenario_weights)
    hazard = aggregator.compute_aggregated_hazard()
    prob = hazard.exceedance_probability_at_level(1.5)
  Reference: QUICK_START.md section "Component 4"

Use Case 3: Test different covariance models
  Code:
    # Compare exponential vs. Matérn
    for model in ['exponential', 'matern']:
        gen = StochasticSlipGenerator(..., covariance_model=model)
  Reference: covariance.py docstring

Use Case 4: Validate moment conservation
  Code:
    python test_framework.py
  Reference: test_framework.py

🔍 KEY CONCEPTS
===============

1. RQMC (Randomized Quasi-Monte Carlo)
   • Systematic sampling of parameter space
   • Sobol sequences: low-discrepancy
   • Owen scrambling: statistical independence
   • Better convergence than standard Monte Carlo

2. Covariance Models
   • Exponential: C(r) = exp(-r/ξ) — simple, fast
   • Matérn: Generalized with smoothness parameter ν
   • Both implemented, switchable

3. Gutenberg-Richter Distribution
   • log₁₀(N) = a - b·M
   • Defines earthquake frequency vs. magnitude
   • a, b values region-specific (~5.0, ~1.0 typical)

4. Moment Scaling
   • Slip generated independently, then scaled
   • Ensures total seismic moment matches Mw target
   • Error: < 1e-6 relative tolerance

5. Tsunami Hazard
   • P(h ≥ height) computed for each scenario magnitude
   • Aggregated with magnitude-frequency weights
   • Return period T = 1 / P (per year)

📋 PARAMETER REFERENCE
======================

Critical Parameters:

SUBFAULT_SIZE_KM
  • Discretization of fault surface [km]
  • Smaller → more realistic but slower
  • Recommended: 20 km (balanced)
  • Range: 10-30 km

CORRELATION_LENGTH_KM
  • Spatial correlation length of slip [km]
  • Smaller → more heterogeneous slip
  • Recommended: 20 km
  • Range: 15-40 km (geophysically motivated)

GR_A_VALUE
  • Gutenberg-Richter intercept
  • Higher → more earthquakes
  • Recommended: 5.0 (subduction zones)
  • Range: 4.0-6.0 (region-dependent)

GR_B_VALUE
  • Gutenberg-Richter slope
  • Typically ≈ 1.0 (±0.2)
  • Interpretation: 10^(-b) = ratio of M+1 to M

See ptha_config.py for all parameters with ranges and defaults.

✅ VALIDATION CHECKLIST
=======================

Before using results in your thesis:

□ Run test_framework.py (all tests pass)
□ Verify moment error < 1e-6 in output
□ Check hazard curve monotonicity in plots
□ Validate against published PTHA studies
□ Document all parameter choices
□ Include methodology references
□ Save configuration and seeds for reproducibility

🎓 FOR YOUR THESIS
==================

1. Methodology Section
   • Cite "README_PTHA.md" methodology section
   • Reference original papers listed there
   • Include mathematical formulations

2. Results Section
   • Include hazard_curve.png
   • Report key return periods (e.g., 100-yr, 500-yr)
   • Compare with published studies

3. Appendix (Optional)
   • Include slip_samples.png examples
   • Document configuration choices
   • List key parameters used

4. Code Availability
   • Provide access to this repository
   • Note: Recommend citing as "PTHA Framework v0.1.0"
   • Include commit hash or version tag

📞 TROUBLESHOOTING
==================

Problem: "Module not found" errors
Solution: Ensure you're in src/slip_generation/ directory
          Check Python path includes the directory

Problem: Covariance matrix singular
Solution: Increase nugget parameter (default: 1e-10)
          Use finer/coarser subfault discretization

Problem: Plots not generating
Solution: Check matplotlib backend (matplotlib.use('Agg'))
          Verify output/plots/ directory exists

Problem: Moment constraint violated
Solution: Check subfault size consistency
          Verify fault geometry matches CSV
          Try different random seed

See QUICK_START.md "Tips & Troubleshooting" for more.

📚 REFERENCES
=============

Key Literature (Comprehensive List in README_PTHA.md):

1. Sobol (1967) - Low-discrepancy sequences
2. Hanks & Kanamori (1979) - Moment magnitude
3. Gutenberg & Richter (1954) - Magnitude-frequency
4. Mai & Beroza (2002) - Spatial slip models
5. Satake (2007) - Tsunami science
6. Goda et al. (2016) - PTHA handbook

🏁 NEXT STEPS
=============

1. Run: python ptha_demo.py
2. Read: QUICK_START.md + README_PTHA.md
3. Modify: ptha_config.py for your scenarios
4. Test: python test_framework.py
5. Integrate into thesis with proper citations

Questions? Check the relevant docstring or README section.
"""


if __name__ == "__main__":
    print(__doc__)
