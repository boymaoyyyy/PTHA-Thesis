"""
README: Physical Fidelity–Based Probabilistic Tsunami Hazard Assessment (PTHA) Framework

This package implements a research-grade PTHA methodology suitable for
undergraduate and graduate research. The framework emphasizes clarity,
reproducibility, and mathematical defensibility over computational optimization.

====================================================================
OVERVIEW
====================================================================

Probabilistic Tsunami Hazard Assessment (PTHA) quantifies the probability of
exceeding specified tsunami wave heights as a function of return period. This
framework implements a complete workflow from earthquake slip modeling through
hazard curve generation.

Main Phases:
  Phase 1: Stochastic Slip Generation (RQMC with Sobol sequences)
  Phase 2: Magnitude-Frequency Relationships (Gutenberg–Richter)
  Phase 3: Tsunami Source Modeling (Elastic dislocation)
  Phase 4: Hazard Aggregation (Exceedance probability)

====================================================================
MODULE STRUCTURE
====================================================================

slip_generation/
├── covariance.py
│   └── Spatial covariance models (exponential, Matérn)
│
├── slip_sampler.py
│   └── RQMC slip generation using Sobol sequences + Owen scrambling
│
├── magnitude_frequency.py
│   └── Gutenberg–Richter distribution and scenario weighting
│
├── tsunami_source.py
│   └── Earthquake slip → seafloor displacement conversion
│
├── hazard_aggregation.py
│   └── Hazard curve computation and exceedance probability
│
└── ptha_demo.py
    └── Complete end-to-end demonstration

data/
└── fault_geometry/
    └── Heidarzadeh (2025) subduction zone parameters

output/
├── slip_samples/
│   └── Generated slip realizations (NPY format)
└── plots/
    └── Diagnostic visualizations

====================================================================
PHASE 1: STOCHASTIC SLIP GENERATION
====================================================================

Objective:
  Generate physically-realistic earthquake slip distributions that:
  • Reproduce observed spatial heterogeneity
  • Satisfy moment-constraint matching target magnitude
  • Are computationally efficient
  • Allow reproducible ensemble sampling

Key Features:
  1. RQMC Sampling:
     - Uses Sobol sequences (low-discrepancy, systematic coverage)
     - Owen scrambling enables independent samples
     - Superior convergence vs. standard Monte Carlo
  
  2. Spatial Correlation:
     - Exponential or Matérn covariance models
     - Realistic fault slip autocorrelation (typically 10–30 km)
     - Built via Cholesky decomposition
  
  3. Moment Scaling:
     - Hanks–Kanamori relation: Mw ↔ M₀
     - Slip scaled post-generation to match target moment
     - Preserves spatial pattern while fixing total seismic moment

Example Usage:
  from slip_sampler import StochasticSlipGenerator
  
  generator = StochasticSlipGenerator(
      magnitude=8.5,
      fault_params_csv='data/fault_geometry/heidarzadeh_2025_table3.csv',
      subfault_size_km=20.0,
      covariance_model='exponential',
      correlation_length_km=20.0,
  )
  
  # Generate ensemble
  samples = generator.generate_ensemble(n_samples=100, seed=42)
  
  # Each sample is shape (n_down, n_along) in meters

References:
  - Sobol (1967): "On the distribution of points in a cube"
  - Owen (1998): "Scrambling Sobol and Niederreiter-Xing points"
  - Mai & Beroza (2002): "A spatial random field model for modulating slip"

====================================================================
PHASE 2: MAGNITUDE-FREQUENCY DISTRIBUTIONS
====================================================================

Objective:
  Characterize earthquake frequency-magnitude relationships and compute
  scenario weights for hazard aggregation.

Key Features:
  1. Gutenberg–Richter Relation:
     log₁₀(N) = a - b·M
     
     where N = annual count of earthquakes with magnitude ≥ M.
  
  2. Scenario Weighting:
     Weight of bin [M_lower, M_upper) = N(M_lower) - N(M_upper)
  
  3. Moment-Magnitude Conversion:
     Hanks–Kanamori: Mw = (log₁₀(M₀) - 4.4) / 1.5
     where M₀ is seismic moment [N·m]

Example Usage:
  from magnitude_frequency import GutenbergRichterRelation
  
  gr = GutenbergRichterRelation(a_value=5.0, b_value=1.0, 
                                M_min=7.0, M_max=9.5)
  
  # Annual rate
  rate = gr.annual_rate(8.5)  # events/year
  
  # Scenario probability
  prob = gr.probability_given_scenario(8.5)

References:
  - Gutenberg & Richter (1954): "Seismicity of the Earth"
  - Hanks & Kanamori (1979): "A moment magnitude scale"
  - Cornell (1968): "Engineering seismic risk analysis"

====================================================================
PHASE 3: TSUNAMI SOURCE CHARACTERIZATION
====================================================================

Objective:
  Convert earthquake slip distributions into seafloor vertical
  displacement (the initial condition for tsunami propagation).

Key Features:
  1. Elastic Dislocation Model:
     Simplified Green's function approach
     ζ(x,y) ≈ (μ/π) * Σ_i s_i * A_i / r_i²
     
     where:
     • ζ = vertical seafloor displacement [m]
     • s_i = slip on subfault i [m]
     • A_i = subfault area [m²]
     • r_i = distance to observation point [m]
  
  2. Interface for External Solvers:
     Provides structured API for coupling to GeoClaw, COMCOT, etc.
  
  3. Metadata Export:
     Supports NetCDF export for standard tsunami solvers

Example Usage:
  from tsunami_source import SimpleDislocationSource
  
  source = SimpleDislocationSource(fault_geometry, slip_distribution)
  
  # Compute displacement at observation points
  obs_points = [[100, 0], [0, 100], [-50, 50]]  # km
  displacement = source.get_seafloor_displacement(obs_points)

References:
  - Okada (1985): "Surface deformation due to shear and tensile faults"
  - Satake (2007): "Tsunamis: Case studies and lessons learned"

====================================================================
PHASE 4: HAZARD AGGREGATION & CURVE GENERATION
====================================================================

Objective:
  Combine slip realizations, magnitude scenarios, and tsunami response
  to produce probabilistic hazard curves.

Methodology:
  For each wave height h:
    P(H ≥ h) = Σ_j w_j · P(H ≥ h | M_j)
  
  where:
  • j indexes magnitude scenarios
  • w_j = annual rate of scenario j
  • P(H ≥ h | M_j) = empirical probability from realizations

Return Period Conversion:
  T [years] = 1 / P(H ≥ h [per year])

Example Usage:
  from hazard_aggregation import HazardAggregator
  
  aggregator = HazardAggregator(scenario_weights={8.5: 0.001, 8.0: 0.01})
  
  # Add tsunami intensity samples for each magnitude
  aggregator.add_scenario_realizations(8.5, max_heights_m)
  
  # Compute hazard curve
  hazard = aggregator.compute_aggregated_hazard()
  
  # Query exceedance probability
  prob = hazard.exceedance_probability_at_level(1.5)  # for 1.5 m wave
  
  # Query return period
  ret_period = hazard.return_period_at_level(1.5)  # years

References:
  - Rikitake & Aida (1988): "Tsunami hazard probability in Japan"
  - Geist & Parsons (2006): "Probabilistic analysis of seismic hazard"
  - Goda (2016): "Importance of rupture variations in seismic hazard"

====================================================================
RUNNING THE DEMONSTRATION
====================================================================

Prerequisites:
  - Python 3.8+
  - NumPy, SciPy, Pandas, Matplotlib
  
  Install: pip install numpy scipy pandas matplotlib

Execute:
  cd src/slip_generation
  python ptha_demo.py

This runs the complete 4-phase workflow and generates:
  • hazard_curve.png: Main PTHA result
  • intensity_distributions.png: Scenario-specific results
  • slip_samples.png: Example slip distributions
  • Console output with detailed statistics

====================================================================
MATHEMATICAL DETAILS
====================================================================

1. RQMC Slip Generation

   Sobol Sequence:
   • n-dimensional low-discrepancy point set
   • Systematic coverage of [0,1)^n with low overlap
   • Owen scrambling: Apply random bit-reversal permutation
   
   Algorithm:
   1. Generate Sobol sample u ~ [0,1)^N (N = number of subfaults)
   2. Map to normal: ξ = Φ⁻¹(u)
   3. Apply covariance: log_s = log(s̄) + L·ξ  [L = Cholesky factor]
   4. Exponentiate: s_raw = exp(log_s)
   5. Scale moment: s_final = s_raw · (M₀_target / M₀_raw)

2. Covariance Models

   Exponential:
   C(r) = exp(-r / ξ)
   
   Matérn (ν = smoothness):
   C(r) = (2^(1-ν) / Γ(ν)) · (√(2ν)r/ξ)^ν · K_ν(√(2ν)r/ξ)

3. Gutenberg–Richter Distribution

   Magnitude sampling (inverse transform):
   M ~ M_min + (M_max - M_min) · (1 - u^(1/(1-10^(-b(M_max-M_min))))) / (1-10^(-b(M_max-M_min)))

4. Hazard Curve Computation

   Log-log interpolation (common in engineering practice):
   log(P(h)) = interp(log(h); data_points)

====================================================================
CODE QUALITY & RESEARCH STANDARDS
====================================================================

This code is designed for:
  ✓ Clarity: Each module has clear purpose and interfaces
  ✓ Reproducibility: Fixed random seeds, explicit parameters
  ✓ Defensibility: Extensive docstrings with methodology references
  ✓ Extensibility: Modular design allows easy modifications
  ✓ Correctness: Moment constraints validated, error checks included

NOT optimized for:
  ✗ Production performance (use purpose-built codes like COMCOT)
  ✗ Very large grids (feasible for ≲500 subfaults)
  ✗ Parallel computing (sequential architecture)

Validation Practices:
  • Moment scaling: rel_err < 1e-6 (error bounded)
  • Covariance matrix: Positive-definite, verified via eigenvalues
  • Hazard aggregation: Probability conservation checked

====================================================================
EXTENDING THE FRAMEWORK
====================================================================

To add new features:

1. Custom covariance model:
   In covariance.py, add function (e.g., von_karman_cov)
   then register in build_covariance_matrix()

2. Alternative magnitude-frequency:
   Extend magnitude_frequency.py with new class inheriting from interface

3. Coupling to external solver:
   Use TsunamiPropagationInterface.export_source_to_netcdf()
   to write seafloor displacement in standard format

4. More realistic tsunami:
   Replace SimpleDislocationSource with full 3D Green's function
   (e.g., call to Okada2003 library)

====================================================================
KEY REFERENCES
====================================================================

Core Methodology:
  1. Goda & Atkinson (2010): "Probabilistic characterization of spatially 
     correlated response spectra for earthquakes in Japan"
  2. Stirling et al. (2012): "National seismic hazard model for New Zealand"

Slip Models:
  3. Mai & Beroza (2002): "A spatial random field model for modulating slip"
  4. Herrero & Bernard (1994): "A kinematic self-similar rupture process"

RQMC Methods:
  5. Dick & Pillichshammer (2010): "Digital Nets and Sequences"
  6. Caflisch (1998): "Monte Carlo and quasi-Monte Carlo methods"

Tsunami Science:
  7. Satake (2007): "Tsunamis: Case studies and lessons learned"
  8. Rikitake & Aida (1988): "Tsunami hazard probability in Japan"

====================================================================
CONTACT & LICENSING
====================================================================

For questions or improvements, contact the thesis advisor.

This code is provided for educational and research purposes.
Modifications should maintain full documentation of changes.

Last updated: February 2026
"""
