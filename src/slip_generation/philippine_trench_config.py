"""
Philippine Trench PTHA Configuration

This module contains domain-specific parameters for the Philippine Trench
Probabilistic Tsunami Hazard Assessment (PTHA) framework based on:

    Thesis: "DEVELOPMENT OF A PHYSICAL FIDELITY-BASED PROBABILISTIC TSUNAMI 
             HAZARD ASSESSMENT FRAMEWORK UTILIZING RANDOMIZED QUASI MONTE CARLO 
             FOR SLIP DISTRIBUTION SAMPLING"
    
    Authors: Henrich Miguel Del Rio Carpio, Carl Vince Callao Dominguez, 
             John Mar Maagad Estimada, Karl Andre Lopez Gutierrez
    Date: December 2025

All parameters are derived from:
  - Heidarzadeh et al. (2025) fault geometry models (Tables 2-3)
  - PHIVOLCS earthquake catalog (1900-2024)
  - Regional seismic studies (Bautista 2019)
  - Global slip inversion studies
  - Bathymetric and exposure data

Physical Fidelity Specifications:
  - Tsunami Simulation: 3rd-order WENO (WENO3) with adaptive mesh refinement (AMR)
  - Slip Model: Matérn covariance with correlation length 20 km, Hurst exponent 0.3
  - Sampling: Randomized Quasi-Monte Carlo (RQMC) with Sobol sequences + Owen scrambling
  - Seismic Recurrence: Gutenberg-Richter model (exponential form, β ≈ 0.9)
  - Validation: Comparison against 2012 and 2023 observed tide gauge records

Study Area: Central to Northern Philippine Trench (7°-12°N)
Focus Location: Dapa coastal zone (inundation depth, building exposure)
Reference Magnitudes: Mw 7.6 (2012, 2023 historical events), Mw 8.0, 8.5, 9.0 (hypothetical)
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple


# ============================================================================
# PHASE 1: SLIP GENERATION PARAMETERS (Section 3.3)
# ============================================================================

@dataclass
class PhilippineTrenchSlipParameters:
    """
    Parameters for stochastic slip generation using RQMC sampling.
    
    References:
      - Section 3.3.2: "Slip Random Field Model"
      - Section 3.3.3: "RQMC Sampling"
      - Section 3.3.1: "RQMC method"
    """
    
    # Gaussian Random Field Parameters (Section 3.3.2)
    covariance_model: str = "matern"  # Options: "exponential", "matern"
    correlation_length_km: float = 20.0  # Global slip inversion standard
    hurst_exponent: float = 0.3  # Roughness/smoothness parameter for Matérn
    
    # RQMC Sampling Configuration (Section 3.3.1)
    sampling_method: str = "rqmc_sobol_owen"  # "rqmc_sobol_owen" or "standard_mc"
    n_subfaults: int = 400  # Approximate grid dimension
    n_samples_per_magnitude: Tuple[int, int] = (2000, 5000)  # Min-max range for ensemble size
    
    # Moment Constraint Parameters (Section 3.3.4)
    moment_error_tolerance: float = 1e-6  # Relative error threshold for moment matching
    moment_scaling_model: str = "hanks_kanamori"  # Moment-magnitude conversion
    
    def __post_init__(self):
        """Validate parameter ranges."""
        assert self.correlation_length_km > 0, "correlation_length_km must be positive"
        assert 0 < self.hurst_exponent < 1, "hurst_exponent must be in (0, 1)"
        assert self.moment_error_tolerance > 0, "moment_error_tolerance must be positive"


# ============================================================================
# PHASE 2: GUTENBERG-RICHTER SEISMICITY PARAMETERS (Section 3.4)
# ============================================================================

@dataclass
class PhilippineTrenchSeismicityParameters:
    """
    Gutenberg-Richter recurrence model parameters for the Philippine Trench.
    
    References:
      - Section 3.4.1: "Selection of the Gutenberg-Richter Recurrence Relation"
      - Section 3.4.2: "Seismicity Parameter Estimation for the Philippine Trench"
      - PHIVOLCS catalog analysis (1900-2024)
    
    Key Findings (Section 3.4.2):
      - Direct MLE on PHIVOLCS catalog (1900-2024): b ≈ 0.95 ± 0.12 (Mw 5.0-7.5 range)
      - Adopted conservative value: β = 0.9 (corresponding to b ≈ 0.39 in log-linear form)
      - Justification: (1) Holocene evidence of M≥8.0 events, (2) Comparison with 
        analogous subduction zones (Nankai, Mediterranean, Java), 
        (3) Sensitivity analysis for conservative hazard estimation
    """
    
    # Gutenberg-Richter Parameters (exponential form: v(m) = β exp(-β m))
    beta: float = 0.9  # Seismicity slope (exponential form)
    b_value: float = 0.39  # Traditional log-linear form (β ≈ 2.303 * b)
    
    # Magnitude Range
    magnitude_minimum: float = 5.0  # Catalog completeness threshold
    magnitude_maximum: float = 9.0  # Physical upper bound (demonstrated Holocene potential)
    
    # Reference Magnitudes (from Heidarzadeh et al. 2025, Table 2-3)
    reference_magnitudes: Dict[str, float] = None
    
    # Annual recurrence rates (derived from beta value)
    annual_rate_at_Mw7p0: float = 0.0  # Estimated from catalog
    
    def __post_init__(self):
        """Initialize reference magnitudes and compute recurrence rates."""
        if self.reference_magnitudes is None:
            self.reference_magnitudes = {
                "Mw7.6_2012": 7.6,
                "Mw7.6_2023": 7.6,
                "Mw8.0_hypothetical": 8.0,
                "Mw8.5_hypothetical": 8.5,
                "Mw9.0_hypothetical": 9.0
            }
        
        # Validate
        assert self.beta > 0, "beta must be positive"
        assert self.magnitude_minimum < self.magnitude_maximum, \
            "magnitude_minimum must be less than magnitude_maximum"


# ============================================================================
# PHASE 3: FAULT GEOMETRY PARAMETERS (Section 3.2.2, Heidarzadeh et al. 2025)
# ============================================================================

@dataclass
class RuptureScenario:
    """Single earthquake rupture scenario with fault geometry."""
    
    name: str  # Scenario identifier
    magnitude: float  # Moment magnitude (Mw)
    length_km: float  # Fault length (L)
    width_km: float  # Fault width (W)
    top_depth_km: float  # Top edge depth
    strike_deg: float  # Strike angle (°)
    dip_deg: float  # Dip angle (°)
    rake_deg: float  # Rake angle (°)
    mean_slip_m: float  # Mean slip on fault (m)
    max_slip_m: float = None  # Maximum slip (for heterogeneity reference)
    event_date: str = None  # For historical events
    
    def __post_init__(self):
        """Validate geometry parameters."""
        assert 5.0 <= self.magnitude <= 9.5, "Magnitude out of physical range"
        assert self.length_km > 0 and self.width_km > 0, "Dimensions must be positive"
        assert 0 <= self.top_depth_km <= 700, "Depth out of subduction range"
        assert 0 <= self.mean_slip_m <= 100, "Slip magnitude unrealistic"


# ============================================================================
# HISTORICAL RUPTURE SCENARIOS (Section 3.2.2, Table 2 from Heidarzadeh et al.)
# ============================================================================

SCENARIO_2012_MW7P6 = RuptureScenario(
    name="2012_Mw7.6_August_31",
    magnitude=7.6,
    length_km=128,
    width_km=90,
    top_depth_km=1,
    strike_deg=345.8,
    dip_deg=45.8,
    rake_deg=61.2,
    mean_slip_m=0.4,
    max_slip_m=3.2,
    event_date="August 31, 2012"
)

SCENARIO_2023_MW7P6 = RuptureScenario(
    name="2023_Mw7.6_December_2",
    magnitude=7.6,
    length_km=144,
    width_km=120,
    top_depth_km=14.5,
    strike_deg=167,
    dip_deg=17,
    rake_deg=62.6,
    mean_slip_m=0.3,
    max_slip_m=1.8,
    event_date="December 2, 2023"
)

# ============================================================================
# HYPOTHETICAL MEGATHRUST SCENARIOS (Section 3.2.2, Table 3 from Heidarzadeh et al.)
# ============================================================================

SCENARIO_MW8P0 = RuptureScenario(
    name="Hypothetical_Mw8.0",
    magnitude=8.0,
    length_km=200,
    width_km=100,
    top_depth_km=7.6,
    strike_deg=160,
    dip_deg=40,
    rake_deg=90,
    mean_slip_m=2.7,
    event_date=None
)

SCENARIO_MW8P5 = RuptureScenario(
    name="Hypothetical_Mw8.5",
    magnitude=8.5,
    length_km=300,
    width_km=100,
    top_depth_km=7.6,
    strike_deg=164,
    dip_deg=39,
    rake_deg=90,
    mean_slip_m=4.7,
    event_date=None
)

SCENARIO_MW9P0 = RuptureScenario(
    name="Hypothetical_Mw9.0",
    magnitude=9.0,
    length_km=600,
    width_km=100,
    top_depth_km=6.1,
    strike_deg=158,
    dip_deg=33,
    rake_deg=90,
    mean_slip_m=8.1,
    event_date=None
)

# Dictionary for convenient access
RUPTURE_SCENARIOS: Dict[str, RuptureScenario] = {
    "2012_Mw7.6": SCENARIO_2012_MW7P6,
    "2023_Mw7.6": SCENARIO_2023_MW7P6,
    "Mw8.0": SCENARIO_MW8P0,
    "Mw8.5": SCENARIO_MW8P5,
    "Mw9.0": SCENARIO_MW9P0,
}


# ============================================================================
# PHASE 3: TSUNAMI SIMULATION PARAMETERS (Section 3.1, Phase 3)
# ============================================================================

@dataclass
class PhilippineTrenchTsunamiParameters:
    """
    High-fidelity tsunami simulation parameters for the Philippine Trench PTHA.
    
    References:
      - Section 3.1 Phase 3: "High-Fidelity Tsunami Simulation"
      - Figures 2.3-2.4: Fully Coupled vs. Superposition methods
      - Section 2.4.1-2.4.2: Validation of fully coupled approach
    
    Physical Fidelity Enhancements:
      - Method 1 (Fully Coupled) vs. Method 4 (Superposition)
      - Key Physics: Non-hydrostatic filtering, acoustic-gravity wave coupling
      - Near-field effects: Rake rotation, shallow fault slip evolution
    """
    
    # Hydrodynamic Solver Configuration
    solver_type: str = "shallow_water_weno"  # Nonlinear shallow-water solver
    weno_order: int = 3  # 3rd-order WENO (WENO3) reconstruction
    time_integration: str = "rk3"  # 3rd-order Runge-Kutta
    
    # Adaptive Mesh Refinement (AMR)
    use_amr: bool = True
    amr_levels: int = 3  # Refinement hierarchy
    base_resolution_m: float = 1000  # 1 km in deep ocean
    coastal_resolution_m: float = 10  # 10 m in Dapa coastal zone
    amr_refine_criterion: str = "gradient_water_elevation"
    
    # Domain Configuration
    domain_extent_lon_deg: Tuple[float, float] = (125.5, 130.5)  # Longitude bounds
    domain_extent_lat_deg: Tuple[float, float] = (6.0, 13.0)    # Latitude bounds
    
    # Bathymetry and Topography
    bathymetry_source: str = "gebco_1km"  # 1 km resolution GEBCO
    include_dynamic_seafloor: bool = True  # Time-dependent seafloor forcing
    
    # Wave Physics Inclusions
    include_nonhydrostatic: bool = False  # Set True for future enhancements
    include_acoustic_coupling: bool = False  # Set True for seismic-acoustic coupling
    
    # Validation Data (Section 3.1, dual-validation strategy)
    validation_tide_gauges: List[str] = None
    
    def __post_init__(self):
        """Initialize validation data."""
        if self.validation_tide_gauges is None:
            # Tide gauges from IOC Sea Level Monitoring Facility + PHIVOLCS
            self.validation_tide_gauges = [
                "Station_Dapa",  # Primary focus (Mw 8.5: 218 building inundation)
                "Station_General_Santos",
                "Station_Zamboanga",
                # Additional stations to be confirmed from thesis
            ]


# ============================================================================
# PHASE 4: HAZARD AGGREGATION PARAMETERS (Section 3.1, Phase 4)
# ============================================================================

@dataclass
class PhilippineTrenchHazardParameters:
    """
    Hazard aggregation and output parameters for Philippine Trench PTHA.
    
    References:
      - Section 3.1 Phase 4: "Hazard Aggregation"
      - Computational outputs: Inundation maps, hazard curves, exposure metrics
    """
    
    # Intensity Metrics
    intensity_metric: str = "inundation_depth_m"  # Primary output
    intensity_thresholds: List[float] = None  # Critical depth levels for analysis
    
    # Return Periods for Hazard Maps
    return_periods: List[int] = None  # Standard return periods in years
    
    # Exposure Analysis
    include_building_exposure: bool = True  # Use OpenStreetMap data
    focus_settlement: str = "Dapa"  # Primary exposure center
    exposure_metric: str = "expected_annual_inundated_buildings"
    
    # Output Products
    output_hazard_curves: bool = True  # Site-specific curves
    output_hazard_maps: bool = True  # Spatial exceedance probability maps
    output_slip_samples: bool = True  # Slip distribution realizations (Phase 1)
    output_inundation_fields: bool = True  # Detailed time series
    
    def __post_init__(self):
        """Initialize default output specifications."""
        if self.intensity_thresholds is None:
            self.intensity_thresholds = [0.5, 1.0, 2.0, 5.0, 10.0, 15.0]  # meters
        
        if self.return_periods is None:
            self.return_periods = [100, 500, 2500]  # years


# ============================================================================
# INTEGRATED CONFIGURATION CLASS
# ============================================================================

@dataclass
class PhilippineTrenchPTHAConfiguration:
    """
    Complete Philippine Trench PTHA configuration integrating all phases.
    
    This configuration class brings together all domain-specific parameters
    for the 4-phase PTHA workflow:
      1. Stochastic Slip Generation (RQMC with Matérn covariance)
      2. Earthquake Probability Model (Gutenberg-Richter)
      3. High-Fidelity Tsunami Simulation (WENO3 + AMR)
      4. Hazard Aggregation (Probabilistic synthesis)
    
    Validation Strategy:
      - Statistical: RQMC vs. Monte Carlo convergence comparison
      - Physical: Simulated vs. observed tide gauge records (2012, 2023)
      - Comparison: Against published hazard curves for other subduction zones
    
    Study Specifications:
      - Study Area: Central to Northern Philippine Trench (7°-12°N)
      - Focus: Dapa coastal zone (building exposure, inundation metrics)
      - Scenarios: 5 reference magnitudes (Mw 7.6 [2x], 8.0, 8.5, 9.0)
      - Ensemble Size: 2000-5000 realizations per magnitude
      - Outputs: Hazard maps, curves, exposure metrics, slip samples
    """
    
    # Phase 1: Slip Generation
    slip_params: PhilippineTrenchSlipParameters = None
    
    # Phase 2: Seismicity Model
    seismicity_params: PhilippineTrenchSeismicityParameters = None
    
    # Phase 3: Tsunami Simulation
    tsunami_params: PhilippineTrenchTsunamiParameters = None
    
    # Phase 4: Hazard Aggregation
    hazard_params: PhilippineTrenchHazardParameters = None
    
    # Study Metadata
    study_name: str = "Philippine_Trench_PTHA_2025"
    study_area: str = "Central to Northern Philippine Trench (7°-12°N)"
    focus_location: str = "Dapa"
    reference_citation: str = "Thesis: DEVELOPMENT OF A PHYSICAL FIDELITY-BASED PTHA"
    authors: List[str] = None
    
    def __post_init__(self):
        """Initialize all parameter groups with defaults."""
        if self.slip_params is None:
            self.slip_params = PhilippineTrenchSlipParameters()
        
        if self.seismicity_params is None:
            self.seismicity_params = PhilippineTrenchSeismicityParameters()
        
        if self.tsunami_params is None:
            self.tsunami_params = PhilippineTrenchTsunamiParameters()
        
        if self.hazard_params is None:
            self.hazard_params = PhilippineTrenchHazardParameters()
        
        if self.authors is None:
            self.authors = [
                "Henrich Miguel Del Rio Carpio",
                "Carl Vince Callao Dominguez",
                "John Mar Maagad Estimada",
                "Karl Andre Lopez Gutierrez"
            ]
    
    def summary(self) -> str:
        """Generate human-readable configuration summary."""
        return f"""
Philippine Trench PTHA Configuration Summary
{'='*60}

STUDY: {self.study_name}
AREA: {self.study_area}
FOCUS: {self.focus_location}

PHASE 1 - SLIP GENERATION:
  Covariance Model: {self.slip_params.covariance_model}
  Correlation Length: {self.slip_params.correlation_length_km} km
  Hurst Exponent: {self.slip_params.hurst_exponent}
  Sampling: {self.slip_params.sampling_method}
  Samples per Magnitude: {self.slip_params.n_samples_per_magnitude[0]}-{self.slip_params.n_samples_per_magnitude[1]}

PHASE 2 - SEISMICITY:
  Beta (Exponential): {self.seismicity_params.beta}
  b-value (Log-linear): {self.seismicity_params.b_value:.2f}
  Magnitude Range: {self.seismicity_params.magnitude_minimum}-{self.seismicity_params.magnitude_maximum}

PHASE 3 - TSUNAMI SIMULATION:
  Solver: {self.tsunami_params.solver_type}
  WENO Order: {self.tsunami_params.weno_order}
  AMR Levels: {self.tsunami_params.amr_levels}
  Base Resolution: {self.tsunami_params.base_resolution_m} m
  Coastal Resolution: {self.tsunami_params.coastal_resolution_m} m

PHASE 4 - HAZARD OUTPUTS:
  Intensity Metric: {self.hazard_params.intensity_metric}
  Return Periods: {self.hazard_params.return_periods}
  Thresholds: {self.hazard_params.intensity_thresholds}

RUPTURE SCENARIOS:
  {len(RUPTURE_SCENARIOS)} scenarios defined (see RUPTURE_SCENARIOS dict)

AUTHORS: {', '.join(self.authors)}
"""


# ============================================================================
# DEFAULT CONFIGURATION INSTANCE
# ============================================================================

PHILIPPINE_TRENCH_PTHA_CONFIG = PhilippineTrenchPTHAConfiguration()


if __name__ == "__main__":
    # Print configuration summary for verification
    print(PHILIPPINE_TRENCH_PTHA_CONFIG.summary())
    
    # List all available scenarios
    print("\nAVAILABLE RUPTURE SCENARIOS:")
    print("="*60)
    for name, scenario in RUPTURE_SCENARIOS.items():
        print(f"{name}: {scenario.magnitude} | L={scenario.length_km} km, "
              f"W={scenario.width_km} km | Mean slip={scenario.mean_slip_m} m")
