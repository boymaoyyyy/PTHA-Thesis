"""
PTHA Framework Configuration Template
=====================================

Copy this file and modify parameters for your specific analysis.

Example usage:
    from ptha_config import PTHAConfig
    
    config = PTHAConfig.from_file('my_config.py')
    generator = StochasticSlipGenerator(**config.slip_generation_params())
"""

class PTHAConfig:
    """Configuration container for PTHA framework parameters."""

    # ====================================================================
    # PHASE 1: STOCHASTIC SLIP GENERATION
    # ====================================================================

    # Target magnitude(s) to simulate
    SLIP_MAGNITUDES = [7.5, 8.0, 8.5]

    # Number of RQMC samples per magnitude
    SLIP_SAMPLES_PER_MAGNITUDE = 50

    # Subfault discretization size [km]
    # Smaller = finer resolution but more samples needed
    # Typical range: 10-30 km
    SUBFAULT_SIZE_KM = 20.0

    # Spatial covariance model
    # Options: 'exponential' (default), 'matern'
    COVARIANCE_MODEL = 'exponential'

    # Correlation length [km]
    # Represents spatial extent of slip heterogeneity
    # Typical range: 15-40 km
    CORRELATION_LENGTH_KM = 20.0

    # Matérn smoothness parameter (only if COVARIANCE_MODEL='matern')
    # Typical range: 0.3-1.5
    # Lower values → rougher slip; higher → smoother
    MATERN_SMOOTHNESS = 0.5

    # Shear modulus (rigidity) [Pa]
    # Typical values:
    #   - 4.0e10 = 40 GPa (standard for subduction zones)
    #   - 3.0e10 = 30 GPa (softer zones)
    #   - 5.0e10 = 50 GPa (stiffer zones)
    RIGIDITY_PA = 4.0e10

    # Random seed for reproducibility
    # Set to fixed value for reproducible runs
    # Set to None for truly random behavior
    SLIP_SEED = 42

    # ====================================================================
    # PHASE 2: MAGNITUDE-FREQUENCY DISTRIBUTION
    # ====================================================================

    # Gutenberg-Richter 'a' parameter
    # Represents activity rate: log₁₀(N) = a - b·M
    # Typical range: 4.0-6.0 (varies by region)
    # Higher a → more earthquakes
    GR_A_VALUE = 5.0

    # Gutenberg-Richter 'b' parameter (slope)
    # Typical range: 0.8-1.2
    # Common approximation: b ≈ 1.0
    # Interpretation: b=1.0 means ~10x fewer earthquakes per +1 magnitude
    GR_B_VALUE = 1.0

    # Minimum magnitude of interest
    # Earthquakes below this are not considered
    GR_M_MIN = 7.5

    # Maximum possible magnitude
    # Practical upper bound for truncated GR distribution
    GR_M_MAX = 9.5

    # ====================================================================
    # PHASE 3: TSUNAMI SOURCE (FAULT GEOMETRY)
    # ====================================================================

    # Fault geometry source
    # Options: 'csv', 'manual'
    FAULT_GEOMETRY_SOURCE = 'csv'

    # Path to fault parameters CSV file
    # Expected columns: "Hypothetical scenario", "L (km)", "W (km)", 
    #                   "Top depth (km)", "Mean strike (°)", "Mean dip (°)", 
    #                   "Rake (°)", "Slip (m)"
    FAULT_PARAMS_CSV = 'data/fault_geometry/heidarzadeh_2025_table3.csv'

    # Manual fault geometry (used if FAULT_GEOMETRY_SOURCE='manual')
    # Uncomment and modify to override CSV loading
    """
    FAULT_GEOMETRY_MANUAL = {
        8.5: {
            'length_km': 300.0,
            'width_km': 100.0,
            'top_depth_km': 7.6,
            'strike_deg': 164.0,
            'dip_deg': 39.0,
            'rake_deg': 90.0,
            'slip_mean_m': 4.7,
        }
    }
    """

    # ====================================================================
    # PHASE 4: HAZARD AGGREGATION
    # ====================================================================

    # Intensity levels at which to compute hazard
    # If None, auto-generated from data
    # Units: meters of wave height
    INTENSITY_LEVELS = None  # Auto-generate
    # INTENSITY_LEVELS = np.linspace(0.1, 3.0, 50)  # Manual specification

    # Observation region for tsunami computation
    # Maximum distance from fault [km]
    OBSERVATION_RANGE_KM = 200.0

    # Number of observation points
    # Higher = more accurate but slower
    N_OBSERVATION_POINTS = 300

    # ====================================================================
    # OUTPUT & VISUALIZATION
    # ====================================================================

    # Output directory
    OUTPUT_DIR = 'output'

    # Plot directory
    PLOT_DIR = 'output/plots'

    # Generate diagnostic plots
    GENERATE_PLOTS = True

    # Plot resolution [DPI]
    PLOT_DPI = 150

    # ====================================================================
    # NUMERICAL PARAMETERS
    # ====================================================================

    # Covariance matrix regularization (nugget)
    # Added to diagonal for numerical stability
    COVARIANCE_NUGGET = 1e-10

    # Moment constraint relative error tolerance
    # Acceptable tolerance: slip_error / target_moment
    MOMENT_ERROR_TOL = 1e-6

    # Log-log interpolation points for hazard curves
    HAZARD_CURVE_POINTS = 100

    # ====================================================================
    # CLASS METHODS
    # ====================================================================

    @classmethod
    def slip_generation_params(cls):
        """Return parameters dict for StochasticSlipGenerator."""
        return {
            'subfault_size_km': cls.SUBFAULT_SIZE_KM,
            'covariance_model': cls.COVARIANCE_MODEL,
            'correlation_length_km': cls.CORRELATION_LENGTH_KM,
            'rigidity_pa': cls.RIGIDITY_PA,
        }

    @classmethod
    def magnitude_frequency_params(cls):
        """Return parameters dict for GutenbergRichterRelation."""
        return {
            'a_value': cls.GR_A_VALUE,
            'b_value': cls.GR_B_VALUE,
            'M_min': cls.GR_M_MIN,
            'M_max': cls.GR_M_MAX,
        }

    @classmethod
    def to_dict(cls):
        """Export all configuration parameters as dictionary."""
        return {
            'slip': {
                'magnitudes': cls.SLIP_MAGNITUDES,
                'samples_per_magnitude': cls.SLIP_SAMPLES_PER_MAGNITUDE,
                'subfault_size_km': cls.SUBFAULT_SIZE_KM,
                'covariance_model': cls.COVARIANCE_MODEL,
                'correlation_length_km': cls.CORRELATION_LENGTH_KM,
                'matern_smoothness': cls.MATERN_SMOOTHNESS,
                'rigidity_pa': cls.RIGIDITY_PA,
                'seed': cls.SLIP_SEED,
            },
            'magnitude_frequency': {
                'a_value': cls.GR_A_VALUE,
                'b_value': cls.GR_B_VALUE,
                'M_min': cls.GR_M_MIN,
                'M_max': cls.GR_M_MAX,
            },
            'tsunami_source': {
                'geometry_source': cls.FAULT_GEOMETRY_SOURCE,
                'params_csv': cls.FAULT_PARAMS_CSV,
            },
            'hazard_aggregation': {
                'intensity_levels': cls.INTENSITY_LEVELS,
                'observation_range_km': cls.OBSERVATION_RANGE_KM,
                'n_observation_points': cls.N_OBSERVATION_POINTS,
            },
            'output': {
                'output_dir': cls.OUTPUT_DIR,
                'plot_dir': cls.PLOT_DIR,
                'generate_plots': cls.GENERATE_PLOTS,
                'plot_dpi': cls.PLOT_DPI,
            },
        }

    @classmethod
    def print_summary(cls):
        """Print configuration summary."""
        print("\nPTHA Framework Configuration")
        print("=" * 60)
        print(f"\nPhase 1: Slip Generation")
        print(f"  Magnitudes: {cls.SLIP_MAGNITUDES}")
        print(f"  Samples/magnitude: {cls.SLIP_SAMPLES_PER_MAGNITUDE}")
        print(f"  Subfault size: {cls.SUBFAULT_SIZE_KM} km")
        print(f"  Covariance: {cls.COVARIANCE_MODEL} (ξ={cls.CORRELATION_LENGTH_KM} km)")
        print(f"  Rigidity: {cls.RIGIDITY_PA:.1e} Pa")

        print(f"\nPhase 2: Magnitude-Frequency")
        print(f"  GR: log(N) = {cls.GR_A_VALUE} - {cls.GR_B_VALUE}·M")
        print(f"  Range: [{cls.GR_M_MIN}, {cls.GR_M_MAX}]")

        print(f"\nPhase 3: Tsunami Source")
        print(f"  Fault geometry: {cls.FAULT_PARAMS_CSV}")

        print(f"\nPhase 4: Hazard Aggregation")
        print(f"  Observation range: {cls.OBSERVATION_RANGE_KM} km")
        print(f"  Observation points: {cls.N_OBSERVATION_POINTS}")

        print(f"\nOutput")
        print(f"  Directory: {cls.OUTPUT_DIR}")
        print(f"  Plots: {cls.GENERATE_PLOTS}")
        print("=" * 60 + "\n")


# ========================================================================
# PRESET CONFIGURATIONS
# ========================================================================

class QuickTestConfig(PTHAConfig):
    """Fast configuration for testing (coarse resolution)."""
    SLIP_MAGNITUDES = [8.5]
    SLIP_SAMPLES_PER_MAGNITUDE = 10
    SUBFAULT_SIZE_KM = 30.0
    N_OBSERVATION_POINTS = 50


class ResearchConfig(PTHAConfig):
    """Research-grade configuration (balanced accuracy)."""
    SLIP_MAGNITUDES = [7.5, 8.0, 8.5]
    SLIP_SAMPLES_PER_MAGNITUDE = 100
    SUBFAULT_SIZE_KM = 20.0
    GR_A_VALUE = 5.2
    GR_B_VALUE = 1.0
    N_OBSERVATION_POINTS = 300


class ProductionConfig(PTHAConfig):
    """High-accuracy configuration (requires more computation)."""
    SLIP_MAGNITUDES = [7.5, 8.0, 8.5, 9.0]
    SLIP_SAMPLES_PER_MAGNITUDE = 500
    SUBFAULT_SIZE_KM = 15.0
    CORRELATION_LENGTH_KM = 25.0
    GR_A_VALUE = 5.1
    GR_B_VALUE = 1.05
    N_OBSERVATION_POINTS = 500
    HAZARD_CURVE_POINTS = 200


# ========================================================================
# EXAMPLE CUSTOMIZATIONS
# ========================================================================

"""
EXAMPLE 1: Japan Trench Subduction Zone
────────────────────────────────────────

class JapanTrenchConfig(PTHAConfig):
    SLIP_MAGNITUDES = [7.5, 8.0, 8.5, 9.0, 9.5]
    GR_A_VALUE = 5.4
    GR_B_VALUE = 0.95
    SUBFAULT_SIZE_KM = 20.0
    CORRELATION_LENGTH_KM = 30.0


EXAMPLE 2: Cascadia Subduction Zone
────────────────────────────────────

class CascadiaConfig(PTHAConfig):
    SLIP_MAGNITUDES = [8.0, 8.5, 9.0]
    GR_A_VALUE = 4.8
    GR_B_VALUE = 1.05
    SUBFAULT_SIZE_KM = 25.0
    CORRELATION_LENGTH_KM = 35.0
    RIGIDITY_PA = 3.5e10


EXAMPLE 3: Mediterranean Subduction
────────────────────────────────────

class MediterraneanConfig(PTHAConfig):
    SLIP_MAGNITUDES = [7.0, 7.5, 8.0, 8.5]
    GR_A_VALUE = 4.5
    GR_B_VALUE = 1.1
    SUBFAULT_SIZE_KM = 15.0
    CORRELATION_LENGTH_KM = 20.0
"""


if __name__ == "__main__":
    # Print default configuration
    PTHAConfig.print_summary()

    # Demonstrate preset configurations
    print("\nAvailable Presets:")
    print("  • PTHAConfig (default)")
    print("  • QuickTestConfig (fast)")
    print("  • ResearchConfig (balanced)")
    print("  • ProductionConfig (accurate)")

    print("\nUsage Example:")
    print("  from ptha_config import ResearchConfig")
    print("  config = ResearchConfig()")
    print("  config.print_summary()")
