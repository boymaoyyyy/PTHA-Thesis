"""
PTHA Framework: Physical Fidelity–Based Probabilistic Tsunami Hazard Assessment

A modular, research-grade Python framework for generating probabilistic tsunami
hazard curves from stochastic earthquake rupture models.

Main components:
  - covariance: Spatial correlation models
  - slip_sampler: RQMC slip generation
  - magnitude_frequency: Gutenberg-Richter relationships
  - tsunami_source: Earthquake-to-seafloor displacement
  - hazard_aggregation: Hazard curve computation

Example:
    from slip_sampler import StochasticSlipGenerator
    from magnitude_frequency import GutenbergRichterRelation
    from tsunami_source import SimpleDislocationSource
    from hazard_aggregation import HazardAggregator
    
    # Phase 1: Generate slip
    gen = StochasticSlipGenerator(magnitude=8.5, ...)
    slip = gen.generate_rqmc_sample(0)
    
    # Phase 2: Set up magnitude-frequency
    gr = GutenbergRichterRelation(a_value=5.0, b_value=1.0)
    
    # Phase 3: Tsunami source
    source = SimpleDislocationSource(fault_geometry, slip)
    
    # Phase 4: Hazard aggregation
    aggregator = HazardAggregator({8.5: 0.001})
    aggregator.add_scenario_realizations(8.5, wave_heights)
    hazard = aggregator.compute_aggregated_hazard()

"""

__version__ = "0.1.0"
__author__ = "Undergraduate Thesis"
__description__ = "Physical Fidelity-Based PTHA Framework"

# Convenience imports
try:
    from .covariance import (
        exponential_covariance,
        matern_covariance,
        build_distance_matrix,
        build_covariance_matrix,
    )
    from .magnitude_frequency import (
        GutenbergRichterRelation,
        HazardScenarioWeights,
        hanks_kanamori_moment,
        hanks_kanamori_magnitude,
    )
    from .slip_sampler import (
        StochasticSlipGenerator,
        generate_slip_sample,
        load_fault_geometry,
    )
    from .tsunami_source import (
        TsunamiSource,
        SimpleDislocationSource,
        TsunamiPropagationInterface,
    )
    from .hazard_aggregation import (
        HazardCurve,
        HazardAggregator,
        compute_empirical_pdf_kde,
    )
except ImportError:
    # Allow partial imports if modules not all available
    pass

__all__ = [
    # Covariance
    'exponential_covariance',
    'matern_covariance',
    'build_distance_matrix',
    'build_covariance_matrix',
    # Magnitude-frequency
    'GutenbergRichterRelation',
    'HazardScenarioWeights',
    'hanks_kanamori_moment',
    'hanks_kanamori_magnitude',
    # Slip generation
    'StochasticSlipGenerator',
    'generate_slip_sample',
    'load_fault_geometry',
    # Tsunami source
    'TsunamiSource',
    'SimpleDislocationSource',
    'TsunamiPropagationInterface',
    # Hazard aggregation
    'HazardCurve',
    'HazardAggregator',
    'compute_empirical_pdf_kde',
]
