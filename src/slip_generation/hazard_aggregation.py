"""
Probabilistic tsunami hazard aggregation and exceedance probability estimation.

This module combines stochastic slip realizations, magnitude-frequency distributions,
and tsunami response metrics to produce hazard curves showing the probability of
exceeding specified tsunami wave heights as a function of return period.

Standard PTHA workflow:
    1. For each magnitude scenario M_i with weight w_i:
    2. Generate multiple stochastic slip realizations
    3. Compute tsunami response (max height) for each realization
    4. Aggregate: P(Wave Height ≥ h) = ∑_i w_i * P(h | M_i)
    5. Convert to return period: T = 1 / P(h)

References:
    - Rikitake & Aida (1988): Tsunami hazard probability in Japan
    - Geist & Parsons (2006): Probabilistic analysis of the seismic hazard
    - Goda (2016): Importance of rupture variations in seismic hazard
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde


class HazardCurve:
    """
    Represents a probabilistic hazard curve: P(intensity ≥ level) vs. intensity level.

    A hazard curve quantifies the relationship between ground motion intensity
    (in this case, tsunami wave height) and annual exceedance probability.
    """

    def __init__(self, intensity_levels, exceedance_probabilities):
        """
        Initialize hazard curve.

        Args:
            intensity_levels (ndarray): Tsunami wave heights or other intensity metric [m].
                                       Should be sorted in ascending order.
            exceedance_probabilities (ndarray): Annual probability of exceeding each level.
                                               Same length as intensity_levels.
        """
        # Sort by intensity level
        order = np.argsort(intensity_levels)
        self.intensity = intensity_levels[order]
        self.probability = exceedance_probabilities[order]

    def exceedance_probability_at_level(self, intensity_level):
        """
        Compute exceedance probability at a specified intensity level.

        Uses linear interpolation on log-log scale (common in hazard analysis).

        Args:
            intensity_level (float): Wave height [m].

        Returns:
            float: Annual exceedance probability [0, 1].
        """
        # Log-log interpolation
        log_intensity = np.log(self.intensity)
        log_prob = np.log(self.probability + 1e-15)
        
        interp_fn = interp1d(log_intensity, log_prob, kind='linear', 
                            fill_value='extrapolate', bounds_error=False)
        
        log_result = interp_fn(np.log(intensity_level))
        return np.exp(log_result)

    def return_period_at_level(self, intensity_level):
        """
        Compute return period (years) at a specified intensity level.

        Return period T (years) = 1 / P(exceeding level per year)

        Args:
            intensity_level (float): Wave height [m].

        Returns:
            float: Return period [years]. Returns inf if probability is very small.
        """
        prob = self.exceedance_probability_at_level(intensity_level)
        if prob < 1e-15:
            return np.inf
        return 1.0 / prob

    def intensity_at_return_period(self, return_period):
        """
        Compute intensity level corresponding to a specified return period.

        Args:
            return_period (float): Return period [years].

        Returns:
            float: Intensity level (wave height) [m].
        """
        target_prob = 1.0 / return_period
        
        # Find intensity corresponding to this probability
        # Use linear interpolation on log-log scale
        log_intensity = np.log(self.intensity)
        log_prob = np.log(self.probability + 1e-15)
        
        interp_fn = interp1d(log_prob, log_intensity, kind='linear',
                            fill_value='extrapolate', bounds_error=False)
        
        log_result = interp_fn(np.log(target_prob))
        return np.exp(log_result)

    def to_dict(self):
        """Export hazard curve as dictionary."""
        return {
            'intensity': self.intensity.copy(),
            'probability': self.probability.copy(),
        }


class HazardAggregator:
    """
    Aggregate tsunami hazard from multiple magnitude scenarios and slip realizations.

    This class implements standard PTHA aggregation:

        P(h) = ∑_j w_j ∑_k p_k(h | s_{jk})

    where:
    - j indexes magnitude scenarios
    - k indexes slip realizations for scenario j
    - w_j is the weight (annual rate) of scenario j
    - s_{jk} is the slip realization
    - p_k(h | s) is the probability of wave height ≥ h given slip s
    """

    def __init__(self, magnitude_scenario_weights, intensity_bins=None):
        """
        Initialize aggregator.

        Args:
            magnitude_scenario_weights (dict): Maps magnitude to annual rate (weight).
                Keys: magnitude values (float)
                Values: annual rates [events/year]

            intensity_bins (ndarray, optional): Intensity levels at which to compute
                                              hazard. If None, will be inferred from data.
        """
        self.scenario_weights = magnitude_scenario_weights
        self.intensity_bins = intensity_bins
        
        # Storage for scenario-specific results
        self.scenario_intensities = {}  # mag -> list of intensities per realization
        self.scenario_labels = []  # For tracking

    def add_scenario_realizations(self, magnitude, intensity_values):
        """
        Add tsunami intensity values from multiple realizations of a single magnitude.

        Args:
            magnitude (float): Scenario magnitude.
            intensity_values (ndarray): Array of max wave heights [m] from each realization.
                                       Shape (n_realizations,).
        """
        if magnitude not in self.scenario_intensities:
            self.scenario_intensities[magnitude] = []
        
        self.scenario_intensities[magnitude].append(intensity_values)

    def compute_aggregated_hazard(self, intensity_levels=None):
        """
        Compute aggregated hazard curve from all scenarios and realizations.

        Method:
            1. For each magnitude scenario, compute empirical CDF of intensities
            2. Weight each scenario by its annual rate
            3. Sum to get total probability

        Args:
            intensity_levels (ndarray, optional): Levels at which to compute hazard.
                                                If None, use provided bins or auto-generate.

        Returns:
            HazardCurve: The aggregated hazard curve object.
        """
        if intensity_levels is None:
            if self.intensity_bins is not None:
                intensity_levels = self.intensity_bins
            else:
                # Auto-generate from data
                all_intensities = []
                for mag_list in self.scenario_intensities.values():
                    for realization_array in mag_list:
                        all_intensities.extend(realization_array)
                
                all_intensities = np.array(all_intensities)
                intensity_levels = np.linspace(np.min(all_intensities) * 0.5,
                                              np.max(all_intensities) * 1.5,
                                              50)
        
        # Compute total exceedance probability
        total_probability = np.zeros_like(intensity_levels, dtype=float)
        
        for magnitude, realization_list in self.scenario_intensities.items():
            scenario_weight = self.scenario_weights.get(magnitude, 0.0)
            
            if scenario_weight <= 0:
                continue
            
            # Combine all realizations for this magnitude
            all_intensities = np.concatenate(realization_list)
            
            # Empirical CDF: P(intensity ≥ h) for each h
            for i, h_level in enumerate(intensity_levels):
                # Count realizations with intensity ≥ h_level
                exceedance_frac = np.mean(all_intensities >= h_level)
                
                # Add to total, weighted by scenario rate
                total_probability[i] += scenario_weight * exceedance_frac
        
        return HazardCurve(intensity_levels, total_probability)

    def compute_scenario_hazard_curve(self, magnitude):
        """
        Compute hazard curve for a single magnitude scenario.

        Args:
            magnitude (float): Scenario magnitude.

        Returns:
            HazardCurve: Hazard curve for this scenario, or None if no data.
        """
        if magnitude not in self.scenario_intensities:
            return None
        
        realization_list = self.scenario_intensities[magnitude]
        all_intensities = np.concatenate(realization_list)
        
        # Auto-generate intensity levels
        h_min, h_max = np.min(all_intensities), np.max(all_intensities)
        intensity_levels = np.linspace(h_min * 0.5, h_max * 1.5, 30)
        
        # Empirical exceedance probabilities
        probs = np.array([np.mean(all_intensities >= h) for h in intensity_levels])
        
        return HazardCurve(intensity_levels, probs)

    def get_summary_statistics(self):
        """
        Compute summary statistics for all scenario realizations.

        Returns:
            dict: Statistics for each magnitude scenario.
        """
        stats = {}
        for mag, realization_list in self.scenario_intensities.items():
            all_int = np.concatenate(realization_list)
            stats[mag] = {
                'n_realizations': len(all_int),
                'mean_intensity': np.mean(all_int),
                'median_intensity': np.median(all_int),
                'std_intensity': np.std(all_int),
                'min_intensity': np.min(all_int),
                'max_intensity': np.max(all_int),
                'weight': self.scenario_weights.get(mag, 0.0),
            }
        return stats


def compute_empirical_pdf_kde(samples, support_range=None, n_points=100):
    """
    Estimate probability density function from samples using kernel density estimation.

    Useful for understanding the distribution of tsunami heights from
    stochastic slip realizations.

    Args:
        samples (ndarray): Data samples (e.g., max wave heights from realizations).
        support_range (tuple, optional): (min, max) range for evaluation.
                                        If None, use data range ± 10%.
        n_points (int, optional): Number of points for PDF evaluation. Default: 100.

    Returns:
        dict: Dictionary with keys:
            'support': x-axis values
            'pdf': probability density function values
            'kde': the fitted KDE object
    """
    if support_range is None:
        s_min, s_max = np.min(samples), np.max(samples)
        margin = 0.1 * (s_max - s_min)
        support_range = (s_min - margin, s_max + margin)
    
    support = np.linspace(support_range[0], support_range[1], n_points)
    
    # Fit KDE
    kde = gaussian_kde(samples)
    pdf = kde(support)
    
    return {
        'support': support,
        'pdf': pdf,
        'kde': kde,
    }


if __name__ == "__main__":
    print("=" * 60)
    print("PTHA Hazard Aggregation Demonstration")
    print("=" * 60)
    
    # Scenario weights (annual rates)
    scenario_weights = {
        7.5: 1.0e-1,  # 0.1 events/year
        8.0: 1.0e-2,  # 0.01 events/year
        8.5: 1.0e-3,  # 0.001 events/year
    }
    
    aggregator = HazardAggregator(scenario_weights)
    
    # Add synthetic realizations for each scenario
    np.random.seed(42)
    
    for mag in [7.5, 8.0, 8.5]:
        # Synthetic tsunami heights (larger for larger magnitude)
        base_height = mag - 6.0
        n_realizations = 100
        
        # Gamma-like distribution (realistic for max wave heights)
        heights = np.random.gamma(shape=2.0, scale=base_height, size=n_realizations)
        aggregator.add_scenario_realizations(mag, heights)
        print(f"Magnitude {mag}: {n_realizations} realizations, "
              f"mean height = {heights.mean():.2f} m, max = {heights.max():.2f} m")
    
    print()
    
    # Compute aggregated hazard
    hazard_curve = aggregator.compute_aggregated_hazard()
    
    print("Aggregated hazard curve (probabilities at selected intensities):")
    test_heights = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
    for h in test_heights:
        prob = hazard_curve.exceedance_probability_at_level(h)
        ret_period = hazard_curve.return_period_at_level(h)
        print(f"  h = {h:.1f} m: P(exceed) = {prob:.4e}/year, Return period = {ret_period:.1f} years")
    
    print()
    print("Wave heights at selected return periods:")
    return_periods = np.array([10, 50, 100, 500])
    for tp in return_periods:
        h = hazard_curve.intensity_at_return_period(tp)
        print(f"  T = {tp:3d} years: h = {h:.2f} m")
    
    print()
    print("Scenario statistics:")
    stats = aggregator.get_summary_statistics()
    for mag, stat in stats.items():
        print(f"  Magnitude {mag}:")
        for key, val in stat.items():
            print(f"    {key}: {val}")
