"""
Magnitude-frequency relationships and seismic hazard probability models.

This module implements the Gutenberg–Richter magnitude-frequency distribution
and related tools for PTHA, including:
    - Annual occurrence rates as a function of magnitude
    - Scenario weights for hazard aggregation
    - Hanks–Kanamori moment-magnitude conversion

References:
    - Gutenberg & Richter (1954): Seismicity of the Earth
    - Hanks & Kanamori (1979): A moment magnitude scale
    - Cornell (1968): Engineering seismic risk analysis
"""

import numpy as np


class GutenbergRichterRelation:
    """
    Gutenberg–Richter magnitude-frequency distribution.

    The log-linear Gutenberg–Richter relation:
        log₁₀(N) = a - b·M
    
    where N is the annual count of earthquakes with magnitude ≥ M,
    a is the activity rate intercept, and b is the slope (typically 0.8–1.2).

    This is the standard model for characterizing seismic activity in PTHA.
    """

    def __init__(self, a_value, b_value, M_min=4.0, M_max=9.5):
        """
        Initialize Gutenberg–Richter relation.

        Args:
            a_value (float): Activity rate parameter (log₁₀ scale).
                            Represents log₁₀(N) for M=0.
            b_value (float): Slope parameter (dimensionless, typically ~1.0).
            M_min (float, optional): Minimum magnitude of completeness. Default: 4.0.
            M_max (float, optional): Maximum possible magnitude. Default: 9.5.
        """
        self.a = a_value
        self.b = b_value
        self.M_min = M_min
        self.M_max = M_max

    def annual_rate(self, magnitude):
        """
        Compute annual occurrence rate for a given magnitude.

        Args:
            magnitude (float or ndarray): Earthquake magnitude(s).

        Returns:
            float or ndarray: Annual occurrence rate(s) for M ≥ magnitude.

        Notes:
            N(M) = 10^(a - b·M) [events per year]
            Applies Gutenberg–Richter below M_max and truncates above.
        """
        if isinstance(magnitude, np.ndarray):
            rate = np.where(magnitude <= self.M_max, 
                           10.0**(self.a - self.b * magnitude),
                           0.0)
            return rate
        else:
            if magnitude <= self.M_max:
                return 10.0**(self.a - self.b * magnitude)
            else:
                return 0.0

    def probability_given_scenario(self, magnitude):
        """
        Conditional probability of scenario magnitude M given occurrence ≥ M_min.

        This is used to weight scenarios when aggregating hazard.
        
        P(M | M ≥ M_min) = (annual rate of M) / (annual rate of M_min)

        Args:
            magnitude (float): Earthquake magnitude.

        Returns:
            float: Conditional probability [0, 1].
        """
        rate_scenario = self.annual_rate(magnitude)
        rate_min = self.annual_rate(self.M_min)
        
        if rate_min == 0:
            return 0.0
        return rate_scenario / rate_min

    def sample_magnitude(self, n_samples, rng=None):
        """
        Draw random magnitude samples from truncated Gutenberg–Richter distribution.

        Uses inverse transform sampling: M = M_min + (M_max - M_min) * u / (1 - u·(1 - b^(M_max - M_min)))
        where u ~ Uniform(0, 1) and b' = 10^(-b).

        Args:
            n_samples (int): Number of magnitude samples.
            rng (np.random.Generator, optional): Random number generator. 
                                               Default: uses global random state.

        Returns:
            ndarray: Array of magnitude samples, shape (n_samples,).
        """
        if rng is None:
            rng = np.random.default_rng()
        
        u = rng.uniform(0, 1, size=n_samples)
        
        # Compute b' = 10^(-b)
        b_prime = 10.0**(-self.b)
        
        # Inverse transform for truncated distribution
        numerator = 1.0 - u * (1.0 - b_prime**(self.M_max - self.M_min))
        M_samples = self.M_min - np.log10(numerator) / self.b
        
        return np.clip(M_samples, self.M_min, self.M_max)


def hanks_kanamori_moment(magnitude):
    """
    Convert moment magnitude to seismic moment via Hanks–Kanamori relation.

    Relation: M_w = (2/3) * log₁₀(M₀) - 10.7
    
    Inverted: M₀ = 10^(1.5 * M_w + 9.1) for units of N·m (SI).
    This constant is widely used in seismic literature and ensures that
    a Mw=8 earthquake corresponds to roughly 2.5e22 N·m of seismic moment.

    Args:
        magnitude (float or ndarray): Moment magnitude M_w.

    Returns:
        float or ndarray: Seismic moment M₀ [N·m].

    Notes:
        Earlier versions mistakenly used a constant of 4.4, which is the
        offset for conversions to dyne·cm rather than N·m; that error
        produced moments five orders of magnitude too small and led to
        unrealistically tiny slip after scaling.
    """
    # Standard SI formula (M₀ in N·m) with correct offset
    # log10(M0) = 1.5*M + 9.1
    return 10.0**(1.5 * magnitude + 9.1)


def hanks_kanamori_magnitude(moment):
    """
    Convert seismic moment to moment magnitude via Hanks–Kanamori relation.

    Args:
        moment (float or ndarray): Seismic moment M₀ [N·m].

    Returns:
        float or ndarray: Moment magnitude M_w.
    """
    return (np.log10(moment) - 4.4) / 1.5


class HazardScenarioWeights:
    """
    Compute and manage scenario weights for hazard aggregation.

    For PTHA, we typically discretize earthquake scenarios by magnitude bins
    and then compute the weight of each scenario as:

        w_i = ΔN_i = N(M_i) - N(M_{i+1})

    where N(M) is the annual rate of earthquakes ≥ M.
    """

    def __init__(self, gr_relation):
        """
        Initialize scenario weighting from a Gutenberg–Richter relation.

        Args:
            gr_relation (GutenbergRichterRelation): The seismic activity model.
        """
        self.gr = gr_relation

    def weight_magnitude_bin(self, M_lower, M_upper):
        """
        Compute weight (rate) for a magnitude bin [M_lower, M_upper).

        Weight = N(M_lower) - N(M_upper)
               = (rate of M ≥ M_lower) - (rate of M ≥ M_upper)

        Args:
            M_lower (float): Lower magnitude bound (inclusive).
            M_upper (float): Upper magnitude bound (exclusive).

        Returns:
            float: Annual occurrence rate for bin, [events/year].
        """
        rate_lower = self.gr.annual_rate(M_lower)
        rate_upper = self.gr.annual_rate(M_upper)
        return rate_lower - rate_upper

    def normalize_weights(self, magnitude_bins):
        """
        Compute and normalize scenario weights for a set of magnitude bins.

        Args:
            magnitude_bins (list or ndarray): Magnitude bin edges [M_min, M_1, M_2, ..., M_max].
                                             Must be sorted.

        Returns:
            dict: Dictionary with keys:
                'bin_centers': Center magnitude of each bin
                'bin_weights': Un-normalized annual rates [events/year]
                'bin_probabilities': Normalized probabilities (sum=1)
        """
        bins = np.sort(magnitude_bins)
        n_bins = len(bins) - 1
        
        centers = 0.5 * (bins[:-1] + bins[1:])
        weights = np.array([self.weight_magnitude_bin(bins[i], bins[i+1]) 
                           for i in range(n_bins)])
        weights = np.maximum(weights, 0.0)  # Ensure non-negative
        
        total_weight = np.sum(weights)
        probs = weights / total_weight if total_weight > 0 else weights
        
        return {
            'bin_centers': centers,
            'bin_weights': weights,
            'bin_probabilities': probs,
            'total_annual_rate': total_weight,
        }


if __name__ == "__main__":
    # Demonstration
    print("=" * 60)
    print("Gutenberg–Richter Magnitude-Frequency Demonstration")
    print("=" * 60)
    
    # Create a GR relation: a=5.0, b=1.0 (typical for subduction zones)
    gr = GutenbergRichterRelation(a_value=5.0, b_value=1.0, 
                                  M_min=7.0, M_max=9.5)
    
    print(f"Parameters: a={gr.a}, b={gr.b}, M_min={gr.M_min}, M_max={gr.M_max}")
    print()
    
    # Compute annual rates
    magnitudes = np.arange(7.0, 9.6, 0.5)
    rates = gr.annual_rate(magnitudes)
    print("Annual occurrence rates:")
    for M, rate in zip(magnitudes, rates):
        print(f"  M ≥ {M:.1f}: {rate:.4e} events/year")
    print()
    
    # Hanks–Kanamori conversion
    print("Hanks–Kanamori moment-magnitude relation:")
    for M in [7.5, 8.0, 8.5, 9.0]:
        M0 = hanks_kanamori_moment(M)
        M_back = hanks_kanamori_magnitude(M0)
        print(f"  M_w={M:.1f} → M₀={M0:.3e} N·m → M_w={M_back:.2f}")
    print()
    
    # Scenario weights
    print("Scenario weights for magnitude bins:")
    mag_bins = [7.0, 7.5, 8.0, 8.5, 9.0, 9.5]
    weighter = HazardScenarioWeights(gr)
    result = weighter.normalize_weights(mag_bins)
    
    print(f"  Total annual rate: {result['total_annual_rate']:.4e} events/year")
    print()
    for i, (center, rate, prob) in enumerate(zip(result['bin_centers'], 
                                                   result['bin_weights'],
                                                   result['bin_probabilities'])):
        print(f"  Bin {i}: M={center:.2f}, Rate={rate:.4e}/yr, Probability={prob:.4f}")
