"""
Phase 4: Hazard Aggregation and Hazard Curve Computation

This script aggregates the tsunami displacement ensembles from Phase 3 with
the seismicity rates from Phase 2 to compute:
  1. Annual exceedance probability (AEP) at each intensity level h
  2. Hazard curves: AEP(h) vs h
  3. Return period estimates (e.g., 475-yr, 2500-yr events)
  4. Probabilistic tsunami hazard map

For each observation point and magnitude scenario:
  - Load displacement samples from Phase 3
  - For each threshold height h:
    * Count fraction of samples that exceed h
    * Weight by annual rate for that magnitude
    * Accumulate to AEP(h)

Usage:
    python phase4_hazard_aggregation.py
"""

import numpy as np
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))


def load_displacement_ensemble(magnitude_str, displacement_root='output/tsunami_sources'):
    """Load all displacement fields for a given magnitude.

    Args:
        magnitude_str (str): 'Mw8.0', 'Mw8.5', 'Mw9.0'
        displacement_root (str): Root directory containing displacement fields

    Returns:
        tuple: (displacements array shape (n_samples, n_obs_pts),
                sample_indices list)
    """
    mag_dir = Path(displacement_root) / magnitude_str
    if not mag_dir.exists():
        raise FileNotFoundError(f"Displacement directory not found: {mag_dir}")

    disp_files = sorted(mag_dir.glob('*.npy'))
    print(f"  Loading {len(disp_files)} displacement fields from {mag_dir}")

    displacements = []
    for f in disp_files:
        d = np.load(f)
        displacements.append(d)

    # Stack into array: shape (n_samples, n_obs_pts)
    displacements = np.array(displacements)
    return displacements


def compute_hazard_curves(displacement_ensembles, rates_dict, magnitude_list,
                          intensity_levels=None):
    """Compute aggregated hazard curves.

    Args:
        displacement_ensembles (dict): {magnitude_str: array of shape (n_samples, n_obs_pts)}
        rates_dict (dict): {magnitude_str: annual rate}
        magnitude_list (list): Ordered list of magnitudes (e.g., [8.0, 8.5, 9.0])
        intensity_levels (ndarray, optional): Displacement levels to evaluate [m].
            If None, auto-compute from data quantiles.

    Returns:
        dict: Results keyed by observation point index:
            'aep_values': array of AEP (annual) for each intensity level
            'intensities': array of intensity levels
            'return_periods': array of return periods in years
    """
    # Concatenate all displacements across magnitudes
    all_displacements = []
    all_weights = []  # weight for each sample

    for mag in magnitude_list:
        mag_str = f"Mw{mag:.1f}"
        disps = displacement_ensembles[mag_str]
        rate = rates_dict[mag_str]

        n_samples = disps.shape[0]
        weight = rate / n_samples  # weight per sample

        all_displacements.append(disps)
        all_weights.extend([weight] * n_samples)

    all_displacements = np.vstack(all_displacements)
    all_weights = np.array(all_weights)

    n_obs_pts = all_displacements.shape[1]

    # Define intensity levels if not provided
    if intensity_levels is None:
        # Use quantiles from the full ensemble
        all_max_disps = np.max(np.abs(all_displacements), axis=1)
        intensity_levels = np.percentile(
            all_max_disps,
            np.linspace(10, 99, 20)
        )

    # Compute AEP at each observation point
    hazard_results = {}

    for obs_idx in range(n_obs_pts):
        obs_displacements = all_displacements[:, obs_idx]

        # For each intensity level, compute P(disp > h)
        aep_values = []
        for h in intensity_levels:
            # Fraction of samples (weighted) that exceed h
            exceeds = obs_displacements > h
            aep = np.sum(all_weights[exceeds])
            aep_values.append(aep)

        aep_values = np.array(aep_values)

        # Convert AEP to return period (years)
        # RT = 1 / AEP, with AEP in annual units
        # Avoid division by zero
        return_periods = np.full_like(aep_values, np.inf)
        valid_aep = aep_values > 0
        return_periods[valid_aep] = 1.0 / aep_values[valid_aep]

        hazard_results[obs_idx] = {
            'aep_values': aep_values,
            'intensities': intensity_levels.copy(),
            'return_periods': return_periods,
            'max_displacement_samples': obs_displacements,
        }

    return hazard_results


def compute_site_statistics(hazard_results):
    """Compute summary stats for each observation point.

    Args:
        hazard_results (dict): Output from compute_hazard_curves

    Returns:
        dict: Summary statistics
    """
    stats = {}
    for obs_idx, result in hazard_results.items():
        aep = result['aep_values']
        h = result['intensities']

        # Find displacement at key return periods
        disps_475yr = np.interp(475.0, result['return_periods'][::-1], h[::-1], left=np.nan, right=np.nan)
        disps_2500yr = np.interp(2500.0, result['return_periods'][::-1], h[::-1], left=np.nan, right=np.nan)

        stats[obs_idx] = {
            '475yr_displacement_m': disps_475yr,
            '2500yr_displacement_m': disps_2500yr,
            'median_displacement_m': np.median(result['max_displacement_samples']),
            'mean_displacement_m': np.mean(result['max_displacement_samples']),
            'max_aep_per_year': np.max(aep),
        }

    return stats


def main():
    print("=" * 70)
    print("Phase 4: Hazard Aggregation")
    print("=" * 70)

    # Load seismicity rates
    rates_file = Path('output') / 'phase2_rates.json'
    with open(rates_file) as f:
        rates_data = json.load(f)
    rates_dict = rates_data['rates']

    # Filter to scenarios we have displacements for
    magnitudes = [8.0, 8.5, 9.0]
    rates_subset = {f"Mw{m:.1f}": rates_dict[f"Mw{m:.1f}"] for m in magnitudes}

    print(f"\nSeismicity rates (annual):")
    for mag_str, rate in rates_subset.items():
        print(f"  {mag_str}: {rate:.4f} events/year")

    # Load displacement ensembles
    print(f"\nLoading displacement ensembles...")
    displacement_ensembles = {}
    for mag in magnitudes:
        mag_str = f"Mw{mag:.1f}"
        displacement_ensembles[mag_str] = load_displacement_ensemble(mag_str)

    # Sample shapes
    for mag_str, disps in displacement_ensembles.items():
        print(f"  {mag_str}: {disps.shape[0]} samples × {disps.shape[1]} observation points")

    # Compute hazard curves
    print(f"\nComputing hazard curves...")
    hazard_results = compute_hazard_curves(
        displacement_ensembles,
        rates_subset,
        magnitudes,
    )
    print(f"  Computed curves for {len(hazard_results)} observation points")

    # Compute site-level statistics
    print(f"\nComputing site-level summary statistics...")
    site_stats = compute_site_statistics(hazard_results)

    # Summary output
    print(f"\n" + "=" * 70)
    print("HAZARD CURVE STATISTICS (Selected Sites)")
    print("=" * 70)

    # Report on a few representative sites
    n_obs = len(site_stats)
    sample_indices = [0, n_obs // 4, n_obs // 2, 3 * n_obs // 4, n_obs - 1]

    for idx in sample_indices:
        if idx not in site_stats:
            continue
        stats = site_stats[idx]
        print(f"\nObservation point {idx}:")
        print(f"  Median displacement: {stats['median_displacement_m']:.3e} m")
        print(f"  Mean displacement: {stats['mean_displacement_m']:.3e} m")
        print(f"  475-yr displacement: {stats['475yr_displacement_m']:.3e} m")
        print(f"  2500-yr displacement: {stats['2500yr_displacement_m']:.3e} m")
        print(f"  Maximum AEP: {stats['max_aep_per_year']:.3e} / year")

    # Save results to JSON
    output_file = Path('output') / 'phase4_hazard_curves.json'
    output_data = {}
    for obs_idx, result in hazard_results.items():
        samples = result['max_displacement_samples']
        output_data[str(obs_idx)] = {
            'intensities_m': result['intensities'].tolist(),
            'aep_per_year': result['aep_values'].tolist(),
            'return_periods_yr': result['return_periods'].tolist(),
            'median_displacement_m': float(np.median(samples)),
            'mean_displacement_m': float(np.mean(samples)),
            'p05_displacement_m': float(np.percentile(samples, 5)),
            'p95_displacement_m': float(np.percentile(samples, 95)),
        }

    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    print(f"\nHazard curves saved to {output_file}")

    # Save site stats
    stats_file = Path('output') / 'phase4_site_statistics.json'
    output_stats = {}
    for obs_idx, stats in site_stats.items():
        output_stats[str(obs_idx)] = {k: float(v) if not np.isnan(v) else None
                                       for k, v in stats.items()}

    with open(stats_file, 'w') as f:
        json.dump(output_stats, f, indent=2)
    print(f"Site statistics saved to {stats_file}")

    print(f"\nPhase 4 complete.")

    return hazard_results, site_stats


if __name__ == '__main__':
    hazard_results, site_stats = main()
