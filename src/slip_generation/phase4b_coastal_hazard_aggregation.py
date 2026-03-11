"""
Phase 4b: Coastal Hazard Aggregation & Curves

This script aggregates coastal inundation ensembles from Phase 3.5 with
seismicity rates from Phase 2 to compute site-specific hazard curves.

Output:
  1. Coastal hazard curves (AEP vs inundation depth) at each tide gauge site
  2. Annual exceedance probabilities (AEP) at key inundation thresholds
  3. Return-period maps and site statistics

Key insight: Coastal hazard is dominated by large-magnitude events (Mw 8.5, 9.0)
due to nonlinear amplification, even though smaller events are more frequent.

Usage:
    python phase4b_coastal_hazard_aggregation.py
"""

import numpy as np
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))


def load_coastal_inundation_ensembles(
    coastal_file='output/coastal_inundation/coastal_inundation_ensembles.json'
):
    """Load coastal inundation data from JSON.

    Returns:
        dict: {magnitude_str: {site_name: {...}}}
    """
    with open(coastal_file) as f:
        return json.load(f)


def load_seismicity_rates(rates_file='output/phase2_rates.json'):
    """Load seismicity rates from Phase 2.

    Returns:
        dict: {magnitude_str: annual_rate}
    """
    with open(rates_file) as f:
        data = json.load(f)
    return data['rates']


def compute_coastal_hazard_curves(coastal_ensembles, rates_dict,
                                   inundation_thresholds=None):
    """
    Compute hazard curves aggregating across magnitudes.

    For each site and inundation threshold h:
        AEP(h) = Σ_m λ_m * P(inundation > h | m)

    Args:
        coastal_ensembles (dict): {Mw: {site: {inundation_samples: [...]}}}
        rates_dict (dict): {Mw: annual_rate}
        inundation_thresholds (ndarray, optional): Levels to evaluate [m]

    Returns:
        dict: {site: {thresholds, aep_values, return_periods, ...}}
    """
    # Extract magnitude strings in order
    magnitudes = sorted([k for k in coastal_ensembles.keys()
                        if 'Mw' in k and k in rates_dict])

    # Collect all samples weighted by rate
    all_sites = list(coastal_ensembles[magnitudes[0]].keys())

    # Build rate-weighted ensemble for each site
    site_ensembles = {}

    for site in all_sites:
        all_samples = []
        all_weights = []

        for mag in magnitudes:
            samples = coastal_ensembles[mag][site]['inundation_samples']
            rate = rates_dict[mag]

            n_samples = len(samples)
            weight = rate / n_samples  # weight per sample

            all_samples.extend(samples)
            all_weights.extend([weight] * n_samples)

        all_samples = np.array(all_samples)
        all_weights = np.array(all_weights)

        # Define thresholds if not provided
        if inundation_thresholds is None:
            thresholds = np.percentile(all_samples, np.linspace(5, 99, 30))
        else:
            thresholds = inundation_thresholds

        # Compute AEP at each threshold
        aep_values = []
        for h in thresholds:
            exceeds = all_samples > h
            aep = np.sum(all_weights[exceeds])
            aep_values.append(aep)

        aep_values = np.array(aep_values)

        # Return periods
        return_periods = np.full_like(aep_values, np.inf)
        valid = aep_values > 0
        return_periods[valid] = 1.0 / aep_values[valid]

        site_ensembles[site] = {
            'inundation_thresholds_m': thresholds.tolist(),
            'aep_per_year': aep_values.tolist(),
            'return_periods_yr': return_periods.tolist(),
            'samples': all_samples.tolist(),
            'sample_weights': all_weights.tolist(),
            'median': float(np.median(all_samples)),
            'mean': float(np.mean(all_samples)),
            'std': float(np.std(all_samples)),
        }

    return site_ensembles


def compute_key_return_periods(site_hazards, return_periods=[475, 2500]):
    """
    Extract inundation depths at key return periods for each site.

    Args:
        site_hazards (dict): Output from compute_coastal_hazard_curves
        return_periods (list): Target return periods [years]

    Returns:
        dict: {site: {return_period: inundation_m}}
    """
    results = {}

    for site, data in site_hazards.items():
        thresholds = np.array(data['inundation_thresholds_m'])
        rps = np.array(data['return_periods_yr'])

        results[site] = {}
        for target_rp in return_periods:
            # Interpolate: find inundation depth at this return period
            # (rps is sorted descending because AEP decreases)
            if target_rp <= rps.max() and target_rp >= rps.min():
                inund_at_rp = np.interp(
                    target_rp,
                    rps[::-1],
                    thresholds[::-1],
                    left=np.nan,
                    right=np.nan
                )
            else:
                inund_at_rp = np.nan

            results[site][f'{int(target_rp)}yr'] = inund_at_rp

    return results


def main():
    print("=" * 70)
    print("Phase 4b: Coastal Hazard Aggregation")
    print("=" * 70)

    # Load data
    print("\nLoading coastal inundation ensembles...")
    coastal_ensembles = load_coastal_inundation_ensembles()

    print("Loading seismicity rates...")
    rates_dict = load_seismicity_rates()

    # Filter rates to magnitudes we have coastal data for
    mags = sorted([m for m in coastal_ensembles.keys() if 'Mw' in m])
    rates_subset = {m: rates_dict[m] for m in mags}

    print(f"\nMagnitudes: {', '.join(mags)}")
    print(f"Seismicity rates (annual):")
    for mag, rate in rates_subset.items():
        print(f"  {mag}: {rate:.4f} events/year")

    # Compute coastal hazard curves
    print("\nComputing coastal hazard curves...")
    site_hazards = compute_coastal_hazard_curves(coastal_ensembles, rates_subset)

    print(f"  Computed curves for {len(site_hazards)} coastal sites:")
    for site in site_hazards.keys():
        print(f"    - {site}")

    # Extract key return periods
    print("\nExtracting inundation depths at key return periods...")
    key_inundations = compute_key_return_periods(site_hazards)

    # Save results
    output_dir = Path('output') / 'coastal_hazard'
    output_dir.mkdir(parents=True, exist_ok=True)

    curves_file = output_dir / 'coastal_hazard_curves.json'
    with open(curves_file, 'w') as f:
        json.dump(site_hazards, f, indent=2)
    print(f"\nCoastal hazard curves saved to {curves_file}")

    returns_file = output_dir / 'key_return_periods.json'
    with open(returns_file, 'w') as f:
        json.dump(key_inundations, f, indent=2)
    print(f"Return period inundations saved to {returns_file}")

    # Summary output
    print("\n" + "=" * 70)
    print("COASTAL HAZARD SUMMARY")
    print("=" * 70)

    for site, data in site_hazards.items():
        print(f"\n{site}:")
        print(f"  Median inundation: {data['median']:.2e} m")
        print(f"  Mean inundation: {data['mean']:.2e} m")
        print(f"  Std. dev: {data['std']:.2e} m")

        print(f"  Inundation at return periods:")
        for rp in [475, 2500]:
            key = f'{rp}yr'
            if key in key_inundations[site]:
                val = key_inundations[site][key]
                if not np.isnan(val):
                    print(f"    {rp}-yr: {val:.2e} m")
                else:
                    print(f"    {rp}-yr: exceeds range")

    print("\n" + "=" * 70)
    print("Phase 4b complete. Coastal hazard assessment ready.")
    print("=" * 70)

    return site_hazards, key_inundations


if __name__ == '__main__':
    site_hazards, key_inundations = main()
