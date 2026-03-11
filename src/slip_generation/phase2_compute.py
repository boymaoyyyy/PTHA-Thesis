"""
Phase 2 Gutenberg-Richter seismicity helper for Philippine Trench PTHA.

This script reads the seismicity parameters from
`philippine_trench_config.PHILIPPINE_TRENCH_PTHA_CONFIG` and evaluates
annual occurrence rates for the prescribed rupture scenarios.  If an
annual rate at the minimum magnitude (M_min) has been estimated from the
catalog, it will be used to infer the GR "a" parameter; otherwise a
reasonable default is chosen.

Usage:
    python phase2_compute.py

The script prints a summary table and optionally writes a JSON file
containing the scenario weights (annual rates and probabilities) which
may be consumed by later phases.
"""

import json
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))

from magnitude_frequency import GutenbergRichterRelation, HazardScenarioWeights
from philippine_trench_config import PHILIPPINE_TRENCH_PTHA_CONFIG, RUPTURE_SCENARIOS


def compute_a_value(b_value, M_min, rate_at_Mmin=None):
    """Compute GR 'a' parameter given slope and optional observed rate.

    If rate_at_Mmin is provided (>0), we set a = log10(rate) + b*M_min.
    Otherwise we use a default of 2.5, which corresponds to about 0.01
    events/yr at M=8 for b~0.4 – a reasonably conservative proxy used in
    the demonstration script.
    """
    if rate_at_Mmin and rate_at_Mmin > 0:
        return np.log10(rate_at_Mmin) + b_value * M_min
    else:
        return 2.5


def main():
    cfg = PHILIPPINE_TRENCH_PTHA_CONFIG
    seis = cfg.seismicity_params

    print("Phase 2: Gutenberg-Richter Seismicity")
    print("========================================")
    print(f"Seismicity parameters from configuration: {seis}")

    # infer a-value
    a_val = compute_a_value(seis.b_value, seis.magnitude_minimum,
                             seis.annual_rate_at_Mw7p0)
    print(f"Using a-value = {a_val:.3f} (b = {seis.b_value:.3f})")

    gr = GutenbergRichterRelation(a_value=a_val,
                                  b_value=seis.b_value,
                                  M_min=seis.magnitude_minimum,
                                  M_max=seis.magnitude_maximum)

    mags = sorted({scenario.magnitude for scenario in RUPTURE_SCENARIOS.values()})
    print("\nAnnual occurrence rates for rupture scenarios:")
    print("  Scenario        Magnitude  Rate (events/yr)  Return period (yr)")
    print("  -------------------------------------------------------------")

    rates = {}
    for name, scenario in RUPTURE_SCENARIOS.items():
        rate = gr.annual_rate(scenario.magnitude)
        rp = 1.0 / rate if rate > 0 else float('inf')
        rates[name] = rate
        print(f"  {name:15} {scenario.magnitude:9.2f}  {rate:16.3e}  {rp:12.1f}")

    # compute bin weights if desired
    bins = sorted(list({scenario.magnitude for scenario in RUPTURE_SCENARIOS.values()}))
    # include extremes
    if cfg.seismicity_params.magnitude_minimum < bins[0]:
        bins.insert(0, cfg.seismicity_params.magnitude_minimum)
    if cfg.seismicity_params.magnitude_maximum > bins[-1]:
        bins.append(cfg.seismicity_params.magnitude_maximum)

    weighter = HazardScenarioWeights(gr)
    normed = weighter.normalize_weights(bins)

    print("\nComputed magnitude bin weights (annual rates & normalized proba):")
    for m, w, p in zip(normed['bin_centers'], normed['bin_weights'], normed['bin_probabilities']):
        print(f"  M = {m:.2f}: weight = {w:.3e} /yr, prob = {p:.3f}")

    # save results
    out = Path('output') / 'phase2_rates.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w') as f:
        json.dump({'a': a_val, 'b': seis.b_value,
                   'rates': rates,
                   'weights': {'centers': normed['bin_centers'].tolist(),
                               'weights': normed['bin_weights'].tolist(),
                               'probabilities': normed['bin_probabilities'].tolist()},
                   }, f, indent=2)
    print(f"\nSaved rates and weights to {out}")


if __name__ == '__main__':
    main()
