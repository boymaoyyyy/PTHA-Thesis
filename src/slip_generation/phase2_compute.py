#!/usr/bin/env python3
"""
Phase 2 Gutenberg-Richter recurrence helper for Philippine Trench PTHA.

This version fixes the earlier issues by:
1. Restricting the operational magnitude range to Mw 7.0–9.0
2. Using 0.5-magnitude bins consistently
3. Treating 2012 and 2023 Mw 7.6 events as validation-only, not hazard bins
4. Supporting optional anchoring to an annual exceedance rate for Mw >= 7.0
5. Outputting clean bin edges, centers, relative weights, probabilities,
   and optionally absolute annual rates per bin

Usage:
    python phase2_compute_fixed.py
    python phase2_compute_fixed.py --annual-rate-at-mw7 0.02
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np


def beta_to_b(beta: float) -> float:
    """
    Convert beta in exponential form to Gutenberg-Richter b-value.

    For natural-log exponential recurrence:
        f(M) ~ beta * exp(-beta * M)

    and classical GR:
        log10 N = a - b M

    the conversion is:
        b = beta / ln(10)
    """
    return float(beta / math.log(10.0))


def build_edges(m_min: float = 7.0, m_max: float = 9.0, delta_m: float = 0.5) -> np.ndarray:
    """
    Build bin edges for the operational range.
    Example: 7.0, 7.5, 8.0, 8.5, 9.0
    """
    n = int(round((m_max - m_min) / delta_m))
    edges = np.array([m_min + i * delta_m for i in range(n + 1)], dtype=float)
    return edges


def build_centers(edges: np.ndarray) -> np.ndarray:
    return 0.5 * (edges[:-1] + edges[1:])


def exponential_bin_probabilities(beta: float, edges: np.ndarray) -> np.ndarray:
    """
    Compute conditional bin probabilities over [M_min, M_max] using
    a truncated exponential recurrence model:

        p(M in [m_i, m_{i+1})) proportional to
            exp(-beta * m_i) - exp(-beta * m_{i+1})

    normalized over [edges[0], edges[-1]].
    """
    m0 = float(edges[0])
    m1 = float(edges[-1])

    denom = math.exp(-beta * m0) - math.exp(-beta * m1)
    if denom <= 0:
        raise ValueError("Invalid beta or magnitude range; denominator must be positive.")

    probs = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        num = math.exp(-beta * float(lo)) - math.exp(-beta * float(hi))
        probs.append(num / denom)

    probs = np.asarray(probs, dtype=float)
    probs /= probs.sum()
    return probs


def annual_bin_rates_from_anchor(
    probabilities: np.ndarray,
    annual_rate_at_mw7: Optional[float],
) -> Optional[np.ndarray]:
    """
    If annual_rate_at_mw7 is supplied, treat it as the annual exceedance rate
    for Mw >= 7.0 over the operational range, and allocate it across bins
    according to the normalized probabilities.
    """
    if annual_rate_at_mw7 is None:
        return None
    if annual_rate_at_mw7 <= 0:
        raise ValueError("--annual-rate-at-mw7 must be positive if supplied.")
    return annual_rate_at_mw7 * probabilities


def incremental_rates_untruncated(beta: float, edges: np.ndarray, annual_rate_at_mw7: float) -> np.ndarray:
    """
    Alternative absolute-rate interpretation:
    Given lambda(M >= 7.0) = annual_rate_at_mw7, compute the absolute
    incremental annual rate in each bin under an untruncated exponential model.

    lambda(M >= m) = lambda(M >= 7.0) * exp[-beta * (m - 7.0)]

    Then each bin rate is:
        lambda([m_i, m_{i+1})) = lambda(M >= m_i) - lambda(M >= m_{i+1})

    This is often more interpretable for hazard calculations.
    """
    mref = float(edges[0])
    rates = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        lam_lo = annual_rate_at_mw7 * math.exp(-beta * (float(lo) - mref))
        lam_hi = annual_rate_at_mw7 * math.exp(-beta * (float(hi) - mref))
        rates.append(max(lam_lo - lam_hi, 0.0))
    return np.asarray(rates, dtype=float)


def validation_only_events() -> Dict[str, Dict[str, float]]:
    """
    Keep historical events separate from operational hazard bins.
    """
    return {
        "2012_Mw7.6": {"magnitude": 7.6},
        "2023_Mw7.6": {"magnitude": 7.6},
    }


def build_output(
    beta: float,
    edges: np.ndarray,
    annual_rate_at_mw7: Optional[float],
) -> Dict:
    centers = build_centers(edges)
    probabilities = exponential_bin_probabilities(beta, edges)
    relative_weights = probabilities.copy()

    anchored_bin_rates = annual_bin_rates_from_anchor(probabilities, annual_rate_at_mw7)

    output: Dict = {
        "model": "truncated_exponential_gutenberg_richter",
        "beta": float(beta),
        "b": beta_to_b(beta),
        "operational_magnitude_range": {
            "min": float(edges[0]),
            "max": float(edges[-1]),
            "delta_m": float(edges[1] - edges[0]),
        },
        "validation_only_events": validation_only_events(),
        "weights": {
            "bin_edges": edges.tolist(),
            "bin_centers": centers.tolist(),
            "relative_weights": relative_weights.tolist(),
            "probabilities": probabilities.tolist(),
        },
    }

    if annual_rate_at_mw7 is not None:
        output["annualization"] = {
            "anchor_definition": "annual exceedance rate for Mw >= 7.0",
            "anchor_value_per_year": float(annual_rate_at_mw7),
            "bin_rates_from_normalized_probabilities_per_year": anchored_bin_rates.tolist(),
            "bin_rates_untruncated_exponential_per_year": incremental_rates_untruncated(
                beta, edges, annual_rate_at_mw7
            ).tolist(),
        }

    return output


def print_summary(result: Dict) -> None:
    print("Phase 2: Gutenberg-Richter Seismicity (Fixed)")
    print("=============================================")
    print(f"Model: {result['model']}")
    print(f"beta = {result['beta']:.3f}")
    print(f"b    = {result['b']:.6f}")
    print(
        "Operational magnitude range: "
        f"Mw {result['operational_magnitude_range']['min']:.1f} "
        f"to {result['operational_magnitude_range']['max']:.1f} "
        f"(ΔM = {result['operational_magnitude_range']['delta_m']:.1f})"
    )

    print("\nValidation-only events:")
    for name, meta in result["validation_only_events"].items():
        print(f"  {name}: Mw {meta['magnitude']:.1f}")

    w = result["weights"]
    print("\nOperational hazard bins:")
    print("  Bin Range      Center    Relative Weight    Probability")
    print("  -------------------------------------------------------")
    edges = w["bin_edges"]
    centers = w["bin_centers"]
    rel = w["relative_weights"]
    probs = w["probabilities"]
    for i in range(len(centers)):
        print(
            f"  {edges[i]:.1f}-{edges[i+1]:.1f}      "
            f"{centers[i]:6.2f}    "
            f"{rel[i]:15.6f}    "
            f"{probs[i]:11.6f}"
        )

    annual = result.get("annualization")
    if annual is not None:
        print("\nAnnualization:")
        print(f"  Anchor: {annual['anchor_definition']}")
        print(f"  Value : {annual['anchor_value_per_year']:.6f} /yr")

        print("\n  Bin rates from normalized probabilities:")
        for i, r in enumerate(annual["bin_rates_from_normalized_probabilities_per_year"]):
            print(f"    {edges[i]:.1f}-{edges[i+1]:.1f}: {r:.8f} /yr")

        print("\n  Bin rates from untruncated exponential:")
        for i, r in enumerate(annual["bin_rates_untruncated_exponential_per_year"]):
            print(f"    {edges[i]:.1f}-{edges[i+1]:.1f}: {r:.8f} /yr")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fixed Phase 2 Gutenberg-Richter recurrence helper")
    p.add_argument("--beta", type=float, default=0.9, help="Central beta value. Default: 0.9")
    p.add_argument("--beta-low", type=float, default=0.7, help="Low sensitivity beta. Default: 0.7")
    p.add_argument("--beta-high", type=float, default=1.1, help="High sensitivity beta. Default: 1.1")
    p.add_argument("--m-min", type=float, default=7.0, help="Operational minimum magnitude. Default: 7.0")
    p.add_argument("--m-max", type=float, default=9.0, help="Operational maximum magnitude. Default: 9.0")
    p.add_argument("--delta-m", type=float, default=0.5, help="Magnitude bin width. Default: 0.5")
    p.add_argument(
        "--annual-rate-at-mw7",
        type=float,
        default=None,
        help="Optional annual exceedance rate for Mw >= 7.0",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("output") / "phase2_rates_fixed.json",
        help="Output JSON path",
    )
    p.add_argument(
        "--write-sensitivity",
        action="store_true",
        help="Also write beta sensitivity outputs for beta-low, beta, beta-high",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    edges = build_edges(args.m_min, args.m_max, args.delta_m)

    result = build_output(
        beta=args.beta,
        edges=edges,
        annual_rate_at_mw7=args.annual_rate_at_mw7,
    )
    print_summary(result)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved fixed Phase 2 output to: {args.output}")

    if args.write_sensitivity:
        sensitivity = {}
        for label, beta in [("low", args.beta_low), ("central", args.beta), ("high", args.beta_high)]:
            sensitivity[label] = build_output(
                beta=beta,
                edges=edges,
                annual_rate_at_mw7=args.annual_rate_at_mw7,
            )

        sens_path = args.output.with_name(args.output.stem + "_sensitivity.json")
        with open(sens_path, "w", encoding="utf-8") as f:
            json.dump(sensitivity, f, indent=2)
        print(f"Saved sensitivity output to: {sens_path}")


if __name__ == "__main__":
    main()