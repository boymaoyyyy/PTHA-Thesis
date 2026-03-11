"""
Phase 1 Alternative: Slip generation using standard Monte Carlo (MC)

This script mirrors phase1_generate.py but uses ordinary Monte Carlo sampling
instead of RQMC. It will allow comparison of RQMC vs MC convergence and quality.

Usage:
    python phase1_generate_mc.py
"""

import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))

from slip_sampler import StochasticSlipGenerator
from philippine_trench_config import PHILIPPINE_TRENCH_PTHA_CONFIG, RUPTURE_SCENARIOS


def compute_subfault_size(length_km, width_km, target_n):
    """Compute a subfault size (km) so that n_along * n_down ~ target_n."""
    return np.sqrt((length_km * width_km) / target_n)


def generate_slip_sample_mc(magnitude, sample_index, subfault_size_km, 
                            fault_params_csv, seed=None):
    """Generate a single slip sample using standard Monte Carlo (Gaussian random sampling).
    
    Instead of using Sobol sequences (RQMC), we generate random normal samples
    directly and apply the same covariance structure.
    """
    from slip_sampler import DATA_DIR
    
    if fault_params_csv is None:
        fault_params_csv = DATA_DIR / "heidarzadeh_2025_table3.csv"
    
    generator = StochasticSlipGenerator(
        magnitude,
        fault_params_csv,
        subfault_size_km=subfault_size_km,
    )
    
    # Instead of using RQMC Sobol, use standard Gaussian random sampling
    if seed is not None:
        np.random.seed(seed + sample_index)
    else:
        np.random.seed(sample_index)
    
    # Generate standard normal random vector (size = number of subfaults)
    xi = np.random.randn(generator.N)
    
    # Apply covariance structure via Cholesky (same as RQMC version)
    log_slip = np.log(generator.slip_mean_m) + generator.L_chol @ xi
    slip_raw = np.exp(log_slip)
    
    # Moment scaling
    M_raw = generator.rigidity * generator.area_subfault_m2 * np.sum(slip_raw)
    scale_factor = generator.M0_target / M_raw
    slip_final = slip_raw * scale_factor
    
    # Validation
    M_final = generator.rigidity * generator.area_subfault_m2 * np.sum(slip_final)
    rel_err = abs(generator.M0_target - M_final) / generator.M0_target
    assert rel_err < 1e-6, f"Moment scaling error: {rel_err:.2e}"
    
    slip_grid = slip_final.reshape(generator.n_down, generator.n_along)
    return slip_grid


def main():
    cfg = PHILIPPINE_TRENCH_PTHA_CONFIG
    slip_cfg = cfg.slip_params

    csv_path = Path("data") / "fault_geometry" / "heidarzadeh_2025_table3.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Fault geometry CSV not found: {csv_path}")

    n_samples = slip_cfg.n_samples_per_magnitude[0]
    seed = 42  # fixed for reproducibility

    print("Phase 1 (Alternative): Slip Generation via Standard Monte Carlo")
    print("================================================================")
    print(f"Using slip parameters: {slip_cfg}")
    print(f"Sampling method: Standard Monte Carlo (NOT RQMC)")
    print(f"Generating {n_samples} samples per scenario, seed={seed}\n")

    # Read available magnitudes from the CSV
    import pandas as pd
    available = pd.read_csv(csv_path)["Hypothetical scenario"].str.replace("Mw ", "").astype(float).tolist()

    for name, scenario in RUPTURE_SCENARIOS.items():
        mag = scenario.magnitude
        if mag not in available:
            print(f"Skipping scenario {name} (Mw {mag}) – not present in CSV")
            continue

        print(f"Processing scenario {name} (Mw {mag})")

        subfault_size = compute_subfault_size(scenario.length_km,
                                              scenario.width_km,
                                              slip_cfg.n_subfaults)
        print(f"  target subfault count: {slip_cfg.n_subfaults},"
              f" computed size: {subfault_size:.2f} km")

        # Create output directory for MC samples
        output_dir = Path("output") / "slip_samples_mc" / f"Mw{mag:.1f}"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate and save each realization using standard MC
        for idx in range(n_samples):
            slip = generate_slip_sample_mc(mag, idx, subfault_size, csv_path, 
                                          seed=seed)
            out_path = output_dir / f"slip_Mw{mag:.1f}_sample{idx}.npy"
            np.save(out_path, slip)
        
        print(f"  generated and saved {n_samples} slip realizations to {output_dir}")
        print()

    print("Phase 1 (MC) generation complete. Slip arrays are stored in output/slip_samples_mc/")


if __name__ == "__main__":
    main()
