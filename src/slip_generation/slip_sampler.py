"""
Randomized Quasi-Monte Carlo (RQMC) slip generation for PTHA.

This module implements stochastic earthquake slip sampling using:
    - Sobol sequences (low-discrepancy sampling)
    - Owen's scrambling (randomization for independent samples)
    - Exponential/Matérn covariance models
    - Moment-scaling constraints via Hanks–Kanamori relation

The RQMC approach provides:
    - Better convergence than standard Monte Carlo
    - Reproducible pseudo-random behavior
    - Efficient exploration of slip parameter space

References:
    - Sobol (1967): On the distribution of points in a cube
    - Owen (1998): Scrambling Sobol and Niederreiter-Xing points
    - Dick & Pillichshammer (2010): Digital Nets and Sequences
"""

import numpy as np
import pandas as pd
from scipy.stats import qmc, norm
from pathlib import Path

from covariance import build_covariance_matrix, build_distance_matrix
from magnitude_frequency import hanks_kanamori_moment

# Resolve repository root relative to this file, then data/output paths
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data" / "fault_geometry"
OUT_DIR = REPO_ROOT / "output" / "slip_samples"

# Ensure output directory exists
OUT_DIR.mkdir(parents=True, exist_ok=True)


class StochasticSlipGenerator:
    """
    RQMC-based generator for stochastic earthquake slip distributions.

    Implements Phase 1 of the PTHA methodology:
    1. Load fault geometry from data
    2. Build spatial covariance matrix
    3. Generate Gaussian random field for log-slip via RQMC
    4. Apply moment-scaling constraint
    5. Export slip realizations
    """

    def __init__(self, magnitude, fault_params_csv, subfault_size_km=20.0,
                 covariance_model='exponential', correlation_length_km=20.0,
                 rigidity_pa=4.0e10):
        """
        Initialize slip generator.

        Args:
            magnitude (float): Target moment magnitude M_w.
            fault_params_csv (Path or str): Path to fault geometry CSV file.
            subfault_size_km (float): Subfault discretization size [km].
            covariance_model (str): 'exponential' or 'matern'.
            correlation_length_km (float): Correlation length ξ [km].
            rigidity_pa (float): Shear modulus μ [Pa]. Default: 40 GPa.
        """
        self.magnitude = magnitude
        self.subfault_size_km = subfault_size_km
        self.covariance_model = covariance_model
        self.correlation_length_km = correlation_length_km
        self.rigidity = rigidity_pa

        # Load fault parameters
        fault_data = pd.read_csv(fault_params_csv)
        row = fault_data[
            fault_data["Hypothetical scenario"] == f"Mw {magnitude:.1f}"
        ].iloc[0]

        self.L_km = row["L (km)"]
        self.W_km = row["W (km)"]
        self.top_depth_km = row["Top depth (km)"]
        self.strike_deg = row["Mean strike ( o )"]
        self.dip_deg = row["Mean dip ( o )"]
        self.rake_deg = row["Rake ( o )"]
        self.slip_mean_m = row["Slip (m)"]

        # Compute subfault grid dimensions
        self.n_along = int(np.round(self.L_km / self.subfault_size_km))
        self.n_down = int(np.round(self.W_km / self.subfault_size_km))
        self.N = self.n_along * self.n_down

        # Recompute grid resolution slightly to ensure exact fit
        self.subfault_length_km = self.L_km / self.n_along
        self.subfault_width_km = self.W_km / self.n_down
        self.area_subfault_m2 = (self.subfault_length_km * 1e3) ** 2

        # Target seismic moment
        self.M0_target = hanks_kanamori_moment(magnitude)

        # Build spatial covariance matrix
        self._build_covariance_matrix()

    def _build_covariance_matrix(self):
        """Build fault-plane coordinate grid and covariance matrix."""
        # Grid of subfault centers
        x_centers = (np.arange(self.n_along) + 0.5) * self.subfault_length_km
        y_centers = (np.arange(self.n_down) + 0.5) * self.subfault_width_km

        X, Y = np.meshgrid(x_centers, y_centers)
        self.coordinates = np.column_stack([X.ravel(), Y.ravel()])

        # Build covariance matrix
        self.cov_matrix = build_covariance_matrix(
            self.coordinates,
            model=self.covariance_model,
            correlation_length=self.correlation_length_km,
        )

        # Compute Cholesky decomposition
        try:
            self.L_chol = np.linalg.cholesky(self.cov_matrix)
        except np.linalg.LinAlgError:
            # Add additional regularization if needed
            eigvals = np.linalg.eigvalsh(self.cov_matrix)
            self.cov_matrix += (1e-8 - np.min(eigvals)) * np.eye(self.N)
            self.L_chol = np.linalg.cholesky(self.cov_matrix)

    def generate_rqmc_sample(self, sample_index, seed=None):
        """
        Generate a single slip realization using scrambled Sobol sequence.

        Method:
            1. Generate Sobol points with Owen scrambling (randomized QMC)
            2. Map to standard normal via inverse CDF
            3. Apply covariance correlation via Cholesky transform
            4. Exponentiate to get slip magnitudes
            5. Moment-scale to match target M_0

        Args:
            sample_index (int): Index in the Sobol sequence (0-indexed).
            seed (int, optional): Random seed for scrambling. Default: None.

        Returns:
            ndarray: Slip distribution, shape (n_down, n_along) [meters].

        Notes:
            - RQMC ensures systematic coverage of parameter space
            - Scrambling breaks lattice structure while preserving low discrepancy
            - Reproducible: same (sample_index, seed) → same result
        """
        # Initialize Sobol sampler with scrambling (Owen)
        sampler = qmc.Sobol(d=self.N, scramble=True, seed=seed)

        # Generate Sobol samples up to requested index
        if sample_index < 0:
            raise ValueError("sample_index must be >= 0")

        u_all = sampler.random(n=sample_index + 1)
        u = u_all[sample_index]  # Take the requested sample

        # Transform to standard normal
        xi = norm.ppf(u)

        # Apply covariance structure
        log_slip = np.log(self.slip_mean_m) + self.L_chol @ xi

        # Exponentiate to get slip
        slip_raw = np.exp(log_slip)

        # Moment scaling: adjust slip so total moment = M0_target
        M_raw = self.rigidity * self.area_subfault_m2 * np.sum(slip_raw)
        scale_factor = self.M0_target / M_raw
        slip_final = slip_raw * scale_factor

        # Validation check
        M_final = self.rigidity * self.area_subfault_m2 * np.sum(slip_final)
        rel_err = abs(self.M0_target - M_final) / self.M0_target
        assert rel_err < 1e-6, f"Moment scaling error: {rel_err:.2e}"

        # Reshape to grid
        slip_grid = slip_final.reshape(self.n_down, self.n_along)

        return slip_grid

    def generate_ensemble(self, n_samples, seed=None):
        """
        Generate an ensemble of RQMC slip samples.

        Args:
            n_samples (int): Number of samples to generate.
            seed (int, optional): Random seed for reproducibility.

        Returns:
            list: List of slip grids, each shape (n_down, n_along).
        """
        samples = []
        for idx in range(n_samples):
            slip = self.generate_rqmc_sample(idx, seed=seed)
            samples.append(slip)
        return samples

    def get_metadata(self):
        """Return metadata describing this generator configuration."""
        return {
            'magnitude': self.magnitude,
            'target_moment_nm': self.M0_target,
            'fault_length_km': self.L_km,
            'fault_width_km': self.W_km,
            'fault_depth_km': self.top_depth_km,
            'strike_deg': self.strike_deg,
            'dip_deg': self.dip_deg,
            'rake_deg': self.rake_deg,
            'n_subfaults_along': self.n_along,
            'n_subfaults_down': self.n_down,
            'subfault_size_km': self.subfault_size_km,
            'covariance_model': self.covariance_model,
            'correlation_length_km': self.correlation_length_km,
            'rigidity_pa': self.rigidity,
        }


# =====================================================================
# Legacy convenience functions (for backward compatibility)
# =====================================================================

def load_fault_geometry(magnitude, fault_params_csv=None):
    """
    Load fault geometry for a given magnitude.

    Args:
        magnitude (float): Target magnitude (e.g., 8.5).
        fault_params_csv (Path, optional): Path to CSV. Auto-resolves if None.

    Returns:
        dict: Fault geometry dictionary.
    """
    if fault_params_csv is None:
        fault_params_csv = DATA_DIR / "heidarzadeh_2025_table3.csv"

    fault_data = pd.read_csv(fault_params_csv)
    row = fault_data[
        fault_data["Hypothetical scenario"] == f"Mw {magnitude:.1f}"
    ].iloc[0]

    return {
        'magnitude': magnitude,
        'L_km': row["L (km)"],
        'W_km': row["W (km)"],
        'top_depth_km': row["Top depth (km)"],
        'strike_deg': row["Mean strike ( o )"],
        'dip_deg': row["Mean dip ( o )"],
        'rake_deg': row["Rake ( o )"],
        'slip_mean_m': row["Slip (m)"],
    }


def generate_slip_sample(magnitude, sample_index=0, subfault_size_km=20.0,
                        fault_params_csv=None, output_dir=None, seed=None):
    """
    Generate a single slip sample and optionally save to file.

    Args:
        magnitude (float): Target magnitude.
        sample_index (int): Index in Sobol sequence.
        subfault_size_km (float): Subfault size [km].
        fault_params_csv (Path, optional): CSV file path.
        output_dir (Path, optional): Output directory for saving.
        seed (int, optional): Random seed.

    Returns:
        ndarray: Slip distribution [m], shape (n_down, n_along).
    """
    if fault_params_csv is None:
        fault_params_csv = DATA_DIR / "heidarzadeh_2025_table3.csv"

    if output_dir is None:
        output_dir = OUT_DIR

    generator = StochasticSlipGenerator(
        magnitude,
        fault_params_csv,
        subfault_size_km=subfault_size_km,
    )

    slip = generator.generate_rqmc_sample(sample_index, seed=seed)

    # Optionally save
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"slip_Mw{magnitude}_sample{sample_index}.npy"
        np.save(output_path, slip)

    return slip


# =====================================================================
# MAIN: Example usage
# =====================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("RQMC Slip Generation Example")
    print("=" * 70)
    print()

    # Create generator for Mw 8.5
    Mw_target = 8.5
    generator = StochasticSlipGenerator(
        magnitude=Mw_target,
        fault_params_csv=DATA_DIR / "heidarzadeh_2025_table3.csv",
        subfault_size_km=20.0,
        covariance_model='exponential',
        correlation_length_km=20.0,
    )

    print("Generator configuration:")
    metadata = generator.get_metadata()
    for key, val in metadata.items():
        if isinstance(val, float):
            print(f"  {key}: {val:.4f}")
        else:
            print(f"  {key}: {val}")
    print()

    # Generate a small ensemble
    n_samples = 3
    print(f"Generating {n_samples} RQMC slip samples...")
    samples = generator.generate_ensemble(n_samples, seed=42)
    print()

    # Save and analyze samples
    for i, slip in enumerate(samples):
        out_path = OUT_DIR / f"slip_Mw{Mw_target}_sample{i}.npy"
        np.save(out_path, slip)

        # Compute moment for validation
        M = generator.rigidity * generator.area_subfault_m2 * np.sum(slip)
        M_err = abs(M - generator.M0_target) / generator.M0_target

        print(f"Sample {i}:")
        print(f"  Saved to: {out_path}")
        print(f"  Mean slip: {np.mean(slip):.3f} m")
        print(f"  Max slip:  {np.max(slip):.3f} m")
        print(f"  Total moment: {M:.3e} N·m")
        print(f"  Moment error: {M_err:.2e}")
        print()

    print("Script completed successfully!")