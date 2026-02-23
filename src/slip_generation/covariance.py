"""
Spatial covariance models for fault slip generation.

This module implements several exponential and Matérn covariance models
commonly used in stochastic rupture mechanics and PTHA frameworks.

References:
    - Mai & Beroza (2002): A spatial random field model for modulating slip
      on a fault surface.
    - Goda et al. (2016): Handbook of seismic risk analysis and management
      of civil infrastructure.
"""

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.special import kv, gamma


def exponential_covariance(distance_matrix, correlation_length, nugget=1e-10):
    """
    Build exponential covariance matrix for slip generation.

    The exponential covariance model is frequently used in PTHA because it:
    (1) represents realistic spatial heterogeneity of slip,
    (2) has a closed-form Cholesky decomposition via FFT for large grids,
    (3) is computationally efficient for moderate grid sizes.

    Covariance: C(r) = σ² * exp(-r / ξ)

    Args:
        distance_matrix (ndarray): Pairwise distance matrix, shape (N, N).
        correlation_length (float): Correlation length ξ (km or user units).
        nugget (float, optional): Small diagonal regularization for numerical stability.
                                 Default: 1e-10.

    Returns:
        ndarray: Positive-definite covariance matrix, shape (N, N).

    Notes:
        - The nugget is added to diagonal after exponentiation to ensure PD.
        - Variance is normalized to σ² = 1 (i.e., C(0) = 1).
    """
    C = np.exp(-distance_matrix / correlation_length)
    C += nugget * np.eye(C.shape[0])
    
    # Ensure positive definiteness via symmetric regularization
    eigvals = np.linalg.eigvalsh(C)
    if np.min(eigvals) < 1e-12:
        C += (1e-10 - np.min(eigvals)) * np.eye(C.shape[0])
    
    return C


def matern_covariance(distance_matrix, correlation_length, smoothness=0.5, nugget=1e-10):
    """
    Build Matérn covariance matrix for slip generation.

    The Matérn family generalizes exponential and Gaussian covariance models
    via the smoothness parameter ν. Smaller ν yields rougher realizations;
    larger ν yields smoother realizations.

    Covariance: C(r) = (2^(1-ν) / Γ(ν)) * (√(2ν) * r / ξ)^ν * K_ν(√(2ν) * r / ξ)

    where K_ν is the modified Bessel function of the second kind.

    Args:
        distance_matrix (ndarray): Pairwise distance matrix, shape (N, N).
        correlation_length (float): Correlation length ξ (km or user units).
        smoothness (float, optional): Smoothness parameter ν. Default: 0.5 (exponential).
                                      ν=0.5 → exponential; ν→∞ → Gaussian.
        nugget (float, optional): Diagonal regularization for stability. Default: 1e-10.

    Returns:
        ndarray: Positive-definite covariance matrix, shape (N, N).

    Notes:
        - At ν=0.5, reduces to exponential covariance.
        - Commonly used: ν ∈ {0.3, 0.5, 1.0, 1.5, 2.5}.
        - For seismic slip, ν ≈ 0.3–0.5 often provides physically realistic roughness.
    """
    def matern_kernel(r, ell, nu):
        """Evaluate Matérn kernel, handling r=0 separately."""
        if isinstance(r, np.ndarray):
            result = np.zeros_like(r, dtype=float)
            nonzero = r > 0
            r_scaled = np.sqrt(2 * nu) * r[nonzero] / ell
            coeff = (2**(1 - nu)) / gamma(nu)
            result[nonzero] = coeff * (r_scaled**nu) * kv(nu, r_scaled + 1e-12)
            result[~nonzero] = 1.0
            return result
        else:
            if r == 0:
                return 1.0
            r_scaled = np.sqrt(2 * smoothness) * r / correlation_length
            coeff = (2**(1 - smoothness)) / gamma(smoothness)
            return coeff * (r_scaled**smoothness) * kv(smoothness, r_scaled + 1e-12)
    
    # Vectorized evaluation for off-diagonal elements
    C = np.zeros_like(distance_matrix, dtype=float)
    for i in range(distance_matrix.shape[0]):
        for j in range(distance_matrix.shape[1]):
            C[i, j] = matern_kernel(distance_matrix[i, j], correlation_length, smoothness)
    
    C += nugget * np.eye(C.shape[0])
    
    # Ensure positive definiteness
    eigvals = np.linalg.eigvalsh(C)
    if np.min(eigvals) < 1e-12:
        C += (1e-10 - np.min(eigvals)) * np.eye(C.shape[0])
    
    return C


def build_distance_matrix(coordinates):
    """
    Compute pairwise Euclidean distance matrix from coordinates.

    Args:
        coordinates (ndarray): Array of shape (N, 2) or (N, 3) representing
                              subfault center locations (e.g., in km).

    Returns:
        ndarray: Pairwise distance matrix of shape (N, N).
    """
    return squareform(pdist(coordinates, metric='euclidean'))


def build_covariance_matrix(coordinates, model='exponential', **kwargs):
    """
    Unified interface to build a covariance matrix given coordinates.

    Args:
        coordinates (ndarray): Subfault locations, shape (N, 2) or (N, 3).
        model (str, optional): Covariance model type: 'exponential' or 'matern'.
                              Default: 'exponential'.
        **kwargs: Keyword arguments passed to the covariance function:
                 - correlation_length (float): Characteristic length scale.
                 - smoothness (float): Matérn smoothness parameter ν.
                 - nugget (float): Diagonal regularization.

    Returns:
        ndarray: Positive-definite covariance matrix, shape (N, N).

    Raises:
        ValueError: If model is not recognized.
    """
    dist_matrix = build_distance_matrix(coordinates)
    
    if model.lower() == 'exponential':
        return exponential_covariance(dist_matrix, **kwargs)
    elif model.lower() == 'matern':
        return matern_covariance(dist_matrix, **kwargs)
    else:
        raise ValueError(f"Unknown covariance model: {model}. "
                        "Choose 'exponential' or 'matern'.")


if __name__ == "__main__":
    # Demonstration: Build and verify covariance matrices
    np.random.seed(42)
    
    # Create a 3x3 grid of subfaults
    x = np.linspace(0, 60, 3)
    y = np.linspace(0, 60, 3)
    X, Y = np.meshgrid(x, y)
    coords = np.column_stack([X.ravel(), Y.ravel()])
    
    print("Demonstrating covariance models for a 3×3 subfault grid:")
    print(f"Coordinates shape: {coords.shape}")
    print(f"Coordinates:\n{coords}\n")
    
    # Exponential covariance
    C_exp = build_covariance_matrix(coords, model='exponential', 
                                     correlation_length=20.0)
    print("Exponential covariance matrix (ξ=20 km):")
    print(C_exp)
    print(f"Eigenvalues: {np.linalg.eigvalsh(C_exp)}\n")
    
    # Matérn covariance
    C_mat = build_covariance_matrix(coords, model='matern', 
                                     correlation_length=20.0, 
                                     smoothness=0.5)
    print("Matérn covariance matrix (ξ=20 km, ν=0.5):")
    print(C_mat)
    print(f"Eigenvalues: {np.linalg.eigvalsh(C_mat)}\n")
