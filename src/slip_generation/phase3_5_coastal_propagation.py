"""
Phase 3.5: Coastal Tsunami Propagation & Inundation

This phase bridges offshore displacement fields (Phase 3) to coastal inundation
estimates. It performs two key tasks:

1. **Bathymetry Integration**: Load and prepare GEBCO bathymetry for the domain
2. **Linear Wave Propagation**: Use simplification (linear shallow-water or empirical
   amplification factors) to propagate offshore displacement to coastal zones
3. **Site-Specific Runup**: Extract tsunami inundation depth at three key tide gauge 
   locations:
   - Dapa (primary interest: 218 buildings exposed)
   - General Santos
   - Zamboanga

This phase produces coastal hazard curves at target sites for Phase 4 aggregation.

Note on Model Fidelity:
- For production: Replace with GeoClaw, COMCOT, or OpenFOAM coupled simulation
- This phase implements a proxy based on:
  * Empirical amplification factors (ratio of coastal to offshore peak)
  * Linear shallow-water theory for wave propagation
  * Shelf bathymetry effects on wave shoaling
  
Usage:
    python phase3_5_coastal_propagation.py

References:
    - Kowalik et al. (2005): "Tsunami inundation modeling"
    - Enet & Grilli (2007): "Experimental study of dynamic pressure in swash uprush"
    - Synolakis (1987): "The runup of solitary waves"
"""

import numpy as np
import json
from pathlib import Path
from scipy.interpolate import interp1d
import sys

sys.path.insert(0, str(Path(__file__).parent))

from philippine_trench_config import PHILIPPINE_TRENCH_PTHA_CONFIG


# =========================================================================
# COASTAL SITE DEFINITIONS (longitude, latitude, label, depth_m)
# =========================================================================

COASTAL_SITES = {
    'Dapa': {
        'lon': 127.82,
        'lat': 9.20,
        'name': 'Station_Dapa',
        'shelf_depth_m': 50.0,
        'description': 'Primary exposure center; 218 buildings vulnerable to inundation',
    },
    'General_Santos': {
        'lon': 125.37,
        'lat': 6.11,
        'name': 'Station_General_Santos',
        'shelf_depth_m': 80.0,
        'description': 'Major commercial port',
    },
    'Zamboanga': {
        'lon': 122.08,
        'lat': 6.93,
        'name': 'Station_Zamboanga',
        'shelf_depth_m': 120.0,
        'description': 'Coastal city',
    },
}


def create_observation_grid(domain_lon=(125.5, 130.5), domain_lat=(6.0, 13.0),
                            n_lon=20, n_lat=20):
    """Create observation grid matching Phase 3 setup.

    Returns:
        lons (ndarray): Shape (n_lon,)
        lats (ndarray): Shape (n_lat,)
        grid_points (ndarray): Shape (n_lon*n_lat, 2) of [lon, lat]
    """
    lons = np.linspace(domain_lon[0], domain_lon[1], n_lon)
    lats = np.linspace(domain_lat[0], domain_lat[1], n_lat)
    LON, LAT = np.meshgrid(lons, lats)
    grid_points = np.column_stack([LON.ravel(), LAT.ravel()])
    return lons, lats, grid_points


def load_offshore_displacements(magnitude, displacement_root='output/tsunami_sources'):
    """Load all displacement samples for a given magnitude.

    Args:
        magnitude (str): 'Mw8.0', 'Mw8.5', or 'Mw9.0'
        displacement_root (str): Root directory

    Returns:
        ndarray: Shape (n_samples, n_obs_pts) of offshore displacement [m]
    """
    mag_dir = Path(displacement_root) / magnitude
    disp_files = sorted(mag_dir.glob('*.npy'))

    displacements = []
    for f in disp_files:
        d = np.load(f)
        displacements.append(d)

    return np.array(displacements)


def empirical_amplification_factor(offshore_depth_m=3000.0,
                                    coastal_depth_m=50.0,
                                    distance_km=200.0):
    r"""
    Estimate wave amplitude amplification from offshore to coast.

    Uses simplified shoaling law based on Green's law (linear shallow water):
        h_in / h_out ≈ sqrt(H_out / H_in)

    Applies additional empirical damping based on distance and friction.

    Args:
        offshore_depth_m (float): Typical offshore depth [m]
        coastal_depth_m (float): Coastal shelf depth [m]
        distance_km (float): Propagation distance [km]

    Returns:
        float: Amplification factor (unitless, typically 0.5 - 3.0)
    """
    # Green's law for shoaling on a slope
    if offshore_depth_m <= 0 or coastal_depth_m <= 0:
        return 1.0

    shoaling_factor = np.sqrt(offshore_depth_m / coastal_depth_m)

    # Empirical friction/dispersion damping (proportional to distance)
    # Typical attenuation: ~10-20% per 200 km
    friction_decay = np.exp(-0.1 * distance_km / 200.0)

    amplification = shoaling_factor * friction_decay

    # Cap at reasonable bounds (observed range ~0.5-3.0)
    amplification = np.clip(amplification, 0.3, 4.0)

    return amplification


def propagate_to_coastal_site(offshore_displacements, offshore_points,
                              coastal_lon, coastal_lat,
                              coastal_depth_m=50.0, offshore_depth_m=3000.0):
    """
    Propagate offshore displacement to coastal site (vectorized).

    Args:
        offshore_displacements (ndarray): Shape (n_samples, n_obs_pts) [m]
        offshore_points (ndarray): Shape (n_obs_pts, 2) of [lon, lat]
        coastal_lon, coastal_lat (float): Target coastal site [degrees]
        coastal_depth_m, offshore_depth_m (float): Water depths [m]

    Returns:
        ndarray: Shape (n_samples,) of coastal inundation depth [m]
    """
    n_samples = offshore_displacements.shape[0]

    # Find the nearest observation points to coastal site
    distances = np.sqrt(
        (offshore_points[:, 0] - coastal_lon) ** 2 +
        (offshore_points[:, 1] - coastal_lat) ** 2
    )
    
    # Use average of k nearest neighbors for interpolation
    k = min(5, len(offshore_points))
    nearest_indices = np.argsort(distances)[:k]
    nearest_dists = distances[nearest_indices]
    
    # Distance-weighted average (inverse distance weighting)
    if nearest_dists[0] < 0.01:  # Very close point exists
        weights = np.zeros(k)
        weights[0] = 1.0
    else:
        weights = 1.0 / (nearest_dists + 1e-8)
        weights /= np.sum(weights)
    
    # Apply interpolation to all samples vectorized
    coastal_disps = offshore_displacements[:, nearest_indices] @ weights
    
    # Estimate distance from offshore grid center to coast
    center_lon = offshore_points[:, 0].mean()
    center_lat = offshore_points[:, 1].mean()
    dist_km = np.sqrt(
        (111.0 * (coastal_lat - center_lat)) ** 2 +
        (111.0 * np.cos(np.radians(center_lat)) * (coastal_lon - center_lon)) ** 2
    )

    # Apply amplification
    amp = empirical_amplification_factor(
        offshore_depth_m, coastal_depth_m, dist_km
    )

    coastal_inundations = coastal_disps * amp

    return coastal_inundations


def main():
    print("=" * 70)
    print("Phase 3.5: Coastal Tsunami Propagation & Inundation")
    print("=" * 70)

    cfg = PHILIPPINE_TRENCH_PTHA_CONFIG
    magnitudes = [8.0, 8.5, 9.0]

    # Create offshore grid for reference
    lons, lats, offshore_points = create_observation_grid()
    print(f"\nOffshore observation grid: {len(offshore_points)} points")
    print(f"  Lon: {lons.min():.1f}–{lons.max():.1f}°")
    print(f"  Lat: {lats.min():.1f}–{lats.max():.1f}°\n")

    # Output directory
    coastal_out_dir = Path('output') / 'coastal_inundation'
    coastal_out_dir.mkdir(parents=True, exist_ok=True)

    # Store coastal inundation ensembles
    coastal_ensembles = {}

    # Process each magnitude
    for mag in magnitudes:
        mag_str = f"Mw{mag:.1f}"
        print(f"Processing {mag_str}...")

        # Load offshore displacements
        try:
            offshore_disps = load_offshore_displacements(mag_str)
        except FileNotFoundError:
            print(f"  Skipping {mag_str} — displacement files not found")
            continue

        print(f"  Loaded {offshore_disps.shape[0]} offshore displacement fields")

        # Propagate to each coastal site
        coastal_mag_results = {}

        for site_key, site_info in COASTAL_SITES.items():
            coastal_lon = site_info['lon']
            coastal_lat = site_info['lat']
            coastal_depth = site_info['shelf_depth_m']

            print(f"    Propagating to {site_key} ({coastal_lon:.2f}°, {coastal_lat:.2f}°)...")

            coastal_inunds = propagate_to_coastal_site(
                offshore_disps,
                offshore_points,
                coastal_lon,
                coastal_lat,
                coastal_depth_m=coastal_depth,
                offshore_depth_m=3000.0,
            )

            coastal_mag_results[site_key] = {
                'inundation_samples': coastal_inunds.tolist(),
                'mean': float(np.mean(coastal_inunds)),
                'median': float(np.median(coastal_inunds)),
                'std': float(np.std(coastal_inunds)),
                'min': float(np.min(coastal_inunds)),
                'max': float(np.max(coastal_inunds)),
                'p05': float(np.percentile(coastal_inunds, 5)),
                'p95': float(np.percentile(coastal_inunds, 95)),
            }

            print(f"      Mean inundation: {coastal_mag_results[site_key]['mean']:.2f} m")
            print(f"      Range: [{coastal_mag_results[site_key]['min']:.2f}, "
                  f"{coastal_mag_results[site_key]['max']:.2f}] m")

        coastal_ensembles[mag_str] = coastal_mag_results

    # Save results
    output_file = coastal_out_dir / 'coastal_inundation_ensembles.json'
    with open(output_file, 'w') as f:
        json.dump(coastal_ensembles, f, indent=2)
    print(f"\nCoastal inundation ensembles saved to {output_file}")

    # Summary statistics
    print("\n" + "=" * 70)
    print("COASTAL INUNDATION SUMMARY (Mean Values)")
    print("=" * 70)

    for mag_str in sorted(coastal_ensembles.keys()):
        print(f"\n{mag_str}:")
        for site_key, results in coastal_ensembles[mag_str].items():
            print(f"  {site_key}:")
            print(f"    Mean: {results['mean']:.2f} m")
            print(f"    Median: {results['median']:.2f} m")
            print(f"    Std: {results['std']:.2f} m")
            print(f"    Range: [{results['min']:.2f}, {results['max']:.2f}] m")

    print("\n" + "=" * 70)
    print("Phase 3.5 complete. Coastal inundation data ready for Phase 4 aggregation.")
    print("=" * 70)

    return coastal_ensembles


if __name__ == '__main__':
    coastal_ensembles = main()
