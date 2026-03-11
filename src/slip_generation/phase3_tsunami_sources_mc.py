"""
Phase 3 (MC): Tsunami Source Generation from MC slip samples

This script mirrors phase3_tsunami_sources.py but uses MC-generated slip
samples from phase1_generate_mc.py instead of RQMC samples.

Usage:
    python phase3_tsunami_sources_mc.py
"""

import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))

from tsunami_source import SimpleDislocationSource
from philippine_trench_config import PHILIPPINE_TRENCH_PTHA_CONFIG, RUPTURE_SCENARIOS


def build_fault_geometry(scenario):
    """Convert RuptureScenario to fault_geometry dict."""
    L = scenario.length_km
    W = scenario.width_km
    target_n = PHILIPPINE_TRENCH_PTHA_CONFIG.slip_params.n_subfaults
    subfault_size = np.sqrt((L * W) / target_n)
    n_along = int(np.round(L / subfault_size))
    n_down = int(np.round(W / subfault_size))
    
    return {
        'length_km': L,
        'width_km': W,
        'top_depth_km': scenario.top_depth_km,
        'strike_deg': scenario.strike_deg,
        'dip_deg': scenario.dip_deg,
        'rake_deg': scenario.rake_deg,
        'n_subfaults_along': n_along,
        'n_subfaults_down': n_down,
    }


def create_observation_grid(domain_lon=(125.5, 130.5), domain_lat=(6.0, 13.0),
                            n_lon=20, n_lat=20):
    """Create a regular grid of observation points."""
    lons = np.linspace(domain_lon[0], domain_lon[1], n_lon)
    lats = np.linspace(domain_lat[0], domain_lat[1], n_lat)
    LON, LAT = np.meshgrid(lons, lats)
    points = np.column_stack([LON.ravel(), LAT.ravel()])
    return points, lons, lats


def main():
    cfg = PHILIPPINE_TRENCH_PTHA_CONFIG
    
    print("Phase 3 (MC): Tsunami Source Generation")
    print("========================================")
    print(f"Using MC slip samples from Phase 1 (MC)\n")
    
    # Where MC slip samples are stored
    slip_root = Path('output') / 'slip_samples_mc'
    
    # Where to save MC displacement results
    disp_root = Path('output') / 'tsunami_sources_mc'
    disp_root.mkdir(parents=True, exist_ok=True)
    
    # Create observation grid
    obs_points, obs_lons, obs_lats = create_observation_grid()
    print(f"Observation grid: {len(obs_points)} points")
    print(f"  Lon: {obs_lons.min():.1f}–{obs_lons.max():.1f}°")
    print(f"  Lat: {obs_lats.min():.1f}–{obs_lats.max():.1f}°\n")
    
    # Process each magnitude scenario
    scenarios_to_process = [s for s in RUPTURE_SCENARIOS.values()
                            if s.magnitude in [8.0, 8.5, 9.0]]
    
    for scenario in scenarios_to_process:
        mag = scenario.magnitude
        mag_str = f"Mw{mag:.1f}"
        
        slip_dir = slip_root / mag_str
        if not slip_dir.exists():
            print(f"Skipping {mag_str} — slip directory not found: {slip_dir}")
            continue
        
        # Load all MC samples
        slip_files = sorted(slip_dir.glob('*.npy'))
        slip_files_subset = slip_files
        
        print(f"Processing {mag_str}: {len(slip_files_subset)} MC samples")
        
        fault_geom = build_fault_geometry(scenario)
        print(f"  Fault geometry: L={fault_geom['length_km']:.0f} km, "
              f"W={fault_geom['width_km']:.0f} km")
        
        # Create output directory
        out_mag_dir = disp_root / mag_str
        out_mag_dir.mkdir(parents=True, exist_ok=True)
        
        max_disps = []
        
        for slip_file in slip_files_subset:
            slip = np.load(slip_file)
            
            source = SimpleDislocationSource(
                fault_geometry=fault_geom,
                slip_distribution=slip,
                rigidity=4.0e10,
            )
            
            displacement = source.get_seafloor_displacement(obs_points)
            max_disp = np.max(np.abs(displacement))
            max_disps.append(max_disp)
            
            out_path = out_mag_dir / slip_file.name.replace('.npy', '_displacement.npy')
            np.save(out_path, displacement)
        
        max_disps = np.array(max_disps)
        print(f"  Generated {len(max_disps)} displacement fields")
        print(f"    Max displacement range: {max_disps.min():.3e}–{max_disps.max():.3e} m")
        print(f"    Mean max displacement: {max_disps.mean():.3e} m\n")
    
    print("Phase 3 (MC) complete. Displacement fields saved to output/tsunami_sources_mc/")


if __name__ == '__main__':
    main()
