"""
Inundation mapper for Philippine tsunami.
Takes Okada displacement field and propagates to determine coastal inundation.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from scipy.interpolate import RegularGridInterpolator
import json
import urllib.request

OUTPUT_DIR = Path(__file__).parent.parent.parent / "output"
DATA_DIR = Path(__file__).parent.parent.parent / "data"
BATHY_DIR = DATA_DIR / "bathymetry"

def load_okada_displacement(magnitude="Mw8.5", n_samples=None):
    """
    Load Okada displacement fields from phase3_okada output.
    Computes mean displacement across samples.
    
    Parameters:
    - magnitude: which magnitude to load (e.g., "Mw8.5")
    - n_samples: how many samples to load (None = load all)
    
    Returns: (mean_displacement, lon, lat)
    """
    try:
        okada_dir = OUTPUT_DIR / "tsunami_sources_okada" / magnitude
        if not okada_dir.exists():
            print(f"ERROR: {okada_dir} not found")
            return None, None, None
        
        # Get list of displacement files
        disp_files = sorted(okada_dir.glob("*_displacement_okada.npy"))
        if not disp_files:
            print(f"No displacement files found in {okada_dir}")
            return None, None, None
        
        if n_samples:
            disp_files = disp_files[:n_samples]
        
        print(f"Loading {len(disp_files)} displacement samples from {magnitude}...")
        
        # Load and stack displacements
        displacements = []
        for f in disp_files:
            data = np.load(f, allow_pickle=True).item()
            displacements.append(data['displacement'])
        
        displacement_array = np.array(displacements)
        mean_disp = displacement_array.mean(axis=0)
        
        # Also load metadata (lon/lat grid)
        if displacements:
            metadata = np.load(disp_files[0], allow_pickle=True).item()
            lon = metadata.get('lon')
            lat = metadata.get('lat')
        else:
            lon, lat = None, None
        
        print(f"✓ Loaded {len(displacements)} displacements")
        print(f"  Mean displacement shape: {mean_disp.shape}")
        print(f"  Max displacement: {mean_disp.max():.2f} m, Min: {mean_disp.min():.2f} m")
        
        return mean_disp, lon, lat
        
    except Exception as e:
        print(f"Could not load Okada displacement: {e}")
        return None, None, None


def fetch_coastline_geojson():
    """Fetch Philippines coastline from Natural Earth."""
    try:
        url = (
            "https://raw.githubusercontent.com/datasets/geo-boundaries-world-110m/"
            "master/countries.geojson"
        )
        print("Fetching Philippines coastline...")
        with urllib.request.urlopen(url, timeout=10) as f:
            gj = json.load(f)
        
        for feat in gj['features']:
            props = feat.get('properties', {})
            name = props.get('admin') or props.get('name') or ''
            if name.lower() == 'philippines':
                geom = feat['geometry']
                return geom
        
        print("WARNING: Philippines not found in geojson")
        return None
    except Exception as e:
        print(f"Could not fetch coastline: {e}")
        return None


def create_synthetic_okada_field(lon, lat):
    """
    Create synthetic Okada displacement field for testing.
    Based on a Gaussian source near Philippine Trench.
    """
    # Trench source location
    source_lon, source_lat = 131.0, 11.5
    source_depth_equiv = 50  # km equivalent
    
    # Create displacement field: Gaussian falloff from source
    LON, LAT = np.meshgrid(lon, lat)
    
    # Distance from source
    dist_km = 111 * np.sqrt((LAT - source_lat)**2 + (LON - source_lon)**2 / np.cos(np.radians(source_lat))**2)
    
    # Displacement amplitude (m)
    # Maximum ~1 m at source, decays with distance
    disp = 1.0 * np.exp(-(dist_km / 50)**2)
    
    # Add directional component (source dip-slip motion)
    # Higher displacement to the west (landward)
    direction = np.cos(np.pi * (LON - source_lon) / 5)
    disp = disp * np.maximum(direction, 0.1)
    
    return disp


def propagate_tsunami_ray_tracing(bathymetry, okada_disp, bathymetry_lon, bathymetry_lat):
    """
    Simplified tsunami propagation using ray-tracing / wavefront approach.
    Much more stable than shallow water equations for this demonstration.
    
    Key idea: tsunami travels at sqrt(g*h), spreading from source with amplitude decay
    based on distance and bathymetry refraction.
    
    Parameters:
    - bathymetry: 2D array of water depths (negative)
    - okada_disp: 2D array of initial displacement (m)
    - bathymetry_lon, bathymetry_lat: coordinate arrays
    
    Returns: water height field over ocean
    """
    print("Propagating tsunami using ray-tracing approach...")
    g = 9.81  # m/s^2
    
    ny, nx = okada_disp.shape
    LON, LAT = np.meshgrid(bathymetry_lon, bathymetry_lat)
    
    # Find source center (maximum displacement)
    src_idx = np.unravel_index(np.argmax(okada_disp), okada_disp.shape)
    src_lat = LAT[src_idx]
    src_lon = LON[src_idx]
    
    print(f"  Source location: ({src_lon:.2f}°, {src_lat:.2f}°)")
    
    # Water depth (convert to positive)
    h = np.abs(bathymetry)
    h[h < 1] = 1  # Minimum depth = 1m
    
    # Wave speed sqrt(g*h)
    c = np.sqrt(g * h)
    
    # Distance from source (in degrees, approximate)
    dx_deg = bathymetry_lon[1] - bathymetry_lon[0] if len(bathymetry_lon) > 1 else 0.1
    dy_deg = bathymetry_lat[1] - bathymetry_lat[0] if len(bathymetry_lat) > 1 else 0.1
    dx_km = dx_deg * 111 * np.cos(np.radians(LAT))
    dy_km = dy_deg * 111
    
    dist_km = 111 * np.sqrt(
        (LAT - src_lat)**2 +
        ((LON - src_lon) * np.cos(np.radians((LAT + src_lat)/2)))**2
    )
    
    # Travel time from source (in hours)
    travel_time = np.zeros_like(dist_km)
    travel_time[c > 0] = dist_km[c > 0] / (c[c > 0] / 3.6)  # convert m/s to km/h
    
    # Water elevation: initial displacement decays with distance and refracts by bathymetry
    # amplitude = amplitude_source * sqrt(c_source / c_point) * exp(-t / tau)
    
    # Use simplified version: amplitude decays as 1/sqrt(distance)
    # and attenuates through shallow water
    water_height = okada_disp.copy()
    
    # Amplitude decay: 1/sqrt(distance) after traveling
    decay = np.ones_like(dist_km)
    decay[dist_km > 5] = 5 / np.sqrt(dist_km[dist_km > 5])
    
    # Shallow water damping: waves dissipate in shallow water faster
    manning = 0.02  # Manning roughness
    time_hours = 1.0  # 1 hour propagation
    damping = np.exp(-manning * time_hours / (h / 1000))
    
    water_height = okada_disp * decay * damping
    
    # Refraction: waves slow and amplify when entering shallow water
    # Use Green's law: amplitude ~ 1 / (4th root of depth)
    h_source = h[src_idx]
    refraction = (h_source / h) ** 0.25
    water_height = water_height * refraction
    
    # Cap at reasonable values
    water_height = np.clip(water_height, -10, 10)
    
    return water_height


def determine_inundation(water_height, dem, dem_lon, dem_lat):
    """
    Determine inundated areas by comparing water height to topography.
    
    Parameters:
    - water_height: 2D array of water surface elevation (m above MSL), already on DEM grid
    - dem: 2D array of topography (m)
    - dem_lon, dem_lat: coordinate arrays for DEM (not used, but kept for API compatibility)
    
    Returns inundation map (0 = dry, 1 = wet)
    """
    print("Determining inundation extent...")
    
    # Ensure same shape
    assert water_height.shape == dem.shape, f"Shape mismatch: water {water_height.shape} vs dem {dem.shape}"
    
    # Inundation occurs where water surface > topography
    inundation = (water_height > dem).astype(float)
    
    # Compute inundation depth (only where inundated)
    inundation_depth = water_height - dem
    inundation_depth[dem >= water_height] = 0
    
    return inundation, inundation_depth


def plot_inundation_map(okada_disp, dem, dem_lon, dem_lat, okada_lon, okada_lat, coastline_geom, 
                        inundation=None, save_path=None):
    """
    Plot inundation map with Philippines coastline.
    
    Parameters:
    - okada_disp: 2D array of displacement (m) on okada grid
    - dem: 2D array of topography (m) on DEM grid
    - dem_lon, dem_lat: coordinate arrays for DEM
    - okada_lon, okada_lat: coordinate arrays for Okada
    - coastline_geom: GeoJSON geometry for coastline
    - inundation: optional inundation mask (on DEM grid)
    - save_path: where to save the figure
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 10))
    
    # --- Panel 1: Displacement field ---
    ax = axes[0]
    im1 = ax.contourf(okada_lon, okada_lat, okada_disp, levels=20, cmap='RdYlBu_r', vmin=-2, vmax=2)
    ax.set_title('Initial Tsunami Displacement (Okada)', fontsize=14, fontweight='bold')
    ax.set_xlabel('Longitude (°)')
    ax.set_ylabel('Latitude (°)')
    
    # Overlay coastline
    if coastline_geom:
        try:
            coords = coastline_geom.get('coordinates', [])
            if coastline_geom.get('type') == 'Polygon':
                coords = [coords]
            elif coastline_geom.get('type') == 'MultiPolygon':
                coords = [c[0] for c in coords]  # Just outer rings
            
            for ring in coords:
                if ring:
                    ring = np.array(ring)
                    ax.plot(ring[:, 0], ring[:, 1], 'k-', linewidth=0.5, alpha=0.7)
        except Exception as e:
            print(f"Could not plot coastline: {e}")
    
    cb1 = plt.colorbar(im1, ax=ax)
    cb1.set_label('Displacement (m)')
    ax.grid(True, alpha=0.3)
    
    # --- Panel 2: Inundation map ---
    ax = axes[1]
    
    # Show topography
    im2 = ax.contourf(dem_lon, dem_lat, dem, levels=20, cmap='terrain', alpha=0.7)
    
    # Overlay inundation
    if inundation is not None:
        # Only show inundated areas
        inundation_masked = np.ma.masked_where(inundation < 0.5, inundation)
        ax.contourf(dem_lon, dem_lat, inundation_masked, levels=[0.5, 1.5], 
                   colors=['cyan'], alpha=0.6)
    
    ax.set_title('Predicted Inundation', fontsize=14, fontweight='bold')
    ax.set_xlabel('Longitude (°)')
    ax.set_ylabel('Latitude (°)')
    
    # Overlay coastline
    if coastline_geom:
        try:
            coords = coastline_geom.get('coordinates', [])
            if coastline_geom.get('type') == 'Polygon':
                coords = [coords]
            elif coastline_geom.get('type') == 'MultiPolygon':
                coords = [c[0] for c in coords]
            
            for ring in coords:
                if ring:
                    ring = np.array(ring)
                    ax.plot(ring[:, 0], ring[:, 1], 'k-', linewidth=1.0, alpha=0.8)
        except Exception as e:
            print(f"Could not plot coastline: {e}")
    
    cb2 = plt.colorbar(im2, ax=ax)
    cb2.set_label('Elevation (m)')
    ax.grid(True, alpha=0.3)
    
    # Add legend
    wet_patch = mpatches.Patch(color='cyan', alpha=0.6, label='Inundated')
    ax.legend(handles=[wet_patch], loc='upper left')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"✓ Map saved: {save_path}")
    else:
        plt.show()


def main():
    print("=" * 70)
    print("INUNDATION MAPPER - PHILIPPINE TSUNAMI")
    print("=" * 70)
    
    # Load or download bathymetry & DEM
    print("\n1. Loading bathymetry and DEM...")
    try:
        from download_bathymetry_dem import load_bathymetry_and_dem
        data = load_bathymetry_and_dem()
        
        bathy_lon = data['bathymetry']['lon']
        bathy_lat = data['bathymetry']['lat']
        bathymetry = data['bathymetry']['elevation']
        
        dem_lon = data['dem']['lon']
        dem_lat = data['dem']['lat']
        dem = data['dem']['elevation']
        
        print(f"✓ Bathymetry: {bathymetry.shape}, lon: {len(bathy_lon)}, lat: {len(bathy_lat)}")
        print(f"✓ DEM: {dem.shape}, lon: {len(dem_lon)}, lat: {len(dem_lat)}")
        
    except Exception as e:
        print(f"ERROR loading data: {e}")
        print("Using synthetic data instead...")
        
        bathy_lon = np.linspace(120, 135, 150)
        bathy_lat = np.linspace(3, 20, 200)
        bathymetry = -np.ones((len(bathy_lat), len(bathy_lon))) * 3000
        
        dem_lon = np.linspace(120, 135, 300)
        dem_lat = np.linspace(3, 20, 400)
        dem = np.ones((len(dem_lat), len(dem_lon))) * 500
    
    # Generate Okada displacement on bathymetry grid
    print("\n2. Computing Okada discharge field...")
    print("  Creating synthetic Okada displacements...")
    okada_lon = bathy_lon
    okada_lat = bathy_lat
    okada_disp = create_synthetic_okada_field(okada_lon, okada_lat)
    print(f"✓ Okada displacement: {okada_disp.shape}, max = {okada_disp.max():.2f} m")
    
    # Propagate tsunami
    print("\n3. Propagating tsunami through ocean...")
    water_height = propagate_tsunami_ray_tracing(bathymetry, okada_disp, bathy_lon, bathy_lat)
    print(f"✓ Final water height: max = {water_height.max():.2f} m, min = {water_height.min():.2f} m")
    
    # Interpolate water height to DEM grid
    print("\n4. Interpolating water field to DEM grid...")
    from scipy.interpolate import RegularGridInterpolator
    
    # Create interpolator for water height
    interp_func = RegularGridInterpolator(
        (bathy_lat, bathy_lon),
        water_height,
        bounds_error=False,
        fill_value=0.0,
        method='linear'
    )
    
    # Evaluate on DEM grid
    DEM_LON, DEM_LAT = np.meshgrid(dem_lon, dem_lat)
    pts = np.column_stack([DEM_LAT.ravel(), DEM_LON.ravel()])
    water_height_on_dem = interp_func(pts).reshape(dem.shape)
    print(f"✓ Interpolated to DEM grid: {water_height_on_dem.shape}")
    
    # Determine inundation
    print("\n5. Computing inundation extent...")
    inundation, inundation_depth = determine_inundation(water_height_on_dem, dem, dem_lon, dem_lat)
    inundated_area_pct = inundation.sum() / inundation.size * 100
    max_inund = np.nanmax(inundation_depth[inundation_depth > 0]) if np.any(inundation_depth > 0) else 0
    print(f"✓ Inundated area: {inundated_area_pct:.2f}% of domain")
    print(f"✓ Max inundation depth: {max_inund:.2f} m")
    
    # Fetch coastline
    print("\n6. Fetching Philippines coastline...")
    coastline_geom = fetch_coastline_geojson()
    if coastline_geom:
        print("✓ Coastline loaded")
    
    # Create plot
    print("\n7. Generating inundation map...")
    output_path = OUTPUT_DIR / "inundation_map_philippines.png"
    plot_inundation_map(okada_disp, dem, dem_lon, dem_lat, okada_lon, okada_lat, coastline_geom, 
                        inundation=inundation, save_path=output_path)
    
    print("\n" + "=" * 70)
    print("INUNDATION MAPPING COMPLETE")
    print("=" * 70)
    print(f"\n✓ Output: {output_path}")
    print(f"\nSummary (Synthetic Okada):")
    print(f"  - Max initial displacement: {okada_disp.max():.2f} m")
    print(f"  - Max water elevation on ocean: {water_height.max():.2f} m")
    print(f"  - Max inundation depth on land: {max_inund:.2f} m")
    print(f"  - Inundated area: {inundated_area_pct:.2f}%")


if __name__ == "__main__":
    main()
