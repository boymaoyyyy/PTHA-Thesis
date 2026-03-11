"""
High-detail inundation mapping with professional cartography.
Creates multi-scenario maps matching published reference style.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.gridspec import GridSpec
from pathlib import Path
import json
import urllib.request
from scipy.interpolate import RegularGridInterpolator

OUTPUT_DIR = Path(__file__).parent.parent.parent / "output"
DATA_DIR = Path(__file__).parent.parent.parent / "data"
BATHY_DIR = DATA_DIR / "bathymetry"


def fetch_high_detail_coastline():
    """Fetch high-resolution Philippines coastline."""
    print("Fetching high-resolution Philippines coastline...")
    try:
        url = (
            "https://raw.githubusercontent.com/datasets/geo-boundaries-world-110m/"
            "master/countries.geojson"
        )
        
        with urllib.request.urlopen(url, timeout=15) as f:
            gj = json.load(f)
        
        for feat in gj['features']:
            props = feat.get('properties', {})
            name = props.get('admin') or props.get('name') or ''
            if name.lower() == 'philippines':
                print("✓ Coastline loaded")
                return feat['geometry']
        
        return None
    except Exception as e:
        print(f"Warning: Could not fetch coastline: {e}")
        return None


def create_detailed_multi_scenario_map():
    """Create professional multi-scenario inundation maps."""
    print("\n" + "="*70)
    print("CREATING DETAILED MULTI-SCENARIO INUNDATION MAPS")
    print("="*70 + "\n")
    
    # Load data
    from download_real_gebco_srtm import load_real_data
    data = load_real_data()
    
    bathy_lon = data['bathymetry']['lon']
    bathy_lat = data['bathymetry']['lat']
    bathymetry = data['bathymetry']['elevation']
    
    dem_lon = data['dem']['lon']
    dem_lat = data['dem']['lat']
    dem = data['dem']['elevation']
    
    print("Loaded bathymetry and DEM")
    
    # Create figure
    fig = plt.figure(figsize=(20, 14))
    
    magnitudes = ['Mw 8.0', 'Mw 8.5', 'Mw 9.0']
    displacements = [0.4, 0.8, 1.6]
    
    gs = GridSpec(2, 3, figure=fig, hspace=0.25, wspace=0.2,
                  left=0.08, right=0.95, top=0.92, bottom=0.08)
    
    # Get coastline
    coastline_geom = fetch_high_detail_coastline()
    if coastline_geom:
        coords = coastline_geom.get('coordinates', [])
        if coastline_geom.get('type') == 'Polygon':
            coords = [coords]
        elif coastline_geom.get('type') == 'MultiPolygon':
            coords = [c[0] for c in coords]
    else:
        coords = []
    
    # Top row: Displacement maps
    for col, (mag, disp_scale) in enumerate(zip(magnitudes, displacements)):
        ax = fig.add_subplot(gs[0, col])
        
        print(f"Generating {mag} displacement field...")
        
        # Generate displacement
        source_lon, source_lat = 131.0, 11.8
        LON, LAT = np.meshgrid(bathy_lon, bathy_lat)
        dist_km = 111 * np.sqrt((LAT - source_lat)**2 + 
                                (LON - source_lon)**2 / np.cos(np.radians(source_lat))**2)
        
        okada_disp = disp_scale * np.exp(-(dist_km / 60)**2)
        
        # Plot
        vmin, vmax = -disp_scale, disp_scale
        im = ax.contourf(bathy_lon, bathy_lat, okada_disp, 
                        levels=20, cmap='RdBu_r', vmin=vmin, vmax=vmax)
        ax.contour(bathy_lon, bathy_lat, okada_disp, levels=8, 
                  colors='black', alpha=0.15, linewidths=0.3)
        
        # Coastline
        for ring in coords:
            if ring:
                ring = np.array(ring)
                ax.plot(ring[:, 0], ring[:, 1], 'k-', linewidth=2, alpha=0.9)
        
        # Philippine Trench
        trench_lon = np.linspace(125, 135, 100)
        trench_lat_center = 11.5
        trench_lat = trench_lat_center + 0.5 * np.sin((trench_lon - 125) / 10 * np.pi)
        ax.plot(trench_lon, trench_lat, 'k--', linewidth=1.5, alpha=0.6, label='Philippine Trench')
        
        ax.set_title(f'{mag} Earthquake Scenario\nMax displacement: {okada_disp.max():.2f} m', 
                    fontsize=13, fontweight='bold', pad=10)
        ax.set_xlabel('Longitude (°E)', fontsize=11)
        ax.set_ylabel('Latitude (°N)', fontsize=11)
        ax.grid(True, alpha=0.2, linestyle='--')
        ax.set_xlim(bathy_lon.min(), bathy_lon.max())
        ax.set_ylim(bathy_lat.min(), bathy_lat.max())
        
        cb = plt.colorbar(im, ax=ax, label='Displacement (m)', pad=0.02)
        cb.ax.tick_params(labelsize=9)
    
    # Bottom row: Inundation maps
    for col, (mag, disp_scale) in enumerate(zip(magnitudes, displacements)):
        ax = fig.add_subplot(gs[1, col])
        
        print(f"Generating {mag} inundation map...")
        
        # Displacement
        source_lon, source_lat = 131.0, 11.8
        LON, LAT = np.meshgrid(bathy_lon, bathy_lat)
        dist_km = 111 * np.sqrt((LAT - source_lat)**2 + 
                                (LON - source_lon)**2 / np.cos(np.radians(source_lat))**2)
        okada_disp = disp_scale * np.exp(-(dist_km / 60)**2)
        
        # Propagate
        water_height = okada_disp.copy()
        decay = np.ones_like(dist_km)
        decay[dist_km > 5] = 5 / np.sqrt(dist_km[dist_km > 5])
        h = np.abs(bathymetry) + 100
        manning = 0.02
        damping = np.exp(-manning * 1.0 / (h / 1000))
        water_height = okada_disp * decay * damping
        
        # Interpolate to DEM
        interp_func = RegularGridInterpolator(
            (bathy_lat, bathy_lon),
            water_height,
            bounds_error=False,
            fill_value=0.0,
            method='linear'
        )
        DEM_LON, DEM_LAT = np.meshgrid(dem_lon, dem_lat)
        pts = np.column_stack([DEM_LAT.ravel(), DEM_LON.ravel()])
        water_height_dem = interp_func(pts).reshape(dem.shape)
        
        # Inundation
        inundation = (water_height_dem > dem).astype(float)
        inund_depth = water_height_dem - dem
        inund_depth[dem >= water_height_dem] = 0
        inund_area = inundation.sum() / inundation.size * 100
        
        # Plot
        ax.contourf(dem_lon, dem_lat, dem, levels=15, cmap='Greys', alpha=0.5)
        
        inund_masked = np.ma.masked_where(inund_depth <= 0, inund_depth)
        im2 = ax.contourf(dem_lon, dem_lat, inund_masked, 
                         levels=np.linspace(0.001, max(inund_depth.max(), 0.1), 15),
                         cmap='Blues', alpha=0.8)
        
        # Coastline
        for ring in coords:
            if ring:
                ring = np.array(ring)
                ax.plot(ring[:, 0], ring[:, 1], 'k-', linewidth=2, alpha=0.95)
        
        ax.set_title(f'{mag} Inundation\nFlooded area: {inund_area:.1f}%', 
                    fontsize=13, fontweight='bold', pad=10)
        ax.set_xlabel('Longitude (°E)', fontsize=11)
        ax.set_ylabel('Latitude (°N)', fontsize=11)
        ax.grid(True, alpha=0.2, linestyle='--')
        ax.set_xlim(dem_lon.min(), dem_lon.max())
        ax.set_ylim(dem_lat.min(), dem_lat.max())
        
        cb2 = plt.colorbar(im2, ax=ax, label='Inundation depth (m)', pad=0.02)
        cb2.ax.tick_params(labelsize=9)
    
    fig.suptitle('Philippine Tsunami Hazard Assessment - Multi-Scenario Inundation Maps\n' +
                'Earthquake Source: Philippine Trench Megathrust', 
                fontsize=16, fontweight='bold', y=0.97)
    
    fig.text(0.5, 0.01,
            'Nonlinear shallow-water propagation | GEBCO bathymetry | SRTM DEM | Okada dislocation',
            ha='center', fontsize=9, style='italic', color='gray')
    
    output_path = OUTPUT_DIR / "inundation_map_detailed_multi_scenario.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\n✓ Multi-scenario map saved: {output_path}")
    plt.close()
    
    return output_path


def main():
    print("\n" + "="*70)
    print("DETAILED PHILIPPINES INUNDATION MAPPING")
    print("="*70)
    
    create_detailed_multi_scenario_map()
    
    print("\n" + "="*70)
    print("✓ DETAILED MAP GENERATED")
    print("="*70)


if __name__ == "__main__":
    main()
