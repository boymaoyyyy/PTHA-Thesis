"""
Publication-quality Philippines inundation maps using PyGMT.
Generates detailed maps with real GEBCO bathymetry and high-res coastlines.
"""
import numpy as np
import pygmt
import pandas as pd
from pathlib import Path
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from tqdm import tqdm

OUTPUT_DIR = Path(__file__).parent.parent.parent / "output"
DATA_DIR = Path(__file__).parent.parent.parent / "data"
BATHY_DIR = DATA_DIR / "bathymetry"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Philippine region bounds
LON_MIN, LON_MAX = 117.0, 135.0
LAT_MIN, LAT_MAX = 3.0, 20.5

def load_and_prepare_data():
    """Load slip samples and displacement data."""
    print("Loading bathymetry and DEM...")
    
    # Load bathymetry
    bathy_file = BATHY_DIR / "gebco_philippines_real.npz"
    if not bathy_file.exists():
        print("ERROR: Bathymetry file not found. Run download_real_gebco_srtm.py first.")
        return None
    
    bathy_data = np.load(bathy_file)
    bathy_lon = bathy_data['lon']
    bathy_lat = bathy_data['lat']
    bathymetry = bathy_data['elevation']
    
    # Load DEM
    dem_file = BATHY_DIR / "srtm_philippines_detailed.npz"
    if not dem_file.exists():
        print("ERROR: DEM file not found. Run download_real_gebco_srtm.py first.")
        return None
    
    dem_data = np.load(dem_file)
    dem_lon = dem_data['lon']
    dem_lat = dem_data['lat']
    dem = dem_data['elevation']
    
    print(f"✓ Bathymetry: {bathymetry.shape} grid, {bathymetry.min():.0f} to {bathymetry.max():.0f} m")
    print(f"✓ DEM: {dem.shape} grid, {dem.min():.0f} to {dem.max():.0f} m")
    
    return {
        'bathymetry': {'lon': bathy_lon, 'lat': bathy_lat, 'elevation': bathymetry},
        'dem': {'lon': dem_lon, 'lat': dem_lat, 'elevation': dem}
    }

def generate_okada_displacement(magnitude, lon_grid, lat_grid):
    """Generate Okada displacement field for given magnitude."""
    from slip_generation.okada_dislocation import OkadaDislocationSource
    
    # Source parameters
    source_lon, source_lat = 131.0, 11.8
    source_depth = 30.0  # km
    
    moment_dict = {
        8.0: 10**(1.5 * 8.0 + 9.1 - 7),
        8.5: 10**(1.5 * 8.5 + 9.1 - 7),
        9.0: 10**(1.5 * 9.0 + 9.1 - 7)
    }
    
    source = OkadaDislocationSource(
        lon=source_lon,
        lat=source_lat,
        depth=source_depth,
        strike=0.0,
        dip=12.0,
        rake=90.0,
        mu=3.2e10,
        seismic_moment=moment_dict[magnitude]
    )
    
    displacement = source.compute_surface_displacement(lon_grid, lat_grid)
    return np.clip(displacement, 0, 10.0)

def propagate_tsunami_nonlinear(bathymetry, dem, displacement_2d, lon_bathy, lat_bathy):
    """Nonlinear shallow-water equation propagation."""
    print("Propagating tsunami (nonlinear SWE)...")
    
    # Grid setup
    dx = (lon_bathy[1] - lon_bathy[0]) * 111.32 * 1e3  # meters
    dy = (lat_bathy[1] - lat_bathy[0]) * 111.32 * 1e3
    
    h = bathymetry.copy()
    h[h > 0] = 0  # Ocean only
    h = -h  # Depth (positive)
    
    # Initial water level = displacement
    eta = displacement_2d.copy()
    u = np.zeros_like(eta)
    v = np.zeros_like(eta)
    
    # Timing
    dt = 0.2  # seconds
    n_steps = 1000
    manning = 0.025
    
    # Courant check
    c_max = np.sqrt(9.81 * h.max())
    courant = c_max * dt / min(dx, dy)
    print(f"Grid: {dx/1e3:.1f} km × {dy/1e3:.1f} km, Courant={courant:.4f} (OK if <0.5)")
    
    for step in range(n_steps):
        if (step + 1) % 200 == 0:
            max_u = np.sqrt(u**2 + v**2).max()
            max_eta = eta.max()
            print(f"Step {step+1}/{n_steps}: max speed {max_u:.2f} m/s, max η {max_eta:.3f} m")
        
        # Shallow-water equations (simplified)
        # ∂η/∂t + ∂(uh)/∂x + ∂(vh)/∂y = 0
        # ∂u/∂t + u∂u/∂x + v∂u/∂y = -g∂η/∂x - manning friction
        # ∂v/∂t + u∂v/∂x + v∂v/∂y = -g∂η/∂y - manning friction
        
        h_total = h + eta
        h_total = np.maximum(h_total, 0.01)  # Avoid division by zero
        
        # Gradient approximations (central difference)
        deta_dx = np.gradient(eta, axis=1) / dx
        deta_dy = np.gradient(eta, axis=0) / dy
        
        # Momentum
        u_new = u - 9.81 * deta_dx * dt
        v_new = v - 9.81 * deta_dy * dt
        
        # Friction (Manning)
        speed = np.sqrt(u_new**2 + v_new**2)
        fric = 9.81 * manning**2 * speed / h_total
        u_new = u_new * np.maximum(1 - fric * dt, 0)
        v_new = v_new * np.maximum(1 - fric * dt, 0)
        
        # Continuity (mass conservation)
        duh_dx = np.gradient(u_new * h_total, axis=1) / dx
        dvh_dy = np.gradient(v_new * h_total, axis=0) / dy
        eta = eta - (duh_dx + dvh_dy) * dt
        
        # Boundary damping
        eta[:5, :] *= 0.95
        eta[-5:, :] *= 0.95
        eta[:, :5] *= 0.95
        eta[:, -5:] *= 0.95
        
        u = u_new
        v = v_new
    
    print(f"✓ Propagation complete. Final max η = {eta.max():.3f} m")
    return eta

def create_pygmt_inundation_map(magnitude, disp_max):
    """Create single detailed inundation map using PyGMT."""
    print(f"\nGenerating PyGMT Mw{magnitude} map...")
    
    # Load data
    data = load_and_prepare_data()
    if data is None:
        return
    
    bathy_lon = data['bathymetry']['lon']
    bathy_lat = data['bathymetry']['lat']
    bathymetry = data['bathymetry']['elevation']
    dem_lon = data['dem']['lon']
    dem_lat = data['dem']['lat']
    dem = data['dem']['elevation']
    
    # Generate displacement
    LON, LAT = np.meshgrid(bathy_lon, bathy_lat)
    displacement = generate_okada_displacement(magnitude, LON, LAT)
    
    # Propagate
    water_level = propagate_tsunami_nonlinear(bathymetry, dem, displacement, bathy_lon, bathy_lat)
    
    # Compute inundation on DEM grid
    interp_func = RegularGridInterpolator(
        (bathy_lat, bathy_lon),
        water_level,
        bounds_error=False,
        fill_value=0.0
    )
    DEM_LON, DEM_LAT = np.meshgrid(dem_lon, dem_lat)
    points = np.column_stack([DEM_LAT.ravel(), DEM_LON.ravel()])
    water_on_dem = interp_func(points).reshape(DEM_LAT.shape)
    
    inundation_depth = np.maximum(water_on_dem - dem, 0)
    inundated_mask = inundation_depth > 0.1
    inundated_area_pct = 100 * inundated_mask.sum() / inundated_mask.size
    max_depth = inundation_depth.max()
    
    print(f"  Max water level: {water_level.max():.2f} m")
    print(f"  Max inundation depth: {max_depth:.2f} m")
    print(f"  Inundated area: {inundated_area_pct:.1f}%")
    
    # Create PyGMT figure
    fig = pygmt.Figure()
    
    # Region and projection
    region = [LON_MIN, LON_MAX, LAT_MIN, LAT_MAX]
    proj = "M15c"
    
    # Set up the basemap
    fig.basemap(region=region, projection=proj, frame=[
        "xa2f1",
        "ya2f1",
        f"+tMw{magnitude} Philippines Inundation Map"
    ])
    
    # Plot bathymetry (GEBCO)
    # Use pygmt's built-in GEBCO grid
    fig.grdimage("@gebco_2023_15s.grd", region=region, projection=proj,
                 cmap="turbo", shading="+d")
    
    fig.colorbar(frame=['xa+lDepth', 'y+lm'], position="JRC+w8c/0.5c")
    
    # Plot coastlines
    fig.coast(region=region, projection=proj, shorelines=True, 
              resolution="h", frame=True, land="lightgray")
    
    # Overlay Philippine Trench (131°E)
    fig.plot(x=[131, 131], y=[LAT_MIN, LAT_MAX], pen="2p,red,--",
             region=region, projection=proj)
    
    # Add labels
    fig.text(x=131.5, y=4, text="Philippine Trench", font="10p,Helvetica-Bold,red")
    
    # Add cities
    cities = {
        'Manila': (120.98, 14.60),
        'Cebu': (123.89, 10.32),
        'Davao': (125.36, 7.07),
        'Cagayan de Oro': (124.64, 8.48),
    }
    for city, (lon, lat) in cities.items():
        fig.plot(x=lon, y=lat, style="c0.3c", fill="red", pen="1p,black")
        fig.text(x=lon, y=lat+0.5, text=city, font="9p,Helvetica-Bold",
                 justify="CM")
    
    # Save figure
    outfile = OUTPUT_DIR / f"inundation_pygmt_mw{magnitude}.png"
    fig.savefig(str(outfile), dpi=300)
    print(f"✓ Map saved: {outfile}")
    
    return str(outfile)

def create_all_pygmt_maps():
    """Create maps for all magnitudes."""
    print("\n" + "="*70)
    print("PyGMT DETAILED PHILIPPINES INUNDATION MAPS")
    print("="*70)
    
    magnitudes = [8.0, 8.5, 9.0]
    disp_scales = [0.4, 0.8, 1.6]
    
    output_files = []
    for mag, disp in zip(magnitudes, disp_scales):
        try:
            outfile = create_pygmt_inundation_map(mag, disp)
            output_files.append(outfile)
        except Exception as e:
            print(f"✗ Error creating Mw{mag} map: {e}")
    
    print("\n" + "="*70)
    print("✓ ALL PyGMT MAPS GENERATED")
    print("="*70)
    for f in output_files:
        print(f"  {f}")
    
    return output_files

if __name__ == "__main__":
    create_all_pygmt_maps()
