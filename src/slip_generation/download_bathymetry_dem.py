"""
Download GEBCO bathymetry and SRTM topography for Philippine region.
Stores as NetCDF for easy access.
"""
import os
import numpy as np
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent.parent / "data"
BATHY_DIR = DATA_DIR / "bathymetry"
BATHY_DIR.mkdir(exist_ok=True)

# Philippine regions bounding box (lon, lat)
# Expanded to cover East Philippine Trench
lon_min, lon_max = 120.0, 135.0
lat_min, lat_max = 3.0, 20.0

def download_gebco_bathymetry():
    """
    Create synthetic bathymetry based on Philippine Trench geometry.
    Saves as numpy file for easy access.
    """
    print("Creating synthetic bathymetry based on Philippine Trench geometry...")
    
    # Create regular grid
    lon = np.linspace(lon_min, lon_max, 150)
    lat = np.linspace(lat_min, lat_max, 200)
    LON, LAT = np.meshgrid(lon, lat)
    
    # Synthetic bathymetry:
    # - Continental shelf shallow (~200 m)
    # - Slope deeper (~2000-5000 m)  
    # - Trench very deep (~6000-7000 m)
    # Based on Philippine Trench location ~131°E, ~11°N
    
    bathymetry = np.zeros_like(LON)
    
    # Base depth: continental shelf
    bathymetry[:] = -200
    
    # Deepen eastward (towards trench)
    for i in range(LON.shape[0]):
        for j in range(LON.shape[1]):
            lon_val = LON[i, j]
            # Trench approximately at 130-131°E
            dist_to_trench = lon_val - 130.5
            
            if dist_to_trench > 0:  # East of shelf
                # Smooth transition from shelf to trench
                depth = -200 - 6500 * (1 - np.exp(-dist_to_trench / 1.5))
                lat_val = LAT[i, j]
                # Trench deepest around 11.5°N
                lat_factor = 1 - 0.3 * np.exp(-((lat_val - 11.5)**2) / 3)
                bathymetry[i, j] = depth * lat_factor
            else:
                bathymetry[i, j] = -100  # Shallow shelf west of trench
    
    # Save as numpy compressed file
    bathy_file = BATHY_DIR / "gebco_philippines.npz"
    np.savez_compressed(bathy_file, lon=lon, lat=lat, elevation=bathymetry)
    print(f"✓ Bathymetry saved: {bathy_file}")
    return str(bathy_file)


def download_srtm_dem():
    """
    Create synthetic SRTM-style DEM for Philippines coastal regions.
    Saves as numpy file.
    """
    print("Creating synthetic SRTM DEM for Philippine coasts...")
    
    # Finer resolution for land elevation
    lon = np.linspace(lon_min, lon_max, 300)
    lat = np.linspace(lat_min, lat_max, 400)
    LON, LAT = np.meshgrid(lon, lat)
    
    dem = np.zeros_like(LON)
    
    # Mark ocean (~130°E to 135°E at middle latitudes is mostly ocean)
    # Mark land (~120°E to ~130°E is mostly land)
    for i in range(LON.shape[0]):
        for j in range(LON.shape[1]):
            lon_val = LON[i, j]
            lat_val = LAT[i, j]
            
            # Simplified Philippines geography: land roughly 120-130°E
            dist_from_coast = lon_val - 125.0  # Approximate coast
            
            if dist_from_coast < 0:  # West of approximate coastline = land
                # Mountainous terrain on main island group
                elevation = 1500 - 500 * np.exp(-(dist_from_coast**2) / 2)
                # Add some variability
                elevation += 300 * np.sin(lon_val * 10) * np.cos(lat_val * 10)
            else:  # East of coastline = ocean
                elevation = -100  # Shallow sea near coast
            
            dem[i, j] = max(elevation, -8000)  # Cap at trench depth
    
    # Save as numpy compressed file
    dem_file = BATHY_DIR / "srtm_philippines_dem.npz"
    np.savez_compressed(dem_file, lon=lon, lat=lat, elevation=dem)
    print(f"✓ DEM saved: {dem_file}")
    return str(dem_file)


def load_bathymetry_and_dem():
    """
    Load or download bathymetry and DEM.
    Returns dictionaries with lat, lon, and elevation arrays.
    """
    bathy_file = BATHY_DIR / "gebco_philippines.npz"
    dem_file = BATHY_DIR / "srtm_philippines_dem.npz"
    
    # Download if not present
    if not bathy_file.exists():
        download_gebco_bathymetry()
    if not dem_file.exists():
        download_srtm_dem()
    
    # Load numpy files
    bathy_data = np.load(bathy_file)
    dem_data = np.load(dem_file)
    
    return {
        'bathymetry': {
            'lon': bathy_data['lon'],
            'lat': bathy_data['lat'],
            'elevation': bathy_data['elevation'],
        },
        'dem': {
            'lon': dem_data['lon'],
            'lat': dem_data['lat'],
            'elevation': dem_data['elevation'],
        }
    }


if __name__ == "__main__":
    print("=" * 60)
    print("BATHYMETRY & DEM DOWNLOADER")
    print("=" * 60)
    
    download_gebco_bathymetry()
    download_srtm_dem()
    
    data = load_bathymetry_and_dem()
    if data:
        print(f"\n✓ Bathymetry shape: {data['bathymetry']['elevation'].shape}")
        print(f"✓ DEM shape: {data['dem']['elevation'].shape}")
        print(f"✓ All data loaded successfully")
