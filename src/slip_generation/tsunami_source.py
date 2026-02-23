"""
Tsunami source characterization from earthquake slip distributions.

This module converts earthquake slip distributions into seafloor vertical 
displacement (bathymetric perturbation), which serves as the initial 
condition for tsunami propagation modeling.

Two approaches are provided:
1. Simplified analytical elastic dislocation (Okada model approximation)
2. Interface for coupling to external tsunami solvers (e.g., GeoClaw)

References:
    - Okada (1985): Surface deformation due to shear and tensile faults
    - Gromov et al. (2021): Towards tsunami early warning systems
    - Satake (2007): Tsunamis: Case studies and lessons learned
"""

import numpy as np
from abc import ABC, abstractmethod


class TsunamiSource(ABC):
    """Abstract base class for tsunami source representations."""

    @abstractmethod
    def get_seafloor_displacement(self, observation_points):
        """
        Compute seafloor vertical displacement at specified locations.

        Args:
            observation_points (ndarray): Array of shape (n_points, 2) or (n_points, 3).
                                         Coordinates where displacement is requested.

        Returns:
            ndarray: Vertical displacement [m] at each observation point.
        """
        pass

    @abstractmethod
    def get_source_metadata(self):
        """
        Return metadata characterizing this tsunami source.

        Returns:
            dict: Dictionary with source parameters (depth, area, moment, etc.)
        """
        pass


class SimpleDislocationSource(TsunamiSource):
    """
    Simplified dislocation-based tsunami source model.

    This model uses an elastic dislocation representation to convert
    fault slip into seafloor vertical displacement. For computational
    efficiency and clarity, we use an approximate Green's function
    appropriate for shallow subduction zones.

    The model is physically-motivated but simplified:
    - Assumes flat fault geometry (rectangular patches)
    - Uses depth-average approximation for Green's function
    - Neglects focal mechanism subtleties (good for vertical slip)

    For production use, couple to Okada's full 3D dislocation or
    boundary element methods.
    """

    def __init__(self, fault_geometry, slip_distribution, rigidity=4.0e10):
        """
        Initialize dislocation source.

        Args:
            fault_geometry (dict): Fault parameters with keys:
                'length_km' (float): Fault length [km]
                'width_km' (float): Fault width (down-dip) [km]
                'top_depth_km' (float): Depth to top edge [km]
                'strike_deg' (float): Strike angle [degrees]
                'dip_deg' (float): Dip angle [degrees]
                'rake_deg' (float): Rake angle [degrees]
                'n_subfaults_along' (int): Number of subfaults along-strike
                'n_subfaults_down' (int): Number of subfaults down-dip
                'origin_lon' (float, optional): Reference longitude [°]
                'origin_lat' (float, optional): Reference latitude [°]

            slip_distribution (ndarray): Slip on each subfault [m].
                                        Shape: (n_down, n_along) matching fault discretization.

            rigidity (float, optional): Shear modulus μ [Pa]. Default: 4.0e10 (40 GPa).
        """
        self.fault_geom = fault_geometry
        self.slip = slip_distribution
        self.rigidity = rigidity
        
        # Extract geometry
        L = fault_geometry['length_km']
        W = fault_geometry['width_km']
        D = fault_geometry['top_depth_km']
        strike = np.radians(fault_geometry.get('strike_deg', 0.0))
        dip = np.radians(fault_geometry.get('dip_deg', 40.0))
        rake = np.radians(fault_geometry.get('rake_deg', 90.0))
        
        n_along = fault_geometry['n_subfaults_along']
        n_down = fault_geometry['n_subfaults_down']
        
        self.length_km = L
        self.width_km = W
        self.depth_km = D
        self.strike_rad = strike
        self.dip_rad = dip
        self.rake_rad = rake
        
        # Subfault dimensions
        self.subfault_length_km = L / n_along
        self.subfault_width_km = W / n_down
        
        # Build subfault center positions (in km, relative to fault origin)
        self._build_subfault_positions()

    def _build_subfault_positions(self):
        """Compute center positions of all subfaults."""
        n_along = self.fault_geom['n_subfaults_along']
        n_down = self.fault_geom['n_subfaults_down']
        
        # Along-strike direction
        x_centers = (np.arange(n_along) + 0.5) * self.subfault_length_km
        
        # Down-dip direction: depth increases with down-dip distance
        downdip_distances = (np.arange(n_down) + 0.5) * self.subfault_width_km
        z_centers = self.depth_km + downdip_distances * np.sin(self.dip_rad)
        
        # Mesh grid
        X, Z = np.meshgrid(x_centers, z_centers)  # shape: (n_down, n_along)
        
        self.subfault_positions = np.column_stack([X.ravel(), Z.ravel()])

    def get_seafloor_displacement(self, observation_points):
        """
        Compute vertical seafloor displacement using simplified Green's function.

        For a subduction zone fault at depth D with slip s, the vertical
        displacement at surface point (x, y) is approximated as:

            ζ(x, y) ≈ (μ / π) * ∑_i s_i * A_i / r_i^2

        where r_i is the distance from subfault i to observation point.

        This is a depth-averaged approximation valid for observation
        points far from the fault and shallow depths.

        Args:
            observation_points (ndarray): Shape (n_obs, 2) or (n_obs, 3).
                If (n_obs, 2): interpreted as (x_km, y_km) relative to fault origin.
                If (n_obs, 3): interpreted as (lon, lat, depth_km) [not yet implemented].

        Returns:
            ndarray: Vertical displacement [m] at each observation point.
        """
        n_obs = len(observation_points)
        
        # Handle 2D or 3D observation points
        if observation_points.shape[1] == 2:
            obs_2d = observation_points  # (x, y) in km
        else:
            # For simplicity, project 3D points to 2D
            obs_2d = observation_points[:, :2]
        
        displacement = np.zeros(n_obs)
        
        # Loop over subfaults
        n_subfaults = len(self.subfault_positions)
        for i in range(n_subfaults):
            # Subfault center position (only x, z components matter for vertical displacement)
            x_sub, z_sub = self.subfault_positions[i, 0], self.subfault_positions[i, 1]
            
            # Slip and area
            slip_i = self.slip.ravel()[i]
            area_i = (self.subfault_length_km * 1e3) * (self.subfault_width_km * 1e3)  # m²
            
            # Horizontal distance from each obs point to subfault center
            dx = obs_2d[:, 0] * 1e3 - x_sub * 1e3  # Convert to meters
            dy = obs_2d[:, 1] * 1e3  # Assume y_sub = 0 (along-strike centered)
            
            r_squared = dx**2 + dy**2 + (z_sub * 1e3)**2  # Distance squared (m)
            r_squared = np.maximum(r_squared, 1.0)  # Avoid division by zero
            
            # Simplified Green's function: vertical displacement ∝ 1/r²
            # Coefficient: μ * sin(rake) / π (sin(rake) ≈ 1 for vertical slip)
            green_coeff = (self.rigidity * np.sin(self.rake_rad)) / np.pi
            
            displacement += (green_coeff * slip_i * area_i) / r_squared
        
        return displacement

    def get_source_metadata(self):
        """Return characterization of this tsunami source."""
        seismic_moment = self.rigidity * np.sum(self.slip) * (
            self.subfault_length_km * 1e3 * self.subfault_width_km * 1e3
        )
        
        return {
            'fault_length_km': self.length_km,
            'fault_width_km': self.width_km,
            'fault_top_depth_km': self.depth_km,
            'strike_deg': np.degrees(self.strike_rad),
            'dip_deg': np.degrees(self.dip_rad),
            'rake_deg': np.degrees(self.rake_rad),
            'total_slip_m': np.sum(self.slip),
            'mean_slip_m': np.mean(self.slip),
            'max_slip_m': np.max(self.slip),
            'seismic_moment_nm': seismic_moment,
            'rigidity_pa': self.rigidity,
        }


class TsunamiPropagationInterface:
    """
    Abstract interface for coupling to external tsunami propagation solvers.

    This provides a standardized API for passing seafloor displacement
    initial conditions from this PTHA framework to solvers like GeoClaw,
    COMCOT, or other depth-integrated shallow water equation solvers.
    """

    def __init__(self, source):
        """
        Initialize interface with a tsunami source.

        Args:
            source (TsunamiSource): Tsunami source object.
        """
        self.source = source

    def export_source_to_netcdf(self, output_path, bathymetry_grid=None):
        """
        Export tsunami source initial condition to NetCDF format.

        Suitable for ingestion by GeoClaw or other standard tsunami solvers.

        Args:
            output_path (str): Path to write NetCDF file.
            bathymetry_grid (dict, optional): Bathymetry grid information:
                {'lon': 1D array, 'lat': 1D array, 'elevation': 2D array}
                If provided, displacement is added to bathymetry.

        Returns:
            None

        Notes:
            Requires scipy.io.netcdf or xarray. Not implemented in prototype.
        """
        raise NotImplementedError("NetCDF export not yet implemented in prototype.")

    def get_initial_condition_on_grid(self, lon_range, lat_range, n_lon, n_lat):
        """
        Compute tsunami source initial condition on a structured grid.

        Args:
            lon_range (tuple): (lon_min, lon_max) in degrees
            lat_range (tuple): (lat_min, lat_max) in degrees
            n_lon (int): Number of grid points in longitude
            n_lat (int): Number of grid points in latitude

        Returns:
            dict: Dictionary with keys:
                'lon': 1D array of longitude values
                'lat': 1D array of latitude values
                'displacement': 2D array of vertical displacement [m]
        """
        # Create grid
        lons = np.linspace(lon_range[0], lon_range[1], n_lon)
        lats = np.linspace(lat_range[0], lat_range[1], n_lat)
        
        # Map from geographic to fault-relative coordinates
        # (This is simplified; in practice, use proper coordinate transformations)
        LON, LAT = np.meshgrid(lons, lats)
        obs_points = np.column_stack([LON.ravel(), LAT.ravel()])
        
        # Compute displacement
        disp = self.source.get_seafloor_displacement(obs_points)
        disp_grid = disp.reshape(LON.shape)
        
        return {
            'lon': lons,
            'lat': lats,
            'displacement': disp_grid,
        }

    def get_source_metadata(self):
        """Get metadata from underlying source."""
        return self.source.get_source_metadata()


if __name__ == "__main__":
    print("=" * 60)
    print("Tsunami Source Representation Demonstration")
    print("=" * 60)
    
    # Example: Create a synthetic fault geometry and slip distribution
    fault_geom = {
        'length_km': 300.0,
        'width_km': 100.0,
        'top_depth_km': 7.6,
        'strike_deg': 164.0,
        'dip_deg': 39.0,
        'rake_deg': 90.0,
        'n_subfaults_along': 5,
        'n_subfaults_down': 3,
    }
    
    # Synthetic slip: increase toward center
    n_down, n_along = 3, 5
    slip_dist = np.ones((n_down, n_along)) * 2.0  # 2 m baseline
    slip_dist[1, 2] = 4.5  # Peak slip at center
    
    # Create source
    source = SimpleDislocationSource(fault_geom, slip_dist, rigidity=4.0e10)
    
    print("Fault geometry:")
    for key, val in fault_geom.items():
        print(f"  {key}: {val}")
    print()
    
    # Observation points (in km, relative to fault origin)
    obs_pts = np.array([
        [0, 0],     # At fault origin
        [100, 0],   # 100 km along-strike
        [0, 100],   # 100 km away perpendicular
        [150, 150], # Diagonal
    ])
    
    displacement = source.get_seafloor_displacement(obs_pts)
    
    print("Vertical seafloor displacement at observation points:")
    for pt, disp in zip(obs_pts, displacement):
        print(f"  Position {pt} km: {disp:.4f} m")
    print()
    
    # Source metadata
    metadata = source.get_source_metadata()
    print("Source metadata:")
    for key, val in metadata.items():
        if isinstance(val, float):
            if 'moment' in key.lower():
                print(f"  {key}: {val:.3e}")
            else:
                print(f"  {key}: {val:.4f}")
        else:
            print(f"  {key}: {val}")
