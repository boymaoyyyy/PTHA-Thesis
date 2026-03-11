"""
Generate spatial maps of Okada mean displacement per magnitude.

Produces a 1x3 figure similar to published deformation figures.  If
network access is available the Philippine coastline is drawn using
Natural Earth geojson; otherwise a simple rectangular boundary is used.

Usage:
    python phase3_okada_spatial_maps.py

Outputs:
    output/phase3_okada_spatial_maps.png
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))

from phase3_okada import create_observation_grid


def fetch_philippines_coastline():
    """Download a minimal geojson and extract Philippines polygon(s)."""
    try:
        import urllib.request, json
        url = (
            "https://raw.githubusercontent.com/datasets/geo-boundaries-world-110m/"
            "master/countries.geojson"
        )
        print("Fetching coastline data from Natural Earth...")
        with urllib.request.urlopen(url, timeout=10) as f:
            gj = json.load(f)
        polys = []
        for feat in gj['features']:
            props = feat.get('properties', {})
            name = props.get('admin') or props.get('name') or ''
            if name.lower() == 'philippines':
                geom = feat['geometry']
                coords = geom['coordinates']
                if geom['type'] == 'Polygon':
                    polys.append(coords)
                elif geom['type'] == 'MultiPolygon':
                    polys.extend(coords)
        if not polys:
            print("WARNING: Philippines feature not found in geojson")
        return polys
    except Exception as e:
        print("Could not download coastline data:", e)
        return []


def plot_maps(mean_displacements, lons, lats, output_file=None):
    n = len(mean_displacements)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 6), constrained_layout=True)
    if n == 1:
        axes = [axes]

    # attempt to get coastline
    coast_polys = fetch_philippines_coastline()
    
    for ax, (mag_str, disp) in zip(axes, mean_displacements.items()):
        disp_grid = disp.reshape(len(lats), len(lons))
        mesh = ax.pcolormesh(lons, lats, disp_grid, cmap='viridis', shading='auto')
        ax.set_title(mag_str + ' mean displacement', fontweight='bold')
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        ax.set_xlim(lons.min(), lons.max())
        ax.set_ylim(lats.min(), lats.max())
        fig.colorbar(mesh, ax=ax, label='Displacement (m)')

        # overlay coastline if available
        if coast_polys:
            for poly in coast_polys:
                # poly may be list of linear rings
                for ring in poly:
                    ring = np.array(ring)
                    ax.plot(ring[:,0], ring[:,1], color='black', linewidth=0.8)
        else:
            # fallback: draw domain rectangle
            ax.plot([lons.min(), lons.max(), lons.max(), lons.min(), lons.min()],
                    [lats.min(), lats.min(), lats.max(), lats.max(), lats.min()],
                    color='black', linewidth=1)
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Spatial maps saved: {output_file}")
    else:
        plt.show()
    plt.close()


def compute_mean_displacement():
    # reuse observation grid from phase3_okada
    _, lons, lats = create_observation_grid()
    disp_root = Path('output') / 'tsunami_sources_okada'
    mean_disp = {}
    for mag in [8.0, 8.5, 9.0]:
        mag_str = f'Mw{mag:.1f}'
        mag_dir = disp_root / mag_str
        if not mag_dir.exists():
            print(f"Warning: directory {mag_dir} missing, skipping")
            continue
        files = sorted(mag_dir.glob('*_displacement_okada.npy'))
        if not files:
            print(f"No displacement files in {mag_dir}")
            continue
        arr = np.vstack([np.load(f).ravel() for f in files])
        mean_disp[mag_str] = arr.mean(axis=0)
    return mean_disp, lons, np.array(lats)


if __name__ == '__main__':
    mean_disp, lons, lats = compute_mean_displacement()
    if mean_disp:
        out = Path('output') / 'phase3_okada_spatial_maps.png'
        plot_maps(mean_disp, lons, lats, output_file=out)
    else:
        print("No displacement data available to plot.")
