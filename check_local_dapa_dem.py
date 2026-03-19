#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

try:
    import rasterio
    from rasterio.transform import xy as rio_xy
except Exception as exc:
    raise SystemExit("rasterio is required for this checker") from exc

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:
    plt = None

DAPA_MARKER = (126.055, 9.745)

@dataclass
class DemReport:
    path: str
    shape_y: int
    shape_x: int
    lon_min: float
    lon_max: float
    lat_min: float
    lat_max: float
    dx_m: float
    dy_m: float
    min_z: float
    max_z: float
    mean_z: float
    nan_count: int
    wet_cells: int
    land_cells: int
    edge_wet_bottom: int
    edge_wet_top: int
    edge_wet_left: int
    edge_wet_right: int
    target_i: int
    target_j: int
    target_lon: float
    target_lat: float
    target_z: float
    target_is_wet: bool
    target_connected_to_boundary: bool
    boundary_connected_wet_cells: int
    wet_cells_total: int
    likely_problem_flags: List[str]


def nearest_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(values - target)))


def read_tif(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    with rasterio.open(path) as ds:
        data = ds.read(1).astype(float)
        rows = np.arange(ds.height)
        cols = np.arange(ds.width)
        xs = np.array([rio_xy(ds.transform, 0, c, offset="center")[0] for c in cols], dtype=float)
        ys = np.array([rio_xy(ds.transform, r, 0, offset="center")[1] for r in rows], dtype=float)
        nodata = ds.nodata
    if nodata is not None:
        data[data == nodata] = np.nan
    ys = ys[::-1]
    data = data[::-1, :]
    return xs, ys, data


def bfs_connectivity(wet: np.ndarray, target_i: int, target_j: int) -> Tuple[int, bool]:
    ny, nx = wet.shape
    visited = np.zeros_like(wet, dtype=np.uint8)
    queue: List[Tuple[int, int]] = []
    for j in range(nx):
        if wet[0, j] and not visited[0, j]:
            visited[0, j] = 1
            queue.append((0, j))
        if wet[ny - 1, j] and not visited[ny - 1, j]:
            visited[ny - 1, j] = 1
            queue.append((ny - 1, j))
    for i in range(ny):
        if wet[i, 0] and not visited[i, 0]:
            visited[i, 0] = 1
            queue.append((i, 0))
        if wet[i, nx - 1] and not visited[i, nx - 1]:
            visited[i, nx - 1] = 1
            queue.append((i, nx - 1))
    head = 0
    nbrs = ((1, 0), (-1, 0), (0, 1), (0, -1))
    while head < len(queue):
        i, j = queue[head]
        head += 1
        for di, dj in nbrs:
            ii, jj = i + di, j + dj
            if 0 <= ii < ny and 0 <= jj < nx and wet[ii, jj] and not visited[ii, jj]:
                visited[ii, jj] = 1
                queue.append((ii, jj))
    connected = bool(visited[target_i, target_j]) if 0 <= target_i < ny and 0 <= target_j < nx else False
    return int(np.count_nonzero(visited)), connected


def analyze_dem(path: Path, wet_threshold_m: float) -> DemReport:
    lon, lat, z = read_tif(path)
    H = np.maximum(-z, 0.0)
    wet = np.isfinite(H) & (H > wet_threshold_m)
    land = np.isfinite(z) & (z > 0.0)
    mean_lat_rad = np.deg2rad(float(np.nanmean(lat)))
    dx_m = abs(float(np.nanmedian(np.diff(lon)))) * 111_320.0 * np.cos(mean_lat_rad)
    dy_m = abs(float(np.nanmedian(np.diff(lat)))) * 110_540.0
    gi = nearest_index(lat, DAPA_MARKER[1])
    gj = nearest_index(lon, DAPA_MARKER[0])
    connected_wet, target_connected = bfs_connectivity(wet.astype(np.uint8), gi, gj)

    flags: List[str] = []
    if np.count_nonzero(wet[0, :]) + np.count_nonzero(wet[-1, :]) + np.count_nonzero(wet[:, 0]) + np.count_nonzero(wet[:, -1]) < 10:
        flags.append("very few wet cells on DEM boundaries; nested boundary forcing may not enter the domain cleanly")
    if not wet[gi, gj]:
        flags.append("Dapa target cell is not wet in the local DEM; gauge and local arrival checks may be misleading")
    if wet[gi, gj] and not target_connected:
        flags.append("Dapa target wet cell is not connected to any wet boundary cells; the local DEM may be clipped too tightly or disconnected")
    if np.nanmin(z) >= 0:
        flags.append("DEM has no negative bathymetry values; this is unlikely for a topo-bathy nesting grid")
    if np.nanmax(z) <= 0:
        flags.append("DEM has no positive land elevations; this is unlikely for a Dapa inundation grid")
    if np.isnan(z).mean() > 0.01:
        flags.append("DEM contains many NaN/nodata cells")

    return DemReport(
        path=str(path),
        shape_y=int(z.shape[0]),
        shape_x=int(z.shape[1]),
        lon_min=float(lon[0]),
        lon_max=float(lon[-1]),
        lat_min=float(lat[0]),
        lat_max=float(lat[-1]),
        dx_m=float(dx_m),
        dy_m=float(dy_m),
        min_z=float(np.nanmin(z)),
        max_z=float(np.nanmax(z)),
        mean_z=float(np.nanmean(z)),
        nan_count=int(np.count_nonzero(~np.isfinite(z))),
        wet_cells=int(np.count_nonzero(wet)),
        land_cells=int(np.count_nonzero(land)),
        edge_wet_bottom=int(np.count_nonzero(wet[0, :])),
        edge_wet_top=int(np.count_nonzero(wet[-1, :])),
        edge_wet_left=int(np.count_nonzero(wet[:, 0])),
        edge_wet_right=int(np.count_nonzero(wet[:, -1])),
        target_i=int(gi),
        target_j=int(gj),
        target_lon=float(lon[gj]),
        target_lat=float(lat[gi]),
        target_z=float(z[gi, gj]),
        target_is_wet=bool(wet[gi, gj]),
        target_connected_to_boundary=bool(target_connected),
        boundary_connected_wet_cells=int(connected_wet),
        wet_cells_total=int(np.count_nonzero(wet)),
        likely_problem_flags=flags,
    )


def save_plot(path: Path, lon: np.ndarray, lat: np.ndarray, z: np.ndarray, report: DemReport) -> None:
    if plt is None:
        return
    fig, ax = plt.subplots(figsize=(7, 6), dpi=180)
    extent = [float(lon[0]), float(lon[-1]), float(lat[0]), float(lat[-1])]
    im = ax.imshow(z, origin="lower", extent=extent, aspect="auto")
    try:
        ax.contour(lon, lat, z, levels=[0.0], colors="black", linewidths=0.5)
    except Exception:
        pass
    ax.plot(DAPA_MARKER[0], DAPA_MARKER[1], marker="*", color="red", ms=8)
    ax.set_title("Local Dapa DEM diagnostic")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    fig.colorbar(im, ax=ax, shrink=0.85, label="Elevation / bathymetry (m)")
    txt = (
        f"edge wet B/T/L/R = {report.edge_wet_bottom}/{report.edge_wet_top}/{report.edge_wet_left}/{report.edge_wet_right}\n"
        f"target wet={report.target_is_wet} | connected={report.target_connected_to_boundary}"
    )
    ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=8, va="bottom", ha="left",
            bbox=dict(boxstyle="round", fc="white", ec="0.5", alpha=0.8))
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description="Check whether a local Dapa DEM is suitable for nested inundation forcing")
    p.add_argument("--dem", type=Path, required=True)
    p.add_argument("--wet-threshold-m", type=float, default=0.10)
    p.add_argument("--out-dir", type=Path, default=Path("output/dem_check"))
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    report = analyze_dem(args.dem, args.wet_threshold_m)
    report_path = args.out_dir / f"{args.dem.stem}_check.json"
    report_path.write_text(json.dumps(asdict(report), indent=2))
    print(json.dumps(asdict(report), indent=2))
    lon, lat, z = read_tif(args.dem)
    save_plot(args.out_dir / f"{args.dem.stem}_check.png", lon, lat, z, report)


if __name__ == "__main__":

    main()
