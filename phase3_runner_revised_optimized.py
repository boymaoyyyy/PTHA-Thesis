#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Optional, Any

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from scipy import interpolate, spatial, special, stats

try:
    from scipy.stats import qmc
except Exception:
    qmc = None

try:
    import xarray as xr
except Exception:
    xr = None

try:
    import rasterio
    from rasterio.transform import from_origin
    from rasterio.transform import xy as rio_xy
except Exception:
    rasterio = None
    from_origin = None
    rio_xy = None

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
except Exception:
    plt = None
    BoundaryNorm = None
    ListedColormap = None

try:
    import pygmt
except Exception:
    pygmt = None

try:
    from tqdm import tqdm
except Exception:
    tqdm = None

try:
    from numba import njit, set_num_threads
    NUMBA_AVAILABLE = True
except Exception:
    NUMBA_AVAILABLE = False
    def njit(*args, **kwargs):
        def deco(fn):
            return fn
        return deco
    def set_num_threads(n):
        return None


PAPER_INUNDATION_BOUNDS = np.array([0.01, 0.25, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50], dtype=float)
PAPER_INUNDATION_COLORS = [
    "#66b3ff",
    "#4aa3ff",
    "#1f7af2",
    "#0d63db",
    "#0a53c1",
    "#083f9c",
    "#082d7a",
    "#051d56",
]
DAPA_MARKER = (126.055, 9.745)


@dataclass
class DomainConfig:
    regional_lon_min: float = 123.5
    regional_lon_max: float = 127.8
    regional_lat_min: float = 6.0
    regional_lat_max: float = 12.0
    dapa_lon_min: float = 125.98
    dapa_lon_max: float = 126.12
    dapa_lat_min: float = 9.69
    dapa_lat_max: float = 9.82
    corridor_lon_min: float = 125.45
    corridor_lon_max: float = 126.95
    corridor_lat_min: float = 7.00
    corridor_lat_max: float = 10.30


@dataclass
class SourceConfig:
    mw: float = 8.5
    shear_modulus_pa: float = 4.0e10
    rupture_length_km: float = 300.0
    rupture_width_km: float = 100.0
    top_depth_km: float = 7.6
    strike_deg: float = 164.0
    dip_deg: float = 39.0
    rake_deg: float = 90.0
    mean_slip_m: float = 4.7
    uplift_fraction_of_slip: float = 0.15
    rise_time_sec_per_km: float = 1.2


@dataclass
class SlipConfig:
    n_along: int = 15
    n_down: int = 5
    corr_length_km: float = 20.0
    hurst_nu: float = 0.3
    rqmc_seed: int = 42


@dataclass
class SolverConfig:
    gravity: float = 9.81
    cfl: float = 0.45
    dt_max: float = 0.15
    t_end_sec: float = 10800.0
    snapshot_stride: int = 60
    diag_stride: int = 250
    min_depth_m: float = 0.05
    eta_guard_m: float = 200.0
    damping: float = 2.0e-4
    viscosity_factor: float = 0.02


@dataclass
class ProbabilityConfig:
    beta: float = 0.9
    thresholds_m: Tuple[float, ...] = (0.5, 1.0, 2.0)


@dataclass
class ExecutionConfig:
    checkpoint_every: int = 25
    downsample_factor: int = 1
    save_per_sample_maps: bool = False
    save_first_gauge_series: bool = True
    paper_style_maps: bool = True
    local_dapa_render: bool = True
    local_dapa_nested: bool = True
    nested_save_stride: int = 10
    write_geotiff: bool = True
    write_npz: bool = True
    write_png: bool = True
    write_catalog: bool = True
    emit_phase4_preview: bool = False
    kernel_backend: str = "auto"
    numba_threads: int = 0


@dataclass
class GaugePoint:
    name: str
    lon: float
    lat: float


DEFAULT_GAUGES = (
    GaugePoint("Davao", 125.6629, 7.1537),
    GaugePoint("Legazpi", 123.7581, 13.1462),
    GaugePoint("Dapa_virtual", DAPA_MARKER[0], DAPA_MARKER[1]),
)


@dataclass
class ProjectConfig:
    project_root: Path
    domain: DomainConfig = field(default_factory=DomainConfig)
    source: SourceConfig = field(default_factory=SourceConfig)
    slip: SlipConfig = field(default_factory=SlipConfig)
    solver: SolverConfig = field(default_factory=SolverConfig)
    probability: ProbabilityConfig = field(default_factory=ProbabilityConfig)
    execution: ExecutionConfig = field(default_factory=ExecutionConfig)
    gauges: Tuple[GaugePoint, ...] = DEFAULT_GAUGES
    dapa_only_mode: bool = False
    local_dem_path_override: Optional[Path] = None
    slip_dir_override: Optional[Path] = None
    results_dirname: str = "phase3_revised_outputs"

    @property
    def bathymetry_dir(self) -> Path:
        return self.project_root / "data" / "bathymetry"

    @property
    def topography_dir(self) -> Path:
        return self.project_root / "data" / "topography"

    @property
    def slip_dir(self) -> Path:
        if self.slip_dir_override is not None:
            p = self.slip_dir_override
            return p if p.is_absolute() else (self.project_root / p)
        return self.project_root / "output" / "slip_samples"

    @property
    def results_dir(self) -> Path:
        return self.project_root / "output" / self.results_dirname

    @property
    def merged_dapa_grid_path(self) -> Path:
        if self.local_dem_path_override is not None:
            p = self.local_dem_path_override
            return p if p.is_absolute() else (self.project_root / p)
        return self.project_root / "output" / "dapa_merged_topobathy_fine.tif"


@dataclass
class Grid2D:
    lon: np.ndarray
    lat: np.ndarray
    z: np.ndarray

    @property
    def shape(self) -> Tuple[int, int]:
        return self.z.shape

    @property
    def dx_m(self) -> float:
        dlon = float(np.nanmedian(np.diff(self.lon if self.lon.ndim == 1 else self.lon[0, :])))
        mean_lat = math.radians(float(np.nanmean(self.lat if self.lat.ndim == 1 else self.lat[:, 0])))
        return abs(dlon) * 111_320.0 * math.cos(mean_lat)

    @property
    def dy_m(self) -> float:
        dlat = float(np.nanmedian(np.diff(self.lat if self.lat.ndim == 1 else self.lat[:, 0])))
        return abs(dlat) * 110_540.0


@dataclass
class PreparedDomain:
    grid: Grid2D
    H: np.ndarray
    z: np.ndarray
    dx: float
    dy: float
    gauge_idx: Dict[str, Tuple[int, int]]
    gauge_interp: Dict[str, Tuple[int, int, float, float]]
    dapa_ys: slice
    dapa_xs: slice
    topo_grid_local: Optional[Grid2D] = None
    nested_local_domain: Optional["PreparedDomain"] = None


@dataclass
class LocalNestedOutput:
    max_eta: np.ndarray
    max_depth: np.ndarray
    max_speed: np.ndarray
    arrival_time: np.ndarray
    wet_mask: np.ndarray
    dapa_p95_inun: float
    meta: Dict[str, Any]


@dataclass
class SampleOutput:
    slip_id: str
    forcing_mode: str
    gauge_time_sec: np.ndarray
    gauge_eta: Dict[str, np.ndarray]
    max_eta: np.ndarray
    max_depth: np.ndarray
    max_speed: np.ndarray
    arrival_time: np.ndarray
    wet_mask: np.ndarray
    gauge_max_abs: Dict[str, float]
    dapa_p95_inun: float
    meta: Dict[str, Any]
    nested_time_sec: Optional[np.ndarray] = None
    nested_eta_box: Optional[np.ndarray] = None
    local_nested: Optional[LocalNestedOutput] = None


@dataclass
class ValidationObservation:
    gauge_name: str
    csv_path: Path
    time_col: str = "time_sec"
    eta_col: str = "eta_m"


@dataclass
class HistoricalEventConfig:
    event_id: str
    mw: float
    output_subdir: str
    source_overrides: Dict[str, float]
    gauges: List[GaugePoint]
    observations: List[ValidationObservation]


# ---------- filesystem / json ----------

def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def save_json(path: Path, payload: Dict[str, Any]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, default=str)


# ---------- IO helpers ----------

def list_files(folder: Path, suffixes: Sequence[str]) -> List[Path]:
    if not folder.exists():
        return []
    out: List[Path] = []
    for s in suffixes:
        out.extend(sorted(folder.rglob(f"*{s}")))
    return out


def nearest_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(values - target)))


def find_bbox_indices(
    lon_1d: np.ndarray,
    lat_1d: np.ndarray,
    lon_min: float,
    lon_max: float,
    lat_min: float,
    lat_max: float,
) -> Tuple[slice, slice]:
    ix0 = nearest_index(lon_1d, lon_min)
    ix1 = nearest_index(lon_1d, lon_max)
    iy0 = nearest_index(lat_1d, lat_min)
    iy1 = nearest_index(lat_1d, lat_max)
    if ix0 > ix1:
        ix0, ix1 = ix1, ix0
    if iy0 > iy1:
        iy0, iy1 = iy1, iy0
    return slice(iy0, iy1 + 1), slice(ix0, ix1 + 1)


def expand_bbox(lon_min: float, lon_max: float, lat_min: float, lat_max: float, factor: float = 1.35) -> Tuple[float, float, float, float]:
    lon_c = 0.5 * (lon_min + lon_max)
    lat_c = 0.5 * (lat_min + lat_max)
    lon_h = 0.5 * (lon_max - lon_min) * factor
    lat_h = 0.5 * (lat_max - lat_min) * factor
    return lon_c - lon_h, lon_c + lon_h, lat_c - lat_h, lat_c + lat_h


def safe_extent(lon: np.ndarray, lat: np.ndarray) -> List[float]:
    x0, x1 = float(lon[0]), float(lon[-1])
    y0, y1 = float(lat[0]), float(lat[-1])
    if np.isclose(x0, x1):
        dx = float(np.median(np.diff(lon))) if len(lon) > 1 else 1.0e-4
        x0 -= 0.5 * abs(dx)
        x1 += 0.5 * abs(dx)
    if np.isclose(y0, y1):
        dy = float(np.median(np.diff(lat))) if len(lat) > 1 else 1.0e-4
        y0 -= 0.5 * abs(dy)
        y1 += 0.5 * abs(dy)
    return [x0, x1, y0, y1]


def describe_subgrid(grid: Grid2D, H: np.ndarray, ys: slice, xs: slice) -> Dict[str, float]:
    subH = H[ys, xs]
    subz = grid.z[ys, xs]
    return {
        "shape_y": int(subH.shape[0]),
        "shape_x": int(subH.shape[1]),
        "wet_cells": int(np.count_nonzero(subH > 0.0)),
        "wet_cells_gt_min": int(np.count_nonzero(subH > 0.05)),
        "max_static_depth": float(np.max(subH)) if subH.size else 0.0,
        "min_z": float(np.min(subz)) if subz.size else 0.0,
        "max_z": float(np.max(subz)) if subz.size else 0.0,
    }


def bilinear_resize(src: np.ndarray, out_shape: Tuple[int, int]) -> np.ndarray:
    sy, sx = src.shape
    oy, ox = out_shape
    y_old = np.linspace(0.0, 1.0, sy)
    x_old = np.linspace(0.0, 1.0, sx)
    y_new = np.linspace(0.0, 1.0, oy)
    x_new = np.linspace(0.0, 1.0, ox)
    f = interpolate.RegularGridInterpolator((y_old, x_old), src, bounds_error=False, fill_value=None)
    yy, xx = np.meshgrid(y_new, x_new, indexing="ij")
    return f(np.column_stack([yy.ravel(), xx.ravel()])).reshape(out_shape)


def choose_best_raster(files: List[Path]) -> Path:
    rank = {".tif": 0, ".tiff": 0, ".nc": 1, ".npy": 2}
    return sorted(files, key=lambda p: (rank.get(p.suffix.lower(), 99), -p.stat().st_size))[0]


def read_grid(path: Path) -> Grid2D:
    suf = path.suffix.lower()
    if suf in {".tif", ".tiff"}:
        if rasterio is None or rio_xy is None:
            raise ImportError("rasterio is required to read GeoTIFF files.")
        with rasterio.open(path) as ds:
            data = ds.read(1).astype(float)
            rows = np.arange(ds.height)
            cols = np.arange(ds.width)
            xs = np.array([rio_xy(ds.transform, 0, c, offset="center")[0] for c in cols], dtype=float)
            ys = np.array([rio_xy(ds.transform, r, 0, offset="center")[1] for r in rows], dtype=float)
        ys = ys[::-1]
        data = data[::-1, :]
        return Grid2D(lon=xs, lat=ys, z=data)
    if suf == ".nc":
        if xr is None:
            raise ImportError("xarray is required to read NetCDF files.")
        ds = xr.open_dataset(path)
        try:
            vname = next(name for name, da in ds.data_vars.items() if da.ndim >= 2)
            da = ds[vname].squeeze()
            lon_name = next(c for c in da.coords if c.lower() in {"lon", "longitude", "x"})
            lat_name = next(c for c in da.coords if c.lower() in {"lat", "latitude", "y"})
            lon = np.asarray(da[lon_name].values, dtype=float)
            lat = np.asarray(da[lat_name].values, dtype=float)
            z = np.asarray(da.values, dtype=float)
            if lat[0] > lat[-1]:
                lat = lat[::-1]
                z = z[::-1, :]
            return Grid2D(lon=lon, lat=lat, z=z)
        finally:
            ds.close()
    if suf == ".npy":
        arr = np.load(path, allow_pickle=True)
        ny, nx = arr.shape
        lon = np.linspace(123.0, 128.0, nx)
        lat = np.linspace(6.0, 12.0, ny)
        return Grid2D(lon=lon, lat=lat, z=np.asarray(arr, dtype=float))
    raise ValueError(f"Unsupported raster format: {path}")


def crop_grid(grid: Grid2D, lon_min: float, lon_max: float, lat_min: float, lat_max: float) -> Grid2D:
    lon_1d = grid.lon if grid.lon.ndim == 1 else grid.lon[0, :]
    lat_1d = grid.lat if grid.lat.ndim == 1 else grid.lat[:, 0]
    ys, xs = find_bbox_indices(lon_1d, lat_1d, lon_min, lon_max, lat_min, lat_max)
    return Grid2D(lon=lon_1d[xs], lat=lat_1d[ys], z=grid.z[ys, xs])


def downsample_grid(grid: Grid2D, factor: int) -> Grid2D:
    if factor <= 1:
        return grid
    return Grid2D(lon=grid.lon[::factor], lat=grid.lat[::factor], z=grid.z[::factor, ::factor])


def load_combined_domain(cfg: ProjectConfig) -> Tuple[Grid2D, Optional[Grid2D], Optional[Grid2D]]:
    bathy_files = list_files(cfg.bathymetry_dir, [".tif", ".tiff", ".nc", ".npy"])
    topo_files = list_files(cfg.topography_dir, [".tif", ".tiff", ".nc", ".npy"])
    if not bathy_files or not topo_files:
        raise FileNotFoundError("Bathymetry/topography files not found.")

    print("[BATHYMETRY] Loading GEBCO...", flush=True)
    bathy_path = choose_best_raster(bathy_files)
    print(f"[BATHYMETRY] Using {bathy_path.name}", flush=True)
    bathy = read_grid(bathy_path)

    print("[TOPOGRAPHY] Loading SRTM...", flush=True)
    topo_path = choose_best_raster(topo_files)
    print(f"[TOPOGRAPHY] Using {topo_path.name}", flush=True)
    topo = read_grid(topo_path)

    bathy = crop_grid(
        bathy,
        cfg.domain.regional_lon_min,
        cfg.domain.regional_lon_max,
        cfg.domain.regional_lat_min,
        cfg.domain.regional_lat_max,
    )
    topo_reg = crop_grid(
        topo,
        cfg.domain.regional_lon_min,
        cfg.domain.regional_lon_max,
        cfg.domain.regional_lat_min,
        cfg.domain.regional_lat_max,
    )

    topo_rs = bilinear_resize(topo_reg.z, bathy.z.shape)
    z = np.where(topo_rs > 0.0, topo_rs, -np.abs(bathy.z))
    grid = Grid2D(lon=bathy.lon.copy(), lat=bathy.lat.copy(), z=z)
    if cfg.execution.downsample_factor > 1:
        grid = downsample_grid(grid, cfg.execution.downsample_factor)

    topo_local = None
    local_merged = None
    need_local = cfg.execution.local_dapa_render or cfg.execution.local_dapa_nested

    if need_local:
        if cfg.execution.local_dapa_render:
            topo_local = crop_grid(
                topo,
                cfg.domain.dapa_lon_min,
                cfg.domain.dapa_lon_max,
                cfg.domain.dapa_lat_min,
                cfg.domain.dapa_lat_max,
            )

        merged_path = cfg.merged_dapa_grid_path
        if merged_path.exists():
            print(f"[LOCAL DEM] Using {merged_path.name}", flush=True)
            local_merged = read_grid(merged_path)
            if cfg.execution.local_dapa_render:
                topo_local = local_merged
        elif cfg.execution.local_dapa_nested:
            print(f"[LOCAL DEM] Not found at {merged_path}; local nested inundation will be skipped.", flush=True)

    return grid, topo_local, local_merged


def sample_gauges(grid: Grid2D, H: np.ndarray, gauges: Sequence[GaugePoint], min_depth: float = 0.50) -> Dict[str, Tuple[int, int]]:
    idx: Dict[str, Tuple[int, int]] = {}
    wet_all = np.argwhere(H > min_depth)
    for g in gauges:
        i0 = nearest_index(grid.lat, g.lat)
        j0 = nearest_index(grid.lon, g.lon)
        if 0 <= i0 < H.shape[0] and 0 <= j0 < H.shape[1] and H[i0, j0] > min_depth:
            idx[g.name] = (i0, j0)
            continue

        best = None
        best_dist = float("inf")
        for radius in range(1, 41):
            i_min, i_max = max(0, i0 - radius), min(H.shape[0], i0 + radius + 1)
            j_min, j_max = max(0, j0 - radius), min(H.shape[1], j0 + radius + 1)
            wet = np.argwhere(H[i_min:i_max, j_min:j_max] > min_depth)
            if wet.size == 0:
                continue
            for wi, wj in wet:
                ii, jj = i_min + int(wi), j_min + int(wj)
                d = math.hypot(float(grid.lat[ii] - g.lat), float(grid.lon[jj] - g.lon))
                if d < best_dist:
                    best = (ii, jj)
                    best_dist = d
            if best is not None:
                break

        if best is None and wet_all.size:
            for wi, wj in wet_all:
                ii, jj = int(wi), int(wj)
                d = math.hypot(float(grid.lat[ii] - g.lat), float(grid.lon[jj] - g.lon))
                if d < best_dist:
                    best = (ii, jj)
                    best_dist = d

        idx[g.name] = best if best is not None else (i0, j0)
    return idx




def build_gauge_interp(grid: Grid2D, gauges: Sequence[GaugePoint]) -> Dict[str, Tuple[int, int, float, float]]:
    lon = np.asarray(grid.lon, dtype=float)
    lat = np.asarray(grid.lat, dtype=float)
    if lon.size < 2 or lat.size < 2:
        return {g.name: (0, 0, 0.0, 0.0) for g in gauges}
    lon0 = float(lon[0])
    lat0 = float(lat[0])
    dlon = float(lon[1] - lon[0])
    dlat = float(lat[1] - lat[0])
    max_j0 = lon.size - 2
    max_i0 = lat.size - 2
    out: Dict[str, Tuple[int, int, float, float]] = {}
    for g in gauges:
        fx = (float(g.lon) - lon0) / dlon
        fy = (float(g.lat) - lat0) / dlat
        j0 = int(np.floor(fx))
        i0 = int(np.floor(fy))
        if j0 < 0:
            j0 = 0
            wx = 0.0
        elif j0 > max_j0:
            j0 = max_j0
            wx = 1.0
        else:
            wx = fx - j0
        if i0 < 0:
            i0 = 0
            wy = 0.0
        elif i0 > max_i0:
            i0 = max_i0
            wy = 1.0
        else:
            wy = fy - i0
        out[g.name] = (i0, j0, float(wy), float(wx))
    return out


def bilinear_sample(field: np.ndarray, spec: Tuple[int, int, float, float]) -> float:
    i0, j0, wy, wx = spec
    v00 = float(field[i0, j0])
    v01 = float(field[i0, j0 + 1])
    v10 = float(field[i0 + 1, j0])
    v11 = float(field[i0 + 1, j0 + 1])
    top = (1.0 - wx) * v00 + wx * v01
    bot = (1.0 - wx) * v10 + wx * v11
    return float((1.0 - wy) * top + wy * bot)

def prepare_domain(cfg: ProjectConfig) -> PreparedDomain:
    grid, topo_local, local_merged = load_combined_domain(cfg)
    z = grid.z.astype(np.float64)
    for _ in range(2):
        z = 0.5 * z + 0.125 * (
            np.roll(z, 1, axis=0) + np.roll(z, -1, axis=0) + np.roll(z, 1, axis=1) + np.roll(z, -1, axis=1)
        )
    H = np.maximum(-z, 0.0)

    if local_merged is not None:
        lon_min, lon_max = float(local_merged.lon[0]), float(local_merged.lon[-1])
        lat_min, lat_max = float(local_merged.lat[0]), float(local_merged.lat[-1])
    else:
        lon_min, lon_max = cfg.domain.dapa_lon_min, cfg.domain.dapa_lon_max
        lat_min, lat_max = cfg.domain.dapa_lat_min, cfg.domain.dapa_lat_max

    dapa_ys, dapa_xs = find_bbox_indices(grid.lon, grid.lat, lon_min, lon_max, lat_min, lat_max)
    dapa_info = describe_subgrid(grid, H, dapa_ys, dapa_xs)
    attempt = 0
    target_corridor_wet = 64
    while dapa_info["wet_cells_gt_min"] < target_corridor_wet and attempt < 6:
        lon_min, lon_max, lat_min, lat_max = expand_bbox(lon_min, lon_max, lat_min, lat_max, factor=1.40)
        dapa_ys, dapa_xs = find_bbox_indices(grid.lon, grid.lat, lon_min, lon_max, lat_min, lat_max)
        dapa_info = describe_subgrid(grid, H, dapa_ys, dapa_xs)
        attempt += 1

    gauges = sample_gauges(grid, H, cfg.gauges, min_depth=max(cfg.solver.min_depth_m, 0.50))
    gauge_interp = build_gauge_interp(grid, cfg.gauges)
    print(
        f"[DAPA] lon=[{grid.lon[dapa_xs][0]:.5f},{grid.lon[dapa_xs][-1]:.5f}] lat=[{grid.lat[dapa_ys][0]:.5f},{grid.lat[dapa_ys][-1]:.5f}] "
        f"shape={H[dapa_ys, dapa_xs].shape} wet_cells={dapa_info['wet_cells']} wet_cells_gt_min={dapa_info['wet_cells_gt_min']} "
        f"max_static_depth={dapa_info['max_static_depth']:.2f}m",
        flush=True,
    )
    for g in cfg.gauges:
        ii, jj = gauges[g.name]
        gi0, gj0, gwy, gwx = gauge_interp[g.name]
        support = [(gi0, gj0), (gi0, gj0 + 1), (gi0 + 1, gj0), (gi0 + 1, gj0 + 1)]
        support_wet = sum(1 for si, sj in support if H[si, sj] > max(cfg.solver.min_depth_m, 0.50))
        print(
            f"[GAUGE] {g.name}: target=({g.lon:.5f},{g.lat:.5f}) interp_support=({grid.lon[gj0]:.5f}..{grid.lon[gj0+1]:.5f},{grid.lat[gi0]:.5f}..{grid.lat[gi0+1]:.5f}) "
            f"nearest_wet=({grid.lon[jj]:.5f},{grid.lat[ii]:.5f}) H={H[ii, jj]:.2f}m support_wet={support_wet}/4",
            flush=True,
        )

    if topo_local is not None:
        topo_local = crop_grid(
            topo_local,
            grid.lon[dapa_xs][0],
            grid.lon[dapa_xs][-1],
            grid.lat[dapa_ys][0],
            grid.lat[dapa_ys][-1],
        )

    nested_local_domain = None
    if cfg.execution.local_dapa_nested and local_merged is not None:
        local_z = local_merged.z.astype(np.float64)
        local_H = np.maximum(-local_z, 0.0)
        local_gauges = sample_gauges(
            local_merged,
            local_H,
            (GaugePoint("Dapa_virtual", DAPA_MARKER[0], DAPA_MARKER[1]),),
            min_depth=max(cfg.solver.min_depth_m, 0.10),
        )
        local_gauge_interp = build_gauge_interp(local_merged, (GaugePoint("Dapa_virtual", DAPA_MARKER[0], DAPA_MARKER[1]),))
        nested_local_domain = PreparedDomain(
            grid=local_merged,
            H=local_H,
            z=local_z,
            dx=local_merged.dx_m,
            dy=local_merged.dy_m,
            gauge_idx=local_gauges,
            gauge_interp=local_gauge_interp,
            dapa_ys=slice(0, local_merged.shape[0]),
            dapa_xs=slice(0, local_merged.shape[1]),
            topo_grid_local=None,
            nested_local_domain=None,
        )
        print(
            f"[LOCAL NEST] shape={local_merged.shape} dx={nested_local_domain.dx:.1f}m dy={nested_local_domain.dy:.1f}m "
            f"wet_cells={int(np.count_nonzero(local_H > cfg.solver.min_depth_m))}",
            flush=True,
        )

    return PreparedDomain(
        grid=grid,
        H=H,
        z=z,
        dx=grid.dx_m,
        dy=grid.dy_m,
        gauge_idx=gauges,
        gauge_interp=gauge_interp,
        dapa_ys=dapa_ys,
        dapa_xs=dapa_xs,
        topo_grid_local=topo_local,
        nested_local_domain=nested_local_domain,
    )


# ---------- stochastic source ----------

def target_seismic_moment(mw: float) -> float:
    return 10.0 ** (1.5 * mw + 9.1)


def estimate_mean_log10_slip(cfg: ProjectConfig) -> float:
    area_m2 = cfg.source.rupture_length_km * 1000.0 * cfg.source.rupture_width_km * 1000.0
    moment_based = target_seismic_moment(cfg.source.mw) / (cfg.source.shear_modulus_pa * area_m2)
    prescribed = float(max(cfg.source.mean_slip_m, 1.0e-4))
    mean_slip_m = 0.5 * moment_based + 0.5 * prescribed
    return float(np.log10(max(mean_slip_m, 1e-4)))


def subfault_centers(cfg: ProjectConfig) -> np.ndarray:
    x = np.linspace(0.5, cfg.slip.n_along - 0.5, cfg.slip.n_along) * (cfg.source.rupture_length_km / cfg.slip.n_along)
    y = np.linspace(0.5, cfg.slip.n_down - 0.5, cfg.slip.n_down) * (cfg.source.rupture_width_km / cfg.slip.n_down)
    xx, yy = np.meshgrid(x, y)
    return np.column_stack([xx.ravel(), yy.ravel()])


def matern_covariance(distance_km: np.ndarray, corr_length_km: float, nu: float) -> np.ndarray:
    r = np.asarray(distance_km, dtype=float)
    scaled = np.sqrt(2.0 * nu) * np.maximum(r, 1e-12) / max(corr_length_km, 1e-12)
    coeff = (2.0 ** (1.0 - nu)) / special.gamma(nu)
    out = coeff * (scaled ** nu) * special.kv(nu, scaled)
    out[r == 0.0] = 1.0
    out[np.isnan(out)] = 1.0
    return out


def generate_slip_ensemble(cfg: ProjectConfig, n_samples: int) -> np.ndarray:
    centers = subfault_centers(cfg)
    dist_km = spatial.distance_matrix(centers, centers)
    cov = matern_covariance(dist_km, cfg.slip.corr_length_km, cfg.slip.hurst_nu) + 1e-8 * np.eye(centers.shape[0])
    chol = np.linalg.cholesky(cov)
    if qmc is None:
        z = np.random.default_rng(cfg.slip.rqmc_seed).standard_normal((n_samples, centers.shape[0]))
    else:
        engine = qmc.Sobol(d=centers.shape[0], scramble=True, seed=cfg.slip.rqmc_seed)
        u = np.clip(engine.random(n_samples), 1e-12, 1 - 1e-12)
        z = stats.norm.ppf(u)
    raw_log10 = estimate_mean_log10_slip(cfg) + z @ chol.T
    raw_slip = np.power(10.0, raw_log10)
    subfault_area_m2 = (cfg.source.rupture_length_km * 1000.0 / cfg.slip.n_along) * (cfg.source.rupture_width_km * 1000.0 / cfg.slip.n_down)
    m0_target = target_seismic_moment(cfg.source.mw)
    out = []
    for i in range(n_samples):
        m0_raw = cfg.source.shear_modulus_pa * subfault_area_m2 * float(np.sum(raw_slip[i]))
        out.append((raw_slip[i] * (m0_target / max(m0_raw, 1e-12))).astype(np.float32))
    return np.asarray(out)


def discover_slip_files_for_mw(cfg: ProjectConfig, mw: float) -> List[Path]:
    files = list_files(cfg.slip_dir, [".npy", ".npz", ".csv"])
    tags = [f"mw{mw:.1f}", f"m{mw:.1f}", f"mw{str(mw).replace('.', '')}"]
    tagged = [p for p in files if any(t in p.as_posix().lower() for t in tags)]
    return tagged if tagged else files


def load_slip_sample(path: Path) -> np.ndarray:
    if path.suffix.lower() == ".npy":
        return np.asarray(np.load(path, allow_pickle=True), dtype=np.float32)
    if path.suffix.lower() == ".npz":
        data = np.load(path)
        key = next((k for k in ("slip", "arr_0", "uplift", "uz") if k in data), None)
        if key is None:
            raise ValueError(f"No usable array found in {path}")
        return np.asarray(data[key], dtype=np.float32)
    return pd.read_csv(path, header=None).values.astype(np.float32)


def reshape_slip(cfg: ProjectConfig, slip: np.ndarray) -> np.ndarray:
    slip = np.asarray(slip, dtype=float)
    if slip.ndim == 1:
        return slip.reshape(cfg.slip.n_down, cfg.slip.n_along)
    if slip.shape == (cfg.slip.n_along, cfg.slip.n_down):
        return slip.T
    return slip


def build_uplift_field(cfg: ProjectConfig, domain: PreparedDomain, slip: np.ndarray) -> np.ndarray:
    slip_grid = reshape_slip(cfg, slip)
    uplift_sub = slip_grid * cfg.source.uplift_fraction_of_slip * math.sin(math.radians(cfg.source.dip_deg))
    field = np.zeros(domain.grid.shape, dtype=np.float64)
    ys, xs = find_bbox_indices(domain.grid.lon, domain.grid.lat, 125.7, 126.8, 7.2, 10.8)
    field[ys, xs] = bilinear_resize(uplift_sub, field[ys, xs].shape)
    for _ in range(3):
        field = 0.4 * field + 0.15 * (
            np.roll(field, 1, 0) + np.roll(field, -1, 0) + np.roll(field, 1, 1) + np.roll(field, -1, 1)
        )
    field = np.where(domain.H > 20.0, field, 0.0)
    return np.clip(field, -10.0, 10.0)


def source_time_derivative(t: float, rise_time: float) -> float:
    if t <= 0 or t >= rise_time:
        return 0.0
    return 0.5 * math.pi / rise_time * math.sin(math.pi * t / rise_time)


# ---------- numerics ----------

def compute_laplacian(a: np.ndarray, dx: float, dy: float) -> np.ndarray:
    lap = np.zeros_like(a)
    lap[1:-1, 1:-1] = ((a[1:-1, 2:] - 2 * a[1:-1, 1:-1] + a[1:-1, :-2]) / (dx * dx) +
                       (a[2:, 1:-1] - 2 * a[1:-1, 1:-1] + a[:-2, 1:-1]) / (dy * dy))
    return lap


def compute_rhs(q: np.ndarray, domain: PreparedDomain, source_rate: np.ndarray, cfg: ProjectConfig) -> np.ndarray:
    eta, hu, hv = q
    dx, dy, g = domain.dx, domain.dy, cfg.solver.gravity
    h_total = np.maximum(eta - domain.z, 0.0)
    active = h_total > cfg.solver.min_depth_m

    rhs = np.zeros_like(q)
    deta_dx = np.zeros_like(eta)
    deta_dy = np.zeros_like(eta)
    deta_dx[:, 1:-1] = (eta[:, 2:] - eta[:, :-2]) / (2.0 * dx)
    deta_dy[1:-1, :] = (eta[2:, :] - eta[:-2, :]) / (2.0 * dy)
    dhu_dx = np.zeros_like(hu)
    dhv_dy = np.zeros_like(hv)
    dhu_dx[:, 1:-1] = (hu[:, 2:] - hu[:, :-2]) / (2.0 * dx)
    dhv_dy[1:-1, :] = (hv[2:, :] - hv[:-2, :]) / (2.0 * dy)

    rhs[0] = -(dhu_dx + dhv_dy) + source_rate
    rhs[1][active] = -(g * h_total * deta_dx)[active] - cfg.solver.damping * hu[active]
    rhs[2][active] = -(g * h_total * deta_dy)[active] - cfg.solver.damping * hv[active]

    nu = cfg.solver.viscosity_factor * min(dx, dy)
    rhs[0] += nu * compute_laplacian(eta, dx, dy)
    rhs[1][active] += nu * compute_laplacian(hu, dx, dy)[active]
    rhs[2][active] += nu * compute_laplacian(hv, dx, dy)[active]
    rhs[1][~active] = 0.0
    rhs[2][~active] = 0.0
    return rhs




def use_numba_backend(cfg: ProjectConfig) -> bool:
    backend = (cfg.execution.kernel_backend or "auto").lower()
    if backend == "numpy":
        return False
    if backend == "numba":
        if not NUMBA_AVAILABLE:
            raise RuntimeError("kernel_backend=numba requested but numba is not available.")
        return True
    return NUMBA_AVAILABLE


@njit(cache=True)
def _rhs_inplace_numba(eta, hu, hv, z, uplift, source_scale, dx, dy, g, min_depth, damping, nu, rhs0, rhs1, rhs2, h_total, wet):
    ny, nx = eta.shape
    inv2dx = 0.5 / dx
    inv2dy = 0.5 / dy
    invdx2 = 1.0 / (dx * dx)
    invdy2 = 1.0 / (dy * dy)

    for i in range(ny):
        for j in range(nx):
            h = eta[i, j] - z[i, j]
            if h < 0.0:
                h = 0.0
            h_total[i, j] = h
            wet[i, j] = h > min_depth
            rhs0[i, j] = source_scale * uplift[i, j]
            rhs1[i, j] = 0.0
            rhs2[i, j] = 0.0

    for i in range(1, ny - 1):
        for j in range(1, nx - 1):
            deta_dx = (eta[i, j + 1] - eta[i, j - 1]) * inv2dx
            deta_dy = (eta[i + 1, j] - eta[i - 1, j]) * inv2dy
            dhu_dx = (hu[i, j + 1] - hu[i, j - 1]) * inv2dx
            dhv_dy = (hv[i + 1, j] - hv[i - 1, j]) * inv2dy

            lap_eta = ((eta[i, j + 1] - 2.0 * eta[i, j] + eta[i, j - 1]) * invdx2 +
                       (eta[i + 1, j] - 2.0 * eta[i, j] + eta[i - 1, j]) * invdy2)
            rhs0[i, j] += -(dhu_dx + dhv_dy) + nu * lap_eta

            if wet[i, j]:
                lap_hu = ((hu[i, j + 1] - 2.0 * hu[i, j] + hu[i, j - 1]) * invdx2 +
                          (hu[i + 1, j] - 2.0 * hu[i, j] + hu[i - 1, j]) * invdy2)
                lap_hv = ((hv[i, j + 1] - 2.0 * hv[i, j] + hv[i, j - 1]) * invdx2 +
                          (hv[i + 1, j] - 2.0 * hv[i, j] + hv[i - 1, j]) * invdy2)
                rhs1[i, j] = -(g * h_total[i, j] * deta_dx) - damping * hu[i, j] + nu * lap_hu
                rhs2[i, j] = -(g * h_total[i, j] * deta_dy) - damping * hv[i, j] + nu * lap_hv


@njit(cache=True)
def _advance_step_numba(eta, hu, hv, z, uplift, src_scale1, src_scale2, dt, t_now, dx, dy, g, min_depth, damping, nu, land_mask,
                        k10, k11, k12, k20, k21, k22, eta1, hu1, hv1, h_total, wet,
                        max_eta, max_depth, max_speed, arrival_time, wet_mask):
    ny, nx = eta.shape
    _rhs_inplace_numba(eta, hu, hv, z, uplift, src_scale1, dx, dy, g, min_depth, damping, nu, k10, k11, k12, h_total, wet)

    for i in range(ny):
        for j in range(nx):
            eta1[i, j] = eta[i, j] + dt * k10[i, j]
            hu1[i, j] = hu[i, j] + dt * k11[i, j]
            hv1[i, j] = hv[i, j] + dt * k12[i, j]

    _rhs_inplace_numba(eta1, hu1, hv1, z, uplift, src_scale2, dx, dy, g, min_depth, damping, nu, k20, k21, k22, h_total, wet)

    max_abs_eta = 0.0
    finite_ok = True
    for i in range(ny):
        for j in range(nx):
            e = eta[i, j] + 0.5 * dt * (k10[i, j] + k20[i, j])
            u_m = hu[i, j] + 0.5 * dt * (k11[i, j] + k21[i, j])
            v_m = hv[i, j] + 0.5 * dt * (k12[i, j] + k22[i, j])
            eta[i, j] = e
            hu[i, j] = u_m
            hv[i, j] = v_m

            h = e - z[i, j]
            if h < 0.0:
                h = 0.0
            is_wet = h > min_depth
            sp = 0.0
            if is_wet:
                sp = (u_m / h) * (u_m / h) + (v_m / h) * (v_m / h)
                if sp > 0.0:
                    sp = np.sqrt(sp)
            else:
                hu[i, j] = 0.0
                hv[i, j] = 0.0

            if not np.isfinite(e) or not np.isfinite(sp):
                finite_ok = False

            ae = abs(e)
            if ae > max_abs_eta:
                max_abs_eta = ae
            if e > max_eta[i, j]:
                max_eta[i, j] = e
            if land_mask[i, j]:
                if h > max_depth[i, j]:
                    max_depth[i, j] = h
                if is_wet and np.isnan(arrival_time[i, j]):
                    arrival_time[i, j] = np.float32(t_now)
            if sp > 1.0e4:
                sp = 1.0e4
            if sp > max_speed[i, j]:
                max_speed[i, j] = sp
            if is_wet:
                wet_mask[i, j] = 1

    return max_abs_eta, finite_ok


@njit(cache=True)
def _precompute_interp_coords(axis_vals, targets):
    n = targets.shape[0]
    idx0 = np.empty(n, dtype=np.int64)
    w = np.empty(n, dtype=np.float64)
    start = axis_vals[0]
    step = axis_vals[1] - axis_vals[0]
    nmax = axis_vals.shape[0] - 2
    for k in range(n):
        fx = (targets[k] - start) / step
        i0 = int(np.floor(fx))
        if i0 < 0:
            i0 = 0
            frac = 0.0
        elif i0 > nmax:
            i0 = nmax
            frac = 1.0
        else:
            frac = fx - i0
        idx0[k] = i0
        w[k] = frac
    return idx0, w


@njit(cache=True)
def _sample_series_bilinear(snaps, iy0, wy, ix0, wx):
    nt = snaps.shape[0]
    npnt = iy0.shape[0]
    out = np.empty((nt, npnt), dtype=np.float32)
    for t in range(nt):
        s = snaps[t]
        for p in range(npnt):
            i0 = iy0[p]
            j0 = ix0[p]
            fy = wy[p]
            fx = wx[p]
            v00 = s[i0, j0]
            v01 = s[i0, j0 + 1]
            v10 = s[i0 + 1, j0]
            v11 = s[i0 + 1, j0 + 1]
            top = (1.0 - fx) * v00 + fx * v01
            bot = (1.0 - fx) * v10 + fx * v11
            out[t, p] = np.float32((1.0 - fy) * top + fy * bot)
    return out

def stable_dt(domain: PreparedDomain, cfg: ProjectConfig) -> float:
    Hwet = domain.H[domain.H > max(cfg.solver.min_depth_m, 1.0)]
    if Hwet.size == 0:
        return min(cfg.solver.dt_max, 0.05)
    cmax = float(np.sqrt(cfg.solver.gravity * np.percentile(Hwet, 99.5)))
    dt = cfg.solver.cfl / max(cmax / domain.dx, cmax / domain.dy)
    return float(min(dt, cfg.solver.dt_max))


def maybe_print(msg: str, pbar=None) -> None:
    if pbar is not None and hasattr(pbar, "write"):
        pbar.write(msg)
    else:
        print(msg, flush=True)


def make_progress(total: int, desc: str):
    if tqdm is not None:
        return tqdm(total=total, desc=desc, leave=True)

    class Dummy:
        def __init__(self, total, desc):
            self.total, self.n, self.desc = total, 0, desc
        def update(self, n=1):
            self.n += n
            print(f"{self.desc}: {self.n}/{self.total}", flush=True)
        def write(self, msg):
            print(msg, flush=True)
        def close(self):
            pass
    return Dummy(total, desc)


# ---------- mapping / raster export ----------

def grid_transform(grid: Grid2D):
    if from_origin is None:
        return None
    dx = abs(float(np.median(np.diff(grid.lon)))) if len(grid.lon) > 1 else 1.0e-4
    dy = abs(float(np.median(np.diff(grid.lat)))) if len(grid.lat) > 1 else 1.0e-4
    west = float(grid.lon[0] - 0.5 * dx)
    north = float(grid.lat[-1] + 0.5 * dy)
    return from_origin(west, north, dx, dy)


def export_geotiff(path: Path, grid: Grid2D, arr: np.ndarray, dtype: str = "float32", nodata: Optional[float] = None) -> None:
    if rasterio is None:
        return
    ensure_dir(path.parent)
    transform = grid_transform(grid)
    data = np.asarray(arr)
    data_write = data[::-1, :].astype(dtype)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=data_write.shape[0],
        width=data_write.shape[1],
        count=1,
        dtype=data_write.dtype,
        crs="EPSG:4326",
        transform=transform,
        nodata=nodata,
        compress="deflate",
    ) as dst:
        dst.write(data_write, 1)


def save_npz_map(path: Path, grid: Grid2D, arr: np.ndarray, variable_name: str) -> None:
    ensure_dir(path.parent)
    np.savez_compressed(path, lon=grid.lon, lat=grid.lat, **{variable_name: np.asarray(arr, dtype=np.float32)})


def classify_inundation_depth(depth: np.ndarray) -> np.ndarray:
    cls = np.zeros(depth.shape, dtype=np.uint8)
    for k in range(len(PAPER_INUNDATION_BOUNDS) - 1):
        lo = PAPER_INUNDATION_BOUNDS[k]
        hi = PAPER_INUNDATION_BOUNDS[k + 1]
        cls[(depth > lo) & (depth <= hi)] = k + 1
    cls[depth > PAPER_INUNDATION_BOUNDS[-1]] = len(PAPER_INUNDATION_BOUNDS) - 1
    return cls


def export_map_bundle(out_dir: Path, prefix: str, grid: Grid2D, depth: np.ndarray, eta: np.ndarray, speed: np.ndarray,
                      arrival_time: np.ndarray, wet_mask: np.ndarray, write_geotiff: bool = True, write_npz: bool = True) -> Dict[str, str]:
    ensure_dir(out_dir)
    depth_cls = classify_inundation_depth(depth)
    outputs: Dict[str, str] = {}
    if write_geotiff:
        export_geotiff(out_dir / f"{prefix}_max_depth.tif", grid, depth, dtype="float32", nodata=np.nan)
        export_geotiff(out_dir / f"{prefix}_max_eta.tif", grid, eta, dtype="float32", nodata=np.nan)
        export_geotiff(out_dir / f"{prefix}_max_speed.tif", grid, speed, dtype="float32", nodata=np.nan)
        export_geotiff(out_dir / f"{prefix}_arrival_time.tif", grid, np.where(np.isfinite(arrival_time), arrival_time, np.nan), dtype="float32", nodata=np.nan)
        export_geotiff(out_dir / f"{prefix}_wet_mask.tif", grid, wet_mask.astype(np.uint8), dtype="uint8", nodata=0)
        export_geotiff(out_dir / f"{prefix}_depth_classes.tif", grid, depth_cls.astype(np.uint8), dtype="uint8", nodata=0)
        outputs.update({
            "max_depth_tif": str(out_dir / f"{prefix}_max_depth.tif"),
            "max_eta_tif": str(out_dir / f"{prefix}_max_eta.tif"),
            "max_speed_tif": str(out_dir / f"{prefix}_max_speed.tif"),
            "arrival_time_tif": str(out_dir / f"{prefix}_arrival_time.tif"),
            "wet_mask_tif": str(out_dir / f"{prefix}_wet_mask.tif"),
            "depth_classes_tif": str(out_dir / f"{prefix}_depth_classes.tif"),
        })
    if write_npz:
        save_npz_map(out_dir / f"{prefix}_max_depth.npz", grid, depth, "depth_max")
        save_npz_map(out_dir / f"{prefix}_max_eta.npz", grid, eta, "eta_max")
        save_npz_map(out_dir / f"{prefix}_max_speed.npz", grid, speed, "speed_max")
        save_npz_map(out_dir / f"{prefix}_arrival_time.npz", grid, np.where(np.isfinite(arrival_time), arrival_time, -1.0), "arrival_time")
        save_npz_map(out_dir / f"{prefix}_wet_mask.npz", grid, wet_mask.astype(np.uint8), "wet_mask")
        save_npz_map(out_dir / f"{prefix}_depth_classes.npz", grid, depth_cls.astype(np.uint8), "depth_class")
    save_json(out_dir / f"{prefix}_fixed_extent.json", {
        "lon_min": float(grid.lon[0]),
        "lon_max": float(grid.lon[-1]),
        "lat_min": float(grid.lat[0]),
        "lat_max": float(grid.lat[-1]),
        "nx": int(grid.shape[1]),
        "ny": int(grid.shape[0]),
        "pixel_dx_deg": float(np.median(np.diff(grid.lon))) if len(grid.lon) > 1 else None,
        "pixel_dy_deg": float(np.median(np.diff(grid.lat))) if len(grid.lat) > 1 else None,
        "paper_inundation_bounds_m": PAPER_INUNDATION_BOUNDS.tolist(),
    })
    return outputs


def _save_paper_style_matplotlib(grid: Grid2D, arr: np.ndarray, title: str, path: Path) -> None:
    if plt is None or ListedColormap is None or BoundaryNorm is None:
        return
    ensure_dir(path.parent)
    fig, ax = plt.subplots(figsize=(8, 5), dpi=180)
    extent = safe_extent(grid.lon, grid.lat)

    background = np.where(grid.z > 0.0, 1.0, 0.0)
    bg_cmap = ListedColormap(["#cfe8f3", "#e6e6e6"])
    ax.imshow(background, origin="lower", extent=extent, aspect="auto", cmap=bg_cmap, interpolation="nearest")

    masked = np.ma.masked_less_equal(np.asarray(arr, dtype=float), PAPER_INUNDATION_BOUNDS[0] - 1.0e-9)
    cmap = ListedColormap(PAPER_INUNDATION_COLORS)
    norm = BoundaryNorm(PAPER_INUNDATION_BOUNDS, cmap.N)
    im = ax.imshow(masked, origin="lower", extent=extent, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest", alpha=0.96)

    try:
        ax.contour(grid.lon, grid.lat, grid.z, levels=[0.0], colors="black", linewidths=0.45)
    except Exception:
        pass

    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    if extent[0] <= DAPA_MARKER[0] <= extent[1] and extent[2] <= DAPA_MARKER[1] <= extent[3]:
        ax.plot(DAPA_MARKER[0], DAPA_MARKER[1], marker="*", ms=7, color="red", mec="black", mew=0.3)
        ax.text(DAPA_MARKER[0] + 0.005, DAPA_MARKER[1] + 0.004, "Dapa", fontsize=8, weight="bold")

    cbar = fig.colorbar(im, ax=ax, shrink=0.92)
    cbar.set_label("Inundation depth (m)")
    cbar.set_ticks(PAPER_INUNDATION_BOUNDS[:-1] + np.diff(PAPER_INUNDATION_BOUNDS) / 2.0)
    cbar.set_ticklabels([
        "0.01–0.25", "0.25–0.50", "0.50–1.00", "1.00–1.50",
        "1.50–2.00", "2.00–2.50", "2.50–3.00", "3.00–3.50"
    ])
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def save_map(grid: Grid2D, arr: np.ndarray, title: str, path: Path, paper_style: bool = False) -> None:
    ensure_dir(path.parent)
    lon, lat = grid.lon, grid.lat
    if paper_style:
        _save_paper_style_matplotlib(grid, arr, title, path)
        return
    if pygmt is not None and xr is not None:
        try:
            da = xr.DataArray(np.asarray(arr, dtype=float), coords={"lat": lat, "lon": lon}, dims=("lat", "lon"))
            fig = pygmt.Figure()
            region = [float(lon[0]), float(lon[-1]), float(lat[0]), float(lat[-1])]
            fig.basemap(region=region, projection="M12c", frame=["af", f"+t{title}"])
            fig.grdimage(grid=da, cmap="turbo")
            fig.coast(shorelines="0.25p,black", borders="1/0.2p,black")
            fig.colorbar(frame="af+lValue")
            if region[0] <= DAPA_MARKER[0] <= region[1] and region[2] <= DAPA_MARKER[1] <= region[3]:
                fig.plot(x=[DAPA_MARKER[0]], y=[DAPA_MARKER[1]], style="c0.18c", fill="red", pen="0.3p,black")
                fig.text(x=DAPA_MARKER[0], y=DAPA_MARKER[1], text="Dapa", font="8p,Helvetica-Bold,black", justify="LM", offset="0.12c/0.12c")
            fig.savefig(str(path))
            return
        except Exception:
            pass
    if plt is None:
        return
    fig, ax = plt.subplots(figsize=(7, 5), dpi=160)
    im = ax.imshow(arr, origin="lower", extent=safe_extent(lon, lat), aspect="auto")
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    if float(lon[0]) <= DAPA_MARKER[0] <= float(lon[-1]) and float(lat[0]) <= DAPA_MARKER[1] <= float(lat[-1]):
        ax.plot(DAPA_MARKER[0], DAPA_MARKER[1], "ro", ms=3)
        ax.text(DAPA_MARKER[0] + 0.005, DAPA_MARKER[1] + 0.004, "Dapa", fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def render_local_dapa_map(topo_local: Optional[Grid2D], coarse_depth: np.ndarray, title: str, path: Path) -> None:
    if topo_local is None:
        return
    positive = np.count_nonzero(np.asarray(coarse_depth) > 0.0)
    if positive <= 0:
        return
    rendered = bilinear_resize(coarse_depth, topo_local.z.shape)
    rendered = np.where(topo_local.z > 0.0, np.maximum(rendered, 0.0), 0.0)
    if np.count_nonzero(rendered > 0.0) <= 0:
        return
    save_map(topo_local, rendered, title, path, paper_style=True)


# ---------- hazard helpers ----------

def exceedance_map(stack: np.ndarray, thr: float) -> np.ndarray:
    return (stack > thr).mean(axis=0).astype(np.float32)


def magnitude_weights(beta: float, magnitudes: Sequence[float]) -> Dict[float, float]:
    mags = np.asarray(magnitudes, dtype=float)
    w = beta * np.exp(-beta * mags)
    w /= np.sum(w)
    return {float(m): float(wi) for m, wi in zip(mags, w)}


def convergence_schedule(target_n: int) -> List[int]:
    base = [8, 16, 32, 64, 128, 256, 512, 1024, 2048]
    return [n for n in base if n <= target_n]


def apply_magnitude_scenario(cfg: ProjectConfig, mw: float) -> None:
    cfg.source.mw = mw
    params = {
        8.0: dict(rupture_length_km=200.0, rupture_width_km=100.0, top_depth_km=7.6, strike_deg=160.0, dip_deg=40.0, rake_deg=90.0, mean_slip_m=2.7),
        8.5: dict(rupture_length_km=300.0, rupture_width_km=100.0, top_depth_km=7.6, strike_deg=164.0, dip_deg=39.0, rake_deg=90.0, mean_slip_m=4.7),
        9.0: dict(rupture_length_km=600.0, rupture_width_km=100.0, top_depth_km=6.1, strike_deg=158.0, dip_deg=33.0, rake_deg=90.0, mean_slip_m=8.1),
        7.6: dict(rupture_length_km=120.0, rupture_width_km=60.0, top_depth_km=10.0, strike_deg=164.0, dip_deg=39.0, rake_deg=90.0, mean_slip_m=1.2),
    }[round(mw, 1)]
    for k, v in params.items():
        setattr(cfg.source, k, v)


def gather_slips(cfg: ProjectConfig, mw: float, n_needed: int) -> List[Tuple[np.ndarray, str]]:
    apply_magnitude_scenario(cfg, mw)
    files = discover_slip_files_for_mw(cfg, mw)
    out: List[Tuple[np.ndarray, str]] = []
    for p in files[:n_needed]:
        out.append((load_slip_sample(p), p.stem))
    if len(out) < n_needed:
        missing = n_needed - len(out)
        print(f"[WARN] Only found {len(out)} slip files for Mw{mw:.1f}. Generating {missing} fallback samples.", flush=True)
        ens = generate_slip_ensemble(cfg, missing)
        for i in range(missing):
            out.append((ens[i], f"generated_Mw{mw:.1f}_{i}"))
    return out


def save_gauge_series(out_dir: Path, sample: SampleOutput) -> Path:
    df = pd.DataFrame({"time_sec": sample.gauge_time_sec})
    for k, v in sample.gauge_eta.items():
        df[k] = v
    out_path = out_dir / f"gauges_{sample.slip_id}_{sample.forcing_mode}.csv"
    df.to_csv(out_path, index=False)
    return out_path


def save_checkpoint(out_dir: Path, mw: float, sample_index: int, gauge_rows: List[Dict[str, Any]], p95_rows: List[Dict[str, Any]]) -> None:
    pd.DataFrame(gauge_rows).to_csv(out_dir / f"Mw{mw:.1f}_gauge_running_summary_checkpoint_{sample_index + 1}.csv", index=False)
    pd.DataFrame(p95_rows).to_csv(out_dir / f"Mw{mw:.1f}_dapa_p95_running_checkpoint_{sample_index + 1}.csv", index=False)
    save_json(out_dir / f"Mw{mw:.1f}_checkpoint_{sample_index + 1}.json", {"mw": mw, "samples_completed": sample_index + 1})


# ---------- simulation cores ----------

def run_tsunami_simulation(
    cfg: ProjectConfig,
    domain: PreparedDomain,
    slip: np.ndarray,
    slip_id: str,
    sample_index: int,
    forcing_mode: str = "dynamic",
    pbar=None,
) -> SampleOutput:
    ny, nx = domain.grid.shape
    q = np.zeros((3, ny, nx), dtype=np.float64)
    eta, hu, hv = q[0], q[1], q[2]
    max_eta = np.full((ny, nx), -np.inf, dtype=np.float32)
    max_depth = np.zeros((ny, nx), dtype=np.float32)
    max_speed = np.zeros((ny, nx), dtype=np.float32)
    arrival_time = np.full((ny, nx), np.nan, dtype=np.float32)
    wet_mask = np.zeros((ny, nx), dtype=np.uint8)

    uplift0 = build_uplift_field(cfg, domain, slip)
    rise_time = cfg.source.rise_time_sec_per_km * cfg.source.rupture_length_km
    dt0 = stable_dt(domain, cfg)

    if forcing_mode == "static":
        eta[:, :] = uplift0

    if sample_index == 0:
        maybe_print(f"Grid: {ny} x {nx}  dx={domain.dx:.0f}m", pbar)
        maybe_print(f"dt={dt0:.4f}s", pbar)
        maybe_print("", pbar)

    dapa_static = domain.H[domain.dapa_ys, domain.dapa_xs]
    maybe_print(
        f"[DIAG] Sample {sample_index}: slip_max={float(np.max(slip)):.4f}, domain_depth_max={float(np.max(domain.H)):.2f}m, "
        f"dapa_static_depth_max={float(np.max(dapa_static)):.2f}m, dapa_static_wet_cells={int(np.count_nonzero(dapa_static > cfg.solver.min_depth_m))}",
        pbar,
    )

    gauge_time: List[float] = []
    gauge_eta: Dict[str, List[float]] = {g.name: [] for g in cfg.gauges}
    nested_times: List[float] = []
    nested_eta_boxes: List[np.ndarray] = []
    t = 0.0
    step = 0
    aborted = False
    abort_reason = ""
    land_mask = np.ascontiguousarray(domain.z > 0.0)
    use_numba = use_numba_backend(cfg)

    if use_numba:
        zero_uplift = np.zeros((ny, nx), dtype=np.float64)
        k10 = np.zeros((ny, nx), dtype=np.float64)
        k11 = np.zeros((ny, nx), dtype=np.float64)
        k12 = np.zeros((ny, nx), dtype=np.float64)
        k20 = np.zeros((ny, nx), dtype=np.float64)
        k21 = np.zeros((ny, nx), dtype=np.float64)
        k22 = np.zeros((ny, nx), dtype=np.float64)
        eta1 = np.zeros((ny, nx), dtype=np.float64)
        hu1 = np.zeros((ny, nx), dtype=np.float64)
        hv1 = np.zeros((ny, nx), dtype=np.float64)
        h_total = np.zeros((ny, nx), dtype=np.float64)
        wet = np.zeros((ny, nx), dtype=np.bool_)
        nu = cfg.solver.viscosity_factor * min(domain.dx, domain.dy)
    else:
        h_total = None
        wet = None

    while t < cfg.solver.t_end_sec:
        dt = min(dt0, cfg.solver.t_end_sec - t)
        if forcing_mode == "dynamic":
            src_scale1 = source_time_derivative(t, rise_time)
            src_scale2 = source_time_derivative(t + dt, rise_time)
        else:
            src_scale1 = 0.0
            src_scale2 = 0.0

        if use_numba:
            maxabseta, finite_ok = _advance_step_numba(
                eta, hu, hv, domain.z, uplift0, src_scale1, src_scale2, dt, t, domain.dx, domain.dy,
                cfg.solver.gravity, cfg.solver.min_depth_m, cfg.solver.damping, nu, land_mask,
                k10, k11, k12, k20, k21, k22, eta1, hu1, hv1, h_total, wet,
                max_eta, max_depth, max_speed, arrival_time, wet_mask,
            )
            if (not finite_ok):
                aborted, abort_reason = True, "non-finite state encountered"
                maybe_print(f"[ABORT] Sample {sample_index}: {abort_reason}", pbar)
                break
            if maxabseta > cfg.solver.eta_guard_m:
                aborted, abort_reason = True, f"|eta| exceeded safety guard ({maxabseta:.2f} m)"
                maybe_print(f"[ABORT] Sample {sample_index}: {abort_reason}", pbar)
                break
        else:
            if forcing_mode == "dynamic":
                src1 = uplift0 * src_scale1
                src2 = uplift0 * src_scale2
            else:
                src1 = np.zeros_like(uplift0)
                src2 = np.zeros_like(uplift0)
            k1 = compute_rhs(q, domain, src1, cfg)
            q1 = q + dt * k1
            k2 = compute_rhs(q1, domain, src2, cfg)
            q = q + 0.5 * dt * (k1 + k2)
            eta, hu, hv = q[0], q[1], q[2]

            h_total = np.maximum(eta - domain.z, 0.0)
            wet = h_total > cfg.solver.min_depth_m
            q[1, ~wet] = 0.0
            q[2, ~wet] = 0.0

            u = np.divide(q[1], h_total, out=np.zeros_like(q[1]), where=wet)
            v = np.divide(q[2], h_total, out=np.zeros_like(q[2]), where=wet)
            speed = np.sqrt(u * u + v * v)

            maxabseta = float(np.max(np.abs(eta)))
            if (not np.all(np.isfinite(eta))) or (not np.all(np.isfinite(speed))):
                aborted, abort_reason = True, "non-finite state encountered"
                maybe_print(f"[ABORT] Sample {sample_index}: {abort_reason}", pbar)
                break
            if maxabseta > cfg.solver.eta_guard_m:
                aborted, abort_reason = True, f"|eta| exceeded safety guard ({maxabseta:.2f} m)"
                maybe_print(f"[ABORT] Sample {sample_index}: {abort_reason}", pbar)
                break

            max_eta = np.maximum(max_eta, eta.astype(np.float32))
            inun = np.where(land_mask, h_total, 0.0)
            max_depth = np.maximum(max_depth, inun.astype(np.float32))
            max_speed = np.maximum(max_speed, np.minimum(speed, 1.0e4).astype(np.float32))
            newly_wet = land_mask & wet & np.isnan(arrival_time)
            arrival_time[newly_wet] = np.float32(t)
            wet_mask |= wet.astype(np.uint8)

        if step % cfg.solver.snapshot_stride == 0:
            gauge_time.append(t)
            for name in domain.gauge_idx.keys():
                if domain.gauge_interp:
                    gauge_eta[name].append(bilinear_sample(eta, domain.gauge_interp[name]))
                else:
                    i, j = domain.gauge_idx[name]
                    gauge_eta[name].append(float(eta[i, j]))

        if domain.nested_local_domain is not None and step % max(1, cfg.execution.nested_save_stride) == 0:
            nested_times.append(t)
            nested_eta_boxes.append(np.asarray(eta[domain.dapa_ys, domain.dapa_xs], dtype=np.float32).copy())

        if step > 0 and step % cfg.solver.diag_stride == 0:
            dapa_eta_now = eta[domain.dapa_ys, domain.dapa_xs]
            dapa_h_now = (np.maximum(eta - domain.z, 0.0))[domain.dapa_ys, domain.dapa_xs]
            dapa_z_now = domain.z[domain.dapa_ys, domain.dapa_xs]
            dapa_offshore_eta = float(np.max(np.abs(dapa_eta_now[dapa_z_now <= 0.0]))) if np.any(dapa_z_now <= 0.0) else 0.0
            dapa_inun = float(np.max(np.where(dapa_z_now > 0.0, dapa_h_now, 0.0))) if dapa_h_now.size else 0.0
            maybe_print(
                f"  Step {step:5d}  t={t:8.1f}s  eta=[{float(np.min(eta)):.2f},{float(np.max(eta)):.2f}]m  "
                f"dapa_offshore_eta={dapa_offshore_eta:.2f}m  dapa_inun={dapa_inun:.2f}m",
                pbar,
            )

        t += dt
        step += 1

    gauge_time_arr = np.asarray(gauge_time, dtype=np.float32)
    gauge_eta_arr = {k: np.asarray(v, dtype=np.float32) for k, v in gauge_eta.items()}
    gauge_max_abs = {k: float(np.max(np.abs(v))) if len(v) else 0.0 for k, v in gauge_eta_arr.items()}
    dapa_depth = max_depth[domain.dapa_ys, domain.dapa_xs]
    dapa_land = domain.z[domain.dapa_ys, domain.dapa_xs] > 0.0
    positive = dapa_depth[(dapa_depth > 0.0) & dapa_land]
    dapa_p95 = float(np.percentile(positive, 95)) if positive.size else 0.0

    nested_time_arr = np.asarray(nested_times, dtype=np.float32) if nested_times else None
    nested_eta_arr = np.stack(nested_eta_boxes, axis=0).astype(np.float32) if nested_eta_boxes else None

    return SampleOutput(
        slip_id=slip_id,
        forcing_mode=forcing_mode,
        gauge_time_sec=gauge_time_arr,
        gauge_eta=gauge_eta_arr,
        max_eta=max_eta,
        max_depth=max_depth,
        max_speed=max_speed,
        arrival_time=arrival_time,
        wet_mask=wet_mask,
        gauge_max_abs=gauge_max_abs,
        dapa_p95_inun=dapa_p95,
        meta={
            "dt0_sec": dt0,
            "rise_time_sec": rise_time,
            "dx_m": domain.dx,
            "dy_m": domain.dy,
            "t_end_sec": cfg.solver.t_end_sec,
            "wetting_drying": "moving shoreline using h=max(eta-z,0)",
            "aborted": aborted,
            "abort_reason": abort_reason,
            "forcing_mode": forcing_mode,
            "kernel_backend": "numba" if use_numba else "numpy",
        },
        nested_time_sec=nested_time_arr,
        nested_eta_box=nested_eta_arr,
    )


def build_local_boundary_driver(regional_domain: PreparedDomain, local_domain: PreparedDomain,
                               nested_time_sec: np.ndarray, nested_eta_box: np.ndarray) -> Dict[str, np.ndarray]:
    reg_lon = np.asarray(regional_domain.grid.lon[regional_domain.dapa_xs], dtype=np.float64)
    reg_lat = np.asarray(regional_domain.grid.lat[regional_domain.dapa_ys], dtype=np.float64)
    loc_lon = np.asarray(local_domain.grid.lon, dtype=np.float64)
    loc_lat = np.asarray(local_domain.grid.lat, dtype=np.float64)

    top_x = loc_lon
    top_y = np.full(loc_lon.shape, loc_lat[-1], dtype=np.float64)
    bottom_x = loc_lon
    bottom_y = np.full(loc_lon.shape, loc_lat[0], dtype=np.float64)
    left_x = np.full(loc_lat.shape, loc_lon[0], dtype=np.float64)
    left_y = loc_lat
    right_x = np.full(loc_lat.shape, loc_lon[-1], dtype=np.float64)
    right_y = loc_lat

    top_iy0, top_wy = _precompute_interp_coords(reg_lat, top_y)
    top_ix0, top_wx = _precompute_interp_coords(reg_lon, top_x)
    bottom_iy0, bottom_wy = _precompute_interp_coords(reg_lat, bottom_y)
    bottom_ix0, bottom_wx = _precompute_interp_coords(reg_lon, bottom_x)
    left_iy0, left_wy = _precompute_interp_coords(reg_lat, left_y)
    left_ix0, left_wx = _precompute_interp_coords(reg_lon, left_x)
    right_iy0, right_wy = _precompute_interp_coords(reg_lat, right_y)
    right_ix0, right_wx = _precompute_interp_coords(reg_lon, right_x)

    snaps = np.asarray(nested_eta_box, dtype=np.float32)
    return {
        "time": np.asarray(nested_time_sec, dtype=np.float32),
        "top": _sample_series_bilinear(snaps, top_iy0, top_wy, top_ix0, top_wx),
        "bottom": _sample_series_bilinear(snaps, bottom_iy0, bottom_wy, bottom_ix0, bottom_wx),
        "left": _sample_series_bilinear(snaps, left_iy0, left_wy, left_ix0, left_wx),
        "right": _sample_series_bilinear(snaps, right_iy0, right_wy, right_ix0, right_wx),
    }


def interpolate_boundary_state(driver: Dict[str, np.ndarray], t: float) -> Dict[str, np.ndarray]:
    tt = driver["time"]
    if t <= float(tt[0]):
        return {k: driver[k][0] for k in ("top", "bottom", "left", "right")}
    if t >= float(tt[-1]):
        return {k: driver[k][-1] for k in ("top", "bottom", "left", "right")}
    hi = int(np.searchsorted(tt, t, side="right"))
    lo = max(0, hi - 1)
    t0 = float(tt[lo])
    t1 = float(tt[hi])
    if t1 <= t0:
        return {k: driver[k][lo] for k in ("top", "bottom", "left", "right")}
    w = (t - t0) / (t1 - t0)
    return {k: ((1.0 - w) * driver[k][lo] + w * driver[k][hi]).astype(np.float32) for k in ("top", "bottom", "left", "right")}


def impose_boundary_state(q: np.ndarray, eta_bc: Dict[str, np.ndarray]) -> None:
    q[0, 0, :] = eta_bc["bottom"]
    q[0, -1, :] = eta_bc["top"]
    q[0, :, 0] = eta_bc["left"]
    q[0, :, -1] = eta_bc["right"]
    q[1, 0, :] = 0.0
    q[1, -1, :] = 0.0
    q[1, :, 0] = 0.0
    q[1, :, -1] = 0.0
    q[2, 0, :] = 0.0
    q[2, -1, :] = 0.0
    q[2, :, 0] = 0.0
    q[2, :, -1] = 0.0


def run_local_nested_simulation(cfg: ProjectConfig, local_domain: PreparedDomain,
                                driver: Dict[str, np.ndarray], pbar=None) -> LocalNestedOutput:
    ny, nx = local_domain.grid.shape
    q = np.zeros((3, ny, nx), dtype=np.float64)
    eta, hu, hv = q[0], q[1], q[2]
    max_eta = np.full((ny, nx), -np.inf, dtype=np.float32)
    max_depth = np.zeros((ny, nx), dtype=np.float32)
    max_speed = np.zeros((ny, nx), dtype=np.float32)
    arrival_time = np.full((ny, nx), np.nan, dtype=np.float32)
    wet_mask = np.zeros((ny, nx), dtype=np.uint8)
    dt0 = min(stable_dt(local_domain, cfg), cfg.solver.dt_max)
    t_end = float(driver["time"][-1])
    land_mask = np.ascontiguousarray(local_domain.z > 0.0)
    use_numba = use_numba_backend(cfg)
    t = 0.0

    if use_numba:
        zero_uplift = np.zeros((ny, nx), dtype=np.float64)
        k10 = np.zeros((ny, nx), dtype=np.float64)
        k11 = np.zeros((ny, nx), dtype=np.float64)
        k12 = np.zeros((ny, nx), dtype=np.float64)
        k20 = np.zeros((ny, nx), dtype=np.float64)
        k21 = np.zeros((ny, nx), dtype=np.float64)
        k22 = np.zeros((ny, nx), dtype=np.float64)
        eta1 = np.zeros((ny, nx), dtype=np.float64)
        hu1 = np.zeros((ny, nx), dtype=np.float64)
        hv1 = np.zeros((ny, nx), dtype=np.float64)
        h_total = np.zeros((ny, nx), dtype=np.float64)
        wet = np.zeros((ny, nx), dtype=np.bool_)
        nu = cfg.solver.viscosity_factor * min(local_domain.dx, local_domain.dy)
    else:
        src = np.zeros((ny, nx), dtype=np.float64)

    while t < t_end:
        dt = min(dt0, t_end - t)
        eta_bc = interpolate_boundary_state(driver, t)
        impose_boundary_state(q, eta_bc)

        if use_numba:
            maxabseta, finite_ok = _advance_step_numba(
                eta, hu, hv, local_domain.z, zero_uplift, 0.0, 0.0, dt, t, local_domain.dx, local_domain.dy,
                cfg.solver.gravity, cfg.solver.min_depth_m, cfg.solver.damping, nu, land_mask,
                k10, k11, k12, k20, k21, k22, eta1, hu1, hv1, h_total, wet,
                max_eta, max_depth, max_speed, arrival_time, wet_mask,
            )
            eta_bc_next = interpolate_boundary_state(driver, t + dt)
            impose_boundary_state(q, eta_bc_next)
            if (not finite_ok) or maxabseta > cfg.solver.eta_guard_m:
                break
        else:
            k1 = compute_rhs(q, local_domain, src, cfg)
            q1 = q + dt * k1
            eta_bc_next = interpolate_boundary_state(driver, t + dt)
            impose_boundary_state(q1, eta_bc_next)
            k2 = compute_rhs(q1, local_domain, src, cfg)
            q = q + 0.5 * dt * (k1 + k2)
            impose_boundary_state(q, eta_bc_next)
            eta, hu, hv = q[0], q[1], q[2]

            h_total = np.maximum(eta - local_domain.z, 0.0)
            wet = h_total > cfg.solver.min_depth_m
            q[1, ~wet] = 0.0
            q[2, ~wet] = 0.0

            u = np.divide(q[1], h_total, out=np.zeros_like(q[1]), where=wet)
            v = np.divide(q[2], h_total, out=np.zeros_like(q[2]), where=wet)
            speed = np.sqrt(u * u + v * v)
            if not np.all(np.isfinite(eta)):
                break
            max_eta = np.maximum(max_eta, eta.astype(np.float32))
            inun = np.where(land_mask, h_total, 0.0)
            max_depth = np.maximum(max_depth, inun.astype(np.float32))
            max_speed = np.maximum(max_speed, speed.astype(np.float32))
            newly_wet = land_mask & wet & np.isnan(arrival_time)
            arrival_time[newly_wet] = np.float32(t)
            wet_mask |= wet.astype(np.uint8)
        t += dt

    positive = max_depth[(max_depth > 0.0) & land_mask]
    dapa_p95 = float(np.percentile(positive, 95)) if positive.size else 0.0
    return LocalNestedOutput(
        max_eta=max_eta,
        max_depth=max_depth,
        max_speed=max_speed,
        arrival_time=arrival_time,
        wet_mask=wet_mask,
        dapa_p95_inun=dapa_p95,
        meta={"dt0_sec": dt0, "t_end_sec": t_end, "nested_boundary_forcing": True, "kernel_backend": "numba" if use_numba else "numpy"},
    )


# ---------- reporting / catalog ----------

def build_viability_report(cfg: ProjectConfig, domain: PreparedDomain, p95_map: np.ndarray, dapa_p95_df: pd.DataFrame,
                           local_p95_map: Optional[np.ndarray] = None, local_domain: Optional[PreparedDomain] = None) -> Dict[str, Any]:
    if local_p95_map is not None and local_domain is not None:
        dapa_map = local_p95_map
        dapa_land = local_domain.grid.z > 0.0
        dx = local_domain.dx
        dy = local_domain.dy
    else:
        dapa_map = p95_map[domain.dapa_ys, domain.dapa_xs]
        dapa_land = domain.grid.z[domain.dapa_ys, domain.dapa_xs] > 0.0
        dx = domain.dx
        dy = domain.dy
    positive = dapa_map[(dapa_map > 0.0) & dapa_land]
    return {
        "dapa_only_mode": cfg.dapa_only_mode,
        "downsample_factor": cfg.execution.downsample_factor,
        "dx_m": float(dx),
        "dy_m": float(dy),
        "t_end_sec": float(cfg.solver.t_end_sec),
        "solver_has_moving_shoreline": True,
        "solver_has_local_nested_grid": bool(local_domain is not None),
        "osm_required_for_inundation_depth": False,
        "osm_required_for_buildings_and_roads_overlays": True,
        "regional_hazard_viable": bool(cfg.solver.t_end_sec >= 1000.0 and domain.dx <= 1000.0),
        "dapa_inundation_viable": bool(cfg.solver.t_end_sec >= 1000.0 and dx <= 60.0 and positive.size > 0),
        "paper_like_risk_overlay_viable_without_osm": False,
        "dapa_nonzero_inundation_cells": int(np.count_nonzero((dapa_map > 0.0) & dapa_land)),
        "dapa_p95_positive_depth_m": float(np.percentile(positive, 95)) if positive.size else 0.0,
        "mean_samplewise_dapa_p95_m": float(dapa_p95_df["dapa_p95_inundation_m"].mean()) if len(dapa_p95_df) else 0.0,
    }


def scenario_summary_row(mw: float, sample: SampleOutput, runtime_sec: float, domain: PreparedDomain, use_local: bool) -> Dict[str, Any]:
    if use_local and sample.local_nested is not None and domain.nested_local_domain is not None:
        depth = sample.local_nested.max_depth
        eta = sample.local_nested.max_eta
        speed = sample.local_nested.max_speed
        arr = sample.local_nested.arrival_time
        land = domain.nested_local_domain.grid.z > 0.0
        dx = domain.nested_local_domain.dx
        dy = domain.nested_local_domain.dy
    else:
        depth = sample.max_depth[domain.dapa_ys, domain.dapa_xs]
        eta = sample.max_eta[domain.dapa_ys, domain.dapa_xs]
        speed = sample.max_speed[domain.dapa_ys, domain.dapa_xs]
        arr = sample.arrival_time[domain.dapa_ys, domain.dapa_xs]
        land = domain.grid.z[domain.dapa_ys, domain.dapa_xs] > 0.0
        dx = domain.dx
        dy = domain.dy

    positive = depth[(depth > 0.0) & land]
    first_arrival = float(np.nanmin(arr[land])) if np.any(np.isfinite(arr[land])) else math.nan
    return {
        "mw": mw,
        "slip_id": sample.slip_id,
        "forcing_mode": sample.forcing_mode,
        "runtime_sec": runtime_sec,
        "max_dapa_depth_m": float(np.max(depth)) if depth.size else 0.0,
        "p95_dapa_depth_m": float(np.percentile(positive, 95)) if positive.size else 0.0,
        "max_dapa_eta_m": float(np.max(np.abs(eta))) if eta.size else 0.0,
        "max_dapa_speed_ms": float(np.max(speed)) if speed.size else 0.0,
        "inundated_area_m2": float(np.count_nonzero((depth > 0.0) & land) * dx * dy),
        "first_arrival_sec": first_arrival,
        **{f"gauge_max_abs_{k}_m": float(v) for k, v in sample.gauge_max_abs.items()},
    }


def write_validation_report(out_dir: Path, event: HistoricalEventConfig, sample: SampleOutput, metrics_rows: List[Dict[str, Any]]) -> None:
    pd.DataFrame(metrics_rows).to_csv(out_dir / f"{event.event_id}_validation_metrics.csv", index=False)
    save_json(out_dir / f"{event.event_id}_validation_summary.json", {
        "event_id": event.event_id,
        "mw": event.mw,
        "nrmse_target": 0.3,
        "arrival_time_target_sec": 120.0,
        "gauges_compared": [r["gauge_name"] for r in metrics_rows],
    })
    save_gauge_series(out_dir, sample)


# ---------- validation ----------

def load_historical_event_config(path: Path) -> List[HistoricalEventConfig]:
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    events: List[HistoricalEventConfig] = []
    for ev in payload.get("events", []):
        gauges = [GaugePoint(**g) for g in ev.get("gauges", [])]
        observations = [
            ValidationObservation(
                gauge_name=o["gauge_name"],
                csv_path=Path(o["csv_path"]),
                time_col=o.get("time_col", "time_sec"),
                eta_col=o.get("eta_col", "eta_m"),
            )
            for o in ev.get("observations", [])
        ]
        events.append(HistoricalEventConfig(
            event_id=ev["event_id"],
            mw=float(ev["mw"]),
            output_subdir=ev.get("output_subdir", ev["event_id"]),
            source_overrides=ev.get("source_overrides", {}),
            gauges=gauges,
            observations=observations,
        ))
    return events


def compute_validation_metrics(sim_time: np.ndarray, sim_eta: np.ndarray, obs_time: np.ndarray, obs_eta: np.ndarray) -> Dict[str, float]:
    if len(sim_time) < 2 or len(obs_time) < 2:
        return {"nrmse": math.nan, "arrival_time_error_sec": math.nan}
    interp = np.interp(obs_time, sim_time, sim_eta)
    denom = float(np.linalg.norm(obs_eta))
    nrmse = float(np.linalg.norm(interp - obs_eta) / denom) if denom > 0 else math.nan
    sim_thr_idx = np.where(np.abs(interp) > 0.01)[0]
    obs_thr_idx = np.where(np.abs(obs_eta) > 0.01)[0]
    if sim_thr_idx.size and obs_thr_idx.size:
        arr_err = float(obs_time[sim_thr_idx[0]] - obs_time[obs_thr_idx[0]])
    else:
        arr_err = math.nan
    return {"nrmse": nrmse, "arrival_time_error_sec": arr_err}


def plot_validation_waveform(path: Path, gauge_name: str, sim_time: np.ndarray, sim_eta: np.ndarray,
                             obs_time: Optional[np.ndarray] = None, obs_eta: Optional[np.ndarray] = None) -> None:
    if plt is None:
        return
    ensure_dir(path.parent)
    fig, ax = plt.subplots(figsize=(8, 4), dpi=180)
    ax.plot(sim_time, sim_eta, label="Simulated")
    if obs_time is not None and obs_eta is not None:
        ax.plot(obs_time, obs_eta, label="Observed", alpha=0.8)
    ax.set_title(f"Waveform comparison | {gauge_name}")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Water level (m)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def apply_source_overrides(cfg: ProjectConfig, overrides: Dict[str, float]) -> None:
    for k, v in overrides.items():
        if hasattr(cfg.source, k):
            setattr(cfg.source, k, v)


# ---------- orchestration modes ----------

def run_single_sample_bundle(cfg: ProjectConfig, domain: PreparedDomain, mw: float, slip: np.ndarray, slip_id: str,
                             out_dir: Path, forcing_mode: str, write_maps: bool = True) -> Tuple[SampleOutput, Dict[str, Any]]:
    start = time.perf_counter()
    sample = run_tsunami_simulation(cfg, domain, slip, slip_id, sample_index=0, forcing_mode=forcing_mode)
    if domain.nested_local_domain is not None and sample.nested_time_sec is not None and sample.nested_eta_box is not None and sample.nested_time_sec.size >= 2:
        driver = build_local_boundary_driver(domain, domain.nested_local_domain, sample.nested_time_sec, sample.nested_eta_box)
        sample.local_nested = run_local_nested_simulation(cfg, domain.nested_local_domain, driver)
    runtime_sec = time.perf_counter() - start

    ensure_dir(out_dir)
    gauge_csv = save_gauge_series(out_dir, sample)
    save_json(out_dir / f"scenario_meta_{slip_id}_{forcing_mode}.json", sample.meta)

    exports: Dict[str, Any] = {"gauge_csv": str(gauge_csv)}
    if write_maps:
        exports["regional"] = export_map_bundle(
            out_dir,
            f"{slip_id}_{forcing_mode}_regional",
            domain.grid,
            sample.max_depth,
            sample.max_eta,
            sample.max_speed,
            sample.arrival_time,
            sample.wet_mask,
            write_geotiff=cfg.execution.write_geotiff,
            write_npz=cfg.execution.write_npz,
        )
        dapa_grid = Grid2D(
            lon=domain.grid.lon[domain.dapa_xs],
            lat=domain.grid.lat[domain.dapa_ys],
            z=domain.grid.z[domain.dapa_ys, domain.dapa_xs],
        )
        exports["dapa_corridor"] = export_map_bundle(
            out_dir,
            f"{slip_id}_{forcing_mode}_dapa_corridor",
            dapa_grid,
            sample.max_depth[domain.dapa_ys, domain.dapa_xs],
            sample.max_eta[domain.dapa_ys, domain.dapa_xs],
            sample.max_speed[domain.dapa_ys, domain.dapa_xs],
            sample.arrival_time[domain.dapa_ys, domain.dapa_xs],
            sample.wet_mask[domain.dapa_ys, domain.dapa_xs],
            write_geotiff=cfg.execution.write_geotiff,
            write_npz=cfg.execution.write_npz,
        )
        if cfg.execution.write_png:
            save_map(domain.grid, sample.max_depth, f"Mw{mw:.1f} {forcing_mode} | max inundation depth", out_dir / f"{slip_id}_{forcing_mode}_regional_max_depth.png", paper_style=cfg.execution.paper_style_maps)
            save_map(dapa_grid, sample.max_depth[domain.dapa_ys, domain.dapa_xs], f"Mw{mw:.1f} {forcing_mode} | Dapa corridor max inundation depth", out_dir / f"{slip_id}_{forcing_mode}_dapa_corridor_max_depth.png", paper_style=True)
        if sample.local_nested is not None and domain.nested_local_domain is not None:
            exports["dapa_local"] = export_map_bundle(
                out_dir,
                f"{slip_id}_{forcing_mode}_dapa_local",
                domain.nested_local_domain.grid,
                sample.local_nested.max_depth,
                sample.local_nested.max_eta,
                sample.local_nested.max_speed,
                sample.local_nested.arrival_time,
                sample.local_nested.wet_mask,
                write_geotiff=cfg.execution.write_geotiff,
                write_npz=cfg.execution.write_npz,
            )
            if cfg.execution.write_png:
                save_map(domain.nested_local_domain.grid, sample.local_nested.max_depth,
                         f"Mw{mw:.1f} {forcing_mode} | Dapa local max inundation depth",
                         out_dir / f"{slip_id}_{forcing_mode}_dapa_local_max_depth.png", paper_style=True)
    exports["runtime_sec"] = runtime_sec
    return sample, exports


def run_proof_of_concept(cfg: ProjectConfig, domain: PreparedDomain) -> None:
    out_dir = cfg.results_dir / "proof_of_concept"
    ensure_dir(out_dir)
    apply_magnitude_scenario(cfg, 8.5)
    slip, slip_id = gather_slips(cfg, 8.5, 1)[0]
    sample, exports = run_single_sample_bundle(cfg, domain, 8.5, slip, slip_id, out_dir / "mw8p5_dynamic", forcing_mode="dynamic")
    row = scenario_summary_row(8.5, sample, exports["runtime_sec"], domain, use_local=True)
    pd.DataFrame([row]).to_csv(out_dir / "proof_of_concept_catalog.csv", index=False)
    save_json(out_dir / "proof_of_concept_summary.json", {
        "mode": "proof-of-concept",
        "deliverable": "One Mw 8.5 synthetic case on the revised dynamic-forcing workflow",
        "runtime_sec": exports["runtime_sec"],
    })


def run_compare_static_dynamic(cfg: ProjectConfig, domain: PreparedDomain, mw: float, slip_index: int = 0) -> None:
    out_dir = cfg.results_dir / "static_vs_dynamic" / f"Mw{mw:.1f}"
    ensure_dir(out_dir)
    slips = gather_slips(cfg, mw, slip_index + 1)
    slip, slip_id = slips[slip_index]
    dyn, dyn_exports = run_single_sample_bundle(cfg, domain, mw, slip, slip_id, out_dir / "dynamic", forcing_mode="dynamic")
    sta, sta_exports = run_single_sample_bundle(cfg, domain, mw, slip, slip_id, out_dir / "static", forcing_mode="static")

    if domain.nested_local_domain is not None and dyn.local_nested is not None and sta.local_nested is not None:
        diff_depth = dyn.local_nested.max_depth - sta.local_nested.max_depth
        diff_grid = domain.nested_local_domain.grid
        prefix = f"{slip_id}_dynamic_minus_static_dapa_local"
    else:
        diff_depth = dyn.max_depth[domain.dapa_ys, domain.dapa_xs] - sta.max_depth[domain.dapa_ys, domain.dapa_xs]
        diff_grid = Grid2D(lon=domain.grid.lon[domain.dapa_xs], lat=domain.grid.lat[domain.dapa_ys], z=domain.grid.z[domain.dapa_ys, domain.dapa_xs])
        prefix = f"{slip_id}_dynamic_minus_static_dapa_corridor"

    if cfg.execution.write_geotiff:
        export_geotiff(out_dir / f"{prefix}_depth_diff.tif", diff_grid, diff_depth, dtype="float32", nodata=np.nan)
    if cfg.execution.write_png:
        save_map(diff_grid, diff_depth, f"Mw{mw:.1f} dynamic - static | max depth difference", out_dir / f"{prefix}_depth_diff.png")

    rows = [
        scenario_summary_row(mw, dyn, dyn_exports["runtime_sec"], domain, use_local=True),
        scenario_summary_row(mw, sta, sta_exports["runtime_sec"], domain, use_local=True),
    ]
    pd.DataFrame(rows).to_csv(out_dir / "static_dynamic_catalog.csv", index=False)
    save_json(out_dir / "static_dynamic_summary.json", {
        "mw": mw,
        "slip_id": slip_id,
        "dynamic_runtime_sec": dyn_exports["runtime_sec"],
        "static_runtime_sec": sta_exports["runtime_sec"],
        "max_depth_diff_m": float(np.nanmax(diff_depth)) if diff_depth.size else 0.0,
        "min_depth_diff_m": float(np.nanmin(diff_depth)) if diff_depth.size else 0.0,
    })


def run_dapa_ensemble(cfg: ProjectConfig, domain: PreparedDomain, mw: float, n_samples: int) -> Dict[str, Any]:
    apply_magnitude_scenario(cfg, mw)
    mag_dir = cfg.results_dir / "dapa_ensemble" / f"Mw{mw:.1f}"
    ensure_dir(mag_dir)
    samples = gather_slips(cfg, mw, n_samples)

    maybe_print(f"\n===== Mw{mw:.1f} Dapa-focused ensemble =====")
    pbar = make_progress(total=len(samples), desc=f"Mw{mw:.1f}")

    gauge_rows: List[Dict[str, Any]] = []
    p95_rows: List[Dict[str, Any]] = []
    catalog_rows: List[Dict[str, Any]] = []
    max_depth_list: List[np.ndarray] = []
    local_max_depth_list: List[np.ndarray] = []

    for i, (slip, slip_name) in enumerate(samples):
        start = time.perf_counter()
        sample = run_tsunami_simulation(cfg, domain, slip, slip_name, i, forcing_mode="dynamic", pbar=pbar)
        if domain.nested_local_domain is not None and sample.nested_time_sec is not None and sample.nested_eta_box is not None and sample.nested_time_sec.size >= 2:
            driver = build_local_boundary_driver(domain, domain.nested_local_domain, sample.nested_time_sec, sample.nested_eta_box)
            sample.local_nested = run_local_nested_simulation(cfg, domain.nested_local_domain, driver, pbar)
        runtime_sec = time.perf_counter() - start

        row = {"sample_index": i, "slip_id": slip_name, "forcing_mode": sample.forcing_mode, **sample.gauge_max_abs}
        gauge_rows.append(row)
        p95_rows.append({
            "sample_index": i,
            "slip_id": slip_name,
            "dapa_p95_inundation_m": sample.local_nested.dapa_p95_inun if sample.local_nested is not None else sample.dapa_p95_inun,
        })
        catalog_rows.append(scenario_summary_row(mw, sample, runtime_sec, domain, use_local=True))
        max_depth_list.append(sample.max_depth.astype(np.float32))
        if sample.local_nested is not None:
            local_max_depth_list.append(sample.local_nested.max_depth.astype(np.float32))

        sample_dir = mag_dir / f"sample_{i:04d}_{slip_name}"
        ensure_dir(sample_dir)
        if i == 0 or cfg.execution.save_per_sample_maps:
            save_gauge_series(sample_dir, sample)
            save_json(sample_dir / "sample_meta.json", sample.meta)
            export_map_bundle(sample_dir, f"{slip_name}_regional", domain.grid, sample.max_depth, sample.max_eta, sample.max_speed, sample.arrival_time, sample.wet_mask,
                              write_geotiff=cfg.execution.write_geotiff, write_npz=cfg.execution.write_npz)
            dapa_grid = Grid2D(lon=domain.grid.lon[domain.dapa_xs], lat=domain.grid.lat[domain.dapa_ys], z=domain.grid.z[domain.dapa_ys, domain.dapa_xs])
            export_map_bundle(sample_dir, f"{slip_name}_dapa_corridor", dapa_grid,
                              sample.max_depth[domain.dapa_ys, domain.dapa_xs],
                              sample.max_eta[domain.dapa_ys, domain.dapa_xs],
                              sample.max_speed[domain.dapa_ys, domain.dapa_xs],
                              sample.arrival_time[domain.dapa_ys, domain.dapa_xs],
                              sample.wet_mask[domain.dapa_ys, domain.dapa_xs],
                              write_geotiff=cfg.execution.write_geotiff, write_npz=cfg.execution.write_npz)
            if sample.local_nested is not None and domain.nested_local_domain is not None:
                export_map_bundle(sample_dir, f"{slip_name}_dapa_local", domain.nested_local_domain.grid,
                                  sample.local_nested.max_depth, sample.local_nested.max_eta, sample.local_nested.max_speed,
                                  sample.local_nested.arrival_time, sample.local_nested.wet_mask,
                                  write_geotiff=cfg.execution.write_geotiff, write_npz=cfg.execution.write_npz)

        pbar.update(1)
        if (i + 1) % cfg.execution.checkpoint_every == 0 or (i + 1) == len(samples):
            maybe_print(f"[CHECKPOINT] Mw{mw:.1f}: completed {i + 1}/{len(samples)} samples", pbar)
            save_checkpoint(mag_dir, mw, i, gauge_rows, p95_rows)
    pbar.close()

    gauge_df = pd.DataFrame(gauge_rows)
    p95_df = pd.DataFrame(p95_rows)
    catalog_df = pd.DataFrame(catalog_rows)
    gauge_df.to_csv(mag_dir / f"Mw{mw:.1f}_gauge_maxima.csv", index=False)
    p95_df.to_csv(mag_dir / f"Mw{mw:.1f}_dapa_p95_inundation.csv", index=False)
    catalog_df.to_csv(mag_dir / f"Mw{mw:.1f}_scenario_catalog.csv", index=False)

    stack = np.stack(max_depth_list, axis=0)
    np.save(mag_dir / f"Mw{mw:.1f}_regional_max_depth_stack.npy", stack)
    regional_p95 = np.percentile(stack, 95, axis=0).astype(np.float32)
    export_map_bundle(mag_dir, f"Mw{mw:.1f}_regional_p95", domain.grid, regional_p95, np.zeros_like(regional_p95), np.zeros_like(regional_p95), np.full_like(regional_p95, np.nan), (regional_p95 > 0).astype(np.uint8),
                      write_geotiff=cfg.execution.write_geotiff, write_npz=cfg.execution.write_npz)
    if cfg.execution.write_png:
        save_map(domain.grid, regional_p95, f"Mw{mw:.1f} | regional p95 inundation depth", mag_dir / f"Mw{mw:.1f}_regional_p95_inundation.png", paper_style=cfg.execution.paper_style_maps)

    dapa_grid = Grid2D(lon=domain.grid.lon[domain.dapa_xs], lat=domain.grid.lat[domain.dapa_ys], z=domain.grid.z[domain.dapa_ys, domain.dapa_xs])
    dapa_corridor_p95 = regional_p95[domain.dapa_ys, domain.dapa_xs]
    export_map_bundle(mag_dir, f"Mw{mw:.1f}_dapa_corridor_p95", dapa_grid, dapa_corridor_p95, np.zeros_like(dapa_corridor_p95), np.zeros_like(dapa_corridor_p95), np.full_like(dapa_corridor_p95, np.nan), (dapa_corridor_p95 > 0).astype(np.uint8),
                      write_geotiff=cfg.execution.write_geotiff, write_npz=cfg.execution.write_npz)
    if cfg.execution.write_png:
        save_map(dapa_grid, dapa_corridor_p95, f"Mw{mw:.1f} | Dapa corridor p95 inundation depth", mag_dir / f"Mw{mw:.1f}_dapa_corridor_p95_inundation.png", paper_style=True)

    local_p95_map = None
    if local_max_depth_list and domain.nested_local_domain is not None:
        local_stack = np.stack(local_max_depth_list, axis=0)
        np.save(mag_dir / f"Mw{mw:.1f}_dapa_local_max_depth_stack.npy", local_stack)
        local_p95_map = np.percentile(local_stack, 95, axis=0).astype(np.float32)
        export_map_bundle(mag_dir, f"Mw{mw:.1f}_dapa_local_p95", domain.nested_local_domain.grid, local_p95_map, np.zeros_like(local_p95_map), np.zeros_like(local_p95_map), np.full_like(local_p95_map, np.nan), (local_p95_map > 0).astype(np.uint8),
                          write_geotiff=cfg.execution.write_geotiff, write_npz=cfg.execution.write_npz)
        if cfg.execution.write_png:
            save_map(domain.nested_local_domain.grid, local_p95_map, f"Mw{mw:.1f} | Dapa local p95 inundation depth", mag_dir / f"Mw{mw:.1f}_dapa_local_p95_inundation.png", paper_style=True)

    conv_rows, diff_rows = [], []
    prev = {}
    for n in convergence_schedule(stack.shape[0]):
        cur_stack = stack[:n]
        vals = cur_stack[cur_stack > 0]
        row = {"mw": mw, "n_samples": n, "regional_p95_inundation_m": float(np.percentile(vals, 95)) if vals.size else 0.0}
        for g in [gp.name for gp in cfg.gauges]:
            row[f"mean_max_wave_{g}_m"] = float(gauge_df[g].iloc[:n].mean()) if g in gauge_df else math.nan
        conv_rows.append(row)
        for thr in cfg.probability.thresholds_m:
            cur = exceedance_map(cur_stack, thr)
            if thr in prev:
                diff = np.abs(cur - prev[thr])
                diff_rows.append({"mw": mw, "threshold_m": thr, "to_n": n, "mean_abs_diff": float(np.mean(diff)), "max_abs_diff": float(np.max(diff))})
            prev[thr] = cur
    pd.DataFrame(conv_rows).to_csv(mag_dir / f"Mw{mw:.1f}_convergence_summary.csv", index=False)
    pd.DataFrame(diff_rows).to_csv(mag_dir / f"Mw{mw:.1f}_ape_successive_differences.csv", index=False)

    if cfg.execution.emit_phase4_preview:
        for thr in cfg.probability.thresholds_m:
            ape = exceedance_map(stack, thr)
            export_map_bundle(mag_dir, f"Mw{mw:.1f}_conditional_exceedance_{thr:.1f}m", domain.grid, ape, np.zeros_like(ape), np.zeros_like(ape), np.full_like(ape, np.nan), (ape > 0).astype(np.uint8),
                              write_geotiff=cfg.execution.write_geotiff, write_npz=cfg.execution.write_npz)
            if cfg.execution.write_png:
                save_map(domain.grid, ape, f"Mw{mw:.1f} | P(depth>{thr:.1f}m)", mag_dir / f"Mw{mw:.1f}_conditional_exceedance_{thr:.1f}m.png")

    viability = build_viability_report(cfg, domain, regional_p95, p95_df, local_p95_map=local_p95_map, local_domain=domain.nested_local_domain)
    save_json(mag_dir / f"Mw{mw:.1f}_viability_report.json", viability)
    save_json(mag_dir / f"Mw{mw:.1f}_summary.json", {
        "mw": mw,
        "n_samples": n_samples,
        "mean_max_wave_by_gauge_m": {g.name: float(gauge_df[g.name].mean()) for g in cfg.gauges if g.name in gauge_df},
        "95th_percentile_inundation_depth_in_dapa_m": float(np.percentile(p95_df["dapa_p95_inundation_m"].values, 95)) if len(p95_df) else 0.0,
        "viability": viability,
    })
    return {"gauge_df": gauge_df, "p95_df": p95_df, "catalog_df": catalog_df, "regional_p95": regional_p95, "local_p95": local_p95_map}


def run_validation_mode(cfg: ProjectConfig, domain: PreparedDomain, event_config_path: Path) -> None:
    events = load_historical_event_config(event_config_path)
    root = cfg.results_dir / "validation"
    ensure_dir(root)
    summary_rows: List[Dict[str, Any]] = []
    for event in events:
        event_dir = root / event.output_subdir
        ensure_dir(event_dir)
        cfg_event = ProjectConfig(project_root=cfg.project_root)
        cfg_event.domain = cfg.domain
        cfg_event.slip = cfg.slip
        cfg_event.probability = cfg.probability
        cfg_event.execution = cfg.execution
        cfg_event.results_dirname = cfg.results_dirname
        cfg_event.solver = cfg.solver
        cfg_event.gauges = tuple(event.gauges) if event.gauges else cfg.gauges
        apply_magnitude_scenario(cfg_event, event.mw)
        apply_source_overrides(cfg_event, event.source_overrides)

        slips = gather_slips(cfg_event, event.mw, 1)
        slip, slip_id = slips[0]
        sample, exports = run_single_sample_bundle(cfg_event, domain, event.mw, slip, slip_id, event_dir, forcing_mode="dynamic")
        metrics_rows: List[Dict[str, Any]] = []
        sim_df = pd.DataFrame({"time_sec": sample.gauge_time_sec, **sample.gauge_eta})
        for obs in event.observations:
            obs_path = obs.csv_path if obs.csv_path.is_absolute() else cfg.project_root / obs.csv_path
            if not obs_path.exists():
                metrics_rows.append({
                    "event_id": event.event_id,
                    "gauge_name": obs.gauge_name,
                    "status": "missing_observation_file",
                    "csv_path": str(obs_path),
                })
                continue
            obs_df = pd.read_csv(obs_path)
            if obs.gauge_name not in sim_df.columns:
                metrics_rows.append({
                    "event_id": event.event_id,
                    "gauge_name": obs.gauge_name,
                    "status": "gauge_not_in_simulation",
                })
                continue
            metrics = compute_validation_metrics(
                sim_df["time_sec"].to_numpy(dtype=float),
                sim_df[obs.gauge_name].to_numpy(dtype=float),
                obs_df[obs.time_col].to_numpy(dtype=float),
                obs_df[obs.eta_col].to_numpy(dtype=float),
            )
            plot_validation_waveform(
                event_dir / f"waveform_{obs.gauge_name}.png",
                obs.gauge_name,
                sim_df["time_sec"].to_numpy(dtype=float),
                sim_df[obs.gauge_name].to_numpy(dtype=float),
                obs_df[obs.time_col].to_numpy(dtype=float),
                obs_df[obs.eta_col].to_numpy(dtype=float),
            )
            metrics_rows.append({
                "event_id": event.event_id,
                "gauge_name": obs.gauge_name,
                "status": "compared",
                **metrics,
            })
        write_validation_report(event_dir, event, sample, metrics_rows)
        summary_rows.extend(metrics_rows)
        save_json(event_dir / "event_run_summary.json", {
            "event_id": event.event_id,
            "mw": event.mw,
            "runtime_sec": exports["runtime_sec"],
            "source_overrides": event.source_overrides,
        })
    pd.DataFrame(summary_rows).to_csv(root / "validation_overview.csv", index=False)


def run_phase3_thesis(cfg: ProjectConfig, domain: PreparedDomain, counts: Dict[float, int]) -> None:
    run_proof_of_concept(cfg, domain)
    run_compare_static_dynamic(cfg, domain, mw=8.5, slip_index=0)
    ensemble_rows = []
    for mw in (8.0, 8.5, 9.0):
        res = run_dapa_ensemble(cfg, domain, mw, counts.get(mw, 1))
        if "catalog_df" in res:
            ensemble_rows.append(res["catalog_df"])

    if ensemble_rows:
        overall_dir = cfg.results_dir / "catalogs"
        ensure_dir(overall_dir)
        full_catalog = pd.concat(ensemble_rows, ignore_index=True)
        full_catalog.to_csv(overall_dir / "phase3_scenario_catalog.csv", index=False)

    if cfg.execution.emit_phase4_preview:
        overall_dir = cfg.results_dir / "phase4_preview"
        ensure_dir(overall_dir)
        weights = magnitude_weights(cfg.probability.beta, [8.0, 8.5, 9.0])
        save_json(overall_dir / "preview_weights.json", {str(k): v for k, v in weights.items()})


# ---------- CLI ----------

def parse_counts(text: str) -> Dict[float, int]:
    out = {8.0: 16, 8.5: 16, 9.0: 16}
    if not text:
        return out
    for item in text.split(","):
        k, v = item.split("=")
        out[float(k.strip())] = int(v.strip())
    return out


def build_config(args: argparse.Namespace) -> ProjectConfig:
    cfg = ProjectConfig(project_root=args.project_root)
    cfg.solver.t_end_sec = args.t_end_sec
    cfg.solver.dt_max = args.dt_max
    cfg.solver.diag_stride = args.diag_stride
    cfg.execution.checkpoint_every = args.checkpoint_every
    cfg.execution.downsample_factor = args.downsample_factor
    cfg.execution.save_per_sample_maps = args.save_per_sample_maps
    cfg.execution.paper_style_maps = not args.no_paper_style_maps
    cfg.execution.local_dapa_render = not args.no_local_dapa_render
    cfg.execution.local_dapa_nested = not args.no_local_dapa_nested
    cfg.execution.nested_save_stride = args.nested_save_stride
    cfg.execution.write_geotiff = not args.no_geotiff
    cfg.execution.write_npz = not args.no_npz
    cfg.execution.write_png = not args.no_png
    cfg.execution.emit_phase4_preview = args.emit_phase4_preview
    cfg.execution.kernel_backend = args.kernel_backend
    cfg.execution.numba_threads = args.numba_threads
    cfg.results_dirname = args.results_dirname
    cfg.dapa_only_mode = bool(args.dapa_only)
    if NUMBA_AVAILABLE and cfg.execution.numba_threads and cfg.execution.numba_threads > 0:
        try:
            set_num_threads(int(cfg.execution.numba_threads))
        except Exception:
            pass
    if args.local_dem is not None:
        cfg.local_dem_path_override = args.local_dem
    if args.slip_dir is not None:
        cfg.slip_dir_override = args.slip_dir

    if args.dapa_only:
        cfg.domain.regional_lon_min = cfg.domain.corridor_lon_min
        cfg.domain.regional_lon_max = cfg.domain.corridor_lon_max
        cfg.domain.regional_lat_min = cfg.domain.corridor_lat_min
        cfg.domain.regional_lat_max = cfg.domain.corridor_lat_max
        cfg.gauges = (GaugePoint("Dapa_virtual", DAPA_MARKER[0], DAPA_MARKER[1]),)
    return cfg


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Revised Phase 3 runner aligned to the thesis memo: dynamic forcing, targeted Dapa refinement, and map-ready exports.")
    p.add_argument("--project-root", type=Path, default=Path("."))
    p.add_argument("--mode", choices=["proof-of-concept", "dapa-ensemble", "compare-static-dynamic", "validation", "thesis"], default="thesis")
    p.add_argument("--counts", type=str, default="8.0=16,8.5=16,9.0=16")
    p.add_argument("--mw", type=float, default=8.5, help="Magnitude for dapa-ensemble or compare-static-dynamic modes")
    p.add_argument("--slip-index", type=int, default=0, help="Slip index for compare-static-dynamic mode")
    p.add_argument("--dapa-only", action="store_true", help="Run a narrower source-to-Dapa corridor and keep only Dapa-focused outputs")
    p.add_argument("--t-end-sec", type=float, default=10800.0, help="Use 10800 for the 3 h synthetic runs in the reference study")
    p.add_argument("--dt-max", type=float, default=0.15, help="Use 0.15 s to align with the reference study")
    p.add_argument("--diag-stride", type=int, default=250)
    p.add_argument("--checkpoint-every", type=int, default=25)
    p.add_argument("--downsample-factor", type=int, default=1)
    p.add_argument("--save-per-sample-maps", action="store_true")
    p.add_argument("--no-paper-style-maps", action="store_true")
    p.add_argument("--no-local-dapa-render", action="store_true")
    p.add_argument("--no-local-dapa-nested", action="store_true")
    p.add_argument("--nested-save-stride", type=int, default=10)
    p.add_argument("--local-dem", type=Path, default=None, help="Path to merged fine Dapa topo-bathy DEM (GeoTIFF)")
    p.add_argument("--validation-config", type=Path, default=None, help="JSON file describing the 2012/2023 validation events and observed gauge CSVs")
    p.add_argument("--slip-dir", type=Path, default=None, help="Override slip-sample directory, e.g. output/RQMC_samples_revised")
    p.add_argument("--kernel-backend", choices=["auto", "numba", "numpy"], default="auto", help="Use numba-compiled kernels when available for large speedups")
    p.add_argument("--numba-threads", type=int, default=0, help="Optional numba thread count; 0 keeps the numba default")
    p.add_argument("--emit-phase4-preview", action="store_true", help="Optional only. Writes conditional exceedance preview maps but keeps Phase 4 separate.")
    p.add_argument("--no-geotiff", action="store_true")
    p.add_argument("--no-npz", action="store_true")
    p.add_argument("--no-png", action="store_true")
    p.add_argument("--results-dirname", type=str, default="phase3_revised_outputs")
    p.add_argument("--export-config", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = build_config(args)
    ensure_dir(cfg.results_dir)
    if args.export_config:
        save_json(cfg.results_dir / "resolved_config.json", asdict(cfg))
    if args.dapa_only:
        print("[NOTE] Dapa-only uses a narrower source-to-Dapa corridor for faster higher-resolution runs.", flush=True)
    print("[NOTE] Phase 3 is now framed as high-fidelity tsunami propagation with dynamic seafloor forcing.", flush=True)
    print("[NOTE] GeoTIFF, NPZ, and PNG map-ready bundles are emitted for Dapa corridor and local Dapa grids when available.", flush=True)
    print("[NOTE] Static-vs-dynamic comparison is a dedicated mode, not mixed into the default ensemble.", flush=True)
    print("[NOTE] Phase 4-style probabilistic aggregation is disabled by default and only available as an optional preview.", flush=True)
    backend = cfg.execution.kernel_backend.lower()
    if backend == "auto":
        resolved = "numba" if NUMBA_AVAILABLE else "numpy"
    else:
        resolved = backend
    print(f"[NOTE] Solver kernel backend: {resolved} (first run may spend time compiling if numba is active).", flush=True)

    domain = prepare_domain(cfg)
    dt0 = stable_dt(domain, cfg)
    print(f"[INFO] Recommended timestep strategy: adaptive CFL with dt_max={cfg.solver.dt_max:.2f}s (initial stable dt≈{dt0:.4f}s)", flush=True)

    if args.mode == "proof-of-concept":
        run_proof_of_concept(cfg, domain)
    elif args.mode == "compare-static-dynamic":
        run_compare_static_dynamic(cfg, domain, mw=args.mw, slip_index=args.slip_index)
    elif args.mode == "dapa-ensemble":
        counts = parse_counts(args.counts)
        run_dapa_ensemble(cfg, domain, mw=args.mw, n_samples=counts.get(args.mw, 1))
    elif args.mode == "validation":
        if args.validation_config is None:
            raise SystemExit("--validation-config is required for validation mode.")
        run_validation_mode(cfg, domain, args.validation_config)
    elif args.mode == "thesis":
        run_phase3_thesis(cfg, domain, parse_counts(args.counts))
    else:
        raise SystemExit(f"Unsupported mode: {args.mode}")


if __name__ == "__main__":
    main()
