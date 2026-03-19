#!/usr/bin/env python3
"""
phase3_dapa_noosm_tools.py

Standalone utilities for Dapa-focused tsunami inundation post-processing,
viability checks, and paper-style map rendering without requiring OSM.

What this script does:
1) Builds a Dapa local DEM from GEBCO and optional SRTM.
2) Computes inundation depth on land from either:
   - a provided max-depth raster, or
   - a provided max-water-surface raster (eta_max) and the DEM.
3) Produces a paper-style map using blue inundation classes similar to the
   Dapa figures in "Tsunami hazards and risks from the PH Trench".
4) Writes a JSON viability report that tells you whether the result is only a
   smoke test or already usable for inundation interpretation.

Important limitation:
This is a post-processor and mapper. It does NOT replace your hydrodynamic
solver. It is meant to sit beside your existing runner and make its Dapa
outputs more trustworthy and easier to interpret.

Inputs accepted for tsunami results:
- NPZ file containing one of the following:
  A) lon, lat, depth_max
  B) lon, lat, eta_max
  C) x, y, depth_max   (aliases accepted if x/y are actually lon/lat)

Typical workflow:
  1. Run your existing solver longer (>= 3600 s, preferably 10800 s for paper-like comparisons).
  2. Export a local or regional max field to NPZ.
  3. Use this script to compute land inundation, map it, and inspect viability.

Example:
python phase3_dapa_noosm_tools.py \
    --gebco /path/to/gebco_2025_n25.0_s1.0_w115.0_e140.0.nc \
    --srtm /path/to/output_SRTMGL1.tif \
    --results /path/to/Mw85_dapa_maxfield.npz \
    --out-prefix /path/to/out/Mw85_dapa \
    --target-res-m 30 \
    --t-end-sec 3600 \
    --dx-m 30

Optional future overlays:
- Buildings and roads are NOT required for viable inundation depth maps.
- If you later acquire building/road vectors, you can add overlays in a
  separate plotting step.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional, Tuple, Dict, Any

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import from_bounds
from rasterio.warp import reproject
from scipy.interpolate import RegularGridInterpolator
import xarray as xr

# -----------------------------
# Dapa defaults
# -----------------------------
DAPA_LON = 126.0538
DAPA_LAT = 9.7592
DEFAULT_HALF_WIDTH_KM = 3.5
DEFAULT_HALF_HEIGHT_KM = 2.5


@dataclass
class ViabilityReport:
    has_srtm: bool
    target_res_m: float
    dx_m_reported: Optional[float]
    t_end_sec: Optional[float]
    inland_wet_cells: int
    inland_wet_area_m2: float
    inland_wet_area_ha: float
    max_land_inundation_m: float
    p95_land_inundation_m: float
    local_dem_min_m: float
    local_dem_max_m: float
    is_viable_for_smoke_test: bool
    is_viable_for_inundation_interpretation: bool
    messages: list[str]


def km_to_deg_lon(km: float, lat_deg: float) -> float:
    return km / (111.32 * max(math.cos(math.radians(lat_deg)), 1e-6))


def km_to_deg_lat(km: float) -> float:
    return km / 110.57


def make_bbox(
    center_lon: float = DAPA_LON,
    center_lat: float = DAPA_LAT,
    half_width_km: float = DEFAULT_HALF_WIDTH_KM,
    half_height_km: float = DEFAULT_HALF_HEIGHT_KM,
) -> Tuple[float, float, float, float]:
    dlon = km_to_deg_lon(half_width_km, center_lat)
    dlat = km_to_deg_lat(half_height_km)
    return center_lon - dlon, center_lon + dlon, center_lat - dlat, center_lat + dlat


def res_m_to_deg(res_m: float, lat_deg: float) -> Tuple[float, float]:
    dlat = res_m / 110570.0
    dlon = res_m / (111320.0 * max(math.cos(math.radians(lat_deg)), 1e-6))
    return dlon, dlat


# -----------------------------
# DEM building
# -----------------------------
def load_gebco_subset(
    gebco_nc: Path,
    bbox: Tuple[float, float, float, float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    minlon, maxlon, minlat, maxlat = bbox
    ds = xr.open_dataset(gebco_nc)
    if "elevation" not in ds:
        raise ValueError("GEBCO netCDF must contain an 'elevation' variable.")

    sub = ds.sel(lon=slice(minlon, maxlon), lat=slice(minlat, maxlat))
    if sub.elevation.size == 0:
        raise ValueError("Requested Dapa bbox is empty in the GEBCO file.")

    lon = np.asarray(sub.lon.values, dtype=float)
    lat = np.asarray(sub.lat.values, dtype=float)
    z = np.asarray(sub.elevation.values, dtype=float)
    return lon, lat, z


def resample_gebco_to_target(
    lon_src: np.ndarray,
    lat_src: np.ndarray,
    z_src: np.ndarray,
    bbox: Tuple[float, float, float, float],
    target_res_m: float,
    center_lat: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    minlon, maxlon, minlat, maxlat = bbox
    dlon, dlat = res_m_to_deg(target_res_m, center_lat)
    nx = max(2, int(np.ceil((maxlon - minlon) / dlon)) + 1)
    ny = max(2, int(np.ceil((maxlat - minlat) / dlat)) + 1)
    lon_t = np.linspace(minlon, maxlon, nx)
    lat_t = np.linspace(minlat, maxlat, ny)

    interp = RegularGridInterpolator(
        (lat_src, lon_src), z_src, bounds_error=False, fill_value=np.nan
    )
    yy, xx = np.meshgrid(lat_t, lon_t, indexing="ij")
    z_t = interp(np.stack([yy, xx], axis=-1))
    return lon_t, lat_t, z_t


def srtm_to_target_grid(
    srtm_tif: Path,
    lon_t: np.ndarray,
    lat_t: np.ndarray,
) -> np.ndarray:
    minlon, maxlon = float(lon_t.min()), float(lon_t.max())
    minlat, maxlat = float(lat_t.min()), float(lat_t.max())
    dst_transform = from_bounds(minlon, minlat, maxlon, maxlat, len(lon_t), len(lat_t))
    dst = np.full((len(lat_t), len(lon_t)), np.nan, dtype=np.float32)

    with rasterio.open(srtm_tif) as src:
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=src.nodata,
            dst_transform=dst_transform,
            dst_crs="EPSG:4326",
            dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )
    return dst.astype(float)


def build_local_dem(
    gebco_nc: Path,
    srtm_tif: Optional[Path],
    bbox: Tuple[float, float, float, float],
    target_res_m: float,
    center_lat: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
    lon_g, lat_g, z_g = load_gebco_subset(gebco_nc, bbox)
    lon_t, lat_t, dem = resample_gebco_to_target(lon_g, lat_g, z_g, bbox, target_res_m, center_lat)
    meta: Dict[str, Any] = {
        "gebco_loaded": True,
        "srtm_loaded": False,
        "target_res_m": float(target_res_m),
    }

    if srtm_tif is not None and srtm_tif.exists():
        srtm = srtm_to_target_grid(srtm_tif, lon_t, lat_t)
        # Prefer SRTM where it exists on land or nearshore. Keep GEBCO offshore.
        # SRTM nodata remains NaN; where valid, it replaces GEBCO.
        valid_srtm = np.isfinite(srtm)
        dem = np.where(valid_srtm, srtm, dem)
        meta["srtm_loaded"] = True
        meta["srtm_path"] = str(srtm_tif)
    return lon_t, lat_t, dem, meta


# -----------------------------
# Results loading and inundation computation
# -----------------------------
def _pick_array(d: Dict[str, Any], keys: Tuple[str, ...]) -> Optional[np.ndarray]:
    for k in keys:
        if k in d:
            return np.asarray(d[k])
    return None


def load_results_npz(results_npz: Path) -> Dict[str, np.ndarray]:
    raw = np.load(results_npz, allow_pickle=True)
    data = {k: raw[k] for k in raw.files}
    lon = _pick_array(data, ("lon", "x", "lons", "longitude"))
    lat = _pick_array(data, ("lat", "y", "lats", "latitude"))
    depth_max = _pick_array(data, ("depth_max", "hmax", "inundation_max"))
    eta_max = _pick_array(data, ("eta_max", "etamax", "zeta_max", "surface_max"))

    if lon is None or lat is None:
        raise ValueError("NPZ must contain lon/lat arrays (aliases x/y also accepted).")
    if depth_max is None and eta_max is None:
        raise ValueError("NPZ must contain either depth_max or eta_max.")

    out = {"lon": lon.astype(float), "lat": lat.astype(float)}
    if depth_max is not None:
        out["depth_max"] = np.asarray(depth_max, dtype=float)
    if eta_max is not None:
        out["eta_max"] = np.asarray(eta_max, dtype=float)
    return out


def resample_field_to_target(
    lon_src: np.ndarray,
    lat_src: np.ndarray,
    field_src: np.ndarray,
    lon_t: np.ndarray,
    lat_t: np.ndarray,
) -> np.ndarray:
    lon_src = np.asarray(lon_src).squeeze()
    lat_src = np.asarray(lat_src).squeeze()
    field_src = np.asarray(field_src)

    if field_src.shape != (len(lat_src), len(lon_src)):
        raise ValueError(
            f"Field shape {field_src.shape} does not match lat/lon lengths {(len(lat_src), len(lon_src))}."
        )

    interp = RegularGridInterpolator(
        (lat_src, lon_src), field_src, bounds_error=False, fill_value=np.nan
    )
    yy, xx = np.meshgrid(lat_t, lon_t, indexing="ij")
    return interp(np.stack([yy, xx], axis=-1))


def compute_land_inundation(
    dem: np.ndarray,
    depth_max: Optional[np.ndarray] = None,
    eta_max: Optional[np.ndarray] = None,
    land_threshold_m: float = 0.0,
) -> np.ndarray:
    land = np.isfinite(dem) & (dem >= land_threshold_m)
    out = np.zeros_like(dem, dtype=float)

    if depth_max is not None:
        # Assume incoming depth_max may contain offshore wet depth too.
        out[land] = np.where(np.isfinite(depth_max[land]), np.maximum(depth_max[land], 0.0), 0.0)
        return out

    if eta_max is None:
        raise ValueError("Either depth_max or eta_max must be provided.")

    # Land inundation depth = peak water surface elevation above MSL minus ground elevation.
    out[land] = np.where(np.isfinite(eta_max[land]), np.maximum(eta_max[land] - dem[land], 0.0), 0.0)
    return out


# -----------------------------
# Viability checks
# -----------------------------
def assess_viability(
    dem: np.ndarray,
    land_inund: np.ndarray,
    target_res_m: float,
    dx_m_reported: Optional[float],
    t_end_sec: Optional[float],
) -> ViabilityReport:
    land = np.isfinite(dem) & (dem >= 0.0)
    wet = land & np.isfinite(land_inund) & (land_inund > 0.01)

    cell_area_m2 = target_res_m * target_res_m
    inland_wet_cells = int(np.count_nonzero(wet))
    inland_wet_area_m2 = float(inland_wet_cells * cell_area_m2)
    land_vals = land_inund[wet]
    max_land_inund = float(np.nanmax(land_vals)) if land_vals.size else 0.0
    p95_land_inund = float(np.nanpercentile(land_vals, 95)) if land_vals.size else 0.0

    messages: list[str] = []

    smoke_ok = True
    inund_ok = True

    if inland_wet_cells == 0:
        messages.append("No inland wet cells were detected in the local Dapa domain.")
        inund_ok = False
    else:
        messages.append(f"Detected {inland_wet_cells} inland wet cells in the local Dapa domain.")

    if target_res_m > 120:
        messages.append(
            f"Local map resolution is {target_res_m:.1f} m; this is too coarse for Dapa inundation interpretation. Aim for <= 60 m, preferably 10-30 m."
        )
        inund_ok = False
    elif target_res_m > 60:
        messages.append(
            f"Local map resolution is {target_res_m:.1f} m; usable for screening, but still coarse for publication-style Dapa inundation footprints."
        )
    else:
        messages.append(f"Local map resolution is {target_res_m:.1f} m, which is appropriate for Dapa-focused inundation mapping.")

    if t_end_sec is not None:
        if t_end_sec < 1000:
            messages.append(
                f"Simulation end time is only {t_end_sec:.0f} s; this is generally too short for a confident Dapa inundation assessment."
            )
            inund_ok = False
        elif t_end_sec < 3600:
            messages.append(
                f"Simulation end time is {t_end_sec:.0f} s; acceptable for intermediate checks, but shorter than the paper-style 3-hour setup."
            )
        else:
            messages.append(f"Simulation end time is {t_end_sec:.0f} s, which is suitable for Dapa-focused inundation interpretation.")

    if dx_m_reported is not None:
        if dx_m_reported > 300:
            messages.append(
                f"Reported solver dx is {dx_m_reported:.1f} m. This is fine for regional propagation, but not for paper-like local inundation mapping."
            )
            inund_ok = False
        else:
            messages.append(f"Reported solver dx is {dx_m_reported:.1f} m.")

    if max_land_inund <= 0.0:
        messages.append("Maximum land inundation depth is zero in the processed local domain.")
        inund_ok = False

    # Smoke-test definition: stable-looking output with finite DEM and finite map products.
    if not np.isfinite(dem).any():
        smoke_ok = False
        messages.append("DEM contains no finite values in the local Dapa domain.")

    return ViabilityReport(
        has_srtm=bool(np.nanmax(dem) > 0),
        target_res_m=float(target_res_m),
        dx_m_reported=float(dx_m_reported) if dx_m_reported is not None else None,
        t_end_sec=float(t_end_sec) if t_end_sec is not None else None,
        inland_wet_cells=inland_wet_cells,
        inland_wet_area_m2=inland_wet_area_m2,
        inland_wet_area_ha=inland_wet_area_m2 / 10000.0,
        max_land_inundation_m=max_land_inund,
        p95_land_inundation_m=p95_land_inund,
        local_dem_min_m=float(np.nanmin(dem)),
        local_dem_max_m=float(np.nanmax(dem)),
        is_viable_for_smoke_test=smoke_ok,
        is_viable_for_inundation_interpretation=inund_ok,
        messages=messages,
    )


# -----------------------------
# Plotting
# -----------------------------
def paper_bins() -> np.ndarray:
    # Similar to the paper legend.
    return np.array([0.01, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])


def paper_cmap() -> Tuple[ListedColormap, BoundaryNorm]:
    colors = [
        "#8fd3ff",  # 0.01–0.25
        "#58afff",  # 0.251–0.5
        "#1f83ff",  # 0.51–1.0
        "#0065d9",  # 1.01–1.5
        "#004cb8",  # 1.51–2.0
        "#003a92",  # 2.01–2.5
        "#002e7a",  # 2.51–3.0
        "#00235e",  # 3.01–3.5+
    ]
    bins = paper_bins()
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(bins, cmap.N, clip=True)
    return cmap, norm


def plot_paper_style_noosm(
    lon: np.ndarray,
    lat: np.ndarray,
    dem: np.ndarray,
    land_inund: np.ndarray,
    out_png: Path,
    title: str,
    subtitle: Optional[str] = None,
) -> None:
    cmap, norm = paper_cmap()
    bins = paper_bins()
    xx, yy = np.meshgrid(lon, lat)

    fig = plt.figure(figsize=(10.5, 7.5), dpi=180)
    ax = plt.axes()

    # Sea/land background close to the paper look.
    sea = np.ma.masked_where(~np.isfinite(dem) | (dem >= 0), dem)
    land = np.ma.masked_where(~np.isfinite(dem) | (dem < 0), dem)

    ax.pcolormesh(xx, yy, sea, shading="auto", cmap=ListedColormap(["#cfe8f7"]))
    ax.pcolormesh(xx, yy, land, shading="auto", cmap=ListedColormap(["#e6e6e6"]))

    inund_masked = np.ma.masked_where(~np.isfinite(land_inund) | (land_inund < bins[0]), land_inund)
    pcm = ax.pcolormesh(xx, yy, inund_masked, shading="auto", cmap=cmap, norm=norm)

    ax.set_title(title, fontsize=13, weight="bold", pad=10)
    if subtitle:
        ax.text(0.01, 0.985, subtitle, transform=ax.transAxes, ha="left", va="top", fontsize=10)

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    cbar = fig.colorbar(pcm, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Inundation depth (m)")
    tick_centers = [
        0.13,
        0.375,
        0.75,
        1.25,
        1.75,
        2.25,
        2.75,
        3.25,
    ]
    cbar.set_ticks(tick_centers)
    cbar.ax.set_yticklabels([
        "0.01–0.25",
        "0.251–0.5",
        "0.51–1.0",
        "1.01–1.5",
        "1.51–2.0",
        "2.01–2.5",
        "2.51–3.0",
        "3.01–3.5+",
    ])

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


# -----------------------------
# CLI
# -----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Dapa no-OSM inundation post-processor and mapper.")
    p.add_argument("--gebco", required=True, type=Path, help="Path to GEBCO netCDF.")
    p.add_argument("--srtm", type=Path, default=None, help="Optional SRTM/DEM GeoTIFF for local land refinement.")
    p.add_argument("--results", type=Path, required=True, help="NPZ with lon/lat and depth_max or eta_max.")
    p.add_argument("--out-prefix", type=Path, required=True, help="Output prefix, e.g. /path/to/out/Mw85_dapa")
    p.add_argument("--center-lon", type=float, default=DAPA_LON)
    p.add_argument("--center-lat", type=float, default=DAPA_LAT)
    p.add_argument("--half-width-km", type=float, default=DEFAULT_HALF_WIDTH_KM)
    p.add_argument("--half-height-km", type=float, default=DEFAULT_HALF_HEIGHT_KM)
    p.add_argument("--target-res-m", type=float, default=30.0)
    p.add_argument("--dx-m", type=float, default=None, help="Reported solver dx in meters for viability reporting.")
    p.add_argument("--t-end-sec", type=float, default=None, help="Simulation end time in seconds for viability reporting.")
    p.add_argument("--title", type=str, default="Dapa tsunami inundation (no OSM overlay)")
    p.add_argument("--subtitle", type=str, default=None)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    bbox = make_bbox(
        center_lon=args.center_lon,
        center_lat=args.center_lat,
        half_width_km=args.half_width_km,
        half_height_km=args.half_height_km,
    )

    lon_t, lat_t, dem, dem_meta = build_local_dem(
        gebco_nc=args.gebco,
        srtm_tif=args.srtm,
        bbox=bbox,
        target_res_m=args.target_res_m,
        center_lat=args.center_lat,
    )

    results = load_results_npz(args.results)
    lon_r = np.asarray(results["lon"], dtype=float)
    lat_r = np.asarray(results["lat"], dtype=float)

    depth_r = results.get("depth_max")
    eta_r = results.get("eta_max")
    depth_t = None
    eta_t = None
    if depth_r is not None:
        depth_t = resample_field_to_target(lon_r, lat_r, depth_r, lon_t, lat_t)
    if eta_r is not None:
        eta_t = resample_field_to_target(lon_r, lat_r, eta_r, lon_t, lat_t)

    land_inund = compute_land_inundation(dem, depth_max=depth_t, eta_max=eta_t)
    report = assess_viability(
        dem=dem,
        land_inund=land_inund,
        target_res_m=args.target_res_m,
        dx_m_reported=args.dx_m,
        t_end_sec=args.t_end_sec,
    )

    out_prefix = args.out_prefix
    out_png = Path(f"{out_prefix}.png")
    out_json = Path(f"{out_prefix}_viability.json")
    out_dem = Path(f"{out_prefix}_local_dem.npz")
    out_inund = Path(f"{out_prefix}_land_inundation.npz")

    subtitle = args.subtitle
    if subtitle is None:
        subtitle = f"Local DEM res={args.target_res_m:.0f} m | SRTM={'yes' if dem_meta['srtm_loaded'] else 'no'}"

    plot_paper_style_noosm(
        lon=lon_t,
        lat=lat_t,
        dem=dem,
        land_inund=land_inund,
        out_png=out_png,
        title=args.title,
        subtitle=subtitle,
    )

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(asdict(report), f, indent=2)

    np.savez_compressed(out_dem, lon=lon_t, lat=lat_t, dem=dem)
    np.savez_compressed(out_inund, lon=lon_t, lat=lat_t, land_inundation=land_inund)

    print(f"[OK] Map saved: {out_png}")
    print(f"[OK] Viability report saved: {out_json}")
    print(f"[OK] Local DEM saved: {out_dem}")
    print(f"[OK] Land inundation saved: {out_inund}")
    print(f"[INFO] Viable for smoke test: {report.is_viable_for_smoke_test}")
    print(f"[INFO] Viable for inundation interpretation: {report.is_viable_for_inundation_interpretation}")
    for msg in report.messages:
        print(f"[NOTE] {msg}")


if __name__ == "__main__":
    main()
