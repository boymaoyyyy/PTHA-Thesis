#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

try:
    import rasterio
    from rasterio.windows import from_bounds
    from rasterio.enums import Resampling
    from rasterio.warp import reproject
except Exception as e:
    raise SystemExit("This script requires rasterio. Install with: pip install rasterio") from e

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:
    plt = None

try:
    from scipy.ndimage import (
        binary_dilation,
        binary_erosion,
        binary_opening,
        binary_closing,
        distance_transform_edt,
        label,
    )
except Exception as e:
    raise SystemExit("This script requires scipy. Install with: pip install scipy") from e


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Create a fine-grid merged Dapa topo-bathy DEM with shoreline cleanup, "
            "QA checks, and Dapa gauge validation."
        )
    )
    p.add_argument("--gebco", type=Path, required=True, help="Path to GEBCO GeoTIFF")
    p.add_argument("--topo", type=Path, required=True, help="Path to Dapa topography GeoTIFF")
    p.add_argument("--xmin", type=float, default=125.95)
    p.add_argument("--ymin", type=float, default=9.65)
    p.add_argument("--xmax", type=float, default=126.15)
    p.add_argument("--ymax", type=float, default=9.85)
    p.add_argument("--blend-cells", type=int, default=8, help="Width of shoreline blending halo in grid cells")
    p.add_argument("--shore-tolerance-m", type=float, default=0.50, help="Treat |z| <= tolerance as shoreline band")
    p.add_argument("--offshore-buffer-m", type=float, default=0.0, help="Optional vertical offset to apply to bathymetry values")
    p.add_argument("--cleanup-min-region-cells", type=int, default=9, help="Remove tiny isolated land/water speckles smaller than this")
    p.add_argument("--gauge-lon", type=float, default=126.055)
    p.add_argument("--gauge-lat", type=float, default=9.745)
    p.add_argument("--gauge-min-depth-m", type=float, default=0.50)
    p.add_argument("--gauge-max-distance-km", type=float, default=2.00)
    p.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("output/dapa_merged_topobathy_fine_qc"),
        help="Output prefix without extension",
    )
    return p.parse_args()


def sanitize_nodata(arr: np.ndarray, nodata_value) -> np.ndarray:
    out = arr.astype(np.float32)
    if nodata_value is not None:
        out[np.isclose(out, nodata_value)] = np.nan
    return out


def crop_topo_to_bounds(
    topo_path: Path,
    bounds: Tuple[float, float, float, float],
) -> Tuple[np.ndarray, dict]:
    xmin, ymin, xmax, ymax = bounds

    with rasterio.open(topo_path) as ds:
        if ds.crs is None:
            raise ValueError("Topography raster has no CRS.")
        if ds.count < 1:
            raise ValueError("Topography raster has no bands.")

        window = from_bounds(xmin, ymin, xmax, ymax, ds.transform)
        window = window.round_offsets().round_lengths()

        if window.width < 2 or window.height < 2:
            raise ValueError("Topography crop is too small. Check the bounds.")

        topo = ds.read(1, window=window).astype(np.float32)
        topo = sanitize_nodata(topo, ds.nodata)
        transform = ds.window_transform(window)

        profile = ds.profile.copy()
        profile.update(
            driver="GTiff",
            dtype="float32",
            count=1,
            height=topo.shape[0],
            width=topo.shape[1],
            transform=transform,
            crs=ds.crs,
            nodata=np.nan,
            compress="deflate",
        )

    return topo, profile


def reproject_gebco_to_topo_grid(
    gebco_path: Path,
    target_profile: dict,
) -> np.ndarray:
    dst = np.full(
        (target_profile["height"], target_profile["width"]),
        np.nan,
        dtype=np.float32,
    )

    with rasterio.open(gebco_path) as src:
        if src.crs is None:
            raise ValueError("GEBCO raster has no CRS.")
        if src.count < 1:
            raise ValueError("GEBCO raster has no bands.")

        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=src.nodata,
            dst_transform=target_profile["transform"],
            dst_crs=target_profile["crs"],
            dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )

    return dst.astype(np.float32)


def remove_small_regions(mask: np.ndarray, min_cells: int) -> np.ndarray:
    if min_cells <= 1:
        return mask.copy()
    labeled, nlab = label(mask)
    if nlab == 0:
        return mask.copy()
    out = np.zeros_like(mask, dtype=bool)
    counts = np.bincount(labeled.ravel())
    for lab_id in range(1, len(counts)):
        if counts[lab_id] >= min_cells:
            out[labeled == lab_id] = True
    return out


def derive_masks(
    topo: np.ndarray,
    bathy: np.ndarray,
    shore_tolerance_m: float,
    cleanup_min_region_cells: int,
) -> Dict[str, np.ndarray]:
    valid_topo = np.isfinite(topo)
    valid_bathy = np.isfinite(bathy)
    valid_any = valid_topo | valid_bathy

    land_seed = valid_topo & (topo > shore_tolerance_m)
    water_seed = valid_bathy & (bathy < -shore_tolerance_m)

    shoreline_band = valid_any & ~(land_seed | water_seed)

    land = binary_closing(binary_opening(land_seed, iterations=1), iterations=1)
    water = binary_closing(binary_opening(water_seed, iterations=1), iterations=1)

    land = remove_small_regions(land, cleanup_min_region_cells)
    water = remove_small_regions(water, cleanup_min_region_cells)

    conflict = land & water
    if np.any(conflict):
        water[conflict] = False

    unresolved = valid_any & ~(land | water)
    shore = unresolved | shoreline_band

    return {
        "valid_any": valid_any,
        "land": land,
        "water": water,
        "shore": shore,
    }


def merge_topobathy(
    topo: np.ndarray,
    gebco_on_topo: np.ndarray,
    blend_cells: int = 8,
    shore_tolerance_m: float = 0.50,
    offshore_buffer_m: float = 0.0,
    cleanup_min_region_cells: int = 9,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    merged = np.full_like(topo, np.nan, dtype=np.float32)

    bathy = gebco_on_topo.copy()
    if offshore_buffer_m != 0.0:
        bathy = bathy + np.float32(offshore_buffer_m)

    masks = derive_masks(
        topo=topo,
        bathy=bathy,
        shore_tolerance_m=shore_tolerance_m,
        cleanup_min_region_cells=cleanup_min_region_cells,
    )

    land = masks["land"]
    water = masks["water"]
    shore = masks["shore"]
    valid_any = masks["valid_any"]

    merged[land] = topo[land]
    merged[water] = bathy[water]

    # Initial fill for unresolved shoreline cells:
    # prefer topo where slightly positive, bathy where slightly negative, otherwise average if both exist
    unresolved = shore & np.isnan(merged)

    topo_pref = unresolved & np.isfinite(topo) & (topo > 0.0)
    bathy_pref = unresolved & np.isfinite(bathy) & (bathy < 0.0)
    both = unresolved & np.isfinite(topo) & np.isfinite(bathy) & ~(topo_pref | bathy_pref)

    merged[topo_pref] = topo[topo_pref]
    merged[bathy_pref] = bathy[bathy_pref]
    merged[both] = 0.5 * (topo[both] + bathy[both])

    fill_from_topo = np.isnan(merged) & np.isfinite(topo)
    merged[fill_from_topo] = topo[fill_from_topo]
    fill_from_bathy = np.isnan(merged) & np.isfinite(bathy)
    merged[fill_from_bathy] = bathy[fill_from_bathy]

    if blend_cells > 0 and np.any(valid_any):
        coastal_band = binary_dilation(shore | binary_dilation(land, iterations=1), iterations=blend_cells)
        coastal_band &= valid_any & np.isfinite(merged) & np.isfinite(topo) & np.isfinite(bathy)

        if np.any(coastal_band):
            dist_to_land = distance_transform_edt(~land).astype(np.float32)
            dist_to_water = distance_transform_edt(~water).astype(np.float32)
            denom = np.maximum(dist_to_land + dist_to_water, 1.0e-6)
            land_weight = np.clip(dist_to_water / denom, 0.0, 1.0)

            blend_mask = coastal_band & ~(land ^ water)
            merged[blend_mask] = (
                land_weight[blend_mask] * topo[blend_mask]
                + (1.0 - land_weight[blend_mask]) * bathy[blend_mask]
            )

    return merged.astype(np.float32), masks


def compute_lon_lat(profile: dict) -> Tuple[np.ndarray, np.ndarray]:
    t = profile["transform"]
    width = profile["width"]
    height = profile["height"]
    lon = np.array([t.c + (i + 0.5) * t.a for i in range(width)], dtype=np.float64)
    lat = np.array([t.f + (j + 0.5) * t.e for j in range(height)], dtype=np.float64)
    return lon, lat


def nearest_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(values - target)))


def estimate_cellsize_m(profile: dict, lat_center_deg: float) -> Tuple[float, float]:
    dx_deg = float(profile["transform"].a)
    dy_deg = float(abs(profile["transform"].e))
    dx_m = abs(dx_deg) * 111_320.0 * np.cos(np.deg2rad(lat_center_deg))
    dy_m = abs(dy_deg) * 110_540.0
    return float(dx_m), float(dy_m)


def validate_gauge(
    merged: np.ndarray,
    profile: dict,
    gauge_lon: float,
    gauge_lat: float,
    min_depth_m: float,
    max_distance_km: float,
) -> Dict[str, object]:
    lon, lat = compute_lon_lat(profile)
    H = np.maximum(-merged, 0.0)

    i0 = nearest_index(lat, gauge_lat)
    j0 = nearest_index(lon, gauge_lon)

    wet = H > min_depth_m
    selected = None
    nearest_km = None

    if wet[i0, j0]:
        selected = (i0, j0)
        nearest_km = 0.0
    else:
        wet_idx = np.argwhere(wet)
        if wet_idx.size > 0:
            best_dist = float("inf")
            best = None
            for wi, wj in wet_idx:
                dlat = float(lat[int(wi)] - gauge_lat)
                dlon = float(lon[int(wj)] - gauge_lon)
                dist_deg = (dlat * dlat + dlon * dlon) ** 0.5
                if dist_deg < best_dist:
                    best_dist = dist_deg
                    best = (int(wi), int(wj))
            if best is not None:
                selected = best
                lat_mid = 0.5 * (float(lat[selected[0]]) + gauge_lat)
                km_per_deg_lon = 111.320 * np.cos(np.deg2rad(lat_mid))
                km_per_deg_lat = 110.540
                nearest_km = (
                    ((float(lat[selected[0]]) - gauge_lat) * km_per_deg_lat) ** 2
                    + ((float(lon[selected[1]]) - gauge_lon) * km_per_deg_lon) ** 2
                ) ** 0.5

    valid = bool(selected is not None and nearest_km is not None and nearest_km <= max_distance_km)

    if selected is None:
        return {
            "target_lon": gauge_lon,
            "target_lat": gauge_lat,
            "selected_lon": None,
            "selected_lat": None,
            "selected_depth_m": None,
            "nearest_wet_km": None,
            "valid": False,
            "reason": "No wet cell found above minimum depth threshold in the DEM.",
        }

    ii, jj = selected
    return {
        "target_lon": gauge_lon,
        "target_lat": gauge_lat,
        "selected_lon": float(lon[jj]),
        "selected_lat": float(lat[ii]),
        "selected_depth_m": float(H[ii, jj]),
        "nearest_wet_km": float(nearest_km),
        "valid": valid,
        "reason": (
            "Gauge location is acceptable for inundation modeling."
            if valid
            else f"Nearest wet cell is too far from target gauge (> {max_distance_km:.2f} km)."
        ),
    }


def build_qa(
    merged: np.ndarray,
    masks: Dict[str, np.ndarray],
    profile: dict,
    gauge_report: Dict[str, object],
    shore_tolerance_m: float,
) -> Dict[str, object]:
    lon, lat = compute_lon_lat(profile)
    lat_center = 0.5 * (float(np.min(lat)) + float(np.max(lat)))
    dx_m, dy_m = estimate_cellsize_m(profile, lat_center)

    H = np.maximum(-merged, 0.0)
    land = masks["land"]
    water = masks["water"]
    shore = masks["shore"]
    valid_any = masks["valid_any"]

    wet_cells = int(np.count_nonzero(H > 0.50))
    land_cells = int(np.count_nonzero(land))
    water_cells = int(np.count_nonzero(water))
    shore_cells = int(np.count_nonzero(shore))

    min_m = float(np.nanmin(merged))
    max_m = float(np.nanmax(merged))
    valid_fraction = float(np.mean(np.isfinite(merged)))

    checks = {
        "has_negative_bathymetry": bool(min_m < 0.0),
        "has_positive_topography": bool(max_m > 0.0),
        "has_shoreline_band": bool(shore_cells > 0),
        "has_minimum_wet_cells": bool(wet_cells >= 100),
        "gauge_is_valid": bool(gauge_report.get("valid", False)),
        "full_valid_coverage": bool(valid_fraction >= 0.99),
    }

    suitable = bool(all(checks.values()))

    return {
        "shape": [int(profile["height"]), int(profile["width"])],
        "xres_deg": float(profile["transform"].a),
        "yres_deg": float(abs(profile["transform"].e)),
        "dx_m": dx_m,
        "dy_m": dy_m,
        "shore_tolerance_m": float(shore_tolerance_m),
        "min_m": min_m,
        "max_m": max_m,
        "mean_m": float(np.nanmean(merged)),
        "valid_fraction": valid_fraction,
        "land_cells": land_cells,
        "water_cells": water_cells,
        "shore_cells": shore_cells,
        "wet_cells_gt_0p5m": wet_cells,
        "suitability_checks": checks,
        "suitable_for_inundation_modeling": suitable,
        "gauge_validation": gauge_report,
    }


def write_geotiff(path: Path, arr: np.ndarray, profile: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    out_profile = profile.copy()
    out_profile.update(driver="GTiff", dtype="float32", count=1, compress="deflate")
    with rasterio.open(path, "w", **out_profile) as dst:
        dst.write(arr.astype(np.float32), 1)


def write_npz(path: Path, arr: np.ndarray, profile: dict) -> None:
    lon, lat = compute_lon_lat(profile)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, lon=lon, lat=lat, z=arr.astype(np.float32))


def save_preview(path: Path, arr: np.ndarray, profile: dict, title: str) -> None:
    if plt is None:
        return

    t = profile["transform"]
    width = profile["width"]
    height = profile["height"]
    xmin = t.c
    ymax = t.f
    xmax = xmin + width * t.a
    ymin = ymax + height * t.e
    extent = [xmin, xmax, ymin, ymax]

    fig, ax = plt.subplots(figsize=(7.5, 7.0), dpi=180)
    im = ax.imshow(arr, origin="upper", extent=extent, aspect="equal")
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def save_mask_png(path: Path, mask: np.ndarray, profile: dict, title: str) -> None:
    if plt is None:
        return
    save_preview(path, mask.astype(np.float32), profile, title)


def save_json(path: Path, payload: Dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def save_summary(
    path: Path,
    gebco_path: Path,
    topo_path: Path,
    bounds: Tuple[float, float, float, float],
    merged: np.ndarray,
    profile: dict,
) -> None:
    summary = {
        "gebco": str(gebco_path),
        "topography": str(topo_path),
        "bounds": {
            "xmin": float(bounds[0]),
            "ymin": float(bounds[1]),
            "xmax": float(bounds[2]),
            "ymax": float(bounds[3]),
        },
        "shape": [int(profile["height"]), int(profile["width"])],
        "crs": str(profile["crs"]),
        "xres_deg": float(profile["transform"].a),
        "yres_deg": float(abs(profile["transform"].e)),
        "min_m": float(np.nanmin(merged)),
        "max_m": float(np.nanmax(merged)),
        "mean_m": float(np.nanmean(merged)),
        "valid_fraction": float(np.mean(np.isfinite(merged))),
    }
    save_json(path, summary)


def main() -> None:
    args = parse_args()
    bounds = (args.xmin, args.ymin, args.xmax, args.ymax)

    gebco_path = args.gebco.resolve()
    topo_path = args.topo.resolve()

    if not gebco_path.exists():
        raise FileNotFoundError(f"GEBCO not found: {gebco_path}")
    if not topo_path.exists():
        raise FileNotFoundError(f"Topography not found: {topo_path}")

    print(f"[GEBCO] {gebco_path}")
    print(f"[TOPO]  {topo_path}")
    print(f"[BBOX]  xmin={args.xmin}, ymin={args.ymin}, xmax={args.xmax}, ymax={args.ymax}")

    topo_crop, target_profile = crop_topo_to_bounds(topo_path, bounds)
    gebco_on_topo = reproject_gebco_to_topo_grid(gebco_path, target_profile)

    merged, masks = merge_topobathy(
        topo=topo_crop,
        gebco_on_topo=gebco_on_topo,
        blend_cells=args.blend_cells,
        shore_tolerance_m=args.shore_tolerance_m,
        offshore_buffer_m=args.offshore_buffer_m,
        cleanup_min_region_cells=args.cleanup_min_region_cells,
    )

    gauge_report = validate_gauge(
        merged=merged,
        profile=target_profile,
        gauge_lon=args.gauge_lon,
        gauge_lat=args.gauge_lat,
        min_depth_m=args.gauge_min_depth_m,
        max_distance_km=args.gauge_max_distance_km,
    )

    qa = build_qa(
        merged=merged,
        masks=masks,
        profile=target_profile,
        gauge_report=gauge_report,
        shore_tolerance_m=args.shore_tolerance_m,
    )

    tif_path = args.out_prefix.with_suffix(".tif")
    npz_path = args.out_prefix.with_suffix(".npz")
    png_path = args.out_prefix.with_name(args.out_prefix.name + "_preview.png")
    land_mask_png = args.out_prefix.with_name(args.out_prefix.name + "_land_mask.png")
    water_mask_png = args.out_prefix.with_name(args.out_prefix.name + "_water_mask.png")
    shoreline_png = args.out_prefix.with_name(args.out_prefix.name + "_shoreline_band.png")
    summary_json = args.out_prefix.with_name(args.out_prefix.name + "_summary.json")
    gauge_json = args.out_prefix.with_name(args.out_prefix.name + "_gauge_validation.json")
    qa_json = args.out_prefix.with_name(args.out_prefix.name + "_qa.json")

    write_geotiff(tif_path, merged, target_profile)
    write_npz(npz_path, merged, target_profile)
    save_preview(png_path, merged, target_profile, "Fine-Grid Merged Dapa Topo-Bathy DEM")
    save_mask_png(land_mask_png, masks["land"], target_profile, "Dapa DEM Land Mask")
    save_mask_png(water_mask_png, masks["water"], target_profile, "Dapa DEM Water Mask")
    save_mask_png(shoreline_png, masks["shore"], target_profile, "Dapa DEM Shoreline Band")
    save_summary(summary_json, gebco_path, topo_path, bounds, merged, target_profile)
    save_json(gauge_json, gauge_report)
    save_json(qa_json, qa)

    print(f"[DONE] Wrote: {tif_path}")
    print(f"[DONE] Wrote: {npz_path}")
    print(f"[DONE] Wrote: {png_path}")
    print(f"[DONE] Wrote: {land_mask_png}")
    print(f"[DONE] Wrote: {water_mask_png}")
    print(f"[DONE] Wrote: {shoreline_png}")
    print(f"[DONE] Wrote: {summary_json}")
    print(f"[DONE] Wrote: {gauge_json}")
    print(f"[DONE] Wrote: {qa_json}")
    print(f"[QA] suitable_for_inundation_modeling = {qa['suitable_for_inundation_modeling']}")
    print(f"[QA] gauge_valid = {gauge_report['valid']} | nearest_wet_km = {gauge_report['nearest_wet_km']}")


if __name__ == "__main__":
    main()