#!/usr/bin/env python3
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse
import copy
import io
import json
import math
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout, redirect_stderr
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except Exception:
    tqdm = None

import phase3_thesis_runner_dapa_nested_patched as base


# Force this runner to write into a dedicated production folder.
base.ProjectConfig.results_dir = property(
    lambda self: Path(
        getattr(self, "_results_dir_override", self.project_root / "output" / "phase3_final_outputs")
    )
)


def set_results_dir(cfg: base.ProjectConfig, path: Path) -> None:
    cfg._results_dir_override = path


_WORKER_CFG = None
_WORKER_DOMAIN = None
_WORKER_MW = None
_WORKER_RETURN_DEPTH = False


@dataclass
class GaugeSupport:
    name: str
    lon: float
    lat: float
    rows: np.ndarray
    cols: np.ndarray
    weights: np.ndarray
    nearest_wet_km: float
    nearest_wet_lon: float
    nearest_wet_lat: float
    valid: bool
    support_radius_km: float


def km_per_deg_lon(lat_deg: float) -> float:
    return 111.320 * math.cos(math.radians(lat_deg))


def build_gauge_support(
    grid: base.Grid2D,
    H: np.ndarray,
    gauges: Sequence[base.GaugePoint],
    min_depth: float = 0.50,
    support_radius_km: float = 8.0,
    max_support_cells: int = 8,
    valid_radius_km: Optional[float] = None,
) -> Dict[str, GaugeSupport]:
    """
    Build a wet-cell support stencil around each virtual gauge.
    This avoids silently snapping a gauge far away and calling it the same point.
    """
    if valid_radius_km is None:
        valid_radius_km = max(5.0, 3.0 * max(grid.dx_m, grid.dy_m) / 1000.0)

    wet_idx = np.argwhere(H > min_depth)
    out: Dict[str, GaugeSupport] = {}

    for g in gauges:
        if wet_idx.size == 0:
            i0 = base.nearest_index(grid.lat, g.lat)
            j0 = base.nearest_index(grid.lon, g.lon)
            out[g.name] = GaugeSupport(
                name=g.name,
                lon=g.lon,
                lat=g.lat,
                rows=np.asarray([i0], dtype=int),
                cols=np.asarray([j0], dtype=int),
                weights=np.asarray([1.0], dtype=float),
                nearest_wet_km=float("inf"),
                nearest_wet_lon=float(grid.lon[j0]),
                nearest_wet_lat=float(grid.lat[i0]),
                valid=False,
                support_radius_km=support_radius_km,
            )
            continue

        wet_rows = wet_idx[:, 0].astype(int)
        wet_cols = wet_idx[:, 1].astype(int)
        dlat_km = (grid.lat[wet_rows] - g.lat) * 110.540
        dlon_km = (grid.lon[wet_cols] - g.lon) * km_per_deg_lon(g.lat)
        dist_km = np.hypot(dlat_km, dlon_km)
        order = np.argsort(dist_km)
        nearest = int(order[0])
        nearest_km = float(dist_km[nearest])

        within = order[dist_km[order] <= support_radius_km]
        chosen = within[:max_support_cells] if within.size else order[:max_support_cells]
        chosen_dist = np.asarray(dist_km[chosen], dtype=float)
        if np.any(chosen_dist < 0.10):
            w = np.zeros_like(chosen_dist)
            w[np.argmin(chosen_dist)] = 1.0
        else:
            w = 1.0 / np.maximum(chosen_dist, 0.25) ** 2
            w /= np.sum(w)

        out[g.name] = GaugeSupport(
            name=g.name,
            lon=g.lon,
            lat=g.lat,
            rows=wet_rows[chosen],
            cols=wet_cols[chosen],
            weights=w.astype(float),
            nearest_wet_km=nearest_km,
            nearest_wet_lon=float(grid.lon[wet_cols[nearest]]),
            nearest_wet_lat=float(grid.lat[wet_rows[nearest]]),
            valid=bool(nearest_km <= valid_radius_km),
            support_radius_km=support_radius_km,
        )
    return out


def sample_gauge_from_support(eta: np.ndarray, support: GaugeSupport) -> float:
    vals = np.asarray(eta[support.rows, support.cols], dtype=float)
    if vals.size == 0:
        return float("nan")
    return float(np.sum(vals * support.weights))


def nanmax_abs(arr: np.ndarray) -> float:
    arr = np.asarray(arr, dtype=float)
    if arr.size == 0 or not np.any(np.isfinite(arr)):
        return float("nan")
    return float(np.nanmax(np.abs(arr)))


def finite_or_none(v: object) -> Optional[float]:
    try:
        vf = float(v)
    except Exception:
        return None
    return vf if np.isfinite(vf) else None


def finite_or_zero(v: object) -> float:
    try:
        vf = float(v)
    except Exception:
        return 0.0
    return vf if np.isfinite(vf) else 0.0


def patch_base_module() -> None:
    def prepare_domain_patched(cfg: base.ProjectConfig) -> base.PreparedDomain:
        grid, topo_local, local_merged = base.load_combined_domain(cfg)
        z = grid.z.astype(np.float64)
        for _ in range(2):
            z = 0.5 * z + 0.125 * (
                np.roll(z, 1, axis=0) + np.roll(z, -1, axis=0) + np.roll(z, 1, axis=1) + np.roll(z, -1, axis=1)
            )
        H = np.maximum(-z, 0.0)

        dapa_nominal_ys, dapa_nominal_xs = base.find_bbox_indices(
            grid.lon,
            grid.lat,
            cfg.domain.dapa_lon_min,
            cfg.domain.dapa_lon_max,
            cfg.domain.dapa_lat_min,
            cfg.domain.dapa_lat_max,
        )
        dapa_nominal_info = base.describe_subgrid(grid, H, dapa_nominal_ys, dapa_nominal_xs)

        if local_merged is not None:
            lon_min, lon_max = float(local_merged.lon[0]), float(local_merged.lon[-1])
            lat_min, lat_max = float(local_merged.lat[0]), float(local_merged.lat[-1])
        else:
            lon_min, lon_max = cfg.domain.dapa_lon_min, cfg.domain.dapa_lon_max
            lat_min, lat_max = cfg.domain.dapa_lat_min, cfg.domain.dapa_lat_max

        dapa_ys, dapa_xs = base.find_bbox_indices(grid.lon, grid.lat, lon_min, lon_max, lat_min, lat_max)
        dapa_info = base.describe_subgrid(grid, H, dapa_ys, dapa_xs)
        attempt = 0
        while dapa_info["wet_cells_gt_min"] <= 0 and attempt < 4:
            lon_min, lon_max, lat_min, lat_max = base.expand_bbox(lon_min, lon_max, lat_min, lat_max, factor=1.40)
            dapa_ys, dapa_xs = base.find_bbox_indices(grid.lon, grid.lat, lon_min, lon_max, lat_min, lat_max)
            dapa_info = base.describe_subgrid(grid, H, dapa_ys, dapa_xs)
            attempt += 1

        gauges = base.sample_gauges(grid, H, cfg.gauges, min_depth=max(cfg.solver.min_depth_m, 0.50))
        gauge_support = build_gauge_support(grid, H, cfg.gauges, min_depth=max(cfg.solver.min_depth_m, 0.50))

        dapa_support = gauge_support.get("Dapa_virtual")
        dapa_window_valid = bool(
            dapa_nominal_info["wet_cells_gt_min"] >= 20
            and dapa_support is not None
            and dapa_support.valid
        )

        print(
            f"[DAPA nominal] lon=[{grid.lon[dapa_nominal_xs][0]:.5f},{grid.lon[dapa_nominal_xs][-1]:.5f}] "
            f"lat=[{grid.lat[dapa_nominal_ys][0]:.5f},{grid.lat[dapa_nominal_ys][-1]:.5f}] "
            f"shape={H[dapa_nominal_ys, dapa_nominal_xs].shape} wet_cells={dapa_nominal_info['wet_cells']} "
            f"wet_cells_gt_min={dapa_nominal_info['wet_cells_gt_min']} max_static_depth={dapa_nominal_info['max_static_depth']:.2f}m",
            flush=True,
        )
        print(
            f"[DAPA effective] lon=[{grid.lon[dapa_xs][0]:.5f},{grid.lon[dapa_xs][-1]:.5f}] "
            f"lat=[{grid.lat[dapa_ys][0]:.5f},{grid.lat[dapa_ys][-1]:.5f}] "
            f"shape={H[dapa_ys, dapa_xs].shape} wet_cells={dapa_info['wet_cells']} "
            f"wet_cells_gt_min={dapa_info['wet_cells_gt_min']} max_static_depth={dapa_info['max_static_depth']:.2f}m",
            flush=True,
        )
        for g in cfg.gauges:
            ii, jj = gauges[g.name]
            support = gauge_support[g.name]
            print(
                f"[GAUGE] {g.name}: target=({g.lon:.5f},{g.lat:.5f}) snapped=({grid.lon[jj]:.5f},{grid.lat[ii]:.5f}) "
                f"H={H[ii, jj]:.2f}m z={z[ii, jj]:.2f}m nearestWet={support.nearest_wet_km:.2f}km "
                f"valid={support.valid}",
                flush=True,
            )
        if dapa_support is not None and not dapa_window_valid:
            print(
                "[DAPA] coarse Dapa outputs will be validity-masked because the nominal window is under-resolved and/or the virtual gauge is too far from Dapa.",
                flush=True,
            )

        if topo_local is not None:
            topo_local = base.crop_grid(
                topo_local,
                grid.lon[dapa_ys][0],
                grid.lon[dapa_ys][-1],
                grid.lat[dapa_ys][0],
                grid.lat[dapa_ys][-1],
            )

        nested_local_domain = None
        if cfg.execution.local_dapa_nested and local_merged is not None:
            local_z = local_merged.z.astype(np.float64)
            local_H = np.maximum(-local_z, 0.0)
            local_gauges = base.sample_gauges(
                local_merged,
                local_H,
                (base.GaugePoint("Dapa_virtual", 126.055, 9.745),),
                min_depth=max(cfg.solver.min_depth_m, 0.10),
            )
            nested_local_domain = base.PreparedDomain(
                grid=local_merged,
                H=local_H,
                z=local_z,
                dx=local_merged.dx_m,
                dy=local_merged.dy_m,
                gauge_idx=local_gauges,
                dapa_ys=slice(0, local_merged.shape[0]),
                dapa_xs=slice(0, local_merged.shape[1]),
                topo_grid_local=None,
                nested_local_domain=None,
            )
            setattr(
                nested_local_domain,
                "gauge_support",
                build_gauge_support(
                    local_merged,
                    local_H,
                    (base.GaugePoint("Dapa_virtual", 126.055, 9.745),),
                    min_depth=max(cfg.solver.min_depth_m, 0.10),
                ),
            )
            setattr(nested_local_domain, "dapa_nominal_ys", slice(0, local_merged.shape[0]))
            setattr(nested_local_domain, "dapa_nominal_xs", slice(0, local_merged.shape[1]))
            setattr(nested_local_domain, "dapa_window_valid", True)
            print(
                f"[LOCAL NEST] shape={local_merged.shape} dx={nested_local_domain.dx:.1f}m dy={nested_local_domain.dy:.1f}m "
                f"wet_cells={int(np.count_nonzero(local_H > cfg.solver.min_depth_m))}",
                flush=True,
            )

        domain = base.PreparedDomain(
            grid=grid,
            H=H,
            z=z,
            dx=grid.dx_m,
            dy=grid.dy_m,
            gauge_idx=gauges,
            dapa_ys=dapa_ys,
            dapa_xs=dapa_xs,
            topo_grid_local=topo_local,
            nested_local_domain=nested_local_domain,
        )
        setattr(domain, "gauge_support", gauge_support)
        setattr(domain, "dapa_nominal_ys", dapa_nominal_ys)
        setattr(domain, "dapa_nominal_xs", dapa_nominal_xs)
        setattr(domain, "dapa_nominal_info", dapa_nominal_info)
        setattr(domain, "dapa_effective_info", dapa_info)
        setattr(domain, "dapa_window_valid", dapa_window_valid)
        return domain

    def run_tsunami_simulation_patched(
        cfg: base.ProjectConfig,
        domain: base.PreparedDomain,
        slip: np.ndarray,
        slip_id: str,
        sample_index: int,
        pbar=None,
    ):
        ny, nx = domain.grid.shape
        q = np.zeros((3, ny, nx), dtype=np.float64)
        max_eta = np.full((ny, nx), -np.inf, dtype=np.float32)
        max_abs_eta = np.zeros((ny, nx), dtype=np.float32)
        max_depth = np.zeros((ny, nx), dtype=np.float32)
        max_speed = np.zeros((ny, nx), dtype=np.float32)
        uplift0 = base.build_uplift_field(cfg, domain, slip)
        rise_time = cfg.source.rise_time_sec_per_km * cfg.source.rupture_length_km
        dt0 = base.stable_dt(domain, cfg)

        if sample_index == 0:
            base.maybe_print(f"Grid: {ny} x {nx}  dx={domain.dx:.0f}m", pbar)
            base.maybe_print(f"dt={dt0:.4f}s", pbar)
            base.maybe_print("", pbar)

        dapa_static = domain.H[domain.dapa_ys, domain.dapa_xs]
        base.maybe_print(
            f"[DIAG] Sample {sample_index}: slip_max={float(np.max(slip)):.4f}, domain_depth_max={float(np.max(domain.H)):.2f}m, "
            f"dapa_static_depth_max={float(np.max(dapa_static)):.2f}m, dapa_static_wet_cells={int(np.count_nonzero(dapa_static > cfg.solver.min_depth_m))}",
            pbar,
        )

        gauge_time: List[float] = []
        gauge_eta: Dict[str, List[float]] = {g.name: [] for g in cfg.gauges}
        gauge_running_max_abs: Dict[str, float] = {g.name: 0.0 for g in cfg.gauges}
        nested_times: List[float] = []
        nested_eta_boxes: List[np.ndarray] = []
        t = 0.0
        step = 0
        aborted = False
        abort_reason = ""
        land_mask = domain.z > 0.0
        gauge_support = getattr(domain, "gauge_support", {})

        while t < cfg.solver.t_end_sec:
            dt = min(dt0, cfg.solver.t_end_sec - t)
            src1 = uplift0 * base.source_time_derivative(t, rise_time)
            k1 = base.compute_rhs(q, domain, src1, cfg)
            q1 = q + dt * k1
            src2 = uplift0 * base.source_time_derivative(t + dt, rise_time)
            k2 = base.compute_rhs(q1, domain, src2, cfg)
            q = q + 0.5 * dt * (k1 + k2)

            eta = q[0]
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
                base.maybe_print(f"[ABORT] Sample {sample_index}: {abort_reason}", pbar)
                break
            if maxabseta > cfg.solver.eta_guard_m:
                aborted, abort_reason = True, f"|eta| exceeded safety guard ({maxabseta:.2f} m)"
                base.maybe_print(f"[ABORT] Sample {sample_index}: {abort_reason}", pbar)
                break

            max_eta = np.maximum(max_eta, eta.astype(np.float32))
            max_abs_eta = np.maximum(max_abs_eta, np.abs(eta).astype(np.float32))
            inun = np.where(land_mask, h_total, 0.0)
            max_depth = np.maximum(max_depth, inun.astype(np.float32))
            max_speed = np.maximum(max_speed, np.minimum(speed, 1.0e4).astype(np.float32))

            for g in cfg.gauges:
                support = gauge_support.get(g.name)
                if support is not None:
                    gval = sample_gauge_from_support(eta, support)
                else:
                    i, j = domain.gauge_idx[g.name]
                    gval = float(eta[i, j])
                if np.isfinite(gval):
                    gauge_running_max_abs[g.name] = max(gauge_running_max_abs[g.name], abs(float(gval)))

            if step % cfg.solver.snapshot_stride == 0:
                gauge_time.append(t)
                for g in cfg.gauges:
                    support = gauge_support.get(g.name)
                    if support is not None:
                        gauge_eta[g.name].append(sample_gauge_from_support(eta, support))
                    else:
                        i, j = domain.gauge_idx[g.name]
                        gauge_eta[g.name].append(float(eta[i, j]))

            if domain.nested_local_domain is not None and step % max(1, cfg.execution.nested_save_stride) == 0:
                nested_times.append(t)
                nested_eta_boxes.append(np.asarray(eta[domain.dapa_ys, domain.dapa_xs], dtype=np.float32).copy())

            if step > 0 and step % cfg.solver.diag_stride == 0:
                dapa_eta_now = eta[domain.dapa_ys, domain.dapa_xs]
                dapa_h_now = h_total[domain.dapa_ys, domain.dapa_xs]
                dapa_z_now = domain.z[domain.dapa_ys, domain.dapa_xs]
                dapa_offshore_eta = float(np.max(np.abs(dapa_eta_now[dapa_z_now <= 0.0]))) if np.any(dapa_z_now <= 0.0) else 0.0
                dapa_inun = float(np.max(np.where(dapa_z_now > 0.0, dapa_h_now, 0.0))) if dapa_h_now.size else 0.0
                base.maybe_print(
                    f"  Step {step:5d}  t={t:8.1f}s  eta=[{float(np.min(eta)):.2f},{float(np.max(eta)):.2f}]m  "
                    f"dapa_offshore_eta={dapa_offshore_eta:.2f}m  dapa_inun={dapa_inun:.2f}m",
                    pbar,
                )

            t += dt
            step += 1

        gauge_time_arr = np.asarray(gauge_time, dtype=np.float32)
        gauge_eta_arr = {k: np.asarray(v, dtype=np.float32) for k, v in gauge_eta.items()}

        gauge_max_abs: Dict[str, float] = {}
        for g in cfg.gauges:
            support = gauge_support.get(g.name)
            raw = gauge_running_max_abs.get(g.name, float("nan"))
            if support is not None and not support.valid:
                gauge_max_abs[g.name] = float("nan")
            else:
                gauge_max_abs[g.name] = float(raw) if np.isfinite(raw) else nanmax_abs(gauge_eta_arr.get(g.name, np.asarray([])))

        nominal_ys = getattr(domain, "dapa_nominal_ys", domain.dapa_ys)
        nominal_xs = getattr(domain, "dapa_nominal_xs", domain.dapa_xs)
        dapa_depth = max_depth[nominal_ys, nominal_xs]
        dapa_z_nominal = domain.z[nominal_ys, nominal_xs]
        dapa_land = dapa_z_nominal > 0.0
        positive_land = dapa_depth[(dapa_depth > 0.0) & dapa_land]
        dapa_p95_raw = float(np.percentile(positive_land, 95)) if positive_land.size else 0.0
        dapa_nonzero_cells = int(np.count_nonzero((dapa_depth > 0.0) & dapa_land))

        dapa_offshore_abs_eta = max_abs_eta[nominal_ys, nominal_xs]
        offshore_mask = dapa_z_nominal <= 0.0
        positive_offshore = dapa_offshore_abs_eta[(dapa_offshore_abs_eta > 0.0) & offshore_mask]
        dapa_offshore_p95_abs_eta = float(np.percentile(positive_offshore, 95)) if positive_offshore.size else 0.0
        dapa_offshore_nonzero_cells = int(np.count_nonzero((dapa_offshore_abs_eta > 0.0) & offshore_mask))

        dapa_window_valid = bool(getattr(domain, "dapa_window_valid", True))
        dapa_p95 = float(dapa_p95_raw) if dapa_window_valid else float("nan")

        nested_time_arr = np.asarray(nested_times, dtype=np.float32) if nested_times else None
        nested_eta_arr = np.stack(nested_eta_boxes, axis=0).astype(np.float32) if nested_eta_boxes else None

        dapa_support = gauge_support.get("Dapa_virtual")
        meta = {
            "dt0_sec": dt0,
            "rise_time_sec": rise_time,
            "dx_m": domain.dx,
            "dy_m": domain.dy,
            "t_end_sec": cfg.solver.t_end_sec,
            "wetting_drying": "moving shoreline using h=max(eta-z,0)",
            "aborted": aborted,
            "abort_reason": abort_reason,
            "dapa_virtual_valid": bool(dapa_support.valid) if dapa_support is not None else None,
            "dapa_virtual_nearest_wet_km": float(dapa_support.nearest_wet_km) if dapa_support is not None else None,
            "dapa_virtual_nearest_wet_lon": float(dapa_support.nearest_wet_lon) if dapa_support is not None else None,
            "dapa_virtual_nearest_wet_lat": float(dapa_support.nearest_wet_lat) if dapa_support is not None else None,
            "dapa_window_valid": dapa_window_valid,
            "coarse_dapa_p95_raw_m": dapa_p95_raw,
            "dapa_nonzero_cells": dapa_nonzero_cells,
            "dapa_offshore_p95_abs_eta_m": dapa_offshore_p95_abs_eta,
            "dapa_offshore_nonzero_cells": dapa_offshore_nonzero_cells,
        }

        sample = base.SampleOutput(
            slip_id=slip_id,
            gauge_time_sec=gauge_time_arr,
            gauge_eta=gauge_eta_arr,
            max_eta=max_eta,
            max_depth=max_depth,
            max_speed=max_speed,
            gauge_max_abs=gauge_max_abs,
            dapa_p95_inun=dapa_p95,
            meta=meta,
            nested_time_sec=nested_time_arr,
            nested_eta_box=nested_eta_arr,
        )
        setattr(sample, "max_abs_eta", max_abs_eta)
        setattr(sample, "dapa_p95_inun_raw", dapa_p95_raw)
        setattr(sample, "dapa_nonzero_cells", dapa_nonzero_cells)
        setattr(sample, "dapa_offshore_p95_abs_eta", dapa_offshore_p95_abs_eta)
        setattr(sample, "dapa_offshore_nonzero_cells", dapa_offshore_nonzero_cells)
        return sample

    base.prepare_domain = prepare_domain_patched
    base.run_tsunami_simulation = run_tsunami_simulation_patched


patch_base_module()


def make_progress(total: int, desc: str):
    if tqdm is not None:
        return tqdm(total=total, desc=desc, leave=True)

    class Dummy:
        def __init__(self, total: int, desc: str):
            self.total = total
            self.n = 0
            self.desc = desc

        def update(self, n: int = 1):
            self.n += n
            print(f"{self.desc}: {self.n}/{self.total}", flush=True)

        def set_postfix_str(self, s: str):
            pass

        def close(self):
            pass

    return Dummy(total, desc)


def parse_counts(text: str) -> Dict[float, int]:
    out = {8.0: 512, 8.5: 512, 9.0: 1024}
    if not text:
        return out
    for item in text.split(","):
        k, v = item.split("=")
        out[float(k.strip())] = int(v.strip())
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Runtime-focused Phase 3 batch runner. "
            "It solves many slip scenarios by magnitude, skips PNG/map rendering, "
            "skips fine local nested inundation, and saves compact machine-readable outputs for Phase 4."
        )
    )
    p.add_argument("--project-root", type=Path, default=Path("."))
    p.add_argument("--counts", type=str, default="8.0=512,8.5=512,9.0=1024")
    p.add_argument(
        "--magnitude",
        type=float,
        default=None,
        help="Optional single magnitude to run, e.g. 8.0. If omitted, runs all magnitudes in --counts.",
    )
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--dapa-only", action="store_true")
    p.add_argument("--t-end-sec", type=float, default=3600.0)
    p.add_argument(
        "--dt-max",
        type=float,
        default=0.5,
        help="Larger dt_max reduces step count. Test 0.5 first; try 1.0 only after a stable single-sample test.",
    )
    p.add_argument(
        "--downsample-factor",
        type=int,
        default=2,
        help="Use 2 for faster production screening. Use 1 for the most faithful regional solve.",
    )
    p.add_argument("--checkpoint-every", type=int, default=32)
    p.add_argument("--diag-stride", type=int, default=10**9)
    p.add_argument("--snapshot-stride", type=int, default=10**9)
    p.add_argument(
        "--save-dapa-stack",
        action="store_true",
        help="Save the coarse nominal Dapa-window max-depth stack per magnitude for Phase 4 aggregation.",
    )
    p.add_argument(
        "--save-regional-gauge-summary",
        action="store_true",
        help="Save per-sample Dapa gauge maxima CSV.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output/phase3_final_outputs"),
    )
    return p.parse_args()


def build_cfg(args: argparse.Namespace) -> base.ProjectConfig:
    ns = argparse.Namespace(
        project_root=args.project_root,
        t_end_sec=args.t_end_sec,
        dt_max=args.dt_max,
        diag_stride=args.diag_stride,
        checkpoint_every=args.checkpoint_every,
        downsample_factor=args.downsample_factor,
        save_per_sample_maps=False,
        no_paper_style_maps=True,
        no_local_dapa_render=True,
        no_local_dapa_nested=True,
        nested_save_stride=10**9,
        dapa_only=args.dapa_only,
        local_dem=None,
    )
    cfg = base.build_config(ns)
    cfg.execution.save_first_gauge_series = False
    cfg.execution.save_per_sample_maps = False
    cfg.execution.paper_style_maps = False
    cfg.execution.local_dapa_render = False
    cfg.execution.local_dapa_nested = False
    cfg.solver.snapshot_stride = max(args.snapshot_stride, 10**9)
    cfg.solver.diag_stride = max(args.diag_stride, 10**9)
    return cfg


def _init_worker(cfg: base.ProjectConfig, mw: float, return_depth: bool) -> None:
    global _WORKER_CFG, _WORKER_DOMAIN, _WORKER_MW, _WORKER_RETURN_DEPTH
    _WORKER_CFG = copy.deepcopy(cfg)
    _WORKER_MW = mw
    _WORKER_RETURN_DEPTH = return_depth

    sink = io.StringIO()
    with redirect_stdout(sink), redirect_stderr(sink):
        base.apply_magnitude_scenario(_WORKER_CFG, mw)
        _WORKER_DOMAIN = base.prepare_domain(_WORKER_CFG)


def _run_one(task: Tuple[int, np.ndarray, str]) -> Dict[str, object]:
    idx, slip, slip_name = task
    sink = io.StringIO()
    with redirect_stdout(sink), redirect_stderr(sink):
        sample = base.run_tsunami_simulation(_WORKER_CFG, _WORKER_DOMAIN, slip, slip_name, idx + 1, pbar=None)

    row = {
        "sample_index": idx,
        "slip_id": slip_name,
        "dapa_gauge_max_m": float(sample.gauge_max_abs.get("Dapa_virtual", float("nan"))),
        "coarse_dapa_p95_m": float(getattr(sample, "dapa_p95_inun", float("nan"))),
        "coarse_dapa_p95_raw_m": float(getattr(sample, "dapa_p95_inun_raw", getattr(sample, "dapa_p95_inun", float("nan")))),
        "dapa_virtual_valid": bool(sample.meta.get("dapa_virtual_valid", False)),
        "dapa_virtual_nearest_wet_km": float(sample.meta.get("dapa_virtual_nearest_wet_km", float("nan"))),
        "dapa_window_valid": bool(sample.meta.get("dapa_window_valid", False)),
        "dapa_offshore_p95_abs_eta_m": float(sample.meta.get("dapa_offshore_p95_abs_eta_m", 0.0)),
        "dapa_offshore_nonzero_cells": int(sample.meta.get("dapa_offshore_nonzero_cells", 0)),
        "aborted": bool(sample.meta.get("aborted", False)),
        "abort_reason": str(sample.meta.get("abort_reason", "")),
    }

    nominal_ys = getattr(_WORKER_DOMAIN, "dapa_nominal_ys", _WORKER_DOMAIN.dapa_ys)
    nominal_xs = getattr(_WORKER_DOMAIN, "dapa_nominal_xs", _WORKER_DOMAIN.dapa_xs)

    if _WORKER_RETURN_DEPTH:
        dapa_depth = sample.max_depth[nominal_ys, nominal_xs].astype(np.float32)
        row["dapa_nonzero_cells"] = int(getattr(sample, "dapa_nonzero_cells", np.count_nonzero(dapa_depth > 0.0)))
        row["dapa_depth"] = dapa_depth
    else:
        row["dapa_nonzero_cells"] = int(getattr(sample, "dapa_nonzero_cells", -1))
        row["dapa_depth"] = None

    return row


def save_dapa_stack(
    out_dir: Path,
    mw: float,
    domain: base.PreparedDomain,
    rows: List[Dict[str, object]],
) -> Path:
    ys = getattr(domain, "dapa_nominal_ys", domain.dapa_ys)
    xs = getattr(domain, "dapa_nominal_xs", domain.dapa_xs)
    dapa_grid = base.Grid2D(
        lon=domain.grid.lon[xs],
        lat=domain.grid.lat[ys],
        z=domain.grid.z[ys, xs],
    )
    stack = np.stack([r["dapa_depth"] for r in rows], axis=0).astype(np.float32)
    slip_ids = np.array([r["slip_id"] for r in rows], dtype=object)
    window_valid = np.array([bool(r.get("dapa_window_valid", False)) for r in rows], dtype=bool)
    gauge_valid = np.array([bool(r.get("dapa_virtual_valid", False)) for r in rows], dtype=bool)
    path = out_dir / f"Mw{mw:.1f}_coarse_dapa_max_depth_stack.npz"
    np.savez_compressed(
        path,
        lon=dapa_grid.lon,
        lat=dapa_grid.lat,
        z=dapa_grid.z.astype(np.float32),
        depth_max_stack=stack,
        slip_ids=slip_ids,
        dapa_window_valid=window_valid,
        dapa_virtual_valid=gauge_valid,
    )
    return path


def write_checkpoint(
    out_dir: Path,
    mw: float,
    rows: List[Dict[str, object]],
    save_gauge_summary: bool,
) -> Tuple[Optional[Path], Path]:
    summary_csv = out_dir / f"Mw{mw:.1f}_sample_summary_checkpoint.csv"
    if save_gauge_summary:
        df = pd.DataFrame([{k: v for k, v in r.items() if k != "dapa_depth"} for r in rows])
        df.to_csv(summary_csv, index=False)
    else:
        summary_csv = None

    gauge_vals = [finite_or_zero(r.get("dapa_gauge_max_m", float("nan"))) for r in rows]
    p95_vals = [finite_or_zero(r.get("coarse_dapa_p95_m", float("nan"))) for r in rows]
    checkpoint = {
        "mw": mw,
        "samples_completed": len(rows),
        "latest_sample_index": int(rows[-1]["sample_index"]) if rows else -1,
        "latest_slip_id": str(rows[-1]["slip_id"]) if rows else None,
        "aborted_samples": int(np.count_nonzero([bool(r["aborted"]) for r in rows])) if rows else 0,
        "valid_dapa_virtual_samples": int(np.count_nonzero([bool(r.get("dapa_virtual_valid", False)) for r in rows])) if rows else 0,
        "valid_dapa_window_samples": int(np.count_nonzero([bool(r.get("dapa_window_valid", False)) for r in rows])) if rows else 0,
        "max_dapa_gauge_max_m": float(max(gauge_vals, default=0.0)),
        "max_coarse_dapa_p95_m": float(max(p95_vals, default=0.0)),
        "sample_summary_csv": str(summary_csv) if summary_csv is not None else None,
    }
    checkpoint_json = out_dir / f"Mw{mw:.1f}_checkpoint.json"
    with open(checkpoint_json, "w", encoding="utf-8") as f:
        json.dump(checkpoint, f, indent=2)
    return summary_csv, checkpoint_json


def run_magnitude(
    cfg: base.ProjectConfig,
    mw: float,
    n_samples: int,
    workers: int,
    out_root: Path,
    save_dapa_stack_flag: bool,
    save_gauge_summary: bool,
    checkpoint_every: int,
) -> Dict[str, object]:
    base.apply_magnitude_scenario(cfg, mw)
    out_dir = out_root / f"Mw{mw:.1f}"
    set_results_dir(cfg, out_root)
    base.ensure_dir(out_dir)

    slips = base.gather_slips(cfg, mw, n_samples)
    print(f"[PHASE3] Mw{mw:.1f}: loaded {len(slips)} slip scenarios", flush=True)

    domain_ref = base.prepare_domain(cfg)
    nominal_ys = getattr(domain_ref, "dapa_nominal_ys", domain_ref.dapa_ys)
    nominal_xs = getattr(domain_ref, "dapa_nominal_xs", domain_ref.dapa_xs)
    nominal_desc = base.describe_subgrid(domain_ref.grid, domain_ref.H, nominal_ys, nominal_xs)
    effective_desc = getattr(domain_ref, "dapa_effective_info", nominal_desc)
    print(
        f"[PHASE3] Mw{mw:.1f}: regional grid={domain_ref.grid.shape[0]}x{domain_ref.grid.shape[1]} dx={domain_ref.dx:.0f}m | "
        f"Dapa nominal={nominal_desc['shape_y']}x{nominal_desc['shape_x']} wet_cells_gt_min={nominal_desc['wet_cells_gt_min']} | "
        f"Dapa effective={effective_desc['shape_y']}x{effective_desc['shape_x']} wet_cells_gt_min={effective_desc['wet_cells_gt_min']} | "
        f"window_valid={getattr(domain_ref, 'dapa_window_valid', False)}",
        flush=True,
    )
    print(f"[PHASE3] Mw{mw:.1f}: workers={workers} | samples={len(slips)} | dt_max={cfg.solver.dt_max}", flush=True)

    rows: List[Dict[str, object]] = []
    tasks = [(i, slip, slip_name) for i, (slip, slip_name) in enumerate(slips)]
    pbar = make_progress(len(tasks), f"Mw{mw:.1f}")

    aborted_count = 0
    max_gauge = 0.0
    max_p95 = 0.0
    max_offshore = 0.0

    if workers <= 1:
        _init_worker(cfg, mw, save_dapa_stack_flag)
        for done, task in enumerate(tasks, start=1):
            row = _run_one(task)
            rows.append(row)
            if row["aborted"]:
                aborted_count += 1
            max_gauge = max(max_gauge, finite_or_zero(row["dapa_gauge_max_m"]))
            max_p95 = max(max_p95, finite_or_zero(row["coarse_dapa_p95_m"]))
            max_offshore = max(max_offshore, finite_or_zero(row["dapa_offshore_p95_abs_eta_m"]))
            pbar.update(1)
            pbar.set_postfix_str(
                f"aborted={aborted_count} validGauge={int(bool(row['dapa_virtual_valid']))} "
                f"maxGauge={max_gauge:.3f}m maxP95={max_p95:.3f}m maxOffshore={max_offshore:.3f}m"
            )
            if done % max(1, checkpoint_every) == 0 or done == len(tasks):
                rows_sorted = sorted(rows, key=lambda r: int(r["sample_index"]))
                write_checkpoint(out_dir, mw, rows_sorted, save_gauge_summary)
    else:
        chunksize = max(1, min(16, len(tasks) // max(1, workers * 8)))
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=_init_worker,
            initargs=(cfg, mw, save_dapa_stack_flag),
        ) as ex:
            for done, row in enumerate(ex.map(_run_one, tasks, chunksize=chunksize), start=1):
                rows.append(row)
                if row["aborted"]:
                    aborted_count += 1
                max_gauge = max(max_gauge, finite_or_zero(row["dapa_gauge_max_m"]))
                max_p95 = max(max_p95, finite_or_zero(row["coarse_dapa_p95_m"]))
                max_offshore = max(max_offshore, finite_or_zero(row["dapa_offshore_p95_abs_eta_m"]))
                pbar.update(1)
                pbar.set_postfix_str(
                    f"aborted={aborted_count} validGauge={int(bool(row['dapa_virtual_valid']))} "
                    f"maxGauge={max_gauge:.3f}m maxP95={max_p95:.3f}m maxOffshore={max_offshore:.3f}m"
                )
                if done % max(1, checkpoint_every) == 0 or done == len(tasks):
                    rows_sorted = sorted(rows, key=lambda r: int(r["sample_index"]))
                    write_checkpoint(out_dir, mw, rows_sorted, save_gauge_summary)

    pbar.close()

    rows.sort(key=lambda r: int(r["sample_index"]))
    df = pd.DataFrame([{k: v for k, v in r.items() if k != "dapa_depth"} for r in rows])

    summary_csv = out_dir / f"Mw{mw:.1f}_sample_summary.csv"
    if save_gauge_summary:
        df.to_csv(summary_csv, index=False)

    dapa_stack_path = None
    if save_dapa_stack_flag:
        dapa_stack_path = save_dapa_stack(out_dir, mw, domain_ref, rows)

    runtime_summary = {
        "mw": mw,
        "n_samples": int(len(rows)),
        "workers": int(workers),
        "downsample_factor": int(cfg.execution.downsample_factor),
        "t_end_sec": float(cfg.solver.t_end_sec),
        "dt_max": float(cfg.solver.dt_max),
        "regional_grid_shape": [int(domain_ref.grid.shape[0]), int(domain_ref.grid.shape[1])],
        "regional_dx_m": float(domain_ref.dx),
        "regional_dy_m": float(domain_ref.dy),
        "dapa_nominal_window_shape": [int(nominal_desc["shape_y"]), int(nominal_desc["shape_x"])],
        "dapa_nominal_max_static_depth_m": float(nominal_desc["max_static_depth"]),
        "dapa_nominal_wet_cells_gt_min": int(nominal_desc["wet_cells_gt_min"]),
        "dapa_effective_window_shape": [int(effective_desc["shape_y"]), int(effective_desc["shape_x"])],
        "dapa_effective_wet_cells_gt_min": int(effective_desc["wet_cells_gt_min"]),
        "dapa_window_valid": bool(getattr(domain_ref, "dapa_window_valid", False)),
        "aborted_samples": int(np.count_nonzero(df["aborted"].values)) if not df.empty else 0,
        "valid_dapa_virtual_samples": int(np.count_nonzero(df["dapa_virtual_valid"].values)) if not df.empty else 0,
        "valid_dapa_window_samples": int(np.count_nonzero(df["dapa_window_valid"].values)) if not df.empty else 0,
        "mean_dapa_gauge_max_m": finite_or_none(df["dapa_gauge_max_m"].mean()) if not df.empty else None,
        "max_dapa_gauge_max_m": finite_or_none(df["dapa_gauge_max_m"].max()) if not df.empty else None,
        "mean_coarse_dapa_p95_m": finite_or_none(df["coarse_dapa_p95_m"].mean()) if not df.empty else None,
        "max_coarse_dapa_p95_m": finite_or_none(df["coarse_dapa_p95_m"].max()) if not df.empty else None,
        "mean_coarse_dapa_p95_raw_m": finite_or_none(df["coarse_dapa_p95_raw_m"].mean()) if not df.empty else None,
        "max_coarse_dapa_p95_raw_m": finite_or_none(df["coarse_dapa_p95_raw_m"].max()) if not df.empty else None,
        "mean_dapa_offshore_p95_abs_eta_m": finite_or_none(df["dapa_offshore_p95_abs_eta_m"].mean()) if not df.empty else None,
        "max_dapa_offshore_p95_abs_eta_m": finite_or_none(df["dapa_offshore_p95_abs_eta_m"].max()) if not df.empty else None,
        "nonzero_inundation_samples": int(np.count_nonzero(df["dapa_nonzero_cells"].values > 0)) if (save_dapa_stack_flag and not df.empty) else None,
        "sample_summary_csv": str(summary_csv) if save_gauge_summary else None,
        "coarse_dapa_stack_npz": str(dapa_stack_path) if dapa_stack_path else None,
        "notes": [
            "Phase 3 runtime-only mode: no PNG maps are generated here.",
            "Fine local nested inundation is disabled to maximize throughput.",
            "BLAS/OpenMP threads are capped at 1 per worker to reduce oversubscription.",
            "Worker initialization and per-sample logs are silenced; only startup and one progress bar are shown.",
            "Dapa gauge maxima are tracked every solver step, so they remain valid even when snapshot_stride is very large.",
            "Coarse Dapa land metrics are validity-masked when the nominal Dapa window is under-resolved or the virtual gauge is too far from Dapa.",
        ],
    }

    with open(out_dir / f"Mw{mw:.1f}_runtime_summary.json", "w", encoding="utf-8") as f:
        json.dump(runtime_summary, f, indent=2)

    return runtime_summary


def main() -> None:
    args = parse_args()
    cfg = build_cfg(args)
    out_root = args.out_dir.resolve()
    base.ensure_dir(out_root)

    counts = parse_counts(args.counts)

    if args.magnitude is not None:
        magnitudes = [float(args.magnitude)]
        if float(args.magnitude) not in counts:
            counts[float(args.magnitude)] = 512 if abs(float(args.magnitude) - 9.0) > 1e-9 else 1024
    else:
        magnitudes = [8.0, 8.5, 9.0]

    print("[PHASE3] Runtime-focused production mode", flush=True)
    print("[PHASE3] No PNG maps will be generated here.", flush=True)
    print("[PHASE3] No fine local nested inundation will be run here.", flush=True)
    print("[PHASE3] BLAS/OpenMP threads capped at 1 per worker.", flush=True)
    print("[PHASE3] Dapa honesty patch enabled: invalid coarse Dapa outputs are masked, and gauge maxima are tracked every time step.", flush=True)
    print(f"[PHASE3] workers={args.workers} | downsample_factor={args.downsample_factor} | dt_max={args.dt_max}", flush=True)
    print(f"[PHASE3] outputs -> {out_root}", flush=True)

    overall_rows: List[Dict[str, object]] = []
    for mw in magnitudes:
        summary = run_magnitude(
            cfg=copy.deepcopy(cfg),
            mw=mw,
            n_samples=counts[mw],
            workers=args.workers,
            out_root=out_root,
            save_dapa_stack_flag=args.save_dapa_stack,
            save_gauge_summary=True,
            checkpoint_every=args.checkpoint_every,
        )
        overall_rows.append(summary)
        overall_df = pd.DataFrame(overall_rows)
        overall_df.to_csv(out_root / "phase3_runtime_overall_summary.csv", index=False)
        with open(out_root / "phase3_runtime_overall_summary.json", "w", encoding="utf-8") as f:
            json.dump(overall_rows, f, indent=2)

    print(f"[DONE] Wrote runtime summaries to {out_root}", flush=True)


if __name__ == "__main__":
    main()
