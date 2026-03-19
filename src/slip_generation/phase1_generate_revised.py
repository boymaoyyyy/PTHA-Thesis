#!/usr/bin/env python3
"""
Revised Phase 1 stochastic slip generator for the Philippine Trench PTHA workflow.

Design goals
------------
- Keep the revised thesis framing: RQMC + Matérn spatial correlation + strict moment conservation.
- Generate slips on the SAME subfault grid that the revised Phase 3 runner expects by default
  (15 along-strike x 5 down-dip = 75 subfaults), avoiding the Phase 1/Phase 3 grid mismatch.
- Reject pathological realizations using transparent quality-control thresholds.
- Save outputs to output/RQMC_samples_revised/
- Emit a per-sample statistics catalog for rapid QC before running expensive Phase 3 simulations.

Examples
--------
# Sensible first pass for the revised workflow (64 per magnitude)
py phase1_generate_revised.py --project-root . --counts 8.0=64,8.5=64,9.0=64

# Quick test with six Mw 9.0 slips only
py phase1_generate_revised.py --project-root . --magnitudes 9.0 --counts 9.0=6

# Stricter QC and overwrite prior output
py phase1_generate_revised.py --project-root . --counts 8.0=64,8.5=64,9.0=128 --overwrite \
    --max-slip-cap-factor 8.0 --p99-to-mean-cap 5.5
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from scipy import special, spatial, stats

try:
    from scipy.stats import qmc
except Exception:
    qmc = None

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:
    plt = None


@dataclass(frozen=True)
class Scenario:
    mw: float
    rupture_length_km: float
    rupture_width_km: float
    top_depth_km: float
    strike_deg: float
    dip_deg: float
    rake_deg: float
    mean_slip_m: float


SCENARIOS: Dict[float, Scenario] = {
    8.0: Scenario(8.0, 200.0, 100.0, 7.6, 160.0, 40.0, 90.0, 2.7),
    8.5: Scenario(8.5, 300.0, 100.0, 7.6, 164.0, 39.0, 90.0, 4.7),
    9.0: Scenario(9.0, 600.0, 100.0, 6.1, 158.0, 33.0, 90.0, 8.1),
}


@dataclass
class Phase1Config:
    project_root: Path
    output_dirname: str = "RQMC_samples_revised"
    seed: int = 42
    n_along: int = 15
    n_down: int = 5
    corr_length_km: float = 20.0
    hurst_nu: float = 0.3
    shear_modulus_pa: float = 4.0e10
    sigma_ln: float = 0.55
    max_attempt_factor: int = 30
    moment_tolerance: float = 1.0e-6
    max_slip_cap_factor: float = 8.0
    max_to_mean_cap: float = 8.0
    p99_to_mean_cap: float = 5.5
    save_preview_png: bool = True

    @property
    def output_root(self) -> Path:
        return self.project_root / "output" / self.output_dirname


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def target_seismic_moment(mw: float) -> float:
    return 10.0 ** (1.5 * mw + 9.1)


def scenario_subfault_centers_km(scenario: Scenario, n_along: int, n_down: int) -> np.ndarray:
    along = np.linspace(0.5, n_along - 0.5, n_along) * (scenario.rupture_length_km / n_along)
    down = np.linspace(0.5, n_down - 0.5, n_down) * (scenario.rupture_width_km / n_down)
    xx, yy = np.meshgrid(along, down)
    return np.column_stack([xx.ravel(), yy.ravel()])


def matern_covariance(distance_km: np.ndarray, corr_length_km: float, nu: float) -> np.ndarray:
    r = np.asarray(distance_km, dtype=float)
    scaled = np.sqrt(2.0 * nu) * np.maximum(r, 1.0e-12) / max(corr_length_km, 1.0e-12)
    coeff = (2.0 ** (1.0 - nu)) / special.gamma(nu)
    out = coeff * (scaled ** nu) * special.kv(nu, scaled)
    out[r == 0.0] = 1.0
    out[np.isnan(out)] = 1.0
    return out


def parse_counts(text: str | None) -> Dict[float, int]:
    default = {8.0: 64, 8.5: 64, 9.0: 64}
    if not text:
        return default
    out = dict(default)
    for item in text.split(","):
        if not item.strip():
            continue
        k, v = item.split("=")
        out[float(k.strip())] = int(v.strip())
    return out


def parse_magnitudes(text: str | None) -> List[float]:
    if not text:
        return [8.0, 8.5, 9.0]
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def subfault_area_m2(scenario: Scenario, n_along: int, n_down: int) -> float:
    return (scenario.rupture_length_km * 1000.0 / n_along) * (scenario.rupture_width_km * 1000.0 / n_down)


def realized_moment(slip_m: np.ndarray, scenario: Scenario, cfg: Phase1Config) -> float:
    area = subfault_area_m2(scenario, cfg.n_along, cfg.n_down)
    return cfg.shear_modulus_pa * area * float(np.sum(slip_m))


def generate_raw_correlated_field(engine, chol: np.ndarray, n_dim: int) -> np.ndarray:
    if qmc is not None and engine is not None:
        u = np.clip(engine.random(1)[0], 1.0e-12, 1.0 - 1.0e-12)
        z = stats.norm.ppf(u)
    else:
        z = np.random.default_rng().standard_normal(n_dim)
    return chol @ z


def build_engines(cfg: Phase1Config, dims: int, n_needed: int) -> object:
    if qmc is None:
        return None
    # Request enough points for likely acceptance/rejection losses.
    max_draws = max(n_needed * cfg.max_attempt_factor, 64)
    return qmc.Sobol(d=dims, scramble=True, seed=cfg.seed)


def sample_one_slip(scenario: Scenario, cfg: Phase1Config, chol: np.ndarray, engine) -> np.ndarray:
    dims = cfg.n_along * cfg.n_down
    corr_norm = generate_raw_correlated_field(engine, chol, dims)

    # Lognormal field with controlled variance.  The unscaled mean is tied to the
    # scenario mean slip so that moment correction is modest rather than extreme.
    mu_ln = math.log(max(scenario.mean_slip_m, 1.0e-6)) - 0.5 * (cfg.sigma_ln ** 2)
    raw = np.exp(mu_ln + cfg.sigma_ln * corr_norm).reshape(cfg.n_down, cfg.n_along)

    # Global moment-conserving scaling.
    m0_target = target_seismic_moment(scenario.mw)
    m0_raw = realized_moment(raw, scenario, cfg)
    scale = m0_target / max(m0_raw, 1.0e-18)
    slip = raw * scale
    return slip.astype(np.float64)


def slip_stats(slip: np.ndarray, scenario: Scenario, cfg: Phase1Config) -> Dict[str, float]:
    flat = np.asarray(slip, dtype=float).ravel()
    mean_slip = float(np.mean(flat))
    max_slip = float(np.max(flat))
    p95 = float(np.percentile(flat, 95))
    p99 = float(np.percentile(flat, 99))
    m0 = realized_moment(slip, scenario, cfg)
    m0_target = target_seismic_moment(scenario.mw)
    rel_err = abs(m0_target - m0) / m0_target
    return {
        "min_slip_m": float(np.min(flat)),
        "mean_slip_m": mean_slip,
        "median_slip_m": float(np.median(flat)),
        "max_slip_m": max_slip,
        "std_slip_m": float(np.std(flat)),
        "p95_slip_m": p95,
        "p99_slip_m": p99,
        "sum_slip_m": float(np.sum(flat)),
        "max_to_mean": max_slip / max(mean_slip, 1.0e-12),
        "p99_to_mean": p99 / max(mean_slip, 1.0e-12),
        "target_moment_Nm": m0_target,
        "realized_moment_Nm": m0,
        "moment_rel_error": rel_err,
    }


def accept_slip(stats_row: Dict[str, float], scenario: Scenario, cfg: Phase1Config) -> Tuple[bool, str]:
    if stats_row["moment_rel_error"] > cfg.moment_tolerance:
        return False, f"moment error {stats_row['moment_rel_error']:.3e} exceeds tolerance"
    if stats_row["max_slip_m"] > cfg.max_slip_cap_factor * scenario.mean_slip_m:
        return False, (
            f"max slip {stats_row['max_slip_m']:.2f}m exceeds "
            f"{cfg.max_slip_cap_factor:.1f}x scenario mean slip"
        )
    if stats_row["max_to_mean"] > cfg.max_to_mean_cap:
        return False, f"max/mean ratio {stats_row['max_to_mean']:.2f} exceeds cap"
    if stats_row["p99_to_mean"] > cfg.p99_to_mean_cap:
        return False, f"p99/mean ratio {stats_row['p99_to_mean']:.2f} exceeds cap"
    if stats_row["mean_slip_m"] <= 0.0:
        return False, "mean slip is non-positive"
    return True, "accepted"


def save_preview_png(out_path: Path, slip: np.ndarray, title: str) -> None:
    if plt is None:
        return
    fig, ax = plt.subplots(figsize=(7, 3.5), dpi=180)
    im = ax.imshow(slip, origin="lower", aspect="auto")
    ax.set_title(title)
    ax.set_xlabel("Along strike index")
    ax.set_ylabel("Down dip index")
    fig.colorbar(im, ax=ax, shrink=0.85, label="Slip (m)")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def generate_for_scenario(scenario: Scenario, n_keep: int, cfg: Phase1Config) -> Tuple[pd.DataFrame, Dict[str, int]]:
    centers = scenario_subfault_centers_km(scenario, cfg.n_along, cfg.n_down)
    dist = spatial.distance_matrix(centers, centers)
    cov = matern_covariance(dist, cfg.corr_length_km, cfg.hurst_nu) + 1.0e-8 * np.eye(dist.shape[0])
    chol = np.linalg.cholesky(cov)
    engine = build_engines(cfg, dims=centers.shape[0], n_needed=n_keep)

    scenario_dir = cfg.output_root / f"Mw{scenario.mw:.1f}"
    ensure_dir(scenario_dir)

    rows: List[Dict[str, float]] = []
    accepted = 0
    attempted = 0
    rejected = 0
    max_attempts = max(n_keep * cfg.max_attempt_factor, n_keep)

    while accepted < n_keep and attempted < max_attempts:
        attempted += 1
        slip = sample_one_slip(scenario, cfg, chol, engine)
        stats_row = slip_stats(slip, scenario, cfg)
        ok, reason = accept_slip(stats_row, scenario, cfg)

        stats_row.update({
            "mw": scenario.mw,
            "attempt_index": attempted - 1,
            "accepted": bool(ok),
            "decision": reason,
            "rupture_length_km": scenario.rupture_length_km,
            "rupture_width_km": scenario.rupture_width_km,
            "top_depth_km": scenario.top_depth_km,
            "strike_deg": scenario.strike_deg,
            "dip_deg": scenario.dip_deg,
            "rake_deg": scenario.rake_deg,
            "reference_mean_slip_m": scenario.mean_slip_m,
            "n_along": cfg.n_along,
            "n_down": cfg.n_down,
            "corr_length_km": cfg.corr_length_km,
            "hurst_nu": cfg.hurst_nu,
            "sigma_ln": cfg.sigma_ln,
        })

        if ok:
            out_name = f"slip_Mw{scenario.mw:.1f}_sample{accepted:04d}.npy"
            np.save(scenario_dir / out_name, slip.astype(np.float32))
            stats_row["sample_index"] = accepted
            stats_row["file_name"] = out_name
            if cfg.save_preview_png and accepted < 3:
                save_preview_png(
                    scenario_dir / f"slip_Mw{scenario.mw:.1f}_sample{accepted:04d}.png",
                    slip,
                    f"Mw {scenario.mw:.1f} slip sample {accepted}",
                )
            accepted += 1
        else:
            rejected += 1
            stats_row["sample_index"] = -1
            stats_row["file_name"] = ""

        rows.append(stats_row)

    df = pd.DataFrame(rows)
    df.to_csv(scenario_dir / f"Mw{scenario.mw:.1f}_generation_log.csv", index=False)

    if accepted < n_keep:
        raise RuntimeError(
            f"Could not obtain {n_keep} accepted samples for Mw {scenario.mw:.1f}. "
            f"Accepted {accepted} out of {attempted} attempts. Relax QC or inspect generator settings."
        )

    accepted_df = df[df["accepted"] == True].copy()  # noqa: E712
    accepted_df.to_csv(scenario_dir / f"Mw{scenario.mw:.1f}_accepted_stats.csv", index=False)

    summary = {
        "mw": scenario.mw,
        "accepted_samples": int(accepted),
        "attempted_samples": int(attempted),
        "rejected_samples": int(rejected),
        "acceptance_rate": float(accepted / max(attempted, 1)),
        "reference_mean_slip_m": scenario.mean_slip_m,
        "accepted_mean_of_means_m": float(accepted_df["mean_slip_m"].mean()),
        "accepted_mean_of_maxima_m": float(accepted_df["max_slip_m"].mean()),
        "accepted_max_of_maxima_m": float(accepted_df["max_slip_m"].max()),
        "accepted_mean_p99_to_mean": float(accepted_df["p99_to_mean"].mean()),
        "accepted_mean_max_to_mean": float(accepted_df["max_to_mean"].mean()),
        "moment_tolerance": cfg.moment_tolerance,
        "max_slip_cap_factor": cfg.max_slip_cap_factor,
        "max_to_mean_cap": cfg.max_to_mean_cap,
        "p99_to_mean_cap": cfg.p99_to_mean_cap,
    }
    with open(scenario_dir / f"Mw{scenario.mw:.1f}_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    return accepted_df, {
        "accepted": accepted,
        "attempted": attempted,
        "rejected": rejected,
    }


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Revised Phase 1 RQMC slip generator with QC and Phase 3-compatible grids")
    p.add_argument("--project-root", type=Path, default=Path("."))
    p.add_argument("--magnitudes", type=str, default="8.0,8.5,9.0", help="Comma-separated magnitudes to generate")
    p.add_argument("--counts", type=str, default="8.0=64,8.5=64,9.0=64", help="Comma-separated Mw=count pairs")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output-dirname", type=str, default="RQMC_samples_revised")
    p.add_argument("--n-along", type=int, default=15)
    p.add_argument("--n-down", type=int, default=5)
    p.add_argument("--corr-length-km", type=float, default=20.0)
    p.add_argument("--hurst-nu", type=float, default=0.3)
    p.add_argument("--sigma-ln", type=float, default=0.55, help="Lognormal sigma in natural-log space")
    p.add_argument("--max-attempt-factor", type=int, default=30)
    p.add_argument("--moment-tolerance", type=float, default=1.0e-6)
    p.add_argument("--max-slip-cap-factor", type=float, default=8.0)
    p.add_argument("--max-to-mean-cap", type=float, default=8.0)
    p.add_argument("--p99-to-mean-cap", type=float, default=5.5)
    p.add_argument("--no-preview-png", action="store_true")
    p.add_argument("--overwrite", action="store_true", help="Delete output/RQMC_samples_revised before regenerating")
    return p


def main() -> None:
    args = build_parser().parse_args()
    counts = parse_counts(args.counts)
    magnitudes = parse_magnitudes(args.magnitudes)

    cfg = Phase1Config(
        project_root=args.project_root,
        output_dirname=args.output_dirname,
        seed=args.seed,
        n_along=args.n_along,
        n_down=args.n_down,
        corr_length_km=args.corr_length_km,
        hurst_nu=args.hurst_nu,
        sigma_ln=args.sigma_ln,
        max_attempt_factor=args.max_attempt_factor,
        moment_tolerance=args.moment_tolerance,
        max_slip_cap_factor=args.max_slip_cap_factor,
        max_to_mean_cap=args.max_to_mean_cap,
        p99_to_mean_cap=args.p99_to_mean_cap,
        save_preview_png=not args.no_preview_png,
    )

    if args.overwrite and cfg.output_root.exists():
        shutil.rmtree(cfg.output_root)
    ensure_dir(cfg.output_root)

    manifest_rows = []
    print("Phase 1 revised: RQMC slip generation")
    print("====================================")
    print(f"Output root: {cfg.output_root}")
    print(f"Grid: {cfg.n_along} along x {cfg.n_down} down = {cfg.n_along * cfg.n_down} subfaults")
    print(f"Matérn correlation: ell={cfg.corr_length_km:.1f} km, nu={cfg.hurst_nu:.2f}")
    print(f"QC caps: max slip <= {cfg.max_slip_cap_factor:.1f}x scenario mean, max/mean <= {cfg.max_to_mean_cap:.1f}, p99/mean <= {cfg.p99_to_mean_cap:.1f}")
    print()

    for mw in magnitudes:
        if mw not in SCENARIOS:
            raise ValueError(f"Unsupported scenario magnitude: {mw}")
        scenario = SCENARIOS[mw]
        n_keep = counts.get(mw, 64)
        print(f"Generating Mw {mw:.1f}: {n_keep} accepted slips")
        accepted_df, tally = generate_for_scenario(scenario, n_keep, cfg)
        print(
            f"  accepted={tally['accepted']} attempted={tally['attempted']} rejected={tally['rejected']} "
            f"max(max_slip)={accepted_df['max_slip_m'].max():.2f}m mean(mean_slip)={accepted_df['mean_slip_m'].mean():.2f}m"
        )
        manifest_rows.append({
            "mw": mw,
            "accepted": tally["accepted"],
            "attempted": tally["attempted"],
            "rejected": tally["rejected"],
            "acceptance_rate": tally["accepted"] / max(tally["attempted"], 1),
            "max_of_max_slip_m": float(accepted_df["max_slip_m"].max()),
            "mean_of_mean_slip_m": float(accepted_df["mean_slip_m"].mean()),
            "output_dir": str(cfg.output_root / f"Mw{mw:.1f}"),
        })

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(cfg.output_root / "generation_manifest.csv", index=False)
    with open(cfg.output_root / "phase1_config_resolved.json", "w", encoding="utf-8") as f:
        json.dump(asdict(cfg), f, indent=2, default=str)

    print()
    print("Done.")
    print(f"Manifest written to: {cfg.output_root / 'generation_manifest.csv'}")


if __name__ == "__main__":
    main()
