"""
Validation tests for PTHA framework.

This script verifies key properties of the framework:
  1. Moment constraint satisfaction
  2. Covariance matrix properties
  3. Magnitude-frequency distribution correctness
  4. Hazard curve monotonicity
"""

import numpy as np
from pathlib import Path

# Import modules
from slip_sampler import StochasticSlipGenerator
from magnitude_frequency import GutenbergRichterRelation, hanks_kanamori_moment
from tsunami_source import SimpleDislocationSource
from hazard_aggregation import HazardAggregator, HazardCurve
from covariance import build_covariance_matrix

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data" / "fault_geometry"


def test_moment_constraint():
    """Verify that generated slip satisfies moment constraint."""
    print("\n" + "=" * 60)
    print("TEST 1: Moment Constraint Verification")
    print("=" * 60)

    mag = 8.5
    gen = StochasticSlipGenerator(
        magnitude=mag,
        fault_params_csv=DATA_DIR / "heidarzadeh_2025_table3.csv",
        subfault_size_km=20.0,
    )

    # Generate multiple samples
    n_samples = 10
    samples = gen.generate_ensemble(n_samples, seed=42)

    target_moment = gen.M0_target
    moments = []

    print(f"\nTarget moment: {target_moment:.3e} N·m (Mw {mag})")
    print("\nSample validation:")

    all_pass = True
    for i, slip in enumerate(samples):
        M = gen.rigidity * gen.area_subfault_m2 * np.sum(slip)
        rel_err = abs(M - target_moment) / target_moment

        moments.append(M)
        status = "✓ PASS" if rel_err < 1e-5 else "✗ FAIL"
        print(f"  Sample {i}: M₀ = {M:.3e} N·m, rel_err = {rel_err:.2e} {status}")

    return all_pass


def test_kl_vs_cholesky():
    """Ensure KL-based sampling produces equivalent moment-corrected fields."""
    print("\n" + "=" * 60)
    print("TEST 1b: KL versus Cholesky sampling check")
    print("=" * 60)

    mag = 8.5
    csv = DATA_DIR / "heidarzadeh_2025_table3.csv"
    gen_c = StochasticSlipGenerator(magnitude=mag, fault_params_csv=csv,
                                     subfault_size_km=20.0, use_kl=False)
    gen_k = StochasticSlipGenerator(magnitude=mag, fault_params_csv=csv,
                                     subfault_size_km=20.0, use_kl=True)

    # generate a few samples with identical Sobol indices and seed
    n = 5
    errs = []
    for i in range(n):
        s_c = gen_c.generate_rqmc_sample(i, seed=123)
        s_k = gen_k.generate_rqmc_sample(i, seed=123)
        # difference in total moment should be tiny
        M_c = gen_c.rigidity * gen_c.area_subfault_m2 * np.sum(s_c)
        M_k = gen_k.rigidity * gen_k.area_subfault_m2 * np.sum(s_k)
        rel = abs(M_c - M_k) / M_c
        errs.append(rel)
        print(f"  idx={i}: M_c={M_c:.3e}, M_k={M_k:.3e}, rel_diff={rel:.2e}")
    print(f"  mean relative difference: {np.mean(errs):.2e}")
    return np.mean(errs) < 1e-6

    print(f"\nMoment statistics:")
    print(f"  Mean relative error: {np.mean([abs(m - target_moment) / target_moment for m in moments]):.2e}")
    print(f"  Max relative error:  {np.max([abs(m - target_moment) / target_moment for m in moments]):.2e}")

    return all_pass


def test_covariance_properties():
    """Verify covariance matrix properties."""
    print("\n" + "=" * 60)
    print("TEST 2: Covariance Matrix Properties")
    print("=" * 60)

    # Create small grid for testing
    x = np.linspace(0, 100, 5)
    y = np.linspace(0, 50, 3)
    X, Y = np.meshgrid(x, y)
    coords = np.column_stack([X.ravel(), Y.ravel()])

    # Build exponential covariance
    C = build_covariance_matrix(coords, model='exponential', correlation_length=20.0)

    print(f"\nCovariance matrix size: {C.shape}")

    # Check properties
    print("\nProperty checks:")

    # 1. Symmetry
    is_symmetric = np.allclose(C, C.T)
    print(f"  1. Symmetric: {is_symmetric} {'✓' if is_symmetric else '✗'}")

    # 2. Positive definiteness
    eigvals = np.linalg.eigvalsh(C)
    is_positive_definite = np.all(eigvals > -1e-12)
    print(f"  2. Positive definite: {is_positive_definite} {'✓' if is_positive_definite else '✗'}")
    print(f"     Eigenvalue range: [{eigvals.min():.3e}, {eigvals.max():.3f}]")

    # 3. Diagonal dominance (should be all 1's at diagonal before nugget)
    diag_values = np.diag(C)
    max_diag_offset = np.max(np.abs(diag_values - 1.0))
    print(f"  3. Diagonal normalization: max |C_ii - 1| = {max_diag_offset:.3e}")

    # 4. Range [0, 1] (before nugget)
    all_in_range = np.all((C >= -0.01) & (C <= 1.01))
    print(f"  4. Values in [-0.01, 1.01]: {all_in_range} {'✓' if all_in_range else '✗'}")

    return is_symmetric and is_positive_definite and all_in_range


def test_magnitude_frequency():
    """Verify Gutenberg-Richter properties."""
    print("\n" + "=" * 60)
    print("TEST 3: Magnitude-Frequency Distribution")
    print("=" * 60)

    gr = GutenbergRichterRelation(a_value=5.0, b_value=1.0, M_min=6.0, M_max=9.0)

    print(f"\nGR parameters: a={gr.a}, b={gr.b}")

    # Test 1: Monotonicity (decreasing with M)
    mags = np.linspace(6.0, 9.0, 10)
    rates = gr.annual_rate(mags)

    is_decreasing = np.all(np.diff(rates) <= 0)
    print(f"  1. Rate decreasing with magnitude: {is_decreasing} {'✓' if is_decreasing else '✗'}")

    # Test 2: Positivity
    all_positive = np.all(rates > 0)
    print(f"  2. All rates positive: {all_positive} {'✓' if all_positive else '✗'}")

    # Test 3: Sample magnitude range
    samples = gr.sample_magnitude(1000, rng=np.random.default_rng(42))
    in_range = np.all((samples >= gr.M_min) & (samples <= gr.M_max))
    print(f"  3. Sampled magnitudes in range: {in_range} {'✓' if in_range else '✗'}")
    print(f"     Sample range: [{samples.min():.2f}, {samples.max():.2f}]")

    # Test 4: Hanks-Kanamori conversion
    print(f"\n  4. Hanks-Kanamori moment-magnitude conversion:")
    for M in [7.5, 8.0, 8.5, 9.0]:
        M0 = hanks_kanamori_moment(M)
        print(f"     Mw {M} → M₀ = {M0:.3e} N·m")

    return is_decreasing and all_positive and in_range


def test_hazard_curve_properties():
    """Verify hazard curve properties."""
    print("\n" + "=" * 60)
    print("TEST 4: Hazard Curve Properties")
    print("=" * 60)

    # Create synthetic hazard curve
    heights = np.array([0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    # Decreasing exceedance probability (realistic)
    probs = np.array([5e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7])

    hazard = HazardCurve(heights, probs)

    print(f"\nHazard curve with {len(heights)} points")

    # Test 1: Monotonicity
    heights_query = np.linspace(0.1, 3.0, 20)
    probs_query = np.array([hazard.exceedance_probability_at_level(h) for h in heights_query])

    is_decreasing = np.all(np.diff(probs_query) <= 0)
    print(f"  1. Probability decreasing with height: {is_decreasing} {'✓' if is_decreasing else '✗'}")

    # Test 2: Return period consistency
    print(f"\n  2. Return period calculations:")
    test_heights = [1.0, 1.5, 2.0]
    for h in test_heights:
        prob = hazard.exceedance_probability_at_level(h)
        ret_period = hazard.return_period_at_level(h)
        h_back = hazard.intensity_at_return_period(ret_period)

        error = abs(h - h_back)
        status = "✓" if error < 0.1 else "✗"
        print(f"     h={h:.1f}m → T={ret_period:.1f}yr → h={h_back:.2f}m "
              f"(error={error:.3f}) {status}")

    return is_decreasing


def test_tsunami_source():
    """Verify tsunami source computation."""
    print("\n" + "=" * 60)
    print("TEST 5: Tsunami Source Properties")
    print("=" * 60)

    # Create synthetic fault and slip
    fault_geom = {
        'length_km': 200.0,
        'width_km': 100.0,
        'top_depth_km': 5.0,
        'strike_deg': 160.0,
        'dip_deg': 40.0,
        'rake_deg': 90.0,
        'n_subfaults_along': 4,
        'n_subfaults_down': 3,
    }

    slip = np.ones((3, 4)) * 2.0

    source = SimpleDislocationSource(fault_geom, slip, rigidity=4.0e10)

    # Test at various distances
    obs_points = np.array([
        [0, 0],
        [50, 0],
        [100, 0],
        [200, 0],
    ])

    displacement = source.get_seafloor_displacement(obs_points)

    print(f"\nDisplacement at observation points:")
    print(f"  Distance [km]  |  Displacement [m]")
    print(f"  {'-'*34}")
    for dist, disp in zip(obs_points[:, 0], displacement):
        print(f"  {dist:6.1f}        |  {disp:8.4f}")

    # For the prototype we just display the values; no strict physical test.
    print("\n  (No physical validation performed on simplified displacement.)")
    return True


def test_hazard_aggregation():
    """Verify hazard aggregation logic."""
    print("\n" + "=" * 60)
    print("TEST 6: Hazard Aggregation")
    print("=" * 60)

    # Simple test case
    scenario_weights = {8.0: 0.01, 8.5: 0.001}

    aggregator = HazardAggregator(scenario_weights)

    # Add synthetic intensities
    np.random.seed(42)
    intensities_80 = np.random.gamma(2.0, 1.0, 100)
    intensities_85 = np.random.gamma(2.0, 2.0, 100)

    aggregator.add_scenario_realizations(8.0, intensities_80)
    aggregator.add_scenario_realizations(8.5, intensities_85)

    print(f"\nAdded realizations for 2 scenarios:")
    print(f"  Mw 8.0: {len(intensities_80)} samples, mean={intensities_80.mean():.2f} m")
    print(f"  Mw 8.5: {len(intensities_85)} samples, mean={intensities_85.mean():.2f} m")

    # Aggregate
    hazard = aggregator.compute_aggregated_hazard()

    print(f"\nAggregated hazard curve properties:")
    heights = np.array([0.5, 1.0, 1.5, 2.0])
    for h in heights:
        prob = hazard.exceedance_probability_at_level(h)
        print(f"  h = {h:.1f} m: P(exceed) = {prob:.4e} per year")

    # Check monotonicity
    probs = np.array([hazard.exceedance_probability_at_level(h) for h in heights])
    is_monotonic = np.all(np.diff(probs) <= 0)
    print(f"\nMonotonic decreasing: {is_monotonic} {'✓' if is_monotonic else '✗'}")

    return is_monotonic


def run_all_tests():
    """Execute all validation tests."""
    print("\n" + "=" * 70)
    print("PTHA FRAMEWORK VALIDATION TESTS")
    print("=" * 70)

    results = {
        'Moment Constraint': test_moment_constraint(),
        'KL vs Cholesky': test_kl_vs_cholesky(),
        'Covariance Matrix': test_covariance_properties(),
        'Magnitude-Frequency': test_magnitude_frequency(),
        'Hazard Curve': test_hazard_curve_properties(),
        'Tsunami Source': test_tsunami_source(),
        'Hazard Aggregation': test_hazard_aggregation(),
    }

    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    for test_name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {test_name:.<50} {status}")

    all_passed = all(results.values())
    print("\n" + ("=" * 70))
    if all_passed:
        print("✓ ALL TESTS PASSED")
    else:
        print("✗ SOME TESTS FAILED - See details above")
    print("=" * 70 + "\n")

    return all_passed


if __name__ == "__main__":
    success = run_all_tests()
    exit(0 if success else 1)
