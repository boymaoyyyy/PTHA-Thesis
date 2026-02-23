# Philippine Trench PTHA Parameters - Thesis Extraction Summary

## Document Overview
This document summarizes all domain-specific parameters extracted from the thesis for the Philippine Trench Probabilistic Tsunami Hazard Assessment (PTHA) framework.

**Thesis Reference:**
- Title: DEVELOPMENT OF A PHYSICAL FIDELITY-BASED PROBABILISTIC TSUNAMI HAZARD ASSESSMENT FRAMEWORK UTILIZING RANDOMIZED QUASI MONTE CARLO FOR SLIP DISTRIBUTION SAMPLING
- Authors: Henrich Miguel Del Rio Carpio, Carl Vince Callao Dominguez, John Mar Maagad Estimada, Karl Andre Lopez Gutierrez
- Date: December 2025

**Data Sources:**
- Heidarzadeh et al. (2025) fault geometry models (Tables 2-3)
- PHIVOLCS earthquake catalog (1900-2024)
- Regional seismic studies (Bautista 2019)
- Global slip inversion studies
- Bathymetric and exposure data

---

## PHASE 1: STOCHASTIC SLIP GENERATION PARAMETERS

### 1.1 Slip Random Field Model (Section 3.3.2)
**Source:** Section 3.3.2 "Slip Random Field Model" + Section 3.3.1 "RQMC method"

Earthquake slip is modeled as a **log-normal Gaussian random field**:
```
log(s) = μₛ + L·ξ
```

| Parameter | Value | Reference | Notes |
|-----------|-------|-----------|-------|
| **Covariance Model** | Matérn | Global slip inversions | Enables spatial correlation with adjustable roughness |
| **Correlation Length** | 20 km | Global subduction zone standards | Controls asperity structure spacing |
| **Hurst Exponent** | 0.3 | Global slip inversions | Roughness/smoothness parameter (0-1) |
| **N_subfaults** | ~400 | Section 3.3.1 | Discretization of fault plane |
| **Slip Distribution** | Log-normal | Section 3.3.2 | Ensures non-negative slip values |

### 1.2 RQMC Sampling Configuration (Section 3.3.1, 3.3.3)
**Source:** Section 3.3.1 "RQMC method" + Section 3.3.3 "RQMC Sampling"

Randomized Quasi-Monte Carlo (RQMC) sampling using Sobol sequences with Owen scrambling:

| Parameter | Value | Reference | Justification |
|-----------|-------|-----------|----------------|
| **Sampling Method** | RQMC Sobol' + Owen | Section 3.3.1 | Superior space-filling in 400-D parameter space |
| **Sequence Type** | Sobol' low-discrepancy | Section 3.3.1 | Convergence rate O(N⁻¹) vs O(N⁻¹/²) for MC |
| **Randomization** | Owen scrambling | Section 3.3.1 | Enables statistical inference (confidence intervals) |
| **Samples/Magnitude** | 2,000-5,000 | Section 3.3.3 | Ensures robust statistical representation |
| **Variance Reduction** | 1-2 orders of magnitude | Section 3.3.1 | vs. standard Monte Carlo |

### 1.3 Moment Constraint (Section 3.3.4)
**Source:** Section 3.3.4 "Seismic Moment Calculation and Generate-and-Scale Correction"

Each slip realization undergoes **generate-and-scale** correction to enforce target moment:

| Parameter | Formula | Constraint |
|-----------|---------|-----------|
| **Seismic Moment (Raw)** | M₀^raw = Σᵢ μᵢ·aᵢ·sᵢ | From raw slip distribution |
| **Target Moment** | M₀ = 10^(1.5·Mw + 4.4) | Hanks-Kanamori relation |
| **Scaling Factor** | α = M₀ / M₀^raw | Applied uniformly to slip field |
| **Corrected Slip** | s = α·s_raw | Final slip for simulation |
| **Moment Error Tolerance** | \|M₀ - M₀_realized\| / M₀ < 10⁻⁶ | Validation threshold |

---

## PHASE 2: GUTENBERG-RICHTER SEISMICITY PARAMETERS

### 2.1 Recurrence Model Selection (Section 3.4.1)
**Source:** Section 3.4.1 "Selection of the Gutenberg-Richter Recurrence Relation"

Classical **Gutenberg-Richter (G-R) law** adopted for earthquake frequency distribution:

**Exponential Form (Continuous):**
```
v(m) = β·exp(-β·m)
```

**Traditional Log-Linear Form:**
```
log₁₀(N) = a - b·M
```

**Relationship:** β = b·ln(10) ≈ 2.303·b

| Reason for Selection | Reference |
|-------------------|-----------|
| Validated globally across diverse tectonic settings (8+ decades) | Section 3.4.1 |
| Recommended in UNESCO-IOC international tsunami hazard guidelines | Section 3.4.1 |
| Addresses fundamental challenge of balancing frequent vs. rare events | Section 3.4.1 |
| Directly integrates multiple magnitude scenarios into single hazard estimate | Section 3.4.1 |

### 2.2 Seismicity Parameter Estimation (Section 3.4.2)
**Source:** Section 3.4.2 "Seismicity Parameter Estimation for the Philippine Trench"

#### 2.2.1 Direct Catalog Analysis
| Analysis Method | Result | Source | Applicability |
|---|---|---|---|
| **Maximum Likelihood Estimation (MLE)** | b ≈ 0.95 ± 0.12 | PHIVOLCS catalog 1900-2024 | Magnitude range Mw 5.0-7.5 |
| **Data Range** | Mw 5.0-7.5 events | Instrumental records | Well-constrained, complete |

#### 2.2.2 Adopted Parameters (Conservative Approach)
**Source:** Section 3.4.2 - Multi-criteria decision process

| Parameter | Value | Reasoning |
|-----------|-------|-----------|
| **β (Exponential Form)** | **0.9** | Conservative choice for megathrust capability |
| **b (Log-linear Form)** | **0.39** | Corresponds to β = 0.9 |
| **Justification Source 1** | Regional Seismic Precedent | Holocene marine terraces (Davao Oriental) show M≥8.0 events (8,080-4,140 cal yr BP) |
| **Justification Source 2** | Comparative Subduction Zones | Nankai (β≈0.85-0.95), Mediterranean (β≈0.90-1.05), Java (β≈0.88-0.98) |
| **Justification Source 3** | Sensitivity Analysis | β<0.85 produces hazard dominated by Mw 8.5-9.0; β≥0.95 dominated by moderate events |

#### 2.2.3 Physical Interpretation
- β=0.9 implies **"flat" magnitude-frequency distribution**
- Large earthquakes not exponentially suppressed compared to moderate events
- Appropriate for subduction zone with demonstrated great-earthquake potential but limited modern records

### 2.3 Magnitude Range
| Parameter | Value | Notes |
|-----------|-------|-------|
| **Magnitude Minimum** | Mw 5.0 | Catalog completeness threshold |
| **Magnitude Maximum** | Mw 9.0 | Physical upper bound (Holocene evidence) |
| **Modeling Scenarios** | Mw 7.6, 8.0, 8.5, 9.0 | Spectrum from common to rare |

---

## PHASE 3: FAULT GEOMETRY PARAMETERS

### 3.1 Historical Ruptures (Heidarzadeh et al. 2025, Table 2)
**Source:** Section 3.2.2 "Collection of Fault Geometry and Seismotectonic Parameters"

#### 3.1.1 August 31, 2012 - Mw 7.6 Earthquake

| Parameter | Value | Unit |
|-----------|-------|------|
| Event Date | August 31, 2012 | - |
| Magnitude | 7.6 | Mw |
| Fault Length (L) | 128 | km |
| Fault Width (W) | 90 | km |
| Top Depth | 1 | km |
| Strike | 345.8 | ° |
| Dip | 45.8 | ° |
| Rake | 61.2 | ° |
| Maximum Slip | 3.2 | m |
| Mean Slip | 0.4 | m |

**Validation Data:**
- Tide gauge amplitude: 3.7 cm (max observed)
- Wave periods: 8.0-18.3 minutes

#### 3.1.2 December 2, 2023 - Mw 7.6 Earthquake

| Parameter | Value | Unit |
|-----------|-------|------|
| Event Date | December 2, 2023 | - |
| Magnitude | 7.6 | Mw |
| Fault Length (L) | 144 | km |
| Fault Width (W) | 120 | km |
| Top Depth | 14.5 | km |
| Strike | 167 | ° |
| Dip | 17 | ° |
| Rake | 62.6 | ° |
| Maximum Slip | 1.8 | m |
| Mean Slip | 0.3 | m |

**Validation Data:**
- Tide gauge amplitude: 12.5 cm (max observed)
- Wave periods: 6.7-28.2 minutes

### 3.2 Hypothetical Megathrust Scenarios (Heidarzadeh et al. 2025, Table 3)
**Source:** Section 3.2.2 "Collection of Fault Geometry and Seismotectonic Parameters"

#### 3.2.1 Mw 8.0 Hypothetical Scenario

| Parameter | Value | Unit |
|-----------|-------|------|
| Scenario Name | Hypothetical Mw 8.0 | - |
| Magnitude | 8.0 | Mw |
| Fault Length (L) | 200 | km |
| Fault Width (W) | 100 | km |
| Top Depth | 7.6 | km |
| Mean Strike | 160 | ° |
| Mean Dip | 40 | ° |
| Rake | 90 | ° |
| Mean Slip | 2.7 | m |

#### 3.2.2 Mw 8.5 Hypothetical Scenario

| Parameter | Value | Unit |
|-----------|-------|------|
| Scenario Name | Hypothetical Mw 8.5 | - |
| Magnitude | 8.5 | Mw |
| Fault Length (L) | 300 | km |
| Fault Width (W) | 100 | km |
| Top Depth | 7.6 | km |
| Mean Strike | 164 | ° |
| Mean Dip | 39 | ° |
| Rake | 90 | ° |
| Mean Slip | 4.7 | m |

**Expected Impact (from thesis intro):**
- Building inundation in Dapa: ~218 buildings

#### 3.2.3 Mw 9.0 Hypothetical Scenario

| Parameter | Value | Unit |
|-----------|-------|------|
| Scenario Name | Hypothetical Mw 9.0 | - |
| Magnitude | 9.0 | Mw |
| Fault Length (L) | 600 | km |
| Fault Width (W) | 100 | km |
| Top Depth | 6.1 | km |
| Mean Strike | 158 | ° |
| Mean Dip | 33 | ° |
| Rake | 90 | ° |
| Mean Slip | 8.1 | m |

**Expected Impact (from thesis intro):**
- Coastal waves: ~17.4 m maximum run-up

### 3.3 Geophysical Parameters

| Parameter | Value | Source | Notes |
|-----------|-------|--------|-------|
| **Depth-dependent Shear Modulus (μ)** | Depth-interpolated | PREM + regional models | Section 3.2.2: "converted seismic velocity profiles" |
| **Fault Geometry Source** | Heidarzadeh et al. (2025) | Seismic tomography, hypocenter distributions, geodetic constraints | Section 3.2.2: "subduction interface model" |
| **Study Area (Latitude)** | 7°-12°N | Philippine Trench segment | Central to Northern region |

---

## PHASE 3: TSUNAMI SIMULATION PARAMETERS

### 4.1 Hydrodynamic Solver Configuration (Section 3.1, Phase 3)
**Source:** Section 3.1 Phase 3: "High-Fidelity Tsunami Simulation"

| Component | Specification | Reference | Justification |
|-----------|---------------|-----------|----------------|
| **Solver Type** | Nonlinear shallow-water equations | Section 3.1 Ph3 | Standard for tsunami propagation |
| **Reconstruction Scheme** | 3rd-order WENO (WENO3) | Section 3.1 Ph3 | Captures steep wave fronts, avoids oscillations |
| **Time Integration** | RK3 (3rd-order Runge-Kutta) | Implicit in WENO scheme | Consistent accuracy order |

### 4.2 Adaptive Mesh Refinement (AMR) (Section 3.1, Phase 3)
**Source:** Section 3.1 Phase 3: "Adaptive Mesh Refinement (AMR)"

| Parameter | Value | Purpose |
|-----------|-------|---------|
| **Use AMR** | Yes | Dynamically allocate computational resolution |
| **Refinement Levels** | 3 (minimum) | Hierarchy depth |
| **Deep Ocean Resolution** | 1 km (1000 m) | Efficient wave propagation |
| **Coastal Resolution (Dapa)** | 10 m (or finer) | Urban-scale detail for inundation |
| **Refine Criterion** | Gradient of water elevation | Automatic shock/discontinuity detection |

### 4.3 Domain Configuration (Section 3.1, Phase 3)
**Source:** Study area: Central to Northern Philippine Trench (7°-12°N), focus on Dapa

| Parameter | Range | Type |
|-----------|-------|------|
| **Longitude Domain** | ~125.5°-130.5°E | Approximate |
| **Latitude Domain** | ~6°-13°N | Approximate |
| **Bathymetry Source** | GEBCO (1 km resolution) | High-resolution global grid |
| **Bathymetry Inclusion** | Yes | Critical for wave shoaling |

### 4.4 Physical Enhancements and Validations (Section 2.4, 3.1)
**Source:** Section 2.4 "Limitations and Validation of Methods" + Section 3.1 "Dual-Validation Strategy"

| Feature | Status | Reference | Importance |
|---------|--------|-----------|-----------|
| **Fully Coupled Approach (Method 1)** | Used | Section 2.4.2 | Captures non-hydrostatic filtering, acoustic-gravity wave coupling |
| **Dynamic Seafloor Forcing** | Used | Section 3.1 Ph3 | Time-dependent boundary condition (vs. static initial conditions) |
| **Non-hydrostatic Filtering** | Noted, Not Required | Section 2.4.1 | Future enhancement; deep-ocean model adequate |
| **Acoustic-Gravity Wave Coupling** | Noted, Not Required | Section 2.4.6 | Future enhancement for near-field details |
| **Poroelastic Effects** | Not required | Section 2.4.6 | Future extension for trench regions |

### 4.5 Validation Data (Section 3.1, "Dual-Validation Strategy")
**Source:** Section 3.1 paragraph on validation

**Physical Validation Against Observations:**
- **2012 Mw 7.6 Event:** Tide gauge time series from IOC Sea Level Monitoring Facility
  - Peak amplitude: 3.7 cm
  - Wave periods: 8.0-18.3 min
- **2023 Mw 7.6 Event:** Tide gauge time series from IOC + PHIVOLCS
  - Peak amplitude: 12.5 cm
  - Wave periods: 6.7-28.2 min

**Metrics Compared:**
- Arrival time (should match observations)
- Amplitude (should match peak heights)
- Waveform shape (should match oscillation pattern)

---

## PHASE 4: HAZARD AGGREGATION PARAMETERS

### 5.1 Output Metrics (Section 3.1, Phase 4)
**Source:** Section 3.1 Phase 4: "Hazard Aggregation"

| Output Product | Specification | Computation Method |
|---|---|---|
| **Primary Metric** | Inundation depth (meters) | Maximum water level above ground |
| **Hazard Curves** | Site-specific | AEP vs. inundation depth at given locations |
| **Hazard Maps** | Spatial | Annual Exceedance Probability maps for return periods |
| **Return Period Maps** | 100-yr, 500-yr, 2500-yr | Standard intervals for risk management |
| **Exposure Metrics** | Buildings inundated per year | Expected annual count (Dapa focus) |

### 5.2 Hazard Aggregation Formula (Section 3.1, Phase 4)
**Source:** Section 3.1 Phase 4

For each location and inundation depth threshold h:
```
AEP(h) = Σᵢ P(rupture i occurs) × P(simulated inundation i > h)
       = Σᵢ λᵢ(Mw) × I(depth_i > h)
```

Where:
- λᵢ(Mw) = annual occurrence rate from Gutenberg-Richter model
- I() = indicator function (1 if depth exceeds threshold, 0 otherwise)
- Summation over all magnitude scenarios and realizations

### 5.3 Output Thresholds and Return Periods
**Source:** Section 3.1 Phase 4 + implicit in hazard analysis

**Inundation Depth Thresholds:**
```
[0.5, 1.0, 2.0, 5.0, 10.0, 15.0] meters
```

**Return Periods for Hazard Maps:**
```
[100, 500, 2500] years
```

**Exposure Metrics (Dapa Focus):**
- Expected annual inundated buildings
- OpenStreetMap infrastructure data
- Comparison benchmark from thesis: ~218 buildings (Mw 8.5)

---

## SUMMARY OF KEY PARAMETERS FOR FRAMEWORK IMPLEMENTATION

### Framework Configuration Class Attributes

```python
PhilippineTrenchSlipParameters:
  ├─ covariance_model = "matern"
  ├─ correlation_length_km = 20.0
  ├─ hurst_exponent = 0.3
  ├─ sampling_method = "rqmc_sobol_owen"
  ├─ n_subfaults = 400
  ├─ n_samples_per_magnitude = (2000, 5000)
  └─ moment_error_tolerance = 1e-6

PhilippineTrenchSeismicityParameters:
  ├─ beta = 0.9
  ├─ b_value = 0.39
  ├─ magnitude_minimum = 5.0
  └─ magnitude_maximum = 9.0

PhilippineTrenchTsunamiParameters:
  ├─ solver_type = "shallow_water_weno"
  ├─ weno_order = 3
  ├─ use_amr = True
  ├─ amr_levels = 3
  ├─ base_resolution_m = 1000
  └─ coastal_resolution_m = 10

PhilippineTrenchHazardParameters:
  ├─ intensity_metric = "inundation_depth_m"
  ├─ intensity_thresholds = [0.5, 1.0, 2.0, 5.0, 10.0, 15.0]
  └─ return_periods = [100, 500, 2500]

RUPTURE_SCENARIOS Dictionary:
  ├─ 2012_Mw7.6 (128×90 km, 0.4 m mean slip)
  ├─ 2023_Mw7.6 (144×120 km, 0.3 m mean slip)
  ├─ Mw8.0 (200×100 km, 2.7 m mean slip)
  ├─ Mw8.5 (300×100 km, 4.7 m mean slip)
  └─ Mw9.0 (600×100 km, 8.1 m mean slip)
```

---

## THESIS VALIDATION BENCHMARKS

**From thesis introduction and results sections:**

| Scenario | Validation Metric | Expected Value | Location |
|----------|-------------------|-----------------|----------|
| **2012 Mw 7.6** | Tide gauge amplitude | 3.7 cm (max) | Tide gauge station(s) |
| **2023 Mw 7.6** | Tide gauge amplitude | 12.5 cm (max) | Tide gauge station(s) |
| **2023 Mw 7.6** | Wave periods | 6.7-28.2 min | Tide gauge records |
| **Mw 8.5** | Building inundation | ~218 buildings | Dapa |
| **Mw 9.0** | Max coastal waves | ~17.4 m | Philippine Trench vicinity |

---

## REFERENCES CITED IN THESIS

**Key References Used for Parameters:**

1. **Heidarzadeh et al. (2025)** - Fault geometry models (Tables 2-3)
2. **Bautista (2019)** - Regional seismicity analysis
3. **Abrahams et al. (2023)** - Fully coupled vs. superposition methods
4. **PHIVOLCS** - Philippine earthquake catalog (1900-2024)
5. **PREM** - Preliminary Reference Earth Model (shear modulus depth profiles)
6. **Ramos et al. (2012)** - Holocene marine terraces evidence of past megathrust events
7. **Nankai/Mediterranean/Java Trench studies** - Comparative subduction zone β-values

---

## CONFIGURATION FILES CREATED

1. **philippine_trench_config.py** (src/slip_generation/)
   - Implements all parameters as Python dataclasses
   - Includes RuptureScenario definitions for all 5 magnitude scenarios
   - Provides PhilippineTrenchPTHAConfiguration for integrated workflow
   - Includes summary() method for verification

---

*Document generated from thesis analysis: December 2025*
*All parameters extracted directly from thesis sections 3.1-3.4, Heidarzadeh et al. (2025) tables, and PHIVOLCS catalog*
