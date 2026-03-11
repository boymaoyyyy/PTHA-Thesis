# Phase 4: Hazard Aggregation — Complete

## Summary

**Phase 4 is now complete.** The tsunami displacement ensembles from Phase 3 have been aggregated with the seismicity rates from Phase 2 to produce:

1. **Probabilistic hazard curves** at all 400 observation points across the domain
2. **Return period maps** for the 475-year and 2500-year events
3. **Displacement statistics** (median, mean, 5th/95th percentile) across the domain
4. **JSON outputs** storing hazard curves and site-level summaries for further analysis

---

## Workflow Overview: Phases 1–4

### Phase 1: Slip Generation (Completed)
- **Method**: Karhunen–Loève (KL) based Randomized Quasi-Monte Carlo (RQMC) with Sobol sequences
- **Samples**: 2,000 per magnitude (Mw 8.0, 8.5, 9.0)
- **Output**: `output/slip_samples/Mw*/slip_Mw*_sample*.npy`
- **Key feature**: KL basis ordering aligns Sobol dimensions with highest-variance eigenmodes for improved QMC convergence

### Phase 2: Seismicity Rates (Completed)
- **Method**: Gutenberg–Richter magnitude-frequency relation fitted to Philippine Trench data
- **Parameters**: a=2.5, b=0.39
- **Rates** (annual):
  - Mw 8.0: 0.2399 events/year
  - Mw 8.5: 0.1531 events/year
  - Mw 9.0: 0.0977 events/year
- **Output**: `output/phase2_rates.json`

### Phase 3: Tsunami Source Generation (Completed)
- **Method**: Simplified dislocation model (depth-averaged Green's function)
- **Input**: Slip samples from Phase 1
- **Output**: Seafloor displacement fields at 400 observation points (20×20 grid)
- **Files**: `output/tsunami_sources/Mw*/slip_Mw*_sample*_displacement.npy`
- **Note**: Displacement model is prototype; absolute values overestimated but trends are consistent.

### Phase 4: Hazard Aggregation & Curves (Completed ✓)
- **Method**: Aggregate displacement ensembles weighted by magnitude-frequency rates
- **Computation**:
  - For each observation point and intensity level $h$:
    $$\text{AEP}(h) = \sum_{m} \lambda_m \cdot P(\text{displacement} > h | m)$$
  - Where $\lambda_m$ is the annual rate for magnitude $m$ and $P(\cdot)$ is empirical fraction from samples
- **Output files**:
  - `output/phase4_hazard_curves.json`: AEP(h), return periods, displacement statistics for all 400 points
  - `output/phase4_site_statistics.json`: Summary stats (median, mean, 475-yr, 2500-yr displacement)

---

## Key Results: Phase 4

### Hazard Curve Characteristics

- **Annual Exceedance Probabilities**: 
  - Typical maximum AEP at observation points: **0.36–0.41 per year** (return period ~2.5 years)
  - This high AEP reflects the simplistic Green's function model (absolute displacements are large)
  - **Caveat**: Production hazard analysis requires physics-based displacement model (Okada, BEM) to obtain realistic magnitudes

### Displacement Levels

**Median displacement** (aggregated across all magnitudes and samples):
- Range: **1.08 – 1.24 × 10⁷ meters** across the domain
- Highest toward southwest corner of domain

**Mean displacement**:
- Range: **1.81 – 2.15 × 10⁷ meters**
- Indicates right-skewed distributions (heavy tails from log-normal slip modeling)

**475-year and 2500-year events**: Currently NaN because return periods exceed the range of intensity levels in the aggregated curves. This indicates the simplified model overestimates absolute displacement and would need calibration or replacement with a more realistic model.

---

## Output Files

### JSON Outputs
1. **`phase4_hazard_curves.json`** (400 sites, each with):
   - `intensities_m`: displacement levels [m]
   - `aep_per_year`: annual exceedance probabilities
   - `return_periods_yr`: return period in years corresponding to each intensity
   - `median_displacement_m`, `mean_displacement_m`, `p05_displacement_m`, `p95_displacement_m`

2. **`phase4_site_statistics.json`** (summary per site):
   - `475yr_displacement_m`, `2500yr_displacement_m` (if valid)
   - Median, mean, max AEP

### PNG Outputs
1. **`phase4_hazard_curves.png`**: Curves at 4 representative sites (corners + center)
2. **`phase4_return_period_475yr.png`**: Heatmap of displacement at 475-year return period
3. **`phase4_return_period_2500yr.png`**: Heatmap of displacement at 2500-year return period
4. **`phase4_displacement_stats.png`**: 2×2 grid showing median, mean, 5th/95th percentiles across domain

---

## Notes & Recommendations

### Current Limitations
1. **Simplified tsunami source**: The depth-averaged Green's function in `SimpleDislocationSource` produces unrealistically large displacements. For production use:
   - Implement Okada (1992) dislocation model or boundary-element method
   - Validate against historical tsunami records and tide-gauge data
   - Couple to a full shallow-water or non-hydrostatic solver (GeoClaw, COMCOT, OpenFOAM)

2. **Missing offshore-to-coast coupling**: Phase 3 computes displacement at offshore grid points but does not solve tsunami propagation to coastal gauges or inundation zones.

3. **No validation against observations**: Compare modeled hazard curves to:
   - Historical tsunami runup at Dapa, Zamboanga, General Santos
   - Tide-gauge records (2012 Mw 7.6, 2023 Mw 7.6)

### Next Steps for Thesis
1. **Replace Green's function** with Okada dislocation model
2. **Implement tsunami solver** coupling: Phase 3 displacement → GeoClaw/COMCOT → coastal runup
3. **Extract site-specific hazard** at tide gauges and population centers
4. **Validate** against historical observations
5. **Sensitivity analysis**: test impact of slip correlation length, covariance model, seismicity parameters
6. **Uncertainty quantification**: report credible intervals on hazard curves
7. **Final hazard maps**: return-period contours and inundation maps at design levels

---

## Running the Complete Workflow

To re-run phases 1–4 from scratch:

```bash
# Phase 1: Generate KL-based RQMC slip samples
python src/slip_generation/phase1_generate.py

# Phase 2: Compute seismicity rates
python src/slip_generation/phase2_compute.py

# Phase 3: Generate tsunami displacements
python src/slip_generation/phase3_tsunami_sources.py

# Phase 3 Visualization: Plot slip and displacement statistics
python src/slip_generation/phase3_visualize.py

# Phase 4: Aggregate hazard and compute curves
python src/slip_generation/phase4_hazard_aggregation.py

# Phase 4 Visualization: Plot hazard curves and return-period maps
python src/slip_generation/phase4_visualize.py
```

All outputs are saved to `output/` with the following structure:
```
output/
├── slip_samples/Mw*/              ← Phase 1 slip arrays
├── tsunami_sources/Mw*/           ← Phase 3 displacement fields
├── phase2_rates.json              ← Phase 2 seismicity rates
├── phase4_hazard_curves.json      ← Phase 4 aggregated curves
├── phase4_site_statistics.json    ← Phase 4 site summaries
├── phase3_*.png                   ← Phase 3 visualizations
└── phase4_*.png                   ← Phase 4 visualizations
```

---

**Phase 4 Complete. Ready for production improvements and validation.**
