# PAGE 1

 
 
DEVELOPMENT OF A PHYSICAL FIDELITY-BASED PROBABILISTIC 
TSUNAMI HAZARD ASSESSMENT FRAMEWORK UTILIZING 
RANDOMIZED QUASI MONTE CARLO FOR  
SLIP DISTRIBUTION SAMPLING 
 
 
 
 
 
 
SUBMITTED TO THE FACULTY OF THE  
DEPARTMENT OF COMPUTER SCIENCE  
COLLEGE OF INFORMATION TECHNOLOGY AND COMPUTING  
BY  
 
 
 
 
 
HENRICH MIGUEL DEL RIO CARPIO  
CARL VINCE CALLAO DOMINGUEZ  
JOHN MAR MAAGAD  ESTIMADA  
KARL ANDRE LOPEZ GUTIERREZ  
 
 
 
 
 
IN PARTIAL FULFILLMENT OF THE 
REQUIREMENTS FOR THE 
DEGREE OF 
 
 
 
 
 
BACHELOR OF SCIENCE IN COMPUTER SCIENCE 
 
 
 
 
 
 
 
 
 
DECEMBER 2025 
 


---

# PAGE 2

 
TABLE OF CONTENTS 
 
 
 
Chapter I....................................................................................................................... 1 
 
INTRODUCTION........................................................................................................1 
 
1.2   Statement of the Problem............................................................................. 4 
1.2.1 Slow Convergence on Earthquake Slip Distributions.......................5 
1.2.2 Uncertainty in Propagation Accuracy.................................................5 
1.3 Research Objectives.........................................................................................6 
1.4 Scope and Limitations.................................................................................... 7 
1.5 Significance of the Study................................................................................ 8 
1.6 Definition of Terms..........................................................................................9 
1.7 Proposed Probabilistic Tsunami Hazard Assessment Framework........ 10 
 
Chapter II....................................................................................................................11 
 
REVIEW OF RELATED LITERATURE...................................................................11 
 
2.1 Tsunami Hazards and Hydrodynamic Complexity in the Philippine 
Trench.................................................................................................................... 11 
2.1.1 Recent Tsunamigenic Events (2012 vs. 2023)....................................13 
2.1.2 Vulnerability of Coastal Infrastructure..............................................14 
2.2 Probabilistic Tsunami Hazard Assessment (PTHA).................................15 
2.2.1 Treatment of Uncertainty.................................................................... 17 
2.2.2 Integration with Seismic Hazard Models......................................... 18 
2.2.3 Case Studies by Region........................................................................18 
2.2.4 Advances in PTHA Modeling............................................................ 21 
2.3 Gutenberg–Richter Law: Earthquake Frequency Distribution...............22 
2.3.1 Application to Tsunami Hazard Analysis.........................................23 
2.3.2 Tsunami Size–Frequency Relations....................................................25 
2.3.3 Shallow vs. Deep Earthquakes in Tsunami Generation..................26 
2.3.4 Critiques and Limitations of Gutenberg–Richter in Tsunami 
Context............................................................................................................27 
2.4 The Physics of Tsunami Generation: Coupled vs. Uncoupled Models. 28 
2.4.1 Limitations of Superposition (Method 4)..........................................29 
2.4.2 Validation of the Fully Coupled Method (Method 1)..................... 30 
2.4.3 The Evolution to Dynamic Rupture-Tsunami Coupling................ 30 
2.4.4 The Evolution to Dynamic Rupture-Tsunami Coupling................ 31 
 


---

# PAGE 3

 
2.4.5 Application in Probabilistic Tsunami Hazard Analysis (PTHA)...32 
2.4.6 Capturing Transient Wave Physics: Seismic and Acoustic 
Coupling......................................................................................................... 32 
2.5 Probabilistic Tsunami Hazard Assessment (PTHA) Paradigms.............33 
2.6 Computational Efficiency in Hazard Modeling........................................34 
2.6.1 The Curse of Dimensionality and Monte Carlo Limitations..........34 
2.6.2 Advanced Variance Reduction Strategies.........................................35 
2.7 Randomized Quasi-Monte Carlo (RQMC) as a Methodological Solution
35 
2.7.1 Mathematical Rigor and Statistical Guarantee.................................36 
2.7.2 Superior Convergence in High-Dimensional Spaces......................36 
 
Chapter III...................................................................................................................38 
 
METHODOLOGY..................................................................................................... 38 
 
3.1 Overall Methodology....................................................................................38 
3.2 Preliminaries.................................................................................................. 41 
3.2.1 Methodological Foundations..............................................................41 
3.2.2 Collection of Fault Geometry and Seismotectonic Parameters......42 
3.3 Phase 1: Stochastic Slip Generation and Distribution sampling............ 43 
3.3.1 RQMC method......................................................................................43 
3.3.2 Slip Random Field Model....................................................................45 
3.3.3 RQMC Sampling...................................................................................46 
3.3.4 Seismic Moment Calculation and Generate-and-Scale Correction...
46 
3.4 Phase 2: Probability Model...........................................................................48 
3.4.1 Selection of the Gutenberg-Richter Recurrence Relation................48 
3.4.2 Seismicity Parameter Estimation for the Philippine Trench.......... 50 
3.4.3 Integration into Hazard Aggregation................................................52 
3.5 Phase 3: High-Fidelity, Fully Coupled Tsunami Simulation...................53 
3.5.1 Numerical Model: Modified GeoClaw..............................................53 
3.5.2 Source Implementation........................................................................54 
3.6 Dual Validation Framework........................................................................ 55 
3.6.1 Statistical validation — convergence efficiency............................... 55 
3.6.2 Physical validation — waveform fidelity..........................................55 
3.7 Hazard Aggregation..................................................................................... 56 
3.7.1 Annual Probability of Exceedance.....................................................56 
3.7.2  Outputs and Evaluation..................................................................... 58 
3.8 Experimental Procedure...............................................................................58 
 


---

# PAGE 4

 
3.8.1 Simulation Campaign Structure.........................................................58 
3.8.2 Computational Workflow................................................................... 59 
3.8.3 Trial Execution and Replication......................................................... 60 
3.8.4 Experimental Setup and Computing Environment.........................60 
3.8.5 Software Implementation....................................................................62 
3.9 Evaluation Metrics........................................................................................ 63 
3.10 Initial Simulation......................................................................................... 63 
 
References:..................................................................................................................74 
 
 


---

# PAGE 5

1 
Chapter I 
INTRODUCTION 
 
 
 
 
The Philippine Trench, a north–south trending subduction zone 
extending over 1,300 kilometers along the eastern margin of the Philippine 
archipelago stands as one of the most seismically active yet critically 
understudied tsunamigenic sources in the western Pacific. Formed by the 
westward subduction of the Philippine Sea Plate beneath the Philippine 
Mobile Belt at rates exceeding 10 cm/year, this trench has repeatedly 
demonstrated its capacity to generate destructive earthquakes and tsunamis. 
Despite its clear hazard potential, systematic tsunami risk assessments for the 
Philippine Trench have been scarce, leaving densely populated coastal 
communities in the eastern Philippines dangerously unprepared. This 
knowledge gap has become increasingly untenable in light of recent seismic 
activity: two Mw 7.6 earthquakes, one in August 2012 and another in 
December 2023, both originated along this interface and produced measurable 
tsunamis with maximum tide gauge amplitudes of 3.7 cm and 12.5 cm, 
respectively (Ocean Engineering, 2025). Notably, the 2023 tsunami exhibited 
longer wave periods (6.7–28.2 minutes) compared to the 2012 event (8.0–18.3 
minutes), a difference attributed to the deeper bathymetry surrounding the 
2023 rupture zone, which allowed longer-wavelength energy to propagate 
more efficiently toward the coast (Ocean Engineering, 2025). These real-world 
events provide critical validation data while simultaneously underscoring the 
urgent need for advanced, probabilistic modeling frameworks that move 
beyond deterministic “what-if” scenarios toward quantifiable risk estimates 
grounded in physical realism and uncertainty quantification. 
 
 
At the heart of tsunami generation lies the spatial distribution of 
coseismic slip along the fault plane, a primary source of uncertainty in hazard 
modeling. The initial sea surface displacement, which serves as the boundary 
condition for hydrodynamic simulations, is governed by the vertical seafloor 
deformation resulting from slip. Using the elastic dislocation theory of Okada 
(1985), this displacement at any surface point x = (x, y) can be expressed as 
 
 
                                                            
 
where Δuz
(i) denotes the vertical displacement due to the i-th subfault 
and θi includes geometric and physical parameters most critically, slip 
magnitude si. Because slip cannot be observed prior to an earthquake and 
exhibits strong spatial heterogeneity due to frictional and stress variations 
 


---

# PAGE 6

2 
along the plate interface, it must be treated as a stochastic field. A common 
and physically motivated approach models log-slip as a Gaussian process: 
 
 
 
 
 
where μs is the mean log-slip, C(⋅;ℓ) is a covariance function (e.g., 
exponential or Matérn), and ℓ is the correlation length that controls the spatial 
smoothness of slip patches (Goda & Song, 2016). This formulation captures 
the physical reality that large slip tends to cluster in asperities, while 
surrounding regions experience little or no displacement. The total seismic 
moment M0​   provides a crucial physical constraint that couples all slip values: 
 
 
 
 
 
 
where Ω is the rupture area,  is the spatial mean slip, A=∣Ω∣ , and μ is 
𝑆
the effective shear modulus. For a target magnitude Mw, M0 is fixed via the 
moment–magnitude relation: Mw = ⅔  log10 (M0) −10.7 (Hanks & Kanamori, 
1979). Thus, while slip heterogeneity is permitted, it must collectively satisfy 
this global moment budget, a nonlinear constraint that renders the sampling 
problem high-dimensional and nontrivial. 
 
 
To explore this uncertainty space efficiently, traditional Monte Carlo 
(MC) methods generate independent slip realizations and compute 
corresponding tsunami metrics, but their slow convergence rate (O(N−1/2) ) 
makes them computationally prohibitive for expensive hydrodynamic 
simulations. Quasi-Monte Carlo (QMC) methods, which use low-discrepancy 
sequences such as Sobol’ (Sobol’, 1967), achieve faster convergence often close 
to O(N−1) in moderate dimensions but lack built-in variance estimation, 
limiting their utility for risk-informed decisions. This limitation is overcome 
by Randomized Quasi-Monte Carlo (RQMC), which applies random 
scrambling to a QMC sequence to produce R independent replicates, enabling  
unbiased estimation and robust error quantification (L’Ecuyer & Lemieux, 
2002). The RQMC estimator for the expected maximum tsunami height at 
location x becomes 
 
 
 
    
 
 
with sample variance that supports confidence intervals essential for 
actionable risk metrics. This framework enables statistically rigorous, 
 


---

# PAGE 7

3 
computationally efficient Probabilistic Tsunami Hazard Analysis (PTHA), a 
capability urgently needed for the Philippine Trench (Ocean Engineering, 
2025). 
 
 
To accurately simulate tsunami generation and propagation, this study 
adopts 
Method 
1: 
the 
fully-coupled 
modeling 
approach, 
which 
simultaneously solves earthquake rupture dynamics, seismic/acoustic wave 
propagation through the solid Earth and ocean, and tsunami hydrodynamics 
across the basin and nearshore zones. This method ensures physical 
consistency by dynamically linking the evolving slip distribution s (r,t) to 
both seismic radiation and seafloor deformation, thereby avoiding the 
artificial decoupling inherent in traditional two-step models (e.g., static 
Okada displacement followed by shallow-water propagation). The system is 
governed by coupled partial differential equations, including the elastic wave 
equation for crustal motion and the nonlinear shallow-water equations for 
tsunami evolution (LeVeque et al., 2011). While this approach offers superior 
physical fidelity, its predictive accuracy particularly in coastal inundation 
depends critically on the numerical treatment of the propagation phase. 
 
 
Recognizing this, this study proposes a key enhancement to 
Fully-coupled Method: the implementation of high-order, adaptive, and 
well-balanced numerical schemes to improve propagation accuracy. Standard 
implementations often rely on second-order finite volume methods on 
uniform grids, which can introduce numerical dispersion and dissipation 
especially in shallow coastal zones with complex bathymetry. To address this, 
this study integrates third-order Weighted Essentially Non-Oscillatory 
(WENO) reconstruction to preserve wavefront sharpness, bathymetry-aware 
adaptive mesh refinement (AMR) to dynamically increase resolution near 
shorelines and steep gradients, and well-balanced fluxes to maintain the 
lake-at-rest condition over irregular topography (George, 2008; LeVeque et al., 
2011). These enhancements ensure that small but critical features such as 
run-up height, inundation extent, and arrival time are captured with greater 
fidelity, directly impacting risk estimates for vulnerable communities. 
 
 
This improved framework will be applied specifically to the Philippine 
Trench, building directly on the scenarios outlined in the 2025 study. That 
work modeled both the 2012 and 2023 Mw 7.6 events and expanded to 
hypothetical Mw 8.0, 8.5, and 9.0 ruptures, revealing that a worst-case Mw 9.0 
earthquake could generate coastal tsunami waves up to 17.4 meters high 
comparable to the 2011 Tohoku and 2004 Indian Ocean tsunamis (Ocean 
Engineering, 2025). Even an Mw 8.5 event could inundate 218 buildings in the 
coastal town of Dapa, including a designated evacuation center, thereby 
compromising the very infrastructure meant to ensure public safety (Ocean 
Engineering, 2025). This study extends this deterministic analysis by 
embedding it within an RQMC-based probabilistic framework that samples 
 


---

# PAGE 8

4 
thousands of physically consistent slip distributions, propagates each through 
the enhanced fully-coupled model, and quantifies exceedance probabilities for 
tsunami height, inundation depth, and building exposure. 
 
 
The annual probability that tsunami height exceeds a threshold h at 
location x is given by 
 
 
 
 
 
where 
ν(m)=β 
exp 
(−βm) 
is 
the 
Gutenberg–Richter 
frequency–magnitude distribution, and the conditional probability is 
estimated via RQMC using an indicator function over the ensemble of 
simulations (González et al., 2009). This formulation produces tsunami hazard 
curves and inundation probability maps that can directly inform land-use 
planning, early warning protocols, and evacuation strategies in the 
Philippines, a nation where over 60% of the population lives in coastal zones 
and institutional capacity for tsunami preparedness remains limited. 
 
 
In conclusion, the Philippine Trench is no longer a hidden hazard but a 
quantifiable, 
high-consequence 
threat 
demanding advanced scientific 
attention. By integrating stochastic slip physics, efficient RQMC sampling, 
and a propagation-enhanced fully-coupled tsunami model, this thesis 
responds directly to the urgent call from the 2025 study for comprehensive 
tsunami risk assessment along this underexplored subduction zone (Ocean 
Engineering, 
2025). 
The 
resulting 
framework 
not 
only 
advances 
methodological rigor but also delivers actionable insights for building true 
resilience in one of the world’s most disaster-prone archipelagos. 
 
 
1.2   Statement of the Problem 
 
 
The Philippine Trench remains significantly understudied in terms of 
tsunami hazard projection, leaving substantial gaps in understanding its 
potential impacts (Heidarzadeh et al., 2025). While existing modeling efforts 
(e.g., Reyes et al., 2020) provide a useful foundation, closer examination 
reveals 
persistent 
limitations: 
uncertainties 
in 
earthquake 
source 
characterization (Qiu et al., 2019) and trade-offs in numerical methods, such 
as simplified wave‐propagation models used to maintain computational 
efficiency (FUNMAP study, 2025). These technical shortcomings can hinder 
realistic simulation, particularly in capturing complex wave behavior and 
achieving numerical convergence. Therefore, this study aims to refine and 
 


---

# PAGE 9

5 
enhance current models to overcome these limitations and provide more 
reliable, physically realistic simulations tailored to the Philippine Trench. 
 
​
1.2.1 Slow Convergence on Earthquake Slip Distributions 
 
Probabilistic tsunami hazard assessment (PTHA) requires the 
generation of a large ensemble of stochastic earthquake slip 
distributions to account for epistemic uncertainties in seismic source 
modeling. These slip realizations are used to simulate tsunami wave 
heights at coastal locations, forming the basis for risk-informed 
planning and mitigation. Traditionally, the Monte Carlo (MC) method 
has been employed to generate these stochastic slip scenarios due to its 
simplicity and general applicability. However, MC is known to suffer 
from slow convergence rates, specifically O(N⁻¹ᐟ²), which significantly 
increases the computational cost required to achieve stable and reliable 
hazard estimates (Hou et al., 2019). 
 
As 
demonstrated 
in 
convergence 
analyses 
of 
tsunami 
simulations, thousands to tens of thousands of samples may be needed 
for each magnitude bin to stabilize maximum coastal tsunami height 
distributions (e.g., via the coefficient of variation). This demand on 
sample size makes MC-based modeling inefficient, especially in 
high-resolution, fully coupled tsunami models. Therefore, there is a 
need to explore alternative sampling methods that can achieve more 
accurate results with fewer simulations.​
 
1.2.2 Uncertainty in Propagation Accuracy  
 
Maintaining Physical Fidelity Against Simplifications: Simpler 
tsunami models (which rely on shortcuts and approximations) fail 
when simulating wave travel in deep water because they neglect 
phenomena like wave dispersion (waves changing speed and shape as 
they travel) and non-hydrostatic filtering (how deep water naturally 
filters out small-scale wave features created by the seafloor movement) 
(Abrahams et al., 2023).When these effects are ignored, predictions 
suffer major errors: the simpler models can "overpredict the tsunami 
wave amplitude by a factor of two" (Abrahams et al., 2023) or 
incorrectly predict a much larger initial wave that "arrives too early" 
(Abrahams et al., 2023). This study must ensure advanced simulation 
which correctly includes these dispersive and filtering effects maintains 
this accuracy as the wave propagates over distance. 
 


---

# PAGE 10

6 
1.3 Research Objectives 
 
This study aims to develop a high-fidelity probabilistic tsunami hazard 
assessment framework specifically for the Philippine Trench. The objectives 
are as follows: 
 
1.​ To implement a Randomized Quasi-Monte Carlo (RQMC) sampling 
model that efficiently generates synthetic earthquake slip distributions. 
 
2.​ To enhance the propagation accuracy of the Fully Coupled Method 
(Method 1) explicitly selected over the Superposition Method (Method 
4) to strictly resolve non-hydrostatic filtering effects that otherwise 
cause 
twofold 
wave amplitude overprediction by integrating 
high-order numerical schemes, and to apply this framework together 
with the RQMC sampling model to the Philippine Trench. 
 
3.​ To design and implement a functional prototype that integrates RQMC 
sampling into the Enhanced Fully-Coupled Method, structured to meet 
established PTHA hazard evaluation frameworks. 
 
4.​ To evaluate the effectiveness of the Enhanced Fully Coupled Method 
with RQMC sampling using statistical and numerical standards 
grounded in NTHMP tsunami validation benchmarks and  established 
PTHA hazard evaluation frameworks. 
 
 
 
 
 
 
 
 


---

# PAGE 11

7 
1.4 Scope and Limitations 
 
Scope 
 
This study focuses exclusively on the Philippine Trench as the primary 
tsunamigenic source, modeling various stochastic rupture scenarios to assess 
hazards along the eastern seaboard of the Philippines. The geographic area of 
interest specifically targets high-risk coastal zones in the Caraga and Davao 
regions, which are geographically exposed to near-field tsunamis generated 
by the subduction interface. The modeling framework integrates a 
Randomized Quasi-Monte Carlo (RQMC) sampling technique to efficiently 
generate synthetic earthquake slip distributions, addressing the need for 
robust uncertainty quantification in probabilistic assessments. These 
stochastic sources drive a fully-coupled tsunami generation model (Method 
1), which simultaneously solves the elastic and hydrodynamic equations to 
capture dynamic rupture effects often missed by static models. To ensure 
physical fidelity during the propagation phase, the study employs enhanced 
numerical 
schemes, 
specifically 
third-order 
Weighted 
Essentially 
Non-Oscillatory (WENO) reconstruction and Adaptive Mesh Refinement 
(AMR), to resolve complex wave interactions over the Philippine Trench 
bathymetry. The resulting hazard assessment quantifies critical metrics 
including tsunami wave height, arrival time, inundation depth, and structural 
exposure, with the physical accuracy of the model validated against 
observational data from the 2012 and 2023 Mw 7.6 tsunamigenic events. 
 
Limitations 
 
Despite the methodological enhancements, this study operates under 
specific technical and physical constraints. First, the accuracy of coastal 
inundation estimates is limited by the resolution of bathymetric and 
topographic data; the model relies partly on global datasets like GEBCO, 
which may not fully resolve fine-scale nearshore features such as coral reefs or 
small man-made structures that influence local run-up behavior. Second, 
regarding source complexity, the study assumes megathrust-only ruptures 
and does not account for secondary tsunamigenic sources such as splay faults, 
outer-rise earthquakes, or submarine landslides, which could act as hazard 
amplifiers in specific local scenarios. Third, while the use of RQMC 
significantly reduces computational cost compared to traditional methods, the 
ensemble size is still constrained by available high-performance computing 
resources, which may limit the precise estimation of tail-risks for extremely 
rare, worst-case Mw 9.0 events. Finally, the assessment of risk is strictly 
geophysical and structural; while building exposure is quantified, complex 
 


---

# PAGE 12

8 
socioeconomic factors such as population dynamics, real-time evacuation 
behavior, and economic loss modeling are beyond the scope of this 
hydrodynamic simulation study. 
 
1.5 Significance of the Study 
 
The primary academic contribution of this research is the development 
of a high-fidelity, fully coupled simulation framework that resolves the 
"computational bottleneck" inherent in Probabilistic Tsunami Hazard Analysis 
(PTHA). Existing studies often rely on simplified static dislocation models to 
save computational time, but recent research indicates that such decouplings 
can lead to significant errors in predicting wave amplitude and arrival times. 
By integrating Randomized Quasi-Monte Carlo (RQMC) sampling, this study 
provides a methodological breakthrough, demonstrating that physically 
realistic slip heterogeneity can be modeled with superior convergence rates 
compared to traditional Monte Carlo methods. This allows for the feasible 
application of the computationally expensive Fully Coupled Method (Method 
1) within a probabilistic framework, setting a new benchmark for accuracy in 
hydrodynamic modeling. Furthermore, this research fills a specific gap in the 
literature by applying these advanced variance-reduction techniques to the 
Philippine Trench, a subduction zone that remains understudied compared to 
the Japan or Nankai trenches despite its potential for generating destructive 
earthquakes. 
 
For the Philippines, specifically the vulnerable communities along the 
eastern seaboard where institutional tsunami preparedness is often limited, 
this study provides actionable, quantitative data to enhance disaster 
resilience. 
Moving 
beyond 
deterministic 
"worst-case" 
scenarios, 
the 
probabilistic outputs such as inundation probability maps and hazard curves 
offer a scientific basis for refining Local Government Unit (LGU) land-use 
plans, such as restricting the construction of critical infrastructure in 
statistically high-risk coastal zones. The enhanced accuracy of the propagation 
model directly supports the calibration of early warning systems by 
establishing site-specific alert thresholds based on physically validated wave 
arrival times rather than generic estimates. Additionally, the study informs 
evacuation planning by identifying areas where designated evacuation 
centers 
may 
be 
compromised 
under specific magnitude scenarios, 
necessitating a re-evaluation of safety protocols for coastal populations. 
Ultimately, this framework contributes technical evidence to support national 
policy updates, including the Philippine Development Plan and the National 
Disaster Risk Reduction and Management Framework, by offering a 
replicable blueprint for science-based risk assessment in data-scarce, 
high-hazard regions. 
 


---

# PAGE 13

9 
 
1.6 Definition of Terms 
 
Philippine Trench: A north–south trending subduction zone along the 
eastern margin of the Philippines, where the Philippine Sea Plate subducts 
westward beneath the Philippine Mobile Belt at rates exceeding 10 cm/year. 
 
Fully Coupled Method: A high-fidelity tsunami modeling approach 
that simultaneously solves the equations of elastodynamics (for seismic wave 
propagation in the solid Earth) and hydrodynamics (for tsunami wave 
propagation in the ocean), capturing the complete multiphysics wavefield. 
 
Stochastic Slip Distribution: A spatially variable model of earthquake 
slip treated as a random field, typically log-normal with spatial correlation, 
used to represent uncertainty in rupture heterogeneity. 
 
Randomized Quasi-Monte Carlo (RQMC): An advanced sampling 
technique that combines low-discrepancy sequences with randomization to 
achieve faster convergence and variance estimation in high-dimensional 
uncertainty quantification. 
 
Propagation Accuracy: The fidelity with which a tsunami model 
simulates wave evolution from deep ocean to coast, particularly in resolving 
inundation extent, run-up height, and arrival time enhanced in this study via 
WENO, AMR, and well-balanced numerics. 
 
Probabilistic Tsunami Hazard Analysis (PTHA): A framework that 
quantifies the likelihood of tsunami impacts (e.g., wave height exceeding a 
threshold) over time, accounting for uncertainties in earthquake source, path, 
and site effects. 
 
Mw (Moment Magnitude): A logarithmic scale measuring the total 
energy released by an earthquake, based on seismic moment M0 ​, defined as 
Mw ​= 3/2 ​ log10 ​(M0​) −10.7. 
 
 


---

# PAGE 14

10 
1.7 Proposed Probabilistic Tsunami Hazard Assessment Framework 
 
Figure 1.1. Proposed Probabilistic Tsunami Hazard Assessment Framework  
 
This study proposes a Physical Fidelity-Based Probabilistic Tsunami 
Hazard Assessment framework designed to rigorously evaluate tsunami 
risks. The methodology commences with the definition of source parameters, 
utilizing Randomized Quasi-Monte Carlo (RQMC) sampling to construct a 
robust probability model that efficiently covers the uncertainty space. These 
probabilistic inputs drive the Generation and Propagation phase, which 
applies "Enhanced Method 1" to simulate wave dynamics. Crucially, a dual 
validation mechanism is integrated at this stage to verify the physical 
accuracy of the simulations before the results are finalized. Ultimately, the 
validated outputs are synthesized through an aggregation process to quantify 
the comprehensive hazard levels.  
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 


---

# PAGE 15

11 
Chapter II 
REVIEW OF RELATED LITERATURE 
 
 
 
 
2.1 Tsunami Hazards and Hydrodynamic Complexity in the Philippine 
Trench 
 
 
The Philippine Trench represents one of the most seismically active yet 
complex subduction zones in the western Pacific. While historical hazard 
assessments often utilized static worst-case scenarios, recent geological events 
have necessitated a shift toward more nuanced, physically validated models. 
 
 
The geodynamic complexity of the Philippine Trench significantly 
influences its seismic and potential tsunami behavior. Lallemand et al. (1998) 
demonstrate that the trench is divided into distinct segments, with the 
deepest portion (~10,100 m) located near 9–10°N, coinciding with the 
maximum depth of the subducting Philippine Sea Plate (~200 km) and the 
oldest arc volcanism (~2.5 Ma). South of 6°N, the trench shallows rapidly, 
interplate seismicity diminishes, and subduction appears to slow down due to 
interaction with the Molucca Sea slab and the Sangihe subduction system. The 
authors interpret this as evidence of decreased plate coupling in the southern 
segment, which may reduce the likelihood of large thrust earthquakes—and 
consequently, major tsunamis—in that region. In contrast, the central segment 
(7–10°N) exhibits stronger mechanical coupling, deeper subduction, and more 
pronounced bending-related faulting, suggesting higher seismic and tsunami 
potential. These tectonic insights are crucial for regional tsunami hazard 
assessments, as they highlight spatial heterogeneity in subduction behavior 
along the Philippine Trench. 
 


---

# PAGE 16

12 
​
 
 
Figure 2.1. Tectonic map of the Philippine and Sangihe Trenches 
 
 
High-resolution P-wave tomography across the southern Philippine 
Trench documents pronounced structural complexity with three subducting 
slabs (Celebes, Molucca, and Philippine Sea) and an overturned geometry of 
the Philippine Sea slab reaching depths of ~450–600 km. Fan & Zhao (2018) 
invert a large dataset of local and teleseismic arrivals to obtain images with 
 


---

# PAGE 17

13 
lateral resolution ≈0.6° and show that the southern segment likely initiated 
subduction at ~20–25 Ma, with a later subduction reorganization about ~5 Ma 
associated with the closure of the Molucca Sea. The resulting slab 
segmentation, depth variations, and lateral heterogeneities inferred from their 
tomographic cross-sections suggest spatially variable seismic coupling and 
fault behavior along the trench conditions that can produce segmented or 
multifocal tsunami sources and complicate near-field hydrodynamic 
responses. Therefore, tomographic constraints such as those presented by Fan 
& Zhao provide essential geophysical context for tsunami hazard models that 
incorporate segmented rupture scenarios, variable slip distributions, and 
site-specific bathymetry.​
 
​
​
The spatiotemporal evolution of subduction along the Philippine 
Trench has direct implications for tsunami potential. Ozawa et al. (2004) used 
K–Ar dating of volcanic rocks along the Philippine Volcanic Arc to 
demonstrate that subduction initiated ~8 Ma in the north (Bicol Peninsula) 
and propagated southward, reaching Leyte by ~3.5 Ma and likely extending 
further south more recently. This north-to-south progression aligns with 
geophysical evidence showing a systematic decrease in trench depth, slab 
length, and interplate coupling toward the south. Consequently, the central 
and northern segments (e.g., 7°–12°N) are tectonically mature and capable of 
generating large megathrust earthquakes, whereas the southern segment 
(south of 6°N) remains immature and weakly coupled, reducing its tsunami 
hazard potential. This segmentation underscores the need for regionally 
differentiated tsunami risk models along the Philippine Trench. 
 
 
2.1.1 Recent Tsunamigenic Events (2012 vs. 2023) 
 
 
The urgency of updating hazard models for this region is 
underscored by recent empirical data. Heidarzadeh et al. (2025) 
performed a comparative analysis of two major earthquakes (Mw 7.6) 
that occurred along the trench in August 2012 and December 2023. 
Their study revealed that despite having identical magnitudes, the 
resulting tsunamis exhibited drastically different hydrodynamic 
behaviors due to the rupture depth and bathymetry. The 2023 event 
produced waves with significantly longer periods (6.7–28.2 minutes) 
compared to the 2012 event (8.0–18.3 minutes), a phenomenon 
attributed to the deeper water column at the 2023 source. This finding 
supports the argument that simple magnitude-based projections are 
insufficient; hazard models must account for the spatial heterogeneity 
of the rupture to capture these complex wave characteristics. 
 
 
 


---

# PAGE 18

14 
2.1.2 Vulnerability of Coastal Infrastructure 
 
 
The gap between scientific understanding and infrastructure 
readiness remains wide. Simulation data from the 2025 study 
highlighted that an Mw 8.5 event well within the trench's capacity 
would likely inundate critical facilities, including the designated 
evacuation center in the coastal town of Dapa. This exposure of "safe 
zones" indicates a failure in current deterministic planning maps, 
validating the need for the probabilistic inundation mapping proposed 
in this study. 
 
 
Geomorphic evidence from uplifted marine terraces along 
eastern Mindanao provides critical insight into the tsunami potential of 
the southern Philippine Trench. Ramos et al. (2012) identified up to 
four Holocene terrace levels (1–12 m above sea level) along the Davao 
Oriental coast, with radiocarbon ages between 8,080 and 4,140 cal yr 
BP. The staircase morphology and 1–3 m terrace risers are interpreted 
as the result of repeated coseismic uplift events, likely caused by large 
to great (M ≥ 8) megathrust earthquakes along the Philippine Trench. 
Notably, historical earthquakes in the same region (e.g., the 1992 Mw 
7.1–7.2 events) generated a 5-meter tsunami but no detectable coastal 
uplift, suggesting that the terraces record significantly larger 
prehistoric ruptures. This implies that the southern Philippine Trench, 
though tectonically younger and less seismically active than its central 
counterparts remains capable of generating tsunamigenic earthquakes, 
challenging assumptions of uniformly low hazard in this segment. 
 
 


---

# PAGE 19

15 
 
 
 
Figure 2. 2. Map of uplifted marine terraces in Davao Oriental 
 
2.2 Probabilistic Tsunami Hazard Assessment (PTHA) 
 
Probabilistic Tsunami Hazard Assessment (PTHA) applies formal 
probabilistic methods (adapted from seismic hazard analysis) to quantify the 
likelihood of tsunami waves of various sizes. In PTHA, all potential tsunami 
sources and their occurrence rates are defined, tsunamis are simulated for 
 


---

# PAGE 20

16 
each scenario, and the results are combined into exceedance‐probability 
curves or maps. The intensity measure (IM) is a physical tsunami metric (e.g. 
coastal wave amplitude or run-up height). PTHA differs from deterministic 
“worst-case” approaches by assigning probabilities to scenarios, yielding 
hazard curves such as the annual probability of exceeding a given run-up. 
Formally, PTHA estimates the mean annual rate λ(IM ≥ im) that a tsunami 
intensity IM at a location exceeds threshold im . This is done by integrating 
tsunami effects (from numerical models) with seismic event rates (from 
magnitude–frequency models). In practice, PTHA proceeds in three core 
stages: 
 
 
Seismic Source Modeling: Define all tsunamigenic sources (e.g. 
subduction megathrusts, crustal faults, submarine slides, volcanoes) and 
assign earthquake occurrence rates. Sources may be zoned by tectonic region 
or fault, and each is given a magnitude–frequency distribution (commonly a 
Gutenberg–Richter (G–R) law or a “characteristic” model). The seismicity 
model (Poissonian or renewal) yields the annual rate of each earthquake 
magnitude. Logic trees or multiple branches are often used to capture 
epistemic uncertainty in parameters (e.g. different b‐values or maximum 
magnitudes). For example, Sorensen et al. (2012) applied truncated G–R 
relations for Mediterranean source zones by fitting catalog data, while Baba et 
al. (2022) used a G–R model for Nankai Trough earthquake rates.  
 
 
Tsunami Simulation: For each seismic scenario (defined by location, 
fault geometry, magnitude, slip distribution, etc.), compute the resulting 
tsunami using numerical wave-propagation models (shallow-water or 
Boussinesq models with nested grids). High-resolution inundation models are 
often used nearshore. Many studies simulate seafloor displacement and wave 
evolution for each event. The output is an intensity measure (wave height, 
flow depth) at coastal sites. For example, Gonzalez et al. (2009) developed the 
first PTHA inundation maps by coupling PSHA with inundation modeling 
for Seaside, Oregon. They computed 100- and 500-year tsunami amplitudes 
(1% and 0.2% annual exceedance) from both far-field (Alaska) and near-field 
(Cascadia) earthquakes. 
 
 
Statistical Hazard Estimation: Combine the tsunami outputs with 
event rates to build hazard curves. Each earthquake scenario has a weight (its 
annual rate) from the seismic model, and produces a maximum IM at the site. 
By summing over all scenarios (e.g. via Monte Carlo or numerical 
integration), one obtains the probability that IM exceeds given levels in a time 
window. The result can be expressed as a curve of excess probability vs. IM, 
or as maps of IM at specified return periods. UNESCO guidelines describe 
PTHA as first “establishing parameters for all possible tsunami sources and 
occurrence rates” then simulating and aggregating results into exceedance 
curves. In practice, PTHA often uses Monte Carlo: e.g. Geist and Parsons 
 


---

# PAGE 21

17 
(2006) 
randomized 
hypocenter 
and 
slip 
to 
build 
a 
tsunami 
amplitude–probability curve.  
 
 
2.2.1 Treatment of Uncertainty 
 
PTHA explicitly quantifies uncertainty. Uncertainties are usually 
divided into epistemic (due to lack of knowledge or competing models) 
and aleatory (intrinsic randomness). Epistemic uncertainties (e.g. 
which seismic source model is correct, how to represent fault 
geometry) are handled via logic trees or ensembles of alternative 
assumptions, each branch weighted by credibility. Aleatory variability 
(e.g. random slip heterogeneity, tide level, simulation error) is treated 
statistically by probability density functions. Common examples 
include:  
 
 
Logic Trees: Many PTHAs use logic trees to capture uncertainty 
in source characterization and magnitude-frequency. For instance, 
alternative Gutenberg–Richter b‐values or corner magnitudes might be 
assigned weights, and the final hazard curve is averaged over 
branches. Geist and Parsons (2005) compared “characteristic” vs. G–R 
magnitude distributions for Cascadia by Monte Carlo, illustrating how 
different recurrence models affect hazard. 
 
 
Slip Variability: The spatial distribution of fault slip in an 
earthquake can strongly affect tsunami size. Modern PTHAs often 
generate multiple slip realizations per rupture to sample aleatory slip 
variability. For example, Mazet‐Roux et al. (2019) modeled each 
rupture area with five stochastic slip distributions based on 
overlapping circular “asperities,” enforcing a k^(-2) decay of slip 
spectrum. This approach captures the natural heterogeneity of slip and 
its effect on tsunami generation.   
 
Tidal and Sea‐Level Uncertainty: Background water level (tide 
or sea level) adds noise to run-up predictions. PTHAs incorporate this 
by sampling over tidal phases or adding a distribution of sealevel 
offsets. Baba et al. (2022) explicitly included astronomical tide 
variations in their hazard curves. The U.S. Crescent City study used a 
“pattern method” to model tidal uncertainty when mapping 
inundation. More generally, stochastic sampling of tide phase and 
projected sea‐level rise is becoming standard in recent assessments. 
 
 
 


---

# PAGE 22

18 
Model and Measurement Uncertainty: Differences between 
modeled and observed tsunami heights introduce additional aleatory 
spread. Davies et al. (2017) noted that model observation mismatches 
were larger than previously assumed, implying that failure to include 
this uncertainty underestimates hazard. They found that explicitly 
accounting for such deviations increases predicted run-ups for a given 
exceedance probability.   
 
 
In all cases, uncertainties are propagated through to the final 
hazard (e.g. via Monte Carlo or convolution). The output is typically 
given with confidence bounds (from logic-tree weighting) on the 
hazard curves or maps. For example, UNESCO guidelines remark that 
PTHA “explicitly addresses different types and sources of uncertainty” 
by developing and weighting alternative models. 
 
 
2.2.2 Integration with Seismic Hazard Models 
 
PTHA is closely linked with seismic hazard assessment. The 
earthquake catalogs and frequency models used in PTHA usually come 
from PSHA-type analysis. Magnitude–frequency parameters (b‐value, 
corner magnitude) and source zonations from seismic studies are 
adopted or adjusted for tsunami application. For example, many 
PTHAs use Gutenberg–Richter fits to regional seismicity as the basis 
for earthquake occurrence rates. In the Mediterranean PTHA, Sørensen 
et al. fitted a truncated G–R law per source zone using historical 
earthquake catalogs. Similarly, Baba et al. based Nankai Trough 
earthquake probabilities on a G–R model. Alternative seismic hazard 
branches (e.g. using different b‐values or a characteristic event 
assumption) can be included in the logic tree. Thus, PTHA builds on 
PSHA methods: earthquake event rates (λ(M)) feed directly into 
tsunami event rates, and Poissonian or time-dependent occurrence 
models from PSHA carry over. The strong coupling means that any 
update in seismic hazard (new sources, revised recurrence) will 
directly affect tsunami hazard estimates.   
 
2.2.3 Case Studies by Region 
 
Probabilistic tsunami hazard studies have been carried out 
worldwide, ranging from local to basin scales. Representative examples 
include: 
 
 
 


---

# PAGE 23

19 
Japan – Nankai Trough: Baba et al. (2022) applied PTHA to 
eastern Shikoku, Japan. They generated 3480 rupture scenarios of past 
Nankai 
earthquakes 
and 
computed 
corresponding 
tsunamis. 
Earthquake probabilities were based on a Gutenberg–Richter model, 
and tidal variation was included. The resulting hazard curves 
(exceedance vs. wave height) agreed well with historical tsunami 
frequencies. Sensitivity tests revealed that slip heterogeneity and 
coastal defenses significantly affect hazard levels. The study also 
produced high-resolution inundation maps with probabilistic contours.  
 
 
Pacific Northwest (Cascadia, USA): Gonzàlez et al. (2009) 
performed a PTHA for Seaside, Oregon (Cascade Range region). They 
combined 
tsunami 
inundation 
modeling 
with 
PSHA-derived 
earthquake rates to produce maps of 100- and 500-year tsunami 
amplitudes. They found that 100- year tsunamis are usually from 
far-field sources (Alaska–Aleutians) with <4 m amplitudes, whereas 
500-year events are dominated by local Cascadia megathrust ruptures 
(>10 m). The main uncertainties were in earthquake recurrence 
intervals, sea level, and bathymetry. Similar PTHAs have been 
conducted for Crescent City, California, with updated source sets and 
tidal uncertainty treatment, confirming that 100- and 500-year 
inundation extents are very sensitive to seismic-source assumptions 
 
 
Mediterranean Sea: Sørensen et al. (2012) presented the first 
basin-wide PTHA for the Mediterranean. Using Monte Carlo 
simulation with zone-specific G–R seismicity models, they estimated 
the annual probability of exceeding a given tsunami amplitude at any 
coast. The highest hazard was along the Hellenic Arc (Greece), but 
almost all coastlines showed non-negligible risk. Notably, they 
estimated that the probability of anywhere in the Mediterranean 
experiencing a tsunami ≥1 m in the next 30 years is nearly 100%. Their 
method identifies the dominant source (region) for hazard at each 
location and can inform warning systems.   
 
 
Global (Earthquake-Generated Tsunamis):Davies et al. (2017) 
performed a global PTHA covering all major subduction zones 
worldwide. They considered only earthquake sources (≈80% of 
historical tsunamis) and generated global maps of tsunami run-up for 
specified exceedance rates. This work highlighted that large 
uncertainties (especially in the frequency of great earthquakes) produce 
large spread in predicted run-ups. They also found that accounting for 
modeling errors increases predicted hazard: models tended to 
under-predict historical run-ups, implying a higher real hazard than 
naive projections. 
 
 


---

# PAGE 24

20 
 
Other Regions: PTHA has also been done for many other areas. 
For example, Tang et al. (2017) assessed PTHA for Southeast 
China/Taiwan, Hokkaido (Japan), Alaska, and several national studies 
in Chile and Indonesia. In each case, local seismicity and source faults 
are combined with tsunami simulation to produce site-specific hazard 
curves. While this study cannot list all studies, the general pattern is 
similar: evaluate a suite of earthquakes (with varied magnitudes and 
slips), run tsunami models, and compute exceedance probabilities.   
 
 
Table 2.1. Tsunami Modeled sources Comparison 
 
 
 
 
 
 
 
 


---

# PAGE 25

21 
2.2.4 Advances in PTHA Modeling 
 
Recent developments are expanding and accelerating PTHA: 
 
 
Machine-Learning 
Surrogates: 
High-fidelity 
tsunami 
simulations are computationally intensive, especially when thousands 
of scenarios are needed. New studies train ML models to predict 
tsunami impacts. For example, Ramalingam et al. (2025) ran hundreds 
of tsunami simulations for Tōhoku subduction earthquakes and 
trained a variational encoder–decoder network to predict coastal 
waveforms and maximum inundation depths from a small number of 
deep-ocean 
input values. The ML surrogate reproduced the 
full-physics results with high accuracy, enabling rapid hazard mapping 
along long coastlines with limited simulation runs. Such surrogate 
models (including neural nets and Gaussian-process emulators) show 
great promise for accelerating PTHA. 
 
  
Stochastic Rupture Models: Better representation of slip 
variability improves realism. Recent PTHA studies use stochastic 
rupture generators (e.g. self-similar spectral models) to sample many 
possible slip patterns. For instance, Mazet‐Roux et al. (2019) generated 
five random slip distributions per rupture by superimposing circular 
sub‐faults (“asperities”) to enforce a k^(-2) spectral decay. This 
captures aleatory slip heterogeneity in each scenario. Similar 
approaches (random kinematic slips) allow PTHA to automatically 
include slip randomness instead of assuming a uniform slip on the 
fault. 
 
 
Multi‐Hazard (“Hybrid”) Frameworks: PTHA is evolving into 
a comprehensive multi-source framework. Grezio et al. (2021) discuss 
“Total PTHA (TotPTHA)” methods that merge multiple tsunamigenic 
sources 
earthquakes, 
submarine 
slides, 
volcanic 
blasts, 
even 
meteorological or anthropogenic triggers into a single probabilistic 
assessment . The TotPTHA formalism accounts for source interactions 
(e.g. an earthquake triggering a landslide) and sums over all source 
types to compute a total tsunami hazard. In essence, each source class 
contributes its own PTHA, and these are combined (with possible 
conditional probabilities) to get the overall hazard. While still 
emerging, this approach reflects the growing recognition that 
non-seismic tsunamis (landslides, volcanoes) can be important locally, 
and that all must be considered for a “complete” hazard analysis.  
 
 
 


---

# PAGE 26

22 
Statistical Emulation and Fast Schemes: Beyond deep learning, 
simpler surrogate approaches are also used. Emulators built from 
regression or machine-learning reduce the number of required full 
simulations. The “Pattern Method” and other statistical schemes 
address specifically the distribution of tidal stages or other factors. 
Such methods are being incorporated into PTHA to speed up 
computations while quantifying uncertainty.  
 
 
Time-Dependent Models: Traditional PTHA assumes Poisson 
(time-invariant) occurrence, but some research is exploring renewal 
models. For subduction zones with quasi-periodic great quakes (e.g. 
Cascadia), 
a 
time-dependent 
hazard 
analysis 
can 
modulate 
probabilities based on elapsed time since the last event. These 
approaches borrow from seismic renewal theory and are an active 
research area. 
 
 
Sea-Level and Land-Motion Effects: Recent PTHAs (especially 
for long-term projections) now include trends in sea-level rise and 
vertical land motion. For instance, Grezio et al. (2021) incorporate 
anticipated sea-level increase into the hazard curves for Mediterranean 
tsunamis, showing higher future inundation probabilities. Such 
extensions ensure PTHA remains relevant under climate change. 
 
 
Inclusion of Non-Seismic Sources: As noted by Geist and 
Lynett (2014), including landslides and volcanic tsunamis greatly 
increases complexity . Nevertheless, some PTHAs now attempt it 
(using slide-volume vs. run-up relationships or empirical distributions 
of collapse events). Volcanic flank collapse tsunamis and seiche-type 
events are increasingly being studied within a probabilistic framework, 
though data scarcity means these remain high-uncertainty components.  
 
 
2.3 Gutenberg–Richter Law: Earthquake Frequency Distribution 
 
The Gutenberg–Richter (G–R) law is a foundational empirical relation 
in seismology. It states that the number of earthquakes N with magnitude ≥  
follows a log-linear distribution:  
 
 


---

# PAGE 27

23 
Where a and b are constants (the “b ‑value” is typically ≈ 1 for tectonic 
earthquakes). Equivalently, the frequency (per unit time) of earthquakes scales 
as 10-bM . Thus, for b≈1 , a one-unit decrease in magnitude produces about ten 
times more events. The G–R law holds over a wide range of magnitudes and 
regions, reflecting the scale-invariance of earthquake processes. (More 
complex variants, like the tapered G– R law, introduce an exponential cutoff 
to limit the occurrence of unrealistically large magnitudes.) 
  
 
Gutenberg and Richter originally formulated this law in 1944 for 
California seismicity. It has since been observed globally, underpinning most 
Probabilistic Seismic Hazard Assessment (PSHA) methods. In these 
assessments, the G–R relation is used to estimate the rate of earthquakes of 
different sizes in a region. In practice, seismologists fit observed catalogs to 
determine a and b , and apply log-likelihood or Bayesian methods to account 
for limited data. For example, a typical modern PTHA will assume a Poisson 
(memoryless) 
occurrence 
process 
with 
a 
G–R 
magnitude-frequency 
distribution for each seismic source zone 
 
 
Key points: The G–R law (log 10 N = a – b M) describes how often large 
vs. small earthquakes occur, with a characteristic b≈1 for global tectonic 
settings. This scale-invariant law is treated as a statistical baseline for 
earthquake rates in hazard models.   
 
 
2.3.1 Application to Tsunami Hazard Analysis  
 
 
In 
Probabilistic 
Tsunami 
Hazard 
Assessment 
(PTHA), 
earthquake occurrence rates are a primary input. Tsunami modeling 
typically adapts the PSHA framework: one defines seismic source 
zones and assigns them magnitude-frequency distributions (often via 
G–R), then simulates tsunami generation and propagation for sampled 
earthquakes . Many PTHA studies explicitly use G–R statistics for 
tsunami sources. For example, Baba et al. (2022) modeled 3480 
candidate ruptures in the Nankai subduction zone and “assumed that 
the probability of earthquake occurrence was based on the 
Gutenberg–Richter law” . They then computed tsunami heights for 
each scenario and found that the resulting hazard curves matched 
historical tsunami frequencies.   
 
 
In practice, hazard analysts often use logic trees to capture 
uncertainty in the magnitude-frequency model. Two common 
hypotheses are: (1) a pure G–R (continuous) distribution of 
 


---

# PAGE 28

24 
magnitudes, or (2) a “characteristic earthquake” model (a few large 
events of fixed size). For example, Geist and Parsons (2005) contrasted 
a characteristic Mw =9 scenario with a G–R ensemble in Cascadia 
PTHA. The choice affects tsunami forecasts: a higher ‑value (steeper 
falloff) reduces the weight of giant earthquakes, whereas a 
characteristic model can inject more weight at the top magnitude.   
 
 
Indeed, Zimmerman et al. (2025) note that “many previous 
PTHA studies… use different forms of on-fault M– F relations (e.g. 
characteristic, Gutenberg–Richter) in logic trees to represent epistemic 
uncertainty”. In earthquake‐tsunami modeling, the G–R law is often 
combined with physical slip models: e.g., Tong et al. (2023) generated 
stochastic fault slips with magnitudes drawn to follow a G–R 
distribution (Mw 7–9) and then computed tsunami heights from each 
slip. Their method identifies which slips are most likely to produce 
large tsunamis, effectively “deriving tsunami height probabilities” 
from an earthquake-driven G–R ensemble. Globally, PTHA efforts 
confirm the primacy of earthquakes as tsunami sources. Davies et al. 
(2017) conducted a global PTHA (80% of recorded damaging tsunamis 
are 
earthquakes‐generated) 
and 
emphasized 
that 
“epistemic 
uncertainties in the exceedance rates of large earthquakes often lead to 
large uncertainties in tsunami run-up” . Likewise, regional PTHAs (e.g. 
Gonzalez et al. 2009 for Seaside, Oregon) integrate tsunami wave 
modeling with PSHA-derived earthquake rates. In all cases, the G– R 
model provides the baseline earthquake occurrence rate, which is then 
used to weight tsunami scenario frequencies. 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 


---

# PAGE 29

25 
Table 
2.2. 
Examples 
of 
PTHA 
studies 
and 
their 
use 
of 
Gutenberg–Richter distributions. 
 
 
(Sources: see text; ERF = earthquake rupture forecast.)​​
 
 
2.3.2 Tsunami Size–Frequency Relations  
 
 
Beyond earthquakes, some studies seek tsunami-specific 
size-frequency laws. Fine et al. (2020) analyzed 32 years of deep-ocean 
pressure data and identified 41 far-field tsunami events. They found 
that tsunami amplitudes themselves follow a power-law (Pareto) 
distribution. In plain language: “the data were used to derive an 
empirical power-law relationship between the ‘size’ (amplitude) of 
tsunami waves and their frequency of occurrence,” establishing a 
tsunami size–frequency law analogous to Gutenberg–Richter. In other 
words, once tsunamis are defined (e.g. by their first-wave height), their 
occurrence also roughly follows a N∝(amplitude)−p law. This is an 
emerging area: if robust, such a law could allow probabilistic tsunami 
hazards to be estimated directly from past tsunami records, similarly to 
how G–R is used for earthquakes. However, unlike earthquakes, 
tsunamis do not have a standardized magnitude measure (so far). Most 
PTHAs still rely on the seismic G–R framework. Tsunami-specific 
power laws (like Fine et al.’s) remain exploratory due to limited event 
counts. Nonetheless, these studies highlight that tsunamis themselves 
 


---

# PAGE 30

26 
exhibit fat-tailed statistics, justifying power-law models in hazard 
analysis. 
 
 
2.3.3 Shallow vs. Deep Earthquakes in Tsunami Generation  
 
Depth matters greatly: tsunamis are overwhelmingly generated 
by shallow submarine earthquakes. For example, NOAA notes that 
“for an earthquake to pose a tsunami hazard it needs to vertically 
move the seafloor; therefore it needs to be large (typically ≧8.0), under 
or near the ocean, and shallow (<100 km)” Likewise, the majority of 
great 
earthquakes 
(Mw ≥ 8) 
occur 
on 
subduction megathrusts 
(convergent plate boundaries) . In practice, PTHAs concentrate on 
shallow interplate thrust events and often ignore deeper intra-slab 
quakes.  
 
 
Deep 
earthquakes 
(hypocenters 
≳ 300–400 km) 
normally 
produce negligible tsunamis because their seafloor displacement fields 
are weakly coupled to the ocean. Okal (2017) rigorously analyzed the 
2013 Sea of Okhotsk earthquake (depth 603 km, Mw 8.3) which 
surprisingly generated a tiny tsunami. He showed theoretically that 
tsunami energy from a deep point-source scales as ~ M2
0/z2
s (where M0 
is seismic moment and depth), versus the classical ~ M4/3 scaling for 
shallow sources. In other words, at great depth the tsunami amplitude 
is severely attenuated (~1/z2), implying only extremely large deep 
quakes (larger than any yet observed) could approach hazard levels. 
This confirms the paradigm that almost all tsunamis come from 
shallow (typically <100 km depth) event. 
 
 
Other 
tsunami‐specific 
nuances: 
Some 
rare 
“tsunami 
earthquakes” defy usual scaling. These are unusually slow ruptures 
(often at shallow, compliant parts of megathrusts) that yield larger 
tsunami waves than their Mw would suggest. Behrens et al. (2021) note 
that “tsunami earthquakes produce excessively large tsunami 
intensities compared to their moment magnitude” and that their 
frequency is “unconstrained” . Standard G–R scaling (based purely on 
magnitude) cannot capture this effect. In practice, PTHAs may include 
a special branch of the logic tree to allow greater slip or efficiency for 
such events, but their statistics remain highly uncertain 
 
 


---

# PAGE 31

27 
2.3.4 Critiques and Limitations of Gutenberg–Richter in Tsunami 
Context  
 
While the Gutenberg–Richter law is a useful starting point, it 
has important caveats in tsunami applications: 
 
●​ Non-Poissonian 
Behavior: 
The classical G–R framework 
assumes earthquake occurrences follow a Poisson process 
(memoryless in time). However, seismicity often clusters in 
space–time (aftershocks, seismic cycles). Kagan and Jackson 
(1991) showed clear evidence of clustering, and most PTHAs 
still use a time-independent Poisson model . Ignoring this may 
under- or overestimate tsunami risk on human timescales (e.g. 
following a megathrust rupture). Only a few studies have 
implemented time‐dependent earthquake models in PTHA  
. 
●​ Upper‐Magnitude Bound: A pure power-law implies an infinite 
tail of large events, which is unphysical. In reality, tectonic 
regions have finite strain energy, and the largest earthquakes 
might fall off faster than G–R predicts. To address this, some 
analyses use a tapered G–R or impose a hard Mmax . Otherwise, a 
simple G–R can overstate the likelihood of truly giant 
tsunamigenic quakes. The characteristic‐earthquake alternative 
(one or few max-size events) stems from this concern. 
 
●​ Spatial Domain Effects: The G–R fit depends on the chosen 
spatial volume. Krushelnitskii et al. (2024) argue that the law’s 
validity can break down if the catalog area is too small or 
includes heterogeneous fault populations. For example, one 
subduction segment might show an apparent deficit of 
moderate quakes (making large events seem too frequent by 
G–R) or vice versa. This means that applying a single global 
‑value to all tsunami sources could be misleading.   
 
●​ Seismic vs. Tsunami Productivity: Gutenberg–Richter addresses 
earthquake counts, but not all large earthquakes produce 
tsunamis. Thrust events on continental plates, strike‐slip 
earthquakes, or very deep shocks are poor tsunami generators. 
Thus, even if G–R correctly predicts earthquake rates, additional 
 


---

# PAGE 32

28 
filters (magnitude ≥ M0, fault type, location) must be applied. 
The law itself doesn’t account for these (they are implemented 
in PTHA source definitions).   
 
●​ Tsunami Earthquake Complexity: As noted, certain events 
(tsunami earthquakes, splay faults, rupture jumps) can produce 
outsized tsunamis. These complex ruptures violate the 
assumptions of simple, self-similar scaling. Behrens et al. 
emphasize that such events “are not well represented in PTHA” 
. In effect, the G–R-based models may under-predict tsunami 
hazard if such anomalous slips occur with non-negligible 
probability.  
 
●​ Parameter Uncertainty: Even if G–R holds on average, the 
parameters a and b have statistical uncertainties from finite 
data. In PTHA, analysts propagate these (e.g. by logic-tree 
branch weights). For tsunami risk, this means a small change in 
b or Mmax or can alter the predicted exceedance rates 
significantly. 
 
●​ Limited Tsunami Records: Unlike seismic catalogs, tsunami 
records are sparse (few tsunamis per century in many locations). 
Using earthquake-based G–R bypasses this, but it also means 
the validation of tsunami probabilities is challenging. As Baba et 
al. (2022) found, PTHA often relies on few historical tsunamis 
for validation .   
 
2.4 The Physics of Tsunami Generation: Coupled vs. Uncoupled Models 
 
A critical methodological choice in modern tsunami science is the 
treatment of the generation phase, the transfer of energy from the solid earth 
to the fluid ocean with the accuracy of the resulting hazard assessment, 
particularly in the near-field, fundamentally depending on the fidelity of this 
physical transfer mechanism (Abrahams et al., 2023). Tsunami generation 
models are classified along a spectrum ranging from uncoupled kinematic 
methods to high-fidelity fully coupled dynamic methods (Abrahams et al., 
2023]; Abrahams et al., 2023). Uncoupled kinematic models, such as Method 4 
(Static Superposition), use a static end-state of seafloor deformation and apply 
it instantaneously to the water surface, neglecting the time-dependent physics 
of the rupture (LeVeque et al., 2016). In contrast, the highest fidelity model, 
 


---

# PAGE 33

29 
the Fully Coupled Dynamic Method (Method 1), simultaneously solves the 
elastic wave equation and hydrodynamic equations, explicitly modeling the 
two-way pressure and stress exchange at the seafloor interface (Abrahams et 
al., 2023).​
​
​
​
​
 
 
 
Figure 2.3. Summary of Methods by  Abrahams et al., 2023 
 
2.4.1 Limitations of Superposition (Method 4) 
 
 
Traditional models often use the Superposition Method (Method 
4), which calculates seafloor deformation statically and applies it 
instantaneously to the water surface. Abrahams et al. (2023) rigorously 
evaluated this approach against fully coupled simulations and 
identified fundamental physical flaws. They demonstrated that 
Method 4 neglects "non-hydrostatic filtering," the natural process 
where the deep ocean column dampens short-wavelength seafloor 
deformations. Consequently, Superposition models can overpredict 
tsunami wave amplitudes by a factor of two, leading to overly 
conservative and expensive engineering requirements. 
 


---

# PAGE 34

30 
​
 
 
            Figure 2.4. Comparison of Tsunami Generation Modeling Approaches 
 
2.4.2 Validation of the Fully Coupled Method (Method 1) 
 
In contrast, Abrahams et al. (2023) found that the Fully Coupled 
Method (Method 1), which simultaneously solves the elastic and 
hydrodynamic equations, captures the correct wave physics, including 
the timing of the initial wave peak. This is further supported by 
Baragamage and Wu (2024), who extended fully coupled modeling to 
three-dimensional environments. Their work on fluid-structure 
interaction for coastal bridges confirmed that coupled models are 
essential not just for wave propagation, but for accurately predicting 
the hydrodynamic loading on infrastructure during the inundation 
phase. 
​
​
​
​
 
​2.4.3 The Evolution to Dynamic Rupture-Tsunami Coupling 
 
The limitations of the Superposition Method (Method 4) are 
severely compounded when applied to modern Probabilistic Tsunami 
Hazard Analysis (PTHA), which traditionally relies on simplified 
kinematic or stochastic slip models that lack the time-dependent 
 


---

# PAGE 35

31 
physics of the earthquake rupture itself (Sara Aniko Wirp et al., 2020). 
The latest generation of fully coupled models addresses this by 
integrating 3D dynamic rupture simulations with tsunami generation 
and propagation (Kutschera et al., 2024; Kutschera et al., 2023). 
Dynamic rupture models are inherently physics-based, relying on 
initial stress conditions and fault friction to spontaneously generate a 
self-consistent, time-dependent slip distribution and the resultant 
vertical seafloor displacement (Ulrich et al., 2024; Gabriel et al., 2023). 
This dynamic approach is critical because it successfully captures 
near-source effects that kinematic models omit, such as near-surface 
rake rotation and large, dynamically evolving shallow fault slip 
(Kutschera et al., 2024; Kutschera et al., 2023). For instance, research on 
the Húsavík–Flatey Fault Zone demonstrated that dynamic rake 
rotation interacting with local bathymetry significantly explains 
higher-than-expected vertical seafloor displacements, thus generating a 
larger tsunami potential, which is a key mechanism in strike-slip 
faulting tsunami generation (Kutschera et al., 2024). Ultimately, 
comparisons of dynamically sourced tsunami waveforms against 
one-way linked (uncoupled) simulations confirm that while arrival 
times are generally similar, complex source effects, such as tsunami 
dispersion and the superposition of seismic and acoustic waves, arise 
only in the fully coupled 3D simulations (Kutschera et al., 2023). This 
provides further physical evidence supporting the argument by 
Abrahams et al. (2023) regarding the superiority of the fully coupled 
approach for accurately modeling complex, near-field hazards 
(Abrahams et al., 2023). 
 
2.4.4 The Evolution to Dynamic Rupture-Tsunami Coupling 
 
The limitations of the Superposition Method (Method 4) are 
severely compounded when applied to modern Probabilistic Tsunami 
Hazard Analysis (PTHA), which traditionally relies on simplified 
kinematic or stochastic slip models that lack the time-dependent 
physics of the earthquake rupture itself (Sara Aniko Wirp et al., 2020). 
The latest generation of fully coupled models addresses this by 
integrating 3D dynamic rupture simulations with tsunami generation 
and propagation (Kutschera et al., 2024; Kutschera et al., 2023). 
Dynamic rupture models are inherently physics-based, relying on 
initial stress conditions and fault friction to spontaneously generate a 
self-consistent, time-dependent slip distribution and the resultant 
vertical seafloor displacement (Ulrich et al., 2024; Gabriel et al., 2023). 
This dynamic approach is critical because it successfully captures 
near-source effects that kinematic models omit, such as near-surface 
rake rotation and large, dynamically evolving shallow fault slip 
 


---

# PAGE 36

32 
(Kutschera et al., 2024; Kutschera et al., 2023). For instance, research on 
the Húsavík–Flatey Fault Zone demonstrated that dynamic rake 
rotation interacting with local bathymetry significantly explains 
higher-than-expected vertical seafloor displacements, thus generating a 
larger tsunami potential, which is a key mechanism in strike-slip 
faulting tsunami generation (Kutschera et al., 2024). Ultimately, 
comparisons of dynamically sourced tsunami waveforms against 
one-way linked (uncoupled) simulations confirm that while arrival 
times are generally similar, complex source effects, such as tsunami 
dispersion and the superposition of seismic and acoustic waves, arise 
only in the fully coupled 3D simulations (Kutschera et al., 2023). This 
provides further physical evidence supporting the argument by 
Abrahams et al. (2023) regarding the superiority of the fully coupled 
approach for accurately modeling complex, near-field hazards 
(Abrahams et al., 2023). 
 
2.4.5 Application in Probabilistic Tsunami Hazard Analysis (PTHA) 
 
The integration of fully coupled, dynamic rupture methods is 
now enabling a paradigm shift in PTHA by incorporating physically 
consistent source variability. Traditional PTHA uses logic trees based 
on pre-defined rupture segments and simplified kinematic slip 
variability (Baoning Wu et al., 2025). This approach is being replaced 
by physics-based dynamic rupture scenarios which provide a 
self-consistent and mechanically viable range of rupture descriptions 
across complex fault networks (Li et al., 2023; Ibrahim et al., 2024). This 
methodology 
effectively 
addresses 
the 
need 
for 
a 
better 
characterization of hazard in tectonically complex regions, thus 
complementing classical seismic hazard assessment methods (Li et al., 
2023). This shift toward generating an ensemble of mechanically sound 
earthquake scenarios aligns with the use of advanced sampling 
techniques like Randomized Quasi Monte Carlo (RQMC), which can be 
applied to the resulting high-fidelity slip distributions to reduce 
numerical uncertainty in hazard forecasts (Nakayama & Tuffin, 2024). 
 
2.4.6 Capturing Transient Wave Physics: Seismic and Acoustic 
Coupling 
 
Beyond the issues of static overprediction and kinematic source 
simplification, the fully coupled methodology is further justified by its 
unique ability to capture complex, time-dependent wave phenomena 
 


---

# PAGE 37

33 
that the Superposition Method entirely neglects (Kutschera et al., 2024). 
The most significant of these effects is the Acoustic-Gravity Wave 
(AGW) coupling, which involves the simultaneous propagation of 
pressure disturbances through the solid earth, the water column, and 
the atmosphere. In a fully coupled 3D solid-fluid environment, seismic 
waves generated by the earthquake propagate into the water as 
acoustic (pressure) waves, resulting in a transient motion of the sea 
surface (Kutschera et al., 2023). On the shallow continental shelf, these 
acoustic wave amplitudes can be unexpectedly high, sometimes 
exceeding the amplitude of the actual tsunami wave, a phenomenon 
which may serve as a rapid indicator of a surface-breaking dynamic 
rupture and is crucial for early warning systems (Kutschera et al., 
2023). Furthermore, advanced coupled models have begun to 
incorporate poroelastic effects the interaction between the solid rock 
skeleton and the pore fluid pressure which is critical for accurately 
modeling dynamic rupture processes in saturated media near the 
trench and for understanding fluid inertia effects in shallow 
subduction environments (Ibrahim et al., 2024; Ulrich et al., 2024). The 
ability of the fully coupled framework to resolve these multi-physics 
processes, including the solid and fluid inertia, is what ultimately 
validates its superior representation of the real-world physics over 
simplified two-step methods (Abrahams et al., 2023).  
 
2.5 Probabilistic Tsunami Hazard Assessment (PTHA) Paradigms 
 
The transition from single-scenario simulations to Probabilistic 
Tsunami Hazard Assessment (PTHA) is driven by the necessity to rigorously 
quantify the large epistemic and aleatory uncertainties inherent in 
tsunamigenic earthquake sources, establishing PTHA as the gold standard for 
regional risk assessment (Baoning Wu et al., 2025). This approach is critical 
because simplified or "average" slip models fail to capture the extreme, 
tail-end risks posed by the heterogeneous rupture patches common in active 
regions like the Philippine Trench (Mulia et al., 2020). Beyond seismic 
uncertainty, the framework's flexibility allows it to accommodate other 
independent sources of hazard, such as the integration of sea-level rise due to 
climate change, which ensures that hazard assessments remain relevant across 
decades (Alhamid et al., 2022). 
 
This comprehensive approach necessitates a fundamental shift in both 
source generation and hydrodynamic modeling to maintain physical fidelity 
across the entire process. While traditional PTHA often relied on simpler 
kinematic source models, the field is rapidly moving toward Physics-Based 
PTHA (PB-PTHA), which mandates the use of dynamic rupture models as 
 


---

# PAGE 38

34 
source inputs because they ensure mechanical consistency across all sampled 
scenarios (Sara Aniko Wirp et al., 2020; Ulrich etal., 2024). By incorporating 
initial stress and friction, dynamic rupture provides a self-consistent, 
mechanically viable range of rupture descriptions, essential for finding hazard 
curves on realistic seismic processes (Li et al., 2023; Gabriel et al., 2023). 
Concurrently, the propagation component must also be high-fidelity, as 
simplified, uncoupled methods (Method 4) introduce significant errors by 
neglecting non-hydrostatic effects in the near-field (Abrahams et al., 2023). 
Therefore, a robust PTHA requires the use of highly resolved, non-linear 
simulations, such as three-dimensional (3D) models or the fully coupled 
approach, to accurately capture phenomena like run-up and hydrodynamic 
loading from the dynamic source to the final inundation metric (Baragamage 
and Wu, 2024) 
 
2.6 Computational Efficiency in Hazard Modeling 
 
 
The necessity of adopting Physics-Based PTHA (PB-PTHA) with fully 
coupled, high-resolution models immediately introduces the computational 
bottleneck as the primary barrier to robust regional hazard assessment 
(Frontiers in Earth Science, 2020). Simulating the thousands of complex events 
required to accurately characterize uncertainty is often computationally 
prohibitive using conventional methodologies. 
 
2.6.1 The Curse of Dimensionality and Monte Carlo Limitations 
 
 
Traditional PTHA has historically relied on standard Monte 
Carlo (MC) methods to sample the uncertainty space, such as those 
used in regional studies by Bayraktar and Sozdinler (2020) and Kotani 
et al. (2020). However, the massive variability inherent in parameters 
like earthquake slip distribution leads to a "curse of dimensionality," 
where the simulation space is vast and highly multi-dimensional 
(Ramalingam et al., 2025). MC simulations suffer from a slow 
convergence rate of O(N^{-1/2}), where N is the number of scenarios, 
meaning that achieving stable, statistically robust results, especially for 
rare, high-consequence events that define the tail-end of the hazard 
curve requires an astronomically large ensemble size (Frontiers in Earth 
Science, 2020). This inherent inefficiency forces a compromise on the 
fidelity of the physical model, often leading researchers to use simpler 
linear wave theory or lower spatial resolution to reduce runtime, a 
trade-off that invalidates the PB-PTHA approach (NHESS, 2025). 
 
 


---

# PAGE 39

35 
2.6.2 Advanced Variance Reduction Strategies 
 
 
To overcome the compromise between accuracy and cost, the 
field has mandated the integration of efficient sampling strategies 
designed to reduce the variance of the hazard estimate without 
increasing the number of simulations. This approach provides the 
conceptual and methodological precedent for this study. The most 
widely adopted strategies fall into two categories. One involves 
creating Machine Learning (ML) surrogate models, such as variational 
encoder-decoders, which are trained on a limited, high-fidelity dataset 
of simulations (Ramalingam et al., 2025). Once trained, these surrogates 
can instantly predict key tsunami metrics (e.g., maximum wave height 
or inundation depth) for millions of untested scenarios, dramatically 
lowering the computational cost by up to 30 times and accelerating 
hazard quantification (Ramalingam et al., 2025). The other strategy 
focuses on variance-reduction techniques that directly improve the 
numerical stability and convergence rate of the sampling process itself, 
as demonstrated by studies using advanced statistical techniques to 
select and weight scenarios to cover the uncertainty space more 
intelligently than random sampling (Davies et al., 2022). This pursuit of 
superior convergence offered by methods beyond standard Monte 
Carlo is essential to make high-fidelity, physics-based PTHA 
computationally feasible. 
 
2.7 Randomized Quasi-Monte Carlo (RQMC) as a Methodological Solution 
 
 
This study proposes Randomized Quasi-Monte Carlo (RQMC) as the 
specific mathematical solution to the efficiency-accuracy trade-off detailed in 
Section 2.4, providing a means to achieve the robust statistical convergence 
required for high-fidelity PB-PTHA. RQMC addresses the critical flaw of 
traditional Monte Carlo by replacing random sampling with quasi-random 
sequences, which are designed to fill the sample space far more uniformly. 
 
2.7.1 Mathematical Rigor and Statistical Guarantee 
 
 
While standard Quasi-Monte Carlo (QMC) methods have long 
been recognized for their capacity to accelerate convergence, their 
application in high-stakes hazard modeling has been historically 
constrained by the inherent difficulty in accurately estimating the error 
of the hazard estimate. The lack of statistical rigor made it challenging 
 


---

# PAGE 40

36 
to define reliable confidence intervals, a necessity for safety-critical 
applications. However, this limitation has been overcome by the 
development of Randomized Quasi-Monte Carlo (RQMC). Nakayama 
and Tuffin (2024) recently provided the theoretical breakthrough by 
establishing sufficient conditions for Central Limit Theorems (CLTs) 
specifically applicable to RQMC estimators. This crucial mathematical 
rigor allows researchers to compute valid confidence intervals for the 
final hazard results, providing the essential statistical guarantee 
needed to deploy these methods in regional tsunami risk assessment. 
 
2.7.2 Superior Convergence in High-Dimensional Spaces 
 
Further supporting this methodological choice is RQMC's 
proven performance in complex, multi-parameter problems. Gobet et 
al. (2022) developed improved mean estimation techniques for RQMC, 
confirming its ability to maintain superior convergence rates. Unlike 
the O(N^-1/2) rate of standard MC, RQMC convergence can approach 
O(N^-1) (where N is the number of simulations), representing a 
substantial gain in efficiency. Crucially, this work validates that RQMC 
maintains this accelerated performance even in the high-dimensional 
uncertainty spaces characteristic of modern tsunami models, such as 
the stochastic slip distributions used in physics-based earthquake 
modeling, ensuring that efficiency is not sacrificed when incorporating 
complex source variability. 
 


---

# PAGE 41

37 
 
 
Figure 2.5. Convergence Rate Comparison  
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 


---

# PAGE 42

38 
 Chapter III 
METHODOLOGY 
 
 
 
3.1 Overall Methodology 
 
 
This study employs a computational–experimental research design 
centered on a Physical Fidelity-Based Probabilistic Tsunami Hazard 
Assessment (PTHA) Framework tailored to the Philippine Trench. The 
framework systematically links stochastic earthquake rupture modeling, 
geophysical probability theory, and high-resolution hydrodynamic simulation 
to produce probabilistic hazard metrics that are both physically realistic and 
statistically robust. The workflow is organized into four sequential phases, 
each serving a distinct purpose in the hazard assessment chain, from source 
characterization to final risk quantification. 
 
 
In Phase 1: Stochastic Slip Generation, this study synthesises a large 
ensemble of heterogeneous slip distributions on the subduction interface of 
the Philippine Trench. Earthquake slip is modeled as a log-normal Gaussian 
random field, with its mean anchored to empirical moment–slip scaling laws 
and its spatial covariance defined by a Matérn kernel parameterized with a 
correlation length of 20 km and a Hurst exponent of 0.3—values consistent 
with global slip inversion studies. To efficiently explore the high-dimensional 
rupture space (approximately 400 subfaults), this study implements 
Randomized Quasi–Monte Carlo (RQMC) sampling using Sobol’ sequences 
with Owen scrambling, which provides superior space-filling properties and 
faster convergence than standard Monte Carlo methods. Each raw slip 
realization 
is 
then 
adjusted 
through 
a 
moment-conserving 
“generate-and-scale” procedure to ensure strict adherence to the target 
seismic moment corresponding to specified magnitudes, including synthetic 
megathrust events (Mw 8.0–9.0) and historical ruptures such as the December 
2023 Mw 7.6 earthquake. 
 
 
Phase 2: Earthquake Probability Model establishes the statistical 
foundation necessary to transform physical scenarios into meaningful hazard 
estimates. This study adopts a Gutenberg–Richter recurrence model, 
parameterized using a and b values derived from the PHIVOLCS earthquake 
catalog and regional seismicity analyses, notably Bautista (2019). This model 
assigns an annual occurrence rate to each rupture scenario based on its 
magnitude, ensuring that frequent moderate events and rare catastrophic 
ruptures are weighted according to their natural likelihood. Without this 
probabilistic weighting, hazard assessments would conflate physical 
plausibility 
with 
frequency, 
potentially 
overemphasizing 
improbable 
megathrust events or underestimating more common but still damaging 
 


---

# PAGE 43

39 
tsunamis. Thus, this phase ensures that the final hazard metrics reflect not 
only what can happen, but also how likely it is to occur. 
 
 
Phase 
3: 
High-Fidelity 
Tsunami 
Simulation 
translates 
each 
moment-conserved slip realization into a dynamic tsunami inundation 
forecast. This study employs a fully nonlinear shallow-water solver enhanced 
with 
third-order 
Weighted 
Essentially 
Non-Oscillatory 
(WENO3) 
reconstruction to accurately capture steep wave fronts and near-shore flow 
complexity. Adaptive Mesh Refinement (AMR) is used to dynamically 
allocate computational resolution—starting at 1 km in the deep ocean and 
refining to 10 meters or finer in the Dapa coastal zone—thereby balancing 
computational efficiency with the need for urban-scale detail. The output for 
each scenario includes time-varying fields of water elevation, flow depth, and 
velocity across the entire domain, providing the physical basis for 
probabilistic synthesis. 
 
 
Phase 4: Hazard Aggregation integrates the ensemble of simulation 
results with the annual occurrence rates from Phase 2 to produce final 
probabilistic hazard metrics. For every location in the study area, this study 
computes the annual exceedance probability (AEP) of inundation depth by 
summing the occurrence rates of all scenarios whose simulated inundation 
exceeds a given threshold. This process yields spatially continuous hazard 
maps for standard return periods (e.g., 100-, 500-, and 2,500-year), site-specific 
hazard curves relating inundation depth to AEP, and exposure-weighted 
metrics such as the expected annual number of inundated buildings, derived 
using OpenStreetMap infrastructure data. Through this aggregation, 
thousands of deterministic simulations are transformed into a coherent 
probabilistic assessment that simultaneously accounts for source variability, 
physical wave dynamics, and seismic recurrence. 
 
 
Throughout this workflow, this study implements a dual-validation 
strategy to ensure methodological integrity. Statistical validation is performed 
by comparing the convergence behavior of RQMC against classical Monte 
Carlo, confirming that the sampling strategy achieves superior variance 
reduction with fewer realizations. Physical validation is conducted by 
comparing simulated sea surface time series against observed tide gauge 
records from the August 2012 and December 2023 (Mw 7.6) earthquakes, 
sourced from the IOC Sea Level Monitoring Facility and PHIVOLCS. 
Agreement in arrival time, amplitude, and waveform shape provides 
confidence in the hydrodynamic model’s fidelity. Together, these four phases 
and their embedded validation constitute a rigorous, end-to-end PTHA 
framework that delivers hazard estimates grounded in both geophysical 
reality and probabilistic reasoning providing a scientifically defensible 
foundation for coastal risk management in Dapa. 
 
 


---

# PAGE 44

40 
Below is the proposed framework of this study; 
 
 
Figure 3.1. Proposed Overall Framework 
 


---

# PAGE 45

41 
3.2 Preliminaries 
 
 
This section outlines the foundational inputs and methodological 
choices that underpin the physics-based probabilistic tsunami hazard 
assessment (PTHA) framework. It begins with a review of key advances in 
stochastic source modeling and sampling efficiency, followed by a summary 
of region-specific geophysical data used to constrain fault geometry and 
seismic recurrence. Together, these elements establish the scientific and 
computational basis for the stochastic slip generation and hazard 
quantification stages described in subsequent phases. 
 
3.2.1 Methodological Foundations 
 
 
This study establishes its methodological foundation through a 
structured 
examination 
of 
prior 
probabilistic 
tsunami hazard 
assessment (PTHA) frameworks, with emphasis placed on documented 
implementation 
procedures 
rather 
than 
reported 
results. 
Methodological elements were identified based on their reproducibility, 
computational feasibility, and relevance to subduction-zone tsunami 
generation, forming the basis for the framework adopted in this work. 
 
 
For stochastic earthquake source representation, established 
procedures for generating heterogeneous slip distributions were 
examined. Previous studies describe the construction of slip fields by 
discretizing the fault plane and sampling slip as a spatially correlated 
random process. Gaussian random field formulations with empirically 
calibrated covariance structures are commonly implemented to control 
spatial continuity and roughness across subfaults. Guided by these 
practices, this study adopts a Matérn-type covariance model, allowing 
correlation length and roughness parameters to be prescribed within 
ranges reported in global finite-fault inversion analyses, thereby 
ensuring physically plausible rupture realizations. 
 
 
The treatment of tsunami generation was informed by 
methodological approaches that couple earthquake rupture with 
hydrodynamic simulation. Prior implementations demonstrate the 
translation 
of 
slip 
distributions 
into 
time-dependent 
seafloor 
deformation directly embedded within nonlinear shallow-water 
solvers. This study follows such approaches by employing dynamic 
seafloor forcing rather than static initial conditions, enabling nonlinear 
wave generation and near-field effects to be explicitly resolved within 
the tsunami model configuration. 
 


---

# PAGE 46

42 
To address the high dimensionality associated with stochastic 
rupture ensembles, sampling strategies reported in the literature were 
evaluated based on convergence efficiency and variance reduction. 
Randomized 
Quasi–Monte 
Carlo 
(RQMC) 
techniques 
using 
low-discrepancy Sobol’ sequences with scrambling have been shown to 
provide uniform coverage of parameter space while maintaining 
unbiased estimators. Consistent with these implementations, this study 
applies RQMC sampling to efficiently explore the discretized rupture 
parameter space, reducing the number of simulations required for 
stable probabilistic estimates. 
 
 
Through the integration of these methodological components, a 
coherent workflow is constructed that links stochastic rupture 
generation, dynamic tsunami simulation, and probabilistic sampling 
within a single framework. The resulting methodology reflects 
established implementation practices while being adapted to the 
tectonic characteristics and data constraints of the Philippine Trench, 
providing the foundation for the phase-based procedures described in 
the succeeding sections. 
 
3.2.2 Collection of Fault Geometry and Seismotectonic Parameters 
 
 
Fault geometry and seismotectonic parameters for the Central to 
Northern Philippine Trench (7°–12°N) were compiled through the 
integration of published geophysical models and instrumental seismic 
datasets. 
The 
rupture 
domain 
was 
delineated 
by extracting 
trench-aligned slab geometry from the subduction interface model of 
Heidarzadeh et al. (2025), which was developed using seismic 
tomography, earthquake hypocenter distributions, and geodetic 
constraints. Strike, dip, and depth profiles were digitized from the slab 
surfaces and interpolated along the trench to construct a continuous 
fault representation. 
 
 
Depth-dependent shear modulus (μ) values were obtained by 
converting seismic velocity profiles from the Preliminary Reference 
Earth Model (PREM) and available regional crustal velocity models 
into elastic parameters. These values were assigned to corresponding 
depth intervals of the fault plane to support seismic moment 
calculations. 
 
 
Seismicity parameters were derived from the PHIVOLCS 
earthquake catalog, which provides instrumental records of regional 
 


---

# PAGE 47

43 
earthquakes. Event magnitudes and occurrence times were processed 
to 
construct 
frequency–magnitude 
distributions, 
from 
which 
Gutenberg–Richter a- and b-values were estimated using standard 
statistical fitting techniques. Together, these data collection and 
processing steps ensured that the stochastic earthquake source model 
was grounded in observed regional tectonics and seismic behavior. 
 
 
3.3 Phase 1: Stochastic Slip Generation and Distribution sampling 
 
 
In this phase, the study generates heterogeneous earthquake slip 
distributions that serve as the primary sources for the ensuing tsunami 
simulations. The goal of this phase is to produce a large ensemble of 
statistically consistent rupture scenarios that incorporate both spatial 
variability and the physical constraints imposed by the target moment 
magnitude. This is accomplished through a structured process involving 
random field modeling, low-discrepancy sampling, and moment-constrained 
scaling. 
 
3.3.1 RQMC method 
 
 
To efficiently explore the high-dimensional uncertainty inherent 
in 
earthquake 
slip 
distributions—typically 
represented 
by 
approximately 
400 
subfaults—this study adopted Randomized 
Quasi–Monte 
Carlo 
(RQMC) 
sampling 
based 
on 
Sobol’ 
low-discrepancy sequences combined with Owen scrambling. This 
choice was motivated by both theoretical properties of the sampling 
method and its demonstrated effectiveness in prior probabilistic 
tsunami and seismic hazard studies. 
 
 
Sobol’ sequences belong to a class of quasi-random point sets 
specifically engineered for high-dimensional numerical integration.  
They were selected primarily because they maintain superior 
space-filling properties even in dimensions well beyond 100, a regime 
where alternative methods such as Halton sequences or Latin 
Hypercube Sampling often degrade in performance. Given that the slip 
field in this study spans a 400-dimensional parameter space, this 
robustness is essential. Moreover, Sobol’ sequences exhibit a 
convergence rate that approaches O(N⁻¹), significantly faster than the 
O(N⁻¹ᐟ²) rate characteristic of classical Monte Carlo sampling. This 
accelerated 
convergence 
translates 
directly 
into 
computational 
efficiency: fewer realizations are required to achieve stable and 
reproducible estimates of tsunami hazard metrics, a critical advantage 
 


---

# PAGE 48

44 
when 
each 
simulation 
incurs 
substantial 
computational 
cost. 
Additionally, Sobol’ sequences are particularly well-suited to modeling 
irregular, rough functions—an attribute that aligns with the Matérn 
covariance structure used to define slip correlations. Because tsunami 
inundation is a highly nonlinear response to slip heterogeneity, the 
smoother estimator behavior and reduced sampling variance afforded 
by Sobol’ sequences enhance the reliability of probabilistic outputs. 
Their efficacy in this domain is further corroborated by a growing body 
of literature in probabilistic tsunami hazard assessment (PTHA) and 
stochastic source inversion, which consistently identifies Sobol’ 
sequences as a preferred tool for high-dimensional stochastic field 
generation due to their stability, reproducibility, and scalability. 
 
 
However, because classical Sobol’ sequences are deterministic, 
they do not permit statistical inference such as variance estimation or 
confidence interval construction—capabilities essential for quantifying 
uncertainty in hazard metrics. To address this limitation, this study 
incorporates Owen scrambling, a randomization technique that applies 
random digit permutations to the Sobol’ sequence while preserving its 
low-discrepancy structure. This transforms the deterministic sequence 
into a randomized quasi–Monte Carlo (RQMC) method that yields 
unbiased 
estimators. 
Crucially, 
Owen scrambling enables the 
computation of statistical error bounds across the ensemble, which is 
indispensable for assessing the convergence and robustness of hazard 
curves and inundation maps. Beyond enabling unbiased estimation, 
scrambling also improves the equidistribution of sample points by 
breaking up potential alignments between the slip field’s dominant 
spatial modes and the eigenvectors of the covariance matrix. Empirical 
studies have shown that RQMC with Owen scrambling can reduce 
sampling variance by one to two orders of magnitude compared to 
standard Monte Carlo, particularly when the underlying model 
exhibits strong nonlinearity—as is the case with tsunami generation 
and propagation, which involve discontinuities, shock formation, and 
steep nearshore gradients. This variance reduction enhances the 
consistency of results across independent ensemble runs and mitigates 
the risk of bias amplification in regions of high model sensitivity. 
 
 
     
Given the study’s objective—to generate physically realistic, 
moment-conserving slip realizations that are statistically representative 
yet computationally feasible—the combination of Sobol’ sequences and 
Owen scrambling emerges as the optimal sampling strategy. The slip 
field 
exhibits 
high 
dimensionality, 
strong 
spatial 
correlation, 
log-normal 
nonlinearity, 
and 
a 
Cholesky-based 
transformation 
structure, all feeding into a highly nonlinear hydrodynamic model. 
Under these conditions, the Sobol’–Owen RQMC approach provides 
the best balance of accuracy, efficiency, and statistical rigor. It ensures 
 


---

# PAGE 49

45 
rapid convergence of the PTHA ensemble while preserving the 
fine-scale variability necessary for moment-constrained rupture 
modeling.  
 
3.3.2 Slip Random Field Model 
 
 
This study begins by modeling earthquake slip as a log-normal 
Gaussian random field to reflect the natural heterogeneity of 
megathrust ruptures while ensuring non-negative slip values. The 
mean of this field is anchored to empirical slip–moment scaling laws, 
providing a physically informed baseline for median slip magnitude. 
Spatial correlations—critical for representing asperity structure—are 
encoded using a Matérn covariance kernel, parameterized with a 
correlation length of 20 km and a Hurst exponent of 0.3, values drawn 
from global subduction zone slip inversions. To enable efficient 
generation of correlated realizations, this study computes the Cholesky 
decomposition of the resulting covariance matrix. This decomposition 
serves as the mathematical bridge that transforms independent random 
inputs into spatially coherent slip patterns mimicking observed rupture 
complexity. 
 
 
Slip s is modeled as a log-normal Gaussian random field: 
 
 
 
 
where: 
            s              is the slip variable  
 
log           denotes the base-10 logarithm 
 
μs​             is the mean of the log-slip distribution 
 
1               is a vector of ones with appropriate dimension 
 
 L              is a linear operator (e.g., Cholesky factor of a  
covariance matrix) 
 
 ξ             is a vector of independent standard normal random  
variables 
 
 


---

# PAGE 50

46 
3.3.3 RQMC Sampling 
 
Building directly on the random field framework established in 
Section 3.3.2, this study employs Randomized Quasi–Monte Carlo 
(RQMC) sampling to generate the high-dimensional input vectors ξ 
required to drive the slip field model. Given that the rupture domain is 
discretized into approximately 400 subfaults, classical Monte Carlo 
sampling would be inefficient due to slow convergence and poor 
coverage in such a high-dimensional space. To overcome this, this 
study uses Sobol’ low-discrepancy sequences, which fill the unit 
hypercube more uniformly than pseudo-random draws. To retain the 
benefits of this uniformity while enabling rigorous statistical 
inference—such as unbiased variance estimation and confidence 
interval 
construction—they 
apply 
Owen 
scrambling, 
which 
randomizes the sequences without degrading their space-filling 
properties. Each resulting RQMC vector is then passed through the 
Cholesky-transformed Gaussian field defined in Section 3.3.2, yielding 
a spatially correlated, geophysically plausible slip realization. This 
study generates between 2,000 and 5,000 such realizations per 
magnitude scenario to ensure a robust statistical representation of 
rupture variability for downstream tsunami modeling. 
 
3.3.4 
Seismic 
Moment 
Calculation 
and 
Generate-and-Scale 
Correction 
 
Although the RQMC-driven Gaussian slip fields capture spatial 
heterogeneity and correlation structures informed by geophysical 
inversions, they do not inherently satisfy the total seismic moment 
required for a specified earthquake magnitude (e.g., Mw 8.0 or 
historical events like the December 2023 Mw 7.6 rupture). To enforce 
this fundamental physical constraint, each raw slip realization 
undergoes a moment-conserving “generate-and-scale” correction. The 
realized seismic moment is first computed as the sum of slip-weighted 
fault area multiplied by depth-dependent shear modulus: 
 
 
 
 
 


---

# PAGE 51

47 
 
 
where: 
 
         is the raw seismic moment estimate for event k 
𝑀𝑟𝑎𝑤
(𝑘)
         μ               is the shear modulus of the medium 
N               is the total number of contributing sub‐sources or elements 
ai               is the area associated with the i-th element 
           
          is the raw slip of the i-th element for event k 
𝑆𝑟𝑎𝑤,𝑖
(𝑘)
 
This is compared to the target moment derived from the 
standard moment–magnitude scaling relation: 
 
​ 
 
where: 
M0          is the seismic moment 
Mω         is the moment magnitude 
N            is the newton, the SI unit of force 
m            is the meter, the SI unit of length 
1.5         is the scaling factor to relate the magnitude to the logarithm 
of the seismic moment. 
9.1           is a constant used to calibrate Mω to be roughly equal to the 
older Richter scale in the magnitude range where they overlap. 
 
A uniform scaling factor 
 is then applied to the 
𝑎
(𝑘) =  
𝑀0
𝑀𝑟𝑒𝑎𝑙𝑖𝑧𝑒𝑑
(𝑘)
entire slip field, yielding a corrected realization 
 ​. This 
𝑠
(𝑘) =  𝑎
(𝑘)𝑠𝑟𝑎𝑤
(𝑘)
 


---

# PAGE 52

48 
step ensures that every synthetic rupture releases precisely the amount 
of energy implied by its target magnitude, preserving the nonlinear 
physics of earthquake source scaling. A final verification confirms that 
the relative moment error satisfies  
, guaranteeing 
𝑀0− 𝑀𝑟𝑒𝑎𝑙𝑖𝑧𝑒𝑑
|
|
𝑀0
 < 10
−6
consistency across the ensemble and producing tsunami initial 
conditions that are both statistically representative and physically 
valid. 
 
3.4 Phase 2: Probability Model 
 
The probability model serves as the statistical backbone of the PTHA 
framework, quantifying the likelihood of earthquake occurrences along the 
Philippine Trench. This model is essential because it transforms the raw 
simulation outputs into meaningful hazard metrics by weighting each 
scenario according to its probability of occurrence in nature. Without this 
probabilistic weighting, all earthquake scenarios—from frequent moderate 
events to extremely rare megathrust ruptures—would be treated as equally 
likely, leading to unrealistic hazard estimates that either underestimate 
common risks or overemphasize catastrophic but improbable events. 
 
3.4.1 Selection of the Gutenberg-Richter Recurrence Relation 
 
This study adopts the classical Gutenberg-Richter (G-R) 
recurrence law as the foundational model for earthquake frequency 
distribution along the Philippine Trench. This choice was driven by 
several compelling factors rooted in both empirical evidence and 
theoretical considerations. 
 
First, the G-R law has been validated globally across diverse 
tectonic 
settings 
over 
the 
past 
eight 
decades, 
consistently 
demonstrating that earthquake occurrence follows an exponential 
decay pattern with magnitude. This universal applicability made it the 
natural choice for a region like the Philippine Trench, where local 
seismic data—though improving—remains limited compared to 
well-instrumented subduction zones such as Japan or Cascadia. 
Second, the G-R framework is explicitly recommended in international 
tsunami hazard assessment guidelines, including those published by 
UNESCO-IOC and adopted in recent regional PTHAs for the 
Mediterranean Sea, Nankai Trough, and Southeast Asia. This study 
 


---

# PAGE 53

49 
therefore ensured methodological consistency with established best 
practices by selecting this model. 
 
Third, and most critically for this study, the G-R law directly 
addresses the fundamental challenge of probabilistic hazard analysis: 
how to rationally balance the contribution of frequent, smaller 
earthquakes against rare, catastrophic ones. In the context of tsunami 
hazard, this balance is particularly important because while magnitude 
7.0–7.5 events occur relatively often along the Philippine Trench (as 
evidenced by the 2012 and 2023 tsunamis), their tsunami impacts are 
localized and manageable. In contrast, a magnitude 9.0 megathrust 
event—though 
potentially 
centuries 
away—could 
generate 
coastline-wide devastation exceeding 17 meters of run-up, as 
demonstrated by the 2011 Tohoku and 2004 Indian Ocean tsunamis. 
The G-R law provides the mathematical structure to integrate both 
types of events into a single, coherent hazard estimate. 
 
The probability density function adopted by this study is 
expressed as: 
 
v(m) = β exp(−β m) 
 
where: 
v(m)   is the relative frequency or occurrence rate associated with 
earthquakes of magnitude m 
 
β    is the seismicity parameter controlling how rapidly 
earthquake frequency decreases with increasing magnitude 
 
This study emphasize that this formulation differs subtly but 
importantly from the original Gutenberg-Richter log-linear form (log₁₀ 
N = a − bM). The exponential form used here is the continuous 
probability density equivalent, where β is directly related to the 
traditional b-value through the conversion β = b · ln(10) ≈ 2.303b. This 
continuous 
formulation 
is 
computationally 
advantageous 
for 
integration into the PTHA framework because it allows for smooth 
 


---

# PAGE 54

50 
weighting across the magnitude range without requiring discrete 
binning, thereby avoiding artificial discontinuities in the hazard curves. 
 
3.4.2 Seismicity Parameter Estimation for the Philippine Trench 
 
The determination of an appropriate β-value for the Philippine 
Trench was a critical methodological decision that required careful 
synthesis of regional seismic data, published literature, and tectonic 
context. This study conducts a systematic analysis of the historical 
earthquake catalog maintained by the Philippine Institute of 
Volcanology and Seismology (PHIVOLCS), covering seismic events 
recorded along the trench system from 1900 to 2024. 
 
This study applies the maximum likelihood estimation (MLE) 
method to fit the Gutenberg-Richter relation to the declustered catalog. 
This analysis revealed a b-value of approximately 0.95 ± 0.12 for 
earthquakes in the Mw 5.0–7.5 range, which is consistent with typical 
subduction 
zone 
seismicity 
worldwide. 
However, 
this 
study 
recognized a critical limitation: the historical catalog for the Philippine 
Trench contains only two instrumentally recorded tsunamigenic 
megathrust events (the 2012 and 2023 Mw 7.6 earthquakes), and no 
observed events exceeding Mw 8.0. This data sparsity at the 
high-magnitude end introduces significant epistemic uncertainty when 
extrapolating recurrence rates for Mw 8.5–9.0 scenarios, which are 
precisely the events most relevant to catastrophic tsunami hazard. 
 
To address this uncertainty while maintaining physical realism, 
this study adopts a conservative β-value of 0.9 for the final probability 
model. This value was selected through a multi-criteria decision 
process that balanced three key considerations: 
 
1. Regional Seismic Precedent: While the Philippine Trench 
itself lacks recent Mw ≥ 8.0 events, geomorphological evidence from 
uplifted marine terraces along Davao Oriental (Ramos et al., 2012) 
indicates that the southern segment has generated M ≥ 8.0 megathrust 
earthquakes during the Holocene. These prehistoric ruptures, dated 
between 8,080 and 4,140 cal yr BP, suggest that the trench is capable of 
producing great earthquakes on millennial timescales. A β-value of 0.9, 
which corresponds to b ≈ 0.39 in the traditional formulation, implies a 
relatively "flat" magnitude-frequency distribution—meaning that large 
 


---

# PAGE 55

51 
earthquakes are not exponentially suppressed compared to moderate 
ones. This parameterization is appropriate for a subduction zone with 
demonstrated 
great-earthquake 
potential 
but 
limited 
modern 
instrumental records. 
 
2. Comparative Analysis with Analogous Subduction Zones: 
This study conducted a literature review of published β-values from 
tectonically similar subduction systems. The Nankai Trough, Japan (β ≈ 
0.85–0.95; Baba et al., 2022), the Mediterranean subduction zones (β ≈ 
0.90–1.05; Sørensen et al., 2012), and the Java Trench, Indonesia (β ≈ 
0.88–0.98) 
all 
exhibit 
β-values 
in 
the 
0.85–1.0 
range 
for 
megathrust-capable segments. The Philippine Trench shares key 
tectonic characteristics with these systems: high convergence rates (>10 
cm/yr), mature subduction geometry, and proximity to major 
population centers. The study therefore judged that β = 0.9 falls within 
the physically plausible range for a subduction margin with high 
seismic coupling and documented tsunamigenic history. 
 
3. Sensitivity Analysis and Conservative Risk Estimation: this 
study performed a sensitivity analysis by computing hypothetical 
hazard curves using β-values ranging from 0.7 to 1.1. This analysis 
revealed that lower β-values (β < 0.85) produce hazard estimates 
dominated by Mw 8.5–9.0 events, resulting in predicted 500-year wave 
heights exceeding 15 meters even in moderately exposed coastal 
locations—a projection that, while physically possible, lacks empirical 
validation from the regional tsunami record. Conversely, higher 
β-values (β > 1.0) suppress the contribution of great earthquakes to the 
extent that the resulting hazard curves fall below the observed impacts 
of the 2023 Mw 7.6 event, which clearly underestimates risk. The value 
β = 0.9 emerged as the optimal balance: it produces hazard estimates 
consistent with historical tsunami observations while appropriately 
incorporating the low-probability, high-consequence tail represented by 
potential Mw 9.0 ruptures. This conservative approach aligns with the 
precautionary principle in disaster risk reduction, particularly given 
the Philippines' high coastal population density and limited tsunami 
preparedness infrastructure. 
 
This study applied the β = 0.9 parameter across a magnitude 
range of Mw 7.0 to Mw 9.0 in discrete bins of ΔM = 0.5. This range was 
chosen based on tsunami physics: earthquakes below Mw 7.0 rarely 
generate significant basin-wide tsunamis due to limited seafloor 
displacement, while events exceeding Mw 9.0 are considered physically 
implausible for the Philippine Trench given its finite rupture length 
 


---

# PAGE 56

52 
(~1,300 km) and observed plate coupling geometry. The study therefore 
focused computational resources on the magnitude window most 
relevant to actionable hazard assessment for Philippine coastal 
communities. 
 
3.4.3 Integration into Hazard Aggregation 
 
The probability model is directly integrated into the final hazard 
quantification 
through 
the 
total 
probability 
theorem, 
which 
mathematically combines the earthquake recurrence rates with the 
conditional tsunami probabilities derived from the 40,000 stochastic 
simulations. The recurrence rate v(m) acts as a weighting function in 
the PTHA integral:      
 
where: 
P(H>h|x)   is the total (unconditional) probability that the tsunami 
intensity measure H exceeds a threshold value h 
 
H                  is the tsunami intensity measure of interest 
  
h​
         is the specified threshold value of the tsunami intensity 
measure H 
 
m ​
          is the earthquake magnitude 
 
dm ​
          is the differential magnitude increment used in the  
continuous formulation of the total probability integral 
 
In applying this formulation, this study relies on the total 
probability theorem to rigorously connect the frequency of earthquake 
occurrences with the tsunami heights they are capable of generating. 
The theorem states that when an outcome, such as tsunami height, can 
result from multiple possible earthquake magnitudes, the total 
probability of that outcome must incorporate both how often each 
magnitude occurs and how strongly it contributes to the resulting 
hazard. This means that the final tsunami hazard cannot be determined 
by simulation outputs alone; it must be weighted by the real-world 
likelihood of the events that produced those outputs 
 
 


---

# PAGE 57

53 
By 
integrating 
recurrence 
rates 
with 
simulation-derived 
exceedance probabilities, this theorem ensures that the hazard estimate 
reflects the true contributions of the full magnitude spectrum. Frequent 
moderate earthquakes (Mw 7.0–7.5) significantly influence the hazard at 
short return periods (50–100 years), while rare but catastrophic Mw 
8.5–9.0 megathrust ruptures dominate the long-period hazard 
(500–2,500 years). Without applying the total probability theorem, the 
assessment would either overemphasize rare extreme events or 
overlook the cumulative effect of more common but smaller 
earthquakes. 
 
Through this integration, the ensemble of modeled tsunami 
heights is transformed into physically meaningful and statistically 
robust exceedance curves. This step provides a scientifically defensible 
foundation for comparing the resulting hazard levels against 
engineering design standards, evacuation planning thresholds, and 
national tsunami risk metrics. 
 
3.5 Phase 3: High-Fidelity, Fully Coupled Tsunami Simulation 
 
To ensure physical realism in tsunami generation and propagation, this 
study adopts the Fully Coupled Dynamic approach (Method 1) as the baseline 
modeling strategy. Unlike static or uncoupled methods that impose an 
instantaneous 
seafloor 
uplift, 
this 
approach 
explicitly 
resolves the 
time-dependent momentum transfer from the rupturing fault into the 
overlying water column—a critical process that governs near-field wave 
amplification and high-frequency content. This study implemented this 
coupling through a modified version of the GeoClaw solver, which they 
extended to simultaneously integrate the elastodynamic equations for 
solid-Earth deformation and the nonlinear shallow-water equations for ocean 
dynamics. 
This 
dual-system 
integration 
enables 
the production of 
dynamically consistent tsunami initial conditions that reflect the true 
kinematics of megathrust rupture. 
 
The successful execution of this high-fidelity simulation framework 
was made possible through three interdependent innovations in numerical 
modeling, source implementation, and high-performance computing—each 
meticulously designed and validated by the study. 
 


---

# PAGE 58

54 
3.5.1 Numerical Model: Modified GeoClaw 
 
This study developed a customized GeoClaw solver that solves 
two core systems of governing equations. First is the elastic wave 
equation, which computes the time-dependent vertical seafloor 
displacement  uz(x,y,t) driven by slip evolution; Second is the nonlinear 
shallow-water equations, which govern the tsunami’s free-surface 
elevation  η(x,y,t) and depth-averaged velocity components  (u,v). 
 
To preserve physical fidelity while maintaining numerical 
stability across disparate spatial scales, four key enhancements were 
implemented: (i) third-order WENO reconstruction to sharpen 
wavefronts 
and 
suppress 
spurious 
oscillations 
near 
steep 
gradients—such 
as during run-up—while minimizing artificial 
diffusion; (ii) bathymetry-aware adaptive mesh refinement (AMR), 
triggered where ∣∇b∣>0.01, enabling automatic resolution refinement 
from 1 km in the open ocean to ≤30 m in coastal zones; (iii) 
well-balanced numerical fluxes that exactly preserve the “lake-at-rest” 
condition over real-world bathymetry (GEBCO/ETOPO), eliminating 
spurious waves in long-distance propagation; and (iv) stable time 
integration using second-order Runge–Kutta (RK2) time stepping with 
a CFL number of 0.5, balancing temporal accuracy and computational 
efficiency. 
 
3.5.2 Source Implementation 
This study designs a kinematically realistic earthquake source 
model that directly feeds into the coupled solver. Each rupture 
realization is represented as a time-dependent slip field(x,y,t), with 
rupture rise time scaled at 1.2 seconds per kilometer of rupture 
length—a value consistent with global kinematic inversions. To convert 
this slip history into dynamic seafloor motion, the study implements an 
approximation of 3D elastodynamic Green’s functions using a 
kinematic time-history approach following Lotto et al. (2018). This 
yields a physically consistent vertical displacement field uz​ (x,y,t) that 
evolves in time and couples directly to the ocean solver. 
 
For every simulation, this study configures the model to output 
three critical products for hazard analysis: 
 
 


---

# PAGE 59

55 
●​ Tide-gauge time series at virtual and real gauge locations; 
 
●​ Maximum offshore and nearshore wave heights; 
 
●​ Spatial maps of maximum inundation depth, written at native AMR 
resolution for downstream probabilistic aggregation. 
 
These outputs form the physical backbone of the hazard metrics 
computed in Phase 4. 
 
3.6 Dual Validation Framework 
 
To demonstrate that this study’s framework is both statistically efficient 
and physically accurate, a comprehensive two-part validation strategy is 
employed. This dual approach ensures that the methodology not only 
converges reliably in a probabilistic sense, producing stable and reproducible 
hazard estimates with minimal computational waste, but also faithfully 
reproduces the observable physical behavior of real tsunami events. 
 
3.6.1 Statistical validation — convergence efficiency. 
 
This study benchmark Randomized Quasi-Monte Carlo (RQMC) 
against baseline Monte Carlo by tracking the Coefficient of Variation 
(CoV) of the maximum wave height at the key site (Dapa). RQMC uses 
scrambled Sobol’ sequences (Owen scrambling) to achieve improved 
space-filling of the stochastic slip parameter space; theoretically and 
empirically, RQMC attains a far faster CoV decay (approaching O(N−1)) 
compared to MC’s O(N−1/2), enabling tighter confidence intervals for a 
fixed computational budget. This study quantifies convergence by 
computing CoV as a function of sample size N across replicates (R=10) 
and compares observed decay rates to the theoretical expectations. The 
RQMC sampling scheme and ensemble sizes are described in Section 
3.2 and the sampling rationale is given in the RQMC discussion. 
3.6.2 Physical validation — waveform fidelity. 
 
 


---

# PAGE 60

56 
This study evaluates physical realism by reproducing the 2012 
and 2023 Mw 7.6 events and comparing simulated to observed 
tide-gauge records. Waveform similarity is quantified using the 
Normalized Root-Mean-Square Error (NRMSE), defined as the L2​ norm 
of the simulated-minus-observed trace normalized by the L2​ norm of 
the observation, and by measuring absolute arrival-time error 
. The study’s target thresholds (NRMSE < 0.3 and 
∆𝑡 =  𝑡𝑠𝑖𝑚 −𝑡𝑜𝑏𝑠 
|
|
) were chosen to ensure realistic timing and amplitude 
∆𝑡 <  2 𝑚𝑖𝑛
reproduction for coastal early-warning and engineering use cases. The 
study 
then 
compares 
these 
metrics 
against 
equivalent 
Okada-plus-shallow-water baseline runs to demonstrate the improved 
fidelity that fully coupled, high-order numerics provide. See Section 3.5 
for formulas and benchmark criteria. 
 
3.7 Hazard Aggregation  
 
After producing the ensemble of coupled tsunami simulations, this 
study aggregates the results into probabilistic hazard products using a 
discrete approximation of the total probability theorem. The aggregation 
integrates stochastic rupture realizations with seismicity rates derived from a 
Gutenberg–Richter parameterization (central estimate β = 0.9) over the 
operational magnitude window Mw ∈ [7.0, 9.0], binned at ΔM = 0.5. This 
framework enables cell-wise computation of annual exceedance probabilities 
across the full raster domain, forming the basis for spatially explicit tsunami 
hazard assessment. 
 
3.7.1 Annual Probability of Exceedance 
 
Concretely, for each grid cell (i.e., each coordinate in the 
computational domain) this study computes the Annual Probability of 
Exceedance (APE) for an intensity threshold by summing the weighted 
contributions of every magnitude bin and simulation realization. 
Operationally this is: 
 
 
 
 


---

# PAGE 61

57 
 
where: 
P(H>h|x)    is the estimated conditional probability that the variable H 
exceeds the threshold h, given x 
H               is the hazard (or response) random variable 
h                is the exceedance threshold 
x                is the conditioning variable or set of observed parameters 
M              is the set of models considered 
m               indexes a specific model in the set M 
v(m)             is the weight or prior probability assigned to model m 
Nm             is the number of realizations (or simulations) generated  
for model 
m 
k                   indexes the realizations for a given model 
H(k)(m)         is the value of H for the k-th realization under model m 
 
1{⋅}               is the indicator function, equal to 1 if the condition inside 
the braces is satisfied and 0 otherwise 
 
This discrete implementation and each component’s role are 
defined and discussed in Section 3.6. 
 
The study derives v(m) from the chosen Gutenberg–Richter 
parameterization 
(central 
estimate 
β=0.9) 
applied 
across 
the 
operational magnitude window Mw ∈ [7.0,9.0] with ΔM=0.5 bins. This 
study emphasises that β=0.9 is treated as a central, defensible estimate 
chosen after sensitivity testing; future work can fold β into a logic-tree 
to capture epistemic uncertainty. The reasoning for the β choice and the 
sensitivity analysis appear in Section 3.3.3. 
 
 


---

# PAGE 62

58 
This study implements the aggregation cell-by-cell over the 
entire raster domain: for each cell the study iterates across the 40,000 
realizations, counts exceedances for each magnitude bin, multiply by 
the bin’s annual rate, and sum to obtain the local APE. The method is 
fully parallelizable (map-reduce style) and was implemented in the 
study’s Python postprocessing stack (xarray, rioxarray, geopandas) to 
scale across the Zarr/NetCDF output shards. The algorithmic steps 
and an illustrative worked example (Site X, h=3.0) are given in Section 
3.6.3–3.6.4 and show how a frequent-but-moderate (e.g., Mw 8.0) and a 
rare-but-extreme (e.g., Mw 8.5) scenario are balanced numerically to 
produce a single annual exceedance probability. 
 
 3.7.2  Outputs and Evaluation 
 
From the aggregated APE fields we produce three primary, 
actionable products: (1) site-specific tsunami hazard curves (APE vs. 
wave height) for engineering return-period lookups, (2) inundation 
probability maps such as P(inundation ≥ 0.5 m) for spatial planning, 
and 
(3) 
building 
exposure 
counts 
obtained 
by 
intersecting 
high-probability footprints with OSM building footprints (Dapa). We 
validate these outputs against our evaluation metrics (CoV decay, 
NRMSE, Δt, moment conservation tolerance <10−6, and CI widths) to 
ensure the results are both statistically robust and physically 
meaningful. See Sections 3.6.5 and 3.7 for full details and target 
thresholds. 
 
3.8 Experimental Procedure 
 
The experimental procedure is conducted in a systematic, phased 
manner to ensure reliable and comparable results across both statistical and 
physical 
validation 
domains. 
The 
prototype 
framework 
integrates 
RQMC-based stochastic slip generation with the Enhanced Fully Coupled 
Method (Method 1) and will be executed across multiple magnitude bins to 
produce a comprehensive ensemble of 40,000 tsunami simulations. 
 
3.8.1 Simulation Campaign Structure 
 
 
Each magnitude bin (Mw 7.0, 7.5, 8.0, 8.5, and 9.0) will undergo 
10,000 stochastic slip realizations, ensuring adequate statistical 
 


---

# PAGE 63

59 
coverage of the uncertainty space. For RQMC validation, 10 
independent 
replicates 
(R=10) 
will 
be 
generated 
using 
Owen-scrambled Sobol' sequences, with sample sizes varying across 
[100, 250, 500, 1000, 2500, 5000, 10000] to enable convergence analysis. 
A parallel baseline Monte Carlo ensemble will be generated using 
identical sample sizes to establish comparative benchmarks for 
convergence efficiency. 
 
For physical validation, historical event reconstructions will be 
performed by configuring the slip model to match the documented 
rupture parameters of the 2012 and 2023 Mw 7.6 Philippine Trench 
earthquakes. Simulated tide gauge waveforms will be extracted at the 
exact coordinates of observational stations (Guiuan, Mati, Surigao) and 
time-aligned with recorded data to compute NRMSE and arrival time 
error metrics. 
 
3.8.2 Computational Workflow 
 
 
Each simulation trial will follow a standardized five-stage pipeline: 
 
1.​ Slip Generation: A moment-constrained slip realization 
is 
𝑆
(𝑘)
synthesized 
using 
the 
Matérn 
covariance 
kernel 
(
) 
and 
scaled 
to 
satisfy 
𝑙 =  20𝑘𝑚,  𝑣 =  0. 3
. 
𝑀0 −𝑀𝑟𝑒𝑎𝑙𝑖𝑧𝑒𝑑
|
| 𝑐𝑖𝑡𝑒𝑠𝑡𝑎𝑟𝑡
[
]/𝑀0 < 10
−6
 
2.​ Probability Model Assignment: Each generated slip realization 
is 
assigned 
a 
probabilistic 
weight 
according 
to 
the 
Gutenberg-Richter recurrence relation 
, 
𝑣(𝑚) = β 𝑒𝑥𝑝(−β𝑚)
where 
=0.9. This weighting ensures that the stochastic 
β
ensemble reflects the relative frequency of earthquakes across 
the magnitude spectrum (Mw 7.0–9.0). 
 
3.​ Coupled Simulation: The Enhanced Fully Coupled Method 
propagates 
the 
tsunami 
using 
WENO3 
reconstruction, 
bathymetry-triggered AMR (
), and well-balanced 
∇𝑏
|
| > 0. 01
fluxes under a CFL=0.5 constraint. 
 


---

# PAGE 64

60 
 
4.​ Output Extraction: Maximum wave heights, arrival times, and 
inundation depths are recorded at 1,247 coastal observation 
points spanning the Caraga and Davao regions. 
 
5.​ Validation Check: For historical events, simulated waveforms 
are compared against tide gauge records; for synthetic events, 
moment conservation and numerical stability are verified. 
 
All simulation outputs will be stored in compressed Zarr/NetCDF 
format 
to 
facilitate 
parallel 
postprocessing 
and 
cloud-based 
aggregation. 
 
3.8.3 Trial Execution and Replication 
 
 
Results from each simulation trial, including maximum coastal 
wave height (
), coefficient of variation (CoV), NRMSE, arrival time 
𝐻𝑚𝑎𝑥
error (
), and computational runtime, will be recorded in structured 
∆𝑡
HDF5 databases for statistical analysis. Comparative evaluation will be 
performed to quantify the relative improvement of RQMC over 
baseline Monte Carlo in terms of convergence rate and confidence 
interval width. For physical validation, the Enhanced Method 1 will be 
benchmarked against equivalent simulations using static Okada 
displacement with standard shallow-water equations to demonstrate 
improvements in waveform fidelity. 
 
Statistical significance testing will be applied to convergence 
metrics using bootstrapped confidence intervals (95% CI) to verify that 
observed CoV decay rates are consistent with theoretical RQMC 
predictions 
(approaching 
). 
For 
waveform 
validation, 
𝑂(𝑁
−1)
acceptance thresholds (NRMSE < 0.3,
< 2 min) will be rigorously 
∆𝑡
|
|
enforced to ensure physical realism before proceeding to hazard 
aggregation. 
 
 
 
 
 


---

# PAGE 65

61 
3.8.4 Experimental Setup and Computing Environment 
 
 
Despite its methodological rigor, this study acknowledges 
several inherent limitations that may affect the generalizability and 
computational scalability of the findings. First, the research team 
conducted all simulations on an ASUS TUF Gaming A15 (2024) laptop 
equipped with an AMD Ryzen 9 8945H processor (8 cores, 16 threads, 
4nm Zen 4 architecture, base clock 4.0 GHz, boost up to 5.2 GHz) and 
16 GB DDR5-5600MHz RAM. While this configuration represents a 
capable 
mobile 
workstation 
with 
modern 
architecture 
and 
high-frequency memory, it operates under significantly different 
constraints compared to dedicated high-performance computing 
(HPC) infrastructure. Consequently, runtime performance is inherently 
limited by thermal throttling under sustained loads, single-node 
memory capacity, and the absence of distributed computing 
capabilities. The research team anticipates that wall-clock time per 
simulation 
will 
be 
substantially 
longer 
than 
HPC-based 
implementations, with estimates ranging from 2–4 hours for Mw 7.5–8.0 
scenarios and potentially 6–10 hours for Mw  9.0 cases due to thermal 
management constraints and extended AMR computations on a 
consumer-grade platform. 
 
Second, 
although 
the 
Philippine 
Trench study domain 
encompasses the primary tsunamigenic source for eastern Philippines, 
the framework does not account for secondary hazard amplifiers such 
as splay faults, submarine landslides triggered by seismic shaking, or 
outer-rise normal faulting events. These mechanisms could locally 
modify tsunami amplitudes but are considered beyond the scope of the 
present megathrust-focused assessment. Additionally, bathymetric and 
topographic data were sourced from GEBCO 2023 (15 arc-second 
resolution globally) and locally enhanced datasets with NAMRIA lidar 
data (5-meter resolution) for critical coastal zones. While this represents 
state-of-the-art coverage for the Philippines, fine-scale features such as 
coral reef structures, seawalls, or small drainage channels are not fully 
resolved and may introduce localized inundation prediction errors. 
 
A further limitation lies in the treatment of uncertainty. The 
Gutenberg-Richter parameter (
) was held constant as a central 
β = 0. 9
estimate, though sensitivity analysis demonstrates that β ϵ [0. 7,  1. 1] 
could shift 500-year wave height estimates by 
. Future 
± 30%
extensions incorporating logic tree frameworks to formally weight 
alternative seismicity models would provide more comprehensive 
epistemic uncertainty quantification. Additionally, while RQMC 
 


---

# PAGE 66

62 
significantly accelerates convergence relative to standard Monte Carlo, 
the 
ensemble 
size 
(N=10,000 
per 
magnitude 
bin) 
remains 
computationally constrained. Tail-risk estimates for extremely rare Mw 
9.0 scenarios may therefore exhibit wider confidence intervals than 
more frequent Mw 7.5–8.0 events. 
 
To mitigate environmental and numerical variability—including 
dynamic load balancing artifacts, minor floating-point rounding 
differences, and AMR grid structure stochasticity—all random number 
generation uses fixed seeds (
 for Sobol', 
 for 
𝑠𝑒𝑒𝑑= 42
𝑠𝑒𝑒𝑑= 2024
NumPy). This ensures bitwise-reproducible slip realizations across 
independent runs despite typical maximum wave height deviations of 
 cm. Lastly, the present implementation prioritizes tsunami 
< 2
hydrodynamics and elastic seafloor response, omitting atmospheric 
acoustic-gravity wave propagation and poroelastic effects in the 
shallow accretionary wedge. These phenomena have second-order 
effects on coastal inundation metrics and are reserved for future 
research. 
 
3.8.5 Software Implementation 
 
 
The software implementation is developed using a hybrid 
Fortran-Python computational stack, chosen for its balance of 
numerical performance and scientific productivity. The core tsunami 
solver is built upon GeoClaw v5.11.1, a finite-volume hydrodynamic 
code written in Fortran 90 that implements the augmented Riemann 
solver framework with built-in support for adaptive mesh refinement 
(AMR). Modifications to GeoClaw's source term routines incorporate 
WENO3 reconstruction for momentum fluxes and extend AMR criteria 
to trigger refinement based on bathymetric gradients (
) in 
∇𝑏
|
| > 0. 01
addition to standard wave height and velocity thresholds. 
 
Stochastic slip generation and moment-constrained scaling are 
implemented in Python 3.11 using NumPy 1.26 for linear algebra 
operations, SciPy 1.11 for Cholesky decomposition of the Matérn 
covariance matrix, and the chaospy 4.3 library for generating 
Owen-scrambled Sobol' sequences. The Matérn kernel computation 
leverages scipy.spatial.distance_matrix for efficient pairwise distance 
calculation across the 400-subfault grid, with correlation length 
 and roughness parameter 
 enforced throughout. 
𝑙= 20 𝑘𝑚
𝑣= 0. 3
Custom Python modules handle the generate-and-scale protocol, 
 


---

# PAGE 67

63 
iteratively adjusting slip fields until seismic moment conservation is 
satisfied to machine precision (relative error 
). 
< 10
−6
 
For benchmarking and validation, the baseline Monte Carlo slip 
generator uses NumPy's default Mersenne Twister pseudo-random 
number generator (np.random.default_rng(seed=2024)) to ensure 
reproducibility. Physical validation against the 2012 and 2023 historical 
events employs Pandas 2.1 for tide gauge data ingestion and 
time-series alignment, with NRMSE and arrival time error (
) 
∆𝑡
computed using vectorized NumPy operations. Postprocessing and 
hazard 
aggregation 
leverage 
xarray 
2023.10 
for 
handling 
multi-dimensional NetCDF outputs, rioxarray 0.15 for spatial 
reprojection and rasterization, and GeoPandas 0.14 for geographic 
overlay operations. Visualization outputs are generated using 
Matplotlib 3.8 and Generic Mapping Tools (GMT) 6.4, with final layout 
composition performed in QGIS 3.34.                                                      
3.9 Evaluation Metrics 
 
The performance and validity of the Physical Fidelity-Based 
Probabilistic Tsunami Hazard Assessment framework will be evaluated using 
a comprehensive dual-validation approach. The primary statistical metric is 
the Coefficient of Variation (CoV), previously defined in Section 3.5.1, which 
integrates convergence rate, confidence interval width, and ensemble stability 
into a single normalized measure of RQMC sampling efficiency. In this 
methodology, CoV will be used to compare all RQMC and Monte Carlo 
experimental runs across varying sample sizes and determine which sampling 
configuration provides the most computationally efficient uncertainty 
quantification. This study will complement this statistical validation with 
physical fidelity metrics, specifically the Normalized Root Mean Square Error 
(NRMSE) introduced in Section 3.5.2, to assess waveform accuracy against the 
2012 and 2023 historical tsunami events. These metrics ensure a balanced 
assessment by capturing both statistical convergence behavior and physical 
realism without isolating any single performance factor, thereby providing a 
rigorous foundation for validating the framework's ability to deliver 
actionable probabilistic hazard estimates for Philippine coastal communities. 
 
3.10 Initial Simulation 
 
The journey toward a physics-based probabilistic tsunami hazard 
assessment began with a single, carefully constructed earthquake slip 
 


---

# PAGE 68

64 
realization. This “initial simulation” served as a proof-of-concept—verifying 
that the full pipeline, from stochastic source generation to moment-conserving 
correction, functioned as intended before scaling to thousands of ensemble 
members. 
 
The scenario targeted a synthetic Mw 8.5 megathrust earthquake on the 
Philippine Trench, based on the hypothetical rupture defined in Heidarzadeh 
et al. (2025, Table 3). Using this reference, the code (slip_sampler.py) loaded 
the following geometric and physical parameters directly from a CSV file: 
 
Target Magnitude: Mw=8.5 
Rupture length: Lkm  = 300.0 km 
Rupture width: Wkm = 100.0 km 
Subfault size: subfault_size_km = 20.0 km 
Grid Size Discretization: nalong = 15, ndown = 5 → N = 75 subfaults 
Shear modulus: μ = 4.0 × 10¹⁰ Pa (40 GPa) 
Subfault area: ai = (20 × 10³)² = 4.0 × 10⁸ m²  
 
From these, the target seismic moment was computed using the 
canonical Hanks and Kanamori (1 979) relation: 
 
 
𝑀0 =  10
1.5𝑀𝑤+9.1 
=  10
21.85 ≈ 7. 08 ×  10
21𝑁𝑚.
 
Next, a stochastic slip field was generated. This study modeled slip as a 
log-normal random field (Equation 3-1), ensuring non-negative values while 
capturing the spatial clustering of high-slip asperities characteristic of real 
earthquake ruptures. Spatial correlation was imposed through a Matérn 
covariance kernel, governed by two key hyperparameters: 
 
Correlation length ℓ = 20 km (controlling the size of slip patches), 
 


---

# PAGE 69

65 
Smoothness (Hurst) parameter ν = 0.3. 
 
The choice of ν = 0.3 was not arbitrary. Empirical studies of global 
megathrust earthquakes—such as those by Mai & Beroza (2002) and Goda & 
Song (2016)—consistently report Hurst exponents in the range 0.2–0.4 for 
inverted slip distributions. A value of ν = 0.3 produces fields with moderate 
roughness: smooth enough to avoid unphysical spikes, yet rough enough to 
generate localized asperities that drive extreme tsunami waves. This 
parameter aligns with the observed k⁻² spectral decay of slip, a hallmark of 
self-similar earthquake rupture physics. 
 
To sample the high-dimensional uncertainty space (75 dimensions, one 
per subfault), the researchers employed Randomized Quasi–Monte Carlo 
(RQMC). Specifically, they used Sobol’ low-discrepancy sequences with Owen 
scrambling—a method that fills the parameter space more uniformly than 
random 
Monte 
Carlo 
while 
retaining 
statistical 
rigor. 
Crucially, 
reproducibility was guaranteed by fixing the randomization seed to 42 in all 
experiments: 
 
This ensures that every slip realization is bitwise reproducible across 
runs, a non-negotiable requirement for scientific validation, debugging, and 
HPC replication. At the same time, Owen scrambling preserves the 
low-discrepancy structure of Sobol’ sequences while enabling unbiased 
variance estimation—a critical feature for constructing confidence intervals in 
hazard curves. 
The raw slip field was then constructed via Cholesky decomposition of 
the Matérn covariance matrix and exponentiation of the correlated Gaussian 
vector. However, this raw field did not satisfy the target moment 
. Its total 
𝑀0
seismic moment, (Equation 3-2)​
 
was typically lower or higher than desired, as stochastic generation does not 
inherently enforce global constraints. 
 
 


---

# PAGE 70

66 
To resolve this, the researchers applied a global moment-conserving 
scaling factor: 
​
 
where: 
F            is the scaling factor applied to the raw slip values 
M0         is the target (desired) seismic moment 
M0,raw     is the raw seismic moment computed from the unscaled  
slip distribution 
 
and scaled every subfault uniformly: 
​
 
where: 
            is the scaled slip (or source parameter) for the i-th element  
𝑠𝑖
(𝑘)
                  in realization k 
 
           is the raw (unscaled) slip value for the i-th element 
𝑠𝑖
𝑟𝑎𝑤
F             is the scaling factor applied to the raw slip values 
I              indexes the spatial element (or sub‐fault) 
k              indexes the realization or simulation 
 
This global scaling is physically justified because it preserves the 
relative spatial heterogeneity of the slip field—maintaining the locations and 
amplitude ratios of asperities—while strictly enforcing the fundamental 
requirement that the total seismic moment matches the target magnitude, and 
it avoids the additional complexity and potential statistical biases associated 
with local or iterative rescaling methods that can distort spatial correlation; 
accordingly, the “generate-and-scale” approach is a standard and validated 
 


---

# PAGE 71

67 
practice in probabilistic tsunami hazard assessment, as demonstrated by Goda 
(2019) and Williamson et al. (2022). 
After scaling, the relative moment error was computed as: 
 
ε =
 
 
𝑀0− 𝑀𝑟𝑒𝑎𝑙𝑖𝑧𝑒𝑑
|
|
𝑀0
 < 10
−6
where: 
ε             is the relative error (dimensionless) 
M0          is the target or reference seismic moment 
M0 realized  is the realized (or computed) seismic moment 
∣⋅∣         denotes the absolute value 
10-6         is the prescribed tolerance indicating the acceptable error threshold 
The final slip field was saved as slip_Mw8.5_sample0.npy, and its 
statistics clearly indicated a successful realization, with a mean slip of 5.90 m 
and a maximum slip of 39.57 m, exhibiting strong spatial heterogeneity and 
well-defined asperities visible in the output plot; these values are consistent 
with the empirical median slip of 4.7 m for Mw 8.5 and reflect the expected 
upward adjustment resulting from global moment-conserving scaling. 
 
This initial simulation thus validated the entire stochastic source 
pipeline: 
from 
RQMC-driven 
spatial 
random 
fields, 
through 
moment-conserving correction, to physically consistent slip outputs. It 
confirmed that the framework could generate geostatistically realistic, 
physically valid, and reproducible earthquake scenarios—laying the essential 
groundwork for large-scale ensemble production and fully coupled tsunami 
simulation. 
 
 


---

# PAGE 72

68 
 
 
Figure 3.2. 
 
 
 
 
 
 
 
 
 
 
 
 


---

# PAGE 73

69 
3.11  Timeline of Research Activities  
 
 
 
PHASE 
TASK 
MONTH 
ACTIVITY 
DELIVERABLE 
PRELIMIN
ARIES 
Concept 
Defense & 
Literature 
Review 
October 
2025  
• Concept Defense  
• Finalize literature 
review  
• Weekly consultation 
with Dr. Landicho 
• Acquire fault geometry 
from Heidarzadeh et al. 
(2025) 
• Download bathymetry 
data (GEBCO 2023, 
NAMRIA lidar) 
• Process PHIVOLCS 
earthquake catalog 
(1900-2024) 
• Approved concept 
paper  
• Annotated 
bibliography 
• Philippine Trench 
database 
• Seismicity 
parameter report 
• Historical 
validation dataset 
Phase 1: 
Stochastic 
Slip 
Generation 
1.1 Finalize 
Slip Random 
Field Model 
November 
2025 (Week 
1-2) 
• Install GeoClaw v5.11.1 
and Python 3.11 
environment 
• Configure ASUS TUF 
Gaming A15 (thermal 
management) 
• Implement Matérn 
kernel (ℓ=20 km, ν=0.3) 
for 400-subfault grid 
• Set log-slip mean using 
moment-magnitude 
scaling 
• Code Cholesky 
decomposition 
 
• Configured 
computational 
environment 
• Validated slip 
covariance model 
• Matérn kernel 
implementation 
• Parameter set 
documentation 
•slip_field_model.py 
module 
 
1.2 
Implement 
RQMC 
Sampling & 
Generate 
Realizations 
November 
2025 (Week 
3-4) - 
December 
2025 (Week 
1-2) 
• Implement 
Owen-scrambled Sobol' 
sequences (chaospy 4.3, 
seed=42) 
• Develop 
"generate-and-scale" 
algorithm for moment 
conservation 
• Generate 10,000 slip 
•moment-conserved 
slip realizations 
stored in compressed 
Zarr format 
• RQMC_sampler.py 
module 
•moment_correction.
py module 
• Validation report 


---

# PAGE 74

70 
 
fields per magnitude bin 
(Mw 7.0–9.0, ΔM=0.5) 
• Enforce moment 
conservation (error < 
10⁻⁶) 
• Create baseline Monte 
Carlo generator 
(seed=2024) 
• Consultations with Dr. 
Landicho 
• Initial results ready 
for proposal 
 
Phase 2: 
Earthquake 
Probability 
Mode 
 
2.1 Calibrate 
Gutenberg-Ri
chter Model 
 
December 
2025 (Week 
2-3) 
 
• Fit β = 0.9 using 
PHIVOLCS catalog (Mw 
5.0–7.5) via MLE 
• Perform sensitivity test 
for β ∈ [0.7, 1.1] 
• Implement v(m) = β 
exp(−β m) 
• Create magnitude bins 
(Mw 7.0–9.0, ΔM = 0.5) 
• Integrate with slip 
generator 
• Prepare Proposal 
Defense 
• Consultations with Dr. 
Landicho 
 
• Finalized annual 
occurrence rates v(m) 
for Mw 7.0–9.0 
•GR_recurrence.py 
module 
• Sensitivity analysis 
report 
• Comparative 
analysis table 
• Proposal Defense 
Presentation 
 
 
 
 
 
 
 
PROPOSAL 
DEFENSE 
 
 
 
 
December 
16, 2025 
 
 
 
 
• Present research 
framework with initial 
results 
• Demonstrate 40,000 
slip realizations with 
RQMC sampling 
• Show moment 
conservation validation 
• Defend Enhanced 
Method 1 selection 
 
 
 
 
• Defense 
documentation 
• Panel 
recommendations 
document 
Phase 3: 
High-Fideli
ty Tsunami 
Simulation 
3.1 Configure 
Enhanced 
Fully 
Coupled 
Model 
December 
2025 (Week 
3-4) - 
January 
2026 (Week 
• Modify GeoClaw: 
integrate WENO3, 
bathymetry-aware AMR 
(∥∇b∥>0.01), 
well-balanced fluxes 
• Validated solver 
configuration on 
ASUS TUF A15 
• Enhanced GeoClaw 
with WENO3 + AMR 


---

# PAGE 75

71 
 
1-2) 
• Set CFL=0.5, RK2 time 
stepping 
• Implement kinematic 
slip time-history (1.2 
s/km) 
• Code elastodynamic 
Green's function for 
uz(x,y,t) 
• Configure Philippine 
Trench domain (7°-12°N) 
• Set up multi-resolution 
grid (1km → 10m) 
• Test OpenMP 
parallelization 
• Weekly consultations 
with Dr. Landicho 
• Fully coupled 
method 
implementation 
• Domain 
configuration files 
• Solver verification 
report 
• Performance 
benchmark 
 
3.2 Run 
Simulation 
Campaign 
January 
2026 (Week 
3-4) - March 
2026 
• Execute 40,000 
RQMC-driven 
simulations using 
OpenMP: 
- Jan W4-Feb: Mw 7.0 
(10,000 runs) 
- Feb-Mar W2: Mw 7.5 
(10,000 runs) 
- Mar W2-4: Mw 8.0 
(10,000 runs) 
- Mar W4-Apr W3: Mw 
8.5 (10,000 runs) 
- Apr W3-4: Mw 9.0 
(10,000 runs) 
• Store as Zarr/NetCDF, 
archive to 4TB external 
SSD 
• Monitor thermal 
performance (CPU < 
85°C) 
• Weekly progress 
consultations with Dr.  
Landicho 
• Complete 
simulation dataset 
(40,000 scenarios) 
• Water elevation, 
flow depth, velocity 
fields 
• Maximum wave 
height catalog 
• Inundation depth 
rasters 
• Tide gauge time 
series (virtual + real 
stations) 
• Computational 
performance log 
 
3.3 Physical 
Validation 
(2012 & 2023 
Events) 
April 2026 
(Week 1-2) 
• Reconstruct Aug 2012 
Mw 7.6 and Dec 2023 
Mw 7.6 events 
• Compare simulated vs. 
observed tide-gauge 
records (Guiuan, Mati, 
• Validation report 
showing NRMSE < 
0.3 and │Δt│ < 2 min 
• 2012 & 2023 
simulation results 
• Tide gauge 


---

# PAGE 76

72 
 
Surigao) 
• Compute NRMSE and 
arrival time error Δt 
• Analyze wave period 
differences 
• Compare against 
baseline Okada model 
• Consultations: May 7, 
14 
• Consultations with Dr. 
Landicho  
comparison plots 
• Waveform fidelity 
analysis 
• Method 1 
superiority 
demonstration 
Phase 4: 
Dual 
Validation 
& Hazard 
Aggregatio
n 
4.1 Statistical 
Validation of 
RQMC 
April 2026 
(Week 2-3) 
• Generate RQMC 
ensembles [100, 250, 500, 
1000, 2500, 5000, 10000] 
• Create R=10 replicates 
with Owen scrambling 
• Generate parallel 
Monte Carlo baseline 
• Compute Coefficient of 
Variation (CoV) at Dapa 
• Confirm RQMC 
convergence ~O(N⁻¹) vs 
MC ~O(N⁻⁰·⁵) 
• Convergence plots 
and performance 
comparison 
• CoV decay analysis 
(RQMC vs MC) 
• Statistical efficiency 
metrics 
• Confidence interval 
comparison 
• Statistical 
validation report 
 
4.2 
Preliminary 
Hazard 
Metrics 
April 2026 
(Week 4) 
• Compute exceedance 
counts per magnitude 
bin at key sites (Dapa, 
Mati, Guiuan) 
• Combine with v(m) for 
preliminary Annual 
Probability of 
Exceedance (APE) 
• Generate site-specific 
hazard curves 
• Extract OpenStreetMap 
building data for Dapa 
• Initial site-specific 
exceedance 
probabilities 
• Preliminary hazard 
curves for key 
locations 
• Building footprint 
intersection data 
• High-risk zone 
identification 
 
4.3 Full 
Hazard 
Aggregation 
May 2026 
(Week 1-2) 
• Implement total 
probability theorem 
• Generate probabilistic 
inundation maps (100-yr, 
500-yr, 2500-yr) 
• Calculate expected 
annual inundated 
buildings in Dapa 
• Perform β-value 
sensitivity analysis 
• Complete 
probabilistic hazard 
maps 
• Tsunami hazard 
curves (all locations) 
• Return period maps 
(1%, 0.2%, 0.04%) 
• Building exposure 
analysis 
•hazard_aggregation.


---

# PAGE 77

73 
 
py 
 
4.4 Document 
Methodology 
& Initial 
Findings 
May 2026 
(Week 2) 
• Compile detailed 
methodology chapter 
• Document algorithms, 
parameters, procedures 
• Write initial results 
chapter with validation 
• Create comprehensive 
figures and tables 
• Draft Chapter III 
(Methodology) - 
complete 
• Draft Chapter IV 
(Results) - initial 
• All methodological 
figures/tables 
• Algorithm 
documentation 
Phase 5: 
Final 
Analysis & 
Defense 
Results 
Analysis & 
Discussion 
May 2026 
(Week 3) 
 Complete analysis of all 
hazard outputs 
• Compare Enhanced 
Method 1 vs baseline 
Okada 
• Quantify RQMC 
efficiency gains 
• Analyze worst-case 
Mw 9.0 scenarios (17.4m 
waves) 
• Finalize all chapters 
(I-V) 
• Develop policy 
recommendations 
• Complete results 
analysis 
• Full thesis 
manuscript (Chapters 
I-V) 
• All figures, tables, 
and appendices 
• Policy 
recommendations 
document 
• Risk interpretation 
narrative 
 
Defense 
Preparation 
May 2026 
(Week 3) 
• Prepare defense 
presentation 
• Create simulation 
demonstration videos 
• Print thesis copies 
• Consultations with Dr. 
Landicho 
• Defense 
presentation 
• Demo videos 
• Printed 
manuscripts 
• Q&A preparation 
 
FINAL 
THESIS 
DEFENSE 
May 2026 
(Week 4) 
• Present 40,000 
simulation ensemble & 
PTHA framework 
• Show RQMC O(N⁻¹) 
convergence 
• Present validation 
(NRMSE<0.3, │Δt│<2 
min) 
• Showcase hazard maps 
& building exposure 
• Address panel 
questions 
• Panel approval and 
signatures 
• Final revision list  


---

# PAGE 78

74 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
Final 
Submission 
May 2026 
(Week 4) 
• Incorporate revisions 
• Archive code/data 
(GitHub) 
• Closing meeting with 
Dr. Landicho 
• Final manuscript 
• Code repository 


---

# PAGE 79

75 
References: 
 
​
[1] Abrahams, L. C., Dunham, E. M., Løvholt, F., & Sladen, A. (2023). 
Comparison of methods for coupled earthquake and tsunami modelling. 
Geophysical 
Journal 
International, 
234(1), 
404–420. 
https://doi.org/10.1093/gji/ggad089 
[2] Alhamid, A. A., Alazmi, H., Alabdan, A., Alfaraj, A. A., & Alshami, M. A. 
(2022). Probabilistic Tsunami Hazard Assessment (PTHA) for the Gulf of 
Aqaba and Red Sea Coast considering Sea Level Rise. Water, 14(12), 1957. 
https://doi.org/10.3390/w14121957 
[3] Alhamid, A. K., Akiyama, M., Ishibashi, H., Aoki, K., Koshimura, S., & 
Frangopol, D. M. (2022). Framework for probabilistic tsunami hazard 
assessment considering the effects of sea-level rise due to climate change. 
Structural Safety, 94, 102152. https://doi.org/10.1016/j.strusafe.2021.102152 
[4] Asai, M., Miyagawa, Y., Idris, N. A., Muhari, A., & Imamura, F. (2016). 
Coupled Tsunami Simulations… Journal of Earthquake and Tsunami. 
https://doi.org/10.1142/S1793431116400194 
[5] Baba, T., Satake, K., Kato, Y., & Matsumoto, H. (2022). Probabilistic 
tsunami hazard assessment along the eastern coast of Shikoku, Japan… 
Natural 
Hazards, 
112(1), 
201–229. 
https://doi.org/10.1007/s11069-021-05165-w 
[6] Baoning Wu, T., Hidayat, R., Kong, T., & Tregoning, P. (2025). Developing a 
Probabilistic Tsunami Hazard Assessment Framework… Environmental & 
Engineering 
Geoscience, 
31(1), 
67–85. 
https://doi.org/10.2113/EEG-2024-0025 
[7] Baragamage, D. S. P. A., & Wu, W. (2024). A three-dimensional 
fully-coupled 
fluid-structure 
model… 
Water, 
16(1), 
189. 
https://doi.org/10.3390/w16010189 
[8] Bayraktar, H. B., & Sozdinler, C. O. (2020). Probabilistic tsunami hazard 
analysis for Tuzla test site… Natural Hazards and Earth System Sciences, 20, 
1741–1764. https://doi.org/10.5194/nhess-20-1741-2020 
[9] Behrens, J., Davies, G., Geist, E., González, F. I., Horspool, N., Mofjeld, H., 
& Thio, H. (2021). Status of probabilistic tsunami hazard assessment 
(PTHA)… 
Pure 
and 
Applied 
Geophysics, 
178(2), 
685–722. 
https://doi.org/10.1007/s00024-020-02597-2 
[10] Cranley, R., & Patterson, T. N. L. (1976). Randomization of number 
theoretic methods for multiple integration. SIAM Journal on Numerical 
Analysis, 13(6), 904–914. 
 


---

# PAGE 80

76 
[11] Davies, G., Griffin, J., & Løvholt, F. (2022). Reducing the computational 
cost of PTHA… Natural Hazards and Earth System Sciences, 22(8), 2685–2701. 
https://doi.org/10.5194/nhess-22-2685-2022 
[12] Davies, G., Griffin, J., Wilson, K., & Burbidge, D. (2017). Probabilistic 
tsunami hazard assessment for earthquake sources in the global oceans. 
Marine Geology, 390, 1–15. https://doi.org/10.1016/j.margeo.2017.06.008 
[13] Davies, G., Weber, R., Wilson, K., & Cummins, P. (2022). From offshore to 
onshore probabilistic tsunami hazard assessment… Geophysical Journal 
International, 230(3), 1630–1651. https://doi.org/10.1093/gji/ggac140 
[14] Entacher, K. (1997). Quasi-Monte Carlo methods for numerical 
integration of multivariate Haar series. BIT, 37(4), 846–861. 
[15] Fan, J., & Zhao, D. (2018). Evolution of the southern segment of the 
Philippine Trench… Geochemistry, Geophysics, Geosystems, 19, 4612–4627. 
https://doi.org/10.1029/2018GC007685 
[16] Fine, I. V., Rabinovich, A. B., Titov, V. V., & Thomson, R. E. (2020). 
Empirical size–frequency relation for distant tsunamis… Journal of 
Geophysical 
Research: 
Oceans, 
125(5), 
e2019JC015647. 
https://doi.org/10.1029/2019JC015647 
[17] Frontiers in Earth Science. (2020). Probabilistic Tsunami Hazard Analysis: 
HPC for Massive Inundation Simulations. Frontiers in Earth Science, 8, 
591549. https://doi.org/10.3389/feart.2020.591549 
[18] Gabriel, A.-A., Ulrich, T., Marchandon, M., Biemiller, J., & Rekoske, J. 
(2023). 3D Dynamic Rupture Modeling of the 2023 Türkiye Mw 7.8–7.7 
Earthquake 
Doublet. 
The 
Seismic 
Record, 
3(4), 
342–358. 
https://doi.org/10.1785/0320230006 
[19] Geist, E. L., & Lynett, P. J. (2014). Source processes for PTHA. 
Oceanography, 27(2), 86–93. https://doi.org/10.5670/oceanog.2014.40 
[20] Geist, E. L., & Parsons, T. (2005). Probabilistic analysis of tsunami 
hazards. 
Natural 
Hazards, 
37(3), 
277–314. 
https://doi.org/10.1007/s11069-005-4646-z 
[21] Genz, A. (1984). Testing multidimensional integration routines. In Ford et 
al. (Eds.). North-Holland. 
[22] Gobet, E., Lemaire, V., & Goudjil, O. (2022). Improved mean estimation 
techniques for RQMC. Journal of Computational and Applied Mathematics, 
411, 114251. https://doi.org/10.1016/j.cam.2021.114251 
[23] Goda, K. (2019). Effects of earthquake source uncertainty on tsunami 
inundation forecasts: A stochastic approach. Journal of Geophysical Research: 
Solid Earth, 124(4), 3804–3826. 
 


---

# PAGE 81

77 
[24] Goda, K., & Song, J. (2016). Uncertainty modeling and sensitivity analysis 
of tsunami inundation due to stochastic slip distributions. Coastal 
Engineering, 117, 86–100. 
[25] Gonzalez, F. I., LeVeque, R. J., Adams, L. M., & Eble, M. C. (2009). 
Probabilistic tsunami hazard assessment for Seaside, Oregon. Journal of 
Geophysical 
Research: 
Oceans, 
114(C11), 
C11023. 
https://doi.org/10.1029/2008JC005132 
[26] Grezio, A., Lorito, S., Maramai, A., & Selva, J. (2021). Towards a unified 
total tsunami hazard model (TotPTHA). Frontiers in Earth Science, 9, 768350. 
https://doi.org/10.3389/feart.2021.768350 
[27] Heidarzadeh, M., Gusman, A. R., Mulia, I. E., Chua, C. T., & Suppasri, A. 
(2025). Tsunami hazards and risks from the Philippine Trench. Ocean 
Engineering, 329, 120985. https://doi.org/10.1016/j.oceaneng.2025.120985 
[28] Hou, Z., Marzouk, Y. M., & Ng, L. W. T. (2019). QMC methods for 
Bayesian inverse problems. MIT Technical Report 19-02. 
[29] Ibrahim, A. A., Zhao, C., & Elbanna, A. E. (2024). Fully Dynamic Rupture 
Modeling of Faults in Saturated Porous Media. SCEC Poster. 
[30] Kotani, M., Mulia, I. E., & Satake, K. (2020). PTHA for the Japan Sea coast 
based on seismic and tsunami records. Pure and Applied Geophysics, 177(11), 
5057–5079. https://doi.org/10.1007/s00024-020-02598-7 
[31] Kotani, T., Tozato, K., Takase, S., Moriguchi, S., Terada, K., Fukutani, Y., 
Otake, Y., Nojima, K., Sakuraba, M., & Choe, Y. (2020). PTHA with 
simulation-based response surfaces. Coastal Engineering, 160, 103719. 
https://doi.org/10.1016/j.coastaleng.2020.103719 
[32] Krause, T. (2011). Probabilistic tsunami hazard assessment for the US East 
Coast [Master’s thesis, University of Rhode Island]. 
[33] Kutschera, F., Gabriel, A.-A., Krenz, L., Heimann, S., Li, M., & 
Halldorsson, B. (2024). Linked and fully coupled 3D earthquake dynamic 
rupture 
and 
tsunami 
modeling… 
Solid 
Earth, 
15(2), 
251–280. 
https://doi.org/10.5194/se-15-251-2024 
[34] Lallemand, S., Popoff, M., Cadet, J.-P., Bader, A.-G., Pubellier, M., Rangin, 
C., & Deffontaines, B. (1998). Genetic relations between the Philippine Trench 
and Sangihe Trench. Journal of Geophysical Research: Solid Earth, 103(B1), 
933–950. https://doi.org/10.1029/97JB02620 
[35] L'Ecuyer, P., & Lemieux, C. (2002). A practical guide to randomized 
quasi-Monte Carlo. INFORMS Journal on Computing, 14(4), 370–381. 
 


---

# PAGE 82

78 
[36] LeVeque, R. J., Waagan, K., González, F. I., Rim, D., & Lin, G. (2016). 
Generating Random Earthquake Events for PTHA. Pure and Applied 
Geophysics, 173, 3671–3692. https://doi.org/10.1007/s00024-016-1357-1 
[37] Leonard, M. (2014). Self-consistent earthquake fault-scaling relations: 
Update and extension to stable continental strike-slip faults. Bulletin of the 
Seismological Society of America, 104(6), 2990–3002. 
[38] Li, M., Gabriel, A.-A., Krenz, L., Heimann, S., & Halldorsson, B. (2023). 
Dynamic Rupture Models… Iceland. Journal of Geophysical Research: Solid 
Earth, 128(7), e2023JB026194. https://doi.org/10.1029/2023JB026194 
[39] Lotto, G. C., Jeppson, T. N., & Dunham, E. M. (2018). Fully-coupled 
simulations of megathrust earthquakes and tsunamis… Manuscript excerpt. 
[40] Mai, P. M., & Beroza, G. C. (2002). A spatial random field model to 
characterize complexity in earthquake slip. Bulletin of the Seismological 
Society of America, 92(8), 2969–2977. 
[41] Mazet-Roux, G., Heinrich, P., & Hébert, H. (2019). Stochastic kinematic 
earthquake rupture modeling for PTHA. Pure and Applied Geophysics, 
176(2), 497–514. https://doi.org/10.1007/s00024-018-1993-3 
[42] Mulia, I. E., Ishibe, T., Satake, K., Gusman, A. R., & Murotani, S. (2020). 
Regional PTHA along eastern margin Sea of Japan. Earth, Planets and Space, 
72, 123. https://doi.org/10.1186/s40623-020-01256-5 
[43] Mulia, I. E., Satake, K., & Guzmán, G. F. (2020). PTHA in the Sea of Japan 
considering fault geometry & slip uncertainty. Pure and Applied Geophysics, 
177(11), 5035–5056. https://doi.org/10.1007/s00024-020-02559-0 
[44] Nakayama, A., & Tuffin, B. (2024). A Central Limit Theorem for RQMC 
Estimators. SIAM Journal on Scientific Computing, 46(2), A1343–A1366. 
https://doi.org/10.1137/23M1577789 
[45] Nakayama, A., & Tuffin, B. (2024). Sufficient conditions for CLTs for 
RQMC. 
Operations 
Research 
Letters, 
52(1), 
101–107. 
https://doi.org/10.1016/j.orl.2023.10.012 
[46] Nakayama, M. K., & Tuffin, B. (2024). Sufficient conditions for CLTs for 
RQMC methods. ACM Transactions on Modeling and Computer Simulation, 
34(3). https://doi.org/10.1145/3643847 
[47] Niederreiter, H. (1992). Random number generation and quasi-Monte 
Carlo methods. SIAM. 
[48] Okal, E. A. (2017). From 3D seismic structure to tsunami generation. Earth 
and 
Planetary 
Science 
Letters, 
479, 
388–396. 
https://doi.org/10.1016/j.epsl.2017.09.034 
 


---

# PAGE 83

79 
[49] Ökten, G. (1998). Error estimation for QMC methods. In Niederreiter et al. 
(Eds.). Springer. 
[50] Owen, A. B. (1997). Monte Carlo variance of scrambled equidistribution 
quadrature. SIAM Journal on Numerical Analysis, 34(5), 1884–1910. 
[51] Owen, A. B. (1997). Scrambled net variance… Annals of Statistics, 25(4), 
1541–1562. 
[52] Owen, A. B. (2025). Error estimation for quasi-Monte Carlo. 
arXiv:2501.00150. 
[53] Ozawa, A., Tagami, T., Listanco, E. L., Arpa, C. B., & Sudo, M. (2004). 
Initiation and propagation of subduction along the Philippine Trench. Journal 
of 
Asian 
Earth 
Sciences, 
23(1), 
105–111. 
https://doi.org/10.1016/S1367-9120(03)00112-3 
[54] Pasmann, S., & Tramm, J. (2025). RQMC Sampling in the Random Ray 
Method for Neutron Transport. M&C 2025, 1851–1858. ANS. 
[55] Ramos, N. T., Tsutsumi, H., Perez, J. S., & Bermas, P. P., Jr. (2012). Uplifted 
marine terraces in Davao Oriental… Journal of Asian Earth Sciences, 45, 
114–125. https://doi.org/10.1016/j.jseaes.2011.07.028 
[56] Ramalingam, N. R., Johnson, K., Pagani, M., & Martina, M. L. V. (2025). 
Machine learning surrogates for nearshore tsunami hazard. Natural Hazards 
and 
Earth 
System 
Sciences, 
25(5), 
1655–1679. 
https://doi.org/10.5194/nhess-25-1655-2025 
[57] Ramalingam, R., Baba, T., Hirata, K., & Nishikawa, T. (2025). 
Deep-learning 
surrogate 
modeling 
for 
tsunami 
hazard prediction… 
Geophysical 
Journal 
International, 
in 
press. 
https://doi.org/10.1093/gji/ggad295 
[58] Razafindrakoto, H. N. T., Cotton, F., & Mai, P. M. (2015). Quantifying 
variability in earthquake rupture models: Lessons from the 2011 Tohoku 
event. Journal of Geophysical Research: Solid Earth, 120(1), 384–405. 
[59] Sánchez-Linares, C., de la Asunción, M., Castro, M. J., González-Vida, J. 
M., Macías, J., & Mishra, S. (2016). Uncertainty quantification in tsunami 
modeling… 
Journal 
of 
Mathematics 
in 
Industry, 
6(5). 
https://doi.org/10.1186/s13362-016-0022-8 
[60] Sara Aniko Wirp, Gabriel, A.-A., Krenz, L., Lorito, S., Selva, J., Romano, 
F., Basili, R., Bader, M., & Dunham, E. M. (2020). Coupled dynamic 
earthquake rupture-tsunami modeling for the Hellenic Arc. AGU Abstract. 
[61] Schorlemmer, D., Marzocchi, W., Lavecchia, G., Wong, K., Boschi, E., & 
Zechar, J. D. (2020). The global tsunami model initiative and its approach to 
 


---

# PAGE 84

80 
probabilistic tsunami hazard and risk assessment. Pure and Applied 
Geophysics, 177, 1405–1425. 
[62] Sørensen, M. B., Spada, M., Babeyko, A., & Sørensen, U. (2012). 
Probabilistic tsunami hazard in the Mediterranean Sea. Journal of Geophysical 
Research: Oceans, 117(C8), C08032. https://doi.org/10.1029/2012JC007710 
[63] Tong, W., Shao, G., & Ji, C. (2023). Optimization-based PTHA using 
stochastic slip models. Geophysical Journal International, 232(1), 45–62. 
https://doi.org/10.1093/gji/ggad030 
[64] Ulrich, T., Madden, E. H., & Gabriel, A.-A. (2022). Stress, rigidity and 
sediment strength control megathrust earthquake and tsunami dynamics. 
Nature Geoscience, 15(1), 1–9. https://doi.org/10.1038/s41561-022-01047-9 
[65] Williamson, A., Dunbar, P. K., Chamberlin, C., & Kong, L. (2022). Efficient 
stochastic modeling of earthquake slip for tsunami hazard applications. 
Natural Hazards and Earth System Sciences, 22(8), 2883–2906. 
[66] Zimmerman, M. B., Powers, J. P., & Geist, E. L. (2025). Logic tree 
approaches 
in 
PTHA: 
An 
overview. 
Natural 
Hazards, 
in 
press. 
https://doi.org/10.1007/s11069-025-06153-9 
 
 
 


---

