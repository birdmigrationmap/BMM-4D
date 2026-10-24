# Model review: what the process looks like and how to model it

Prepared 23 October 2026 for the prod (dev) reconstruction. Sources: a review of the radar-aeroecology literature, a code review of BMM-US (`birdmigrationmap/BMM-US`), the references and model description of our own manuscript, and diagnostics on the 21–22 October 2026 eastern-US night (101 radars, 26 frames). Nothing here is decided; section 6 lists the choices.

## 1. What the process looks like

### Vertical structure
- **Most migration is low and above the ground, not above sea level.**
  - Eastern-US flight heights average about 420–560 m above ground (AGL) at six NEXRADs (Horton et al. 2016, *Biol Lett*).
  - Across 143 NEXRADs, the 50% activity envelope starts at about 355 m AGL and is about 520 m wide, higher and wider in the west (Giuntini et al. 2025, *Ecology*).
  - In the northeast, densities concentrate at roughly 500–2000 m above sea level (La Sorte et al. 2015).
  - Our manuscript states that roughly half of nocturnal migration is below 700 m AGL.
- **Shape changes slowly and smoothly, with occasional layering.**
  - Birds climb during the first ~30% of the night, then descend (Horton 2016).
  - Altitude alone explains 50–70% of density variation; birds go higher when winds aloft are better (Kemp et al. 2013).
  - Two layers can occur, e.g. about 0.8 and 1.9 km (Dokter et al. 2011).
- **Height reference is unresolved in the literature.** Kranstauber et al. 2020 overestimated densities near the Alps, probably because birds rise over terrain. Dokter et al. 2018 bin heights above sea level. For the eastern US only the Appalachians and the western plains edge make this matter.

### Horizontal and temporal structure
- **Large-scale correlation.**
  - European column densities correlate 50% at 300 km and 25% at 500 km, with a 16% nugget; this is not identified below about 100 km (Nussbaumer et al. 2019).
  - BMM-US fitted about 500–600 km spatial range in migration seasons (October: about 600 km) and about 100–250 km in winter or summer.
- **Within a night.**
  - Correlation is 67% at 1 h (Nussbaumer 2019). BMM-US intra-night fits give about 1–2.4 h and 100–450 km.
  - Take-off is synchronous: 90% of tagged songbirds leave within 69 min of civil dusk (Cooper et al. 2023). Density is highest in the first half of the night.
- **Skewed and heteroscedastic.**
  - In Europe the median column density is 8 and the mean 19 birds km⁻² (Nussbaumer 2019).
  - A strong power transform (p ≈ 0.13) was needed for normality, and uncertainty grows with intensity.
  - Intra-night variance is larger at dusk and dawn (BMM-US `nnt_std`).

### Radar artefacts that shape the reconstruction
- **Beam overshoot.**
  - The 0.5° beam is about 1.5 km above the radar at 100 km.
  - Range corrections (Buler & Diehl 2009; Buler et al. 2012; Kranstauber 2020) all assume one vertical profile per radar domain. Kranstauber's adjusted column is literally *column total × fixed shape*.
  - Buler 2012 excludes gates that sample 5% or less of the profile to avoid false zeros.
- **Detection limit.**
  - The WSR-88D minimum detectable signal is about −19.5 dBZe at 50 km in clear-air mode (−9.5 in precipitation mode).
  - It degrades by 6 dB per doubling of range. A "no echo" gate means *below the limit*, not *zero*.
- **Near-range suppression.** vol2bird uses 5–35 km and excludes gates within 5 km. Our data show a much larger effect, see section 3.
- **Calibration.**
  - The specification is ±1 dB, but clear-air comparisons between neighbouring WSR-88Ds can differ by several dB.
  - Nussbaumer 2019 cross-validation found constant radar-specific biases with no spatial pattern.
- **Other contaminants.**
  - Wind turbines (blade Doppler, not removed by clutter filters).
  - Anomalous propagation under nocturnal inversions (|v| < 1 m/s masks).
  - Insects (low airspeed, high ZDR, depolarization ratio; mostly the lowest layer).
- **Body orientation.** It changes the cross-section by up to a factor of 2, and the cross-section itself drifts within seasons (11 ± 6 cm²). These are multiplicative errors.

### Closest methodological analogues
- **Vertical profile of reflectivity (VPR) as an inverse problem.** In hydrology, Andrieu & Creutin 1995 and Germann & Joss 2002 estimate local profiles from volume scans: column × shape with a beam operator.
- **3D variational assimilation** with a background covariance built from recursive filters and mass continuity as a weak constraint (Gao et al. 2004). This is almost our objective, with the covariance made explicit.
- **Gaussian-process tomography** (Li, Svensson et al. 2013) and **Gaussian-process velocity retrieval for bird radar** (Angell & Sheldon 2018), both with a linear line-integral operator and a GP prior.
- **What doesn't exist yet:** no published INLA/SPDE, hurdle or Tweedie model for radar bird density. Nussbaumer 2019 (power transform plus GP) is the closest.

## 2. Lessons from BMM-US

| BMM-US choice | Keep, adapt or drop |
| --- | --- |
| Work on the column total (VID), 2D | Keep for the horizontal part: the column total is where the horizontal structure is. Add the vertical shape as a separate, smooth field. |
| Below 10 birds km⁻² is censored ("nothing matters"), imputed from a truncated Gaussian (~75% of night samples) | Adapt: don't impute; *down-weight* by density. The structure is learned and fitted only where birds are; empty areas just have to stay near zero. |
| Box-Cox on the upper tail, kriging the transformed values, back-transform without bias correction | Adapt: a multiplicative (log-scale) model, but solved in density space by iterative reweighting, so the beam operator stays linear and there is no back-transform bias. |
| Trends removed first (day-of-year mean/std, DEM, LOWESS mean/std maps, night-time template), stationary residual | Adapt: replace with hierarchical levels estimated jointly: a large-scale field (the "trend") and a small-scale field whose variance scales with it. That's the heteroscedastic structure BMM-US got from dividing by `sf_std` and `nnt_std`, without a separate detrending step. |
| Separable space × time covariance: exponential in space, Cauchy in time; nuggets forced to zero | Keep separability and the ranges as starting values. Allow a nugget: our data show a radar-level nugget, see section 3. |
| Night-time (NNT) mean template and variance | Keep the idea: take-off and landing are tied to sun angle. We already relax advection there; the variance (how freely the field may change) should also depend on night time. |
| Covariance fitted by weighted least squares on empirical covariances per 14-day window | Keep the empirical approach for starting values, fitted only where birds are. Later, possibly maximum likelihood. |
| Daily + intra-night components | Here the window is one night at 30 min: the "daily" level becomes the large-scale (inter-radar) field, the intra-night one the fine (intra-radar) field. |

The code review also found:
- debug restrictions in the committed scripts (2010 only, one night, no simulations);
- a dropped elevation term on the grid (`f_heightdem`);
- swapped legend labels in the parameter plots;
- a parameter reuse bug in the speed fit.

These matter if BMM-US is rerun, not for the design here.

## 3. Diagnostics from our own night (us_east, 21–22 October 2026)

From the 101 radar vertical profiles (100 m layers) and their column totals (VID), frames 01:30–10:30 UTC:

| Quantity | Value | Consequence |
| --- | --- | --- |
| VID distribution | quantiles 10/25/50/75/90/99%: 5 / 15 / 36 / 74 / 127 / 273 birds km⁻²; 10% below 5, 31% below 20 | Skewed but not dominated by zeros on a migration night. A density weight or threshold around 10–20 matches BMM-US. |
| Vertical shape, above ground, 0–2.5 km | 1, 2, 3 EOFs explain 76%, 91%, 96% of the shape variance (1808 profiles) | 3 shape patterns are enough. Mean shape: median height 750 m AGL, 90% below 1850 m. |
| VID semivariogram (fraction of variance) | 75 km: 0.15 (log: 0.19); 300 km: 0.64; 425 km: 0.86; sill about 500–600 km (log: about 800 km) | Large-scale range about 500–800 km, consistent with BMM-US October (about 600 km). A nugget of about 15–19% between neighbouring radars: station effects. |
| Shape semivariogram (EOF1 / EOF2) | 75 km: 0.28 / 0.23; 425 km: 0.41 / 0.70; 850 km: 0.74 / 0.93 | The shape varies at the regional scale, slower than VID, with a radar-level nugget of about 25%. |
| Heteroscedasticity, neighbours (< 250 km) | median ΔVID by level 5–20 / 20–50 / 50–100 / 100–200 / > 200: 4.9 / 14.8 / 31.5 / 50 / 116, about 0.4 × the level everywhere | A multiplicative process (constant coefficient of variation about 0.4), not additive. |
| Heteroscedasticity, 30-min changes | 5–12% of the level | Same. |
| 30-min autocorrelation of log VID | 0.99 | Strong temporal persistence: advection (or a temporal GP) carries most of the information between frames. |

Near-range suppression, measured earlier (lowest sweep, median over 59 radars with birds, relative to 30–50 km):

| Range (km) | 0–5 | 5–10 | 10–15 | 15–20 |
| --- | ---: | ---: | ---: | ---: |
| Median ratio | 0.06 | 0.27 | 0.89 | 1.19 |
| 25% quantile | — | — | 0.36 | 0.64 |

In the sandbox, every regularization choice left the radar imprint unless these readings were down-weighted (`studies/sandbox/README.md`).

## 4. What our current model does, in this light

Objective: data misfit + VP prior (radar profiles spread horizontally) + 3D Laplacian smoothness + advection, all quadratic in linear density, with uniform data weights.

**Strengths:**
- The beam operator is exact.
- Advection carries information through time, which matches the 0.99 autocorrelation.
- The no-echo gates are used.

**Mismatches:**
- **Additive regularization on a multiplicative process.** The same absolute deviation costs the same in empty and dense areas, so dense areas are over-smoothed relative to their natural variability, and empty areas get structure they don't need.
- **No explicit large/small-scale split.** The prior is a fixed, hand-built large-scale field (radar profiles plus inverse-distance or Gaussian weights). Its deviations are uncorrelated (white), so data corrections stop at the beam edge. That's the radar imprint.
- **Vertical shape is not separated from amount.** Height dependence is carried by 30 free layers plus a prior per sea-level height. 3 shape patterns above ground would carry 96% of it.
- **Every observation is treated as an exact, unbiased linear measurement.**
  - The near-range suppression, the detection limit at long range (no echo ≠ zero), station calibration offsets and multiplicative cross-section errors are all ignored.
  - We have been compensating with hard cuts or weights, which you'd rather avoid.

## 5–7. Two-level separable model (set aside, 24 October 2026)
The proposed model (density = column amount × a few vertical shape patterns, coarse level V_L plus fine deviation, near-range recovery curve instead of data cuts), the decisions taken with Raphaël, the prototype plan and its results are now in `studies/two_level/README.md`. It removed the imprint in the sandbox but was worse than production on the full us_east night, so work went back to the production (lean) model.

## 8. Learning from many nights (later, a separate task)

One night cannot validate the model or set its parameters: nights differ strongly, and the peak nights that matter most happen only a few times per season. Plan for later:
- **A test library of nights,** stratified by intensity (peak, typical, weak), season (spring, autumn), weather (clear, partly rainy, frontal), and region (plains, east, Gulf coast, terrain), for example 20–40 nights. Regional cases at sandbox size (about 30 radars) cost minutes per night after cleaning; full eastern-US nights cost about 1.5 h of cleaning each.
- **Peak nights selected from historical column-density time series** (BMM-US VID data 1995–2021, or profiles from cleaned volumes), not by hand.
- **Hyperparameters learned across nights:** correlation lengths and variances per season from radar profiles alone (cheap: no inversion needed), as BMM-US did for 1995–2021 by day of year. Its fitted tables are a starting point.
- **Validation:** leave-one-radar-out on every test night, scored on predicted vs withheld observations, density-weighted, and pooled by night type, so that a choice is judged on the whole library and not on one night.

## References
- Andrieu H, Creutin JD (1995) *J Appl Meteorol* 34:240–258.
- Angell R, Sheldon D (2018) Inferring latent velocities from weather radar data using Gaussian processes. NeurIPS.
- Buler JJ, Diehl RH (2009) *IEEE Trans Geosci Remote Sens* 47:2741–2751.
- Buler JJ et al. (2012) *PLoS ONE* 7:e41571.
- Cooper NW et al. (2023) Songbirds initiate migratory flights synchronously relative to civil dusk. *Mov Ecol*.
- Dokter AM et al. (2011) *J R Soc Interface* 8:30–43.
- Dokter AM et al. (2018) *Nat Ecol Evol* 2:1603–1609.
- Dokter AM et al. (2019) bioRad. *Ecography* 42:852–860.
- Gao J, Xue M, Brewster K, Droegemeier KK (2004) *J Atmos Ocean Technol* 21:457–466.
- Germann U, Joss J (2002) *J Appl Meteorol* 41:542–557.
- Giuntini S et al. (2025) Structuring the skies. *Ecology* 106:e70247.
- Horton KG et al. (2016) Where in the air? *Biol Lett* 12:20160591.
- Kemp MU et al. (2013) *Ibis* 155:734–749.
- Kranstauber B et al. (2020) *Remote Sens* 12:635.
- La Sorte FA et al. (2015) *R Soc Open Sci* 2:150347.
- Li D, Svensson J et al. (2013) *Rev Sci Instrum* 84:083506.
- Lin T-Y et al. (2019) MistNet. *Methods Ecol Evol*.
- Nussbaumer R et al. (2019) *Remote Sens* 11:2233.
- Nussbaumer R et al. (2021) *J R Soc Interface* 18:20210194.
- BMM-US code: github.com/birdmigrationmap/BMM-US (density/ and speed/ inference live scripts).
