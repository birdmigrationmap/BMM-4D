# Bird aspect correction (dev)

This study looks for the best way to remove the aspect bias in radar bird density without removing real gradients. The bias comes from birds seen side-on reflecting more than birds seen head-on or tail-on.

The data are the us_east night of 21–22 October 2026: 107 radars and 26 frames from 23:30 to 12:00 UTC, every 30 min. The starting point is the cleaned volumes with no aspect correction. These are recovered from `prod/volumes/us_east_fix` by undoing its symmetric fit; the recovery matches the no-shape cache `322b6cd8` to within 0.04%.

This is dev only: nothing under `toolbox/`, `examples/` or `studies/artifacts/` was changed.

```matlab
project_root = setup_paths();
run(fullfile(project_root, "studies", "aspect", "scripts", "AS_A_bins.m"))         % overlap cells, ~25 min (26 frames)
run(fullfile(project_root, "studies", "aspect", "scripts", "AS_B_overlap.m"))      % fits and scores, ~4 min
run(fullfile(project_root, "studies", "aspect", "scripts", "AS_C_polarimetric.m")) % raw Z_H / ZDR / rhoHV, 2 frames, ~7 min
run(fullfile(project_root, "studies", "aspect", "scripts", "AS_D_figures.m"))      % figures, ~1 min
```

Data are in `<local_data_dir>/aspect/`:
- `bins/` (9 MB per frame);
- `polar/`;
- `overlap_results.mat`, `polar_results.mat` and `aspect_models.mat`;
- `test_clean/`, which holds one volume from `as_clean_volume`.

Nothing was downloaded.

## Notation

ψ = azimuth from the radar − direction of the birds' motion.
- ψ = 0 is tail-on: the bird flies away from the radar.
- ψ = 180° is head-on.
- ψ = ±90° is side-on.

The aspect function is log T(ψ) = b2·(−cos 2ψ) + b1·cos ψ.
- b2 > 0 means side-on birds are brighter.
- b1 < 0 means tail-on birds are dimmer than head-on birds.

All amplitudes are in natural-log units. For small values they equal the linear amplitude a in 1 + a cos.

## What the literature says

- **Single birds in the lab, X band, wings closed.** Broadside is 14–90 × head-on (11–20 dB), and head-on is 1.1–1.4 × tail-on (Edwards & Houghton 1959, *Nature* 184:1059).
- **Standing ducks at UHF.** Breast toward the radar is 14 dB above tail toward (Blacksmith & Mack 1965, *Proc. IEEE* 53:1125).
- **Size and wingbeat.** Passerines are in the resonance region at S band, and the wingbeat modulates the cross-section by 10–20 dB (Vaughn 1985, *Proc. IEEE* 73:205).
- **Spread wings in a model.** In the electromagnetic model of a bat, spread wings reverse the HH pattern: side/head is 0.54 (Mirkovic et al. 2016, *Sci. Rep.* 6:35637). Flying birds are therefore much less anisotropic than lab bodies.
- **Birds in flight, X band.** Terns show broadside 3–4 dB above head and tail, i.e. 2–2.5× (Urmy & Warren 2017, *Methods Ecol. Evol.*).
- **No S-band value for passerine migration.** We found no published S-band Z_H side/head ratio for passerine migration.
- **Weather-radar mosaics.** Kranstauber et al. (2020, *Remote Sens.* 12:635) see higher density perpendicular to the flight direction, at most about a factor of 2.
- **Correction in standard tools.** vol2bird and bioRad (Dokter et al. 2011, 2019) do not correct for aspect. They rely on averaging over azimuth.
- **Expected amplitude.** Heading spread σ damps the population amplitude by exp(−2σ²), so a plausible population b2 at S band is 0.05–0.3. A front/back term is physically expected.
- **Polarimetric signatures.**
  - rhoHV is lowest tail-on and highest head-on. This is used to get the heading (Stepanian & Horton 2015, *IEEE TGRS*; Stepanian et al. 2016, *Ecosphere* 7:e01539; Horton et al. 2016, *Sci. Rep.*).
  - NEXRAD ZDR and rhoHV azimuth patterns are distorted by the transmit differential phase (Melnikov et al. 2015, *JTECH*; Melnikov & Bridge 2024, *IEEE GRSL*), so they are unreliable for the amplitude.
  - Insects have stronger and sharper side-on peaks, and their axis can differ from the birds'.
- **Not checked in full text:** Konrad 1968, Bruderer & Joss 1969, Zrnić & Ryzhkov 1998, and Diehl & Larkin 2005.

## What the data say

### Radar overlap: the gradient-free measurement (`AS_B_overlap.m`, `figures/AS_aspect_function.png`)

When two radars see the same 4 × 4 km × 100 m cell at the same time, the log ratio of their densities depends only on their two viewing angles, not on the density field. The method:
- Cells are moved to the frame time with each radar's velocity profile.
- Fixed effects per pair, frame and 500 m layer absorb calibration, per-volume scale and beam smoothing.
- A log range-ratio slope per layer absorbs beam broadening.
- Errors are clustered by radar pair.

The data: 2.96 M cell pairs from 197 radar pairs, at 0.8–3.5 km above the antennas.

**Measured aspect function** (24 bins, no shape imposed): tail-on 0.70, side-on 1.29, head-on 1.01, relative to the azimuthal mean.
- Side/head is 1.28.
- Side/tail is 1.84.
- Head/tail is 1.44.

**Harmonics**
- b2 = 0.200 ± 0.012 and b1 = −0.208 ± 0.019.
- The 4th harmonic is 0.020 ± 0.009.
- The terms that are zero by symmetry are zero: sin 2ψ is 0.014 ± 0.012 and sin ψ is −0.017 ± 0.020.
- Shape: cos 2ψ plus cos ψ; no 4th harmonic is needed.

**Robustness**

| Subset | b2 | b1 |
| --- | ---: | ---: |
| equal ranges, matched beams | 0.209 | −0.19 |
| 23:30–03:00 UTC | 0.198 | −0.19 to −0.22 |
| 03:30–07:00 UTC | 0.184 | −0.19 to −0.22 |
| 07:30–12:00 UTC | 0.195 | −0.19 to −0.22 |
| Z_H | 0.225 | −0.32 |
| Z_V | 0.227 | −0.15 |

- b1 is −0.19 to −0.22 in all three periods (the range is shared across the time rows).
- Z_V is modulated as much as Z_H, so Z_V is no aspect-free proxy.
- b1 depends on polarization, which points to the scatterers rather than radar processing.

**Varies with height, speed and region**

| Height above antennas (km) | 0.75–1 | 1–1.5 | 1.5–2 | 2–2.5 | 2.5–3.5 |
| --- | ---: | ---: | ---: | ---: | ---: |
| b2 | 0.07 | 0.13 | 0.16 | 0.19 | 0.24 |
| b1 | −0.39 (22 pairs) | −0.17 to −0.25 | −0.22 | −0.24 | −0.16 to −0.21 |

| Speed (m/s) | < 6 | 6–10 | 10–14 | > 14 (82 pairs) |
| --- | ---: | ---: | ---: | ---: |
| b2 | 0.06 | 0.17 | 0.24 | 0.12 |

| Longitude | < −92° | −92° to −84° | > −84° |
| --- | ---: | ---: | ---: |
| b2 | 0.28 | 0.13 | 0.12 |

**Axis** (2 frames, `AS_C_polarimetric.m`)

| Flight direction used | b2 |
| --- | ---: |
| track at the cell height | 0.203 |
| rhoHV heading | 0.201 |
| fixed 190° | 0.195 |
| volume track over 200–1500 m (as `ar_clean_volume`) | 0.094 |

- The low-level track is up to 23° off the flight direction aloft and halves the amplitude.
- With track and fixed-axis terms fitted together (all frames): b2 is 0.11 + 0.11 and b1 is −0.09 − 0.16.
- So the axis lies between the track and a preferred heading, consistent with partial drift compensation.

### Single radar (population means; real gradients average out only if they have no preferred direction)

**Every volume shows the same amplitude.** At 10–60 km the per-volume sd of the side-on harmonic equals that of the 45° harmonic (0.145 and 0.145). Per volume, the amplitude cannot be told apart from real structure. A per-volume fit (`ar_clean_volume`) therefore mostly fits structure.

**Smaller than the overlap at all heights:**
- b2 ≈ 0.00 below 750 m and 0.05–0.07 above.
- b1 is −0.07 to −0.14.
- b2 increases with range: −0.03 at 5–10 km, 0.03 at 10–30 km, 0.09 at 30–60 km and 0.11 at 60–150 km.

**ZDR does not vary with aspect.**
- The side-on harmonic of ZDR (= log Z_H/Z_V, gradient-free) is −0.011 ± 0.007.
- Its correlation with the Z_H harmonic is 0.02.
- So ZDR cannot predict the amplitude. rhoHV gives the axis about as well as the track (heading − track: median 3°, IQR 45°).

**Track turning in time (`as_time_fit`).** At each radar the track turns about ±20° during the night. A radar × sweep × band fixed effect removes the radar's persistent pattern. The results:
- b2 is −0.04 to +0.02, i.e. the side-on pattern does not turn with the track.
- b1 is −0.03 to −0.10.

**Per-radar amplitudes** (on the overlap, around the shared model):
- sd 0.24 against a standard error of 0.035.
- Reproducible between even and odd frames: r = 0.98.

With an almost uniform track in this night, these take every fixed azimuthal pattern of a radar onto the side-on axis: partial blockage, sector calibration and lasting gradients. They are not aspect. Per-radar fitting lowers the overlap disagreement further (0.589) for that reason, so it is not recommended.

## Methods compared

The overlap scores are out of fold: shared parameters are fitted on even frames and scored on odd frames, and the other way round.

| Column | Meaning |
| --- | --- |
| Disagreement | weighted RMS log ratio after the fixed effects and beam slopes; the reduction is relative to no correction |
| b2 / b1 left | aspect harmonics still on the overlap, ± se |
| Single-radar columns | 10–60 km, local track; first = median first harmonic (gradient across the disc plus front/back) |
| Wedge | deepest 10° sector against ±30° sectors on polar data: median, share below 0.8, and KDMX's deepest at 04:30–07:30 |

On polar data the wedge score is noise-dominated (it is about 0.80 for every method). Only the ellipse stands out.

| Method | Disagreement (reduction) | b2 left | b1 left | First | b2 | b1 | Wedge median / < 0.8 | KDMX deepest |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| no correction | 0.632 | 0.200 ± 0.016 | −0.220 ± 0.033 | 0.273 | 0.037 | −0.110 | 0.800 / 50% | 0.87 |
| ellipse (paper) | 0.613 (2.9%) | 0.073 ± 0.016 | −0.116 ± 0.040 | 0.202 | −0.025 | −0.037 | 0.774 / 60% | 0.72 |
| per-volume symmetric (`ar_clean_volume`) | 0.617 (2.3%) | 0.133 ± 0.019 | −0.217 ± 0.035 | 0.273 | −0.030 | −0.114 | 0.802 / 49% | 0.87 |
| shared symmetric, b2 = 0.20 | 0.599 (5.2%) | 0.009 ± 0.016 | −0.213 ± 0.034 | 0.279 | −0.159 | −0.121 | 0.803 / 48% | 0.85 |
| shared b2 + b1 (0.20, −0.21) | 0.597 (5.5%) | 0.011 ± 0.016 | −0.016 ± 0.033 | 0.237 | −0.157 | +0.097 | 0.802 / 49% | 0.81 |
| shared b2(h) + b1 | 0.597 (5.6%) | 0.016 ± 0.014 | −0.016 ± 0.034 | 0.240 | −0.028 | +0.093 | 0.804 / 49% | 0.83 |
| shared b2(h, speed) + b1 | 0.596 (5.7%) | 0.008 ± 0.014 | −0.016 ± 0.033 | 0.237 | −0.095 | +0.091 | 0.803 / 49% | 0.81 |
| **shared b2(h) + b1 ramp (recommended)** | **0.597 (5.6%)** | **0.016 ± 0.014** | **−0.018 ± 0.033** | **0.225** | **−0.026** | **+0.010** | **0.805 / 48%** | **0.83** |

### Reading the table

**The ellipse removes both aspect and gradients.** It removes part of the aspect: b2 left is 0.07 and b1 left is −0.12. It also removes gradients: the first harmonic falls from 0.27 to 0.20 at 10–60 km. At KDMX 06:30 its factor runs from 0.66 to 2.9 by azimuth, which is the bow-tie (`figures/AS_kdmx_maps.png`).

**The per-volume symmetric fit removes a third of the side-on effect aloft and none of the front/back term.**
- Its amplitude comes from the near range, where the effect is small: 0.02 at KDMX 06:30.
- Its track comes from 200–1500 m.

**Any shared function does better on the overlap than per-volume fits** (5.2–5.7% against 2.3–2.9%), because it does not fit the field.

**A constant amplitude over-corrects the near range** (single-radar b2 −0.16). The near range needs b2 to grow with height.

**Keeping b1 at its full value down to the ground** turns the near-range first harmonic around (+0.09). The ramp below 1 km brings it to +0.01. This leaves the near-range single-radar means close to zero without losing anything on the overlap.

**At KDMX 06:30 the recommended factor is 0.82–1.49 (1st–99th percentile of the gates within 90 km, 0.5–1.3° sweeps), smooth with azimuth and height.** It is largest downstream, where birds are seen tail-on, and aloft. The dense NW band and the low S–SE are kept; there is no bow-tie. KDMX's deepest polar sector goes from 0.87 to 0.83, against 0.72 with the ellipse.

## Recommendation

Use one shared aspect function, not fitted to each volume (`functions/as_clean_volume.m` → `as_apply_aspect` → `as_aspect_template`):

```
log T = b2(h) (−cos 2ψ) + b1(h) cos ψ − log N,   ψ = az − track(h)
b2(h) = 0.092 h_km,   b1(h) = −0.21 min(1, h / 1 km)
```

- **Height and track.** h is the beam-centre height above the antenna. track(h) is the direction of the volume's own EP-VVP profile at that height.
- **Normalization.** N gives T an azimuthal mean of 1, so a uniform field keeps its mean and no per-sweep rescaling is needed.
- **Applying it.** Each gate and its sd are divided by T before thinning.
- **Parameters.**
  - From this night's overlap: b2_per_km = 0.092 ± 0.005 and b1 = −0.21 ± 0.02 (fold values 0.092 / 0.092 and −0.209 / −0.213).
  - They are passed in `filter.aspect`: a struct with b2, b2_per_km, b2_per_ms, h0_m, v0_ms, b1, b1_ramp_m.
  - The defaults are written in `as_clean_volume`. `filter.aspect = []` leaves the volume uncorrected.
  - The parameters used are saved in `volume.aspect`, and `as_read_volume` undoes them.
- **Swap-in.** `as_clean_volume` was tested on KEAX 06:30: it equals the no-shape cache divided by T to machine precision. Use it with `read_options` `"shape", false`.

### Estimating the parameters for another night

1. Clean once without any shape or aspect correction (`"shape", false`, and `filter.aspect = []` or `bmm4d.radar.clean_volume`).
2. Run `as_frame_bins` with a manifest whose `volume_files` point to that cache and `fix_dir` = the same cache, frame by frame (about 1 min each).
3. Run `as_overlap_pairs`, then `as_overlap_design(pairs, "raw")`, then `as_fit_aspect_model(D, true(size(D.y)), radar_height, "height")`, and set `b1_ramp_m = 1000`.
4. Apply the parameters to the cached volumes with `as_apply_aspect`, which is cheap (no re-cleaning). Rerun `as_clean_volume` only if a fresh cleaning is wanted anyway.

```matlab
volume = load(file).volume;
volume = as_apply_aspect(volume, aspect);
save(new_file, "volume")
```

## What remains uncertain

- **Below about 1 km**, where most birds fly, there is no gradient-free measurement. The overlap starts at 0.8 km and has few pairs below 1.25 km.
  - The model there is an extrapolation: b2 ≈ 0.05 at 0.5 km, and b1 ramps to zero at the ground.
  - Single-radar means agree with it, but they include real gradients.
  - The 1 km ramp height is a choice, not a fit.
- **Axis: heading or track.**
  - The side-on pattern does not follow the track when the track turns during the night.
  - On the overlap, the axis lies between the track and a fixed preferred heading.
  - Without wind data, the heading is only available from rhoHV, which needs the raw volume. It did as well as the track at the cell height in 2 frames, but was not tested on all frames.
- **One autumn night, tracks mostly toward S.**
  - The amplitude varies by region (0.12–0.28) and with ground speed.
  - The b2(h, speed) model is only marginally better on the overlap and over-corrects the near range.
  - No spring night was tested: every other cache has the shape correction applied, and re-cleaning was limited.
  - The parameters should be refit for each night or season. A night with a different track direction would also separate aspect from fixed geographic patterns.
- **Above 3 km**, b2 linear in height overshoots: 0.30 against 0.24 measured.
- **No inversion was run.** Wedges and gradients were scored on polar volumes, not on reconstructed columns as in `studies/artifacts` (`AR_B_shape.m`). A solve of `ar_midwest` with `as_clean_volume` (or `as_apply_aspect` on the no-shape cache `322b6cd8`) is the next check.
  - Done since in `studies/artifacts` (`AR_I_aspect.m`, full night `us_east_shared`).
- **Aspect or fixed radar patterns.** We fitted the shared function together with a fixed azimuthal pattern per radar in the geographic frame (`as_overlap_fit_radar`; 1st and 2nd harmonics, summing to zero over the radars).
  - On this night, b2_per_km falls from 0.092 to 0.017 and b1 from −0.21 to −0.06.
  - With one track direction, and each radar seeing its neighbours in a few directions only, the two cannot be told apart. Only nights with other tracks can decide.

## Multi-night plan

`PLAN_multi_night.md` is a plan, with code, for 16 more nights (2021–2025, chosen from Dark Ecology for different track directions). It asks:
- whether the shape holds;
- whether the function is aspect or radar patterns;
- whether one fixed function is enough, or a fit per night is needed;
- which covariates matter;
- what happens below 0.8 km.

It needs the user's OK for about 90 GB of downloads and about 20 h of machine time. The scripts are below; the parts tested so far are listed in the plan.

## Files

- **Scripts:**
  - `scripts/AS_A_bins.m`: overlap cells and sector sums.
  - `scripts/AS_B_overlap.m`: fits, strata, scores.
  - `scripts/AS_C_polarimetric.m`: Z_H / Z_V / ZDR / rhoHV.
  - `scripts/AS_D_figures.m`: figures.
  - `scripts/AS_E_select_nights.m`: the nights of the multi-night plan (Dark Ecology archive, no download).
  - `scripts/AS_F_multi_night_run.m`: download, clean and bin them in a batch job (pilot, then main run).
  - `scripts/AS_G_multi_night_fit.m`: the multi-night analysis and decision.
- **Functions:**
  - `as_night_bins` and `as_multi_night_job`: binning and the batch job of the multi-night plan.
  - `as_overlap_fit_radar`: the shared function with fixed radar patterns (sparse).
  - `as_clean_volume`, `as_apply_aspect` and `as_aspect_template`: the correction.
  - `as_read_volume`: reads a volume uncorrected.
  - `as_frame_bins` and `as_frame_polar`: binning.
  - `as_overlap_pairs`, `as_overlap_design`, `as_overlap_fit`, `as_fit_strata`, `as_aspect_curve` and `as_overlap_score`: the overlap regression.
  - `as_fit_aspect_model`, `as_fit_radar_offsets` and `as_pair_templates`: models.
  - `as_volume_harmonics`, `as_cell_harmonics`, `as_sector_corrected` and `as_time_fit`: single-radar analyses.
  - `as_polar_summary`: rhoHV heading and ZDR.
  - `as_volume_map`: KDMX maps.
- **Figures:**
  - `figures/AS_aspect_function.png`: the measured aspect function and its amplitudes against height.
  - `figures/AS_kdmx_maps.png`: KDMX 06:30 UTC before and after, with the factor of each correction.
  - `figures/AS_overlap_validation.png`: overlap disagreement, aspect left, the function before and after, and near-range harmonics.
- **Tables:** `figures/AS_*.csv`.
