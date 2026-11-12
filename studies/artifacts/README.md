# Artefacts of the us_east night (dev)

Five places in the us_east video (21–22 October 2026, round 0) where the reconstruction looks wrong, what makes each, and fixes tested one at a time on small cases cut from the same night, before a full us_east run. Dev only, except the beam-blockage lookup, now fixed in the toolbox (`8cabaa1` on `v2-prep`, MIGRATION.md "Beam-blockage maps read at compass azimuth"); `dev/` holds the study override the test cleaning jobs used before (now the same as the toolbox).

```matlab
project_root = setup_paths();
run(fullfile(project_root, "studies", "artifacts", "scripts", "AR_A_cases.m"))   % test cases (observations ~1–3 min each)
run(fullfile(project_root, "studies", "artifacts", "scripts", "AR_B_shape.m"))   % fix 1
run(fullfile(project_root, "studies", "artifacts", "scripts", "AR_C_terrain.m")) % fix 2
run(fullfile(project_root, "studies", "artifacts", "scripts", "AR_D_edges.m"))   % fix 3
runtests(fullfile(project_root, "studies", "artifacts", "scripts", "test_ar_terrain.m"))
```

**Test cases** (`ar_test_case`): radars in a box, 04:30–07:30 UTC every 30 min (7 frames), grid = radar extent + 160 km, us_east cleaning, solve grid and domain; the us_east manifest and static clutter masks are reused (no download, no clutter stage). `ar_midwest` (28 radars, KDMX), `ar_appalachians` (16), `ar_plains` (12, to the −104° domain edge). Solves use the production inversion settings with plain PCG (`ar_solve`).

## What each screenshot shows

| Screenshot | Where | Cause |
| --- | --- | --- |
| 1 | High Plains: dark discs with circular edges at KLNX, KGLD, KUDX | The prior does not honour each radar's own level: its 150 km Gaussian blend mixes in denser neighbours (at KLNX the prior aloft is 20–40 birds km⁻³ where KLNX reads 0–5). The white prior lets the field fall back to it within ~1 cell where the rows stop (here ~134 km: farther beams leave the 3000 m grid top). KUDX also has the beam-blockage bug. |
| 2 | Blue Ridge / Smokies (KMRX, KGSP): column over the ridges 0.38 of the valleys, streaks downwind | (a) The beam-blockage maps are indexed by math angle (counter-clockwise from east) but read with compass azimuths: blocked gates are kept and their "no echo" becomes zero density (KMRX 0.5° at 70–100 km toward the Smokies reads 0.2 birds km⁻³; valley rows at the same height above ground 45–100). (b) Prior and advection are layered by height above sea level: the ridges cut off the dense low layer (prior alone: ridge / valley 0.40). |
| 3 | Missouri (KEAX, KSGF): arcs between radars | As 1: level steps where the data of one radar end. |
| 4 | Central Texas: bright NE–SW line through KEWX–KGRK at 02:30, moving west ~13 m/s | In the data: a fine line (convergence boundary), 20–25 dBZ, ρHV 0.4–0.6 (biological, not rain), high ZDR. Kept, as all biological echo; no fix. |
| 5 | KDMX 05:30–07:30: wedge NE of the radar between dense areas | The anisotropy (shape) correction. Without it KDMX is densest NW and lowest S–SE (a real gradient; NE about the median); the elliptical template (different front and back axes, oriented by ρHV) absorbs that gradient and dividing by it makes a bow-tie (NE 0.65–0.8 of the median, SE 1.3–1.4). |

The beam-blockage bug: `studies/beam_blockage/beam_blockage_script.m` builds the polar grid with `pol2cart(deg2rad(AZ), R)` (x = R cos az: angle from east, counter-clockwise), and `bmm4d.radar.filter.beam_blockage` looks the map up at the gate's compass azimuth. Map column a holds compass bearing 90° − a (correlation of the map with the DEM horizon: 0.93–1.00 that way, −0.68 to 0.23 as used, for KMRX, KCXX, KFCX, KRLX, KUDX). At 90 km and 0.5°, 55–113 of 360 azimuths of these radars keep blocked gates. `dev/+bmm4d/+radar/+filter/beam_blockage.m` looks up 90° − az (no new maps needed).

## Fixes and results

### Fix 1: aspect correction (screenshot 5) — `AR_B_shape.m`

Harmonics of each radar's observations at 10–60 km per frame (medians): the first is a gradient across the disc; the second is taken along the side-on axis (track + 90°, track from the radar's profile at 200–1500 m), i.e. the aspect effect.

| Cleaning (ar_midwest) | Wedge, median | Radar-frames with a wedge < 0.8 | KDMX deepest | 1st harmonic | Side-on 2nd harmonic (share > 0) | Misfit (30–120 km) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| elliptical shape correction (production) | 0.896 | 3.1% | 0.85 | 0.15 | −0.01 (43%) | 0.210 (0.195) |
| no correction | 0.896 | 0.5% | 0.89 | 0.27 | 0.12 (86%) | 0.232 (0.215) |
| symmetric aspect correction (+ blockage fix) | 0.904 | 1.5% | 0.91 | 0.28 | 0.02 (73%) | 0.230 (0.212) |

- The aspect effect is real but smaller than the gradients: without a correction the side-on harmonic is positive in 86% of radar-frames, median 0.12 (side-on ≈ 1.27 × head-on), while the gradients across the discs are 0.27. The elliptical template removes both, the gradients included (0.27 → 0.15), and this is how it makes the bow-ties.
- The symmetric correction, T(az) = 1 + a cos 2(az − track − 90°), a from the volume's own second harmonic on that axis, within [0, 0.5] (`ar_clean_volume`), removes the side-on pattern and keeps the gradients; it cannot make one-sided wedges. On KDMX 06:30 a = 0.02. The maps around KDMX (`figures/AR_B_shape.png`) show one dense band NW of the radar, joined to KOAX's, instead of the bow-tie. The three wedges left below 0.8 are KFSD 04:30 and KUEX 07:00–07:30.
- The misfits are against each case's own observations: the production ones are azimuthally flatter and easier to fit, so the misfit does not rank the cleanings. The symmetric correction fits its data slightly better than no correction (0.230 against 0.232).
- **Chosen:** the symmetric aspect correction; superseded by the shared aspect function below.

**Follow-up: the shared aspect function of `studies/aspect`** (`AR_I_aspect.m`). That study measured the aspect effect where two radars see the same cell at the same time (no gradient can enter): side-on 1.29, head-on 1.01, tail-on 0.70 of the mean, growing with height, and proposed one function shared by all volumes instead of a fit per volume. The same cleaned volumes (beam-blockage fix) with no correction, the per-volume symmetric fit, and the shared function, solved with the production settings; held-out radars as in fix 2b (prior rebuilt without the radar).

| Cleaning (ar_midwest) | Wedge median / < 0.8 | KDMX deepest | First | Side-on | Front/back | Misfit | Held-out: relative error | median \|log ratio\| |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ellipse (round 0) | 0.896 / 3.1% | 0.85 | 0.15 | −0.01 | −0.04 | 0.210 | 0.487 | 0.430 |
| no correction | 0.896 / 0.5% | 0.89 | 0.27 | 0.12 | −0.18 | 0.232 | 0.544 | 0.477 |
| per-volume symmetric | 0.904 / 1.5% | 0.91 | 0.28 | 0.02 | −0.19 | 0.230 | 0.546 | 0.479 |
| shared function | 0.899 / 0.5% | 0.88 | 0.19 | 0.06 | −0.05 | 0.221 | 0.510 | 0.414 |

Front/back: the first harmonic on the track axis (> 0 denser where the birds fly away from the radar). Held out: KDMX, KDVN, KOAX, KEAX, KARX.
- The shared function predicts the held-out radars better than no correction (0.510 against 0.544; median log ratio 0.414, the best of the four), the per-volume fit does not (0.546): it fits the field, not the aspect (as `studies/aspect` found on the polar data).
- The ellipse has the lowest relative error (0.487), but it flattens each radar's own data (first harmonic 0.15), which makes a held-out radar easier to predict without being more right; the median log ratio, less driven by the dense areas, puts the shared function first.
- The shared function removes most of the front/back pattern (−0.18 → −0.05) and half of the side-on one at 10–60 km (its amplitude grows with height, near the radar it is small), and keeps the dense band NW of KDMX without a bow-tie (`figures/AR_I_aspect.png`).
- **Chosen:** the shared function. The full-night run used the per-volume fit (it ran before this test).

### Fix 2: terrain (screenshot 2) — `AR_C_terrain.m`, `test_ar_terrain.m`

- `ar_steps(..., "agl")`: an advection link keeps the height above the lowest active cell of its column (birds carried over a ridge rise with it), `ar_profile_prior(..., "agl")`: the radar profiles are laid from each column's ground up. Tests: with "asl" both equal the toolbox exactly; AGL links carry a field that depends only on height above ground without residual; the AGL prior's ridge / valley column is 0.98 (production 0.40).

| Solve (ar_appalachians) | Ridge / valley column | Column ratio per km of ground | Misfit (30–120 km) |
| --- | ---: | ---: | ---: |
| production (ASL) | 0.41 | 0.35 | 0.223 (0.201) |
| AGL links | 0.45 | 0.43 | 0.220 (0.197) |
| AGL prior | 0.59 | 0.50 | 0.220 (0.198) |
| AGL links + prior | 0.63 | 0.61 | 0.223 (0.201) |
| blockage fixed (ASL) | 0.46 | 0.39 | 0.217 (0.192) |
| blockage fixed + AGL links + prior | 0.76 | 0.72 | 0.217 (0.193) |

- With the production cleaning, the few rows over the ridges (smoothed ground > 1000 m) read 2–8 birds km⁻³ at 0–1200 m above ground, against 45–100 in the valleys: KMRX and KGSP 0.5° and 0.9° at 80–100 km, beams blocked by the Smokies that the buggy lookup keeps (KMRX 0.5° toward 165° is 100% "no echo" at 70–100 km and becomes 0.2 birds km⁻³). With the lookup fixed these gates are missing, and the misfit drops.
- The two parts need each other: the blockage fix alone gives 0.46, AGL alone 0.63, both 0.76 (`figures/AR_C_terrain.png`: the dark Smokies band is mostly gone).
- What is left is in the data: the rows that now see the ridges (800–2000 m above ground) read 0.55–0.75 of the valley rows at the same height above ground, but more than the valley rows at the same height above sea level: the birds follow the terrain only in part. Columns over high ground also end at the 3000 m grid top (1500–1800 m of air over the ridges).
- **Chosen:** the blockage lookup fixed, and AGL links and prior.

### Fix 2b: advection over terrain — `AR_G_terrain_advection.m`

Does the advection carry birds over the ridges? Three rules for the advection links, each with its prior, on `ar_appalachians_bb` (blockage fixed):
- ASL (production): a link keeps its height above sea level. A link into a ridge is dropped, and the lee gets no inflow from below the ridge top.
- AGL: a link keeps its height above the ground. The column rises over a ridge, but the layers pushed above the 3000 m grid top are dropped.
- sigma (`ar_steps` / `ar_profile_prior` "sigma"): a link keeps its relative height between the ground and the top of its column (source cell j of n layers links to every destination layer overlapping [j, j+1)/n of the destination column) and carries the ratio of the two columns' air thickness, x_dest = (H_src / H_dst) x_src, so a column squeezed over a ridge keeps its birds. The production solver only links equal densities, so this needs a dev copy with link ratios (`ar_solve_advection_lean`; without ratios it is the toolbox solver). The sigma prior lays each radar's profile as birds per unit relative height, interpolated between radars and divided by each column's air thickness. Tests (`test_ar_terrain.m`): sigma links carry a field with the same birds in every column exactly; sigma prior ridge / valley column 1.01 (ASL 0.40, AGL 0.98).

| Rule | Uphill source cells that keep a link | Lee cells that receive one |
| --- | ---: | ---: |
| ASL (production) | 0.85 | 0.89 |
| AGL | 0.83 | 0.79 |
| sigma | 1.00 | 0.98 |

Uphill: destination ground > 300 m above the source's; lee: in a column reached from ground > 300 m higher.

| Solve (white prior) | Ridge / valley column | Per km of ground | Misfit (30–120 km) | Held-out radars: relative error | median \|log ratio\| |
| --- | ---: | ---: | ---: | ---: | ---: |
| ASL (production) | 0.46 | 0.39 | 0.217 (0.192) | 0.508 | 0.568 |
| AGL | 0.76 | 0.72 | 0.217 (0.193) | 0.484 | 0.510 |
| sigma | 0.85 | 0.81 | 0.219 (0.193) | 0.482 | 0.508 |
| sigma links, AGL prior | 0.78 | 0.76 | 0.218 (0.193) | 0.480 | 0.499 |

Held out: KMRX, KGSP, KFCX, KRLX, KJKL, KHTX in turn (rows weight 0, prior rebuilt from the other radars' profiles), their rows at 30–120 km predicted from the other radars (`ar_holdout_score`, `AR_G_holdout.csv`). A first version kept the left-out radar's profile in the prior (the prior falls back on it where the data stop), which favoured the solves that stay closest to the prior: ASL then came out best (0.442 against 0.457–0.462).

- The production advection does not handle the terrain: 15% of the uphill links go into the ground and 12% of the lee cells get no inflow. AGL loses about as much through the grid top. Sigma loses nothing.
- Following the terrain predicts the held-out mountain radars better: relative error 0.48 against 0.51, median |log ratio| 0.50–0.51 against 0.57; most at KFCX (on a ridge: 0.63 → 0.51) and KHTX (0.54 → 0.48), equal at KJKL and KRLX, slightly worse at KMRX (0.46 → 0.50). AGL and sigma are equal on this test.
- The section across the Smokies at 06:00 (`figures/AR_G_section.png`): ASL stops the dense low layer at the ridge; AGL lifts it over; sigma squeezes it into a dense band over the ridges (× H_valley / H_ridge, about 1.5). That squeeze comes from the 3000 m grid top acting as a lid: real birds can climb above it. A higher grid top over mountains (or a hybrid coordinate that follows the terrain near the ground and is flat higher up, under a raised top) would remove the lid; not tested.
- AGL and sigma predict equally; AGL (fix 2) stays in the full-night run. Sigma is the one that conserves birds and would be the choice with a higher grid top.

### Fix 2c: how birds cross terrain — `AR_K_terrain_hybrid.m`

What terrain the grid uses, what the radar data say about flight over the ridges, and the simplest advection that fits them. Case `ar_appalachians_bb`, production settings, white prior; solves and held-out radars in a batch job (`ar_terrain_solves`). New modes of `ar_steps` / `ar_profile_prior`: `"smooth"`, `"hybrid"`, option `detour` (tests in `test_ar_terrain.m`).

**The grid and the terrain** (`figures/AR_K_grid.png`)

- `data/reference/DEM/dem.mat` is GTOPO30: 30″ (0.93 × 0.75 km at 36° N), int16, not smoothed (the Clingmans Dome pixel is 1981 m, the summit 2025 m).
- `create_grid` samples it at the 2 km grid cell centres (nearest). The smoothing line, `imgaussfilt(dem, [100 100])` (a Gaussian of 100 pixels, about 75–90 km), is commented out. `coarsen_grid` averages 2 × 2 cells to the 4 km solve grid (`g.dem`). Over the Smokies the 4 km ground is within 7 m (median) of the true 4 km mean; peaks stand up to 714 m above their cell (median 72 m).
- Heights are above sea level (`g.z` 50–2950 m, top 3000 m). `build_active_mask` switches off the cells whose top is below `g.dem`; `air_fraction` gives the lowest cell its air share. The prior and the advection links are layered by height above sea level.
- So there is no smoothed topography and no terrain-following coordinate. The only terrain in the model is the cells switched off under the unsmoothed 4 km ground; the smoothed topography is the commented-out line.
- The paper code (`BMM-within`) has the same `create_grid`, but nothing reads `g.dem`: no cell is switched off and the grid starts at sea level everywhere. Its `construct_F` builds a blockage height map (`bb_alt`) and does not use it either. The Midwest case (ground 150–450 m) hardly needed either.
- Rockies edge (39.75° N, same sampling): at the us_east edge (−104°) the 4 km ground is 1540 m, so the column has 1460 m of air under the 3000 m top; west of −105.7° the ground is above the top. A domain further west needs a higher top.

**What the radar data show** (`AR_K_data.png`, `AR_K_coordinate.png`; `AR_K_profiles.csv`, `AR_K_same_beam.csv`, `AR_K_coordinate.csv`)

Each observation row gets its beam centre (F-weighted height over the active cells) and the ground under its footprint. Rows at 30–150 km, 7 frames, 5.6 M observations, 0.37% of them over ground above 1000 m. Density is relative to the radar-frame mean.

| Relative density | ground < 300 m | 300–600 m | 600–1000 m | > 1000 m |
| --- | ---: | ---: | ---: | ---: |
| 1100 m above ground | 1.99 | 1.89 | 1.19 | 0.45 |
| 1500 m above ground | 1.22 | 0.86 | 0.44 | 0.31 |
| 1300 m above sea level | 1.94 | 2.46 | 2.54 | 2.14 |
| 1500 m above sea level | 1.54 | 1.89 | 1.94 | 1.90 |
| 1700 m above sea level | 1.18 | 1.32 | 1.50 | 1.45 |
| 1900 m above sea level | 0.84 | 0.90 | 1.07 | 0.90 |

- By height above ground, the ground classes differ by 2–4 times. By height above sea level they agree within about 25%, with the high ground 1.1–1.3 times the lowlands. Aloft, the density is set by the height above sea level, whatever the ground below.
- **Same beam, different ground.** Rows of one radar, frame, sweep and 5 km range bin share the beam height above sea level (within ~70 m) and differ only in azimuth. If birds followed the terrain, rows over higher ground would be denser, because they sit closer to the dense low layer: × 1.48 for ground 400 m higher at 600–1200 m above ground, × 3.8 for 1200 m higher at 1200–1800 m (AGL expectation from the lowland profile). They read 0.91 and 0.84. Above 1800 m they are a little denser (1.08–1.51, against 2.4–7.2 expected): at most a fifth of the AGL lift, in log terms. The lowest beams (< 600 m above ground) read less over higher ground (0.64 at +200 m; few rows).
- **Which height coordinate collapses each radar-frame's rows onto one profile?** The score is the median |log residual| from the radar-frame's profile, on rows over ground ≥ 300 m from the radar's.
  - AGL scores 0.61 and ASL 0.40.
  - The hybrid (terrain following near the ground, fading to level at H) equals ASL for H ≤ 1000 m and is worse above.
  - Height above the ground smoothed over L is worse than ASL below 150 km and equal at 150–300 km. On all rows it is the best (0.343 at 150–200 km, against 0.348 for ASL and 0.357 for AGL).
- Within one radar's view (≤ 150 km) the birds fly level over the ridges, at every height the radars see there. What follows the terrain is the regional level: a radar on high ground has its own profile shifted up with its ground. This is the KFCX case (Blue Ridge plateau), where AGL predicted better in fix 2b.

**What the radars cannot see** (`AR_K_trust.png`, `AR_K_blockage.csv`, `AR_K_section.png` top)

- Over the ridges (ground above 1000 m, 679 columns), 62% of the columns have a row at all. Their lowest beam centre is 1045 m above the ground (median; 10th percentile 540 m), and the rows over them are at 0.6–1.7 km above the ground. The section at 35.6° N has no row below 1.3 km asl over the Smokies at 06:00.
- Share of azimuths blocked at 60–120 km (0.5° / 0.9°): KMRX 83% / 40%, KRLX 35% / 0%, KGSP 26% / 22%, KFCX 25% / 3%, KJKL 9% / 0%, KHTX 0%. The kept gates have the lower half-power edge above the terrain horizon (`8cabaa1`), so they lose at most ~5% of the two-way beam to the terrain. Ground clutter would make ridge rows denser; they are lighter or equal.
- **Limit:** the lowest ~500 m over the ridges is not in the data. Whether the low layer climbs, goes around or is not there comes from the model (prior and links). This is also one night, with flight to the SW, along the ridges.

**Rules tested**

- **smooth L** (links and prior): keep the height above the ground smoothed with a Gaussian of L. A fraction-of-a-layer shift is split between two layers (coefficient √(share of the cell's air covered)) rather than rounded away. Infinite L gives the ASL steps and prior exactly (tested). One constant.
- **hybrid H**: the lift fades linearly from the ground rise, at the ground, to 0 at H above the ground. Links are level above H, and everywhere the rise is ≥ H (the low layer cannot follow a ridge higher than H). The prior uses the same mapping from each radar's column.
- **detour**: a link whose destination layer is under the ground of an active column goes sideways instead of being dropped. It goes to the nearest column, perpendicular to the displacement, where that layer is air, at most the displacement length away (≤ 45°). Equidistant columns on both sides share it.

| Links (`AR_K_links.csv`) | Uphill sources that keep a link | Lee cells reached | Links (M) | Detour links |
| --- | ---: | ---: | ---: | ---: |
| ASL (production) | 0.85 | 0.89 | 12.4 | 0 |
| AGL | 0.83 | 0.79 | 12.4 | 0 |
| sigma | 1.00 | 0.98 | 15.9 | 0 |
| hybrid H 1000 m | 1.00 | 0.98 | 16.8 | 0 |
| smooth L 100 km | 0.85 | 0.92 | 24.4 | 0 |
| smooth L 100 km + detour | 0.99 | 0.93 | 24.4 | 77 k |
| ASL + detour | 0.99 | 0.89 | 12.5 | 81 k |

**Solves and held-out radars** (`AR_K_solves.csv`, `AR_K_holdout.csv`, `AR_K_section.png`, `AR_K_maps.png`). Each mountain radar is left out in turn, as in 2b, with the prior rebuilt without it. ASL and AGL reproduce 2b exactly. The 25 and 50 km rows ran in a second batch job with the same settings.

| Rule (links + prior) | Ridge / valley column | Misfit (30–120 km) | Held out: relative error | median \|log ratio\| |
| --- | ---: | ---: | ---: | ---: |
| ASL (production) | 0.46 | 0.192 | 0.508 | 0.568 |
| AGL | 0.76 | 0.193 | 0.484 | 0.510 |
| sigma (2b) | 0.85 | 0.193 | 0.482 | 0.508 |
| hybrid H 1000 m | 0.50 | 0.192 | 0.504 | 0.565 |
| smooth L 25 km | 0.64 | 0.192 | **0.482** | **0.502** |
| smooth L 50 km | 0.58 | 0.193 | 0.488 | 0.513 |
| smooth L 100 km | 0.51 | 0.193 | 0.497 | 0.535 |
| smooth L 100 km + detour | 0.51 | 0.193 | 0.497 | 0.535 |
| smooth L 150 km + detour | 0.48 | 0.192 | 0.501 | 0.553 |
| ASL + detour | 0.47 | 0.192 | 0.508 | 0.568 |

| Held out: relative error | KMRX | KGSP | KFCX | KRLX | KJKL | KHTX |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ASL (production) | 0.465 | 0.564 | 0.634 | 0.459 | 0.383 | 0.543 |
| AGL | 0.508 | 0.569 | 0.511 | 0.453 | 0.380 | 0.485 |
| hybrid H 1000 m | 0.469 | 0.564 | 0.626 | 0.459 | 0.377 | 0.528 |
| smooth L 25 km | 0.488 | 0.560 | 0.515 | 0.461 | 0.378 | 0.492 |
| smooth L 50 km | 0.471 | 0.557 | 0.536 | 0.465 | 0.390 | 0.511 |
| smooth L 100 km | 0.455 | 0.557 | 0.578 | 0.475 | 0.389 | 0.526 |

- **Detour changes nothing measurable** (ASL 0.5077 with and without). The 77–81 k links it adds are the lowest layers beside the ridges, where no radar sees (above). It is ~60 lines and harmless, but these data cannot test it.
- **The hybrid ("only the low layer follows") is ASL on this test** (0.504). Over the ridges the low layer is never seen, and between radars it does not move the profiles.
- **The held-out gain of AGL is regional, not ridge-scale.** It comes from KFCX (on the Blue Ridge plateau: 0.63 → 0.51) and KHTX. Yet KFCX's own rows are level within its view: same-beam ratios are 0.99–1.06 over ground 200–700 m lower or higher. What matters is that each radar's profile sits on its own regional ground. Smoothing over L keeps that and flies level over anything narrower than L:
  - 25 km matches AGL and sigma on the held-out radars (0.482, median log 0.502, the best);
  - 50 km is close (0.488);
  - 100–150 km gives back half the gain.
  - Smoothing also lifts the column over the ridges less than AGL (ridge / valley 0.58–0.64 against 0.76) and does not stack a lifted layer above the crests (`AR_K_section.png`).
- **The two tests pull L apart.** Rows within one radar's view want L ≥ 150 km (level flight, coordinate score above); profiles carried between radars want L ≤ 50 km. One night, six radars: the difference between 25 and 50 km (0.482 against 0.488) is within the radar-to-radar spread.

**Chosen, without further fitting: height above the ground smoothed over 50 km, with climb** (in the toolbox since `7e6eddf`: `inversion.ground_smoothing_m`, MIGRATION.md "Prior and advection by height above the smoothed ground")

- **Prior:** each radar's profile is laid by height above the ground smoothed with a Gaussian of 50 km (`ar_profile_prior(..., "smooth", smooth_length_m=50e3)`).
- **Links:** keep the height above that smoothed ground (`ar_steps(..., "smooth", smooth_length_m=50e3, climb=true)`). A shift of a fraction of a layer is split between the two nearest layers. A link whose destination is under the real 4 km ground goes to the lowest air cell of that column (the birds climb) instead of being dropped; the sideways detour changed nothing and is left out.
- **Why this rule:**
  - Birds follow plateaus and basins, and fly level over ridges and valleys narrower than ~50 km. It is the same rule on small ridges and on the Rockies, whatever the ridge orientation.
  - The radars keep measuring at their true altitude: the grid and the forward model are unchanged.
  - Mass: the links tie equal densities, and the smoothed ground varies slowly. A link therefore shifts a column without squeezing it, so it carries the same birds.
  - The solver is unchanged, and an infinite length is exactly today's ASL.
- **Losses:** source cells left without a link fall from 2.6% to 2.1% (mostly the open domain boundary). Those lost through the 3000 m top are 0.2%, so the top need not rise for us_east; it must rise further west anyway.

**Check for artefacts** (section "Check of the chosen rule"; `AR_K_check_<case>_sections.png`, `_ratio.png`, `_maps.png`, `.csv`)

| Case | Rule | Misfit (30–120 km) | Lowest cell / 300 m above, ground > 1000 m | ground < 450 m | Edge (fix 3) |
| --- | --- | ---: | ---: | ---: | ---: |
| Appalachians | ASL (production) | 0.192 | 1.36 | 1.07 | 1.49 |
| | AGL | 0.193 | 1.06 | 1.05 | 1.26 |
| | smooth 50 km + climb | 0.193 | 1.17 | 1.05 | 1.43 |
| Plains | ASL (production) | 0.149 | 1.28 | 1.08 | 2.86 |
| | AGL | 0.149 | 2.05 | 1.30 | 2.71 |
| | smooth 50 km + climb | 0.149 | 1.65 | 1.26 | 2.76 |

- **Appalachians:**
  - The column ratio to ASL is smooth everywhere: up to ~1.3 over the uplands, slightly below 1 on the coastal plain. There are no edges or stripes.
  - In the sections the low layer is gently raised over the uplands. There is no layer stacked above the crests and no spike over single tall cells (AGL has one at 35.0° N, −84.5°).
  - Climbing does not pile birds into the lowest cell over the ridges (1.17, below ASL's 1.36).
- **Plains (ground up to 1900 m at −104°):**
  - Over ground above 1000 m the column moves by about as much as with AGL (median 13%, 95th percentile ×1.6; AGL ×1.6), and the low layer sits on the High Plains.
  - The rings in the ratio map are the end-of-data discs of fix 3: the edge score is unchanged (2.76 against 2.86), and the GP prior removes them.
- Same misfit and iterations as ASL in both cases. Plain PCG, white prior, one night.

### Fix 3: coverage edges (screenshots 1, 3) — `AR_D_edges.m`

Edge score (`ar_edge_score`): on every seen column next to an unseen one, the mean log column within 12 km on each side; the step where the data end against the same step 40 km inside the data (1 = the end of the data does not show).

| Solve (ar_plains) | Step at the edge | Step inside | Edge | Misfit (30–120 km) | Time | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| production (white prior) | ×1.22 | ×1.07 | 2.82 | 0.165 (0.149) | 22 s | 39 |
| prior 75 km | ×1.20 | ×1.07 | 2.58 | 0.163 (0.147) | 19 s | 37 |
| smoothed prior 75 km (second solve) | ×1.15 | ×1.08 | 1.99 | 0.164 (0.148) | 49 s | 39 + 27 |
| GP 50 km | ×1.12 | ×1.08 | 1.41 | 0.160 (0.144) | 143 s | 214 |
| GP 100 km | ×1.11 | ×1.08 | 1.26 | 0.160 (0.144) | 193 s | 283 |
| GP 150 km | ×1.10 | ×1.08 | 1.21 | 0.161 (0.145) | 158 s | 315 |

Times of the final run, with a cleaning job beside it (the first runs alone: 22 / 92 / 69 / 124 s for production / GP 50 / 100 / 150).

- Why the discs: in the west–east section through KLNX (41.96°N), KLNX's rows end ~134 km out (farther beams leave the 3000 m grid top) and read 0–5 birds km⁻³ west of the radar, while the prior is 20–40: the 150 km Gaussian blend of radar profiles pulls in KUEX's dense profile even at KLNX. A 75 km blend honours each radar more but changes the edges little.
- GP: the existing correlated prior option of `solve_advection_lean` (`prior_correlation_m`, same marginal variance as the white prior): a deviation set by the data fades over that length, so the field is continuous where the data end. The discs around KLNX, KGLD, KUDX and KUEX's circular edge go (`figures/AR_D_edges.png`); it also fits the data slightly better. Cost: 5–8 × the iterations (a warm start saves 10%). Seen with it: GP 150 km spreads KLNX's low far rows into KUEX's dense area (a dark lens at 1.3–2 km), and at 06:00–07:30 a thin arc of KUEX's far rows shows west of KUEX (inside the disc before).
- With the blockage fix (`ar_plains_bb`) the edge ratios are the same (2.86 white, 1.23 GP 100 km).
- In the Appalachians the radars are close and the edges are weaker (1.21); AGL 1.13, AGL + GP 100 km 0.99.
- Smoothing the first solve's deviation into a new prior (`ar_smooth_prior`, column support, one or two passes, 75 or 150 km) stops at an edge ratio of about 2: the step between the local data level and the smoothed one remains.
- **Chosen:** GP 100 km (edges nearly as smooth as 150 km, less spreading), with the smoothness term kept (below).

**Smoothness term and held-out radars** (`AR_H_gp_smoothness.m`, `ar_plains_bb`). Roughness: mean |second difference| along z and x over the mean density. Held out: KLNX, KUEX, KGLD, KOAX, KTWX (prior rebuilt without the radar).

| Solve | Edge | Misfit | Roughness z | Roughness x | Iterations | Held-out: relative error | median \|log ratio\| |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| white, smoothness 0.3 (production) | 2.86 | 0.165 | 0.079 | 0.091 | 38 | 0.802 | 0.862 |
| GP 100 km, smoothness 0.3 | 1.23 | 0.160 | 0.124 | 0.066 | 281 | 0.807 | 0.838 |
| GP 100 km, smoothness 0 | 1.18 | 0.159 | 0.448 | 0.071 | 320 | 0.805 | 0.837 |
| GP 100 km / 300 m vertical, smoothness 0 | 1.23 | 0.158 | 0.246 | 0.080 | 284 | 0.829 | 0.851 |

- Horizontally the GP replaces the smoothness term: without it the edges, the misfit and the horizontal roughness are the same. Vertically it does not: the GP correlates only horizontally, and the profiles become 3.6 × rougher; a 300 m vertical correlation halves that and predicts worse. The smoothness term stays.
- The GP does not predict the held-out radars better on average: the same relative error as the white prior, a slightly lower log error, better at KLNX (the disc: 1.22 → 0.98) and worse at KTWX (0.85 → 1.00) and KUEX (0.51 → 0.60). It removes the visible edges; it does not make the field between radars more right. (A first version kept the left-out radar's profile in the prior and showed the GP worse, 0.71 against 0.60: the white prior falls back on that profile.)

### Fix 4: the Texas line (screenshot 4)

No fix. The line is in the reflectivity (a fine line along a moving boundary) and its echo is biological (low ρHV), so it is kept like all biological echo. Its density is real echo read with the bird cross-section; it moves with the boundary, not with the bird winds, so the advection term does not carry it.

## Full night

`AR_E_us_east.m` (cleaning in a batch job, cache `prod/volumes/us_east_fix`, then the observations and one solve with AGL links and prior and GP 100 km, `ar_invert_fix` as a batch job; section 4 repeats it with the shared aspect function, `prod/us_east_shared`), `AR_F_video.m` (video with the round-0 colour scale and the five spots, round 0 against the fixes: `figures/us_east_night_<fix|shared>.mp4`, `AR_F_spots_<fix|shared>.png`; set `fix_folder` and `output_name`).

| Step | Time | Note |
| --- | ---: | --- |
| Cleaning (2392 volumes, 5 workers) | 2.6 h | 1869 cleaned, 523 from the test cases, 1 failed; each MistNet worker grows to 3.5–6 GB, 8 workers thrash on this 39 GB machine |
| Observations | 7.8 min | 5.83 M rows, as round 0 |
| Solve (batch, 8 threads) | 100 min | 352 PCG iterations; the worker held up to ~36 GB (compressed and swapped; resident 5–14 GB) |
| Video (batch) | 10 min | |
| Shared aspect function: volumes from the `us_east_fix` cache | 23 min | no new MistNet run: the aspect correction is undone and the shared function applied (`as_apply_aspect`) |
| Shared: observations, solve, video | 7 min, 64 min, 7 min | 354 PCG iterations |

- Misfit 0.361 (30–120 km 0.338) on its own observations; round 0 0.362 (0.333) on its own; the shared aspect function 0.361 (0.339).
- The five spots (`AR_F_spots_fix.png`): the High Plains discs (KLNX, KGLD) and the Missouri arcs are gone; the Smokies band is gone; KDMX has no wedge (one dense band NW of the radar); the Texas line is still there, as it should.
- The shared aspect function changes little at this scale (`AR_F_spots_shared.png`; log columns at 04:00 UTC correlate 0.96 with the per-volume fit).
- Patches around many radars, stretched NNE–SSW, clearest in the Appalachians (KJKL, KFCX, KHTX, KGSP). I first read them as more birds near the radars than between them; `AR_J_near_radar.m` says otherwise: averaged around each radar the column is flat from 40 to 150 km and *lower* within 20 km, and it is not denser along the flight direction (below). The patches look like the field between the radar-centred dips rather than an excess.

### Near the radars — `AR_J_near_radar.m`

Column against the distance to the nearest radar, each radar and frame relative to its own 60–100 km ring, median over radar-frames (`figures/AR_J_distance.csv`, `AR_J_near_radar.png`):

| Distance (km) | 0–10 | 10–20 | 20–30 | 30–40 | 40–60 | 100–150 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| us_east, per-volume fit | 0.71 | 0.83 | 0.93 | 0.98 | 1.00 | 1.01–1.05 |
| us_east, shared function | 0.72 | 0.83 | 0.93 | 0.98 | 1.00 | 1.01–1.04 |
| Midwest (white prior), reconstruction | 0.86 | 0.91 | 1.01 | 1.05 | 1.03–1.06 | 0.94–0.98 |
| Midwest, prior | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00–1.03 |

- Every run has a dip within 20 km of the radars (us_east 0.72 at 0–10 km); the prior has none, so it comes from the data. Along the track (within 30° of it) against across it at 30–100 km the column is 0.97 (per-volume fit) and 0.96 (shared): no plume along the flight direction on average.
- The data near the radars disagree with the model by elevation. Observed over predicted density of the Midwest solve (the forward model holds the beam, so at a fixed range an elevation dependence is in the data):

| Sweep at 5–15 km | 0.5° | 0.9° | 1.3° | 1.8° | 2.4° | 3.1° | 4.0° | 5.1° | 6.4° |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Observed / predicted | 1.04 | 1.07 | 1.06 | 0.92 | 0.93 | 0.93 | 0.91 | 0.88 | 0.82 |

- The step is between 1.3° and 1.8°, the boundary between the NEXRAD split cuts (separate surveillance and Doppler scans) and the batch cuts. It is ~0.87 at every predicted density (0.82 against 0.94 in the sparsest quarter, 0.94 against 1.07 in the densest), strongest at 5–10 km (0.76–0.88 above 1.6° against 1.02–1.04 below) and gone beyond 20 km (`AR_J_residuals.csv`). The same cells seen by two radars (`studies/aspect` overlap pairs) are 0.4–0.6 lower from the radar at 3–7° than from the one at 0.5°, but most of that is the far radar's 2 km wide beam, which the forward model accounts for.
- It is the batch cuts within 20 km. Solved without them, the column near the radars is flat (0.997 at 0–10 km, 1.02 at 10–20 km), and the dropped rows read 0.72 of what the other data predict. Solved without the split cuts within 20 km instead, the dip stays (0.86, 0.89) and those rows read 1.23 of the prediction (`AR_J_drop.csv`). The step is at 25 of the 28 radars (batch over split cuts at 5–15 km: median 0.86, 0.78–0.98).
- The raw reflectivity at 5–15 km (12 Midwest radars, one volume an hour 04–07 UTC, `AR_J_raw_5_15km.csv`):

| Cut | 0.5° | 0.9° | 1.3° | 1.8° | 2.4° | 3.1° | 4.0° | 5.1° | 6.4° |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Radials | 720 | 720 | 720 | 360 | 360 | 360 | 360 | 360 | 360 |
| Gates with echo | 0.96 | 0.99 | 1.00 | 0.93 | 0.91 | 0.93 | 0.94 | 0.95 | 0.97 |
| 1st percentile, dBZ at 1 km | −32.5 | −36.7 | −39.1 | −38.4 | −38.7 | −39.5 | −39.7 | −40.2 | −39.7 |
| Estimated density per unit Z, / 1.3° | 0.94 | 0.99 | 1 | 1.20 | 1.27 | 1.27 | 1.24 | 1.18 | 1.20 |
| Gates the clutter filter took power from | 0.96 | 0.89 | 0.83 | 0.84 | 0.83 | 0.78 | 0.73 | 0.69 | 0.66 |
| Their mean Z / that of the other gates | 0.32 | 0.22 | 0.15 | 0.13 | 0.12 | 0.13 | 0.13 | 0.14 | 0.14 |

- Not the detection limit: almost every gate has echo and the batch cuts are as sensitive as the split cuts. Not the gate-density estimator: it gives the batch cuts 20–27% *more* density per unit reflectivity than the split cuts, so their reflectivity itself is lower near the radar (by roughly 0.72 / 1.2 ≈ 0.6).
- A candidate: the clutter filter. Within 15 km it takes power from 66–96% of the gates of every cut up to 6.4° (from the CFP moment; on KDMX at 30–60 km the median removed above 1.2° is 0 dB), where the main lobe is 0.5–1.7 km above the ground. The gates it filters have 0.12–0.32 of the reflectivity of the others at the same range and cut: either it removes bird echo, or the clutter decision leaves the strongest bird echo unfiltered (selection); this test cannot tell the two apart. The cleaning reads the surveillance scan of each split pair (720 radials, 460 km range); in the batch cuts reflectivity comes from fewer pulses, and a filter notch might take a larger share of the bird spectrum there. Not settled.
- What it means for the maps: the dark centre of every radar (0.72 at 0–10 km in us_east). Options, none applied (method changes): drop or down-weight the batch cuts within ~20 km; keep only gates the filter did not touch (CFP code 2 or 8, read in the cleaning); or check the filtered fraction against the residual per gate first. Next test: the residual of near-range rows by CFP, which needs the CFP code per gate carried through the cleaning.
