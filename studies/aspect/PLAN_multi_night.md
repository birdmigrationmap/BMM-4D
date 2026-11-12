# Multi-night aspect fit: plan for a later run

> **6 Nov 2026: the data part moves to the cluster.** The nights now come from the shared night set (`studies/nights/PLAN.md`): 17 full nights plus the reference night, downloaded and cleaned once on merula, with the same cleaning settings. That replaces `AS_E` and `AS_F` below. `AS_G_multi_night_fit.m` still applies; point it at the night-set products and bin them with `as_night_bins`.

This plan is for an agent, or the user, running the analysis later. It needs the user's go-ahead for the download and about 20 h of machine time.

## Why one night is not enough

The 21–22 Oct 2026 night measures the aspect function precisely: b2_per_km = 0.092 ± 0.005 and b1 = −0.21 ± 0.02. But that night cannot say four things.

1. **Whether it holds on other nights.** It has one season and one flight direction (mostly toward S). The amplitude already varies by region (0.12–0.28) and with speed.
2. **Whether it is aspect at all, or fixed patterns of the radars.** Partial blockage, calibration by sector and lasting gradients stay in the geographic frame of each radar.
   - The joint fit with a fixed azimuthal pattern per radar (`as_overlap_fit_radar`, patterns summing to zero over the radars) brings b2_per_km down from 0.092 to 0.017 and b1 from −0.21 to −0.06 on this night.
   - With one flight direction and each radar seeing its neighbours only in a few directions, the two cannot be told apart. This is not evidence that the aspect is small; it shows that one night cannot decide.
3. **The shape and size below 0.8 km** above the antennas, where most birds fly. No two radars see the same cell there.
4. **Whether a fixed function is enough,** or whether it must be fitted again for every night.

**Why more nights help.** The aspect turns with the birds' track, while a radar's own pattern does not. So nights with different tracks separate them, both in the overlap and in single-radar harmonics.

The side-on term repeats every 180°. Spring and autumn (opposite tracks) therefore do not help it; tracks turned by 30–90° do. The front/back term flips between spring and autumn, so the seasons help that one.

## Questions and decision rules

Each question is answered by one section of `AS_G_multi_night_fit.m`.

**1. Shape: are cos 2ψ and cos ψ enough?**
- Analysis: section 2, the harmonics s2, s1 and c4 and the 24-bin curves, pooled and by season.
- Rule: keep the form when |s2|, |s1| and |c4| are below 0.03 and the curves are alike.

**2. Aspect or fixed radar patterns?**
- Analysis: section 3, the fit with and without radar patterns.
- Rule:
  - b2_per_km and b1 stay within 2 se → they are aspect.
  - They fall → the one-night function was partly radar artefacts. The per-radar patterns (`AS_G_radar_patterns.csv`) then become a topic of their own: an azimuthal calibration of each radar.

**3. One fixed function, or a fit per night?**
- Analysis: section 5, leave one night out.
- Rule:
  - Fixed function when the climatological fit gets at least 90% of the per-night fit's reduction of the disagreement, and leaves |c2| and |c1| below 0.03 on at least 75% of the nights.
  - Otherwise, fit per night in the pipeline: clean without shape, bin the overlap, fit, then apply `as_apply_aspect`. This costs about 1 min per frame.

**4. Which covariates matter?** (height, speed, season, longitude)
- Analysis: section 4, a two-stage regression; τ is the between-stratum sd that no covariate explains.
- Rule: keep a covariate when |estimate| > 2 se and it lowers τ clearly.
  - Test on two pseudo-nights cut from the one night: τ for c2 fell from 0.085 (height only) to 0.027 (all covariates).

**5. Below 0.8 km**
- Analysis: section 6, single-radar harmonics with a fixed effect per radar, sweep and range band, over all nights (`as_time_fit`).
- Rule: where single radars and the overlap both have values (≥ 0.8 km) and agree within 2 se, use the single-radar values below 0.8 km. This replaces the 1 km ramp, which is a choice now.

**6. Axis: heading or track?**
- Analysis: the `b2_turned` and `b1_turned` columns of section 6 (an axis turned from the track).
- Rule: only a first hint. A real test needs the winds (heading = track − wind), for example ERA5, which needs the user's own CDS account. This is not part of this run.

## Data and budget

**The nights.** `AS_E_select_nights.m` reads the Dark Ecology daily archive (2021–2025, 107 us_east radars; no download):
- 192 dry, strong nights qualify;
- 16 are chosen in turn from each season and 30° track sector:
  - 6 spring nights, with tracks 354–38°;
  - 10 autumn nights, with tracks 169–244°;
- their track axes span about 75°, and several nights have very different tracks across the domain (axis coherence 0.2–0.5);
- the 21–22 Oct 2026 night is added from `aspect/bins` (6 frames, 01–06 UTC).

| Item | Estimate |
| --- | --- |
| Volumes | 16 nights × 6 hourly frames × 107 radars = 10,272 |
| Download | ≈ 90 GB (8.8 MB per volume); the pilot ≈ 5.7 GB |
| Cleaned volumes | ≈ 21 GB (2.1 MB each), cache `prod/volumes/493fd773` |
| Overlap bins | ≈ 0.5 GB |
| Cleaning | 6.3 s per volume on 4 MistNet workers: ≈ 70 min per night, ≈ 18 h in all |
| Binning | 19 s per frame (measured), ≈ 2 min per night |
| Download time | unknown; the pilot measures it |
| `AS_G_multi_night_fit.m` | ≈ 1 h and ≈ 10 GB of memory for ~12 M cell pairs (101 s for 1.4 M pairs on the test) |

**Disk.** 151 GB were free on 5 Nov 2026. Keeping all the raw files leaves about 40 GB. `delete_raw = true` deletes, after each night is binned, only the raw files the job downloaded. It needs the user's OK.

**Memory.** 4 MistNet workers at 3.5–6 GB each, plus the job's own worker, stay at or below about 27 GB of the 36 GB. While the job runs, run nothing else large: no other parallel pool, no inversion, no other cleaning.

## Steps

1. **Ask the user** for:
   - the OK to download (the pilot about 5.7 GB, then about 85 GB);
   - whether to delete the downloaded raw files after binning;
   - whether the machine is free for about 20 h.
2. **Choose the nights.** Run `scripts/AS_E_select_nights.m` (seconds) and look at `figures/AS_E_nights.png`.
3. **Pilot.** Run section 1 of `scripts/AS_F_multi_night_run.m`, a batch job on the strongest night (about 1.5 h).
   - Do not `wait()` in the MATLAB MCP: the tool times out.
   - Check `pilot_job.State` and `run_log.csv` from time to time.
4. **Check the pilot** (section 2 of AS_F). Go on only if all of these hold:
   - cleaning takes at most about 10 s per volume;
   - at most 5% of the volumes failed;
   - there are at least 300 k cell pairs;
   - b2_per_km and b1 are within a factor of 2 of the one-night values.

   Otherwise stop and report.
5. **Main run.** Section 3 of AS_F starts a batch job on the other nights (about 20 h). It can be restarted: a night with `done.txt` is skipped. Check `run_log.csv` every few hours.
6. **Fit.** Run `scripts/AS_G_multi_night_fit.m` after the job has finished. Through the MCP, run it with `cd(project_root)` and `run(...)`. If the tool times out, check that `multi_night_results.mat` was written.
7. **Write up.** Add a section "Multi-night fit" to `README.md`: the tables of AS_G, the figure `figures/AS_G_multi_night.png`, and the decision of each rule above.
8. **Report** to the user: the decision, and what it means for the production cleaning. Nothing goes into the toolbox without the user's approval.

## Rules

- Nothing in `toolbox/`, `examples/` or the manuscript changes, and nothing is committed.
- MATLAB runs through the MCP tools, from the project root; never with `matlab -batch`.
- The pool stays at 4 workers or fewer. Never call `restoredefaultpath`.
- Delete only raw files this job downloaded, and only with the user's OK.
- Other sessions share the MATLAB session. Keep interactive calls short while the job runs.

## What was tested (5 Nov 2026)

- **Binning:** `as_night_bins` reproduces the existing bins of the us_east night exactly (frame 10: 239,880 cells, maximum difference 0) in 19 s.
- **Settings:** the night settings of `as_multi_night_job` give production cleaning without the shape correction: MistNet, the gate estimator, cleaning version 2 (with the beam-blockage fix), and cache `493fd773`. On 2026-10-22 04–05 UTC, 200 of 214 night slots were local.
- **Fit script:** `AS_G_multi_night_fit.m` runs end to end on two pseudo-nights cut from the 21–22 Oct night, plus that night: 1.38 M cell pairs, 101 s.
  - Its numbers mean nothing: the pseudo-nights are the same night.
  - The run did show the radar-pattern confounding described above.
- **Not tested:** the batch job end to end (download and cleaning). The pilot is that test.

## Sharing with other studies

The cleaned volumes (cache `493fd773`) and each night's manifest (`prod/aspect_YYYYMMDD/manifest.mat`) are ordinary production-style caches. The near-radar study (`studies/near_radar/PLAN.md`) can compute its own per-volume diagnostics on them after the job, with no new download or cleaning.
