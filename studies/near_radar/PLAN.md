# Near-radar dip: drop the steep cuts near each radar

Plan for an agent or the user. Everything here is dev, in `studies/near_radar/`. Nothing goes into `toolbox/`, `examples/` or the manuscript, and nothing is committed, without the user's approval.

The goal is to remove the near-radar dip (column 0.72 of the background within 10 km on us_east) by **dropping** the observations that read low, not correcting them, if the data around them in space and time can fill the gap. That needs a rule that works across many radars and nights.

## Do we have enough data?

**To choose the criterion: yes, for a first pass.** The data on disk cover:
- the us_east night of 21–22 Oct 2026: ~101 radars × 24 frames, raw files plus the production cleaning (cache `9c2b4ef2`);
- the night of 13 May 2021: 28 Midwest radars × 16 frames, raw files plus cleaning (cache `fad5e1e1`, an older estimator cleaning with the density cap and patch filter);
- 9 sandbox nights (2021–2025) × 28 radars, one volume each, raw files plus `fad5e1e1`.

In all: ~3,100 cleaned volumes, ~130 radars, 11 nights in both seasons, VCPs 12/31/32/34/35/212/215, and staggered-PRT cuts at 7 radars. That is enough to tell which cuts read low and how the deficit falls off with range.

**To show the rule holds across nights and seasons: no.**
- Only two nights have full time series, one per season, and they already differ: the adaptive R_min has a median of 20 km on the autumn 2026 night and 5 km on the spring 2021 night.
- The 9 sandbox nights have one volume per radar. That is too few for a per-radar-night R_min, and no use for a solve.
- The gap-filling test (a solve) exists only for these two nights.

**The shared night set on the cluster** (`studies/nights/PLAN.md`) is the data for that: 18 full nights (17 Dark Ecology nights 2021–2025 in six migration phases, weak and strong, plus the 21–22 Oct 2026 reference), the 107 us_east radars every 30 min, ~38,000 volumes cleaned with production settings without the shape correction (cache `493fd773`), raw files kept. It needs no download and no cleaning.

## What the light analyses found (5 Nov 2026)

### 1. Scan strategies — `NR_A_inventory.m`
1,005 volumes, from the Level II header (Message 5, waveform of every cut) by `nr_read_vcp`. Output: `figures/NR_A_vcp_summary.csv`.

| VCP | Split cuts (surveillance reflectivity) | Other cuts below 7° | Volumes |
| --- | --- | --- | ---: |
| 35, 212, 215, 12 | 0.5, 0.9, 1.3° | batch from 1.8°; staggered pulse pair at 7 radars in 2026 | 950 |
| 34 | 0.5, 0.9, 1.3° | batch and contiguous Doppler | 21 |
| 32 | 0.5, 1.5° | batch from 2.5° | 84 |
| 31 | 0.5, 1.5, 2.5° | contiguous Doppler 3.5, 4.5° | 14 |

- **"Elevation ≥ 1.8°" is the batch boundary only for VCP 12/34/35/212/215.**
  - For VCP 31 it would drop a split cut (2.5°) and two Doppler cuts.
  - Staggered-PRT radars (KBIS, KBUF, KCXX, KFFC, KICT, KMHX, KMPX) have no batch cuts above 1.8°.
- **23 of 129 radar-nights change VCP during the night,** so the cut type must be read per volume, not per radar.

### 2. Each cut against the split-cut profile, no solve — `NR_B_volume_ratio.m`
Light run: 723 volumes (every 6th of the two full nights, plus all the sandbox nights). `nr_volume_ratio` fits a beam-weighted vertical profile to the split cuts at 20–60 km of the same volume. It then compares every sweep and range ring with what that profile predicts there.

Observed / profile, pooled, by cut type:

| Ring (km) | 3–5 | 5–7.5 | 7.5–10 | 10–12.5 | 12.5–15 | 15–20 | 20–25 | 25–30 | 30–40 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| batch | 0.57 | 0.66 | 0.74 | 0.80 | 0.83 | 0.85 | 0.87 | 0.89 | 0.91 |
| staggered | 0.46 | 0.60 | 0.73 | 0.83 | 0.90 | 0.93 | 0.97 | 0.97 | 0.96 |
| Doppler (VCP 31, 34) | 0.91 | 1.06 | 1.10 | 1.12 | 1.11 | 1.11 | 1.13 | 1.11 | 1.04 |
| split | 0.78 | 0.89 | 0.95 | 1.03 | 1.07 | 1.07 | 1.06 | 1.04 | 1.01 |

**It is range, not height.** At a fixed height the batch ratio still rises with range: at 600–900 m it is 0.51, 0.64, 0.77, 0.89, 0.95 and 0.96 at 5, 7.5, 10, 15, 20 and 25 km. So the right shape for the rule is a minimum range.

**Two levels are mixed in the batch row above** (`nr_radar_night`, per radar-night):
- **Near-radar deficit:** each ring relative to the radar-night's own batch level at 25–50 km. Median over radar-nights:

  | Ring (km) | 3–5 | 5–7.5 | 7.5–10 | 10–12.5 | 12.5–15 | 15–20 | 20–25 | 25–30 | 30–40 |
  | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
  | Deficit | 0.52 | 0.62 | 0.71 | 0.81 | 0.85 | 0.90 | 0.96 | 0.98 | 0.99 |

  The spread over radar-nights is wide (10–12.5 km: 0.46–1.21 for the 10–90%), partly noise.
- **Far batch level:** batch cuts at 25–50 km over split cuts at 10–40 km, median 0.79 (10–90%: 0.52–1.03; pooled sums ~0.9). This offset does not depend on range, so dropping near rows does not touch it. It is a separate topic (see `SB_cut_levels`) and is flagged here, not fixed.

**Adaptive R_min.** Outer edge of the outermost ring within 25 km whose deficit is below 0.95 by more than 2 se (jackknife over volumes), on radar-nights with ≥ 3 volumes:
- us_east (autumn 2026, 93 radar-nights): median 20 km (0–25);
- 13 May 2021 (spring, 18 radar-nights): median 5 km (0–25);
- VCP 34: 0–5 km.

Without the 5% tolerance (`r_min_significant_km`) the R_min grows with the amount of data, so the tolerance stays.

**Staggered cuts** (7 radar-nights) are mixed: a strong deficit at KCXX and KBIS, none beyond 10 km at KMPX, KICT and KMHX.

**Caveat.** Above about 1.2 km the reference profile rests on wide far beams, so absolute levels at the top are uncertain by about ±10%. The range shape at a fixed height is the robust part.

### 3. The clutter-filter flag per gate — `NR_C_clutter_filter.m`
Smoke test only: 3 volumes (KDMX, KMPX, KDVN), so not conclusive.
- **Batch gates read low in every CFP class,** including "filter applied, 0 dB removed": 0.54–0.69 at 7.5–15 km.
- **The filter acts on almost every batch gate at every range:** "not applied" (code 2) is below 10% of gates.
- **Staggered gates read low only where ≥ 4 dB were removed.**

If the full run confirms this, a per-gate clutter flag cannot isolate the batch deficit.

### 4. Drop rules in the solve — `NR_D_drop_solves.m`
Smoke test: `ar_midwest_shared`, 7 frames, with production settings, except the prior and weights are rebuilt without the dropped rows. The split entries within 20 km of every second radar are held out of every solve as an independent check.

| Rule | Entries dropped | Column 0–10 km | Column 10–20 km | Dropped rows obs/pred | Held-out split obs/pred | Misfit 30–120 km |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| all rows | 0 | 0.83 | 0.89 | – | 1.30 | 0.165 |
| batch < 20 km | 2.7% | 0.98 | 0.99 | 0.70 | 1.20 | 0.165 |
| batch, adaptive R_min | 3.6% | 0.99 | 1.01 | 0.72 | 1.21 | 0.165 |

- AR_J's drop gave 0.997 at 0–10 km because its prior still contained the dropped rows. Rebuilt without them, it is 0.98.
- The held-out split entries read 1.2–1.3 under every rule. The lowest ~300 m near a radar is denser than the other data predict, and the hold-out itself removes that information. Compare the rules here, not the level.

**Synthetic test** (same case; truth T1 = the "batch < 20 km" field; observations F × truth, no noise), with the batch rows within 20 km dropped:
- **T1:** reconstructed / true column 1.008 at 0–10 km and 1.006 at 10–20 km; median error 1.5%, the same as far from the radars.
- **T1 with the measured batch deficit** on the batch entries: the same (1.002, 1.003).
- **T2 = T1 × ±30% structure** at a 20 km scale, new each frame: no bias (1.02, 1.01), but a median column error of 11% at 0–10 km and 9% at 10–20 km, against 3% beyond 30 km.

So **smooth fields fill the gap.** Local structure smaller than the gap and not carried from frame to frame does not; dropping trades it for the removal of a 30–50% bias. The full run adds the "all rows" and "batch < 30 km" rules to this test.

## Candidate criteria

| Criterion | For | Against | Status |
| --- | --- | --- | --- |
| Elevation ≥ 1.8° | Simple; no header needed | Wrong for VCP 31/32 and for staggered-PRT radars | Fallback only |
| **Cut type from the header** (batch; possibly also staggered) | Where the step is; correct for every VCP; per volume | Needs the Level II header (0.1 s per file, metadata record only) | **Recommended** |
| **Fixed minimum range** | One number; the deficit is a function of range at a fixed height | One setting for nights that differ (20 vs 5 km) | Candidate: 20 or 25 km |
| Minimum range by elevation or beam height | – | At a fixed height the deficit still depends on range | Rejected by the data |
| Clutter flag per gate (CFP) | Physical, no range threshold | Early evidence: batch gates are low with 0 dB removed; needs CFP in the cleaning (dev patch: `studies/sandbox/dev/cfp_error_toolbox.patch`) | Test in run 1 |
| **Adaptive R_min per radar-night** | Follows the night differences; cheap (no solve, ~20 volumes a night) | One tolerance (5%); needs ≥ 3 volumes per radar-night | Candidate |

## Run 1: about 1 h 15 min, ready to launch

Uses only data on disk: no cleaning, no download, no pool. Run from the project root through the MATLAB MCP. The scripts are set for the full run (`every_kth = 1` in NR_B; the light results above used 6).

```matlab
project_root = setup_paths();
cd(project_root); addpath(project_root)
run(fullfile(project_root, "studies", "near_radar", "scripts", "NR_B_volume_ratio.m"))    % ~15 min, ~3,100 volumes
run(fullfile(project_root, "studies", "near_radar", "scripts", "NR_D_drop_solves.m"))     % ~50 min, 22 solves
run(fullfile(project_root, "studies", "near_radar", "scripts", "NR_C_clutter_filter.m"))  % ~10 min, ~180 volumes
```

- **NR_B first:** NR_D reads its `NR_B_radar_night.csv` (adaptive R_min) and `NR_B_pooled.csv` (the deficit put on the synthetic data).
- **NR_D** runs the rules below:

  | Case | Rules |
  | --- | --- |
  | ar_midwest_shared | all rows, elevation < 20 km, batch < 10/20/30 km, batch + staggered < 20 km, batch adaptive, batch at all ranges |
  | night_20210512 | all rows, elevation < 20 km, batch < 20/30 km, batch adaptive |
  | synthetic (T1, T1 + deficit, T2) | all rows, batch < 20 km, batch < 30 km |

- **Timeouts:** MCP calls time out before MATLAB finishes. Check that `figures/NR_*_*.csv` and `<local_data_dir>/near_radar/NR_D_results.mat` were updated before assuming a failure.
- **After run 1,** write the results into a `README.md` here and set `cut_set` and `fixed_range_km` at the top of `NR_E_multi_night.m` (decision rule, steps 1–3).

## Run 2: on the cluster, ~1 h + ~3 h, on the shared night set

Runs as cluster jobs, submitted from the laptop in the repository with `hpc/submit.sh` (SSH to merula; see `hpc/README.md`). The heavy parts are two functions; each writes one file per night (and region) to `products/near_radar/NR_E/` and skips files that exist, so a cut-off job is finished by submitting it again.

Set `cut_set` and `fixed_range_km` from run 1 in the commands below and at the top of `NR_E_multi_night.m`.

**Job 1, the per-volume diagnostic** (`nr_diagnose_nights`, ~38,000 volumes on 8 workers, ~1–1.5 h):

```bash
hpc/submit.sh -c 10 -m 32G -t 0-06:00:00 near-radar-diag 'addpath(fullfile(project_root, "studies", "near_radar", "functions")); n = readtable(fullfile(project_root, "studies", "nights", "figures", "NS_A_nights.csv"), TextType="string"); nr_diagnose_nights(n.name(startsWith(n.name, "night_us_")), 8)'
```

**Jobs 2–3, gap filling,** one per region, after job 1 (`nr_gap_fill`: per night, a regional case from the cleaned volumes, `build_observations` with its 6 workers, then 4 solves; ~10 min per night, ~3 h per job):

```bash
hpc/submit.sh -c 8 -m 48G -t 1-00:00:00 near-radar-midwest 'addpath(fullfile(project_root, "studies", "near_radar", "functions")); n = readtable(fullfile(project_root, "studies", "nights", "figures", "NS_A_nights.csv"), TextType="string"); nr_gap_fill(transpose(n.name(startsWith(n.name, "night_us_"))), "midwest", "batch", 20)'
```

and the same with `near-radar-east` and `"east"`.

**Then on the laptop** (share mounted), `NR_E_multi_night.m` finds every file done and only collects, analyses and plots (Helvetica):
1. the deficit by ring and cut type, with its spread over radar-nights;
2. a two-way radar + night model of the 5–15 km deficit (`NR_E_variance.csv`);
3. R_min by season, phase, year and VCP, and its Spearman correlation with track speed and density;
4. the share of radar-nights where the fixed range leaves a low ring;
5. the gap-filling table over 18 nights × 2 regions (Midwest: the 28 sandbox radars; East: Appalachians to the coast), four rules: all rows, the fixed range, adaptive, the whole cut set;
6. the decision table, `figures/NR_E_decision.csv`.

Each regional `case.mat` (1–3 GB) is deleted after scoring unless `keep_cases = true`. Jobs 2 and 3 open 6 workers each: run at most these two at once, and check `seff` (or the job log) after the first night.

## Decision rule: from the results to a production setting

1. **Which cuts.**
   - Drop the **batch** waveform, read from the header.
   - Add **staggered** cuts if their deficit at 5–15 km (pooled per radar-night, `nr_radar_night(rows, "staggered")`) is below 0.95 by more than 2 se in at least half of their radar-nights.
   - Doppler cuts stay unless they show the same.
   - Elevation ≥ 1.8° is used only where no header can be read.
2. **Clutter flag.** Keep it as a candidate only if, in NR_C, batch gates with 0 dB removed read ≥ 0.95 at 5–15 km while the gates with removed power carry the deficit. Otherwise drop it from consideration.
3. **Fixed range R\*.** The smallest of 10/20/30 km beyond which a larger range no longer changes the near-radar column (0–20 km over 60–100 km) by more than its standard error over night-regions. The tolerance is the noise, not a chosen number.
4. **Fixed or adaptive.** Adaptive replaces the fixed R\* only if all three hold:
   - over the NR_E night-regions, it brings the near-radar column of the worst 10% closer to 1 than R\* does;
   - it drops no more entries in median;
   - the radar + night model shows the deficit varies mostly by night.

   If it varies mostly by radar, a static per-radar range table from the multi-night run is the alternative.
5. **Gap filling holds when:**
   - the misfit at 30–120 km is unchanged (within 1%);
   - the held-out split ratio moves toward 1 compared with "all rows";
   - in the synthetic test, the column within the gap is unbiased (T1 and T2 ratio within the T2 error) and its T1 error is no more than 1.5× the error beyond 30 km.
6. **Production (needs the user's approval).** A proposal, not a decision:
   - store each sweep's waveform (Message 5) with the cleaned volume;
   - drop the chosen cuts within R\* (or the adaptive R_min) when the observation rows are built (`build_radar_rows`), not in the cleaning, so the setting can change without cleaning again;
   - the cleaned cache keeps all data.

## What the shared night set provides

Nothing more is needed from it: the 18 nights are full (every 30 min), the raw files are kept (waveforms from the headers, and CFP per gate for NR_C on two nights if run 1 asks for it), and the cleaning matches production without the shape correction.

## Files

| File | What it does |
| --- | --- |
| `functions/nr_read_vcp.m` | VCP and waveform of every cut from the Level II header (metadata record only) |
| `functions/nr_volume_ratio.m` | The cheap diagnostic: every sweep and ring against the beam-weighted split-cut profile of the same volume (optionally per gate) |
| `functions/nr_radar_night.m` | Per radar-night: deficit by ring (jackknife se), far level, adaptive R_min |
| `functions/nr_cut_types.m` | Waveform of every observation entry (per time, from the raw files) |
| `functions/nr_drop_rows.m` | Entries a rule drops (elevation / batch / batch + staggered / not split; fixed, per-radar or infinite range) |
| `functions/nr_apply_drop.m` | Removes entries, then rebuilds the prior and weights without them |
| `functions/nr_score.m` | Column against distance, obs/pred by cut, ring and drop, misfit; against a truth |
| `functions/nr_night_case.m` | Regional case from any cleaned night's folder (as `ar_test_case`) |
| `functions/nr_shared_path.m` | Maps share paths between the cluster and the laptop's mount |
| `functions/nr_diagnose_nights.m` | Cluster job: the diagnostic of every volume of shared nights |
| `functions/nr_gap_fill.m` | Cluster job: the four rules solved on regional cases of shared nights |
| `scripts/NR_A_inventory.m` | VCP inventory of the raw files on disk (~2 min) |
| `scripts/NR_B_volume_ratio.m` | Diagnostic over the cleaned volumes on disk (~15 min) |
| `scripts/NR_C_clutter_filter.m` | Deficit by clutter-filter class per gate (~10 min) |
| `scripts/NR_D_drop_solves.m` | Drop rules in the solve, real and synthetic (~50 min) |
| `scripts/NR_E_multi_night.m` | Run 2: drives or collects the two jobs, then the analysis and decision table |

**Checks done:**
- `check_matlab_code`: no issues in any file.
- Smoke tests:
  - NR_A and NR_B in full (light thinning);
  - NR_C on 3 volumes;
  - NR_D reduced to 3 rules on `ar_midwest_shared` plus 3 synthetic solves (654 s);
  - NR_E sections 1–2 (the earlier version) on 2 nights × 6 radars;
  - the cluster job functions on the cluster (job 58735, 12 Nov 2026, 5 min, exit 0; pilot files `*_pilot.mat`, which the full run ignores):
    - `nr_diagnose_nights` on 2 nights × 3 frames: 404 of 404 volumes, 22 s per night on 4 workers (~0.1 s per volume, so ~1 h for all ~38,000);
    - `nr_gap_fill` on 2023-05-05, Midwest, first 4 frames (3 after the sun limit): case built (26 of 28 radars, 42 s), four solves of 36–45 s.
    - Pilot numbers, only to show the chain runs (3 frames; the adaptive rule fell back to the 20 km median because no Midwest radar had 3 volumes yet):

      | Rule | Dropped | Column 0–10 km | Column 10–20 km | Held-out split | Misfit 30–120 km |
      | --- | ---: | ---: | ---: | ---: | ---: |
      | all rows | 0 | 0.89 | 0.93 | 0.97 | 0.215 |
      | batch < 20 km | 3.1% | 0.96 | 0.97 | 0.95 | 0.224 |
      | whole batch cuts | 21% | 0.95 | 0.96 | 0.94 | 0.206 |

    - Peak memory 25 GB (MaxRSS, with the 6 `build_observations` workers): book at least 32 GB for the gap-fill jobs.
  - A first cluster test found and fixed a bug: `bmm4d.config` drops frames before the sun limit, so `nr_night_case` now keeps only the manifest columns of the case's times.

## Rules

- Dev only. Nothing in `toolbox/`, `examples/` or the manuscript changes, and nothing is committed without the user's OK.
- Run 1: MATLAB through the MCP from the project root, never `matlab -batch`; never `restoredefaultpath`; no pool. Run 2: cluster jobs through `hpc/submit.sh` only, at most two at once, and only with the user's go.
- **Keep all biological echo.** The rule drops a measurement artefact by cut type and range, the same whatever the echo is (birds, insects, local movements).
