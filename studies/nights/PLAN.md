# Night set on the cluster: which data, and how to process it

One set of full nights, downloaded and cleaned once on the cluster (merula), serves every open question of the v2 studies. This plan says which nights, why, and how to run it. Nothing has been downloaded yet: the user's go-ahead is needed first.

Dev only. Nothing in `toolbox/`, `examples/` or the manuscript changes, and nothing is committed.

## What the studies need

Each line below is a finding from the dev studies (5–6 Nov 2026), what blocks it, and the data that unblocks it.

**Aspect** (`studies/aspect`: `README.md`, `PLAN_multi_night.md`, `PLAN_robust_aspect.md`)
- **The blocker:** one night cannot separate the birds' aspect from fixed patterns of each radar. With a fixed pattern per radar in the fit, c2 falls from 0.092 to 0.017 per km.
  - The axis (heading or track) is unknown.
  - Below 0.8 km there is no gradient-free measurement.
- **The data:**
  - nights whose tracks are turned 30–90° against each other (the side-on term);
  - both seasons (the front/back term);
  - winds for the heading and airspeed;
  - a few days of widespread rain, where the radar patterns show without birds.

**Near-radar dip** (`studies/near_radar/PLAN.md`)
- **The blocker:** the batch and staggered cuts read 0.57–0.85 of the profile within 20 km. The deficit depends on range, not height. The cut type must be read per volume, because 23 of 129 radar-nights change their scan strategy during the night. An adaptive minimum range was 20 km in autumn 2026 and 5 km in spring 2021.
- **The data:**
  - many full radar-nights in both seasons, weak and strong (clear-air and precipitation scan strategies);
  - the raw files kept, for the headers and the clutter-filter flag;
  - full nights for the gap-filling solves.

**Terrain** (`studies/artifacts`, `AR_K_terrain_hybrid.m`; that session stopped at a spend limit, so its results are saved but not written up)
- **The blocker:** the signals disagree.
  - The radar rows line up best with height above sea level, or above terrain smoothed over ≥ 150 km. They line up worst with height above the ground.
  - But held-out mountain radars are predicted best with height above the ground:

    | Coordinate | Relative error | Log error |
    | --- | ---: | ---: |
    | Above sea level | 0.508 | 0.568 |
    | Smoothed terrain, 100 km | 0.497 | 0.535 |
    | Above ground | 0.484 | 0.510 |

  - Detours change nothing.
- **The data:** strong, dry Appalachian nights with tracks across the ridges and along them.

**Prior** (`studies/artifacts/PRIOR_OPTIONS.md`)
- **Where it stands:**
  - In the paper, "prior" is the mean field; the GP changes the covariance.
  - White prior, GP and smoothness are the three terms of one Matérn (SPDE) square. Option B, one Matérn term around the profile mean, is the recommendation.
  - The 7× cost comes from the diagonal preconditioner.
- **The data:**
  - stage 0: profile anomalies over many radars and nights, for the horizontal range;
  - edge tests: strong, dry nights over the sparse Plains radars.

**Paper rerun with v2** (the paper reruns with v2; no ports to 1.0.x)
- **The data:** the Key West and Midwest cases on the cluster.

## The set — `scripts/NS_A_select_nights.m`

**Source.** The Dark Ecology daily archive, 2021–2025: 930 nights in spring and autumn over the 107 us_east radars. No download.

**Qualifying nights: 234.** A night qualifies when:
- at least 90 radars report;
- rain falls on less than 10% of the night at 70% or more of them;
- the median night-mean reflectivity is at least 100 cm² km⁻².

**Choosing.** Three nights in each of six migration phases:
- spring: 15 Mar–20 Apr, 21 Apr–20 May, 21 May–10 Jun;
- autumn: 15 Aug–14 Sep, 15 Sep–20 Oct, 21 Oct–20 Nov.

The first night of each phase is the strongest. Each next one is the night farthest from all chosen nights in track axis, track direction, speed and reflectivity, with at least 4 days between chosen nights and at most 5 per year. Then nights would be added until the terrain and edge questions have enough. They already did, so none was added.

| Night (local date) | Phase | Track (deg) | Speed (m/s) | Median refl. | Also serves |
| --- | --- | ---: | ---: | ---: | --- |
| 2021-05-12 | spring peak | 345 | 10.4 | 344 | terrain across, edges |
| 2021-05-25 | spring late | 33 | 12.4 | 157 | |
| 2021-08-27 | autumn early | 258 | 7.3 | 208 | |
| 2021-10-16 | autumn peak | 178 | 17.2 | 567 | terrain across |
| 2022-09-07 | autumn early | 198 | 9.6 | 971 | |
| 2022-10-10 | autumn peak | 210 | 7.8 | 128 | |
| 2022-11-07 | autumn late | 213 | 11.6 | 245 | terrain along |
| 2023-04-10 | spring early | 10 | 11.0 | 169 | |
| 2023-05-05 | spring peak | 4 | 14.9 | 1025 | |
| 2023-05-22 | spring late | 3 | 11.1 | 445 | edges |
| 2024-04-14 | spring early | 17 | 13.7 | 289 | terrain along, edges |
| 2024-05-10 | spring peak | 15 | 7.8 | 342 | edges |
| 2024-08-27 | autumn early | 163 | 5.4 | 274 | terrain across |
| 2024-10-23 | autumn late | 185 | 9.9 | 524 | |
| 2025-04-17 | spring early | 14 | 18.4 | 306 | |
| 2025-10-08 | autumn peak | 200 | 12.9 | 1306 | terrain along |
| 2025-10-30 | autumn late | 170 | 11.9 | 167 | terrain across |
| 2026-10-21 | reference (every dev study so far) | | | | |

**How the set covers the questions**
- **Track axes** span about 95° (−17° to 78° from N–S), with speeds of 5–18 m/s and reflectivity of 128–1306 cm² km⁻².
- **Varied tracks within a night:** several nights have tracks that differ across the domain (axis coherence 0.08–0.35). These help the aspect fit too.
- **Every full night serves** the aspect, near-radar and prior studies.
- **Terrain across the ridges:** the track is at least 70° (|sin| ≥ 0.7) from the ridge axis (NE–SW), on nights that are dry and strong over the Appalachians.
- **Terrain along the ridges:** the track is within about 17° (|sin| ≤ 0.3) of the ridge axis, on nights that are dry and strong there.
- **Edges:** nights that are dry and strong over the Plains radars.

**Rain days** (radar patterns without birds). The three wettest April–May and September–October days, in different years: 11 Apr 2024, 5 Apr 2025 and 28 Oct 2021. Six hourly frames, 15–20 UTC. On these days 28–32% of the radars have rain over half the day or more.

**Paper cases.** Key West (23 Apr 2018) and Midwest (16 May 2020).

**Frames.** Every 30 min, from 1 h after sunset at the east edge to 30 min before sunrise at the west edge (21–28 frames). Each radar keeps only its night volumes.

**Output.** `figures/NS_A_nights.csv`, which drives the cluster jobs, and `figures/NS_A_nights.png`.

## Processing

- **Download once, keep the raw files** in the cluster's raw archive (`weather_radar_us/nexrad_level2`). There are 116 TB free. The headers, the clutter-filter flag and any later cleaning variant then need no second download.
- **Clean the full nights once,** as production but without the shape correction (read option `shape` false). The result is one shared cache, the same settings as the multi-night aspect plan (`volume_cache_dir` hash `493fd773` on the laptop), plus the static clutter masks.
  - The shared aspect function is applied afterwards (`as_apply_aspect`, minutes per night), without cleaning again.
  - The paper's ellipse cleaning is not repeated unless a comparison needs it.
- **Rain days and paper cases are only downloaded.**
  - The radar-pattern analysis reads the rain reflectivity directly; the cleaning would remove it.
  - The paper cases are cleaned once the v2 method is chosen.
- **Winds** (not part of the jobs below): ERA5 on pressure levels, 1000–600 hPa (u, v, t, z), hourly, 22–52°N, 106–64°W, for each night, about 30 MB a night. Two ways to get them:
  - the CDS API with the `cds_token` already in the environment (by variable name only);
  - NCAR's open ERA5 copy on AWS.

  Either way the cluster needs a small Python environment with `cdsapi` (a container or `pip --user`).
- **Where the products land:** `BMM-4D/products/prod/<name>` (manifest, static clutter) and `BMM-4D/products/prod/volumes/<hash>` on Ceph. The laptop sees the same files under `/Volumes/UNIT_Vogelzug`.

## Cluster check (6 Nov 2026)

**Access and nodes**
- `ssh merula` works (key login).
- Partitions: `short` (1 day), `normal` (7 days) and `long`.
- merops-02 has 352 CPUs and 469 GB. At the check, others used 128 CPUs.

**Storage**
- The UNIT_Vogelzug share (Ceph) has 159 TB, of which 116 TB are free.
- The raw archive holds 3.3 GB (33 radar-days).

**Tooling**
- `hpc/submit.sh` snapshots the laptop's working tree, the untracked studies included, into a job folder with provenance.
- The cluster's `local_settings.json` points the raw files and products to the shared storage.

**Today's sandbox job**
- It downloaded 168 volumes in 30 s (about 50 MB/s).
- It cleaned them with MistNet in 1084 s on 4 workers: 26 CPU-seconds per volume, as on the laptop.
- Observations, the solve (37 iterations) and the maps also ran. The first try failed only at the MP4 export, which is not available on Linux.

**Downloads.** The NEXRAD bucket answered from the cluster (HTTP 200, 0.34 s).

**Not checked:** whether two MATLAB jobs can run at the same time under the individual licence. The pilot tests it.

## Budget

| | Volumes | Raw | Cleaning |
| --- | ---: | ---: | ---: |
| 17 nights + reference | 39,700 | 350 GB | 286 CPU-hours |
| 3 rain days | 1,930 | 17 GB | – |
| 2 paper cases | 350 | 3 GB | – |
| **Total** | **41,900** | **370 GB** | **286 CPU-hours** |

- **Cleaned volumes:** about 84 GB (2.1 MB each).
- **Cleaning time:** `clean_volumes` uses 4 MistNet workers; this is fixed in the code, so no toolbox change is needed. That makes about 6.5 s per volume per job.
- **Jobs:** six jobs of three nights each (balanced to about 6,600 volumes) take about 12 h each, run side by side.
  - Together they use 48 CPUs and about 240 GB (40 GB each; the pilot measures the real peak).
  - With four jobs at once instead, the whole run takes about 18 h.
- **Downloads** run inside the jobs, at about 7 min per night.

## Steps

**0. The user's go.**
- A download of about 370 GB into the cluster's raw archive.
- Up to 48 CPUs and 240 GB on merops-02 for about 12 h. Fewer jobs at once if the node is busy.

**1. Pilot, about 1 h.** From the laptop, in the repository: the first 6 frames of the strongest spring night, and at the same time one rain day (download only), to test that two jobs can run under the licence.

```bash
hpc/submit.sh -c 8 -m 40G -t 0-04:00:00 nights-pilot 'addpath(fullfile(project_root, "studies", "nights", "functions"), fullfile(project_root, "examples")); ns_process_nights("night_us_20230505", 6)'
```

```bash
hpc/submit.sh -c 2 -m 8G -t 0-02:00:00 nights-rain1 'addpath(fullfile(project_root, "studies", "nights", "functions"), fullfile(project_root, "examples")); ns_process_nights("rain_us_20240411")'
```

Check:
- `products/nights/run_log.csv`;
- the job logs;
- `sacct -j <id> --format=JobID,Elapsed,MaxRSS,State` (peak memory per job);
- seconds per volume, and failed volumes.

Go on if both jobs ran, cleaning takes at most about 10 s per volume per job, and no more than 2% of the volumes fail.

**2. Main run.** Six jobs. The pilot's volumes are reused from the cache.

| Job | Nights |
| --- | --- |
| nights-1 | 20221107, 20230410, 20210525 |
| nights-2 | 20241023, 20240414, 20230505 |
| nights-3 | 20251030, 20240827, 20230522 |
| nights-4 | 20211016, 20261021, 20240510 |
| nights-5 | 20221010, 20210827, 20250417 |
| nights-6 | 20251008, 20220907, 20210512 |
| nights-extra | the other two rain days and the two paper cases, download only (`-c 2 -m 8G`) |

Each night job is submitted like the pilot, with a time limit of 1 day:

```bash
hpc/submit.sh -c 8 -m 40G -t 1-00:00:00 nights-1 'addpath(fullfile(project_root, "studies", "nights", "functions"), fullfile(project_root, "examples")); ns_process_nights(["night_us_20221107" "night_us_20230410" "night_us_20210525"])'
```

A night that was cut off is finished by submitting it again: cached volumes are skipped.

**3. Then each study**, with its own plan: about 1 h first, then the long part. They read the shared products; nothing is cleaned again.

- **Aspect:**
  - bin the overlap of every night (`as_night_bins`; about 20 s per frame);
  - run the multi-night fit (`AS_G_multi_night_fit.m`, set to read this night set instead of the AS_E list; hourly or all frames);
  - the heading from the winds (`PLAN_robust_aspect.md` A2);
  - radar patterns from the rain days (B1).
- **Near-radar:**
  - the per-volume diagnostic on every night (`NR_B`);
  - the gap-filling solves, Midwest and East (`NR_E`);
  - the clutter-filter flag on two nights from the kept raw files.
- **Terrain:** the Appalachian cases on the across and along nights, with held-out mountain radars (AR_K settings).
- **Prior:**
  - stage 0, the profile anomalies over all nights;
  - stage 1, the test cases;
  - stage 2, one us_east night with the best variant.

## Rules

- Jobs only through `hpc/submit.sh` from the laptop. Nothing is cloned, pulled or edited on the cluster (`hpc/README.md`).
- The snapshot carries the laptop's uncommitted toolbox changes too. `provenance.txt` of each job records them.
- Keep the raw files. Delete nothing on the cluster without the user's OK.
- No more than six of these jobs at once, and fewer if others need the node.
- Dev only: no commit, and no toolbox or manuscript change.

## Files

| File | What it does |
| --- | --- |
| `scripts/NS_A_select_nights.m` | Chooses the set from the Dark Ecology archive (seconds, no download) |
| `functions/ns_process_nights.m` | The cluster job: downloads, cleans and builds the static clutter masks of a list of nights; logs to `products/nights/run_log.csv` |
| `figures/NS_A_nights.csv` | The set: frames, purpose, what each night serves, budget |
| `figures/NS_A_nights.png` | Track against day of year, and speed against reflectivity, chosen nights marked |

**Checked:**
- `check_matlab_code` reports no issues.
- `NS_A_select_nights.m` runs.
- The night settings of `ns_process_nights` give the cleaning cache of the multi-night aspect plan.
- The rain-day settings keep all 642 daytime radar-frames.

**Not run:** the cluster jobs.
