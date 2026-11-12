# Robust aspect correction: test plan

Goal: one aspect correction that holds for every night and every radar, without fitting anything to a single volume. This plan orders the tests from cheap checks on data we already have to fits over many nights. Each test states the question, the method, the cost and what decides.

This is a proposal for discussion. Nothing here has been run. The multi-night fit of `PLAN_multi_night.md` is one stage of it (C1).

## Data available (9 Nov 2026)

The shared night set of `studies/nights/PLAN.md` (list in `studies/nights/figures/NS_A_nights.csv`) is on the cluster:
- **18 full nights, cleaned:** 17 Dark Ecology nights (2021–2025, both seasons) plus 21–22 Oct 2026. They use 107 us_east radars every 30 min, production cleaning with `"shape", false`, in cache `493fd773`, and their raw files are kept. Track axes span about 95°, and several nights have tracks that vary across the domain (axis coherence 0.08–0.35).
- **3 rain days, downloaded but not cleaned:** 11 Apr 2024, 5 Apr 2025 and 28 Oct 2021, hourly from 15 to 20 UTC.
- **Not yet available:** winds. ERA5 on pressure levels is planned there, about 30 MB a night.

So stages B and C need no new radar downloads, only binning and fits. Runs go through `hpc/submit.sh` and need the user's go-ahead.

## What we know and what blocks us

From the 21–22 Oct 2026 night (`README.md`):
- The overlap gives a clean shape: c2 (side-on) = 0.20 and c1 (front/back) = −0.21 pooled, with c2 growing with height (0.092 per km).
- One night cannot separate bird aspect from fixed radar patterns. With a fixed azimuthal pattern per radar in the fit, c2 falls from 0.092 to 0.017 per km and c1 from −0.21 to −0.06.
- The axis is uncertain: the side-on pattern does not turn with the track during the night, and on the overlap it lies between the track and a fixed heading.
- Below 0.8 km there is no gradient-free measurement.
- The amplitude varies with region (0.12–0.28) and ground speed (0.06 below 6 m/s, 0.24 at 10–14 m/s).

Three things must be identified before a general function makes sense:
1. **Radar patterns.** What part of the overlap signal belongs to the radars (blockage, calibration by sector) and stays in their geographic frame.
2. **Axis.** Whether the pattern follows the heading (body axis) or the track.
3. **Covariates.** What the amplitude depends on (height, airspeed, insects, transmit frequency, season, region), so that one function transfers.

Principles:
- Separate effects by their frame. Aspect turns with the birds; radar patterns stay with the radar. Any data where the two frames move apart (other tracks, other headings, precipitation) separates them.
- Treat the amplitude as empirical. Single-bird values (11–20 dB side/head with wings closed at X band) do not transfer to populations in flight at S band, where passerines are near resonance and wing flapping alone moves the cross-section by at least 10 dB (Gong et al. 2022; Melnikov et al. 2012). No literature prior on c2.
- Score out of sample: other frames, other nights, other radars.

## Stage A: existing night, no new radar data (hours each)

### A1. Can this night separate aspect from radar patterns at all?

- **Question.** Is the collapse of c2 in the joint fit (`as_overlap_fit_radar`) real, or a design that cannot tell the two apart?
- **Method.** Simulate log ratios on the real overlap design (same cells, pairs, frames, fixed effects, noise from the residuals):
  - (a) pure aspect with the shared function;
  - (b) pure radar patterns with the per-radar harmonics of the joint fit;
  - (c) both.
  Refit the joint model on each. If (a) also collapses, the night cannot decide and the joint fit carries no information.
- **Also.** Compute the variance inflation of c2 under the joint model as a function of how much the track varies per radar. This gives the track spread new nights must bring (input to B2 and C1).
- **Cost.** Minutes; uses the saved design.

### A2. Heading or track: add winds

- **Why.** Headings differ from tracks by 34° on average at inland radars and 42° at coastal ones (Horton et al. 2016, ρHV headings and NARR winds). If the axis is the heading, a track-based fit sees the side-on term reduced by the mean of cos 2Δ over the heading–track angles Δ, and a constant offset shows up as a sin 2ψ term. A KNMI report found the axis of the reflectivity pattern at the heading (210°) rather than the track (230°) at one radar (Dokter et al. 2009; not checked in full text).
- **Method.**
  - Get model winds for the night from an open copy: ERA5 on AWS (`nsf-ncar-era5`, no account) or HRRR on AWS (`noaa-hrrr-bdp-pds`, from 2014). The night set plans ERA5 for all 18 nights; reuse it once it exists. Interpolate to each radar's 100 m layers and frames.
  - Heading = direction of (ground velocity − wind); also airspeed.
  - Refit the overlap with ψ from the heading; compare c2, the disagreement and the residual sin 2ψ term with the track-based fit.
  - Redo the track-turning test (`as_time_fit`) with the heading. If the pattern turns with the heading during the night but not with the track, the axis is the heading.
- **Decides.** Which axis `as_aspect_template` uses. If the heading wins, production needs winds too (or the ρHV heading, which needs the raw volume).
- **Cost.** A small wind download for one night; refits in minutes.

### A3. What c2 depends on

- **Question.** Is the growth of c2 with height a property of the birds, or of insects at low levels, slow targets, viewing geometry or the radar?
- **Candidates and why.**
  - Airspeed (from A2) instead of ground speed: at low ground speed the heading–track angle grows, so a track-based c2 shrinks.
  - An insect indicator (airspeed below about 5 m/s, ZDR, ρHV). Insects are Rayleigh scatterers at S band with their own side-on maximum on their body axis (Riley 1985; Stepanian et al. 2016), and they often head downwind (Chapman et al. 2010), so their axis differs from the birds'.
  - Drift behaviour with height: compensation at low levels and more drift aloft (Liechti 1993, via Bruderer et al. 2025; secondary). More uniform headings aloft would make c2 grow with height.
  - Elevation angle and range of each radar, and the beam-width ratio of the pair. A KNMI report saw a stronger contrast at 0.5° than at 2.1° (Dokter et al. 2009; not checked).
  - Transmit frequency of each radar: within S band, bird reflectivity can change by more than 10 dB with frequency (Melnikov et al. 2012), so the aspect contrast may differ between radars for a physical reason.
- **Method.** One regression of the overlap harmonics on these covariates together, by stratum (as `as_fit_strata`), checked on the equal-range, matched-beam subset. Keep the sin ψ and sin 2ψ terms per radar pair as diagnostics: an axis offset or a transmit-phase effect would show there.
- **Decides.** The covariates kept in T, and whether gates likely dominated by insects get T = 1 or their own pattern.
- **Cost.** Minutes to an hour (airspeed needs A2).

### A4. How precise must T be?

- **Question.** How much do reconstructed columns and maps change between plausible functions?
- **Method.** Reconstruct the Midwest night with T scaled by 0, 0.5, 1 and 1.5, and with the axis from the track or the heading. Score held-out radars and the first harmonic of each radar's own data, as in `studies/artifacts`.
- **Decides.** The tolerance for the later fits. If columns change little between 0.5 and 1.5, a rough climatological function is enough and later stages can be lighter.
- **Cost.** A few reconstructions (the longest step of Stage A).

## Stage B: small data sets, already downloaded (hours each)

### B1. Radar patterns without birds

- **Question.** What are the fixed azimuthal patterns of each radar, measured without birds?
- **Method.**
  - Precipitation: the three rain days of the night set. Read the raw reflectivity (the cleaning would remove the rain), keep stratiform areas and low sweeps below the melting layer. Run the same overlap regression on rain reflectivity, as in neighbour-radar comparisons (Gourley et al. 2003). Rain has no aspect, so the per-radar 1st and 2nd harmonics in the geographic frame are radar patterns.
  - Ground clutter: relative calibration from clutter (Silberstein et al. 2008; Hunzinger et al. 2020), binned by azimuth, as a cheaper check of calibration by sector.
  - Spaceborne radar (GPM) gives an absolute reference of about 1 dB (Schwaller & Morris 2011; Louf et al. 2019), probably too coarse for patterns of a few tens of percent, but a check for the worst radars.
  - Compare with the per-radar patterns of the joint bird fit. If they agree, subtract them from the bird overlap (or use them as a prior) and refit c2 and c1 on the bird night.
- **Decides.** Whether the one-night amplitude survives once radar patterns are fixed independently. This attacks the main confound without needing many bird nights.
- **Caveats.** Rain fills the beam differently from birds and has its own vertical structure, so patterns from partial blockage may depend on the target's profile. Bird cross-sections also depend on frequency, which rain does not share. Check on radars with known blockage.
- **Cost.** 3 days × 6 frames × 107 radars of raw reflectivity, already downloaded; reading and binning take hours on the cluster. On these days only 28–32% of the radars have rain over half the day or more, so some radar pairs will have no rain overlap.

### B2. Two or three nights with turned tracks

- **Question.** Does the shared function transfer to nights whose track differs by 45–90°, and does the joint fit keep the amplitude there?
- **Method.** Take cleaned nights from the night set whose track is turned relative to 21–22 Oct (axis about 190°, or 10° modulo 180°). The side-on term repeats every 180°, so a spring night with the opposite track does not help it; turned tracks do. Candidates:
  - 2021-08-27, track 258° (axis turned about 68°);
  - the nights with tracks that vary most across the domain: 2024-05-10, 2023-05-22 and 2023-04-10 (axis coherence 0.08–0.26);
  - one spring night (tracks 345–33°), which flips the front/back term and so tests c1.

  Bin the overlap with `as_night_bins` and fit:
  - with the one-night function fixed (transfer score);
  - free (per-night fit);
  - jointly with radar patterns (from B1, or free) over all nights together.
- **Decides.** Whether to scale up to C1, and which nights matter most. This is the pilot of `PLAN_multi_night.md`.
- **Cost.** The nights are already cleaned; binning takes minutes per night (19 s per frame measured).

## Stage C: fit through many scans

### C1. Multi-night overlap fit

- **Question.** One climatological function, or a fit per night?
- **Method.** As in `PLAN_multi_night.md`, with three additions from Stages A and B:
  - the axis chosen in A2 (heading with winds, if it wins);
  - radar patterns fixed from B1, or fitted jointly now that tracks vary;
  - the covariates kept in A3, with a random effect per night for what they do not explain.
  Leave one night out. The climatological function is enough when it gets most of the per-night fit's reduction of the disagreement, within the tolerance set in A4.
- **Cost.** The 18 nights are cleaned. Binning them all takes about 2–3 h, and the fit about 1 h with about 10 GB of memory (see `PLAN_multi_night.md`).

### C2. Below 0.8 km from single radars

- **Question.** What is T where most birds fly and no two radars see the same cell?
- **Method.** Stack single-radar harmonics over many nights in the flight frame (as `as_time_fit`, but across nights). Real gradients and radar patterns have no preferred direction relative to the flight once tracks vary between nights, so they average out; aspect does not. Anchor to the overlap where both exist (0.8–1.5 km) and use the single-radar values below.
- **Decides.** Replaces the 1 km ramp of c1, which is a choice now.

## Stage D: operational checks

- **Fallback.** For nights or radars with few overlap pairs (isolated radars, the coast, the west), use the climatological function from C1.
- **Reconstruction.** Held-out radars over several nights and regions, with the wedge and first-harmonic scores of `studies/artifacts`. A correction is accepted when it improves the median held-out error without flattening each radar's own data.
- **Production.** Only after this: a proposal for `toolbox/` and Appendix B, for the user's approval.

## Order and decision points

| Step | Needs | Cost | Go on if |
| --- | --- | --- | --- |
| A1 identifiability | saved design | minutes | always; sets the track spread needed |
| A2 winds and heading | open wind data | hours | always |
| A3 covariates | A2 | hours | always |
| A4 tolerance | reconstructions | about a day | always; sets how precise later fits must be |
| B1 radar patterns | 3 rain days (downloaded) | hours | A1 shows the night cannot decide alone |
| B2 turned-track nights | 2–3 cleaned nights | hours | B1 done or skipped |
| C1 multi-night fit | 18 cleaned nights, `PLAN_multi_night.md` | a few hours | B2 shows the function moves between nights by more than the A4 tolerance |
| C2 below 0.8 km | C1 data | hours | C1 done |
| D operational | all above | days | user's approval |

## Suggested changes to the dev code

None is applied.

- `as_aspect_template`: an option for the axis from the heading (ground velocity − wind) instead of the track; cap c2 above about 3 km (linear growth overshoots: 0.30 against 0.24 measured); the covariates kept in A3.
- `as_overlap_fit_radar`: accept fixed per-radar patterns (from B1) as offsets or a prior.
- `as_overlap_design`: carry airspeed, heading, transmit frequency and an insect indicator per cell, so A2 and A3 need no new binning.

## Literature notes for README.md

A literature check on 5 Nov 2026 suggests these corrections to "What the literature says". Some items were checked only in abstracts or secondary sources, as marked.

- **Edwards & Houghton 1959.** The 14–90× and 1.1–1.4× values come from a table in Konrad, Hicks & Dobson 1968 (Science 159:274): suspended birds, wings closed, X band, 18° elevation (secondary source).
- **Vaughn 1985.** The "10–20 dB wingbeat" figure was not found. Gong et al. 2022 (IEEE Access 10:85137) report at least 10 dB fluctuation from flapping at S, X and Ku band.
- **Mirkovic et al. 2016.** A measured bat with spread wings at X band: HH head 2.78, side 1.49, tail 1.88 cm² (side/head 0.54); a wingless ellipsoid model gives about 4.2.
- **Urmy & Warren 2017** (MEE 8:860): a commercial marine radar; X band not confirmed.
- **Kranstauber et al. 2020.** No measured factor of 2: the paper notes lobes perpendicular to the flight and expects the error not to exceed a factor of two, citing Dokter et al. 2009. It applies no correction.
- **Melnikov et al. 2015** (Melnikov, Istok & Westbrook, JTECH 32:659) and **Melnikov & Bridge 2024** are about insects, not birds. Stepanian et al. 2016 also show insect ρHV/ΦDP patterns flipping between radars because of the transmit differential phase.
- **Bruderer & Joss 1969:** no aspect numbers found.
- **To add:**
  - Hoekstra et al. 2025 (J Environ Manage 395:127753): remove a 180°-period sine fitted to each scan's dBZ. This is prior art for the per-scan approach, with the same weakness (it fits real gradients and radar patterns).
  - Dokter, Liechti & Holleman 2009 (KNMI WR 2009-06): azimuthal correction not implemented; pattern axis at the heading; stronger at 0.5° than 2.1° (not checked in full text).
  - Horton et al. 2016 (Sci Rep 6:21249): heading − track 34° inland, 42° coastal.
  - Melnikov, Lee & Langlieb 2012 (IEEE GRSL 9:413): more than 10 dB change of bird Z with frequency within S band.
  - Schmid et al. 2019 (Ecography 42:931): a vertical view largely removes the aspect dependence.
  - Riley 1985 (Proc IEEE 73:228) and Chapman et al. 2010 (Science 327:682): insect scattering and orientation.
  - Calibration: Gourley et al. 2003; Silberstein et al. 2008 (JTECH 25:2037); Hunzinger et al. 2020 (AMT 13:3147); Schwaller & Morris 2011 (JTECH 28:301); Louf et al. 2019 (JTECH 36:17).
