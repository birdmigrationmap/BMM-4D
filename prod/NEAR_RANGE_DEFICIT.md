# Near-range echo deficit

Working note, started 24 October 2026. Data: sandbox case (Midwest, 28 radars, 6 frames 00:30–03:00 UTC on 22 October 2026, VCP 35 on 25 radars). Scripts: `sandbox/SB_near_range_raw.m` (raw data), `sandbox/SB_near_range.m` (inversion), and `functions/read_level2_moments.m` (reads REF and CFP from the Level II files). Figures: `sandbox/figures/near_range_*.png`.

## The problem, defined without a model

**Near-range echo deficit:** within one volume, the bird reflectivity a radar reports for a given height above ground falls the closer the gate is to the radar.

For radar *k*, height layer *h* and range *r*, let D_k(r, h) = η_k(r, h) / η_k(30–60 km, h), where η is the mean density over all gates whose beam centre lies in layer *h*, and no-echo gates count as 0. If birds are roughly uniform over about 50 km around the radar, D should be about 1 at every range. Measured on the radar's own output (REF, before any of our processing), the median over radars of D at 300–500 m above ground is:

| Range | 3 km | 5 km | 7 km | 9 km | 11 km | 13 km | 15 km | 17 km | 19 km | 21 km | 24 km | 26 km | ≥ 29 km |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| D | 0.03 | 0.05 | 0.10 | 0.17 | 0.29 | 0.40 | 0.52 | 0.65 | 0.72 | 0.87 | 0.90 | 0.94 | ~1 |

The deficit reaches farther for higher layers: at 800–1200 m, D only reaches 1 at about 35 km. It is present at every radar. Per radar, D at 10 km and 300–500 m ranges from about 0.03 to 0.3.

The same deficit shows when observations are compared with an inversion that never saw the near rows (`near_range_observations.png`): observed / predicted at 4 km is 0.10 at 0.5° and 0.01 at 3.1°.

**How it relates to what we saw in the maps**
- The **radar-distance artefact** in the reconstructed maps (a dark disc at each radar, then a ring) is a *consequence* of this deficit. The inversion reads the suppressed near-range values as true densities. Data cuts and weights only hide them.
- The **far-range column deficit** (columns 100–150 km from a radar at 0.85 of their regional mean) is a *separate* artefact: vertical sampling. Beyond about 60 km the beams pass above the dense low layer, and the flatter prior profile fills the low layers there (`near_range_vertical.png`).

## What the raw data show

1. **The deficit is in the radar's output, not created by our processing.** Our cleaning makes it worse: at 15 km and 300–500 m, D is 0.52 in REF and 0.22 in our cleaned volumes; at 11 km, 0.29 against 0.10.
2. **The clutter filter removes power almost everywhere close to the radar.** Within about 20–30 km, 75–88% of gates have power removed (CFP code ≥ 8), at every elevation.
   - The amount removed falls with range. At 1.3° and above it is 20–25 dB at 3 km, about 8 dB at 10 km and 2–3 dB at 25 km. At 0.5° it is 37 dB at 3 km and still 10 dB at 40–60 km (ground clutter).
   - Each elevation cut has a sharp outer limit to the filtered region, between 15 and over 60 km.
3. **Adding the removed power back (REF + CFP) recovers only part of it.** At 3–5 km, D goes from 0.03–0.05 to 0.12–0.15. From about 9 km out, the unfiltered D is the same as the output (0.24 against 0.17 at 9 km; 0.55 against 0.52 at 15 km). So the clutter filter's notch explains the first few kilometres, not the deficit at 9–30 km.
4. **The bird flight direction matters only within about 10 km.** Across the flight, where the radial velocity is near 0, D is about half of D along the flight within 10 km, as a velocity notch filter would cause. Beyond 10 km there is no difference.
5. **The near-range reflectivity is not spiky.** At 5–10 km and 300–500 m, 92% of gates have echo, and the values form a broad distribution with a mode near −10 dBZ. At 30–60 km the mode is +14 dBZ.
   - A gate at 5–10 km holds about 0.001–0.005 km³, i.e. about 0.05–0.25 birds at 50 birds km⁻³. Most gates should be empty, and a few should hold one bird each at about +20–25 dBZ.
   - Those spikes are almost absent: 1.4% of gates are above 20 dBZ.
6. **The point clutter filter (CFP code 1) flags 17% of gates at 3–7 km.** That share falls to 10% at 21 km, 5% at 36 km and 2% at 60 km, the same range as the deficit. Its decline with range matches the expected number of birds per gate rising above 1 (about 1 at 20 km, 2–3 at 30 km). Flagged gates keep an ordinary, low reflectivity.

## Explanations, ranked by the current evidence

| # | Explanation | Status |
| --- | --- | --- |
| H1 | **Discrete targets removed as point clutter.** Close to the radar, a gate holds less than one bird, so each bird is an isolated strong point. The radar's point clutter filter (and our outlier, continuity and max-density filters) replaces these spikes by the background. Linear averaging over the gates then misses most of the birds. The deficit should end where gates hold several birds, so it should depend on density and gate volume (range² × beam width²), not on elevation as such. | **Leading.** It fits the range of the deficit, the missing spikes and the code-1 share. Not proven: CFP gives no removed power for code-1 gates. |
| H2 | **Velocity notch of the clutter filter (GMAP/CMD).** It removes bird power near zero radial velocity, more where the clutter is stronger (wider notch). | **Real within about 8–10 km** (unfiltered recovery, across-flight deficit). Too small beyond that. |
| H3 | **Our cleaning** (MAD outliers, 5 × 5 continuity filter, 3000 birds km⁻³ cap, 25 km² patch filter). These are spike filters too. | **Adds to the deficit** (D halved again at 10–20 km). Each step not yet isolated. |
| H4 | Fewer birds near radar sites (siting near towns and airports). | Unlikely as the main cause: factor 3–30 at all 28 radars and in every layer up to 1.2 km. |
| H5 | Geometry: beam broadening makes the far reference average over a thicker layer, including the dense low layer. | Can bias the far reference by a factor of about 2 in the lowest layers, not by 10–30. To check with the forward model on a synthetic profile. |

## Tests to decide between them

1. **Volumes from before the point clutter filter existed.** CFP and point clutter filtering came with RDA Build 19 (2019–2020). The same check on volumes from 2018 (we have some, e.g. KBYX 2018) would show whether the deficit at 9–30 km was smaller then. If so, H1 is the radar's processing.
2. **Expected spikes.** Simulate Poisson-distributed birds through the gate volumes (range, beam width, gate length) at the density of the 30–60 km ring, and compare the expected near-range reflectivity distribution with the observed one, with and without removing isolated gates.
3. **Density dependence.** If H1 holds, the deficit should reach less far on dense frames and radars (more birds per gate). Compare D across radars and frames against the 30–60 km density.
4. **Our cleaning, step by step.** Run read_rd with each filter switched off in turn (outliers, continuity, cap, patch filter, smoothing, shape correction) on the same volumes.
5. **Unfiltered reflectivity.** The Level II files have no unfiltered REF. If a site provides it (Level I, or a research radar), compare it directly.

## Sources
- CFP codes (0 filter not applied, 1 point clutter filter applied, 2 dual-pol-only filter applied, 3–7 reserved, ≥ 8 power removed = code − 8 dB): open-source NEXRAD decoder [danielway/nexrad](https://github.com/danielway/nexrad) (`nexrad-model/src/data/moment.rs`), which follows the RDA/RPG ICD ([2620002Y](https://www.roc.noaa.gov/public-documents/icds/2620002Y.pdf); the NOAA server refused the connection, so the ICD table itself is not yet checked). CFP was added with RDA Build 19 ([pyart issue 887](https://github.com/ARM-DOE/pyart/issues/887)).

## Follow-up (24 October 2026)
`prod/gate_density/` (README) explains the radar's clutter processing and CFP in more detail. It replaces outlier removal and median smoothing by a homogeneous-density estimator, tests it on simulated sweeps, and applies it to the sandbox sweeps.
- Our moving-median smoothing alone accounts for a factor of 2–4 at 10–20 km.
- With the estimator, the deficit at 300–500 m is gone beyond 13 km and reduced within 10 km, where the radar filters almost every gate.
