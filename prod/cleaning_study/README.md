# Scan Cleaning Study

Purpose: decide which radar gates are birds before reconstruction, comparing the current ρHV rain mask with MistNet on real migration volumes. The study is separate from the manuscript workflows. Its result will set the `read_rd` options used by `us_night/` and later runs.

## Run

```matlab
project_root = setup_paths();
run(fullfile(project_root, "prod", "cleaning_study", "scripts", "SC_A_screen_volumes.m"))
run(fullfile(project_root, "prod", "cleaning_study", "scripts", "SC_B_select_cases.m"))
run(fullfile(project_root, "prod", "cleaning_study", "scripts", "SC_C_plot_gallery.m"))
```

| Script | Does | Output |
| --- | --- | --- |
| `SC_A_screen_volumes.m` | Draws 300 random volumes: 107 radars east of the Rockies, spring (10 Apr–31 May) and fall (25 Aug–31 Oct) nights in 2022–2025, within ±4 h of local solar midnight. Each volume is cleaned with `removeRain` set to `none`, `rho` and `mistnet`. | `data/screening.mat`, `data/screening.csv` |
| `SC_B_select_cases.m` | Selects a review set of 24 volumes: 6 where MistNet removes much more than ρHV, 6 where ρHV removes much more, 4 with heavy removal by both methods, 4 clean controls with strong echo, and 4 at random. No radar appears more than twice. | `data/cases.mat`, `data/cases_review.csv` |
| `SC_C_plot_gallery.m` | Exports one review figure per case. | `figures/gallery/*.png` |
| `SC_E_rain_methods.m` | Cleans 60 screening volumes (30 with rain, 30 others) with ρHV > 0.95, DR < −12 dB and DR < −15 dB, and compares each rain mask with MistNet. | `data/rain_methods.mat`, `data/rain_methods.csv` |
| `SC_D_filter_variants.m` | Cleans the gallery cases and the volumes behind the prod map aberrations with the current prod settings and candidate settings (MistNet + ρHV cap at 0.95 / 0.90, 2 km rain dilation, no-echo zeros). Results are in `prod/README.md`. | `figures/variants/*.png` |

## Metrics (per volume)

- **`removed_<method>`:** the share of the `none` signal (density-weighted) that the method removes.
- **`disagreement`:** the density-weighted share removed by one method but not the other.
- **`vid_<method>`:** the vertically integrated density (birds km⁻²) from the mean density profile, 0–3000 m.

All other `read_rd` options stay at the `us_night` settings. Only the rain removal changes.

## Review figure

- **Row 1:** raw reflectivity, radial velocity, ρHV and ZDR at 0.5°.
- **Row 2:** bird density after each method, and a map of the gates each method removed.
- **Row 3:** removed-gate maps at about 1.5° and 3.5°, the density profiles, and the screening metrics.

Cleaning in the figures runs at full gate resolution. The screening metrics use the 2 km decimation applied in the reconstruction.

## Review

Fill in `data/cases_review.csv`:

- **`verdict_rho` and `verdict_mistnet`:** one of `ok`, `contaminated` (non-bird echo kept) or `over_removed` (birds removed).
- **`contamination`:** the type of echo left, for example rain, insects, clutter, sea clutter, wind farm or anomalous propagation (AP).
- **`notes`:** free text, for example which sector is affected.

Neither method separates birds from insects, and neither removes ground or sea clutter. The clean-control group is there to show whether these need their own filter.

## Expected limitations

- MistNet classifies only precipitation versus biology. It was trained on legacy single-polarization data, and the wrapper interpolates its Cartesian predictions back to the polar gates.
- The ρHV threshold (0.95) misses low-ρHV rain: storm edges, the melting layer and hail.
- Raw volumes are cached outside OneDrive (see `local_settings.json` and `nexrad_raw_dir`).

## First results (22 October 2026)

- Screening: 300 volumes drawn, 287 with files, all 287 cleaned after fixing a `read_rd` bug. When ρHV and reflectivity sweep lists differed (for example an extra 0.48° ρHV sweep at KCLE), the shape correction paired sweeps by position. It now pairs them by elevation. Before the fix, 18 volumes failed and others could silently use the wrong ρHV sweep.
- Median removed signal: ρHV 2.8%, MistNet 0.5%. In the top 10% of volumes, both remove more than 50%.
- Disagreement is strongly asymmetric. Where MistNet removes much more, it is usually widespread precipitation that the 0.95 ρHV threshold misses (stratiform rain at 0.90–0.97 and the melting-layer ring). Where ρHV removes more, the volumes are almost empty (VID < 15 birds km⁻²).
- MistNet can over-remove biology next to rain (C01, possibly C14).
- Draft per-case readings are in `data/cases_draft_notes.csv` (by Claude, not expert-validated). Fill the review columns in `data/cases_review.csv`.
- Candidate next methods to tune: a lower ρHV threshold (0.90), and MistNet restricted to gates with ρHV above about 0.85, so that clearly biological gates next to rain are kept.

## Depolarization ratio against MistNet (23 October 2026)

`read_rd` has a `removeRain = "dr"` option: a gate is rain when its depolarization ratio, estimated from ZDR and ρHV for simultaneous H/V transmission (Kilambi, Fabry and Meunier 2018), DR = 10 log10[(1 + Z − 2ρ√Z) / (1 + Z + 2ρ√Z)] with Z the linear ZDR, is below `rainDrThreshold` (default −12 dB). Rain (ρHV ≈ 1, ZDR ≈ 0) has low DR; birds, insects and clutter have high DR. Unlike a ρHV threshold, DR also uses ZDR, so biology with high ρHV but high ZDR is kept. The default stays MistNet.

`SC_E_rain_methods` (rough estimate, 60 volumes, same `read_rd` options, no rain dilation), with MistNet as the reference. "Extra" is signal removed by the rule but kept by MistNet; "missed" is signal MistNet removes but the rule keeps (both weighted by density, means over volumes):

| Rule | Rainy volumes (n = 30): extra / missed | Other volumes (n = 30): extra / missed |
| --- | ---: | ---: |
| ρHV > 0.95 | 2.6% / 22.2% | 2.6% / 0.5% |
| DR < −12 dB | 2.5% / 11.7% | 1.8% / 0.5% |
| DR < −15 dB | 1.5% / 21.5% | 1.1% / 0.5% |

DR at −12 dB is closer to MistNet: it misses about half as much rain as the ρHV rule, with the same extra removal in rainy volumes and less in rain-free ones. DR at −15 dB removes the least biology but misses as much rain as ρHV. MistNet is the reference here, not the truth. The DR rule costs about 0.2 s per volume on one core (ZDR is already read), against about 5.6 s for MistNet.
