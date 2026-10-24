# BMM-4D: bird-migration reconstruction in space and time from weather radars

Code for reconstructing animated 3D bird-density fields from a network of NEXRAD radars: radar processing, forward model, density prior, advection-regularized inversion, and the operational pipeline built on them.

This repository was split from the paper repository `BMM-within` on 24 October 2026. The paper keeps its own frozen copy of the code; tag `paper-v1` marks the code state behind the submitted manuscript figures (paper commit `e52a81b`). Changes made here reach the paper only when they are ported on purpose and the affected workflows rerun.

## Layout

| Path | Content |
| --- | --- |
| `toolbox/` | The package, in the `bmm4d` namespace (`toolbox/+bmm4d/`): radar reading and cleaning (`bmm4d.radar`), beam kernels and forward operator (`bmm4d.forward`), grid and masks (`bmm4d.grid`), priors (`bmm4d.prior`), advection (`bmm4d.advect`), solvers (`bmm4d.solve`), plotting (`bmm4d.plot`), helpers (`bmm4d.util`). Calls are written in full, e.g. `bmm4d.radar.read_rd(...)`; a script can `import bmm4d.radar.*` to use short names. |
| `pipeline/` | Operational pipeline P1–P5 with its case settings, caches and logging (`pipeline/functions/`). See `pipeline/README.md`. |
| `synthetic/` | Synthetic truths, advection and observations for tests and studies; not part of the package. |
| `studies/` | Research beside the package: gate density near the radar, cleaning study, sandbox, two-level model, beam-blockage maps, legacy US night, notes. |
| `tests/` | Unit tests of the shared code, and `paperRegressionTest` (tag `Regression`): today's code against outputs of tag `paper-v1` on three Key West volumes. Its reference is local (`data/regression/`), built by `tests/regression/make_regression_reference.m`. |

## Setup

Copy `local_settings.example.json` to `local_settings.json` (not tracked) and set the folders of this machine, then:

```matlab
project_root = setup_paths();
```

`setup_paths` adds `toolbox/` and the external code, whose folders come from `local_settings.json` (absolute, or relative to the repository root; defaults under `external/`, which is not tracked):

| Setting | Code | Needed for |
| --- | --- | --- |
| `wsrlib_dir` | [wsrlib](https://github.com/darkecology/wsrlib): NEXRAD reading (rsl2mat), VVP, MatConvNet for MistNet. On Apple Silicon, MatConvNet must be built without finite-math assumptions (branch `fix/maca64-finite-math` of [Rafnuss/wsrlib](https://github.com/Rafnuss/wsrlib)); `mistnet_polar_radar` stops otherwise. | Required; `setup_paths` stops with an error without it |
| `regtools_dir` | Hansen's Regularization Tools (`cgsvd`, `tikhonov`, L-curve corner) | `compute_vp` and `compute_scan` |
| `fast_gaussian_simulation_dir` | FastGaussianSimulation | Synthetic Gaussian truth fields (`synthetic/`) |

MathWorks toolboxes used: Optimization, Mapping, Image Processing, Statistics and Machine Learning, Curve Fitting, Parallel Computing.

## Data

Code reads data through `bmm4d.util.data_path(...)`, from the data folder: `data_dir` in `local_settings.json`, or `data/` in the repository root. Nothing in it is tracked.

| Path in the data folder | Content | Source |
| --- | --- | --- |
| `reference/radar.csv` | Radar sites | Copied from the paper repository |
| `reference/DEM/`, `onwater/`, `scan_pattern_new.mat` | Terrain, land–water mask, scan patterns | Copied from the paper repository |
| `reference/beam_blockage/` | Terrain-horizon maps | `studies/beam_blockage/beam_blockage_script.m` |
| `derived/F_elev/150km/` | Beam elevation kernels | Copied; `pipeline/scripts/generate_F_elev.m` rebuilds them |
| `raw_scans/nexrad/` | Older raw volumes | Symbolic link to the paper repository's folder (read only) |
| `regression/` | Reference of the paper-v1 regression test | `tests/regression/make_regression_reference.m` |

Two more folders are set in `local_settings.json`: `nexrad_raw_dir`, the cache of downloaded raw volumes (`bmm4d.radar.nexrad_raw_dir()`), and `local_data_dir`, the large generated products such as radar caches, cases and inversions (`bmm4d.util.local_data_dir()`).
