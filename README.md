# bmm: 3D bird-migration reconstruction from weather radars

Code for reconstructing animated 3D bird-density fields from a network of NEXRAD radars: radar processing, forward model, density prior, advection-regularized inversion, and the operational pipeline built on them.

This repository was split from the paper repository `BMM-within` on 24 October 2026. The paper keeps its own frozen copy of the code; tag `paper-v1` marks the code state behind the submitted manuscript figures (paper commit `e52a81b`). Changes made here reach the paper only when they are ported on purpose and the affected workflows rerun.

## Layout

| Path | Content |
| --- | --- |
| `functions/` | Shared code: radar reading and cleaning (`radar/read_rd.m`), beam kernels (`F_elev/`), forward operator (`construct_F.m`), grid, priors, advection, synthetic-case tools, plotting helpers. |
| `prod/` | Operational pipeline (P1–P5), lean solver, sandbox, and studies: gate density near the radar, rain cleaning, the two-level model. See `prod/README.md`. |
| `preprocessing/beam_blockage/` | Builds the terrain beam-blockage maps used by `read_rd`. |
| `tests/` | Unit tests of the shared code. |

## Setup

```matlab
project_root = setup_paths();
```

`setup_paths` adds `functions/` and two external libraries from fixed local paths: [wsrlib](https://github.com/darkecology/wsrlib) (required) and FastGaussianSimulation (optional, synthetic cases). `functions/regtools/regu/` (Hansen's Regularization Tools, needed for `get_l`) is not tracked; copy it locally.

## Data

Nothing under `data/` is tracked. The code reads it by relative path from the project root:

| Path | Content | Source |
| --- | --- | --- |
| `data/reference/radar.csv` | Radar sites | Copied from the paper repository |
| `data/reference/DEM/`, `onwater/`, `scan_pattern_new.mat` | Terrain, land–water mask, scan patterns | Copied from the paper repository |
| `data/reference/beam_blockage/` | Terrain-horizon maps | `preprocessing/beam_blockage/beam_blockage_script.m` |
| `data/derived/F_elev/150km/` | Beam elevation kernels | Copied; `functions/F_elev/generate_F_elev.m` rebuilds them |
| `data/raw_scans/nexrad/` | Older raw volumes | Symbolic link to the paper repository's folder (read only) |

`local_settings.json` (not tracked) sets machine-specific folders: `nexrad_raw_dir` for the raw NEXRAD cache and `local_data_dir` for large generated products (radar caches, cases, inversions).
