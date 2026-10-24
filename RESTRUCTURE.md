# BMM-4D restructure: scope and target layout (draft for review)

## 1. Scope

**BMM-4D is a development package for reconstructing 4D bird density from real weather-radar networks.** It covers the whole applied path: get volumes, clean them, build the observation system, run the inversion, make maps. It is still a research tool, and cleaning, observation model, prior and solver will keep changing. The structure must make trying a new method easy, not freeze one.

**What came over from the paper is already the right amount.**
- **Came over:** the shared code (`functions/`) and its unit tests.
- **Stayed in the paper repository:**
  - the manuscript and figures;
  - the eleven numbered workflows and their scripts;
  - the synthetic case files and the paper data products.

Two things still need sorting:
- **Synthetic tools (19 functions):** move to a separate `synthetic/` folder beside the package. Tests and studies use them, but the package does not depend on them.
- **Paper-figure helpers (about 10 functions):** drop them. The paper keeps its own copies.

**Principles**
1. **Keep alternatives, remove duplicates.**
   - A different method (ρHV vs MistNet rain removal, IDW vs Gaussian prior, with or without non-negativity) stays, selected by name in the settings, with a documented default.
   - A second implementation of the same method goes. Example: the full-grid solver, which the lean solver reproduces to 2e-5 with every cell active.
2. **Steps with a common interface.** Each stage of the path (cleaning filters, prior, solver terms) is a function with the same inputs and outputs as its alternatives. A new method is then one new file plus one settings value.
3. **Settings in one struct, with defaults.** Today's `prod_case` becomes `bmm4d.config`. A case is a few lines that override defaults.
4. **Studies graduate into the package as methods.**
   - Gate density, the cleaning study and the sandbox options stay outside until they work.
   - Each then enters as a new selectable method, and becomes the default only after it is checked on the sandbox.

## 2. Stages and their methods

| Stage | Methods kept (default first) | Restructure |
| --- | --- | --- |
| Get volumes | list, local cache, HTTPS download | – |
| Clean | rain: MistNet, ρHV, DR, MistNet+ρHV, none; beam blockage; static clutter; no-echo gates; outliers; patch filter; thinning | **Split `read_rd` (790 lines) into one function per filter**, applied in the order given in settings. Each filter takes and returns the same sweep struct and a mask of what it removed. Each filter is then testable and replaceable on its own, and gate density fits in as one more step. |
| Density from gates | outlier + median smoothing (default); gate-density estimator later | Same filter interface |
| Grid and domain | radar-range domain, geographic domain, aggregation, water mask | – |
| Forward operator | `construct_F` (with aggregation); kernel loading and generation | `construct_F_compact` moves to `synthetic/` if only synthetic cases use it |
| Velocity and advection | radar velocity profiles → advection steps; relaxed links | – |
| Density prior | Gaussian weights (default), IDW; land and water profiles; temporal smoothness | `build_density_prior_fast` already covers both; the older `build_radar_density_prior` is a duplicate |
| Vertical profile | `compute_vp` (Tikhonov + L-curve) | Kept as a radar-level product |
| Inversion | lean solver: data, prior, 3D smoothness, advection; options PCG or projected FISTA (non-negative), row weights, correlated prior, shape term | Remove duplicate implementations (`compute_scans_advection(_active)`, `prepare_scans_advection(_active)`, `solve_advection_fast`) once the lean solver is shown to match them. Group solver options in one sub-struct of the settings. `compute_scan` (static, one radar) stays only if it is not the lean solver with one time step. |
| Products | maps, video, column density, takeoff and landing balance, radar-distance diagnostic | – |

## 3. Layout

```
BMM-4D/
├── README.md  LICENSE  CITATION.cff  AGENTS.md  buildfile.m  setup_paths.m
├── toolbox/                         THE PACKAGE (resources/mpackage.json)
│   ├── +bmm4d/                      config, get_volumes, clean_volumes, build_observations, invert, plot_map, export_video
│   │   ├── +radar/                  reading, download, cache
│   │   │   └── +filter/             one function per cleaning step (rain_mistnet, rain_rho, blockage, clutter, outliers, patches, smooth_median, …)
│   │   ├── +forward/  +grid/  +prior/  +advect/  +solve/  +plot/  +util/
│   │   └── +internal/
│   ├── data/                        small reference files
│   └── examples/                    keywest.m, midwest.m, us_east.m
├── synthetic/                       truth fields, advection, observations through bmm4d.forward
├── studies/                         gate_density, cleaning_study, sandbox, two_level, beam_blockage, notes
├── tests/                           unit, synthetic recovery, paper-v1 regression
└── data/                            local, not tracked
```

The five stage functions are the short way through. Scripts can also call any module directly to inspect intermediate objects, as P1–P5 do today.

## 4. Changes beyond moving files

1. **Data paths.**
   - Problem: code reads `data/...` relative to the current folder, which works only because `setup_paths` changes into the project root.
   - Fix: add one `bmm4d.util.data_dir()`. Small files come from `toolbox/data`; large layers come from a folder set in `local_settings.json`.
2. **External code.**
   - wsrlib: its folder is set in `local_settings.json` and checked with a clear error.
   - regtools: only `get_l` is used (4 functions). Replace it with our own difference operator.
   - FastGaussianSimulation: used only by `synthetic/`.
3. **Required MathWorks toolboxes** go in the README and package metadata.
4. **License.**

## 5. Order of work

Each step is one commit, with tests green after it.

1. **paper-v1 regression test.** Store outputs from today's code on a small Key West case (cleaned sweeps, forward rows, prior, solution). It must keep passing through every pure restructuring step, and it reports any change once methods change.
2. **Move files** with `git mv`: package, `synthetic/`, `studies/`; drop the paper-figure helpers. Record old → new in `MIGRATION.md`.
3. **Remove duplicate implementations.** Do it after an equality check against the kept one (solver paths, prior builder).
4. **Data paths and external code.**
5. **Namespace** `+bmm4d/+module`.
6. **Split the cleaning into filter steps.** The regression test must still pass.
7. **Settings and stage functions:** `bmm4d.config`, the five functions, and examples that replace P1–P5.
8. **Package and build files:** `mpackage.json`, `buildfile.m`, README, AGENTS.md.

## 6. Decisions for you

- **A. Synthetic tool:** a folder here (recommended), or a separate package later.
- **B. The cleaning split (step 6):** worth doing now, since it is the stage you expect to change most?
- **C. `compute_scan` and `compute_vp`:** keep as applied single-radar tools, or move them to the paper-only side?
