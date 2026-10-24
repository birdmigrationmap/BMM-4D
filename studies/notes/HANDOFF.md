# Handoff: prod model development (24 October 2026)

State of the `prod/` dev work on branch `dev/production` (nothing committed since PR #10; do not commit without Raphaël's explicit approval). Read in this order:

1. `pipeline/README.md`: the pipeline (P1–P5, P2b, P4a, P4b), current us_east defaults, radar-imprint history.
2. `studies/sandbox/README.md`: the fast sandbox (Midwest part of the night, 28 radars, 6 frames, 4 km × 100 m), its scores and the lean-solver option tests (data weights, prior strength, GP prior).
3. `studies/notes/MODEL_REVIEW.md`: literature, BMM-US review, diagnostics on the 21–22 October 2026 night, and what the current model gets wrong (sections 1–4).
4. `studies/two_level/README.md`: the two-level separable model, set aside (see below).

## Environment
- MATLAB **R2026b** through the MCP tool, run from the project root: `cd(project_root); addpath(project_root); project_root = setup_paths(); addpath(fullfile(project_root,"prod","functions"), fullfile(project_root,"prod","sandbox"))`.
- Data and caches live outside OneDrive in `/Users/rafnuss/Data/BMM-within/prod/` (cases `us_east`, `sandbox`, `us_east_near0_zero150`, cleaned volume caches). Disk space is short (~20 GB free).
- Long runs time out the MCP tool after 30 min but keep running; watch `stage_log.csv` or the output file with a background shell loop.

## Where we are
- **Production** (`scripts/P4_invert.m`, `solve_advection_lean`): 3D density on active cells (4 km × 100 m, ground to 3 km AGL), radar-profile prior (Gaussian 150 km, weight 0.3), 3D Laplacian (0.3), advection (1, relaxed at takeoff and landing), PCG then projected FISTA for non-negativity. Data cuts instead of an observation model: rows within 25 km dropped, no-echo zeros only within 75 km. us_east night: 13 min solve, 7–8 GB, radar-distance bias 0.036.
- **Radar imprint** (`functions/radar_distance_bias.m`) is the main open issue. It comes mostly from suppressed near-range readings (elevation dependent, to 30–40 km at 2° and above) and from how the low layers far from radars are filled. Regularization weights alone trade imprint against fit and detail.
- **Two-level separable model** (`experiments/two_level/`): removed the imprint in the sandbox and predicted withheld radars better, but was worse than production on the full night (bias 0.098, dark spot at each radar from an under-corrected near-range curve) and 3× slower. Set aside on 24 October 2026; its code runs from that folder.

## Direction (Raphaël, 24 October 2026)
The main goal is to understand why the radar observations are not consistent with what birds are known to do, i.e. why the reconstruction shows a radar-distance artefact.
1. **Fast iteration first.** Keep the original 3D model (no separable or Gaussian structure) and drop the non-negativity step that makes it slow (plain PCG); add it back later. *Done:* plain PCG is the sandbox default; same scores, 1.6× faster (`sandbox/SB_fast_solver.m`).
2. **Work on the observation model and the prior/smoothing** in that fast model.
3. **A video for every test,** so progress can be judged visually. *Done:* every `sandbox_run` writes `sandbox/figures/runs/<name>.mp4` and appends to `sandbox/figures/sandbox_log.csv`.
4. **Clear metrics of success:** radar-distance bias, trusted misfit, detail, radar cross-validation (`sandbox/sandbox_fold.m`). `sandbox_summary` compares a set of runs.
5. **A synthetic case** with the expected structure (smooth large-scale field, vertical profiles, takeoff and landing), observed through the forward model with the suspected observation errors (near-range suppression, detection limit, multiplicative error), to check whether the synthetic observations reproduce the patterns seen in the real ones.
