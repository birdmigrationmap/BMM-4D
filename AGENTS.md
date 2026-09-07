# AGENTS.md

## Project Style

This project is an early-stage research codebase for bird migration model reconstruction.

Prioritize:

- linear scripts
- interactive execution
- transparent intermediate objects
- easy debugging
- reproducible workflows

Do not optimize for backward compatibility unless explicitly requested.

When replacing an old implementation, update the callers instead of keeping legacy aliases or deprecated parameter names.

## Canonical Grid Convention

Use one Cartesian grid convention across the project.

- Use the variable name `g` for the reconstruction grid. Do not use a local variable named `grid` in active MATLAB code because it conflicts with MATLAB's `grid` plotting command and creates avoidable Code Analyzer issues.
- Arrays are ordered as `[x y z]`
- `g.sz = [numel(g.x) numel(g.y) numel(g.z)]`
- `reshape(vector, g.sz)` is the standard vector-to-grid conversion
- `sub2ind(g.sz, x_id, y_id, z_id)` is the matching linear indexing rule
- If a figure needs a transpose or `permute(...)`, do it only in the local plotting code and comment why

## Script Structure

Prefer top-to-bottom scripts that can be run section by section in MATLAB.

For project setup in scripts, use:

```matlab
project_root = setup_paths();
```

Capture the return value when the script needs to construct paths or call `run(...)` on other scripts. Do not call `setup_paths();` without assigning `project_root` in new or edited numbered workflow scripts.

Do not use manual root discovery such as:

```matlab
fileparts(fileparts(fileparts(mfilename("fullpath"))))
```

Use explicit section headers such as:

```matlab
%% Load packages and paths
%% Read radar data
%% Build forward model
%% Reconstruct vertical profile
%% Reconstruct 3D field
%% Plot results
```

Avoid hiding the main workflow inside nested functions or abstractions.

## MATLAB Execution Lessons

For MATLAB MCP and Codex runs:

- Treat the project root as the working directory for workflow execution.
- Before running a numbered script through MATLAB MCP, set the current folder to the project root and add the root to the MATLAB path.
- Prefer:

```matlab
cd(project_root)
addpath(project_root)
run(fullfile(project_root, "03_simulator_fit", "scripts", "S03_B_prepare_scans.m"))
```

- Do not rely on `run_matlab_file(...)` for numbered workflow scripts when they call `setup_paths()`, because that tool starts MATLAB in the script folder rather than the project root.
- If a script saves generated artifacts that are consumed downstream, keep the saved variable names aligned with the active codebase and rerun the upstream script after a naming-contract change.
- When a long MATLAB job is launched through MCP, expect tool timeouts before MATLAB necessarily finishes. Check whether output files were updated before assuming the run failed.

## Comments

Comments should explain what is being done and why.

Prefer short comments above logical blocks, for example:

```matlab
% Remove rain-contaminated gates using rhoHV
% Build the vertical profile operator
% Smooth density along range and azimuth
```

Do not add verbose or defensive comments that restate obvious syntax.

## Functions

Keep functions short and explicit.

Only create helpers when the same logic is reused or when the function is a clear conceptual unit.

Prefer simple argument lists over compatibility wrappers.

Do not extract one-off script logic into separate helper files unless that logic is reused across multiple files or multiple numbered sections.

If logic is reused only within one script, prefer inline code first. If repetition becomes large enough to hurt readability, use a local function at the end of that same script instead of creating a separate `.m` file.

For numbered workflow scripts such as `A_*.m`, `B_*.m`, and `C_*.m`, keep the script body linear and move reusable helper functions into shared files under `functions/` instead of defining local functions at the end of the script.

Before adding any new helper for plotting, exporting, formatting labels, or figure layout, first check `functions/plot/` and reuse an existing helper when possible.

Do not introduce local helper functions inside numbered workflow scripts when the logic is already covered by an existing shared helper. Prefer calling the shared helper and applying a few explicit script-level axis adjustments afterward.

## Naming

Use readable variable names.

Fix typos in names when touching code instead of preserving them for compatibility.

## Scientific Workflow

Make assumptions explicit in code.

Keep intermediate objects accessible for plotting, checking, and debugging.

Prefer practical code that is easy to modify during method development.

For synthetic generation and inversion scripts, do not introduce `selected_cases` or similar saved case-selection control variables. Loop over all available cases explicitly in the script. If you need to run one case temporarily, edit the loop locally for that run instead of adding a persistent selection layer to the workflow.

For synthetic case and inversion products, keep saved `.mat` files compact and intentional. Before saving a new variable, ask whether it is needed to rerun downstream analysis without regenerating the expensive synthetic case. If it is only useful for debugging, plotting convenience, or can be rebuilt cheaply from saved controls and derived inputs, do not put it in the main case file.

- Save large dense scientific fields as `single`: `truth_fields`, `truth_sequences`, reconstructed fields, dense prior fields, and cached baseline fields.
- Keep solver operators, grids, scan geometry, metadata tables, controls, and compact observation structs in their natural precision unless there is a measured reason to change them.
- Save observation data in compact form: sparse vectors/matrices for single-radar synthetic observations and compact `observations` structs for mosaic workflows.
- Do not save repeated full radar scan structs, full template radar volumes, elevation kernels, forward operators, or true gridded velocity fields inside synthetic case files.
- Save compact descriptors instead of bulky duplicates: velocity/profile arrays, radar/time identifiers, generation controls, seeds, and lightweight scan metadata. Rebuild minimal scan structs locally when a downstream function needs scan metadata.
- Split heavy products by role: `case*.mat` for truth/observations/controls, `inversion*.mat` for reconstruction outputs, `validation_products.mat` for derived metrics, and figures under `figures/`.
- If a saved file grows unexpectedly large, inspect variable sizes before adding more output and prefer a one-time compaction script over rerunning expensive simulations only to change storage format.

## Manuscript Writing

Keep Results sections concise and focused on the new finding shown by each analysis or figure. Avoid repeating methodological details or conclusions already established in earlier sections, and avoid broad statements that are not directly supported by the presented result.

Treat manuscript and appendix files as publication text, not as code documentation. Do not include function calls, filenames, script variable names, literal assignments, or other source-code syntax unless the user explicitly requests them. Describe methods with scientific notation and publication-appropriate terminology.

Do not add tables, new subsections, or other manuscript structure unless explicitly requested or approved. If the appropriate level of implementation detail or its placement is uncertain, ask before editing the manuscript.

## MATLAB Figure Style

For manuscript and appendix figures exported from MATLAB:

- Assign colormaps by scientific quantity so the same color has the same meaning across workflows:
  - bird density and other primary field values: `parula`
  - error magnitudes such as RMSE and NRMSE: `hot(256)`
  - pattern correlation: `winter(256)`
  - observation support or data sensitivity: `flipud(gray(256))`, with unsupported bins shown as white/blank
  - secondary case descriptors used to color scatter plots: `copper(256)`
- Use light figure backgrounds consistently. Standard non-map figures should export on white or near-white backgrounds and should never rely on transparency that turns dark in viewers or manuscripts.
- For maps with masked regions or land-water separation, render the background explicitly instead of leaving transparent holes:
  - use a light water color such as a pale blue for water or non-target masked areas
  - keep land visually lighter than the data overlay so density and error colors remain dominant
  - avoid dark satellite or dark-theme basemaps for normal manuscript figures unless a figure explicitly needs nighttime imagery
- When a panel uses `AlphaData` to hide unsupported or masked cells, also set the axes background color intentionally so exported PNGs remain light in all viewers.
- For 1D summaries of those same quantities such as lines, histograms, shaded bands, and markers, use colors that match the assigned colormap family rather than unrelated defaults:
  - RMSE and NRMSE summaries: warm orange-red tones consistent with `hot(256)`
  - pattern-correlation summaries: cool green-teal tones consistent with `winter(256)`
  - reference or comparison distributions: neutral gray
- Apply colormaps to individual axes, not the whole figure, when panels show different quantities.
- Reuse shared color limits when directly comparing the same metric between methods or panels.
- Do not use `turbo` unless a manuscript-specific exception is explicitly requested.
- prefer white figure backgrounds
- create exported figures with `create_export_figure(width_px, height_px, ...)` so output size is independent of the display and window manager
- export raster figures with `export_fixed_figure(fig, file, Resolution=300)` rather than calling `exportgraphics` directly; `exportgraphics` crops to content and can produce inconsistent pixel dimensions
- keep exported figures offscreen with `Visible="off"` unless interactive inspection is the purpose of the script section
- treat `width_px` and `height_px` in `create_export_figure` as the intended final pixel dimensions at `Resolution=300`
- size figures explicitly for paper export rather than relying on default screen dimensions or visible figure windows
- use concise titles and axis labels so multi-panel figures remain legible after export
- use `tiledlayout(..., "TileSpacing", "tight", "Padding", "tight")` for figure layouts unless a script has a specific documented reason to use another spacing
- export manuscript figures as high-resolution `.png` by default, using `export_fixed_figure(..., ".png", Resolution=300)` at minimum and `Resolution=600` when finer detail matters
- use vector export only when it is explicitly needed and the figure content is lightweight enough to stay responsive
- when XY and XZ panels are compared in the same figure, keep the horizontal physical span visually consistent across panels

## MATLAB MCP

Use the MATLAB MCP tools for MATLAB work in Codex instead of trying to launch MATLAB from the shell.

Prefer these MCP actions:

- inspect the MATLAB installation and toolboxes
- evaluate MATLAB code directly
- run `.m` scripts
- run MATLAB test files
- use MATLAB code analysis when reviewing or debugging MATLAB code

As of `2026-06-06`, the MATLAB MCP connection was verified from Codex and reported:

- `MATLAB R2026a Update 2`
- `MATLAB MCP Core Server Toolbox 0.1.1`

Do not add instructions that rely on `matlab -batch` or other command-line MATLAB startup from the Codex shell unless that workflow is re-verified and explicitly needed.
