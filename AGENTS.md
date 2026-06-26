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

## Manuscript Writing

Keep Results sections concise and focused on the new finding shown by each analysis or figure. Avoid repeating methodological details or conclusions already established in earlier sections, and avoid broad statements that are not directly supported by the presented result.

## MATLAB Figure Style

For manuscript and appendix figures exported from MATLAB:

- use the standard MATLAB colormap `parula`, not `turbo`, unless a manuscript-specific exception is explicitly requested
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
