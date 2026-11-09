# Web map viewer

The one place to look at inversion results: a local, zoomable, animated map of any inversion, on the laptop or the cluster, production or test. It shows the column density on a MapLibre basemap, morphed between the 30-min frames along the birds' wind, with particles (moving at the birds' speed) or arrows (length ∝ speed) showing the flight direction.

Use it instead of writing a viewer, a video or a map script to check a run. When it lacks a view you need, extend it here.

| Piece | What it does |
|---|---|
| `bmm4d.export_webmap` (`toolbox/+bmm4d/export_webmap.m`) | Writes one run as `<c.data_dir>/webmap/<run>/` (frames, birds' wind, `meta.json`) and prints its link |
| `serve.py` | Lists every export on the laptop and the cluster share, newest first, and serves the viewer on http://127.0.0.1:8642/ |
| `app/` | The viewer (HTML, JavaScript, MapLibre) |
| `share/` | One night as a self-contained page for sharing (claude.ai Artifact or any static host), with the offline basemap |

## Look at an inversion

1. **Export the run** (MATLAB, seconds to a minute, ~10 MB for a us_east night), right after the inversion, on the laptop or in a cluster job:

   ```matlab
   [active_field, active_id] = bmm4d.invert(c);
   bmm4d.export_webmap(c, active_field, active_id=active_id, run="shared_prior", note="shared prior, 2 relaxation rounds")
   ```

   It writes `<c.data_dir>/webmap/<run>/` next to the case's products and prints the link:

   ```text
   Web map: http://127.0.0.1:8642/?data=local/prod/keywest/webmap/inversion
   ```

   `examples/run_case.m` does this as its last stage (run `"inversion"`), so every pipeline run, on the laptop or through `hpc/submit.sh`, ends with its link (in the job's `slurm-*.out` on the cluster). For an inversion saved earlier, `bmm4d.export_webmap(c)` reads `c.data_dir/inversion.mat`.

2. **Start the viewer** on the laptop (from the repository root; it runs until you stop it with Ctrl-C):

   ```bash
   python3 webmap/serve.py
   ```

3. **Open the link**, or <http://127.0.0.1:8642/> for the newest run. The menu at the top lists every run, newest first, grouped by root and case; ↻ looks for new ones. The panel shows the run's note, when and where it was exported (host and commit, or cluster job) and its id. `?data=<id>` opens a run; the URL hash keeps the view (`#zoom/lat/lon`).

**Cluster runs:** a job exports into the cluster's products (`/mnt/ceph/.../BMM-4D/products/<...>/webmap/<run>/`) and prints a `cluster/...` link. The viewer reads them through the Vogelzug share, also from home over VPN: mount it first (`open "smb://rn@pallidus.vogelwarte.ch/UNIT_Vogelzug"`). If the share isn't mounted, its runs are simply not listed.

### `bmm4d.export_webmap(c, field, ...)`

- `c`: the case settings (`case_settings`, `bmm4d.config`); reads `c.data_dir/case.mat` (grid, radar scans, radars) unless `case_data` gives them (e.g. the output of `bmm4d.build_observations`).
- `field`: `[n_active n_time]` density in birds km⁻³ on the active cells, as `bmm4d.invert` returns it; `active_id` the cells (rebuilt from the case masks if omitted). Without `field`, the case's `inversion.mat` (`active_field`, `active_id`).
- `run`: folder name and label (default `"inversion"`). A rerun with the same name overwrites the export. Give each variant its own name.
- `note`: one line shown in the viewer: what the run tests.
- `out_dir`: another folder (the link then needs `serve.py --root`); `meta_extra`: fields added to `meta.json` (e.g. solve time); `column_floor`, `wind_limit`: encodings.

### Where the viewer looks (`serve.py`)

| Root | Folder | Link prefix |
|---|---|---|
| `local` | `local_data_dir` of `local_settings.json` (this laptop's products) | `local/...` |
| `cluster` | `cluster_products_dir` of `local_settings.json`: the cluster's `local_data_dir` as mounted here, `/Volumes/UNIT_Vogelzug/BMM-4D/products` | `cluster/...` |
| more | `--root name=/some/folder` (repeatable) | `name/...` |

An export is any folder with `meta.json` and `col_00.png` up to five levels below a root (folders `volumes` and `code` are skipped). The list is cached for a minute; ↻ rescans. The first scan after mounting the share takes ~10 s, later ones under a second. `--port` changes the port (default 8642). The server only listens on 127.0.0.1.

## Export format

One folder per run:

- `col_XX.png` (one per frame): column density in birds km⁻² (cells weighted by their air share, as the MATLAB maps), on a grid with the model's longitudes and rows regular in Web Mercator y, so the map places it without distortion. Red: `1 + round(254 log(1 + c / floor) / log(1 + c_max / floor))` (floor 5 birds km⁻², 2–3% steps); green: data mask. Outside the mask the red channel repeats the nearest data pixel, so bilinear filtering does not darken the edge.
- `wind_XX.png` (one per interval): the birds' mean wind, the layer winds of the interval (radar profiles of both frames, spread to the grid with the case's velocity weights, Gaussian by default, as the advection) averaged over the layers with the birds of both frames as weights. Red u, green v, ±50 m/s on 0–255.
- `meta.json`: `format` (`bmm4d-webmap-1`), `case`, `run`, `note`, `created` (UTC), `host`, `commit` (laptop; `+` if the tree had uncommitted changes), `job` (cluster job folder); times, bounds, encodings, the 99.5% column (default colour limit), the sun limit, parula, the radars.

The us_east night is 26 frames of 900 × 760 pixels (4 km), ~10 MB. 
## Viewer (`app/`)

- `main.js`: a WebGL custom layer draws the night under the basemap labels. Between frames i and i+1 at fraction τ, the fragment shader takes frame i at the departure point of a parcel that moved τ·dt with the wind, and frame i+1 at the point it reaches in (1−τ)·dt, wind taken at the half-way point (two stages), and blends them (1−τ, τ). This is the morph of the MATLAB video (`studies/sandbox/SB_us_east_video.m`), with the column-mean bird wind instead of each layer's wind. "Advect between frames" off gives a plain cross-fade.
- Night only: data fade out where the sun is above the case's limit (−6°), as in the video (`(limit + 1 − elevation) / 2`); sun position from a low-precision ephemeris in the shader.
- Flow (Settings → Flow, or `?flow=particles|arrows|off`), on a canvas redrawn after every map frame:
  - Particles move at the birds' own speed. When the clock advances by s simulated seconds (real time × playback speed, or a step of the slider or the arrow keys), a particle moves by wind × s, in Mercator units as the shader's morph (1 / (EARTH cos lat) per metre, wind at the half-way point). Particles and the advected density therefore move together at any zoom and playback speed; when paused they stand still. Each particle is a streak back along its wind over an hour of flight (the streak length also shows the speed), fading in and out over its first and last 20 minutes. Particles start where there are birds in the view (probability by column density, night only), live 2–5 simulated hours, and restart when they leave the birds, the night or the view. Their number follows the birds in the view: the Particles setting (default 4 000) is the count at the night's maximum, so dusk and dawn show fewer. Whenever the nearest frame or the view changes, the view's 8 × 8 blocks holding more particles than their share of the birds restart the excess in the blocks short of particles, so the particles keep following the density (night spreading west, a zoom out) instead of staying where they started.
  - Arrows: a 44-px screen grid; each arrow points along the birds' mean wind of the interval and its length is proportional to the speed (2.6 px per m/s at every zoom; the legend shows 10 m/s; 19 m/s, the fastest of the us_east night, is one grid step). Drawn where the displayed (morphed) column is at least 1 bird km⁻², more opaque where denser, faded with the sun as the density.
- Colour: parula, linear from 0 to the colour limit (slider, log scale; kept when switching runs, so runs compare on one scale); opacity; optional fade of low densities; radar sites; under the cursor, the density as drawn (morphed) and the birds' speed and heading.
- Keys: space plays / pauses, ← → step 5 min.
- Basemap: OpenFreeMap Positron (no key). `?mapbox_token=<token>` uses Mapbox GL JS and the Mapbox Light style instead (same layer code, Mercator projection). `?basemap=local` uses the offline Natural Earth style `share/basemap.json` (`share/build_basemap.R`).
- The map cannot be pitched (`maxPitch: 0`): the flow overlay maps the image to the screen by an affine transform.

Check of the particle speed (us_east_shared, frame 8 → 9, browser console). Particles were advanced in 30-s steps by the app's own code. After one 30-min interval they travelled 14.7 km on average, and ended 0.12 km on average (max 1.4 km) from the point whose frame-8 density the shader draws there. In the rendered GPU image, a density feature moved 93 px against 92 px for a particle at zoom 8 (21 against 23 px at zoom 6). Playback at 10 min/s, 30 min/s and 2 h/s gives the same 0.463 km per simulated minute.

## Sharing one night

`share/build_share.py` builds one exported run as a self-contained page for a claude.ai Artifact (private link, shared from its Share menu). The Artifact's security policy allows scripts from a few CDNs only, and no external stylesheet, tile or fetch. The build therefore inlines the CSS (with `share/maplibre-lite.css`, the parts of the MapLibre stylesheet the viewer uses) and `main.js`, loads MapLibre 4.7.1 from unpkg, and uses the offline basemap. One night is ~10.5 MB (50 kB page, 0.5 MB basemap, 51 PNG). Artifact limits are 64 MB per publish and 256 MB per version, so a handful of nights fit in one artifact.

```bash
python3 webmap/share/build_share.py --export <local_data_dir>/prod/<case>/webmap/<run> --name <case> --out <folder>
```

`share/basemap.json` is committed; `Rscript webmap/share/build_basemap.R` rebuilds it (R packages sf, rnaturalearth, rnaturalearthhires, jsonlite). Publish `<folder>/index.html` with the files under `<folder>/data/` at the same paths. The same folder works on any static host. Do not use Git LFS (Pages serves the pointer files) or release assets (no CORS headers, so the PNGs cannot be read into the canvas).

## Limits

- The morph uses one wind per column (the birds' mean); the video advects each layer with its own wind. Where layers fly different ways, the in-between frames differ.
- Column only: no altitude view yet (one more image per layer band and frame).
- All frames are held in memory (~10 MB of PNG, ~35 MB decoded for us_east); a whole season would need tiles or streaming.
- The flow uses the wind of the current interval (as the morph), so particle directions and arrows change at the 30-min frame times.
- The export needs a grid in longitude and latitude and the case's radar scans (for the wind).
