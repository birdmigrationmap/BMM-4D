# Terrain beam-blockage maps

`beam_blockage_script.m` builds the static terrain-horizon maps used by
`toolbox/radar/read_rd.m` to exclude blocked radar gates.

## Inputs

- `data/reference/radar.csv`: radar latitude, longitude, antenna height, and station name.
- `data/reference/DEM/dem.mat`: the project DEM (`dem`) and its latitude and longitude coordinates (`dem_lat`, `dem_lon`).

## Method

For each radar, the script samples the DEM onto a polar grid with 250 m range
spacing from 1 to 500 km and 1 degree azimuth spacing. It converts the terrain
height at every polar cell to the corresponding elevation angle from the radar
antenna with `slant2range`. Taking the cumulative maximum along each azimuth
produces the terrain-horizon elevation: the minimum beam-centre elevation that
clears every nearer terrain cell.

The output is therefore an elevation-angle map, not an estimate of the fraction
of beam power occulted. `slant2range` uses its default curved-Earth, 4/3
effective-radius propagation model. The maps include terrain represented by the
DEM but not buildings, vegetation, diffraction, side lobes, or
scene-specific atmospheric-refraction variability.

## Outputs and use

The script writes `data/reference/beam_blockage/beam_blockage.mat`, containing
`bb(range, azimuth, radar)`, `bb_range`, `bb_az`, and `bb_radar`, plus one
station-specific map per radar. These files are reference data and are ignored
by Git.

During scan cleaning, `read_rd.m` nearest-neighbour interpolates a station map
to each sweep's native polar gates. A gate is excluded when the terrain-horizon
elevation is greater than the lower half-power edge of that sweep:

```matlab
terrain_horizon_elevation > sweep_elevation - vertical_half_beamwidth
```

The sparse forward operator does not apply a second blockage mask; it only uses
the observations that remain after this cleaning step.

## Relation to published occultation maps

This approach shares the DEM-based terrain-screening principle of Kucera,
Krajewski, and Young (2004) and Krajewski, Ntelekos, and Goska (2006), but it
is intentionally simpler: those studies model elevation-specific partial and
complete beam blockage/power loss, whereas this workflow conservatively drops
gates whose lower half-power beam edge is not terrain-clear. Do not describe
these maps as a 10% occultation mask.
