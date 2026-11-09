# Offline basemap for the web map viewer (webmap/)
#
# Natural Earth 1:10m (public domain, R packages rnaturalearth and
# rnaturalearthhires) around the eastern US: land with the lakes cut out,
# coastline, country borders and state/province lines. Written as a MapLibre
# style with inline GeoJSON and no external source, sprite or glyphs, so the
# map runs where external tiles are blocked (the shared Artifact build) and
# offline (?basemap=local in dev).
#
#   Rscript webmap/share/build_basemap.R            # from the repository root
#
# Output: webmap/share/basemap.json (~0.5 MB, committed; rerun only to change it)

suppressPackageStartupMessages({
  library(sf)
  library(rnaturalearth)
})
sf_use_s2(FALSE)                         # planar operations on lon/lat are enough at this scale
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else file.path("webmap", "share")

# Region: the us_east domain (-106..-65 E, 23..50 N) with a wide margin
bbox <- st_bbox(c(xmin = -125, ymin = 12, xmax = -50, ymax = 60), crs = st_crs(4326))
crop <- function(x) suppressWarnings(suppressMessages(st_crop(st_make_valid(x), bbox)))
tolerance <- 0.015                       # deg (~1.5 km), line simplification
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# Land: admin-1 polygons of the three countries (lakes cut out, unlike the
# admin-0 polygons) plus the other countries (Cuba, Bahamas, Central America)
big3 <- c("United States of America", "Canada", "Mexico")
states <- crop(ne_states(country = big3, returnclass = "sf"))
others <- crop(ne_countries(scale = 10, returnclass = "sf"))
others <- others[!others$admin %in% big3, ]
by_country <- lapply(big3, function(name) quiet(st_union(states[states$admin == name, ])))
land <- quiet(st_union(do.call(c, c(by_country, list(st_union(others))))))

# Country borders: country outlines minus the coast (and the crop edge)
coast <- quiet(st_boundary(land))
coast_zone <- quiet(st_buffer(coast, 0.02))
outlines <- quiet(st_union(do.call(c, c(lapply(by_country, st_boundary), list(st_boundary(st_union(others)))))))
borders <- quiet(st_difference(outlines, coast_zone))

# State and province lines: state outlines minus the country outlines
state_lines <- quiet(st_union(st_boundary(states)))
state_lines <- quiet(st_difference(state_lines, st_buffer(outlines, 0.02)))

# Simplify for size and write compact GeoJSON (3 decimals, ~100 m)
geojson <- function(geometry, type) {
  # One MultiPolygon or MultiLineString feature collection. Crops and
  # differences can return collections, so the wanted type is extracted;
  # line pieces are merged first, as two-point segments cannot be simplified.
  # Written directly (GDAL's GeoJSON writer turned the land into a collection).
  geometry <- quiet(st_union(st_collection_extract(st_sfc(geometry, crs = 4326), type)))
  if (type == "LINESTRING") geometry <- quiet(st_line_merge(geometry))
  geometry <- quiet(st_simplify(geometry, preserveTopology = TRUE, dTolerance = tolerance))
  multi <- if (type == "POLYGON") "MULTIPOLYGON" else "MULTILINESTRING"
  x <- quiet(st_cast(geometry, multi))[[1]]
  ring <- function(m) round(m[, 1:2, drop = FALSE], 3)
  coordinates <- if (type == "POLYGON") lapply(x, function(polygon) lapply(polygon, ring)) else lapply(x, ring)
  feature <- list(type = "Feature", properties = setNames(list(), character(0)),
    geometry = list(type = if (type == "POLYGON") "MultiPolygon" else "MultiLineString", coordinates = coordinates))
  as.character(jsonlite::toJSON(list(type = "FeatureCollection", features = list(feature)), auto_unbox = TRUE, digits = NA))
}
land_json <- geojson(land, "POLYGON")
coast_json <- geojson(st_boundary(st_simplify(land, preserveTopology = TRUE, dTolerance = tolerance)), "LINESTRING")
borders_json <- geojson(borders, "LINESTRING")
states_json <- geojson(state_lines, "LINESTRING")

# MapLibre style: water background, land, then (above the birds, marked
# "bmm:above-birds" so the viewer inserts the density layer below them)
# coast, state lines and country borders
source <- function(json) sprintf('{"type": "geojson", "data": %s, "attribution": "Natural Earth"}', json)
style <- sprintf('{
  "version": 8,
  "name": "BMM-4D offline basemap (Natural Earth 1:10m)",
  "sources": {"land": %s, "coast": %s, "states": %s, "borders": %s},
  "layers": [
    {"id": "water", "type": "background", "paint": {"background-color": "#dae6ee"}},
    {"id": "land", "type": "fill", "source": "land", "paint": {"fill-color": "#f4f4f1", "fill-antialias": true}},
    {"id": "coast", "type": "line", "source": "coast", "metadata": {"bmm:above-birds": true},
      "paint": {"line-color": "#7f95a5", "line-width": ["interpolate", ["linear"], ["zoom"], 3, 0.5, 8, 1.2]}},
    {"id": "states", "type": "line", "source": "states", "metadata": {"bmm:above-birds": true},
      "paint": {"line-color": "#55636e", "line-opacity": 0.55, "line-width": ["interpolate", ["linear"], ["zoom"], 3, 0.4, 8, 1]}},
    {"id": "borders", "type": "line", "source": "borders", "metadata": {"bmm:above-birds": true},
      "paint": {"line-color": "#3c4852", "line-opacity": 0.8, "line-width": ["interpolate", ["linear"], ["zoom"], 3, 0.8, 8, 1.8]}}
  ]
}', source(land_json), source(coast_json), source(states_json), source(borders_json))
out_file <- file.path(out_dir, "basemap.json")
writeLines(style, out_file)
cat(sprintf("Wrote %s (%.0f kB)\n", out_file, file.size(out_file) / 1024))
