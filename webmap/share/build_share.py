#!/usr/bin/env python3
"""Build the shareable web map (Artifact page) of one exported night.

The Artifact page may load scripts from a few CDNs only, and no external
stylesheet, tile or fetch. The build therefore inlines app/style.css, a
minimal MapLibre stylesheet (share/maplibre-lite.css) and app/main.js, loads
MapLibre from unpkg, and uses the offline Natural Earth basemap
(share/basemap.json, from share/build_basemap.R).

    python3 webmap/share/build_share.py \
        --export /Users/rafnuss/Data/BMM-within/prod/<case>/webmap/<run> --name <case> --out <folder>

<folder>/index.html is the page (no <html>/<head>/<body>: the Artifact host
wraps it); <folder>/data/ holds basemap.json and the night as data/<name>/
(meta.json, col_XX.png, wind_XX.png, copied from the export folder that
bmm4d.export_webmap wrote). Publish the page with the data files at the same
relative paths. The same folder also works on any static host.
"""
import argparse
import datetime
import json
import pathlib
import shutil

WEBMAP = pathlib.Path(__file__).resolve().parent.parent      # webmap/
MAPLIBRE = "https://unpkg.com/maplibre-gl@4.7.1/dist/maplibre-gl.js"


def night_title(meta):
    # "Night of 21–22 October 2026" from the first and last frame (UTC)
    first, last = (datetime.datetime.fromisoformat(t.replace("Z", "+00:00")) for t in (meta["times"][0], meta["times"][-1]))
    if first.date() == last.date():
        return f"Night of {first:%-d %B %Y}"
    return f"Night of {first:%-d}–{last:%-d %B %Y}" if first.month == last.month else f"Night of {first:%-d %B}–{last:%-d %B %Y}"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--export", required=True, help="export folder of the night (bmm4d.export_webmap)")
    parser.add_argument("--name", required=True, help="short name of the night in the page (data/<name>/)")
    parser.add_argument("--basemap", default=str(WEBMAP / "share" / "basemap.json"),
                        help="offline basemap (share/build_basemap.R)")
    parser.add_argument("--out", required=True, help="output folder (replaced)")
    args = parser.parse_args()
    export, out = pathlib.Path(args.export), pathlib.Path(args.out)

    # Data files at the paths the viewer reads: data/basemap.json, data/<run>/...
    if out.exists():
        shutil.rmtree(out)
    (out / "data").mkdir(parents=True)
    shutil.copy2(args.basemap, out / "data" / "basemap.json")
    shutil.copytree(export, out / "data" / args.name, ignore=shutil.ignore_patterns(".*"))
    meta = json.loads((export / "meta.json").read_text())
    run = f"{meta['case']}, run {meta['run']}" if "case" in meta else args.name

    # Viewer settings for the shared build (read by main.js as window.WEBMAP)
    config = {
        "runs": [args.name], "basemap": "local", "hash": False, "attribution": False,
        "maxBounds": [[-128, 10], [-47, 60]],
        "title": night_title(meta),
        "credit": f"Column bird density reconstructed with BMM-4D from NEXRAD weather radars "
                  f"(development run: {run}). Basemap: Natural Earth.",
    }

    # Page: the app's body without its script tags, with styles and scripts inlined
    html = (WEBMAP / "app" / "index.html").read_text()
    body = html.split("<body>", 1)[1].split('<script src="https://unpkg.com/maplibre-gl', 1)[0].strip()
    css = (WEBMAP / "share" / "maplibre-lite.css").read_text() + (WEBMAP / "app" / "style.css").read_text()
    js = (WEBMAP / "app" / "main.js").read_text()
    page = "\n".join([
        '<meta charset="utf-8">',               # the Artifact host adds it; plain static hosts need it
        "<title>Eastern US Migration Map</title>",
        f"<style>\n{css}</style>",
        body,
        f'<script src="{MAPLIBRE}"></script>',
        f"<script>window.WEBMAP = {json.dumps(config, ensure_ascii=False)};</script>",
        f"<script>\n{js}</script>",
    ])
    (out / "index.html").write_text(page)
    files = sorted(str(p.relative_to(out)) for p in (out / "data").rglob("*") if p.is_file())
    size = sum((out / f).stat().st_size for f in files)
    print(f"{out / 'index.html'}: {len(page) / 1024:.0f} kB page, {len(files)} data files, {size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
