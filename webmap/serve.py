#!/usr/bin/env python3
"""Web map viewer of BMM-4D inversions, on localhost.

Lists every web map export (a folder with meta.json and col_00.png, written
by bmm4d.export_webmap) found under the data roots, newest first, and serves
the viewer (app/) at / and the exports at /data/<root>/<path>/.

    python3 webmap/serve.py                            # http://127.0.0.1:8642/
    python3 webmap/serve.py --root tests=/some/folder  # an extra root (repeatable)

Default roots, from local_settings.json:
  local    local_data_dir (this laptop's products)
  cluster  cluster_products_dir (the cluster's local_data_dir as mounted here,
           /Volumes/UNIT_Vogelzug/BMM-4D/products; listed when the share is
           mounted, also over VPN)

A run opens with http://127.0.0.1:8642/?data=<root>/<path>, the link
bmm4d.export_webmap prints. The list is read again after RESCAN_SECONDS or with
/data/index.json?refresh=1 (the viewer's Refresh button).
"""
import argparse
import http.server
import json
import os
import pathlib
import time
import urllib.parse

WEBMAP = pathlib.Path(__file__).resolve().parent
APP = WEBMAP / "app"
BASEMAP = WEBMAP / "share" / "basemap.json"        # offline basemap (?basemap=local)
PROJECT = WEBMAP.parent
MAX_DEPTH = 5                               # <root>/prod/<case>/webmap/<run>/ is depth 4
SKIP = {"volumes", "code", "__pycache__"}   # large folders without exports
RESCAN_SECONDS = 60


def default_roots():
    # local_data_dir and cluster_products_dir of local_settings.json
    settings_file = PROJECT / "local_settings.json"
    settings = json.loads(settings_file.read_text()) if settings_file.is_file() else {}
    roots = {}
    for name, key in [("local", "local_data_dir"), ("cluster", "cluster_products_dir")]:
        if settings.get(key):
            roots[name] = pathlib.Path(settings[key]).expanduser()
    return roots


def find_exports(root, depth=0):
    # Folders with meta.json and col_00.png, without looking inside them or into SKIP
    try:
        entries = list(os.scandir(root))
    except OSError:
        return []
    names = {e.name for e in entries}
    if "meta.json" in names and "col_00.png" in names:
        return [pathlib.Path(root)]
    found = []
    if depth < MAX_DEPTH:
        for e in entries:
            if e.is_dir(follow_symlinks=False) and e.name not in SKIP and not e.name.startswith("."):
                found += find_exports(e.path, depth + 1)
    return found


def describe(root_name, root, folder):
    # Index entry of one export: id for the URL, case, run, note, provenance
    meta = json.loads((folder / "meta.json").read_text())
    relative = folder.relative_to(root).as_posix()
    parts = relative.split("/")
    in_case = len(parts) >= 3 and parts[-2] == "webmap"      # <...>/<case>/webmap/<run>
    mtime = (folder / "meta.json").stat().st_mtime
    return {
        "id": f"{root_name}/{relative}", "root": root_name,
        "case": meta.get("case") or (parts[-3] if in_case else parts[-1]),
        "run": meta.get("run") or parts[-1],
        "note": meta.get("note", ""),
        "created": meta.get("created") or time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(mtime)),
        "host": meta.get("host", ""), "commit": meta.get("commit", ""), "job": meta.get("job", ""),
        "frames": len(meta.get("times", [])),
        "start": meta["times"][0] if meta.get("times") else "",
    }


class Index:
    def __init__(self, roots):
        self.roots, self.entries, self.time = roots, [], 0.0

    def get(self, refresh=False):
        if refresh or time.time() - self.time > RESCAN_SECONDS:
            entries = []
            for name, root in self.roots.items():
                for folder in find_exports(root):
                    try:
                        entries.append(describe(name, root, folder))
                    except (OSError, ValueError, KeyError) as error:
                        print(f"skipped {folder}: {error}")
            self.entries = sorted(entries, key=lambda e: e["created"], reverse=True)
            self.time = time.time()
        return self.entries


def make_handler(roots, index):
    class Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(APP), **kwargs)

        def end_headers(self):
            # Development server: always read the latest export
            self.send_header("Cache-Control", "no-cache")
            super().end_headers()

        def translate_path(self, path):
            # /data/<root>/<path> -> file under that root; /data/basemap.json -> share/basemap.json
            route = urllib.parse.unquote(path.split("?", 1)[0].split("#", 1)[0])
            if route == "/data/basemap.json":
                return str(BASEMAP)
            if route.startswith("/data/"):
                name, _, rest = route[len("/data/"):].partition("/")
                root = roots.get(name)
                if root is not None:
                    target = (root / rest).resolve()
                    if root.resolve() in target.parents:
                        return str(target)
                return str(APP / "__missing__")
            return super().translate_path(path)

        def do_GET(self):
            route, _, query = self.path.partition("?")
            if route == "/data/index.json":
                body = json.dumps(index.get("refresh=1" in query)).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            super().do_GET()

        def log_message(self, format, *args):
            pass                                    # quiet: one line per image otherwise

    return Handler


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", action="append", default=[], metavar="NAME=FOLDER", help="extra data root")
    parser.add_argument("--port", type=int, default=8642)
    args = parser.parse_args()
    roots = default_roots()
    for item in args.root:
        name, _, folder = item.partition("=")
        roots[name] = pathlib.Path(folder).expanduser()
    roots = {name: root for name, root in roots.items() if root.is_dir()}
    index = Index(roots)
    start = time.time()
    count = len(index.get(True))
    print(f"{count} exports in {time.time() - start:.1f} s under " +
          ", ".join(f"{name} = {root}" for name, root in roots.items()))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", args.port), make_handler(roots, index))
    print(f"Web map on http://127.0.0.1:{args.port}/")
    server.serve_forever()
