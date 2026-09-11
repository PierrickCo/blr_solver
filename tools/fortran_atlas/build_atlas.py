#!/usr/bin/env python3
"""
build_atlas.py -- generate a self-contained, interactive HTML "atlas" of a
Fortran codebase's Makefile + source dependency graph: executables -> source
files -> subroutines/functions/programs, with a per-routine call graph,
clickable source code, and search.

Reuses analyze_deps.py (same folder) for the actual parsing, then injects
its result as JSON into template.html (also same folder) to produce ONE
standalone HTML file -- no server, no companion files, works from a plain
file:// open.

Usage (run from the target codebase's own root -- the directory the given
Makefile paths are relative to):

    python3 /path/to/fortran_atlas/build_atlas.py OUTPUT.html Makefile [Makefile2 ...]
    python3 /path/to/fortran_atlas/build_atlas.py OUTPUT.html Makefile --title "My Project Atlas"

Re-run any time the source changes; purely read-only on the inspected repo.
"""
import argparse
import json
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import analyze_deps as ad  # noqa: E402


def build_atlas(output_path, seed_makefiles, title=None):
    result = ad.analyze(seed_makefiles)
    n_exec = len(result['executables'])
    n_sources = result.pop('n_sources')
    n_routines = result.pop('n_routines')

    template_path = os.path.join(_HERE, 'template.html')
    with open(template_path, 'r', encoding='utf-8') as f:
        html = f.read()

    display_title = title or 'Fortran Dependency Atlas'
    mark = (title or 'atlas').lower().replace(' ', '.')
    subtitle = 'genere depuis ' + ' + '.join(seed_makefiles)

    html = html.replace('__ATLAS_TITLE__', display_title, 1)
    html = html.replace('__ATLAS_MARK__', mark, 1)
    html = html.replace('__ATLAS_SUBTITLE__', subtitle, 1)

    data = json.dumps(result, separators=(',', ':'), sort_keys=True)
    data = data.replace('</', '<\\/')  # keep </script> out of the JSON text
    marker = '<script>\n"use strict";'
    assert marker in html, 'marker not found in template.html'
    inject = '<script type="application/json" id="depgraph-data">' + data + '</script>\n' + marker
    html = html.replace(marker, inject, 1)

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"Wrote {output_path}: {n_exec} executables, {n_sources} source "
          f"files, {n_routines} routines. Open it directly in a browser "
          f"(no server needed).")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('output', help='Path of the HTML file to write')
    parser.add_argument('makefiles', nargs='+', help='One or more seed Makefiles')
    parser.add_argument('--title', help='Atlas title, shown as the page title and header (default: "Fortran Dependency Atlas")')
    args = parser.parse_args()
    build_atlas(args.output, args.makefiles, args.title)


if __name__ == '__main__':
    main()
