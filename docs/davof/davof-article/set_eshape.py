#!/usr/bin/env python3
"""Replace the E_shape column of the transport table bodies by the re-scored values.

The shape error of version 1 (tag davof-v1.0) was scored against the cell fractions of the
isosurface cut of the vertex values, whose construction floor on the exact state is
3.6e-4 m^3 at N = 32. The written end fields were re-scored after version 1 against the exact
sphere--box cell fractions (fork davof-swept 1f3f238, volumeFractionError exactState). This
script takes the new N = 128 values from data/eshape/preprint_eshape_map.csv and writes them
into the E_shape column of table_advection_sphere.tex and table_advection_curved.tex; every
other cell is unchanged. Run it before prepare_tables.py (on Lichtenberg, not on the laptop):
    python3 set_eshape.py && python3 prepare_tables.py
"""
import csv

MAP = "data/eshape/preprint_eshape_map.csv"
TABLES = {"table_advection_sphere.tex": ("table_advection_sphere", 8),
          "table_advection_curved.tex": ("table_advection_curved", 9)}
NOTE = ("% E_shape column re-scored against the exact sphere (fork davof-swept 1f3f238) and "
        "written by set_eshape.py from " + MAP)

new = {}
with open(MAP) as fh:
    for r in csv.DictReader(fh):
        if r["N"] == "128" and r["where"].startswith(("table_advection_sphere", "table_advection_curved")):
            key = (r["where"].split()[0], r["row"], float(r["Co"]))
            new[key] = float(r["new"])

for path, (where, col) in TABLES.items():
    lines = open(path).read().splitlines()
    lines = [l for l in lines if not l.startswith("% E_shape column re-scored")]
    name, done = None, 0
    for i, line in enumerate(lines):
        body = line.rstrip()
        if not body.endswith(r"\\") or "&" not in body:
            continue
        cells = [c.strip() for c in body[:-2].split("&")]
        if cells[0]:
            name = cells[0]
        if len(cells) <= col:
            continue
        try:
            co = float(cells[1])
        except ValueError:
            continue
        key = (where, name, co)
        if key in new:
            cells[col] = "%.2e" % new[key]
            lines[i] = " & ".join(cells) + r" \\"
            done += 1
    lines.insert(1, NOTE)
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(path, "cells replaced:", done)
