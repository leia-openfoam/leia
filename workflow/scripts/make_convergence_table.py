#!/usr/bin/env python3
"""Consolidate the semi-Lagrangian convergence studies into publication CSVs.

The two semi-Lagrangian METHOD LINES are handled as independent methods, each
with its own study set, docs theme and output prefix:

    --method quadratic   uncachedConv* studies  -> semi-lagrangian-level-set (sl_*)
    --method linear      linearConv*   studies  -> linear-semi-lagrangian-level-set (lsl_*)

Reads the curated ``<study>_errors.csv`` of the uniform-refinement studies
(2D vortex, 3D shear/deformation, hex + polyhedral) and emits tidy CSVs the
paper and slides consume directly:

    <prefix>_convergence.csv
        one row per resolution with ALL error metrics:
        case, mesh, reconstruction, cfl, h, maxCellSize,
        gradientError, shapeError, volumeError, gradientErrorBand,
        gradientErrorHalf, gradientErrorBandHalf, volumeErrorHalf,
        minGradPsiBand, minGradPsiBandHalf
    <prefix>_convergence_orders.csv
        one row per (case, mesh, reconstruction, cfl) with the least-squares
        convergence order (slope of log error vs log h) of EVERY error metric.

plus two LaTeX tabulars (\\input by the articles):
    convergence_orders.tex           the three primary metrics (shape, volume,
                                     gradient) -- unchanged layout
    convergence_orders_extended.tex  all seven convergent metrics, incl. the
                                     band-restricted and half-time (t=T/2,
                                     maximum deformation) variants

Least-squares order only (the uniform refinement ladder makes the slope clean);
no per-triplet Richardson / GCI. Missing studies are skipped with a note.

Usage (from repo root):
    python3 workflow/scripts/make_convergence_table.py [--method quadratic|linear]
"""
import argparse
import csv
import os
import re
import sys

import numpy as np

import method_label  # shared, LaTeX-safe method label
import paths  # thematic docs layout (single source of truth for output dirs)

REPO = paths.REPO

# Per method line: docs theme, output-file prefix, and the study set
# (study dir -> (case label, mesh label)) for the publication tables.
METHODS = {
    "quadratic": {
        "theme": "semi-lagrangian-level-set",
        "prefix": "sl",
        "studies": [
            ("uncachedConv2Dvortex",          "2Dvortex",      "hex"),
            ("uncachedConv3Dshear",           "3Dshear",       "hex"),
            ("uncachedConv3Ddeformation",     "3Ddeformation", "hex"),
            ("uncachedConv3DshearPoly",       "3Dshear",       "poly"),
            ("uncachedConv3DdeformationPoly", "3Ddeformation", "poly"),
        ],
    },
    "linear": {
        "theme": "linear-semi-lagrangian-level-set",
        "prefix": "lsl",
        "studies": [
            ("linearConv2Dvortex",          "2Dvortex",      "hex"),
            ("linearConv3Dshear",           "3Dshear",       "hex"),
            ("linearConv3Ddeformation",     "3Ddeformation", "hex"),
            ("linearConv3DshearPoly",       "3Dshear",       "poly"),
            ("linearConv3DdeformationPoly", "3Ddeformation", "poly"),
        ],
    },
    # The SDPLS source line. Its arms are distinguished by the METHOD label
    # (euler, euler+SDPLS:R/simpleImp, euler+SDPLS:beta/strictNegSp, ...), NOT
    # by the semi-Lagrangian `reconstruction` column, which is blank for every
    # Eulerian run -- hence the groupby override. Without it all six arms of a
    # study collapse into a single series and no fitted order means anything.
    #
    # This entry is the wiring that finally produces the two numbers the SDPLS
    # record has never had: gradientBandOrder (t = T, a pure residual meter on
    # a reversible flow) and gradientBandHalfOrder (t = T/2, maximal stretch).
    "sdpls": {
        "theme": "sdpls-level-set",
        "prefix": "sdpls",
        "groupby": "method",
        "studies": [
            # benchVortexEulerT2 / benchVortexEulerT8 are DELIBERATELY ABSENT.
            # They are reversed (oscillating) runs, and a reversed t=T error is
            # not an accuracy measurement -- see _drop_reversed() below. They
            # stay on disk as evidence of the artifact itself; they contribute
            # no published order.
            #
            # The beta-target sweep. Its arms are separable here ONLY because
            # method_label._beta() renders an off-default beta into the label --
            # without that they collapse to one series per h and the
            # least-squares slope below gets fitted straight through them.
            ("sdplsBetaSweep",         "2Dvortex",      "hex"),
            # Non-reversing vortex: the only 2D study whose t=T column is a
            # measurement rather than a time-reversal symmetry check.
            ("sdplsConv2Dvortex",      "2Dvortex",      "hex"),
            # Order ablation. Its arms differ in the linearUpwind reconstruction
            # gradient, which method_label now renders (`...grad(psi):lsq`), so
            # the four variants stay separate series here.
            ("sdplsOrderAblation",     "2Dvortex",      "hex"),
            ("sdplsConv3Dshear",       "3Dshear",       "hex"),
            ("sdplsConv3Ddeformation", "3Ddeformation", "hex"),
        ],
    },
}

# All convergent error metrics: (CSV column, orders column, tex header).
# The half-time (T/2) metrics probe the maximally deformed interface, where the
# level set is furthest from a signed distance field; the band variants restrict
# |grad(psi)|-1 to the narrow band that the numerics actually use.
METRICS = [
    ("shapeError",           "shapeOrder",            r"shape"),
    ("volumeError",          "volumeOrder",           r"volume"),
    ("gradientError",        "gradientOrder",         r"$E_{\nabla\psi}$"),
    ("gradientErrorBand",    "gradientBandOrder",     r"$E_{\nabla\psi}^{\mathrm{band}}$"),
    ("gradientErrorHalf",    "gradientHalfOrder",     r"$E_{\nabla\psi}(T/2)$"),
    ("gradientErrorBandHalf", "gradientBandHalfOrder", r"$E_{\nabla\psi}^{\mathrm{band}}(T/2)$"),
    ("volumeErrorHalf",      "volumeHalfOrder",       r"volume$(T/2)$"),
]
# Health (not convergence) metrics carried through the per-resolution CSV.
# The band MEAN and MAX join the min because a source relaxing toward a wrong
# TARGET (sdplsBeta drives |grad psi| -> beta - a) shifts the centre of the band
# distribution, which the min barely registers; max - min is the spread, which a
# pure target offset should NOT change. meanStrainBand is `a` itself, so
# `beta - meanStrainBand` can be checked against meanGradPsiBand directly.
HEALTH = ["minGradPsiBand", "minGradPsiBandHalf",
          "meanGradPsiBand", "meanGradPsiBandHalf",
          "maxGradPsiBand", "maxGradPsiBandHalf",
          "meanStrainBand", "meanStrainBandHalf"]


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


# REVERSED FLOW IS NOT AN ACCURACY MEASUREMENT.
# On the oscillating vortex/shear the velocity carries a cos(pi t/tau) factor,
# so the interfacial stretching integrates to zero over the period:
#     \int_0^T n.grad(U).n dt = 0.
# Error accumulated on the forward half is undone on the return half, so the t=T
# state flatters whichever scheme accumulated the most RECOVERABLE error. Not
# hypothetical: at t=T the reversed 2D vortex made sourceless transport look
# 477x more accurate than SDPLS-R, while the SAME case run non-reversing has
# sourceless DIVERGING (order -0.26) and only R converging (+0.74). The reversal
# inverted the ranking outright.
#
# So a row that ran with oscillation ON never contributes a published order.
# Keyed on the RECORDED value rather than on a study name, so a config that
# forgets `OSCILLATION: "off"` is caught by its own data instead of by whoever
# remembers to update a list here.
#
# Rows predating the column (value absent) are deliberately left alone. The
# semi-Lagrangian shape-error benchmark is reversed BY CONSTRUCTION -- shape
# error at t=T is defined against the initial condition, and there is no t=T
# shape reference without the return leg -- so voiding those published tables is
# a separate decision from this one, and not one to take silently.
def _drop_reversed(rows):
    """Split rows into (publishable, excluded-because-reversed)."""
    rev = [r for r in rows if (r.get("oscillation") or "").strip() == "on"]
    keep = [r for r in rows if (r.get("oscillation") or "").strip() != "on"]
    return keep, rev


# STALE-LADDER DETECTION.
# `studies/` is git-ignored, survives config edits, and case directories are
# named by ORDINAL index over the axis product -- so changing N_CELLS renumbers
# them and a stale tree can masquerade as a fresh one. Measured: after the 3D
# ladder became 64/81/102/128/161/203, studies/sdplsConv3Dshear still held a
# complete 3-rung 32/64/128 table from the previous ladder AND a previous code
# version, and the only reason it was not published as the new result is that
# someone happened to read the row count. The missing-study guard cannot see
# this -- the study is present, just wrong.
#
# So compare the resolutions IN THE DATA against the ladder the config declares.
# h is written as 1/N by aggregate.py, so N comes back exactly.
def _config_ladder(study):
    """Declared N_CELLS for a study, or None when it cannot be determined."""
    path = os.path.join(REPO, "config", f"{study}.yaml")
    if not os.path.isfile(path):
        return None
    try:
        import yaml
        with open(path) as fh:
            cfg = yaml.safe_load(fh) or {}
    except Exception:
        return None
    # POLYHEDRAL studies are sized by cfMesh's MAX_CELL_SIZE, not by N_CELLS,
    # and aggregate.py sets h = maxCellSize for them -- so 1/h is not a cell
    # count and comparing it to a ladder is meaningless. Worse, those configs
    # pin `N_CELLS: [64]` as an explicit dummy ("unused by pMesh; pin so the
    # base [32,64] does not multiply runs"), which reads as a one-rung ladder
    # and made this guard reject all four healthy poly studies on first run.
    if (cfg.get("mesh") or "").strip() == "poly":
        return None

    # N_CELLS lives under `axes_override` today, but has moved before. Search
    # the whole tree rather than a guessed key list: a lookup that silently
    # returns None disables this guard, which is exactly how the stale 3D tree
    # got past it on the first attempt.
    def _find(node):
        if isinstance(node, dict):
            for key, val in node.items():
                if key == "N_CELLS" and isinstance(val, (list, tuple)):
                    return {int(v) for v in val}
                found = _find(val)
                if found:
                    return found
        elif isinstance(node, list):
            for item in node:
                found = _find(item)
                if found:
                    return found
        return None
    return _find(cfg)


def _ladder_complaints(study, rows):
    """Resolutions present in the data but absent from the declared ladder."""
    declared = _config_ladder(study)
    if not declared:
        return []
    found = set()
    for r in rows:
        h = _f(r.get("h"))
        if h and h > 0:
            found.add(int(round(1.0 / h)))
    return sorted(found - declared)


def _collapse_duplicates(rows, groupby):
    """Collapse rows identical in (arm, h, CFL) AND in every metric.

    A study that crosses an INERT token against an arm runs that arm twice:
    sdplsConv2Dvortex crosses SOURCE_SCHEME (the source LINEARIZATION) with
    SDPLS_SOURCE=noSource, which has no source to linearize, so the pair is
    bit-identical (band gradient error 1.59990364597 from both) and the arm
    appears at 14 rungs on a 7-rung ladder. The least-squares slope is unchanged
    by exact duplicates, but the published per-resolution table would print
    every rung of that arm twice.

    Collapsing is safe ONLY where the duplicates agree exactly. Differing values
    under one (arm, h) mean the label is hiding a real degree of freedom -- a
    labelling bug -- so those stay visible and are reported, never averaged.
    """
    seen, out, collapsed, conflicts = {}, [], 0, []
    for r in rows:
        key = (r.get(groupby, ""), r.get("h", ""), r.get("cfl", ""))
        sig = tuple(r.get(c, "") for c, _o, _l in METRICS)
        if key not in seen:
            seen[key] = sig
            out.append(r)
        elif seen[key] == sig:
            collapsed += 1
        else:
            conflicts.append(key)
            out.append(r)
    return out, collapsed, conflicts


def _order(pts):
    """Least-squares order (slope of log |e| vs log h) over (h, e) pairs.
    Magnitudes are used so signed residuals (volume error) fit cleanly."""
    pts = [(h, abs(e)) for h, e in pts if h and e is not None and h > 0 and abs(e) > 0]
    if len(pts) < 2:
        return None
    h, e = zip(*pts)
    return float(np.polyfit(np.log(h), np.log(e), 1)[0])


# Stable-envelope detection -- same rule as workflow/scripts/make_sl_3d_fig.py, so
# the order table and the convergence figures agree on which resolutions count.
# The value-fitting reconstructions are stable by construction (the full ladder
# counts); the guard is a safety net that truncates the field-differentiation
# Taylor family at the first resolution where the gradient defect runs away
# (|grad psi|-1 > GRAD_CEILING; the elevated-but-bounded polyhedral baseline stays
# under it) or the shape error turns upward (> running min x SHAPE_TOL). Orders
# are fitted over the stable prefix only.
GRAD_CEILING = 10.0
SHAPE_TOL    = 1.30

# UNDER-RESOLVED COARSE PREFIX. The envelope above drops UNSTABLE FINE rungs --
# the semi-Lagrangian failure mode, where refinement destabilizes. Non-reversing
# flows fail the other way round: the steady vortex winds the interface into a
# filament the coarse meshes cannot represent, so the COARSE end is invalid and
# refinement cures it. Measured on sdplsConv2Dvortex, relative volume error of
# the sourceless baseline: 0.245, 0.079, 0.019, 0.0064 at N = 32/45/64/90. A run
# that has lost a quarter of its phase volume is not a convergence data point
# whatever its gradient error reads.
#
# SCOPED DELIBERATELY: applied only when `shapeError` is blank, which
# aggregate.py does exactly for oscillation=off rows. On a reversed flow the
# shape error is a valid interface-integrity check and the existing envelope
# already covers it, so the semi-Lagrangian and reversed tables cannot be
# touched by this gate.
VOL_CEILING = 0.05


def _stable_hset(grp):
    """Set of stable h for one (reconstruction, CFL) group, and the first
    destabilized h (or None)."""
    pts = sorted(
        ((_f(r.get("h")), _f(r.get("shapeError")), _f(r.get("gradientError")))
         for r in grp),
        key=lambda t: (t[0] if t[0] else 0.0), reverse=True)
    stable, limit, shape_min = set(), None, None
    for h, shape, grad in pts:
        if h is None:
            continue
        ok = True
        if grad is not None and abs(grad) > GRAD_CEILING:
            ok = False
        if shape is not None:
            if shape_min is not None and shape > shape_min * SHAPE_TOL:
                ok = False
            shape_min = shape if shape_min is None else min(shape_min, shape)
        if ok:
            stable.add(h)
        else:
            limit = h
            break

    # Drop the under-resolved COARSE prefix (non-reversing runs only -- see
    # VOL_CEILING). Walk coarse -> fine and discard until the volume error comes
    # under the ceiling; everything finer is kept, since volume error decreases
    # monotonically once the interface is represented at all.
    vol = sorted(((_f(r.get("h")), _f(r.get("volumeError")),
                   (r.get("shapeError") or "").strip())
                  for r in grp),
                 key=lambda t: (t[0] if t[0] else 0.0), reverse=True)
    for h, v, shape_str in vol:
        if h is None or shape_str != "":
            break          # reversed / SL rows: leave the existing behaviour
        if v is not None and abs(v) > VOL_CEILING:
            stable.discard(h)
        else:
            break
    return stable, limit


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=sorted(METHODS), default="quadratic",
                    help="semi-Lagrangian method line (default: quadratic)")
    ap.add_argument("--allow-partial", action="store_true",
                    help="publish even when studies are missing or an arm fits "
                         "no order at all (overwrites the durable docs/ record)")
    args = ap.parse_args(argv)
    m = METHODS[args.method]
    outdir, prefix = paths.tables_dir(m["theme"]), m["prefix"]
    # Which column separates the arms of a study. The semi-Lagrangian lines are
    # separated by their reconstruction; the Eulerian SDPLS line by its method
    # label (see METHODS["sdpls"]).
    groupby = m.get("groupby", "reconstruction")
    arm_col = (groupby != "reconstruction")

    conv_rows, order_rows = [], []
    missing, stale = [], []
    for study, case, mesh in m["studies"]:
        path = os.path.join(REPO, "studies", study, f"{study}_errors.csv")
        if not os.path.isfile(path):
            print(f"[convtable] no {path}; skip {case}/{mesh}")
            missing.append(study)
            continue
        rows = _rows(path)
        foreign = _ladder_complaints(study, rows)
        if foreign:
            stale.append((study, foreign))
        rows, reversed_rows = _drop_reversed(rows)
        if reversed_rows:
            print(f"[convtable] {study}: EXCLUDED {len(reversed_rows)} reversed "
                  f"(oscillation=on) row(s) -- a reversed t=T error is a "
                  f"time-reversal symmetry check, not an accuracy measurement")
        if not rows:
            print(f"[convtable] {study}: every row was reversed; no order fitted")
            continue
        rows, n_dup, dup_conflicts = _collapse_duplicates(rows, groupby)
        if n_dup:
            print(f"[convtable] {study}: collapsed {n_dup} exactly-duplicated "
                  f"row(s) (an inert axis crossed against an arm re-runs it)")
        for key in dup_conflicts:
            print(f"[convtable] {study}: WARNING two DIFFERENT results share "
                  f"arm/h/CFL {key} -- the method label is hiding a real "
                  f"degree of freedom; both rows kept")
        # per (reconstruction, cfl) group + its stable prefix
        groups = {}
        for r in rows:
            key = (r.get(groupby, ""), r.get("cfl", ""))
            groups.setdefault(key, []).append(r)
        stable = {key: _stable_hset(grp) for key, grp in groups.items()}
        # tidy per-resolution rows (all metrics + health columns + stable flag)
        for r in rows:
            key = (r.get(groupby, ""), r.get("cfl", ""))
            stable_h = stable.get(key, (set(), None))[0]
            out = {"case": case, "mesh": mesh,
                   "study": study,
                   "arm": r.get(groupby, ""),
                   "reconstruction": r.get("reconstruction", ""),
                   "cfl": r.get("cfl", ""),
                   "h": r.get("h", ""),
                   "maxCellSize": r.get("maxCellSize", ""),
                   "stable": int(_f(r.get("h")) in stable_h)}
            for col, _oc, _lbl in METRICS:
                out[col] = r.get(col, "")
            for col in HEALTH:
                out[col] = r.get(col, "")
            conv_rows.append(out)
        # per group -> LSQ order of each error metric over the STABLE prefix only
        for (rec, cfl), grp in sorted(groups.items()):
            stable_h, limit = stable[(rec, cfl)]
            sgrp = [r for r in grp if _f(r.get("h")) in stable_h]
            hs = [_f(r.get("h")) for r in sgrp if _f(r.get("h"))]
            row = {"case": case, "mesh": mesh, "arm": rec,
                   # End time distinguishes otherwise identical study rows: the
                   # T=2 and T=8 reversed-vortex studies share case AND mesh, so
                   # without this the reader cannot tell them apart.
                   "T": (grp[0].get("T", "") if grp else ""),
                   # ...and T is NOT enough once two studies run the same arm on
                   # the same case at the same T. sdplsBetaSweep's beta = 1.0
                   # rung is deliberately identical in configuration to
                   # benchVortexEulerT2's beta/simpleImp arm (it is the
                   # cross-check), and method_label suppresses the default beta,
                   # so both label identically. Their ladders differ (5 rungs to
                   # N=512 vs 7 rungs to N=256), hence so do their fitted
                   # orders, and the published table would have shown two rows
                   # with the same visible key and different numbers.
                   "study": study,
                   "reconstruction": rec, "cfl": cfl,
                   "nLevels": len(grp), "stableLevels": len(sgrp),
                   "hMin": f"{min(hs):.6g}" if hs else "",
                   "hMax": f"{max(hs):.6g}" if hs else "",
                   "hLimit": f"{limit:.6g}" if limit is not None else ""}
            for col, oc, _lbl in METRICS:
                p = _order([(_f(r.get("h")), _f(r.get(col))) for r in sgrp])
                row[oc] = f"{p:.3f}" if p is not None else ""
            order_rows.append(row)
        n_unstable = sum(1 for k in groups if stable[k][1] is not None)
        print(f"[convtable] {study}: {len(rows)} rows, {len(groups)} group(s), "
              f"{n_unstable} with a stability limit")

    if not conv_rows:
        print(f"[convtable] no {args.method} studies found; nothing written")
        return 1

    # PUBLISHED-TABLE GUARD. These outputs live under docs/**/data and ARE the
    # durable record: the article and deck build from them without re-running
    # anything. Regenerating them from a partial `studies/` tree silently
    # replaces good published orders with rows of "--" -- measured, not
    # hypothetical: a local run with an incomplete LSL study set rewrote a
    # 7-level 2D vortex row as `1$^\dagger$ & -- & -- & --`. Refuse unless the
    # whole configured study set is present, or the caller says it means it.
    # Stale data is worse than missing data: it publishes silently. Refuse hard,
    # and do NOT let --allow-partial through -- "publish anyway" is a statement
    # about COVERAGE, never a licence to publish rungs from a superseded ladder.
    if stale:
        for study, foreign in stale:
            print(f"[convtable] REFUSING: studies/{study} contains resolution(s) "
                  f"{foreign} that config/{study}.yaml does not declare. That "
                  f"tree is left over from a previous ladder; its rungs and the "
                  f"current ones would be published as one table.")
        print("[convtable] Move or delete the stale study tree(s) and re-pull.")
        return 3
    if missing and not args.allow_partial:
        print(f"[convtable] REFUSING to overwrite the published {args.method} "
              f"tables: {len(missing)} of {len(m['studies'])} studies are "
              f"absent from studies/ ({', '.join(missing)}).")
        print("[convtable] These files are the durable record that docs/ builds "
              "from. Pull the missing studies (make pull-study STUDY=<name>) or "
              "re-run them, then retry. Pass --allow-partial to override.")
        return 2
    # Even with a complete set, never publish an arm whose fit collapsed to
    # nothing -- that is the same failure wearing a different hat.
    empty_arms = [r for r in order_rows
                  if all(not r.get(oc) for _c, oc, _l in METRICS)]
    if empty_arms and not args.allow_partial:
        print(f"[convtable] REFUSING: {len(empty_arms)} arm(s) fitted NO order "
              f"at all (e.g. {empty_arms[0].get('study','?')} / "
              f"{empty_arms[0].get('arm','?')} at CFL {empty_arms[0].get('cfl','?')}) "
              f"-- the study data is incomplete. Pass --allow-partial to override.")
        return 2

    os.makedirs(outdir, exist_ok=True)
    conv_path = os.path.join(outdir, f"{prefix}_convergence.csv")
    with open(conv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(conv_rows[0].keys()))
        w.writeheader(); w.writerows(conv_rows)
    order_path = os.path.join(outdir, f"{prefix}_convergence_orders.csv")
    with open(order_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(order_rows[0].keys()))
        w.writeheader(); w.writerows(order_rows)
    print(f"[convtable] wrote {conv_path} ({len(conv_rows)} rows)")
    print(f"[convtable] wrote {order_path} ({len(order_rows)} rows)")

    # Booktabs LaTeX tabulars so the article tables stay a single source
    # (regenerated here; \input'd by the .tex). Body only -- the article supplies
    # the surrounding table float + caption.
    pretty = {"2Dvortex": "2D reversed vortex", "3Dshear": "3D shear",
              "3Ddeformation": "3D deformation"}

    # Every row of a given table shares one discretization (that is enforced by
    # workflow/scripts/check_discretization.py), so repeating it in all 20 arm
    # cells only costs width. Strip it into a footnote instead.
    _disc = sorted({m.group(0) for r in order_rows
                    for m in [re.search(r"\+div:.*$", r["arm"])] if m})
    def _arm(a):
        return re.sub(r"\+div:.*$", "", a).replace("euler+", "").replace("euler", "no source")

    # "Levels" cell: stable levels used for the fit; a dagger flags a series
    # trimmed at a stability limit (orders are over the stable envelope only).
    def _levels(r):
        s, n = r.get("stableLevels", r["nLevels"]), r["nLevels"]
        return f"{s}$^\\dagger$" if r.get("hLimit") else str(s)
    any_limit = any(r.get("hLimit") for r in order_rows)
    # An extra leading column when the arms are not the reconstruction (SDPLS).
    arm_hdr = "Arm & " if arm_col else ""
    arm_spec = "l" if arm_col else ""
    # The T column disambiguates otherwise identical rows (the T=2 and T=8
    # reversed-vortex studies share case AND mesh). Shown only alongside the arm
    # column, i.e. only for the SDPLS suite -- the semi-Lagrangian tables are
    # published artifacts of other method lines and must stay byte-identical.
    t_hdr = "$T$ & " if arm_col else ""
    t_spec = "c" if arm_col else ""
    # Study column, same reasoning as T and under the same arm_col guard: two
    # studies can contribute the SAME arm at the same case/mesh/T/CFL (the
    # sdplsBetaSweep beta = 1.0 cross-check against benchVortexEulerT2), and
    # without this the two rows are visually identical while their fitted orders
    # differ because their ladders do.
    s_hdr = "Study & " if arm_col else ""
    s_spec = "l" if arm_col else ""
    note = ("\\multicolumn{%d}{l}{\\footnotesize $^\\dagger$stable envelope only; " % (7 + 2*int(arm_col))
            + "finer resolutions destabilized (excluded from the fit).}\\\\\n"
            if any_limit else "")

    # 1. primary table (layout: shape, volume, gradient; Levels = stable levels)
    tex_path = os.path.join(outdir, "convergence_orders.tex")
    with open(tex_path, "w") as fh:
        fh.write("% Auto-generated by workflow/scripts/make_convergence_table.py -- do not edit.\n")
        fh.write("\\begin{tabular}{ll" + arm_spec + s_spec + t_spec + "ccccc}\n\\toprule\n")
        fh.write("Case & Mesh & " + arm_hdr + s_hdr + t_hdr + "CFL & Levels & shape & volume & "
                 "$\\bigl\\|\\,|\\nabla\\psi|-1\\bigr\\|$ \\\\\n\\midrule\n")
        for r in order_rows:
            fh.write(" & ".join([
                pretty.get(r["case"], r["case"]), r["mesh"]]
                + ([method_label.latex_escape(_arm(r["arm"]))] if arm_col else [])
                + ([method_label.latex_escape(r.get("study", ""))] if arm_col else [])
                + ([str(r.get("T", ""))] if arm_col else [])
                + [str(r["cfl"]),
                _levels(r),
                r["shapeOrder"] or "--", r["volumeOrder"] or "--",
                r["gradientOrder"] or "--"]) + " \\\\\n")
        fh.write("\\midrule\n" + note if note else "")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(f"[convtable] wrote {tex_path}")

    # 2. extended table: ALL convergent metrics incl. band + half-time variants
    ncol = 3 + 3*int(arm_col) + len(METRICS)
    ext_note = ("\\multicolumn{%d}{l}{\\footnotesize $^\\dagger$orders over the stable "
                "envelope only; finer resolutions destabilized.}\\\\\n" % ncol) if any_limit else ""
    ext_path = os.path.join(outdir, "convergence_orders_extended.tex")
    with open(ext_path, "w") as fh:
        fh.write("% Auto-generated by workflow/scripts/make_convergence_table.py -- do not edit.\n")
        fh.write("\\begin{tabular}{ll" + arm_spec + s_spec + t_spec + "c" + "c" * len(METRICS) + "}\n\\toprule\n")
        fh.write("Case & Mesh & " + arm_hdr + s_hdr + t_hdr + "CFL & "
                 + " & ".join(lbl for _c, _o, lbl in METRICS) + " \\\\\n\\midrule\n")
        for r in order_rows:
            mesh_cell = r["mesh"] + ("$^\\dagger$" if r.get("hLimit") else "")
            fh.write(" & ".join(
                [pretty.get(r["case"], r["case"]), mesh_cell]
                + ([method_label.latex_escape(_arm(r["arm"]))] if arm_col else [])
                + ([method_label.latex_escape(r.get("study", ""))] if arm_col else [])
                + ([str(r.get("T", ""))] if arm_col else [])
                + [str(r["cfl"])]
                + [(r[oc] or "--") for _c, oc, _lbl in METRICS]) + " \\\\\n")
        fh.write("\\midrule\n" + ext_note if ext_note else "")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(f"[convtable] wrote {ext_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
