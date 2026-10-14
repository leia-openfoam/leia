#!/usr/bin/env python3
"""The static source test (STATUS 11.27; config/staticSource/): metrics, null checks, verdicts.

Every case of cases/2DstaticCircle integrates ddt(psi) == S(psi) at u = 0 (leiaLevelSetFoam) from a
psi0 whose zero set is the exact circle (centre CIRCLE_X0/Y0, radius CIRCLE_R; system/initFields).
This script reads the reconstructed fields at the write times and the per-step CSVs, and measures
everything against the EXACT geometry and against the sampled signed distance d (the target):

  zcL1, zcL2       error of the linearly interpolated zero crossing of psi on every centre-to-centre
                   edge against the exact crossing of the circle, in units of h
  shift1, shift2   the crossing shift against tau = 0, in h (edges that cross at both times)
  spuriousNear/Far edges where psi changes sign and d does not (|d| <= 3h / farther); missing: the
                   reverse
  shapeEx1         sum V |alpha - alpha_exact| / V_exact (alpha_exact: the closed-form area of the
                   disk in each cell); shapeEx2: RMS of alpha - alpha_exact over the cut cells
  shapeEnd1        the same against alphaEnd = alpha_DA(d); shapeChg1 against alpha(0)
  VDA, VDAend, Vex the Detrixhe-Aslam volume, that of d, pi R^2 * depth; volEx, volDh, volChg
  dist3/6/G_1/_2   L1 and L2 of (psi - d)/h in |d| <= 3h, 6h, everywhere
  eik3/6/G_1/_2    L1 and L2 of |grad psi| - 1 with the centred difference (interior cells)
  saw_1/_2         L1 and L2 of D2(D2(psi - psi0))/64 in |d| <= 3h (D2 the 5-point difference)
  ratio            mean of psi/d in h/4 < |d| < h
  centreChange, wallChange   L2 of psi - psi0 in r < 3h of the centre / in the wall cells, over
                   the L2 of psi0 there
Norms are L1 and L2 only (never L_inf). Cell order: i + N j (blockMesh, x fastest).

Commands:
  selftest
  check-init <case>... [--out <dir>]   (P0.2; with --out also staticSource_init.csv)
  analyse --studies <study dir>... --out <dir> [--jobs N]
  null --studies <study dir>... --out <dir>
  seam --serial <study dir>... --parallel <study dir>... --out <dir>
  pace --studies <study dir>...
  curate --metrics <dir> [--theme gradient-controlled-level-set] [--stamp <libleiaCore stamp>]
numpy and the standard library only; runs on the laptop and under snakemake's interpreter on the
cluster.
"""
import argparse
import concurrent.futures as cf
import csv
import glob
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEPTH = 0.1          # z extent of cases/2DstaticCircle/system/blockMeshDict.template
TAUS = (0, 10, 20, 30, 40, 50, 60, 70, 80)
HIST_TAUS = (0, 0.1, 0.2, 0.5, 1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 50, 60, 70, 80)
SKIP_COLS = {"ELAPSED_CPU_TIME", "ELAPSED_CLOCK_TIME"}
LINEAR_LAWS = ("linearQ", "linearZ", "saturatedLinearZ", "twoThirdsZReg")
FLAT_LAWS = ("cubicQ", "cubicZ", "softWall")


# --------------------------------------------------------------------------- reading
def read_scalar_field(path, n=None):
    text = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", text)
    if m:
        a = m.end()
        b = text.index(")", a)
        v = np.array(text[a:b].split(), dtype=float)
        if v.size != int(m.group(1)):
            raise ValueError(f"{path}: {v.size} values, header says {m.group(1)}")
        return v
    m = re.search(r"internalField\s+uniform\s+([-+0-9.eEinfa]+)\s*;", text)
    if m and n:
        return np.full(n, float(m.group(1)))
    raise ValueError(f"cannot read the internalField of {path}")


def read_csv(path):
    """Rows of a solver or function-object CSV. A function-object CSV (gradPsiError.csv) is
    APPENDED to when a case is re-run without its preprocess rule (MEASURED 2026-10-14: a solve
    stopped by the time limit and re-run left the stopped run's rows in front of the new ones), so
    the last row of every TIME is kept, which is the row of the latest run."""
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if rows and any(v in (None, "") for v in rows[-1].values()):
        rows = rows[:-1]          # an incomplete last row of a run that stopped
    if rows and "TIME" in rows[0]:
        rows = [r for r in rows if r["TIME"] != "TIME"]      # a repeated header line
        last = {}
        for r in rows:
            last[r["TIME"]] = r
        if len(last) != len(rows):
            rows = sorted(last.values(), key=lambda r: float(r["TIME"]))
    return rows


def time_dirs(case):
    out = []
    for name in os.listdir(case):
        if re.fullmatch(r"[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?", name) and os.path.isdir(os.path.join(case, name)):
            out.append((float(name), name))
    return sorted(out)


def read_kv(path):
    d = {}
    if os.path.isfile(path):
        for line in open(path):
            if "=" in line:
                k, v = line.rstrip("\n").split("=", 1)
                d[k.strip()] = v.strip()
    return d


def case_info(case):
    p = json.load(open(os.path.join(case, "case_params.json")))
    tok = p["tokens"]
    N = int(float(tok["N_CELLS"]))
    info = dict(case=case, name=os.path.basename(case), study=os.path.basename(os.path.dirname(case)),
                N=N, h=1.0 / N, R=float(tok["CIRCLE_R"]), cx=float(tok["CIRCLE_X0"]),
                cy=float(tok["CIRCLE_Y0"]), init=tok["INIT_FIELD"], source=tok["SDPLS_SOURCE"],
                law=tok.get("GC_LAW", "none"), weight=tok.get("GC_STRAIN_WEIGHT", "none"),
                moll=tok.get("MOLLIFIER", "none"), endTime=float(tok["END_TIME"]),
                git=p.get("gitCommit", ""), mode=p.get("mode", ""), np=p.get("np", ""), tokens=tok)
    info["arm"] = arm_name(info)
    return info


def arm_name(c):
    if c["source"] == "gradientControl":
        return f"gc-{c['law']}-{c['weight']}-{c['moll']}"
    return f"{c['source']}-{c['moll']}"


def case_dirs(study_dirs):
    out = []
    for s in study_dirs:
        out += sorted(d for d in glob.glob(os.path.join(s, "2DstaticCircle_*"))
                      if os.path.isfile(os.path.join(d, "case_params.json")))
    return out


def log_state(case):
    log = os.path.join(case, "log.leiaLevelSetFoam")
    try:
        r = subprocess.run(["bash", os.path.join(HERE, "foam_log_state.sh"), log],
                           capture_output=True, text=True, timeout=120)
        return (r.stdout.strip().splitlines() or [""])[0]
    except Exception as e:  # the classifier is a report, never a failure of this script
        return f"UNKNOWN ({e})"


# --------------------------------------------------------------------------- exact geometry
def _F(u, R):
    """Antiderivative of sqrt(R^2 - u^2)."""
    u = max(-R, min(R, u))
    return 0.5 * (u * math.sqrt(max(R * R - u * u, 0.0)) + R * R * math.asin(u / R))


def _breakpoints(x0, x1, y0, y1, cx, cy, R):
    xs = [x0, x1, cx - R, cx + R]
    for yk in (y0, y1):
        dy = yk - cy
        if abs(dy) < R:
            w = math.sqrt(R * R - dy * dy)
            xs += [cx - w, cx + w]
    return sorted(x for x in xs if x0 <= x <= x1)


def disk_rect_area(x0, x1, y0, y1, cx, cy, R):
    """Exact area of the disk (cx, cy, R) inside the rectangle [x0, x1] x [y0, y1]."""
    xs = _breakpoints(x0, x1, y0, y1, cx, cy, R)
    area = 0.0
    for a, b in zip(xs[:-1], xs[1:]):
        if b <= a:
            continue
        u = 0.5 * (a + b) - cx
        if abs(u) >= R:
            continue
        s = math.sqrt(R * R - u * u)
        top_c = cy + s < y1
        bot_c = cy - s > y0
        top = cy + s if top_c else y1
        bot = cy - s if bot_c else y0
        if top <= bot:
            continue
        Is = _F(b - cx, R) - _F(a - cx, R)
        it = cy * (b - a) + Is if top_c else y1 * (b - a)
        ib = cy * (b - a) - Is if bot_c else y0 * (b - a)
        area += it - ib
    return area


def _chord_length(x, y0, y1, cx, cy, R):
    u = x - cx
    if abs(u) >= R:
        return 0.0
    s = math.sqrt(R * R - u * u)
    return max(0.0, min(y1, cy + s) - max(y0, cy - s))


def disk_rect_area_gauss(x0, x1, y0, y1, cx, cy, R, npts=64):
    """The same area by Gauss-Legendre quadrature between the same breakpoints (selftest), in
    the variable phi with x = cx + R sin(phi): the chord length has a square-root edge at
    x = cx +- R, which the substitution removes, so the quadrature converges exponentially."""
    xg, wg = np.polynomial.legendre.leggauss(npts)
    xs = _breakpoints(x0, x1, y0, y1, cx, cy, R)
    area = 0.0
    for a, b in zip(xs[:-1], xs[1:]):
        a, b = max(a, cx - R), min(b, cx + R)
        if b <= a:
            continue
        pa, pb = math.asin((a - cx) / R), math.asin((b - cx) / R)
        pm, pr = 0.5 * (pa + pb), 0.5 * (pb - pa)
        area += pr * sum(w * _chord_length(cx + R * math.sin(pm + pr * t), y0, y1, cx, cy, R)
                         * R * math.cos(pm + pr * t) for t, w in zip(xg, wg))
    return area


class Geometry:
    """Cell centres, d, the exact alpha and the exact edge crossings of one mesh."""

    _cache = {}

    @classmethod
    def get(cls, N, cx, cy, R):
        key = (N, cx, cy, R)
        if key not in cls._cache:
            cls._cache[key] = cls(N, cx, cy, R)
        return cls._cache[key]

    def __init__(self, N, cx, cy, R):
        self.N, self.cx, self.cy, self.R = N, cx, cy, R
        h = self.h = 1.0 / N
        x = (np.arange(N) + 0.5) * h
        self.X, self.Y = np.meshgrid(x, x)
        self.r = np.sqrt((self.X - cx) ** 2 + (self.Y - cy) ** 2)
        self.d = self.r - R
        self.theta = np.arctan2(self.Y - cy, self.X - cx)
        ae = (self.d < 0).astype(float)
        for j, i in zip(*np.nonzero(np.abs(self.d) < h)):
            ae[j, i] = disk_rect_area(i * h, (i + 1) * h, j * h, (j + 1) * h, cx, cy, R) / (h * h)
        self.alpha_exact = ae
        self.Vcell = h * h * DEPTH
        self.Vex = math.pi * R * R * DEPTH
        self.edges = {}
        d = self.d
        for axis in ("x", "y"):
            if axis == "x":
                da, db, xa, ya = d[:, :-1], d[:, 1:], self.X[:, :-1], self.Y[:, :-1]
                w = np.sqrt(np.maximum(R * R - (ya - cy) ** 2, 0.0))
                r1, r2 = (cx - w - xa) / h, (cx + w - xa) / h
            else:
                da, db, xa, ya = d[:-1, :], d[1:, :], self.X[:-1, :], self.Y[:-1, :]
                w = np.sqrt(np.maximum(R * R - (xa - cx) ** 2, 0.0))
                r1, r2 = (cy - w - ya) / h, (cy + w - ya) / h
            cross = da * db < 0
            s = np.where(cross & (r1 >= 0) & (r1 <= 1), r1, np.where(cross & (r2 >= 0) & (r2 <= 1), r2, np.nan))
            self.edges[axis] = (cross, s, np.minimum(np.abs(da), np.abs(db)))
        interior = np.zeros((N, N), bool)
        interior[1:-1, 1:-1] = True
        self.interior = interior
        inner2 = np.zeros((N, N), bool)
        inner2[2:-2, 2:-2] = True
        self.inner2 = inner2
        self.wall = ~interior


def analytic_psi0(case, geo, init):
    """psi0 of the variant from the coefficients in the case's system/initFields/<init>."""
    text = open(os.path.join(case, "system", "initFields", init)).read()
    coef = dict(re.findall(r'"(c[ABEM])\s*=\s*([^"]+)"', text))
    env = {"sqr": lambda v: v * v, "rad": geo.R, "hc": geo.h, "__builtins__": {}}
    c = {k: float(eval(v, env)) for k, v in coef.items()}
    sd = geo.d
    return c["cA"] * (sd + c["cB"] * sd * sd / (2 * geo.R)) * (1 + c["cE"] * np.cos(c["cM"] * geo.theta))


# --------------------------------------------------------------------------- the zero-set ODE
def law_G(tok, source, law):
    """G(q) at u = 0 (sigma = 0) for the source and law of a case (gradientControlLaws.H)."""
    f = lambda key, dflt: float(tok.get(key, dflt))
    mu = f("GC_MU", 1)
    if source == "beta":
        beta = f("SDPLS_BETA", 1)
        return lambda q: beta - q
    if source != "gradientControl" or law == "none":
        return lambda q: 0.0 * q
    zeta = lambda q: 0.5 * (q * q - 1)
    if law == "linearQ":
        return lambda q: -mu * (q - 1)
    if law == "linearZ":
        return lambda q: -mu * zeta(q)
    if law == "cubicQ":
        return lambda q: -mu * (q - 1) ** 3
    if law == "cubicZ":
        return lambda q: -mu * zeta(q) ** 3
    if law == "twoThirdsZReg":
        eps = f("GC_EPS", 0.02)
        return lambda q: -mu * zeta(q) / (zeta(q) ** 2 + eps * eps) ** (1.0 / 6.0)
    if law == "saturatedLinearZ":
        c = f("GC_SAT_C", 1)
        return lambda q: -mu * zeta(q) / (1 + c * np.abs(zeta(q)))
    if law == "softWall":
        K = f("SW_C_KAPPA", 1.25) * f("SW_EPS_D", 0.0)      # sigma = 0 at u = 0
        gam, ds, p = f("SW_GAMMA", 1.4722194895832), f("SW_DELTA_S", 0.08), int(f("SW_P", 5))
        return lambda q: -K * np.tanh(gam * ((q - 1) / ds) ** p)
    raise ValueError(f"unknown law {law}")


def ode_q(G, q0, taus, dtau=1e-3):
    """RK4 of dq/dtau = q G(q) from q0 (scalar or array); values at the sorted taus."""
    q = np.array(q0, dtype=float)
    out, t, k = [], 0.0, 0
    rhs = lambda v: v * G(v)
    for target in taus:
        while t < target - 1e-12:
            st = min(dtau, target - t)
            k1 = rhs(q)
            k2 = rhs(q + 0.5 * st * k1)
            k3 = rhs(q + 0.5 * st * k2)
            k4 = rhs(q + st * k3)
            q = q + st / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += st
        out.append(q.copy())
    return out


def q0_of(init, R, theta=None):
    if init == "alg":
        return 2 * R
    if init in ("modM2", "modM4"):
        m = 2 if init == "modM2" else 4
        th = np.linspace(-math.pi, math.pi, 720, endpoint=False) if theta is None else theta
        return 1 + 0.1 * np.cos(m * th)
    return 1.0


# --------------------------------------------------------------------------- metrics of one field
def _norms(v):
    v = np.asarray(v, dtype=float)
    if v.size == 0:
        return float("nan"), float("nan")
    return float(np.mean(np.abs(v))), float(np.sqrt(np.mean(v * v)))


def _d2(f):
    g = np.zeros_like(f)
    g[1:-1, 1:-1] = f[2:, 1:-1] + f[:-2, 1:-1] + f[1:-1, 2:] + f[1:-1, :-2] - 4 * f[1:-1, 1:-1]
    return g


def field_metrics(geo, P, P0, A, A0, Aend):
    N, h, d = geo.N, geo.h, geo.d
    m = {}
    errs, shifts = [], []
    sp_near = sp_far = miss = 0
    for axis in ("x", "y"):
        cross_d, s_ex, dmin = geo.edges[axis]
        if axis == "x":
            pa, pb, qa, qb = P[:, :-1], P[:, 1:], P0[:, :-1], P0[:, 1:]
        else:
            pa, pb, qa, qb = P[:-1, :], P[1:, :], P0[:-1, :], P0[1:, :]
        cross_p = pa * pb < 0
        both = cross_p & cross_d
        errs.append(pa[both] / (pa[both] - pb[both]) - s_ex[both])
        b2 = cross_p & (qa * qb < 0)
        shifts.append(pa[b2] / (pa[b2] - pb[b2]) - qa[b2] / (qa[b2] - qb[b2]))
        spur = cross_p & ~cross_d
        sp_near += int(np.sum(spur & (dmin <= 3 * h)))
        sp_far += int(np.sum(spur & (dmin > 3 * h)))
        miss += int(np.sum(cross_d & ~cross_p))
    m["zcL1"], m["zcL2"] = _norms(np.concatenate(errs))
    m["shift1"], m["shift2"] = _norms(np.concatenate(shifts))
    m["spuriousNear"], m["spuriousFar"], m["missing"] = sp_near, sp_far, miss
    V, Vex = geo.Vcell, geo.Vex
    dA = A - geo.alpha_exact
    m["shapeEx1"] = float(np.sum(np.abs(dA)) * V / Vex)
    cut = ((geo.alpha_exact > 0) & (geo.alpha_exact < 1)) | ((A > 0) & (A < 1))
    m["shapeEx2"] = float(np.sqrt(np.mean(dA[cut] ** 2))) if cut.any() else float("nan")
    m["shapeEnd1"] = float(np.sum(np.abs(A - Aend)) * V / Vex)
    m["shapeChg1"] = float(np.sum(np.abs(A - A0)) * V / Vex)
    m["VDA"] = float(np.sum(A) * V)
    m["VDAend"] = float(np.sum(Aend) * V)
    m["VDA0"] = float(np.sum(A0) * V)
    m["Vex"] = Vex
    m["volEx"] = (m["VDA"] - Vex) / Vex
    m["volDh"] = (m["VDA"] - m["VDAend"]) / Vex
    m["volChg"] = (m["VDA"] - m["VDA0"]) / Vex
    regions = {"3": np.abs(d) <= 3 * h, "6": np.abs(d) <= 6 * h, "G": np.ones_like(d, bool)}
    for k, mask in regions.items():
        m[f"dist{k}_1"], m[f"dist{k}_2"] = _norms(((P - d) / h)[mask])
    gx = (P[1:-1, 2:] - P[1:-1, :-2]) / (2 * h)
    gy = (P[2:, 1:-1] - P[:-2, 1:-1]) / (2 * h)
    e = np.sqrt(gx * gx + gy * gy) - 1
    di = d[1:-1, 1:-1]
    for k, mask in {"3": np.abs(di) <= 3 * h, "6": np.abs(di) <= 6 * h, "G": np.ones_like(di, bool)}.items():
        m[f"eik{k}_1"], m[f"eik{k}_2"] = _norms(e[mask])
    cb = _d2(_d2(P - P0)) / 64.0
    m["saw_1"], m["saw_2"] = _norms(cb[(np.abs(d) <= 3 * h) & geo.inner2])
    band = (np.abs(d) > h / 4) & (np.abs(d) < h)
    m["ratio"] = float(np.mean(P[band] / d[band])) if band.any() else float("nan")
    centre = geo.r < 3 * h
    for k, mask in (("centreChange", centre), ("wallChange", geo.wall)):
        den = float(np.sqrt(np.mean(P0[mask] ** 2))) if mask.any() else float("nan")
        num = float(np.sqrt(np.mean((P - P0)[mask] ** 2))) if mask.any() else float("nan")
        m[k] = num / den if den and den > 0 else float("nan")
    return m


# --------------------------------------------------------------------------- one case
def analyse_case(case):
    c = case_info(case)
    geo = Geometry.get(c["N"], c["cx"], c["cy"], c["R"])
    N = c["N"]
    times = time_dirs(case)
    launch = read_kv(os.path.join(case, ".leia_launch"))
    state = log_state(case)
    summary = dict(study=c["study"], case=c["name"], arm=c["arm"], source=c["source"], law=c["law"],
                   weight=c["weight"], moll=c["moll"], init=c["init"], N=N, mode=c["mode"], np=c["np"],
                   gitCommit=c["git"], libStamp=open(os.path.join(case, "leia.version")).read().strip().replace("\n", " ")
                   if os.path.isfile(os.path.join(case, "leia.version")) else "",
                   state=state.split(" ")[0] if state else "", stateLine=state,
                   wallSeconds=(float(launch["end"]) - float(launch["start"])) if launch.get("end") and launch.get("start") else "",
                   slurmJobId=launch.get("slurmJobId", ""), lastTime=times[-1][0] if times else "",
                   endTime=c["endTime"])
    rows, hist = [], []
    if not times or times[0][0] != 0.0:
        summary["note"] = "no time 0"
        return summary, rows, hist
    rd = lambda t, f: read_scalar_field(os.path.join(case, t, f), N * N).reshape(N, N)
    try:
        P0 = rd(times[0][1], "psi")
        A0 = rd(times[0][1], "alpha")
        Aend = rd("0", "alphaEnd")
    except Exception as e:
        summary["note"] = f"unreadable fields at 0: {e}"
        return summary, rows, hist
    for tv, tn in times:
        try:
            P = rd(tn, "psi")
            A = rd(tn, "alpha")
        except Exception as e:
            summary.setdefault("note", f"unreadable fields at {tn}: {e}")
            continue
        r = dict(study=c["study"], case=c["name"], arm=c["arm"], init=c["init"], N=N, tau=tv)
        r.update(field_metrics(geo, P, P0, A, A0, Aend))
        r["psiBitEqual0"] = int(np.array_equal(P, P0))
        rows.append(r)
    adv = read_csv(os.path.join(case, "leiaLevelSetFoam.csv"))
    gp = read_csv(os.path.join(case, "gradPsiError.csv"))
    G = law_G(c["tokens"], c["source"], c["law"])
    q0 = q0_of(c["init"], c["R"])
    qode = ode_q(G, q0, HIST_TAUS)
    for k, tau in enumerate(HIST_TAUS):
        if tau > c["endTime"] + 1e-12:
            continue
        h = dict(study=c["study"], case=c["name"], arm=c["arm"], init=c["init"], N=N, tau=tau)
        for rows_, cols in ((adv, ("E_VOL_ALPHA_REL", "E_VOL_ALPHA_SIGNED_REL", "E_GEOM_ALPHA_REL")),
                            (gp, ("NARROW_MEAN_MAG_GRAD_PSI", "E_NARROW_L1_GRAD_PSI", "E_NARROW_L2_GRAD_PSI",
                                  "E_L2_GRAD_PSI"))):
            if not rows_:
                continue
            best = min(rows_, key=lambda r: abs(float(r["TIME"]) - tau))
            if abs(float(best["TIME"]) - tau) > 0.5 * geo.h + 1e-12:
                continue
            for col in cols:
                h[col] = best.get(col, "")
        qv = qode[k]
        if np.ndim(qv) == 0:
            h["qODE"] = float(qv)
        else:
            h["qODErms"] = float(np.sqrt(np.mean((qv - 1) ** 2)))
        hist.append(h)
    return summary, rows, hist


# --------------------------------------------------------------------------- orders and verdicts
def lsq_order(hs, vs):
    pts = [(math.log(h), math.log(abs(v))) for h, v in zip(hs, vs)
           if v is not None and isinstance(v, (int, float)) and math.isfinite(v) and v != 0]
    if len(pts) < 2:
        return float("nan")
    x = np.array([p[0] for p in pts])
    y = np.array([p[1] for p in pts])
    return float(np.polyfit(x, y, 1)[0])


ORDER_METRICS = ("zcL1", "zcL2", "shift2", "shapeEx1", "shapeEx2", "volEx", "volDh", "dist3_1", "dist3_2",
                 "eik3_2", "saw_2")


def orders_table(fields):
    groups = {}
    for r in fields:
        groups.setdefault((r["study"], r["arm"], r["init"], r["tau"]), []).append(r)
    out = []
    for (study, arm, init, tau), rs in sorted(groups.items()):
        rs = sorted(rs, key=lambda r: r["N"])
        if len(rs) < 2:
            continue
        hs = [1.0 / r["N"] for r in rs]
        for met in ORDER_METRICS:
            vs = [r.get(met) for r in rs]
            signs = {np.sign(v) for v in vs if isinstance(v, float) and math.isfinite(v) and v != 0}
            pair = []
            for a, b in zip(rs[:-1], rs[1:]):
                va, vb = a.get(met), b.get(met)
                pair.append(math.log(abs(va / vb)) / math.log(b["N"] / a["N"])
                            if va and vb and math.isfinite(va) and math.isfinite(vb) and va != 0 and vb != 0 else float("nan"))
            out.append(dict(study=study, arm=arm, init=init, tau=tau, metric=met,
                            N=" ".join(str(r["N"]) for r in rs),
                            values=" ".join(f"{r.get(met):.4e}" if isinstance(r.get(met), float) else str(r.get(met)) for r in rs),
                            lsqOrder=lsq_order(hs, vs), pairOrders=" ".join(f"{p:.2f}" for p in pair),
                            signChange=int(len(signs) > 1)))
    return out


def _get(fields, **kw):
    return [r for r in fields if all(r.get(k) == v for k, v in kw.items())]


def verdicts(fields, hist, cases, null_rows, init_rows=()):
    V = []

    def add(pid, scope, value, threshold, verdict, note=""):
        V.append(dict(item=pid, scope=scope, value=value, threshold=threshold, verdict=verdict, note=note))

    ctrl = [r for r in fields if r["arm"] == "noSource-none"]
    # P0.1: noSource psi bit-identical at every write
    if ctrl:
        bad = [r for r in ctrl if not r["psiBitEqual0"]]
        add("P0.1", "noSource, every write", f"{len(ctrl) - len(bad)}/{len(ctrl)} bit-identical", "all",
            "PASS" if not bad else "FAIL", "; ".join(f"{r['case']} tau {r['tau']}" for r in bad[:5]))
    else:
        add("P0.1", "noSource", "", "", "NOT EVALUATED", "no control cases")
    # P1.1, P1.3: representation orders at tau = 0 from the control
    for met, pid, lo in (("shapeEx1", "P1.1", 1.5), ("shapeEx2", "P1.1", 1.5), ("zcL2", "P1.1", 1.5)):
        for init in sorted({r["init"] for r in ctrl}):
            rs = sorted(_get(ctrl, init=init, tau=0.0), key=lambda r: r["N"])
            if len(rs) < 3:
                continue
            o = lsq_order([1 / r["N"] for r in rs], [r[met] * (1 / r["N"] if met == "zcL2" else 1) for r in rs])
            add(pid, f"{met} {init} tau 0 (absolute)", f"{o:.2f}", "2 +- 0.3 (FAIL < 1.5)",
                "PASS" if abs(o - 2) <= 0.3 else ("FAIL" if o < lo else "OUT OF BAND"))
    for init in sorted({r["init"] for r in ctrl}):
        rs = sorted(_get(ctrl, init=init, tau=0.0), key=lambda r: r["N"])
        if len(rs) >= 3:
            o = lsq_order([1 / r["N"] for r in rs], [r["volEx"] for r in rs])
            add("P1.1", f"|volEx| {init} tau 0", f"{o:.2f}", "2 +- 0.3", "PASS" if abs(o - 2) <= 0.3 else "FAIL")
    # P1.2: blend volume differences are eps times the algNorm difference
    for N in sorted({r["N"] for r in ctrl}):
        if N < 128:
            continue
        get = lambda init: next((r for r in _get(ctrl, init=init, N=N, tau=0.0)), None)
        sdf, an = get("sdf"), get("algNorm")
        if not (sdf and an):
            continue
        base = an["VDA"] - sdf["VDA"]
        h, R = 1.0 / N, 0.25
        for init, eps in (("blend01", 0.1), ("blendH1", h / R), ("blendH2", (h / R) ** 2)):
            b = get(init)
            if not b or base == 0:
                continue
            ratio = (b["VDA"] - sdf["VDA"]) / (eps * base)
            add("P1.2", f"{init} N {N}", f"{ratio:.3f}", "[0.85, 1.15]", "PASS" if 0.85 <= ratio <= 1.15 else "FAIL")
    # P1.3: |V_DA(alg) - V_DA(sdf)| order 2
    rs_alg = sorted(_get(ctrl, init="alg", tau=0.0), key=lambda r: r["N"])
    rs_sdf = {r["N"]: r for r in _get(ctrl, init="sdf", tau=0.0)}
    pts = [(1 / r["N"], abs(r["VDA"] - rs_sdf[r["N"]]["VDA"]) / r["Vex"]) for r in rs_alg if r["N"] in rs_sdf]
    if len(pts) >= 3:
        o = lsq_order([p[0] for p in pts], [p[1] for p in pts])
        add("P1.3", "|V_DA(alg) - V_DA(sdf)|/V_exact", f"{o:.2f}", "2 +- 0.3", "PASS" if abs(o - 2) <= 0.3 else "FAIL")
    # P1.4: the modulated starts have no derived value; their representation orders are recorded
    for init in ("modM2", "modM4"):
        rs = sorted(_get(ctrl, init=init, tau=0.0), key=lambda r: r["N"])
        if len(rs) >= 3:
            o = {met: lsq_order([1 / r["N"] for r in rs], [r[met] for r in rs]) for met in ("shapeEx1", "shapeEx2")}
            add("P1.4", f"{init} tau 0", f"shapeEx1 order {o['shapeEx1']:.2f}, shapeEx2 order {o['shapeEx2']:.2f}",
                "no derived value", "RECORDED")
    # P0.2 from check-init --out (staticSource_init.csv); P0.3, P6.1, P7.1, P8.1 from null
    if init_rows:
        bad = [r for r in init_rows if r["verdict"] != "PASS"]
        add("P0.2", "every case (check-init)", f"{len(init_rows) - len(bad)}/{len(init_rows)} PASS",
            "L1, L2 <= 1e-14; signs of d; sdf = psiEnd", "PASS" if not bad else "FAIL",
            "; ".join(f"{r['case']}: {r['result']}" for r in bad[:5]))
    else:
        add("P0.2", "every case", "", "", "NOT EVALUATED", "no staticSource_init.csv (check-init --out)")
    for r in null_rows:
        add(r["item"], r["pair"], r["result"], "bit-identical", r["verdict"], r.get("note", ""))
    # P2.1: qband against the zero-set ODE; one row per arm (band: tau 1, 2, 5, 7, 10; none: tau 1, 2)
    gc_arms = sorted({h["arm"] for h in hist if h["arm"].startswith("gc-") or h["arm"].startswith("beta-")})
    for arm in gc_arms:
        taus = (1, 2, 5, 7, 10) if arm.endswith("-band") else (1, 2)
        fails, checked, worst = [], 0, 0.0
        for h in hist:
            if h["arm"] != arm or h["N"] < 128 or h["tau"] not in taus or "qODE" not in h:
                continue
            if h["init"] not in ("sdf", "algNorm", "alg", "blend01", "blendH1", "blendH2"):
                continue
            if h.get("NARROW_MEAN_MAG_GRAD_PSI") in (None, ""):
                continue
            checked += 1
            dev = abs(float(h["NARROW_MEAN_MAG_GRAD_PSI"]) - h["qODE"])
            worst = max(worst, dev)
            if dev > 0.02:
                fails.append(f"{h['init']} N {h['N']} tau {h['tau']:g}: {dev:.3f}")
        if checked:
            add("P2.1", f"{arm}: qband vs dq/dtau = q G(q), N >= 128, tau {','.join(map(str, taus))}",
                f"{checked - len(fails)}/{checked} within 0.02, worst {worst:.2e}", "0.02",
                "PASS" if not fails else "FAIL", "; ".join(fails[:8]))
        else:
            add("P2.1", arm, "", "", "NOT EVALUATED", "no histories at N >= 128")
    # P2.2: modM2, modM4: eikB follows the RMS over theta of qODE - 1 within 0.01 above the noSource
    # floor, tau <= 10. Applied as |eikB - qODErms| <= floor + 0.01 with floor = eikB of noSource on
    # sdf at the same N and tau 0 (the estimator's own error), at the history samples 0 < tau <= 10;
    # every rung, both mollifiers (the header names no restriction). Fixed 2026-10-14 before any
    # cluster data was read.
    floor = {h["N"]: float(h["E_NARROW_L2_GRAD_PSI"]) for h in hist
             if h["arm"] == "noSource-none" and h["init"] == "sdf" and h["tau"] == 0
             and h.get("E_NARROW_L2_GRAD_PSI") not in (None, "")}
    for arm in gc_arms:
        fails, checked, worst = [], 0, -1.0
        for h in hist:
            if h["arm"] != arm or h["init"] not in ("modM2", "modM4") or not (0 < h["tau"] <= 10):
                continue
            if "qODErms" not in h or h.get("E_NARROW_L2_GRAD_PSI") in (None, "") or h["N"] not in floor:
                continue
            checked += 1
            dev = abs(float(h["E_NARROW_L2_GRAD_PSI"]) - h["qODErms"])
            tol = floor[h["N"]] + 0.01
            worst = max(worst, dev - tol)
            if dev > tol:
                fails.append(f"{h['init']} N {h['N']} tau {h['tau']:g}: {dev:.4f} > {tol:.4f}")
        if checked:
            add("P2.2", f"{arm}: modM eikB vs RMS_theta(qODE - 1), 0 < tau <= 10",
                f"{checked - len(fails)}/{checked} within floor + 0.01, worst margin {worst:+.2e}",
                "|eikB - qODErms| <= floor + 0.01", "PASS" if not fails else "FAIL", "; ".join(fails[:8]))
        else:
            add("P2.2", arm, "", "", "NOT EVALUATED", "no modM histories or no noSource floor")
    # P3.1: the DA volume excess of algNorm over sdf
    for arm in sorted({r["arm"] for r in fields if r["arm"].startswith("gc-")}):
        law = arm.split("-")[1]
        for N in sorted({r["N"] for r in fields if r["arm"] == arm}):
            an = {r["tau"]: r for r in _get(fields, arm=arm, init="algNorm", N=N)}
            sd = {r["tau"]: r for r in _get(fields, arm=arm, init="sdf", N=N)}
            if 0.0 not in an or 0.0 not in sd:
                continue
            dv0 = an[0.0]["VDA"] - sd[0.0]["VDA"]
            if dv0 == 0:
                continue
            for tau in (10.0, 40.0, 80.0):
                if tau not in an or tau not in sd:
                    continue
                if max(an[tau]["shift2"], sd[tau]["shift2"]) > 0.02:
                    add("P3.1", f"{arm} N {N} tau {tau:g}", "", "", "NOT EVALUATED", "shift2 > 0.02h")
                    continue
                ratio = (an[tau]["VDA"] - sd[tau]["VDA"]) / dv0
                if law in LINEAR_LAWS and tau == 10.0:
                    add("P3.1", f"{arm} N {N} tau 10", f"{ratio:.3g}", "|ratio| <= 0.05",
                        "PASS" if abs(ratio) <= 0.05 else "FAIL")
                elif law in FLAT_LAWS:
                    add("P3.1", f"{arm} N {N} tau {tau:g}", f"{ratio:.3g}", "|ratio - 1| <= 0.1",
                        "PASS" if abs(ratio - 1) <= 0.1 else "FAIL")
    # P4: band arm on sdf
    for arm in sorted({r["arm"] for r in fields if r["arm"].startswith("gc-") and r["arm"].endswith("-band")}):
        law = arm.split("-")[1]
        rs = {r["N"]: r for r in _get(fields, arm=arm, init="sdf", tau=80.0)}
        r10 = {r["N"]: r for r in _get(fields, arm=arm, init="sdf", tau=10.0)}
        if not rs:
            continue
        Ns = sorted(rs)
        if law in FLAT_LAWS:
            worst = max(rs[N]["shift2"] for N in Ns)
            add("P4.1", f"{arm} sdf tau 80", f"max over rungs {worst:.4f} h", "<= 0.01h (FAIL > 0.05h)",
                "PASS" if worst <= 0.01 else ("FAIL" if worst > 0.05 else "OUT OF BAND"))
        if law in ("linearQ", "linearZ", "saturatedLinearZ"):
            ok10 = all(r10[N]["shift2"] <= 0.02 for N in Ns if N in r10)
            ok80 = all(rs[N]["shift2"] >= 0.05 for N in Ns)
            fin = rs[Ns[-1]]
            o = lsq_order([1 / N for N in Ns], [rs[N]["shift2"] / N for N in Ns])
            falsified = fin["shift2"] < 0.02 and o >= 1.8
            add("P4.2", f"{arm} sdf", f"shift2(10) <= 0.02h: {ok10}; shift2(80) >= 0.05h: {ok80}; "
                f"finest {fin['shift2']:.4f} h, abs order {o:.2f}", "see header",
                "FALSIFIED" if falsified else ("PASS" if ok10 and ok80 else "FAIL"))
    lq = {r["N"]: r for r in _get(fields, arm="gc-linearQ-none-band", init="sdf", tau=80.0)}
    tt = {r["N"]: r for r in _get(fields, arm="gc-twoThirdsZReg-none-band", init="sdf", tau=80.0)}
    common = sorted(set(lq) & set(tt))
    if common:
        ok = all(tt[N]["shift2"] >= lq[N]["shift2"] for N in common)
        add("P4.3", "twoThirdsZReg vs linearQ, band, sdf, tau 80",
            " ".join(f"N{N}: {tt[N]['shift2']:.4f}/{lq[N]['shift2']:.4f}" for N in common), ">=", "PASS" if ok else "FAIL")
    # P5: none arm
    none_arms = sorted({r["arm"] for r in fields if r["arm"].startswith("gc-") and r["arm"].endswith("-none-none")})
    for arm in none_arms:
        law = arm.split("-")[1]
        r10 = _get(fields, arm=arm, tau=10.0)
        if r10:
            changed = all((r["centreChange"] > 0) and (r["wallChange"] > 0) for r in r10)
            add("P5.1", f"{arm} tau 10, all starts and rungs", f"changed in {sum((r['centreChange'] > 0) and (r['wallChange'] > 0) for r in r10)}/{len(r10)}",
                "all", "PASS" if changed else "FAIL")
            far = [r for r in r10 if r["spuriousFar"] > 0]
            if law in ("softWall", "saturatedLinearZ"):
                add("P5.3", f"{arm} tau 10", f"{len(far)} case(s) with far spurious crossings", "none",
                    "PASS" if not far else "FAIL")
            if law == "linearQ":
                sd10 = [r for r in r10 if r["init"] == "sdf"]
                ok = all(r["spuriousFar"] > 0 for r in sd10)
                add("P5.3", f"{arm} sdf tau 10", f"far spurious in {sum(r['spuriousFar'] > 0 for r in sd10)}/{len(sd10)}",
                    "all", "PASS" if ok else "FAIL")
            if law in ("linearQ", "linearZ", "twoThirdsZReg"):
                # predicted to appear; outside the FAIL condition of P5.3, so recorded
                add("P5.3", f"{arm} tau 10, all starts and rungs", f"{len(far)}/{len(r10)} case(s) with far spurious crossings",
                    "predicted: appear", "RECORDED")
    states = [c for c in cases if c.get("state") not in ("COMPLETED",)]
    add("P5.2", "every case", f"{len(cases) - len(states)}/{len(cases)} COMPLETED", "all",
        "PASS" if not states else "FAIL", "; ".join(f"{c['study']}/{c['case']} {c['state']}" for c in states[:8]))
    for arm in none_arms:
        band = arm.replace("-none-none", "-none-band")
        for init in sorted({r["init"] for r in fields if r["arm"] == arm}):
            a = {r["N"]: r for r in _get(fields, arm=arm, init=init, tau=80.0)}
            b = {r["N"]: r for r in _get(fields, arm=band, init=init, tau=80.0)}
            for N in sorted(set(a) & set(b)):
                if N < 128:
                    continue
                ok = a[N]["shift2"] >= 0.5 * b[N]["shift2"]
                add("P5.4", f"{arm} {init} N {N} tau 80", f"{a[N]['shift2']:.4f} vs band {b[N]['shift2']:.4f}",
                    "none >= band/2", "PASS" if ok else "FAIL")
    return V


# --------------------------------------------------------------------------- null and seam
def _fields_equal(ca, cb):
    ta = [t for _, t in time_dirs(ca)]
    tb = [t for _, t in time_dirs(cb)]
    if ta != tb:
        return False, f"write times differ ({len(ta)} vs {len(tb)})"
    for t in ta:
        for f in ("psi", "alpha"):
            pa, pb = os.path.join(ca, t, f), os.path.join(cb, t, f)
            if not (os.path.isfile(pa) and os.path.isfile(pb)):
                return False, f"missing {f} at {t}"
            if not np.array_equal(read_scalar_field(pa), read_scalar_field(pb)):
                return False, f"{f} differs at {t}"
    return True, f"{len(ta)} writes"


def _csv_equal(fa, fb, tol=0.0):
    ra, rb = read_csv(fa), read_csv(fb)
    if len(ra) != len(rb):
        return False, f"rows {len(ra)} vs {len(rb)}"
    cols = [c for c in (ra[0].keys() if ra else []) if c not in SKIP_COLS]
    worst = 0.0
    for a, b in zip(ra, rb):
        for col in cols:
            if a.get(col) == b.get(col):
                continue
            if tol == 0.0:
                return False, f"{col} differs at TIME {a.get('TIME')}"
            try:
                va, vb = float(a[col]), float(b[col])
            except (TypeError, ValueError):
                return False, f"{col} not numeric"
            rel = abs(va - vb) / max(abs(va), abs(vb), 1e-14)
            worst = max(worst, rel)
            if rel > tol:
                return False, f"{col} rel {rel:.2e} at TIME {a.get('TIME')}"
    return True, f"{len(ra)} rows" + (f", worst rel {worst:.1e}" if tol else "")


def null_checks(study_dirs):
    cases = [case_info(c) for c in case_dirs(study_dirs)]
    index = {}
    for c in cases:
        index[(c["source"], c["law"], c["weight"], c["moll"], c["init"], c["N"])] = c
    out = []

    def pair(item, a, b):
        okf, nf = _fields_equal(a["case"], b["case"])
        oka, na = _csv_equal(os.path.join(a["case"], "leiaLevelSetFoam.csv"), os.path.join(b["case"], "leiaLevelSetFoam.csv"))
        okg, ng = _csv_equal(os.path.join(a["case"], "gradPsiError.csv"), os.path.join(b["case"], "gradPsiError.csv"))
        ok = okf and oka and okg
        out.append(dict(item=item, pair=f"{a['study']}/{a['name']} ({a['arm']}, {a['init']}, N {a['N']}) vs "
                        f"{b['study']}/{b['name']} ({b['arm']})", result=f"fields: {nf}; csv: {na}; grad csv: {ng}",
                        verdict="PASS" if ok else "FAIL"))

    for c in cases:
        key = (c["init"], c["N"])
        if c["source"] in ("R", "Rdiv", "RdivStrictSp"):
            ref = index.get(("noSource", "none", "none", "none") + key)
            if ref:
                pair("P7.1", c, ref)
        elif c["source"] == "gradientControl" and c["law"] == "linearQ" and c["weight"] == "full":
            ref = index.get(("gradientControl", "linearQ", "none", c["moll"]) + key)
            if ref:
                pair("P8.1", c, ref)
        elif c["source"] == "beta":
            ref = index.get(("gradientControl", "linearQ", "none", c["moll"]) + key)
            if ref:
                pair("P6.1", c, ref)
        elif c["source"] == "noSource" and c["init"] == "alg":
            ref = index.get(("noSource", "none", "none", "none", "algNorm", c["N"]))
            if ref:   # P0.3: alpha(alg) = alpha(algNorm) at every write
                ta = [t for _, t in time_dirs(c["case"])]
                ok = all(np.array_equal(read_scalar_field(os.path.join(c["case"], t, "alpha")),
                                        read_scalar_field(os.path.join(ref["case"], t, "alpha"))) for t in ta)
                out.append(dict(item="P0.3", pair=f"{c['name']} alg vs {ref['name']} algNorm, N {c['N']}",
                                result=f"alpha at {len(ta)} writes", verdict="PASS" if ok else "FAIL"))
    return out


SEAM_SKIP = ("ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME,L_INF_E_PSI,E_BOUND_ALPHA,MAX_CFL,MAX_ABS_SUM_VOL_FLUXES,"
             "E_MAX_GRAD_PSI,E_NARROW_MAX_GRAD_PSI,MAX_MAG_GRAD_PSI,NARROW_MAX_MAG_GRAD_PSI,MIN_MAG_GRAD_PSI,"
             "NARROW_MIN_MAG_GRAD_PSI,NARROW_MIN_R,NARROW_MAX_R")


def _compare(fa, fb):
    """The committed comparator: relative 1e-10 with an absolute floor of 1e-7 in the column's
    units (the volume columns are m^3 of size 1e-17 to 1e-8, so the floor makes it an absolute
    1e-17 m^3, the round-off of the volume sums). The extremum columns (minima, maxima, L_INF_E_PSI)
    are L_inf-type statistics and are not compared."""
    r = subprocess.run([sys.executable, os.path.join(HERE, "compare_metrics_csv.py"), fa, fb, "--tol", "1e-10",
                        "--floor", "1e-7", "--skip", SEAM_SKIP], capture_output=True, text=True)
    last = (r.stdout.strip().splitlines() or [""])[-1]
    return r.returncode == 0, last


def seam_checks(serial_dirs, parallel_dirs):
    ser = {(c["arm"], c["init"], c["N"]): c for c in map(case_info, case_dirs(serial_dirs))}
    par = {(c["arm"], c["init"], c["N"]): c for c in map(case_info, case_dirs(parallel_dirs))}
    out = []
    for key in sorted(set(ser) & set(par)):
        a, b = ser[key], par[key]
        ok1, n1 = _compare(os.path.join(a["case"], "leiaLevelSetFoam.csv"), os.path.join(b["case"], "leiaLevelSetFoam.csv"))
        ok2, n2 = _compare(os.path.join(a["case"], "gradPsiError.csv"), os.path.join(b["case"], "gradPsiError.csv"))
        worst = 0.0
        for (_, t) in time_dirs(a["case"]):
            pb = os.path.join(b["case"], t, "psi")
            if os.path.isfile(pb):
                da = read_scalar_field(os.path.join(a["case"], t, "psi")) - read_scalar_field(pb)
                worst = max(worst, float(np.sqrt(np.mean(da * da))))
        ok = ok1 and ok2 and worst <= 1e-12
        judged = key[0] == "noSource-none" or key[0].endswith("-band")   # P9.1: noSource and the band arms
        out.append(dict(item="P9.1", pair=f"{key[0]} {key[1]} N {key[2]}: serial {a['name']} vs np {b['np']} {b['name']}",
                        result=f"csv: {n1}; grad csv: {n2}; L2(dpsi) max over writes {worst:.1e}",
                        verdict=("PASS" if ok else "FAIL") if judged else ("REPORTED (within)" if ok else "REPORTED (outside)")))
    return out


# --------------------------------------------------------------------------- output
def write_rows(path, rows):
    if not rows:
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.10g}" if isinstance(v, float) else v) for k, v in r.items()})


def cmd_selftest(_args):
    ok = True
    R, cx, cy = 0.25, 0.5013, 0.4987
    for N in (64, 128, 256, 512):
        g = Geometry(N, cx, cy, R)
        tot = float(np.sum(g.alpha_exact)) * g.h * g.h
        rel = abs(tot - math.pi * R * R) / (math.pi * R * R)
        print(f"selftest: N {N}: sum of exact cell areas / (pi R^2) - 1 = {rel:.2e}")
        ok &= rel < 1e-13
    rng = np.random.default_rng(20261014)
    worst_mid = worst_gl = worst_refl = 0.0
    for k in range(200):
        h = 10 ** rng.uniform(-3, -1)
        ang = rng.uniform(0, 2 * math.pi)
        rr = R + rng.uniform(-1.5, 1.5) * h
        x0 = cx + rr * math.cos(ang) + rng.uniform(-h, 0)
        y0 = cy + rr * math.sin(ang) + rng.uniform(-h, 0)
        A = disk_rect_area(x0, x0 + h, y0, y0 + h, cx, cy, R)
        n = 2000
        t = (np.arange(n) + 0.5) / n * h
        XX, YY = np.meshgrid(x0 + t, y0 + t)
        Amid = float(np.mean((XX - cx) ** 2 + (YY - cy) ** 2 < R * R)) * h * h
        Agl = disk_rect_area_gauss(x0, x0 + h, y0, y0 + h, cx, cy, R)
        worst_mid = max(worst_mid, abs(A - Amid) / (h * h))
        worst_gl = max(worst_gl, abs(A - Agl) / (h * h))
        Arefl = disk_rect_area(2 * cx - (x0 + h), 2 * cx - x0, y0, y0 + h, cx, cy, R)
        worst_refl = max(worst_refl, abs(A - Arefl) / (h * h))
    print(f"selftest: 200 random cells: |A - midpoint 2000^2| <= {worst_mid:.1e} h^2, "
          f"|A - Gauss-Legendre 64| <= {worst_gl:.1e} h^2, |A - A(reflected about x = cx)| <= {worst_refl:.1e} h^2")
    # F(b) - F(a) subtracts two numbers of size R^2 for a result of size h^2: the closed form has a
    # round-off of about 1e-12 h^2 at h = 1e-3. A formula error would be O(1) h^2.
    ok &= worst_mid <= 1e-4 and worst_gl <= 1e-9 and worst_refl <= 1e-9
    G = law_G({"GC_MU": "1"}, "gradientControl", "linearQ")
    taus = (1, 2, 5, 10)
    qs = ode_q(G, 0.5, taus)
    dev = max(abs(float(q) - 1 / (1 + math.exp(-t))) for q, t in zip(qs, taus))
    print(f"selftest: linearQ zero-set ODE from 0.5 against 1/(1 + e^-tau): {dev:.1e}")
    ok &= dev < 1e-10
    print("selftest:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def cmd_check_init(args):
    ok_all = True
    rows = []
    for case in args.cases:
        c = case_info(case)
        geo = Geometry.get(c["N"], c["cx"], c["cy"], c["R"])
        n = c["N"] * c["N"]
        psi0 = read_scalar_field(os.path.join(case, "0", "psi"), n).reshape(c["N"], c["N"])
        ana = analytic_psi0(case, geo, c["init"])
        diff = psi0 - ana
        l1, l2 = float(np.mean(np.abs(diff))), float(np.sqrt(np.mean(diff * diff)))
        sign_ok = bool(np.all(np.sign(psi0) == np.sign(geo.d)))
        msgs = [f"L1 {l1:.1e} L2 {l2:.1e}", f"sign {'ok' if sign_ok else 'DIFFERS'}"]
        ok = l1 <= 1e-14 and l2 <= 1e-14 and sign_ok
        pe = os.path.join(case, "0", "psiEnd")
        if os.path.isfile(pe) and c["init"] == "sdf":
            eq = np.array_equal(read_scalar_field(pe, n), psi0.ravel())
            msgs.append(f"sdf = psiEnd bit for bit: {eq}")
            ok &= eq
        procs = sorted(glob.glob(os.path.join(case, "processor*")))
        if procs:
            miss = [p for p in procs for f in ("psiEnd", "alphaEnd") if not os.path.isfile(os.path.join(p, "0", f))]
            msgs.append(f"psiEnd/alphaEnd in {len(procs)} processor dirs: {'ok' if not miss else 'MISSING'}")
            ok &= not miss
        else:
            ok &= all(os.path.isfile(os.path.join(case, "0", f)) for f in ("psiEnd", "alphaEnd"))
        print(f"check-init {c['study']}/{c['name']} {c['init']} N {c['N']}: {'PASS' if ok else 'FAIL'} ({'; '.join(msgs)})")
        rows.append(dict(item="P0.2", case=f"{c['study']}/{c['name']}", init=c["init"], N=c["N"],
                         result="; ".join(msgs), verdict="PASS" if ok else "FAIL"))
        ok_all &= ok
    if args.out:
        os.makedirs(args.out, exist_ok=True)
        write_rows(os.path.join(args.out, "staticSource_init.csv"), rows)
    print(f"check-init: {sum(r['verdict'] == 'PASS' for r in rows)}/{len(rows)} PASS")
    return 0 if ok_all else 1


def cmd_analyse(args):
    cases = case_dirs(args.studies)
    if not cases:
        print("analyse: no cases")
        return 1
    commits = {case_info(c)["git"] for c in cases}
    if any(g.endswith("-dirty") for g in commits) and not args.allow_dirty:
        print(f"analyse: refusing -dirty commits {sorted(commits)} (use --allow-dirty for a smoke)")
        return 2
    os.makedirs(args.out, exist_ok=True)
    summaries, fields, hist = [], [], []
    with cf.ProcessPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        for s, rows, h in ex.map(analyse_case, cases, chunksize=4):
            summaries.append(s)
            fields += rows
            hist += h
    write_rows(os.path.join(args.out, "staticSource_cases.csv"), summaries)
    write_rows(os.path.join(args.out, "staticSource_fields.csv"), fields)
    write_rows(os.path.join(args.out, "staticSource_histories.csv"), hist)
    rep = [r for r in fields if r["arm"] == "noSource-none" and r["tau"] == 0.0]
    write_rows(os.path.join(args.out, "staticSource_representation.csv"),
               [{k: r[k] for k in ("init", "N", "zcL1", "zcL2", "shapeEx1", "shapeEx2", "VDA", "VDAend", "Vex",
                                   "volEx", "volDh", "dist3_2", "eik3_2")} for r in sorted(rep, key=lambda r: (r["init"], r["N"]))])
    write_rows(os.path.join(args.out, "staticSource_orders.csv"), orders_table(fields))
    nulls = []
    npath = os.path.join(args.out, "staticSource_null.csv")
    if os.path.isfile(npath):
        nulls = read_csv(npath)
    inits = []
    ipath = os.path.join(args.out, "staticSource_init.csv")
    if os.path.isfile(ipath):
        inits = read_csv(ipath)
    V = verdicts(fields, hist, summaries, nulls, inits)
    write_rows(os.path.join(args.out, "staticSource_verdicts.csv"), V)
    print(f"analyse: {len(cases)} cases, {len(fields)} field rows, {len(hist)} history rows, {len(V)} verdict lines; commits {sorted(commits)}")
    for v in V:
        if v["verdict"] not in ("PASS",):
            print(f"  {v['item']} {v['scope']}: {v['verdict']} {v['value']} {v.get('note', '')}")
    return 0


def cmd_null(args):
    os.makedirs(args.out, exist_ok=True)
    rows = null_checks(args.studies)
    write_rows(os.path.join(args.out, "staticSource_null.csv"), rows)
    bad = [r for r in rows if r["verdict"] != "PASS"]
    print(f"null: {len(rows)} pairs, {len(bad)} not PASS")
    for r in bad[:20]:
        print(f"  {r['item']} {r['pair']}: {r['result']}")
    return 0 if not bad else 1


def cmd_seam(args):
    os.makedirs(args.out, exist_ok=True)
    rows = seam_checks(args.serial, args.parallel)
    write_rows(os.path.join(args.out, "staticSource_seam.csv"), rows)
    bad = [r for r in rows if r["verdict"] == "FAIL"]
    rep = [r for r in rows if r["verdict"].startswith("REPORTED (outside")]
    print(f"seam: {len(rows)} pairs, {len(bad)} FAIL (noSource and band arms), {len(rep)} reported outside (none arms)")
    for r in bad[:20]:
        print(f"  {r['pair']}: {r['result']}")
    return 0 if not bad else 1


def cmd_pace(args):
    for case in case_dirs(args.studies):
        c = case_info(case)
        L = read_kv(os.path.join(case, ".leia_launch"))
        m = re.search(r"steps=(\d+)", L.get("state", ""))
        if L.get("start") and L.get("end") and m and int(m.group(1)) > 0:
            sps = (float(L["end"]) - float(L["start"])) / int(m.group(1))
            full = sps * 160 * c["N"]
            print(f"pace {c['study']}/{c['name']} {c['arm']} {c['init']} N {c['N']} np {c['np']}: "
                  f"{sps:.2e} s/step; a full run (160 N = {160 * c['N']} steps) about {full / 60:.1f} min")
    return 0


def cmd_curate(args):
    theme = args.theme
    tables = os.path.join(REPO, "docs", theme, "gcls-level-set-article", "data", "tables")
    names = ("representation", "orders", "null", "verdicts", "cases", "seam")
    os.makedirs(tables, exist_ok=True)
    copied = []
    for n in names:
        src = os.path.join(args.metrics, f"staticSource_{n}.csv")
        if os.path.isfile(src):
            dst = os.path.join(tables, f"staticSource_{n}.csv")
            shutil.copyfile(src, dst)
            copied.append(dst)
    if args.stamp:
        arch = os.path.join(REPO, "docs", theme, "gcls-level-set-article", "data", "archive", args.stamp, "staticSource")
        os.makedirs(arch, exist_ok=True)
        man = []
        for src in sorted(glob.glob(os.path.join(args.metrics, "staticSource_*.csv"))):
            dst = os.path.join(arch, os.path.basename(src))
            shutil.copyfile(src, dst)
            man.append(dict(file=os.path.basename(src), sha256=hashlib.sha256(open(dst, "rb").read()).hexdigest()))
        write_rows(os.path.join(arch, "MANIFEST.csv"), man)
        copied.append(arch)
    print("curate:", *copied, sep="\n  ")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    p = sub.add_parser("check-init")
    p.add_argument("cases", nargs="+")
    p.add_argument("--out", default="")
    p = sub.add_parser("analyse")
    p.add_argument("--studies", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--allow-dirty", action="store_true")
    p = sub.add_parser("null")
    p.add_argument("--studies", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("seam")
    p.add_argument("--serial", nargs="+", required=True)
    p.add_argument("--parallel", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("pace")
    p.add_argument("--studies", nargs="+", required=True)
    p = sub.add_parser("curate")
    p.add_argument("--metrics", required=True)
    p.add_argument("--theme", default="gradient-controlled-level-set")
    p.add_argument("--stamp", default="")
    a = ap.parse_args()
    return {"selftest": cmd_selftest, "check-init": cmd_check_init, "analyse": cmd_analyse, "null": cmd_null,
            "seam": cmd_seam, "pace": cmd_pace, "curate": cmd_curate}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
