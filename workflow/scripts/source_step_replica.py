#!/usr/bin/env python3
"""A numpy copy of the gradient-control source step at zero velocity (STATUS 11.26). Research instrument.

slGradientControlSource::apply with law linearQ and no strain weight, on the uniform 2D hex mesh of
cases/2Dtranslation (unit square, circle centre (0.25, 0.5), R = 0.15). At u = 0 the SL step is the
identity: the quadratic fit is anchored at the cell value (quadraticWeightedLeastSquares::evaluateRaw)
and is evaluated at the cell centre. Per step, with all cells updated from the same old field:
    q_c = |g_c|,  band: |psi_c| <= bandCells h q_c,  psi_c <- psi_c exp(dt (-mu (q_c - 1))).
The estimators of g_c:
    centred      the centred difference (OpenFOAM leastSquares on a uniform hex mesh; at a
                 zeroGradient wall the face value equals the cell value, at distance h/2)
    godunov      the first-order Godunov upwind difference (Cartesian, for reference only)
    upwindFace   the least-squares gradient over the face neighbours upwind of c (s_c (psi_u - psi_c) < 0)
    upwindPoint  the same over the point neighbours
    upwindQuad   an anchored quadratic least-squares fit over the upwind cells of the 5x5 ring
    ownFit       the gradient of the cell's own anchored quadratic fit (point stencil, weights 1/|d|)
    upwindFit    the upwind-weighted mean of the face neighbours' fit gradients at x_c:
                 g_c = sum_u w_u (g_u + H_u (x_c - x_u)) / sum_u w_u, w_u = s_c (psi_c - psi_u) > 0
                 (qScheme upwindFit of slGradientControlSource)
Metrics at the last time (as source_checkerboard_diag.py "final:"): the zero-crossing shift in h,
rms(psi - d)/h in |d| < 3h, the sawtooth amplitude cb = D2(D2(psi - d))/64, the median psi/d in |d| < h.

Usage:
  source_step_replica.py validate <E1 case dir>     the copy against an OpenFOAM run of gate E1
  source_step_replica.py run <q> <N> [--band B] [--mudt X] [--mut T] [--scale S] [--anchor] [--taper]
  source_step_replica.py table                      the tables of STATUS 11.26 (minutes)
"""
import argparse, os, re, sys, time
import numpy as np

XC, YC, R = 0.25, 0.5, 0.15
FACE = [(1, 0), (-1, 0), (0, 1), (0, -1)]
POINT = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1) if (i, j) != (0, 0)]
RING2 = [(i, j) for i in range(-2, 3) for j in range(-2, 3) if (i, j) != (0, 0)]


def setup(n):
    h = 1.0/n
    x = (np.arange(n) + 0.5)*h
    X, Y = np.meshgrid(x, x)                       # [j][i], x fastest as in OpenFOAM
    return h, np.sqrt((X - XC)**2 + (Y - YC)**2) - R


def shifted(psi, oi, oj):
    p = np.pad(psi, 2, mode="edge")
    n = psi.shape[0]
    return p[2 + oj:2 + oj + n, 2 + oi:2 + oi + n]


def q_centred(psi, h):
    n = psi.shape[0]
    gx = np.empty_like(psi); gy = np.empty_like(psi)
    gx[:, 1:-1] = (psi[:, 2:] - psi[:, :-2])/(2*h); gy[1:-1, :] = (psi[2:, :] - psi[:-2, :])/(2*h)
    gx[:, 0] = (psi[:, 1] - psi[:, 0])/(2*h); gx[:, -1] = (psi[:, -1] - psi[:, -2])/(2*h)
    gy[0, :] = (psi[1, :] - psi[0, :])/(2*h); gy[-1, :] = (psi[-1, :] - psi[-2, :])/(2*h)
    return np.sqrt(gx**2 + gy**2), np.ones(psi.shape, bool)


def q_godunov(psi, h):
    p = np.pad(psi, 1, mode="edge")
    dmx = (p[1:-1, 1:-1] - p[1:-1, :-2])/h; dpx = (p[1:-1, 2:] - p[1:-1, 1:-1])/h
    dmy = (p[1:-1, 1:-1] - p[:-2, 1:-1])/h; dpy = (p[2:, 1:-1] - p[1:-1, 1:-1])/h
    pos = psi > 0
    ax = np.where(pos, np.maximum(np.maximum(dmx, 0)**2, np.minimum(dpx, 0)**2),
                  np.maximum(np.minimum(dmx, 0)**2, np.maximum(dpx, 0)**2))
    ay = np.where(pos, np.maximum(np.maximum(dmy, 0)**2, np.minimum(dpy, 0)**2),
                  np.maximum(np.minimum(dmy, 0)**2, np.maximum(dpy, 0)**2))
    return np.sqrt(ax + ay), np.ones(psi.shape, bool)


def q_upwind_linear(psi, h, offs):
    s = np.sign(psi)
    A11 = np.zeros_like(psi); A12 = np.zeros_like(psi); A22 = np.zeros_like(psi)
    b1 = np.zeros_like(psi); b2 = np.zeros_like(psi); cnt = np.zeros(psi.shape, int)
    best = np.zeros_like(psi)
    for oi, oj in offs:
        dx, dy = oi*h, oj*h; r2 = dx*dx + dy*dy
        dpsi = shifted(psi, oi, oj) - psi
        up = s*dpsi < 0
        wu = np.where(up, 1.0/r2, 0.0)
        A11 += wu*dx*dx; A12 += wu*dx*dy; A22 += wu*dy*dy
        b1 += wu*dx*dpsi; b2 += wu*dy*dpsi; cnt += up
        best = np.maximum(best, np.where(up, np.abs(dpsi)/np.sqrt(r2), 0.0))
    det = A11*A22 - A12*A12
    ok = det > 1e-12*(A11 + A22)**2
    gx = np.where(ok, (A22*b1 - A12*b2)/np.where(ok, det, 1), 0.0)
    gy = np.where(ok, (A11*b2 - A12*b1)/np.where(ok, det, 1), 0.0)
    return np.where(ok, np.sqrt(gx**2 + gy**2), best), cnt > 0


def q_upwind_quad(psi, h, band):
    q, valid = q_upwind_linear(psi, h, POINT)
    idx = np.argwhere(np.abs(psi) <= (band + 0.5)*h*1.1)
    if len(idx) == 0:
        return q, valid
    jj, ii = idx[:, 0], idx[:, 1]
    s = np.sign(psi[jj, ii])
    P = np.pad(psi, 2, mode="edge")
    dps = np.stack([P[jj + 2 + oj, ii + 2 + oi] - psi[jj, ii] for oi, oj in RING2], 1)
    up = s[:, None]*dps < 0
    DX = np.array([oi*h for oi, oj in RING2]); DY = np.array([oj*h for oi, oj in RING2])
    B = np.stack([DX, DY, 0.5*DX*DX, DX*DY, 0.5*DY*DY], 1)/np.array([h, h, h*h, h*h, h*h])
    w = np.where(up, (1.0/(DX**2 + DY**2))[None, :], 0.0)
    M = np.einsum("cm,mk,ml->ckl", w, B, B)
    r = np.einsum("cm,mk,cm->ck", w, B, dps)
    ev = np.linalg.eigvalsh(M + 1e-30*np.eye(5))
    good = (up.sum(1) >= 7) & (ev[:, 0] > 1e-6*ev[:, -1])
    sol = np.zeros((len(jj), 5))
    if good.any():
        sol[good] = np.linalg.solve(M[good], r[good][..., None])[..., 0]
    q[jj[good], ii[good]] = np.sqrt((sol[good, 0]/h)**2 + (sol[good, 1]/h)**2)
    return q, valid


def fit_point(psi, h):
    """leia's quadraticWeightedLeastSquares on the point stencil: rows w_i b_k(d_i), w_i = 1/|d_i|,
    anchored at psi_c; returns gx, gy, hxx, hxy, hyy per cell."""
    DX = np.array([oi*h for oi, oj in POINT]); DY = np.array([oj*h for oi, oj in POINT])
    W = 1.0/np.sqrt(DX**2 + DY**2)
    B = np.stack([DX, DY, 0.5*DX*DX, DX*DY, 0.5*DY*DY], 1)
    Minv = np.linalg.pinv(W[:, None]*B)
    dps = np.stack([shifted(psi, oi, oj) - psi for oi, oj in POINT], -1)
    return np.einsum("km,ijm->ijk", Minv*W[None, :], dps)


def q_own_fit(psi, h):
    cf = fit_point(psi, h)
    return np.sqrt(cf[..., 0]**2 + cf[..., 1]**2), np.ones(psi.shape, bool)


def q_upwind_fit(psi, h):
    cf = fit_point(psi, h)
    s = np.sign(psi)
    acc = np.zeros(psi.shape + (2,)); wsum = np.zeros_like(psi)
    for oi, oj in FACE:
        pu = shifted(psi, oi, oj)
        cu = np.stack([shifted(cf[..., k], oi, oj) for k in range(5)], -1)
        dx, dy = -oi*h, -oj*h                      # x_c - x_u
        w = np.where(s*(pu - psi) < 0, s*(psi - pu), 0.0)
        acc[..., 0] += w*(cu[..., 0] + cu[..., 2]*dx + cu[..., 3]*dy)
        acc[..., 1] += w*(cu[..., 1] + cu[..., 3]*dx + cu[..., 4]*dy)
        wsum += w
    valid = wsum > 0
    g = acc/np.where(valid, wsum, 1)[..., None]
    return np.sqrt(g[..., 0]**2 + g[..., 1]**2), valid


def q_of(name, psi, h, band):
    if name == "centred":
        return q_centred(psi, h)
    if name == "godunov":
        return q_godunov(psi, h)
    if name == "upwindFace":
        return q_upwind_linear(psi, h, FACE)
    if name == "upwindPoint":
        return q_upwind_linear(psi, h, POINT)
    if name == "upwindQuad":
        return q_upwind_quad(psi, h, band)
    if name == "ownFit":
        return q_own_fit(psi, h)
    if name == "upwindFit":
        return q_upwind_fit(psi, h)
    raise SystemExit(f"unknown q estimator {name}")


def taper(r, band):
    """1 for r <= band - 2, a C1 smoothstep to 0 at r = band (gate E2.2's taper, research only)."""
    s = np.clip((band - r)/2.0, 0, 1)
    return s*s*(3 - 2*s)


def step(psi, h, mu, dt, band, qname, anchor=False, use_taper=False):
    q, valid = q_of(qname, psi, h, band)
    inb = valid & (q > 1e-150) & (np.abs(psi) <= band*h*q)
    if anchor:                                     # interface cells keep their values
        sc = np.zeros(psi.shape, bool)
        for oi, oj in FACE:
            sc |= shifted(psi, oi, oj)*psi < 0
        inb &= ~sc
    x = np.clip(dt*(-mu*(q - 1.0)), -30, 30)
    if use_taper:
        x = x*taper(np.abs(psi)/np.maximum(q, 1e-300)/h, band)
    return np.where(inb, psi*np.exp(x), psi)


def run(qname, n, band=3.0, mudt=0.01, mut=80.0, mu=27.0, scale=1.0, anchor=False, use_taper=False,
        marks=()):
    h, d = setup(n)
    psi = scale*d
    dt = mudt/mu
    nsteps = int(round(mut/mudt))
    snaps = {}
    want = {int(round(m/mudt)): m for m in marks}
    for k in range(nsteps + 1):
        if k in want:
            snaps[want[k]] = psi.copy()
        if k < nsteps:
            psi = step(psi, h, mu, dt, band, qname, anchor, use_taper)
    return h, d, psi, snaps


def d2(f):
    g = np.zeros_like(f)
    g[1:-1, 1:-1] = f[2:, 1:-1] + f[:-2, 1:-1] + f[1:-1, 2:] + f[1:-1, :-2] - 4*f[1:-1, 1:-1]
    return g


def metrics(h, d, psi):
    n = psi.shape[0]
    inner = np.zeros((n, n), bool); inner[3:-3, 3:-3] = True
    sh = []
    for a, b in (((slice(3, n-3), slice(3, n-4)), (slice(3, n-3), slice(4, n-3))),
                 ((slice(3, n-4), slice(3, n-3)), (slice(4, n-3), slice(3, n-3)))):
        pa, pb, da, db = psi[a], psi[b], d[a], d[b]
        m = (pa*pb < 0) & (da*db < 0)
        sh.append(np.abs(pa[m]/(pa[m] - pb[m]) - da[m]/(da[m] - db[m])))
    sh = np.concatenate(sh)
    b3 = inner & (np.abs(d) < 3*h)
    fi = inner & (np.abs(d) < h)
    cb = d2(d2(psi - d))/64.0
    return dict(ratio=float(np.median((psi/d)[fi])), shift_max=float(sh.max()),
                shift_rms=float(np.sqrt(np.mean(sh**2))), cbd=float(np.sqrt(np.mean(cb[b3]**2))),
                dev=float(np.sqrt(np.mean(((psi - d)/h)[b3]**2))))


def fmt(m):
    return (f"median psi/d {m['ratio']:.4f} | zero shift/h max {m['shift_max']:.4f} rms {m['shift_rms']:.4f} "
            f"| rms(psi-d)/h {m['dev']:.2e} | cb(psi-d) {m['cbd']:.2e}")


def internal(path):
    t = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", t)
    a = m.end(); b = t.index(")", a)
    return np.array([float(x) for x in t[a:b].split()])


def cmd_validate(case):
    import json
    tok = json.load(open(os.path.join(case, "case_params.json")))["tokens"]
    n = int(tok["N_CELLS"]); mu = float(tok["GC_MU"]); band = float(tok["SL_SOURCE_BAND_CELLS"])
    qname = tok.get("GC_Q_SCHEME", "centred")
    h, d = setup(n)
    psi = internal(os.path.join(case, "0", "psi")).reshape(n, n)
    dt = 1e-2/mu
    times = sorted((float(t), t) for t in os.listdir(case)
                   if re.fullmatch(r"[0-9.eE+-]+", t) and os.path.isfile(os.path.join(case, t, "psi")))
    bystep = {int(round(tv/dt)): tn for tv, tn in times}
    last = max(bystep)
    for k in range(1, last + 1):
        psi = step(psi, h, mu, dt, band, qname)
        if k in bystep and (k % 1000 == 0 or k == last or k == 50):
            of = internal(os.path.join(case, bystep[k], "psi")).reshape(n, n)
            print(f"step {k}: max|copy - OpenFOAM| = {np.abs(psi - of).max():.3e}, "
                  f"max|OpenFOAM - d| = {np.abs(of - d).max():.3e}", flush=True)


def cmd_table():
    print("== 1. the band width and the taper, centred q, psi0 = d, N = 100, mu dt = 0.01, mu t = 80")
    for band, tp in ((3, False), (6, False), (12, False), (3, True)):
        h, d, psi, _ = run("centred", 100, band=band, use_taper=tp)
        print(f"band {band:2d} taper {tp!s:5s}: " + fmt(metrics(h, d, psi)), flush=True)
    print("== 2. mu dt, centred and godunov, band 3 and 12, N = 100: sawtooth amplitude at mu t = 40 and 80")
    for qn in ("centred", "godunov"):
        for band in (3, 12):
            for mudt in (0.001, 0.01, 0.1):
                h, d, psi, snaps = run(qn, 100, band=band, mudt=mudt, marks=(40,))
                c40 = metrics(h, d, snaps[40])["cbd"]; c80 = metrics(h, d, psi)["cbd"]
                print(f"{qn:8s} band {band:2d} mu dt {mudt:5.3f}: cb {c40:.2e} -> {c80:.2e}, "
                      f"rate {np.log(c80/c40)/40:.3f} mu", flush=True)
    print("== 3. the estimators, psi0 = d and 1.5 d, band 3, mu dt = 0.01, mu t = 80, N = 50, 100, 200")
    for qn, anchor in (("centred", False), ("upwindFace", False), ("upwindFace", True), ("ownFit", False),
                       ("upwindQuad", False), ("upwindFit", False)):
        for scale in (1.0, 1.5):
            for n in (50, 100, 200):
                t0 = time.time()
                h, d, psi, _ = run(qn, n, scale=scale, anchor=anchor)
                print(f"{qn:10s} anchor {anchor!s:5s} psi0 {scale}d N {n:3d}: " + fmt(metrics(h, d, psi))
                      + f" ({time.time() - t0:.0f} s)", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate"); v.add_argument("case")
    r = sub.add_parser("run"); r.add_argument("q"); r.add_argument("n", type=int)
    r.add_argument("--band", type=float, default=3.0); r.add_argument("--mudt", type=float, default=0.01)
    r.add_argument("--mut", type=float, default=80.0); r.add_argument("--scale", type=float, default=1.0)
    r.add_argument("--anchor", action="store_true"); r.add_argument("--taper", action="store_true")
    sub.add_parser("table")
    a = ap.parse_args()
    if a.cmd == "validate":
        cmd_validate(a.case)
    elif a.cmd == "run":
        h, d, psi, _ = run(a.q, a.n, band=a.band, mudt=a.mudt, mut=a.mut, scale=a.scale, anchor=a.anchor,
                           use_taper=a.taper)
        print(f"{a.q} N {a.n}: " + fmt(metrics(h, d, psi)))
    else:
        cmd_table()


if __name__ == "__main__":
    main()
