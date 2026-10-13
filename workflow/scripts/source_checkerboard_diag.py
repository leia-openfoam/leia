#!/usr/bin/env python3
"""The zero-flow source gate (config/gates/sourceZeroFlow2D.yaml): the checkerboard amplitude of psi.

For every written time of a case of cases/2Dtranslation (single block, N x N cells, x fastest):
  cb    = D2(D2 psi)/64 in the band |d| <= 3h, D2 the 5-point second difference. A mode
          eps (-1)^(i+j) g gives D2 = -8 eps (-1)^(i+j) g, so D2(D2) = 64 times the mode; the smooth
          SDF of the circle leaves about h^4/R^3 (3e-6 at N = 100).
  q     = the centred |grad psi| (the q of the source step on a uniform mesh), band mean
  dev   = RMS of psi - d in the band (d the exact signed distance of the circle)
The growth rate of cb is a least-squares fit of log(cb_rms) over the samples between 30 times the
first sample and 1e-3. The zero-set displacement is E_GEOM_ALPHA / (perimeter * depth) at T.

The line "final:" (gate E2.1, config/gates/sourceZeroFlowUpwind2D.yaml) is read at the last time,
on psi - d, with the definitions of the numpy copy of the source step:
  shift  the zero-crossing shift along x and y between neighbours with a sign change in psi and
         in d: |psi_a/(psi_a - psi_b) - d_a/(d_a - d_b)|, max and RMS, in units of h
  dev    rms(psi - d)/h in |d| < 3h
  cbd    the RMS of D2(D2 (psi - d))/64 in |d| < 3h, at the last time and at half of it
  ratio  the median of psi/d in |d| < h
All of them exclude a border of 3 cells.

Usage: source_checkerboard_diag.py <case dir> [<case dir> ...]
"""
import csv, json, math, os, re, sys
import numpy as np

XC, YC, R, DEPTH = 0.25, 0.5, 0.15, 0.1

def internal(path):
    t = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", t)
    a = m.end(); b = t.index(")", a)
    return np.array([float(x) for x in t[a:b].split()])

def analyse(case):
    tok = json.load(open(os.path.join(case, "case_params.json")))["tokens"]
    n = int(tok["N_CELLS"]); h = 1.0/n; mu = float(tok["GC_MU"]); band = tok["SL_SOURCE_BAND_CELLS"]
    i = np.arange(n); x = (i + 0.5)*h
    X, Y = np.meshgrid(x, x)                       # [j][i]
    d = np.sqrt((X - XC)**2 + (Y - YC)**2) - R
    inner = np.zeros((n, n), bool); inner[2:-2, 2:-2] = True
    bandmask = (np.abs(d) <= 3*h) & inner
    times = sorted((float(t), t) for t in os.listdir(case)
                   if re.fullmatch(r"[0-9.eE+-]+", t) and os.path.isfile(os.path.join(case, t, "psi")))
    rows = []
    inner3 = np.zeros((n, n), bool); inner3[3:-3, 3:-3] = True
    b3 = inner3 & (np.abs(d) < 3*h)
    def d2z(f):
        g = np.zeros_like(f)
        g[1:-1, 1:-1] = f[2:, 1:-1] + f[:-2, 1:-1] + f[1:-1, 2:] + f[1:-1, :-2] - 4*f[1:-1, 1:-1]
        return g
    def cbd_of(psi):
        return float(np.sqrt(np.mean((d2z(d2z(psi - d))/64.0)[b3]**2)))
    half = times[len(times)//2] if times else None
    for tv, tn in times:
        psi = internal(os.path.join(case, tn, "psi")).reshape(n, n)
        def d2(f):
            g = np.zeros_like(f)
            g[1:-1, 1:-1] = f[2:, 1:-1] + f[:-2, 1:-1] + f[1:-1, 2:] + f[1:-1, :-2] - 4*f[1:-1, 1:-1]
            return g
        cb = d2(d2(psi))/64.0
        gx = np.zeros_like(psi); gy = np.zeros_like(psi)
        gx[:, 1:-1] = (psi[:, 2:] - psi[:, :-2])/(2*h); gy[1:-1, :] = (psi[2:, :] - psi[:-2, :])/(2*h)
        q = np.sqrt(gx**2 + gy**2)
        rows.append((tv, float(np.sqrt(np.mean(cb[bandmask]**2))), float(np.max(np.abs(cb[bandmask]))),
                     float(np.mean(q[bandmask])), float(np.sqrt(np.mean((psi - d)[bandmask]**2)))))
    rate = float("nan")
    if len(rows) > 3:
        c0 = rows[0][1] if rows[0][1] > 0 else rows[1][1]
        sel = [(t, c) for t, c, *_ in rows if 30*c0 < c < 1e-3]
        if len(sel) >= 3:
            ts = np.array([s[0] for s in sel]); ls = np.log([s[1] for s in sel])
            rate = float(np.polyfit(ts, ls, 1)[0])
    csvp = os.path.join(case, "leiaSemiLagrangeLevelSetFoam.csv")
    disp = float("nan")
    if os.path.isfile(csvp):
        r = list(csv.DictReader(open(csvp)))
        disp = float(r[-1]["E_GEOM_ALPHA"]) / (2*math.pi*R*DEPTH) / h
    final = None
    if times:
        psi = internal(os.path.join(case, times[-1][1], "psi")).reshape(n, n)
        psih = internal(os.path.join(case, half[1], "psi")).reshape(n, n)
        sh = []
        for a, b in (((slice(3, n-3), slice(3, n-4)), (slice(3, n-3), slice(4, n-3))),
                     ((slice(3, n-4), slice(3, n-3)), (slice(4, n-3), slice(3, n-3)))):
            pa, pb, da, db = psi[a], psi[b], d[a], d[b]
            m = (pa*pb < 0) & (da*db < 0)
            sh.append(np.abs(pa[m]/(pa[m] - pb[m]) - da[m]/(da[m] - db[m])))
        sh = np.concatenate(sh)
        fi = inner3 & (np.abs(d) < h)
        final = dict(t=times[-1][0], mut=times[-1][0]*mu, shift_max=float(sh.max()),
                     shift_rms=float(np.sqrt(np.mean(sh**2))),
                     dev=float(np.sqrt(np.mean(((psi - d)/h)[b3]**2))),
                     cbd=cbd_of(psi), cbd_half=cbd_of(psih), mut_half=half[0]*mu,
                     ratio=float(np.median((psi/d)[fi])))
    return dict(case=os.path.basename(case), mu=mu, band=band, n=n, qs=tok.get("GC_Q_SCHEME", "centred"),
                rows=rows, rate=rate, disp=disp, final=final)

for case in sys.argv[1:]:
    r = analyse(case)
    first, last = r["rows"][0], r["rows"][-1]
    print(f"{r['case']}: mu={r['mu']:g} bandCells={r['band']} | cb_rms {first[1]:.2e} -> {last[1]:.2e} "
          f"(max {last[2]:.2e}), fitted rate {r['rate']:.1f} 1/s = {r['rate']/r['mu']:.3f} mu | "
          f"q band mean {first[3]:.4f} -> {last[3]:.4f} | rms(psi-d) {first[4]:.2e} -> {last[4]:.2e} | "
          f"zero-set displacement at T {r['disp']:.3f} h")
    f = r["final"]
    if f:
        print(f"  final: N={r['n']} q={r['qs']} mu t={f['mut']:.0f} | zero shift/h max {f['shift_max']:.4f} "
              f"rms {f['shift_rms']:.4f} | rms(psi-d)/h {f['dev']:.3e} | cb(psi-d) {f['cbd']:.2e} "
              f"(mu t={f['mut_half']:.0f}: {f['cbd_half']:.2e}) | median psi/d {f['ratio']:.4f}")
