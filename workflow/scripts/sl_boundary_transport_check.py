#!/usr/bin/env python3
"""Far-field level-set transport check for a translating droplet (Popinet cases).

Compares the written psi against the exact signed distance of the translated sphere
(centre c0 + U t e_x, radius R; psi positive outside) and reports, by region, the error
and the TRANSPORT FRACTION of the first cell layer at each boundary:

    (psi(T) - psi(0)) / (psi_exact(T) - psi_exact(0))      median over the layer

1.00 = the level set moves with the flow; ~0 = the layer is frozen. MEASURED 2026-09-05
with the boundary faces in the stencil (`stencilBoundaryFaces include`): side walls
1.002 on every mesh, inlet 0.009 (2D) / 0.20 (3D hex) / 0.08 (3D poly), outlet 0.96 /
0.49 / 0.24. With `inflowOnly` the outlet reads 1.000; the inlet stays pinned by its
zeroGradient condition (stable, harmless while the incoming level set stays away from 0).

    python3 workflow/scripts/sl_boundary_transport_check.py <case> [--time T] [--centre X Y Z]
                                                            [--radius R] [--U 1] [--xlen 2]
"""
import argparse, math, os, re, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import leia_refine as lr


def read_vectors(path, n):
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*(\d+)\s*\(", txt)
    body = txt[m.end():]; body = body[:body.index("\n)\n")]
    return [[float(v) for v in s.split()] for s in re.findall(r"\(([^)]*)\)", body)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case"); ap.add_argument("--time", default=None)
    ap.add_argument("--centre", nargs=3, type=float, default=[0.5, 0.5, 0.5])
    ap.add_argument("--radius", type=float, default=0.2); ap.add_argument("--U", type=float, default=1.0)
    ap.add_argument("--xlen", type=float, default=2.0); ap.add_argument("--h", type=float, default=None)
    a = ap.parse_args(); os.chdir(a.case); n = lr.n_cells(".")
    times = sorted([d for d in os.listdir(".") if re.match(r"^[0-9]*\.?[0-9]+(e-?[0-9]+)?$", d) and float(d) > 0], key=float)
    t = a.time or times[-1]; T = float(t)
    if not os.path.exists(f"{t}/C"):
        subprocess.run(f"postProcess -func writeCellCentres -time {t} > log.writeCellCentres 2>&1", shell=True)
    C = read_vectors(f"{t}/C", n); psi = lr.read_scalar_field(f"{t}/psi", n); psi0 = lr.read_scalar_field("0/psi", n)
    twoD = max(c[2] for c in C) - min(c[2] for c in C) < 1e-9
    cx, cy, cz = a.centre; R = a.radius
    def ex(p, tt):
        dz = 0.0 if twoD else (p[2] - cz) ** 2
        return math.sqrt((p[0] - cx - a.U * tt) ** 2 + (p[1] - cy) ** 2 + dz) - R
    if a.h:
        h = a.h
    else:
        # median cube-root cell volume (2D: sqrt of volume / thickness) -- mesh-agnostic
        if not os.path.exists("0/V"):
            subprocess.run("postProcess -func writeCellVolumes -time 0 > log.writeCellVolumes 2>&1", shell=True)
        V = sorted(lr.read_scalar_field("0/V", n)); vm = V[len(V) // 2]
        if twoD:
            zs = sorted(set(round(c[2], 9) for c in C[:2000])); thick = 2 * abs(zs[0]) if len(zs) == 1 else 1.0
            h = math.sqrt(vm / thick)
        else:
            h = vm ** (1 / 3)
    err = [psi[i] - ex(C[i], T) for i in range(n)]
    def region(c):
        bx = min(c[0], a.xlen - c[0]); wall = min(c[1], 1 - c[1]) if twoD else min(c[1], 1 - c[1], c[2], 1 - c[2])
        if bx < 1.5 * h and wall < 1.5 * h: return "edge/corner"
        if c[0] < 1.5 * h: return "inlet layer"
        if a.xlen - c[0] < 1.5 * h: return "outlet layer"
        if wall < 1.5 * h: return "wall layer"
        return "interior"
    groups = {}
    for i in range(n): groups.setdefault(region(C[i]), []).append(i)
    print(f"{a.case}: t = {t}, {n} cells, h = {h:.4e}, exact centre x = {cx + a.U * T:.3f}")
    for k in ("interior", "wall layer", "inlet layer", "outlet layer", "edge/corner"):
        idx = groups.get(k, [])
        if not idx: continue
        e = sorted(abs(err[i]) for i in idx); se = sorted(err[i] for i in idx)
        print(f"  {k:13s} n={len(idx):7d} |err| max={e[-1]:.3e} median={e[len(e)//2]:.3e}  signed min={se[0]:+.3e} max={se[-1]:+.3e}")
    far = [i for i in range(n) if ex(C[i], T) > 5 * h]
    print(f"  far cells with psi < 0 (a fake zero set): {sum(1 for i in far if psi[i] < 0)}")
    mid = (lambda c: 0.2 < c[1] < 0.8) if twoD else (lambda c: 0.2 < c[1] < 0.8 and 0.2 < c[2] < 0.8)
    for name, sel in (("inlet first layer", lambda c: c[0] < 0.8 * h and mid(c)),
                      ("outlet first layer", lambda c: a.xlen - c[0] < 0.8 * h and mid(c)),
                      ("side-wall first layer", lambda c: c[1] < 0.8 * h and 0.3 < c[0] < a.xlen - 0.3)):
        idx = [i for i in range(n) if sel(C[i])]
        fr = sorted((psi[i] - psi0[i]) / (ex(C[i], T) - ex(C[i], 0)) for i in idx if abs(ex(C[i], T) - ex(C[i], 0)) > 1e-6)
        if fr: print(f"  transport fraction {name:22s} n={len(idx):5d} median {fr[len(fr)//2]:.3f} (p10 {fr[len(fr)//10]:.3f}, p90 {fr[9*len(fr)//10]:.3f})")


if __name__ == "__main__":
    main()
