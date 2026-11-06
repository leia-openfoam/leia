#!/usr/bin/env python3
"""The independent truth for the dumped swept pieces (2026-11-06): the pieces on which
the fork's closure volume (davof::sweptPieceVolume, the cap system from the closure of
the clipped wet boundary) and the exact volume (davof::sweptPieceExactVolume, signed
tetrahedra clipped as convex polyhedra) disagree, dumped by isoAdvection's key
sweptPieceDumpTol as the raw closed triangle list of the swept volume, the cell's
half-spaces and the PLIC half-space (postProcessing/sweptPieceDump/*.txt).

Two independent volumes of the piece P = {swept volume} ∩ {cell} ∩ {PLIC half-space}:
  tets  the signed tetrahedra of the closed list about a point, each clipped by the
        half-spaces as a convex polyhedron (Sutherland-Hodgman face by face, the cap
        of each clip from the crossing points; I. E. Sutherland, G. W. Hodgman,
        Reentrant polygon clipping, Commun. ACM 17 (1974) 32-42,
        doi:10.1145/360767.360802), the volume by the divergence theorem; written
        here independently of the fork's C++ (numpy, a different apex and a
        different cap construction);
  MC    Monte Carlo over the bounding box of the list: a sample counts with its
        generalized winding number w(p) with respect to the closed swept list
        (the sum of the signed solid angles of the triangles, A. Van Oosterom,
        J. Strackee, The solid angle of a plane triangle, IEEE Trans. Biomed. Eng.
        BME-30 (1983) 125-126, doi:10.1109/TBME.1983.325207; A. Jacobson, L. Kavan,
        O. Sorkine-Hornung, Robust inside-outside segmentation using generalized
        winding numbers, ACM Trans. Graph. 32 (2013) 33, doi:10.1145/2461912.2461916)
        where it lies inside every half-space; V = box volume x mean(w), the error
        box volume x std(w)/sqrt(N). A self-overlapping list shows |w| = 2 samples.
Also the topology of the PLIC section: the segments of the PLIC plane with the
swept list clipped to the cell, chained into loops (a cap polygon that is not simply
connected shows two or more loops).
    python3 piece_truth.py <case dir> [N samples, default 1000000] [n worst, default 10]
"""
import glob, os, re, sys, math
import numpy as np

VEC = re.compile(r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)")


def parse(fn):
    L = open(fn).read().split("\n")
    kv = L[0].split()
    head = {kv[i]: kv[i + 1] for i in range(0, len(kv) - 1, 2)}
    v = lambda s: np.array([float(x) for x in VEC.search(s).groups()])
    i = 1
    xref = v(L[i]); i += 1
    xdon = v(L[i]); i += 1
    m = VEC.findall(L[i]); tail = L[i].split(")")[-1].split()
    plic = dict(n=np.array([float(x) for x in m[0]]), x0=np.array([float(x) for x in m[1]]), axis=int(tail[0]), sgn=float(tail[1])); i += 1
    np_ = int(L[i].split()[1]); i += 1
    planes = []
    for _ in range(np_):
        m = VEC.findall(L[i]); tail = L[i].split(")")[-1].split()
        planes.append(dict(n=np.array([float(x) for x in m[0]]), x0=np.array([float(x) for x in m[1]]), axis=int(tail[0]), sgn=float(tail[1]))); i += 1
    nt = int(L[i].split()[1]); i += 1
    T = np.array([[[float(x) for x in t] for t in VEC.findall(L[i + k])] for k in range(nt)])
    return head, xref, xdon, plic, planes, T


def dist(hs, p):
    """signed distance(s); inside where <= 0 (the fork's convention)"""
    p = np.asarray(p)
    if hs["axis"] < 0:
        return (p - hs["x0"]) @ hs["n"]
    return hs["sgn"]*(p[..., hs["axis"]] - hs["x0"][hs["axis"]])


def normal_of(hs):
    if hs["axis"] < 0:
        n = hs["n"]
    else:
        n = np.zeros(3); n[hs["axis"]] = hs["sgn"]
    return n/np.linalg.norm(n)


def clip_poly(poly, hs):
    """Sutherland-Hodgman of one polygon (list of points) by the half-space; returns the
    kept polygon and the crossing points"""
    out, cross = [], []
    n = len(poly)
    for j in range(n):
        p0, p1 = poly[j], poly[(j + 1) % n]
        s0, s1 = dist(hs, p0), dist(hs, p1)
        if s0 <= 0: out.append(p0)
        if (s0 <= 0) != (s1 <= 0):
            x = p0 + (p1 - p0)*(s0/(s0 - s1))
            out.append(x); cross.append(x)
    return (out if len(out) >= 3 else []), cross


def clip_polyhedron(polys, hs):
    """the convex polyhedron (outward polygons) clipped; the cap closed from the crossings"""
    new, cross = [], []
    for q in polys:
        k, c = clip_poly(q, hs)
        if k: new.append(k)
        cross += c
    if len(cross) >= 3:
        nrm = normal_of(hs)
        cen = np.mean(cross, axis=0)
        u = np.array([1.0, 0, 0]) if abs(nrm[0]) < 0.9 else np.array([0, 1.0, 0])
        u -= (u @ nrm)*nrm; u /= np.linalg.norm(u); vv = np.cross(nrm, u)
        ang = [math.atan2((x - cen) @ vv, (x - cen) @ u) for x in cross]
        cap = [cross[i] for i in np.argsort(ang)]
        a = np.zeros(3)
        for j in range(len(cap)):
            a += np.cross(cap[j], cap[(j + 1) % len(cap)])
        if a @ nrm < 0: cap = cap[::-1]
        new.append(cap)
    return new


def poly_volume_term(polys, o):
    v3 = 0.0
    for q in polys:
        a = np.zeros(3)
        for j in range(len(q)):
            a += np.cross(q[j] - o, q[(j + 1) % len(q)] - o)
        v3 += (q[0] - o) @ (0.5*a)
    return v3/3


def tets_volume(T, planes, plic, wet, apex):
    """the signed tetrahedra about apex, each clipped as a convex polyhedron"""
    vol = 0.0
    for t in T:
        a = 0.5*np.cross(t[1] - t[0], t[2] - t[0])
        vs = ((t[0] - apex) @ a)/3
        if abs(vs) < 1e-30: continue
        sgn = 1.0 if vs > 0 else -1.0
        A, B, C = t[0], t[1], t[2]
        if sgn < 0: B, C = C, B
        polys = [[A, B, C], [apex, B, A], [apex, C, B], [apex, A, C]]
        for hs in planes:
            polys = clip_polyhedron(polys, hs)
            if not polys: break
        if polys and wet:
            polys = clip_polyhedron(polys, plic)
        if polys:
            vol += sgn*poly_volume_term(polys, apex)
    return vol


def winding(P, T):
    """generalized winding numbers of the points P (n x 3) w.r.t. the closed list T"""
    w = np.zeros(len(P))
    for t in T:
        a = t[0] - P; b = t[1] - P; c = t[2] - P
        la = np.linalg.norm(a, axis=1); lb = np.linalg.norm(b, axis=1); lc = np.linalg.norm(c, axis=1)
        num = np.einsum("ij,ij->i", a, np.cross(b, c))
        den = la*lb*lc + np.einsum("ij,ij->i", a, b)*lc + np.einsum("ij,ij->i", a, c)*lb + np.einsum("ij,ij->i", b, c)*la
        w += 2*np.arctan2(num, den)
    return w/(4*math.pi)


def mc_volume(T, planes, plic, wet, N, seed=1):
    lo = T.reshape(-1, 3).min(axis=0); hi = T.reshape(-1, 3).max(axis=0)
    # the box can be cut by the cell's axis planes
    for hs in planes:
        if hs["axis"] >= 0:
            x = hs["x0"][hs["axis"]]
            if hs["sgn"] > 0: hi[hs["axis"]] = min(hi[hs["axis"]], x)
            else: lo[hs["axis"]] = max(lo[hs["axis"]], x)
    if np.any(hi <= lo): return 0.0, 0.0, 0.0, 0
    box = float(np.prod(hi - lo))
    rng = np.random.default_rng(seed)
    vals = np.zeros(N); done = 0; nTwo = 0
    while done < N:
        m = min(100000, N - done)
        P = lo + rng.random((m, 3))*(hi - lo)
        inside = np.ones(m, dtype=bool)
        for hs in planes:
            inside &= dist(hs, P) <= 0
        if wet:
            inside &= dist(plic, P) <= 0
        w = np.zeros(m)
        if inside.any():
            w[inside] = winding(P[inside], T)
        vals[done:done + m] = w
        nTwo += int(np.sum(np.abs(w) > 1.5))
        done += m
    return box*vals.mean(), box*vals.std()/math.sqrt(N), box, nTwo


def plic_section_loops(T, planes, plic):
    """the segments of the PLIC plane with the swept list clipped to the cell, chained"""
    segs = []
    for t in T:
        poly = [t[0], t[1], t[2]]
        for hs in planes:
            poly, _ = clip_poly(poly, hs)
            if not poly: break
        if not poly: continue
        _, cross = clip_poly(poly, plic)
        if len(cross) == 2:
            segs.append((cross[0], cross[1]))
    if not segs: return 0, 0
    L = max(np.linalg.norm(s[1] - s[0]) for s in segs)
    tol = 1e-7*max(L, 1e-30)
    used = [False]*len(segs); loops = 0; open_chains = 0
    for i in range(len(segs)):
        if used[i]: continue
        used[i] = True; start = segs[i][0]; cur = segs[i][1]; closed = False
        while True:
            if np.linalg.norm(cur - start) <= tol: closed = True; break
            found = False
            for j in range(len(segs)):
                if used[j]: continue
                for a, b in ((segs[j][0], segs[j][1]), (segs[j][1], segs[j][0])):
                    if np.linalg.norm(a - cur) <= tol:
                        used[j] = True; cur = b; found = True; break
                if found: break
            if not found: break
        if closed: loops += 1
        else: open_chains += 1
    return loops, open_chains


def main():
    case = sys.argv[1]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 1000000
    nworst = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    files = glob.glob(os.path.join(case, "postProcessing", "sweptPieceDump", "*.txt"))
    rows = []
    for fn in files:
        head, xref, xdon, plic, planes, T = parse(fn)
        Vc = float(head["Vc"]); cl = float(head["closureVolume"]); ex = float(head["exactVolume"])
        rows.append((abs(cl - ex)/Vc, fn, head, xref, xdon, plic, planes, T))
    rows.sort(key=lambda r: -r[0])
    print("dumped pieces: %d; the %d worst by |closure - exact|/V_c" % (len(rows), nworst))
    print("file | step cell nCaps resid | closure/Vc exact(C++)/Vc | tets(py)/Vc | MC/Vc +- sigma (|w|=2 samples) | PLIC-section loops/open | verdict")
    for d, fn, head, xref, xdon, plic, planes, T in rows[:nworst]:
        Vc = float(head["Vc"]); cl = float(head["closureVolume"]); ex = float(head["exactVolume"])
        wet = head["wetPlane"] in ("1", "true")
        apex = T.reshape(-1, 3).mean(axis=0) + np.array([0.37, 0.11, 0.23])*1e-3*(T.reshape(-1, 3).max() - T.reshape(-1, 3).min())
        vt = tets_volume(T, planes, plic, wet, apex)
        vm, sm, box, nTwo = mc_volume(T, planes, plic, wet, N)
        loops, opens = plic_section_loops(T, planes, plic) if wet else (0, 0)
        okE = abs(ex - vm) <= 3*sm + 1e-12*Vc; okC = abs(cl - vm) <= 3*sm + 1e-12*Vc
        verdict = ("exact agrees with MC" if okE else "exact DISAGREES with MC") + ", " + ("closure agrees with MC" if okC else "closure DISAGREES with MC")
        print("%s | %s %s %s %s | %.6e %.6e | %.6e | %.6e +- %.1e (%d) | %d/%d | %s" % (
            os.path.basename(fn), head["step"], head["cell"], head["nCaps"], head["closureResidual"][:9], cl/Vc, ex/Vc, vt/Vc, vm/Vc, sm/Vc, nTwo, loops, opens, verdict))
    # the whole set: how many agree
    nE = nC = 0; nTot = 0
    for d, fn, head, xref, xdon, plic, planes, T in rows:
        Vc = float(head["Vc"]); cl = float(head["closureVolume"]); ex = float(head["exactVolume"]); wet = head["wetPlane"] in ("1", "true")
        vm, sm, box, nTwo = mc_volume(T, planes, plic, wet, min(N, 200000), seed=2)
        nTot += 1
        if abs(ex - vm) <= 3*sm + 1e-12*Vc: nE += 1
        if abs(cl - vm) <= 3*sm + 1e-12*Vc: nC += 1
    print("all %d dumped pieces (2e5 samples each): exact agrees with MC on %d, closure on %d" % (nTot, nE, nC))


if __name__ == "__main__":
    main()
