#!/usr/bin/env python3
"""Geometry of the swept-face figure (fig_sweptflux.tex, panels a and b).

The face f is the unit square x = 1, 0 <= y, z <= 1, with the area vector S_f = (1, 0, 0)
from the donor (x < 1) to the neighbour (x > 1); a uniform velocity with the displacement
U dt = d per step traces it back to the pre-image f - d. The swept volume is the closed list
of 4F = 16 triangles (face fan, reversed pre-image fan, two triangles per lateral
quadrilateral), outward for phi_f > 0. Panel b cuts it by the plane (x - x0).n = 0, liquid
where (x - x0).n <= 0: the wetted part of every triangle, the cap, and the liquid volume by
V^- = (1/3) sum_T (p_T - x0).a_T^- (checked against the cap's polygon area for the closure and
against a midpoint quadrature of the liquid indicator for the volume). Writes fig_sweptflux_a.tex and fig_sweptflux_b.tex (TikZ drawing
commands in the oblique projection (u, v) = (y + 0.7 x, z + 0.35 x)).
Run on Lichtenberg, not on the laptop:  python3 make_fig_sweptflux.py
"""
import math

D = (0.55, 0.30, 0.22)                      # U dt
N_RAW = (0.30, 0.20, 1.00)                  # the donor's plane
X0 = (0.75, 0.35, 0.50)


def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def mul(s, a): return tuple(s * x for x in a)
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b): return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])
def norm(a): return math.sqrt(dot(a, a))


NRM = mul(1.0 / norm(N_RAW), N_RAW)
P = [(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)]
C = (1, 0.5, 0.5)
Pp = [sub(p, D) for p in P]
Cp = sub(C, D)

tris = []                                     # (kind, a, b, c)
for i in range(4):
    j = (i + 1) % 4
    tris.append(("face", C, P[i], P[j]))
    tris.append(("pre", Cp, Pp[j], Pp[i]))
    tris.append(("lat", P[i], Pp[i], Pp[j]))
    tris.append(("lat", P[i], Pp[j], P[j]))


def area_vec(poly):
    s = (0.0, 0.0, 0.0)
    for k in range(1, len(poly) - 1):
        s = add(s, mul(0.5, cross(sub(poly[k], poly[0]), sub(poly[k + 1], poly[0]))))
    return s


def clip(poly):
    out = []
    for k in range(len(poly)):
        a, b = poly[k], poly[(k + 1) % len(poly)]
        fa, fb = dot(sub(a, X0), NRM), dot(sub(b, X0), NRM)
        if fa <= 0:
            out.append(a)
        if fa * fb < 0:
            out.append(add(a, mul(fa / (fa - fb), sub(b, a))))
    return out


vol = sum(dot(t[1], cross(t[2], t[3])) for t in tris) / 6.0
closure = norm(tuple(sum(area_vec(t[1:])[k] for t in tris) for k in range(3)))
wet = [(t[0], clip(list(t[1:]))) for t in tris]
aw = [area_vec(w) if len(w) >= 3 else (0.0, 0.0, 0.0) for _, w in wet]
vminus = sum(dot(sub(w[0], X0), a) for (_, w), a in zip(wet, aw) if len(w) >= 3) / 3.0
cap_closure = -dot(NRM, tuple(sum(a[k] for a in aw) for k in range(3)))

# the cap: the plane's crossings of the triangle edges, ordered about their mean
pts = []
for t in tris:
    for k in range(3):
        a, b = t[1 + k], t[1 + (k + 1) % 3]
        fa, fb = dot(sub(a, X0), NRM), dot(sub(b, X0), NRM)
        if fa * fb < 0:
            q = add(a, mul(fa / (fa - fb), sub(b, a)))
            if all(norm(sub(q, r)) > 1e-9 for r in pts):
                pts.append(q)
cm = mul(1.0 / len(pts), tuple(sum(p[k] for p in pts) for k in range(3)))
e1 = sub(pts[0], cm); e1 = mul(1.0 / norm(e1), e1)
e2 = cross(NRM, e1)
pts.sort(key=lambda p: math.atan2(dot(sub(p, cm), e2), dot(sub(p, cm), e1)))
cap_area = dot(area_vec(pts), NRM)

# independent reference: the swept volume is {(1, eta, zeta) - s d : s, eta, zeta in [0, 1]}
# (Jacobian d_x); the midpoint rule on the liquid indicator in (s, eta, zeta), n^3 points
NQ = 160
acc = 0
for i in range(NQ):
    s = (i + 0.5) / NQ
    for j in range(NQ):
        eta = (j + 0.5) / NQ
        base = (1 - s * D[0] - X0[0]) * NRM[0] + (eta - s * D[1] - X0[1]) * NRM[1]
        for k in range(NQ):
            if base + ((k + 0.5) / NQ - s * D[2] - X0[2]) * NRM[2] <= 0:
                acc += 1
vref = D[0] * acc / NQ**3


def uv(p): return (p[1] + 0.7 * p[0], p[2] + 0.35 * p[0])
def path(poly): return " -- ".join("(%.4f,%.4f)" % uv(p) for p in poly) + " -- cycle"
def seg(a, b): return "(%.4f,%.4f) -- (%.4f,%.4f)" % (uv(a) + uv(b))


hdr = ("%% written by make_fig_sweptflux.py -- do not edit. d = %s, n = (%.5f, %.5f, %.5f), x0 = %s\n"
       "%% signed volume %.12f (d.S_f = %.12f), closure |sum a_T| = %.1e, V^- = %.12f, midpoint quadrature (160^3) %.6f,\n"
       "%% cap area by closure %.12f, by its polygon %.12f\n"
       % (D, NRM[0], NRM[1], NRM[2], X0, vol, D[0], closure, vminus, vref, cap_closure, cap_area))

with open("fig_sweptflux_a.tex", "w") as fh:
    fh.write(hdr)
    for kind, a, b, c in tris:
        col = {"face": "facecol!60", "pre": "precol!70", "lat": "black!35"}[kind]
        fh.write("\\draw[%s, line width=0.3pt] %s;\n" % (col, path([a, b, c])))
    fh.write("\\draw[facecol, line width=1.4pt] %s;\n" % path(P))
    fh.write("\\draw[precol, line width=1.2pt, dashed] %s;\n" % path(Pp))
    for p, q in zip(P + [C], Pp + [Cp]):
        fh.write("\\draw[->, black!65, line width=0.6pt] %s;\n" % seg(p, add(q, mul(0.06, sub(p, q)))))
    fh.write("\\fill[facecol] (%.4f,%.4f) circle[radius=0.018];\n" % uv(C))
    fh.write("\\fill[precol] (%.4f,%.4f) circle[radius=0.018];\n" % uv(Cp))
    fh.write("\\coordinate (fC) at (%.4f,%.4f);\n\\coordinate (pC) at (%.4f,%.4f);\n" % (uv(C) + uv(Cp)))
    fh.write("\\coordinate (fS) at (%.4f,%.4f);\n" % uv(add(C, (0.45, 0, 0))))
    fh.write("\\coordinate (fP) at (%.4f,%.4f);\n\\coordinate (pP) at (%.4f,%.4f);\n" % (uv(P[2]) + uv(Pp[0])))

with open("fig_sweptflux_b.tex", "w") as fh:
    fh.write(hdr)
    for kind, a, b, c in tris:
        fh.write("\\draw[black!30, line width=0.3pt] %s;\n" % path([a, b, c]))
    for (kind, w), a in zip(wet, aw):
        if len(w) >= 3 and norm(a) > 1e-12:
            col = "wetcol" if kind != "pre" else "wetcol!70!orange"
            fh.write("\\fill[%s, opacity=0.45] %s;\n" % (col, path(w)))
    fh.write("\\fill[capcol, opacity=0.55] %s;\n" % path(pts))
    fh.write("\\draw[ifacecol, line width=1.0pt] %s;\n" % path(pts))
    fh.write("\\draw[facecol, line width=1.4pt] %s;\n" % path(P))
    fh.write("\\draw[precol, line width=1.2pt, dashed] %s;\n" % path(Pp))
    # the cone of one wetted pre-image triangle from x0 on the plane
    k = next(i for i, (kind, w) in enumerate(wet) if kind == "pre" and len(w) == 3)
    w = wet[k][1]
    fh.write("\\fill[black, opacity=0.10] %s;\n" % path(w))
    for p in w:
        fh.write("\\draw[black!70, dashed, line width=0.5pt] %s;\n" % seg(X0, p))
    fh.write("\\fill[black] (%.4f,%.4f) circle[radius=0.02];\n" % uv(X0))
    fh.write("\\coordinate (xz) at (%.4f,%.4f);\n" % uv(X0))
    fh.write("\\coordinate (capC) at (%.4f,%.4f);\n" % uv(cm))
    fh.write("\\coordinate (capN) at (%.4f,%.4f);\n" % uv(add(cm, mul(0.35, NRM))))
    tc = mul(1.0 / 3, add(add(w[0], w[1]), w[2]))
    fh.write("\\coordinate (tC) at (%.4f,%.4f);\n" % uv(tc))

print(hdr)
