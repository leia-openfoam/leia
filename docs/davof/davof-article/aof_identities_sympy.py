#!/usr/bin/env python3
"""Exact verification of the two Area-of-Fluid identities of the article
(sections/area-of-fluid-method-interface-recon.tex, corrected 2026-09-29).

For the liquid sub-volume Omega_k^- of a cell Omega_k cut by ONE plane with
unit normal n (pointing out of the liquid), the Gauss theorem gives

    normal:    n |A|  = - sum_f n_f |S_f^-|   = - sum_f alpha_f S_f        (eq:nsigmasurfsum)
    position:  p_k    = ( d |Omega_k| alpha_k - sum_f alpha_f (x_f - x_k).S_f ) / |A|
                                                                          (eq:positiondiscrete)

with |S_f^-| the wetted area of face f, alpha_f = |S_f^-|/|S_f|, S_f = |S_f| n_f
the outward face area vector, x_f the face centroid, x_k the cell centroid,
|A| the area of the plane polygon inside the cell, d the space dimension, and
p_k = (x - x_k).n the signed distance of the plane from x_k.

Part A  symbolic, 3D: a general plane a x' + b y' + c z' = t cutting the corner
        of the unit cube (x' = x + 1/2, ...): both identities simplify to zero;
        the 2020 sign convention returns -p_k; a volume defect dV moves the
        plane by 3 dV/|A|.
Part B  symbolic, 2D: the same for the unit square (d = 2).
Part C  exact rationals, 3D: random rational planes cutting the cube in every
        configuration (corner, edge, slab, through vertices). The wetted face
        polygons come from exact half-space clipping, the plane polygon from
        chaining the clip edges, the liquid volume from an INDEPENDENT
        tetrahedral decomposition about an interior point. Residuals of both
        identities are exactly 0 (rational arithmetic, no tolerance).

    python3 verification/aof_identities_sympy.py        # prints PASS/FAIL
"""
import random
import sys

import sympy as sp
from sympy import Rational as R

# ----------------------------------------------------------------------------
# Part A: symbolic corner cut of the unit cube centred at the origin
# ----------------------------------------------------------------------------
def part_a():
    a, b, c, t, dV = sp.symbols("a b c t deltaV", positive=True)
    half = R(1, 2)
    # liquid: a(x+1/2) + b(y+1/2) + c(z+1/2) <= t with t < min(a,b,c): the
    # tetrahedron at the corner (-1/2,-1/2,-1/2) with legs t/a, t/b, t/c.
    V = t**3/(6*a*b*c)
    wetted = {  # face -> (outward unit normal, wetted area, face centroid)
        "x-": (sp.Matrix([-1, 0, 0]), t**2/(2*b*c), sp.Matrix([-half, 0, 0])),
        "y-": (sp.Matrix([0, -1, 0]), t**2/(2*a*c), sp.Matrix([0, -half, 0])),
        "z-": (sp.Matrix([0, 0, -1]), t**2/(2*a*b), sp.Matrix([0, 0, -half])),
    }
    xk = sp.zeros(3, 1)
    m = sp.zeros(3, 1)
    for n, A, _ in wetted.values():
        m -= n*A
    n_exact = sp.Matrix([a, b, c])/sp.sqrt(a**2 + b**2 + c**2)
    A_exact = t**2*sp.sqrt(a**2 + b**2 + c**2)/(2*a*b*c)   # the cap triangle
    ok = sp.simplify(m - n_exact*A_exact) == sp.zeros(3, 1)
    print(f"[A] 3D corner cut: normal identity n|A| = -sum n_f|S_f^-|  ... {'ok' if ok else 'FAIL'}")

    num = 3*V - sum(((xf - xk).dot(n))*A for n, A, xf in wetted.values())
    p = sp.simplify(num/A_exact)
    # the plane a x + b y + c z = t - (a+b+c)/2 has signed distance from the origin
    p_exact = (t - (a + b + c)/2)/sp.sqrt(a**2 + b**2 + c**2)
    ok2 = sp.simplify(p - p_exact) == 0
    print(f"[A] 3D corner cut: position p_k = (3V - sum alpha_f (x_f-x_k).S_f)/|A| ... {'ok' if ok2 else 'FAIL'}")
    p2020 = sp.simplify(-num/A_exact)
    ok3 = sp.simplify(p2020 + p_exact) == 0
    print(f"[A] the 2020 numerator (sum - 3V) returns -p_k ... {'ok' if ok3 else 'FAIL'}")
    num_def = 3*(V + dV) - sum(((xf - xk).dot(n))*A for n, A, xf in wetted.values())
    ok4 = sp.simplify(num_def/A_exact - p - 3*dV/A_exact) == 0
    print(f"[A] a volume defect dV moves the plane by 3 dV/|A| ... {'ok' if ok4 else 'FAIL'}")
    return ok and ok2 and ok3 and ok4


# ----------------------------------------------------------------------------
# Part B: symbolic corner cut of the unit square (d = 2)
# ----------------------------------------------------------------------------
def part_b():
    a, b, t = sp.symbols("a b t", positive=True)
    half = R(1, 2)
    # liquid: a(x+1/2) + b(y+1/2) <= t, t < min(a,b): the triangle at (-1/2,-1/2)
    V = t**2/(2*a*b)
    wetted = {
        "x-": (sp.Matrix([-1, 0]), t/b, sp.Matrix([-half, 0])),
        "y-": (sp.Matrix([0, -1]), t/a, sp.Matrix([0, -half])),
    }
    xk = sp.zeros(2, 1)
    m = sp.zeros(2, 1)
    for n, A, _ in wetted.values():
        m -= n*A
    n_exact = sp.Matrix([a, b])/sp.sqrt(a**2 + b**2)
    A_exact = t*sp.sqrt(a**2 + b**2)/(a*b)       # the hypotenuse
    ok = sp.simplify(m - n_exact*A_exact) == sp.zeros(2, 1)
    num = 2*V - sum(((xf - xk).dot(n))*A for n, A, xf in wetted.values())
    p_exact = (t - (a + b)/2)/sp.sqrt(a**2 + b**2)
    ok2 = sp.simplify(num/A_exact - p_exact) == 0
    print(f"[B] 2D corner cut (d = 2): normal {'ok' if ok else 'FAIL'}, position {'ok' if ok2 else 'FAIL'}")
    return ok and ok2


# ----------------------------------------------------------------------------
# Part C: exact rational arithmetic, random planes through the cube
# ----------------------------------------------------------------------------
def cross(u, v):
    return sp.Matrix([u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]])


def clip_polygon(poly, n, s):
    """Sutherland-Hodgman: the part of the polygon with n.x <= s (exact)."""
    out = []
    new_edge = []
    npts = len(poly)
    for i in range(npts):
        P, Q = poly[i], poly[(i + 1) % npts]
        fP, fQ = n.dot(P) - s, n.dot(Q) - s
        if fP <= 0:
            out.append(P)
        if (fP < 0 and fQ > 0) or (fP > 0 and fQ < 0):
            X = P + (fP/(fP - fQ))*(Q - P)
            out.append(X)
            new_edge.append(X)
        elif fP == 0 and i > 0:
            pass
    # vertices exactly on the plane also bound the cap
    on_plane = [P for P in out if n.dot(P) - s == 0]
    return out, on_plane


def polygon_vector_area(poly):
    va = sp.zeros(3, 1)
    for i in range(len(poly)):
        va += cross(poly[i], poly[(i + 1) % len(poly)])
    return va/2


def cube_faces():
    h = R(1, 2)
    V = {}
    for i, (x, y, z) in enumerate([(sx*h, sy*h, sz*h) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]):
        V[(x, y, z)] = sp.Matrix([x, y, z])
    def P(x, y, z):
        return V[(x*h, y*h, z*h)]
    # counter-clockwise seen from outside, with the outward normal
    return [
        ([P(1, -1, -1), P(1, 1, -1), P(1, 1, 1), P(1, -1, 1)], sp.Matrix([1, 0, 0])),
        ([P(-1, -1, -1), P(-1, -1, 1), P(-1, 1, 1), P(-1, 1, -1)], sp.Matrix([-1, 0, 0])),
        ([P(-1, 1, -1), P(-1, 1, 1), P(1, 1, 1), P(1, 1, -1)], sp.Matrix([0, 1, 0])),
        ([P(-1, -1, -1), P(1, -1, -1), P(1, -1, 1), P(-1, -1, 1)], sp.Matrix([0, -1, 0])),
        ([P(-1, -1, 1), P(1, -1, 1), P(1, 1, 1), P(-1, 1, 1)], sp.Matrix([0, 0, 1])),
        ([P(-1, -1, -1), P(-1, 1, -1), P(1, 1, -1), P(1, -1, -1)], sp.Matrix([0, 0, -1])),
    ]


def liquid_polyhedron_volume(face_polys):
    """Volume of the convex polyhedron bounded by the given face polygons, by
    tetrahedral decomposition about the mean of its vertices (an interior
    point of a convex body). Independent of the Gauss identities."""
    verts = []
    for poly in face_polys:
        for v in poly:
            if not any(v == w for w in verts):
                verts.append(v)
    if len(verts) < 4:
        return R(0)
    o = sum(verts, sp.zeros(3, 1))/len(verts)
    vol = R(0)
    for poly in face_polys:
        for i in range(1, len(poly) - 1):
            M = sp.Matrix.hstack(poly[0] - o, poly[i] - o, poly[i + 1] - o)
            vol += abs(M.det())/6
    return vol


def check_plane(n, s):
    """Both identities for the liquid {x : n.x <= s} in the unit cube; exact."""
    faces = cube_faces()
    wet_polys = []          # wetted polygons (for the independent volume)
    m = sp.zeros(3, 1)      # - sum_f alpha_f S_f
    face_moment = R(0)      # sum_f alpha_f (x_f - x_k).S_f, x_k = 0
    cap_edges = []
    for poly, nf in faces:
        wet, on_plane = clip_polygon(poly, n, s)
        if len(wet) >= 3:
            va = polygon_vector_area(wet)
            A_wet = va.dot(nf)                     # exact, >= 0 for ccw polygons
            assert A_wet >= 0
            m -= nf*A_wet
            xf = sum(poly, sp.zeros(3, 1))/len(poly)   # centroid of the square face
            face_moment += A_wet*(xf.dot(nf))
            wet_polys.append(wet)
            if len(on_plane) == 2:
                cap_edges.append(on_plane)
    # the cap polygon: chain the clip edges
    cap = []
    if cap_edges:
        chain = list(cap_edges[0])
        edges = cap_edges[1:]
        while edges:
            last = chain[-1]
            for k, e in enumerate(edges):
                if e[0] == last:
                    chain.append(e[1]); edges.pop(k); break
                if e[1] == last:
                    chain.append(e[0]); edges.pop(k); break
            else:
                break
        if len(chain) >= 4 and chain[0] == chain[-1]:
            chain.pop()
        cap = chain
    if len(cap) < 3:
        return None     # the plane misses the cube or touches it in a point/edge
    cap_va = polygon_vector_area(cap)
    # (i) normal identity: m equals the cap's vector area (up to the chain orientation)
    ok_normal = (m - cap_va == sp.zeros(3, 1)) or (m + cap_va == sp.zeros(3, 1))
    ok_parallel = cross(m, n) == sp.zeros(3, 1) and m.dot(n) > 0
    # (ii) position identity: with |n| irrational in general, verify it in the
    # rational form  (3V - face_moment) (n.n) == s (m.n),  since p = s/|n| and
    # |A| = (m.n)/|n|.
    V = liquid_polyhedron_volume(wet_polys + [cap])
    lhs = (3*V - face_moment)*n.dot(n)
    rhs = s*m.dot(n)
    ok_position = sp.simplify(lhs - rhs) == 0
    # the 2020 sign gives -p
    ok_sign = sp.simplify((face_moment - 3*V)*n.dot(n) + s*m.dot(n)) == 0
    return ok_normal and ok_parallel, ok_position, ok_sign, len(cap)


def part_c(n_random=150, seed=20260929):
    rng = random.Random(seed)
    cases = [
        (sp.Matrix([1, 0, 0]), R(1, 3)),               # slab, axis-aligned
        (sp.Matrix([1, 1, 0]), R(1, 4)),               # edge cut
        (sp.Matrix([1, 1, 1]), R(-3, 4)),              # corner cut
        (sp.Matrix([1, 1, 1]), R(0)),                  # hexagonal section through the centre
        (sp.Matrix([1, 2, 3]), R(1, 2)),               # general
        (sp.Matrix([2, -3, 1]), R(-1, 2)),             # general, mixed signs
        (sp.Matrix([1, 1, 0]), R(1, 2)),               # through two vertices
    ]
    for _ in range(n_random):
        n = sp.Matrix([R(rng.randint(-7, 7)), R(rng.randint(-7, 7)), R(rng.randint(-7, 7))])
        if n == sp.zeros(3, 1):
            continue
        smax = sum(abs(v) for v in n)/2
        s = R(rng.randint(-20, 20), 21)*smax
        cases.append((n, s))
    n_ok = n_bad = n_skip = 0
    shapes = {}
    for n, s in cases:
        res = check_plane(n, s)
        if res is None:
            n_skip += 1; continue
        okn, okp, oks, ncap = res
        shapes[ncap] = shapes.get(ncap, 0) + 1
        if okn and okp and oks:
            n_ok += 1
        else:
            n_bad += 1
            print(f"[C] FAIL for n = {list(n)}, s = {s}: normal {okn}, position {okp}, sign {oks}")
    print(f"[C] exact rational planes through the cube: {n_ok} ok, {n_bad} failed, {n_skip} skipped "
          f"(cap polygons by vertex count: {dict(sorted(shapes.items()))})")
    return n_bad == 0 and n_ok > 0


if __name__ == "__main__":
    ok = part_a() & part_b() & part_c()
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
