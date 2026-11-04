#!/usr/bin/env python3
"""
Checks of joint_prototype.py's evaluators on synthetic geometry with known
answers (no solver data):

  1. the signed-tetrahedra clip with no plane returns the list's volume;
  2. one plane: the clip equals the closed form of the davof-swept-volume note,
     V = (1/3)[sum_T p_T.a_T^- - (x0.n)(n.sum_T a_T^-)], on a twisted swept
     volume (a warped pre-image);
  3. two planes: the intersection volume of an axis-aligned box with two
     half-spaces against the analytic value;
  4. the max as the complement: V_max + V(∩ complements) = V_total;
  5. the disc area of a triangle: a triangle inside the disc (its area), a
     disc inside the triangle (pi r^2), a half-disc.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from joint_prototype import (add, sub, mul, dot, cross, norm, newell, clip_polygon, list_clipped_volume,
                             fan_clipped_area, triangle_ball_area)


def swept_list(face, pre, xf, xq):
    T = []
    F = len(face)
    for i in range(F):
        a = face[i]; b = face[(i + 1) % F]; aq = pre[i]; bq = pre[(i + 1) % F]
        T.append((xf, a, b)); T.append((xq, bq, aq)); T.append((a, aq, bq)); T.append((a, bq, b))
    # outward for a face at the top (z = 1) with the pre-image below: the face fan
    # (xf, a, b) with a, b counter-clockwise seen from above points up: outward
    return T


def list_volume(T):
    return sum(dot(t[0], newell(list(t))) for t in T)/3.0


def closed_form(T, n, x0):
    total = (0.0, 0.0, 0.0); moment = 0.0
    for (a, b, c) in T:
        poly, _ = clip_polygon([a, b, c], n, x0)
        if not poly:
            continue
        aw = newell(poly)
        total = add(total, aw); moment += dot(a, aw)
    cap = -dot(n, total)
    return (moment + dot(x0, n)*cap)/3.0


def main():
    fails = 0
    def check(name, got, want, tol):
        nonlocal fails
        ok = abs(got - want) <= tol*max(1.0, abs(want))
        print("%-52s got %.15g want %.15g %s" % (name, got, want, "ok" if ok else "FAIL"))
        if not ok: fails += 1

    # a unit square face at z = 1, the pre-image warped and shifted (a twisted swept volume)
    face = [(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    pre = [(0.05, -0.02, 0.72), (1.03, 0.01, 0.66), (0.98, 1.04, 0.74), (-0.03, 0.97, 0.69)]
    xf = (0.5, 0.5, 1.0)
    xq = (0.51, 0.49, 0.70)
    T = swept_list(face, pre, xf, xq)
    V = list_volume(T)
    v0, caps = list_clipped_volume(T, xf, [], want_caps=True)
    check("1 no plane: clip = list volume", v0, V, 1e-13)
    # one tilted plane
    n = (0.3, -0.2, 0.9); L = norm(n); n = mul(n, 1.0/L); x0 = (0.4, 0.6, 0.85)
    v1, caps = list_clipped_volume(T, xf, [(n, x0)], want_caps=True)
    check("2 one plane: clip = closed form", v1, closed_form(T, n, x0), 1e-12)
    print("   caps", caps)
    # an axis-aligned box [0,1]^2 x [0.5,1] as a 'swept volume' (face at z = 1, pre-image at z = 0.5)
    pre2 = [(0, 0, 0.5), (1, 0, 0.5), (1, 1, 0.5), (0, 1, 0.5)]
    Tb = swept_list(face, pre2, xf, (0.5, 0.5, 0.5))
    check("3a box volume", list_volume(Tb), 0.5, 1e-14)
    # two planes: x <= 0.3 and y <= 0.6 -> 0.3 * 0.6 * 0.5
    p1 = ((1.0, 0.0, 0.0), (0.3, 0.0, 0.0)); p2 = ((0.0, 1.0, 0.0), (0.0, 0.6, 0.0))
    v2, caps = list_clipped_volume(Tb, xf, [p1, p2], want_caps=True)
    check("3b two planes: box ∩ {x<=0.3} ∩ {y<=0.6}", v2, 0.3*0.6*0.5, 1e-13)
    print("   caps", caps)
    # a tilted plane through the box: x + y <= 1 -> half of the box
    p3 = ((1/math.sqrt(2), 1/math.sqrt(2), 0.0), (0.5, 0.5, 0.0))
    v3, _ = list_clipped_volume(Tb, xf, [p3])
    check("3c tilted plane x + y <= 1: half the box", v3, 0.25, 1e-13)
    # 4 the max by complement on the twisted volume with two planes
    m = (-0.1, 0.4, 0.9); m = mul(m, 1.0/norm(m)); y0 = (0.6, 0.3, 0.8)
    vmin, _ = list_clipped_volume(T, xf, [(n, x0), (m, y0)])
    vdry, _ = list_clipped_volume(T, xf, [(mul(n, -1.0), x0), (mul(m, -1.0), y0)])
    vA, _ = list_clipped_volume(T, xf, [(n, x0)]); vB, _ = list_clipped_volume(T, xf, [(m, y0)])
    # inclusion-exclusion: V(A ∪ B) = V(A) + V(B) - V(A ∩ B) = V - V(A^c ∩ B^c)
    check("4 max by complement = inclusion-exclusion", V - vdry, vA + vB - vmin, 1e-12)
    # 5 the disc areas
    o = (0.5, 0.5, 0.0); R = 1.0
    tri = ((0.4, 0.4, 0.0), (0.6, 0.4, 0.0), (0.5, 0.6, 0.0))
    a, _ = triangle_ball_area(tri, o, R)
    check("5a triangle inside the disc", a, 0.5*0.2*0.2, 1e-14)
    tri = ((-5.0, -5.0, 0.0), (5.0, -5.0, 0.0), (0.0, 8.0, 0.0))
    a, _ = triangle_ball_area(tri, o, R)
    check("5b disc inside the triangle", a, math.pi, 1e-12)
    tri = ((0.5, -5.0, 0.0), (6.0, -5.0, 0.0), (0.5, 6.0, 0.0))   # the half-plane x >= 0.5 (as a big triangle)
    a, _ = triangle_ball_area(tri, o, R)
    check("5c half-disc", a, 0.5*math.pi, 1e-12)
    # a tilted triangle: the plane z = 0.6 (distance 0.6 from the centre at z = 0): r = 0.8
    tri = ((-5.0, -5.0, 0.6), (5.0, -5.0, 0.6), (0.0, 8.0, 0.6))
    a, mvec = triangle_ball_area(tri, o, R)
    check("5d disc of the plane z = 0.6 inside the triangle", a, math.pi*0.64, 1e-12)
    # the fan clip: the pre-image fan of the box by x <= 0.3: area 0.3
    Tpre = [((0.5, 0.5, 0.5), pre2[i], pre2[(i + 1) % 4]) for i in range(4)]
    aw = fan_clipped_area(Tpre, [p1])
    check("6 fan clipped by x <= 0.3", abs(aw[2]), 0.3, 1e-14)
    print("FAILS", fails)
    return fails


if __name__ == "__main__":
    sys.exit(main())
