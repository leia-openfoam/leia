#!/usr/bin/env python3
"""
The joint clipping function of the swept-face flux: a prototype on the solver's
own data (Tomislav Maric, 2026-11-04; fork branch davof-swept-joint, key
sweptJointDumpSteps; the working note leia/docs/davof/davof-swept-joint).

The pieces algorithm (fluxIntegration sweptFaceCells) cuts the swept volume
Omega_f of a face into the pieces inside the donor's vertex-ring cells and
clips the piece of an interface cell by THAT cell's PLIC plane. The joint
function replaces the split by cell faces with ONE function over Omega_f,

    F(x) = min_c phi_c(x),   phi_c(x) = -(x - x0_c).n_c   (liquid where phi_c > 0),

over the PLIC planes of the interface cells the swept volume touches (the
intersection of their liquid half-spaces, exact for the supporting planes of a
convex liquid body), or the max (the union, a concave body). The liquid
volume is Omega_f clipped by F > 0: sequential half-space clips of the whole
volume by k planes, no split by cell faces.

Evaluated per face from the dump (postProcessing/sweptJointDump/step_*.txt):

  V_exact   the liquid volume of the traced swept volume under the exact sphere
            at t^n (in-code, the exact sphere-tetrahedron overlap of Strobl,
            Formella, Poeschel 2016, doi:10.1016/j.jcp.2016.02.003)
  V_pieces  the pieces' liquid volume (in-code, unscaled: the geometric volume)
  V_min     Omega_f ∩ (∩_c H_c), the joint min        (this script)
  V_max     Omega_f ∩ (∪_c H_c), the joint max        (this script)
  V_curv    V_min or V_max by the sign of the planes' dihedral opening
  V_band    the joint min inside the interface cells' pieces, the full cells'
            pieces wet, the empty ones dry (the state-guarded joint function:
            Tomislav's third guard)
  the same for the face fraction of the pre-image (the fan clipped by F),
  the jump volume V_max - V_min (the region where the planes disagree),
  k (the planes that cut Omega_f) and k_caps (the planes that bound the
  final intersection: the cap closure of the pieces holds for k_caps <= 3),
  per donor the sum of the outgoing joint liquid against alpha_c V_c.

The geometric evaluator is independent of the fork's cap closure: every
boundary triangle T spans the signed tetrahedron (x_f, T); each tetrahedron is
clipped as a convex polyhedron by the half-spaces in turn (Sutherland-Hodgman
on its faces, I. E. Sutherland, G. W. Hodgman, Reentrant polygon clipping,
Commun. ACM 17 (1974) 32-42, doi:10.1145/360767.360802; the cap polygon from
the crossing points, ordered by angle, convex) and its volume is the
divergence-theorem sum (1/3) sum_faces p.a over the clipped faces. The signed
sum over the tetrahedra is the volume of Omega_f ∩ (∩ H) for any closed
oriented list, as in swept_volume.py of the davof-swept-volume note. The
exact face fraction of a planar fan triangle is the area of its intersection
with the ball, the disc of the sphere's trace in the triangle's plane, by the
signed circle-triangle areas about the disc's centre.

Pure Python (no numpy): tiny vectors are faster as tuples, and the cluster's
system python has no dependency to check. Usage:

    joint_prototype.py --dump step_1_0.txt --out prefix [--max-faces N]
                       [--sample-check 200] [--co 1] [--N 32] [--row exact]

Writes prefix_faces.csv (one row per sweeping face), prefix_donors.csv (the
boundedness check per donor), prefix_summary.csv (one row), and prints the
summary.
"""
import argparse
import csv
import math
import sys
import time

# --------------------------------------------------------------------------- vectors

def sub(a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
def add(a, b): return (a[0] + b[0], a[1] + b[1], a[2] + b[2])
def mul(a, s): return (a[0]*s, a[1]*s, a[2]*s)
def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def cross(a, b): return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])
def norm(a): return math.sqrt(dot(a, a))


def newell(poly):
    """The area vector (1/2) sum p_j x p_{j+1} of a polygon."""
    ax = ay = az = 0.0
    m = len(poly)
    for j in range(m):
        p = poly[j]; q = poly[(j + 1) % m]
        ax += p[1]*q[2] - p[2]*q[1]
        ay += p[2]*q[0] - p[0]*q[2]
        az += p[0]*q[1] - p[1]*q[0]
    return (0.5*ax, 0.5*ay, 0.5*az)


# --------------------------------------------------------------------------- polygon and polyhedron clips

def clip_polygon(poly, n, x0):
    """The part of a convex polygon with (x - x0).n <= 0 (Sutherland-Hodgman, one plane).
    Returns (polygon, crossed)."""
    m = len(poly)
    s = [dot(sub(p, x0), n) for p in poly]
    anyIn = False; anyOut = False
    for v in s:
        if v <= 0.0: anyIn = True
        else: anyOut = True
    if not anyOut: return poly, False
    if not anyIn: return [], False
    out = []
    for j in range(m):
        p = poly[j]; q = poly[(j + 1) % m]; sp = s[j]; sq = s[(j + 1) % m]
        if sp <= 0.0:
            out.append(p)
            if sq > 0.0:
                t = sp/(sp - sq); out.append(add(p, mul(sub(q, p), t)))
        elif sq <= 0.0:
            t = sp/(sp - sq); out.append(add(p, mul(sub(q, p), t)))
    return (out if len(out) >= 3 else []), True


def clip_polyhedron(faces, n, x0, pid):
    """Clip a convex polyhedron (a list of (polygon, label)) by (x - x0).n <= 0.
    The cap polygon (label pid) is built from the crossing points, ordered by
    angle about their mean in the plane (convex). Returns the new face list."""
    new_faces = []
    cap = []
    for poly, lab in faces:
        s = [dot(sub(p, x0), n) for p in poly]
        anyIn = False; anyOut = False
        for v in s:
            if v <= 0.0: anyIn = True
            else: anyOut = True
        if not anyOut:
            new_faces.append((poly, lab)); continue
        if not anyIn:
            continue
        m = len(poly); out = []
        for j in range(m):
            p = poly[j]; q = poly[(j + 1) % m]; sp = s[j]; sq = s[(j + 1) % m]
            if sp <= 0.0:
                out.append(p)
                if sq > 0.0:
                    t = sp/(sp - sq); x = add(p, mul(sub(q, p), t)); out.append(x); cap.append(x)
            elif sq <= 0.0:
                t = sp/(sp - sq); x = add(p, mul(sub(q, p), t)); out.append(x); cap.append(x)
        if len(out) >= 3:
            new_faces.append((out, lab))
    if len(cap) >= 3:
        cx = sum(p[0] for p in cap)/len(cap); cy = sum(p[1] for p in cap)/len(cap); cz = sum(p[2] for p in cap)/len(cap)
        cen = (cx, cy, cz)
        e1 = sub(cap[0], cen); e1 = sub(e1, mul(n, dot(e1, n)))
        L = norm(e1)
        if L > 1e-300:
            e1 = mul(e1, 1.0/L); e2 = cross(n, e1)
            ang = sorted((math.atan2(dot(sub(p, cen), e2), dot(sub(p, cen), e1)), i) for i, p in enumerate(cap))
            ordered = [cap[i] for _, i in ang]
            dedup = [ordered[0]]
            for p in ordered[1:]:
                if norm(sub(p, dedup[-1])) > 1e-15:
                    dedup.append(p)
            if len(dedup) >= 3 and norm(sub(dedup[0], dedup[-1])) <= 1e-15:
                dedup.pop()
            if len(dedup) >= 3:
                a = newell(dedup)
                if dot(a, n) < 0.0:
                    dedup.reverse()
                new_faces.append((dedup, pid))
    return new_faces


def polyhedron_volume(faces):
    v = 0.0
    for poly, _ in faces:
        v += dot(poly[0], newell(poly))
    return v/3.0


def tet_faces(o, a, b, c):
    """The outward faces of the positively oriented tetrahedron (o, a, b, c)."""
    return [([o, b, a], -1), ([o, c, b], -1), ([o, a, c], -1), ([a, b, c], -1)]


def list_clipped_volume(T, xref, planes, want_caps=False):
    """The volume of Omega (the closed outward list T) ∩ (∩ half-spaces), by the
    signed tetrahedra (xref, T_k) clipped in turn. planes: list of (n, x0).
    Returns (volume, set of plane ids with a final cap of non-zero area)."""
    vol = 0.0
    caps = set()
    for (a, b, c) in T:
        tp = dot(cross(sub(a, xref), sub(b, xref)), sub(c, xref))
        if abs(tp) <= 1e-300:
            continue
        if tp > 0:
            sign = 1.0; faces = tet_faces(xref, a, b, c)
        else:
            sign = -1.0; faces = tet_faces(xref, a, c, b)
        L2 = max(dot(sub(a, xref), sub(a, xref)), dot(sub(b, xref), sub(b, xref)), dot(sub(c, xref), sub(c, xref)))
        for pid, (n, x0) in enumerate(planes):
            faces = clip_polyhedron(faces, n, x0, pid)
            if not faces:
                break
        if faces:
            vol += sign*polyhedron_volume(faces)
            if want_caps:
                for poly, lab in faces:
                    if lab >= 0 and norm(newell(poly)) > 1e-14*L2:
                        caps.add(lab)
    return vol, caps


def fan_clipped_area(Tpre, planes):
    """The area vector of the pre-image fan clipped by the half-spaces."""
    a = (0.0, 0.0, 0.0)
    for tri in Tpre:
        poly = list(tri)
        for n, x0 in planes:
            poly, _ = clip_polygon(poly, n, x0)
            if not poly:
                break
        if poly:
            a = add(a, newell(poly))
    return a


# --------------------------------------------------------------------------- the exact disc area of a planar triangle

def circle_segment_area(c, r, p, q):
    """The signed area of the intersection of the disc (c, r) with the triangle
    (c, p, q) in 2D: the segment p -> q split at the circle crossings; a part
    inside the disc contributes its triangle, a part outside its sector."""
    def tri(a, b):
        return 0.5*((a[0] - c[0])*(b[1] - c[1]) - (a[1] - c[1])*(b[0] - c[0]))

    def sector(a, b):
        th = math.atan2((a[0] - c[0])*(b[1] - c[1]) - (a[1] - c[1])*(b[0] - c[0]),
                        (a[0] - c[0])*(b[0] - c[0]) + (a[1] - c[1])*(b[1] - c[1]))
        return 0.5*r*r*th

    d = (q[0] - p[0], q[1] - p[1])
    f = (p[0] - c[0], p[1] - c[1])
    A = d[0]*d[0] + d[1]*d[1]
    if A <= 1e-300:
        return 0.0
    B = 2.0*(f[0]*d[0] + f[1]*d[1])
    C = f[0]*f[0] + f[1]*f[1] - r*r
    disc = B*B - 4.0*A*C
    ts = [0.0]
    if disc > 0.0:
        sq = math.sqrt(disc)
        for t in ((-B - sq)/(2.0*A), (-B + sq)/(2.0*A)):
            if 0.0 < t < 1.0:
                ts.append(t)
    ts.append(1.0)
    ts.sort()
    area = 0.0
    for i in range(len(ts) - 1):
        t0, t1 = ts[i], ts[i + 1]
        if t1 - t0 <= 0.0:
            continue
        a = (p[0] + t0*d[0], p[1] + t0*d[1]); b = (p[0] + t1*d[0], p[1] + t1*d[1])
        tm = 0.5*(t0 + t1); m = (p[0] + tm*d[0] - c[0], p[1] + tm*d[1] - c[1])
        if m[0]*m[0] + m[1]*m[1] <= r*r:
            area += tri(a, b)
        else:
            area += sector(a, b)
    return area


def triangle_ball_area(tri, o, R):
    """The area of the planar triangle inside the ball (o, R), and the unit
    normal of the triangle (the area vector's direction)."""
    a, b, c = tri
    av = cross(sub(b, a), sub(c, a)); A = 0.5*norm(av)
    if A <= 1e-300:
        return 0.0, (0.0, 0.0, 0.0)
    m = mul(av, 1.0/(2.0*A))
    d = dot(sub(o, a), m)
    if abs(d) >= R:
        return 0.0, m
    r = math.sqrt(R*R - d*d)
    cen = sub(o, mul(m, d))
    e1 = sub(b, a); e1 = mul(e1, 1.0/norm(e1)); e2 = cross(m, e1)
    def to2(p):
        v = sub(p, a); return (dot(v, e1), dot(v, e2))
    c2 = to2(cen); p2 = [to2(a), to2(b), to2(c)]
    s = 0.0
    for j in range(3):
        s += circle_segment_area(c2, r, p2[j], p2[(j + 1) % 3])
    # the polygon's orientation in (e1, e2) is positive by construction (b - a along e1,
    # c on the e2 side), so s is the inside area
    return max(0.0, min(abs(s), A)), m


# --------------------------------------------------------------------------- the dump

class Face:
    __slots__ = ("facei", "donor", "phi", "sweeps", "vTotal", "vPieces", "vUsed", "afPieces", "afDonor",
                 "fallback", "vExact", "h", "F", "preSign", "Apre", "xf", "xq", "p", "q", "cells")


class Cell:
    __slots__ = ("cell", "cls", "isDonor", "alpha", "wetPlane", "wet", "n", "x0", "vol", "caps", "ok",
                 "cov", "covPre", "aw", "planes")


def read_dump(path, max_faces=None):
    header = {}
    faces = []
    with open(path) as fh:
        line = fh.readline()
        tok = line.split()
        # "# sweptJointDump step idx tn t dt dt sphere s o x y z R R rank r"
        header["step"] = int(tok[3]); header["tn"] = float(tok[5]); header["dt"] = float(tok[7])
        header["sphere"] = int(tok[9]); header["o"] = (float(tok[11]), float(tok[12]), float(tok[13]))
        header["R"] = float(tok[15])
        for line in fh:
            if line.startswith("#"):
                continue
            tok = line.split()
            if tok[0] != "face":
                raise RuntimeError("unexpected line: " + line[:60])
            f = Face()
            f.facei = int(tok[1]); f.donor = int(tok[2]); f.phi = float(tok[3]); f.sweeps = int(tok[4])
            f.vTotal = float(tok[5]); f.vPieces = float(tok[6]); f.vUsed = float(tok[7])
            f.afPieces = float(tok[8]); f.afDonor = float(tok[9]); f.fallback = int(tok[10])
            f.vExact = float(tok[11]); f.h = float(tok[12]); f.F = int(tok[13]); f.preSign = int(tok[14])
            f.Apre = (float(tok[15]), float(tok[16]), float(tok[17]))
            t = fh.readline().split(); f.xf = (float(t[1]), float(t[2]), float(t[3]))
            t = fh.readline().split(); f.xq = (float(t[1]), float(t[2]), float(t[3]))
            f.p = []
            for _ in range(f.F):
                t = fh.readline().split(); f.p.append((float(t[1]), float(t[2]), float(t[3])))
            f.q = []
            for _ in range(f.F):
                t = fh.readline().split(); f.q.append((float(t[1]), float(t[2]), float(t[3])))
            t = fh.readline().split(); nc = int(t[1])
            f.cells = []
            for _ in range(nc):
                t = fh.readline().split()
                c = Cell()
                c.cell = int(t[1]); c.cls = t[2]; c.isDonor = int(t[3]); c.alpha = float(t[4])
                c.wetPlane = int(t[5]); c.wet = int(t[6])
                c.n = (float(t[7]), float(t[8]), float(t[9])); c.x0 = (float(t[10]), float(t[11]), float(t[12]))
                c.vol = float(t[13]); c.caps = int(t[14]); c.ok = int(t[15]); c.cov = float(t[16]); c.covPre = float(t[17])
                c.aw = (float(t[18]), float(t[19]), float(t[20]))
                npl = int(t[21]); c.planes = []
                for k in range(npl):
                    b = 22 + 6*k
                    c.planes.append(((float(t[b]), float(t[b + 1]), float(t[b + 2])), (float(t[b + 3]), float(t[b + 4]), float(t[b + 5]))))
                f.cells.append(c)
            faces.append(f)
            if max_faces and len(faces) >= max_faces:
                break
    return header, faces


def triangle_list(f):
    """The closed list as sweptTriangles builds it (the file header of the dump)."""
    T = []
    F = f.F
    for i in range(F):
        a = f.p[i]; b = f.p[(i + 1) % F]; aq = f.q[i]; bq = f.q[(i + 1) % F]
        T.append((f.xf, a, b)); T.append((f.xq, bq, aq)); T.append((a, aq, bq)); T.append((a, bq, b))
    if f.preSign == 1:
        T = [(t[0], t[2], t[1]) for t in T]
    return T


def preimage_fan(f):
    F = f.F
    return [(f.xq, f.q[i], f.q[(i + 1) % F]) for i in range(F)]


# --------------------------------------------------------------------------- the evaluation of one face

def dihedral_sign(planes):
    """The opening of the plane pairs toward the liquid: with n OUT of the liquid
    and the plane points x0, a convex liquid body has (n_j - n_i).(x0_j - x0_i) > 0
    for every pair (the outward normal turns with the position). Returns
    (+1 convex, -1 concave, 0 saddle/mixed, None for one plane)."""
    k = len(planes)
    if k < 2:
        return None
    pos = neg = 0
    for i in range(k):
        for j in range(i + 1, k):
            s = dot(sub(planes[j][0], planes[i][0]), sub(planes[j][1], planes[i][1]))
            if s > 1e-14: pos += 1
            elif s < -1e-14: neg += 1
    if pos and not neg: return 1
    if neg and not pos: return -1
    if pos and neg: return 0
    return 1   # parallel planes (no opening): take the min


def evaluate_face(f, hdr, check_pieces=False):
    T = triangle_list(f)
    Tpre = preimage_fan(f)
    res = {}
    res["facei"] = f.facei; res["donor"] = f.donor; res["phi"] = f.phi; res["fallback"] = f.fallback
    res["vTotal"] = f.vTotal; res["vExact"] = f.vExact; res["h"] = f.h
    vol_list = sum(dot(t[0], newell(list(t))) for t in T)/3.0
    res["vList"] = vol_list
    # the pieces' geometric liquid volume: in-code; a fallback face carries the
    # one-donor clip scaled, undone here
    scale = abs(f.phi)*hdr["dt"]/f.vTotal if f.vTotal > 0 else 1.0
    res["scale"] = scale
    res["vPieces"] = f.vPieces if not f.fallback else (f.vUsed/scale if scale > 0 else f.vUsed)
    # the plane set: the interface cells (a PLIC plane, or a wisp's matched
    # plane) whose piece has coverage (the cells the swept volume touches)
    planes = []
    plane_cells = []
    nFull = nEmpty = nWisp = nInt = nIntAll = 0
    for c in f.cells:
        if c.wetPlane:
            nIntAll += 1
        if c.cov <= 0.0:
            continue
        if c.wetPlane:
            planes.append((c.n, c.x0)); plane_cells.append(c); nInt += 1
        elif c.cls == "full": nFull += 1
        elif c.cls == "empty": nEmpty += 1
        else: nWisp += 1
    res["k"] = len(planes); res["kAll"] = nIntAll
    res["nFull"] = nFull; res["nEmpty"] = nEmpty; res["nWisp"] = nWisp
    # the planes that cut the swept volume at all (both signs among its vertices)
    verts = [f.xf, f.xq] + f.p + f.q
    kcut = 0
    for n, x0 in planes:
        s = [dot(sub(v, x0), n) for v in verts]
        if min(s) < 0.0 < max(s):
            kcut += 1
    res["kCut"] = kcut
    if not planes:
        # no plane at all: the swept volume lies in full / empty / wisp cells
        res["vMin"] = res["vMax"] = res["vBand"] = res["vCurv"] = res["vPieces"]
        res["kCaps"] = 0; res["curv"] = ""; res["afMin"] = res["afMax"] = res["afBand"] = f.afPieces
        res["afPieces"] = f.afPieces; res["afDonor"] = f.afDonor
        res["afExact"] = exact_face_fraction(f, Tpre, hdr)
        res["vJump"] = 0.0
        return res
    vMin, caps = list_clipped_volume(T, f.xf, planes, want_caps=True)
    res["vMin"] = vMin; res["kCaps"] = len(caps)
    # the max: the complement of the intersection of the complements
    planesC = [(mul(n, -1.0), x0) for n, x0 in planes]
    vDry, _ = list_clipped_volume(T, f.xf, planesC)
    res["vMax"] = vol_list - vDry
    res["vJump"] = res["vMax"] - vMin
    sgn = dihedral_sign(planes)
    res["curv"] = {1: "convex", -1: "concave", 0: "saddle", None: "one"}[sgn]
    res["vCurv"] = vMin if (sgn is None or sgn >= 0) else res["vMax"]
    # the state-guarded joint function: inside every interface cell's piece
    # the joint min, the full cells' pieces wet, the empty ones dry, a wisp by
    # its rule (the in-code piece volume times its wet flag)
    vBand = 0.0
    for c in f.cells:
        if c.cov <= 0.0:
            continue
        if c.wetPlane:
            v, _ = list_clipped_volume(T, f.xf, c.planes + planes)
            vBand += v
        elif c.wet:
            vBand += c.vol
    res["vBand"] = vBand
    if check_pieces:
        # the pieces reproduced here (each interface piece by its own plane)
        vp = 0.0
        for c in f.cells:
            if c.cov <= 0.0:
                continue
            if c.wetPlane:
                v, _ = list_clipped_volume(T, f.xf, c.planes + [(c.n, c.x0)])
                vp += v
            elif c.wet:
                v, _ = list_clipped_volume(T, f.xf, c.planes)
                vp += v
        res["vPiecesPy"] = vp
    # the face fraction: the pre-image fan clipped by the joint function
    A2 = dot(f.Apre, f.Apre)
    def frac(a):
        return dot(a, f.Apre)/A2 if A2 > 0 else 0.0
    res["afPieces"] = f.afPieces; res["afDonor"] = f.afDonor
    res["afMin"] = frac(fan_clipped_area(Tpre, planes))
    aDry = fan_clipped_area(Tpre, planesC)
    res["afMax"] = frac(sub(f.Apre, aDry))
    # the band variant of the face update: per cell as the pieces, the joint min
    # inside the interface cells
    aBand = (0.0, 0.0, 0.0)
    for c in f.cells:
        if c.covPre <= 0.0:
            continue
        if c.wetPlane:
            aBand = add(aBand, fan_clipped_area(Tpre, c.planes + planes))
        elif c.wet:
            aBand = add(aBand, fan_clipped_area(Tpre, c.planes))
    res["afBand"] = frac(aBand)
    res["afExact"] = exact_face_fraction(f, Tpre, hdr)
    return res


def exact_face_fraction(f, Tpre, hdr):
    if not hdr["sphere"]:
        return float("nan")
    A2 = dot(f.Apre, f.Apre)
    if A2 <= 0:
        return 0.0
    a = (0.0, 0.0, 0.0)
    for tri in Tpre:
        area, m = triangle_ball_area(tri, hdr["o"], hdr["R"])
        a = add(a, mul(m, area))
    return dot(a, f.Apre)/A2


# --------------------------------------------------------------------------- the run

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-faces", type=int, default=0)
    ap.add_argument("--sample-check", type=int, default=200, help="faces on which the pieces are reproduced here")
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    t0 = time.time()
    hdr, faces = read_dump(args.dump, args.max_faces or None)
    print("read %d faces, step %d, tn %.6g, dt %.6g, sphere %d o %s R %g" % (len(faces), hdr["step"], hdr["tn"], hdr["dt"], hdr["sphere"], hdr["o"], hdr["R"]))
    sweeping = [f for f in faces if f.sweeps == 1]
    rows = []
    check_every = max(1, len(sweeping)//max(1, args.sample_check))
    for i, f in enumerate(sweeping):
        rows.append(evaluate_face(f, hdr, check_pieces=(i % check_every == 0)))
        if (i + 1) % 1000 == 0:
            print("  %d / %d faces, %.0f s" % (i + 1, len(sweeping), time.time() - t0)); sys.stdout.flush()
    cols = ["facei", "donor", "phi", "fallback", "h", "vTotal", "vList", "scale", "vExact", "vPieces", "vPiecesPy",
            "vMin", "vMax", "vCurv", "vBand", "vJump", "k", "kAll", "kCut", "kCaps", "nFull", "nEmpty", "nWisp", "curv",
            "afExact", "afPieces", "afDonor", "afMin", "afMax", "afBand"]
    with open(args.out + "_faces.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])

    # the summary: L1 over faces against the truth, the worst faces, the jump
    # region, k, the boundedness per donor
    def l1(key):
        return sum(abs(r[key] - r["vExact"]) for r in rows if r["vExact"] >= 0)
    def lmax(key):
        return max((abs(r[key] - r["vExact"]) for r in rows if r["vExact"] >= 0), default=0.0)
    nE = sum(1 for r in rows if r["vExact"] >= 0)
    sumExact = sum(r["vExact"] for r in rows if r["vExact"] >= 0)
    sumTotal = sum(r["vTotal"] for r in rows)
    summary = {}
    summary["label"] = args.label; summary["dump"] = args.dump; summary["step"] = hdr["step"]; summary["dt"] = hdr["dt"]
    summary["h"] = rows[0]["h"] if rows else 0
    summary["faces"] = len(sweeping); summary["facesExact"] = nE
    summary["fallbacks"] = sum(r["fallback"] for r in rows)
    summary["sumVexact"] = sumExact; summary["sumVtotal"] = sumTotal
    for key in ("vPieces", "vMin", "vMax", "vCurv", "vBand"):
        summary["L1_" + key] = l1(key); summary["max_" + key] = lmax(key)
        summary["net_" + key] = sum(r[key] - r["vExact"] for r in rows if r["vExact"] >= 0)
    summary["L1_vMin_vPieces"] = sum(abs(r["vMin"] - r["vPieces"]) for r in rows)
    summary["L1_vBand_vPieces"] = sum(abs(r["vBand"] - r["vPieces"]) for r in rows)
    summary["sumJump"] = sum(r["vJump"] for r in rows)
    summary["sumJumpRel"] = summary["sumJump"]/sumExact if sumExact > 0 else float("nan")
    summary["maxJumpRel"] = max((r["vJump"]/r["vTotal"] for r in rows if r["vTotal"] > 0), default=0.0)
    chk = [(abs(r["vPiecesPy"] - r["vPieces"]), r["vTotal"]) for r in rows if "vPiecesPy" in r]
    summary["piecesCheckFaces"] = len(chk)
    summary["piecesCheckMaxRel"] = max((d/v for d, v in chk if v > 0), default=0.0)
    chkL = [abs(r["vList"] - r["vTotal"])/r["vTotal"] for r in rows if r["vTotal"] > 0]
    summary["listCheckMaxRel"] = max(chkL, default=0.0)
    for kk in range(0, 8):
        summary["k%d" % kk] = sum(1 for r in rows if r["k"] == kk)
    summary["kGt3"] = sum(1 for r in rows if r["k"] > 3)
    summary["kCapsGt3"] = sum(1 for r in rows if r["kCaps"] > 3)
    summary["kCutGt3"] = sum(1 for r in rows if r["kCut"] > 3)
    summary["kCapsLe3Share"] = (1.0 - summary["kCapsGt3"]/len(rows)) if rows else float("nan")
    for cv in ("convex", "concave", "saddle", "one", ""):
        summary["curv_" + (cv or "none")] = sum(1 for r in rows if r["curv"] == cv)
    summary["facesWithFull"] = sum(1 for r in rows if r["nFull"] > 0)
    summary["facesWithEmpty"] = sum(1 for r in rows if r["nEmpty"] > 0)
    summary["facesWithWisp"] = sum(1 for r in rows if r["nWisp"] > 0)
    # the face fractions
    ok = [r for r in rows if r["afExact"] == r["afExact"]]
    for key in ("afPieces", "afDonor", "afMin", "afMax", "afBand"):
        summary["L1_" + key] = sum(abs(r[key] - r["afExact"]) for r in ok)/max(1, len(ok))
        summary["max_" + key] = max((abs(r[key] - r["afExact"]) for r in ok), default=0.0)
    # the clip counts: the pieces today (every ring cell: 7 planes x 4F
    # triangles... counted as plane-triangle clips), the pieces after the
    # pre-test (the cells with coverage), the joint (k planes on the whole list)
    nRing = sum(len(f.cells) for f in sweeping)
    nCov = sum(1 for f in sweeping for c in f.cells if c.cov > 0)
    nCovPlane = sum(1 for f in sweeping for c in f.cells if c.cov > 0 and c.wetPlane)
    nTri = sum(4*f.F for f in sweeping)
    summary["ringCellsPerFace"] = nRing/len(sweeping) if sweeping else 0
    summary["piecesPerFace"] = nCov/len(sweeping) if sweeping else 0
    summary["clipsPiecesNow"] = sum(4*f.F*(len(c.planes) + (1 if c.wetPlane else 0)) for f in sweeping for c in f.cells)
    summary["clipsPiecesPretest"] = sum(4*f.F*(len(c.planes) + (1 if c.wetPlane else 0)) for f in sweeping for c in f.cells if c.cov > 0)
    summary["clipsJoint"] = sum(4*f.F*r["k"] for f, r in zip(sweeping, rows))
    summary["clipsBand"] = sum(4*f.F*(len(c.planes) + r["k"]) for f, r in zip(sweeping, rows) for c in f.cells if c.cov > 0 and c.wetPlane) \
        + sum(4*f.F*len(c.planes) for f, r in zip(sweeping, rows) for c in f.cells if c.cov > 0 and not c.wetPlane)
    summary["seconds"] = time.time() - t0

    # the boundedness per donor: the outgoing liquid (scaled as the flux) of
    # the joint variants against alpha_c V_c
    donors = {}
    for f, r in zip(sweeping, rows):
        alpha = None
        for c in f.cells:
            if c.isDonor:
                alpha = c.alpha
        if alpha is None:
            continue
        d = donors.setdefault(f.donor, {"alpha": alpha, "V": f.h**3, "n": 0, "pieces": 0.0, "min": 0.0, "max": 0.0, "band": 0.0, "exact": 0.0})
        d["n"] += 1
        for key, src in (("pieces", "vPieces"), ("min", "vMin"), ("max", "vMax"), ("band", "vBand"), ("exact", "vExact")):
            d[key] += r["scale"]*max(r[src], 0.0)
    with open(args.out + "_donors.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["donor", "alpha", "V", "nFaces", "outPieces", "outMin", "outMax", "outBand", "outExact",
                    "overPieces", "overMin", "overMax", "overBand"])
        for c, d in sorted(donors.items()):
            aV = d["alpha"]*d["V"]
            w.writerow([c, d["alpha"], d["V"], d["n"], d["pieces"], d["min"], d["max"], d["band"], d["exact"],
                        (d["pieces"] - aV)/d["V"], (d["min"] - aV)/d["V"], (d["max"] - aV)/d["V"], (d["band"] - aV)/d["V"]])
    for key in ("pieces", "min", "max", "band"):
        over = [(d[key] - d["alpha"]*d["V"])/d["V"] for d in donors.values()]
        summary["overdrawMax_" + key] = max(over, default=0.0)
        summary["overdrawCount_" + key] = sum(1 for o in over if o > 1e-12)
    summary["donors"] = len(donors)

    with open(args.out + "_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(list(summary.keys())); w.writerow(list(summary.values()))
    # the worst faces of the joint min against the truth
    worst = sorted(rows, key=lambda r: -abs(r["vMin"] - r["vExact"]))[:20]
    with open(args.out + "_worst.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in worst:
            w.writerow([r.get(c, "") for c in cols])
    for k, v in summary.items():
        print("%-24s %s" % (k, v))


if __name__ == "__main__":
    main()
