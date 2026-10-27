#!/usr/bin/env python3
"""
One-ring quadric (Monge) fit on analytically sampled surfaces: what converges
under which error model.  Replacement study for Table 3 of the DAVOF proposal
(``\\label{tab:svd-convergence}``).

Diagnosis of the lost Table 3
-----------------------------
The proposal's Table 3 reports a local quadratic fit on a randomly oriented,
analytically sampled ellipsoid with additive Gaussian noise sigma = h^2 over
five levels h = 1.25e-2 ... 7.81e-4, with RMSE_fit 1.18e-4 ... 5.72e-7 and
RMSE_kappa 2.60e-2 ... 3.12e-5, rates about 2.  The script that produced it is
lost; its patch size and sample count were never recorded.  Second order of the
curvature under i.i.d. noise of standard deviation sigma is only possible if
the patch has a FIXED PHYSICAL half-width L: the standard deviation of a fitted
second derivative on n samples spread over L scales like
sigma * sqrt(45/n) / L^2, so with sigma = h^2 and L fixed the curvature error is
O(h^2) (and even O(h^2) * sqrt(h^2/L^2) when n grows like L^2/h^2), whereas with
L = 1.5 h (one ring of cells) it is O(1).  The method described in WP2 fits the
quadric through ONE RING of neighbouring cells, so its patch width is O(h); the
method's own one-ring fit on the ellipsoid converged at order 1.55 with
pairwise rates 1.83, 1.49, 1.26 (falling toward first order).  The lost table
therefore measured a fixed-size patch, not the method.  This script measures
the one-ring fit in cells, on a quadric (ellipsoid) and a non-quadric (torus)
surface, under explicit error models.

Design
------
* Surfaces (mm): ellipsoid with half-axes 1.0, 0.8, 0.6; torus with major
  radius R = 1 and minor radius r = 0.4 (kappa = 1/r + cos(phi)/(R + r cos(phi)),
  the Gaussian curvature changes sign).  Implicit functions with exact normal,
  exact total curvature kappa = kappa_1 + kappa_2 = div(n) (2/R on a sphere)
  and the closest-point projection (closed form on the torus, iterative on the
  ellipsoid).
* Sampling = one ring in cells.  Cartesian grid of spacing h = 4 mm / N,
  N in (20, 40, 80, 160, 320), in a 4 mm box centred on the surface.  A cell is
  cut when the implicit function changes sign among its eight vertices.  One
  sample per cut cell: the closest-point projection onto the surface of the
  cell centre displaced by a uniform random offset of up to 0.3 h in each
  coordinate (a PLIC-piece centroid is not the cell centre; a symmetric
  Cartesian stencil cancels the cubic term and would fake second order on the
  torus; ``--symmetric`` uses the undisplaced centre as a demonstration).  Each
  realisation applies a random rotation (uniform in SO(3)) and a random grid
  offset (uniform in [-h/2, h/2]^3) to the surface.  A patch around a centre
  cut cell = the samples of the cut cells in the 3x3x3 block of cells around it
  ("one ring"; the actual count is recorded) and, with ``--rings 2``, the 5x5x5
  block.  20 realisations; patch centres: a random subsample of 4000 cut cells
  per (surface, level, realisation), all of them with ``--all-centres``.
* Fit: a Monge patch z = a0 + a1 x + a2 y + a3 x^2 + a4 x y + a5 y^2 in the
  local frame whose z axis is the EXACT normal at the patch centre (the role of
  the method's DAVOF normal) and whose origin is the exact (unperturbed) centre
  sample.  Points only by default; ``--normals`` adds, per sample, the two
  tangential components of the exact sample normal as equations for the
  gradient of f (weighted by h, so that the scales match).  The coordinates are
  scaled by h before the normal-equation solve.  From the fit: the distance of
  the fitted point (0, 0, a0) to the exact surface (closest-point projection),
  and the total curvature at the centre from the general Monge formula with the
  fitted gradient, kappa = -[(1+q^2) r - 2 p q s + (1+p^2) t] / (1+p^2+q^2)^{3/2},
  p = a1, q = a2, r = 2 a3, s = a4, t = 2 a5 (sign such that a sphere seen from
  outside has kappa = +2/R, as the implicit formula).  Compared with the exact
  kappa at the exact centre point.
* Error metrics per (surface, level, error model, rings): RMSE, L1 (mean
  absolute) and Linf of the position error [mm]; RMSE, L1, Linf of
  |kappa_fit - kappa_exact| / RMS(kappa_exact over the centres) (the
  normalisation of the proposal's Table 1).  Observed orders: least-squares
  slope of log(error) against log(h) over all levels and over the three finest
  levels, and the pairwise rates between consecutive levels.
* Error models, applied to the sample positions along the exact normal at
  each sample (units: h in mm, sigma = (h / 1 mm)^2 mm, so sigma = 0.04 mm at
  h = 0.2 mm and 1.56e-4 mm at h = 0.0125 mm):
  E0 none;
  E1 i.i.d. Gaussian with standard deviation sigma = h^2;
  E2 a smooth correlated field: 64 random Fourier modes in 3D with
     |k| = 2 pi / L_c, L_c = 4 h, random directions and phases, scaled to
     RMS h^2, evaluated at the sample points;
  E3 a smooth field with FIXED correlation length L_c = 0.5 mm and RMS h^2 plus
     i.i.d. jitter of standard deviation c h^3 / (1 mm)^2 with c = 1 (no
     calibration; reported in its own column set).

PRE-REGISTERED PREDICTIONS (written before the first run; one ring, Monge,
points only, torus unless stated)
---------------------------------
P1 (E0, torus and ellipsoid): curvature error first order (aliasing of the
   cubic term on an asymmetric stencil); position error third order.
P2 (E1): curvature error O(1), about 2 sigma/h^2 = 2 1/mm, i.e. about 60 % of
   the torus curvature RMS, order about 0; two rings about six times smaller,
   same order; position error second order.
P3 (E2, L_c = 4h): curvature error constant about h^2/L_c^2 = 1/16 1/mm, order
   about 0; position error second order.
P4 (E3): curvature order 2 while the smooth part dominates, falling to 1 at the
   finest levels; position order 2.
P5 (``--symmetric``, E0, torus): curvature second order (the cubic term
   cancels), demonstrating why the asymmetric stencil is the honest case.
P6 (``--normals``, E0, ellipsoid): errors reduced by a constant factor, orders
   unchanged.

Usage
-----
    python3 qsvd_onering_convergence.py                      # default study
    python3 qsvd_onering_convergence.py --symmetric --surfaces torus \
        --models E0 --tag table3_symmetric                  # P5
    python3 qsvd_onering_convergence.py --normals --surfaces ellipsoid \
        --models E0 --tag table3_normals                    # P6
    python3 qsvd_onering_convergence.py --selftest

Outputs (in ``--outdir``, prefixed by ``--tag``): ``<tag>.csv`` (long format),
``<tag>_orders.csv``, ``<tag>.tex`` (tabular only, default tag), and
``<tag>_convergence.png``.
"""

import argparse
import csv
import itertools
import math
import os
import sys
import time

import numpy as np

SURFACE_IDS = {"torus": 0, "ellipsoid": 1}
MODEL_IDS = {"E0": 0, "E1": 1, "E2": 2, "E3": 3}
BOX = 4.0           # mm, edge of the Cartesian box
N_MODES = 64        # Fourier modes of the smooth error fields
LC_E3 = 0.5         # mm, fixed correlation length of E3
JITTER_C = 1.0      # E3 jitter: std = c h^3 / (1 mm)^2


# --------------------------------------------------------------------------
# surfaces (in their own frame, lengths in mm)
# --------------------------------------------------------------------------
def unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def implicit_total_curvature(g, H):
    """kappa = div(grad F / |grad F|) = (|g|^2 tr H - g^T H g) / |g|^3."""
    g2 = np.einsum("mi,mi->m", g, g)
    trH = np.einsum("mii->m", H)
    gHg = np.einsum("mi,mij,mj->m", g, H, g)
    return (g2 * trH - gHg) / g2 ** 1.5


class Ellipsoid:
    name = "ellipsoid"

    def __init__(self, a=1.0, b=0.8, c=0.6):
        self.axes = np.array([a, b, c])
        self.inv2 = 1.0 / self.axes ** 2

    def f(self, q):
        return np.einsum("mi,i->m", q * q, self.inv2) - 1.0

    def grad(self, q):
        return 2.0 * q * self.inv2

    def hess(self, q):
        H = np.zeros(q.shape + (3,))
        H[:, 0, 0] = 2.0 * self.inv2[0]
        H[:, 1, 1] = 2.0 * self.inv2[1]
        H[:, 2, 2] = 2.0 * self.inv2[2]
        return H

    def kappa(self, q):
        return implicit_total_curvature(self.grad(q), self.hess(q))

    def closest(self, p, steps=100):
        """Exact foot point.  Lagrange condition q_i = p_i a_i^2 / (a_i^2 + mu)
        with G(mu) = sum p_i^2 a_i^2 / (a_i^2 + mu)^2 - 1 = 0; G is monotone on
        (-c^2, inf) (c the smallest half-axis) and its root there is the closest
        point, inside and outside (Eberly).  Bisection on a bracket that holds
        the root by construction, then a Newton polish."""
        a2 = self.axes ** 2
        imin = int(np.argmin(a2))
        c2 = a2[imin]
        pa2 = p * p * a2
        lo = -c2 + np.abs(p[:, imin]) * math.sqrt(c2) / math.sqrt(2.0)   # G(lo) >= 1
        hi = np.sqrt(pa2.sum(-1)) - c2                                  # G(hi) <= 0

        def G(mu):
            return (pa2 / (a2 + mu[:, None]) ** 2).sum(-1) - 1.0
        for _ in range(steps):
            mid = 0.5 * (lo + hi)
            up = G(mid) > 0
            lo = np.where(up, mid, lo)
            hi = np.where(up, hi, mid)
        mu = 0.5 * (lo + hi)
        for _ in range(3):
            dG = -2.0 * (pa2 / (a2 + mu[:, None]) ** 3).sum(-1)
            mu = mu - G(mu) / dG
        return p * a2 / (a2 + mu[:, None])


class Torus:
    name = "torus"

    def __init__(self, R=1.0, r=0.4):
        self.R = R
        self.r = r

    def f(self, q):
        rho = np.hypot(q[:, 0], q[:, 1])
        return (rho - self.R) ** 2 + q[:, 2] ** 2 - self.r ** 2

    def grad(self, q):
        rho = np.hypot(q[:, 0], q[:, 1])
        g = rho - self.R
        return 2.0 * np.stack([g * q[:, 0] / rho, g * q[:, 1] / rho, q[:, 2]], -1)

    def hess(self, q):
        x, y = q[:, 0], q[:, 1]
        rho = np.hypot(x, y)
        g = rho - self.R
        H = np.zeros(q.shape + (3,))
        H[:, 0, 0] = 2.0 * (x * x / rho ** 2 + g * y * y / rho ** 3)
        H[:, 1, 1] = 2.0 * (y * y / rho ** 2 + g * x * x / rho ** 3)
        H[:, 0, 1] = H[:, 1, 0] = 2.0 * (x * y / rho ** 2 - g * x * y / rho ** 3)
        H[:, 2, 2] = 2.0
        return H

    def kappa(self, q):
        return implicit_total_curvature(self.grad(q), self.hess(q))

    def kappa_closed_form(self, q):
        """1/r + cos(phi)/(R + r cos(phi)) at surface points, outward normal."""
        rho = np.hypot(q[:, 0], q[:, 1])
        cphi = (rho - self.R) / self.r
        return 1.0 / self.r + cphi / (self.R + self.r * cphi)

    def closest(self, p):
        rho = np.hypot(p[:, 0], p[:, 1])
        c = np.stack([self.R * p[:, 0] / rho, self.R * p[:, 1] / rho,
                      np.zeros_like(rho)], -1)
        d = p - c
        return c + self.r * unit(d)


def random_rotation(rng):
    """Uniform random rotation from a normalised Gaussian quaternion."""
    w, x, y, z = rng.normal(size=4) / 1.0
    s = 1.0 / math.sqrt(w * w + x * x + y * y + z * z)
    w, x, y, z = w * s, x * s, y * s, z * s
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


class Placed:
    """A surface rotated by R and shifted by t: world p = R q + t."""

    def __init__(self, surf, R, t):
        self.surf, self.R, self.t = surf, R, t

    def to_local(self, p):
        return (p - self.t) @ self.R

    def to_world(self, q):
        return q @ self.R.T + self.t

    def f(self, p):
        return self.surf.f(self.to_local(p))

    def normal(self, p):
        return unit(self.surf.grad(self.to_local(p)) @ self.R.T)

    def kappa(self, p):
        return self.surf.kappa(self.to_local(p))

    def closest(self, p):
        return self.to_world(self.surf.closest(self.to_local(p)))

    def distance(self, p):
        return np.linalg.norm(p - self.closest(p), axis=-1)


# --------------------------------------------------------------------------
# cut cells and samples
# --------------------------------------------------------------------------
def cut_cells(placed, N, h):
    """Indices (i, j, k) of the cells whose eight vertices do not share the
    sign of the implicit function; lexicographically sorted."""
    x0 = -BOX / 2
    v = x0 + h * np.arange(N + 1)
    sgn = np.empty((N + 1,) * 3, np.int8)
    slab = max(1, int(2.0e6 // ((N + 1) * (N + 1))))
    for i0 in range(0, N + 1, slab):
        i1 = min(N + 1, i0 + slab)
        X, Y, Z = np.meshgrid(v[i0:i1], v, v, indexing="ij")
        P = np.stack([X.ravel(), Y.ravel(), Z.ravel()], -1)
        sgn[i0:i1] = np.sign(placed.f(P)).reshape(i1 - i0, N + 1, N + 1)
    pos = np.zeros((N, N, N), bool)
    neg = np.zeros((N, N, N), bool)
    for dx, dy, dz in itertools.product((0, 1), repeat=3):
        s = sgn[dx:dx + N, dy:dy + N, dz:dz + N]
        pos |= s >= 0
        neg |= s <= 0
    return np.argwhere(pos & neg).astype(np.int64)


def tangent_basis(n):
    """Orthonormal (t1, t2) perpendicular to the unit normals n (M, 3)."""
    a = np.zeros_like(n)
    idx = np.argmin(np.abs(n), axis=1)
    a[np.arange(len(n)), idx] = 1.0
    t1 = unit(np.cross(n, a))
    t2 = np.cross(n, t1)
    return t1, t2


def fourier_field(rng, k, rms):
    """Superposition of N_MODES random plane waves with |k| fixed, RMS = rms."""
    u = unit(rng.normal(size=(N_MODES, 3)))
    kv = k * u
    theta = rng.uniform(0.0, 2.0 * np.pi, N_MODES)
    amp = rms / math.sqrt(N_MODES / 2.0)

    def phi(x):
        return amp * np.cos(x @ kv.T + theta).sum(-1)
    return phi


def perturb(model, pts, nrm, h, rng):
    """Displace the samples along their exact normals (units of mm, h in mm)."""
    sigma = h * h
    if model == "E0":
        return pts
    if model == "E1":
        d = rng.normal(size=len(pts)) * sigma
    elif model == "E2":
        d = fourier_field(rng, 2.0 * np.pi / (4.0 * h), sigma)(pts)
    elif model == "E3":
        d = fourier_field(rng, 2.0 * np.pi / LC_E3, sigma)(pts)
        d = d + rng.normal(size=len(pts)) * JITTER_C * h ** 3
    else:
        raise ValueError(model)
    return pts + d[:, None] * nrm


def patch_neighbours(ijk, lin, N, cidx, rings):
    """For each centre, the sample index of every cut cell in the
    (2 rings + 1)^3 block and a mask of which block cells are cut."""
    r = np.arange(-rings, rings + 1)
    off = np.array(list(itertools.product(r, r, r)), dtype=np.int64)
    nb = ijk[cidx][:, None, :] + off[None]
    valid = ((nb >= 0) & (nb < N)).all(-1)
    nbl = (nb[..., 0] * N + nb[..., 1]) * N + nb[..., 2]
    pos = np.searchsorted(lin, nbl)
    pos = np.minimum(pos, len(lin) - 1)
    found = valid & (lin[pos] == nbl)
    return pos, found


def fit_patches(pts, nrm, h, x0, n0, nbidx, mask, use_normals, cond_max=1e10):
    """Batched Monge fit.  pts: (possibly perturbed) sample positions;
    nrm: exact sample normals; x0, n0: exact centre point and normal.
    Returns the fitted point, the fitted curvature, the ok flag and the
    number of samples per patch."""
    t1, t2 = tangent_basis(n0)
    d = pts[nbidx] - x0[:, None, :]
    xi = np.einsum("pmi,pi->pm", d, t1) / h
    eta = np.einsum("pmi,pi->pm", d, t2) / h
    zeta = np.einsum("pmi,pi->pm", d, n0) / h
    w = mask.astype(float)
    one = np.ones_like(xi)
    zero = np.zeros_like(xi)
    A = np.stack([one, xi, eta, xi * xi, xi * eta, eta * eta], -1) * w[..., None]
    b = zeta * w
    if use_normals:
        m = nrm[nbidx]
        m1 = np.einsum("pmi,pi->pm", m, t1)
        m2 = np.einsum("pmi,pi->pm", m, t2)
        m3 = np.einsum("pmi,pi->pm", m, n0)
        Ax = np.stack([zero, one, zero, 2 * xi, eta, zero], -1) * w[..., None]
        Ay = np.stack([zero, zero, one, zero, xi, 2 * eta], -1) * w[..., None]
        A = np.concatenate([A, Ax, Ay], axis=1)
        b = np.concatenate([b, -m1 / m3 * w, -m2 / m3 * w], axis=1)
    M = np.einsum("pmi,pmj->pij", A, A)
    rhs = np.einsum("pmi,pm->pi", A, b)
    n_s = mask.sum(1)
    with np.errstate(all="ignore"):
        cond = np.linalg.cond(M)
    ok = (n_s >= 6) & np.isfinite(cond) & (cond < cond_max)
    coef = np.full((len(x0), 6), np.nan)
    if ok.any():
        coef[ok] = np.linalg.solve(M[ok], rhs[ok][..., None])[..., 0]
    a0 = h * coef[:, 0]
    a1, a2 = coef[:, 1], coef[:, 2]
    a3, a4, a5 = coef[:, 3] / h, coef[:, 4] / h, coef[:, 5] / h
    W2 = 1.0 + a1 * a1 + a2 * a2
    kappa_fit = -((1 + a2 * a2) * 2 * a3 - 2 * a1 * a2 * a4
                  + (1 + a1 * a1) * 2 * a5) / W2 ** 1.5
    pfit = x0 + a0[:, None] * n0
    return pfit, kappa_fit, ok, n_s


# --------------------------------------------------------------------------
# study
# --------------------------------------------------------------------------
def make_surface(name):
    return Torus() if name == "torus" else Ellipsoid()


def run_study(args):
    levels = list(args.levels)
    rings_list = list(args.rings)
    models = list(args.models)
    results = {}   # (surface, N, rings, model) -> dict of lists
    t_start = time.time()
    for sname in args.surfaces:
        surf = make_surface(sname)
        sid = SURFACE_IDS[sname]
        for N in levels:
            h = BOX / N
            t_lev = time.time()
            for rr in range(args.realisations):
                rng = np.random.default_rng([args.seed, sid, N, rr])
                R = random_rotation(rng)
                t = rng.uniform(-h / 2, h / 2, 3)
                placed = Placed(surf, R, t)
                ijk = cut_cells(placed, N, h)
                lin = (ijk[:, 0] * N + ijk[:, 1]) * N + ijk[:, 2]
                centres = -BOX / 2 + h * (ijk + 0.5)
                if not args.symmetric:
                    centres = centres + rng.uniform(-0.3 * h, 0.3 * h, centres.shape)
                pts = placed.closest(centres)
                nrm = placed.normal(pts)
                kap = placed.kappa(pts)
                ncut = len(ijk)
                if args.all_centres or ncut <= args.centres:
                    cidx = np.arange(ncut)
                else:
                    cidx = np.sort(rng.choice(ncut, args.centres, replace=False))
                x0 = pts[cidx]
                n0 = nrm[cidx]
                k0 = kap[cidx]
                neigh = {rg: patch_neighbours(ijk, lin, N, cidx, rg) for rg in rings_list}
                for model in models:
                    rng_m = np.random.default_rng([args.seed, sid, N, rr, MODEL_IDS[model]])
                    ppts = perturb(model, pts, nrm, h, rng_m)
                    for rg in rings_list:
                        nbidx, mask = neigh[rg]
                        pfit, kfit, ok, n_s = fit_patches(
                            ppts, nrm, h, x0, n0, nbidx, mask, args.normals)
                        epos = np.full(len(cidx), np.nan)
                        epos[ok] = placed.distance(pfit[ok])
                        key = (sname, N, rg, model)
                        res = results.setdefault(key, dict(
                            epos=[], ekap=[], kex=[], n_s=[], ok=[], ncut=[]))
                        res["epos"].append(epos[ok])
                        res["ekap"].append(np.abs(kfit[ok] - k0[ok]))
                        res["kex"].append(k0)
                        res["n_s"].append(n_s)
                        res["ok"].append(ok)
                        res["ncut"].append(ncut)
            print(f"  {sname:9s} N={N:4d} h={h:.4e}  cut cells={ncut:7d}  "
                  f"{time.time() - t_lev:6.1f} s  (total {time.time() - t_start:6.1f} s)",
                  flush=True)
    return results


METRIC_NAMES = ("RMSE", "L1", "Linf")


def norms(e):
    return {"RMSE": float(np.sqrt(np.mean(e * e))),
            "L1": float(np.mean(np.abs(e))),
            "Linf": float(np.max(np.abs(e)))}


def summarise(results):
    """Per key: metrics of the position and the normalised curvature error."""
    out = {}
    for key, res in results.items():
        epos = np.concatenate(res["epos"])
        ekap = np.concatenate(res["ekap"])
        kex = np.concatenate(res["kex"])
        n_s = np.concatenate(res["n_s"])
        ok = np.concatenate(res["ok"])
        krms = float(np.sqrt(np.mean(kex * kex)))
        out[key] = dict(
            position=norms(epos),
            curvature=norms(ekap / krms),
            kappa_rms=krms,
            n_patches=int(ok.size),
            frac_excluded=float(1.0 - ok.mean()),
            n_samples_mean=float(n_s.mean()),
            n_samples_min=int(n_s.min()),
            n_samples_max=int(n_s.max()),
            n_cut_mean=float(np.mean(res["ncut"])),
        )
    return out


def slope(h, e):
    h, e = np.asarray(h, float), np.asarray(e, float)
    good = np.isfinite(e) & (e > 0)
    if good.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(h[good]), np.log(e[good]), 1)[0])


def orders(summary, surfaces, levels, rings_list, models):
    rows = []
    hs = [BOX / N for N in levels]
    for sname in surfaces:
        for rg in rings_list:
            for model in models:
                keys = [(sname, N, rg, model) for N in levels]
                if not all(k in summary for k in keys):
                    continue
                for qty in ("position", "curvature"):
                    for nm in METRIC_NAMES:
                        e = [summary[k][qty][nm] for k in keys]
                        rates = [math.log(e[i] / e[i + 1]) / math.log(hs[i] / hs[i + 1])
                                 if e[i] > 0 and e[i + 1] > 0 else float("nan")
                                 for i in range(len(e) - 1)]
                        rows.append(dict(
                            surface=sname, rings=rg, model=model, quantity=qty, norm=nm,
                            order_all=slope(hs, e), order_fine3=slope(hs[-3:], e[-3:]),
                            rates=rates, errors=e))
    return rows


# --------------------------------------------------------------------------
# outputs
# --------------------------------------------------------------------------
def write_long_csv(path, summary, args):
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["surface", "N", "h_mm", "rings", "model", "normals", "symmetric",
                     "n_patches", "n_cut_cells_mean", "n_samples_mean", "n_samples_min",
                     "n_samples_max", "frac_excluded", "kappa_rms_exact", "quantity",
                     "norm", "value"])
        for (sname, N, rg, model), s in sorted(summary.items(), key=lambda kv: (
                kv[0][0], kv[0][1], kv[0][2], kv[0][3])):
            for qty in ("position", "curvature"):
                for nm in METRIC_NAMES:
                    wr.writerow([sname, N, f"{BOX / N:.6e}", rg, model, int(args.normals),
                                 int(args.symmetric), s["n_patches"],
                                 f"{s['n_cut_mean']:.1f}", f"{s['n_samples_mean']:.3f}",
                                 s["n_samples_min"], s["n_samples_max"],
                                 f"{s['frac_excluded']:.3e}", f"{s['kappa_rms']:.6e}",
                                 qty, nm, f"{s[qty][nm]:.6e}"])


def load_long_csv(path):
    """Inverse of write_long_csv, for --replot."""
    summary = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            key = (row["surface"], int(row["N"]), int(row["rings"]), row["model"])
            s = summary.setdefault(key, dict(
                position={}, curvature={}, kappa_rms=float(row["kappa_rms_exact"]),
                n_patches=int(row["n_patches"]), frac_excluded=float(row["frac_excluded"]),
                n_samples_mean=float(row["n_samples_mean"]),
                n_samples_min=int(row["n_samples_min"]),
                n_samples_max=int(row["n_samples_max"]),
                n_cut_mean=float(row["n_cut_cells_mean"]),
                normals=bool(int(row["normals"])), symmetric=bool(int(row["symmetric"]))))
            s[row["quantity"]][row["norm"]] = float(row["value"])
    return summary


def write_orders_csv(path, rows, levels, args):
    pair_names = [f"rate_{levels[i]}_{levels[i + 1]}" for i in range(len(levels) - 1)]
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["surface", "rings", "model", "normals", "symmetric", "quantity",
                     "norm", "order_all", "order_fine3"] + pair_names
                    + [f"err_N{N}" for N in levels])
        for r in rows:
            wr.writerow([r["surface"], r["rings"], r["model"], int(args.normals),
                         int(args.symmetric), r["quantity"], r["norm"],
                         f"{r['order_all']:.3f}", f"{r['order_fine3']:.3f}"]
                        + [f"{x:.3f}" for x in r["rates"]]
                        + [f"{x:.6e}" for x in r["errors"]])


def sci(x):
    """$1.25\\times10^{-2}$ in the style of the existing Table 3."""
    if not np.isfinite(x) or x == 0:
        return "--"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(m, 2) >= 10.0:
        m /= 10.0
        e += 1
    return f"${m:.2f}\\times10^{{{e}}}$"


TEX_COLUMNS = [  # (surface, model, quantity, header)
    ("torus", "E0", "curvature", r"torus E0 $\kappa$"),
    ("torus", "E1", "curvature", r"torus E1 $\kappa$"),
    ("torus", "E2", "position", r"torus E2 pos."),
    ("torus", "E2", "curvature", r"torus E2 $\kappa$"),
    ("ellipsoid", "E2", "curvature", r"ellipsoid E2 $\kappa$"),
]


def write_tex(path, summary, order_rows, levels):
    """A complete tabular (one ring, points only, RMSE) to \\input in the float."""
    need = [(s, N, 1, m) for (s, m, _, _) in TEX_COLUMNS for N in levels]
    if not all(k in summary for k in need):
        return False
    ordmap = {(r["surface"], r["rings"], r["model"], r["quantity"], r["norm"]): r
              for r in order_rows}
    ns = {s: np.mean([summary[(s, N, 1, "E0")]["n_samples_mean"] for N in levels])
          for s in ("torus", "ellipsoid")}
    lines = []
    lines.append("% Generated by figures/qsvd/qsvd_onering_convergence.py; do not edit.")
    lines.append("% One ring (3x3x3 block of cells), points only, 20 realisations,")
    lines.append("% 4000 patch centres per realisation and level. RMSE over all patches:")
    lines.append("% position = distance of the fitted centre point to the exact surface [mm];")
    lines.append("% kappa = |kappa_fit - kappa_exact| / RMS(kappa_exact). E0 exact samples;")
    lines.append("% E1 i.i.d. Gaussian sigma = h^2; E2 smooth field, L_c = 4h, RMS h^2.")
    lines.append(f"% Mean samples per patch: torus {ns['torus']:.1f}, ellipsoid {ns['ellipsoid']:.1f}.")
    lines.append(r"\begin{tabular}{cccccc}")
    lines.append(r"\toprule")
    lines.append("$h$ [mm] & " + " & ".join(c[3] for c in TEX_COLUMNS) + r" \\")
    lines.append(r"\midrule")
    for N in levels:
        h = BOX / N
        cells = [sci(summary[(s, N, 1, m)][q]["RMSE"]) for (s, m, q, _) in TEX_COLUMNS]
        lines.append(sci(h) + " & " + " & ".join(cells) + r" \\")
    lines.append(r"\midrule")
    cells = [f"${ordmap[(s, 1, m, q, 'RMSE')]['order_fine3']:.2f}$".replace("$-0.00$", "$0.00$")
             for (s, m, q, _) in TEX_COLUMNS]
    lines.append("order (three finest) & " + " & ".join(cells) + r" \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    with open(path, "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return True


def write_plot(path, summary, surfaces, levels, rings_list, models, title_tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter

    colours = {"E0": "#2a78d6", "E1": "#eb6834", "E2": "#1baf7a", "E3": "#eda100"}
    labels = {"E0": "E0 exact samples", "E1": r"E1 i.i.d. $\sigma = h^2$",
              "E2": r"E2 smooth, $L_c = 4h$", "E3": r"E3 smooth $L_c$ = 0.5 mm + $h^3$ jitter"}
    styles = {1: dict(ls="-", marker="o"), 2: dict(ls="--", marker="s")}
    hs = np.array([BOX / N for N in levels])
    ink, grid = "#0b0b0b", "#dcdcd8"
    fig, axes = plt.subplots(2, len(surfaces), figsize=(4.6 * len(surfaces) + 0.6, 7.8),
                             squeeze=False, sharex=True)
    for j, sname in enumerate(surfaces):
        for i, (qty, ylab) in enumerate((
                ("position", "position RMSE [mm]"),
                ("curvature", r"curvature RMSE / RMS($\kappa$)"))):
            ax = axes[i, j]
            emin, emax = np.inf, 0.0
            for model in models:
                for rg in rings_list:
                    keys = [(sname, N, rg, model) for N in levels]
                    if not all(k in summary for k in keys):
                        continue
                    e = np.array([summary[k][qty]["RMSE"] for k in keys])
                    ax.plot(hs, e, color=colours[model], lw=1.6, ms=5,
                            label=f"{labels[model]}, {rg} ring" + ("s" if rg > 1 else ""),
                            **styles[rg])
                    emin, emax = min(emin, e.min()), max(emax, e.max())
            # reference slopes anchored just below the smallest error at the
            # finest level, labelled above their mid point; the y range is set
            # by the data, so the reference lines are clipped
            anchor = emin * 0.5
            for p, txt in ((1, "$h$"), (2, "$h^2$"), (3, "$h^3$")):
                ref = anchor * (hs / hs[-1]) ** p
                ax.plot(hs, ref, color="#9a9a9a", lw=0.9, ls=":", zorder=0)
                kmid = len(hs) // 2
                ax.annotate(txt, (hs[kmid], ref[kmid]), textcoords="offset points",
                            xytext=(0, 2), ha="center", va="bottom", fontsize=8,
                            color="#52514e", annotation_clip=True,
                            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none",
                                      alpha=0.85))
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_ylim(emin * 0.3, emax * 3.0)
            ax.set_xlim(hs[-1] * 0.8, hs[0] * 1.25)
            ax.set_xticks(hs)
            ax.set_xticklabels([f"{x:g}" for x in hs])
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.grid(True, which="major", color=grid, lw=0.6)
            ax.tick_params(colors=ink, labelsize=8)
            for sp in ax.spines.values():
                sp.set_color(grid)
            ax.set_ylabel(ylab, fontsize=9, color=ink)
            if i == 0:
                ax.set_title(sname, fontsize=10, color=ink)
            if i == 1:
                ax.set_xlabel("$h$ [mm]", fontsize=9, color=ink)
    handles, labs = axes[0, 0].get_legend_handles_labels()
    if not handles:
        handles, labs = axes[1, 0].get_legend_handles_labels()
    fig.legend(handles, labs, loc="lower center", ncol=2, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(f"One-ring Monge fit, {title_tag}", fontsize=10, color=ink)
    fig.tight_layout(rect=(0, 0.11, 1, 0.97))
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)


def markdown_orders(order_rows, summary, levels, rings_list):
    """A markdown table of the three-finest orders and pairwise rates (RMSE)."""
    lines = ["| surface | model | rings | samples/patch | position order (rates) | "
             "curvature order (rates) | curvature RMSE at finest |",
             "|---|---|---|---|---|---|---|"]
    bykey = {(r["surface"], r["rings"], r["model"], r["quantity"], r["norm"]): r
             for r in order_rows}
    seen = []
    for r in order_rows:
        k = (r["surface"], r["rings"], r["model"])
        if k in seen:
            continue
        seen.append(k)
        s, rg, m = k
        rp = bykey[(s, rg, m, "position", "RMSE")]
        rc = bykey[(s, rg, m, "curvature", "RMSE")]
        ns = np.mean([summary[(s, N, rg, m)]["n_samples_mean"] for N in levels])
        fmt = lambda r_: f"{r_['order_fine3']:.2f} ({', '.join(f'{x:.2f}' for x in r_['rates'])})"
        lines.append(f"| {s} | {m} | {rg} | {ns:.1f} | {fmt(rp)} | {fmt(rc)} | "
                     f"{rc['errors'][-1]:.2e} |")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------
def selftest():
    rng = np.random.default_rng(1)
    tor = Torus()
    # torus curvature: implicit formula against the closed form
    th = rng.uniform(0, 2 * np.pi, 1000)
    ph = rng.uniform(0, 2 * np.pi, 1000)
    q = np.stack([(tor.R + tor.r * np.cos(ph)) * np.cos(th),
                  (tor.R + tor.r * np.cos(ph)) * np.sin(th), tor.r * np.sin(ph)], -1)
    d1 = np.abs(tor.kappa(q) - tor.kappa_closed_form(q)).max()
    print(f"torus kappa implicit vs closed form: max diff {d1:.2e}")
    # torus closest point: foot condition
    p = q + rng.normal(size=q.shape) * 0.1
    c = tor.closest(p)
    n = unit(tor.grad(c))
    d2 = np.abs(np.cross(p - c, n)).max()
    print(f"torus closest point: |F| {np.abs(tor.f(c)).max():.1e}, tangential residual {d2:.1e}")
    # ellipsoid closest point: on the surface and (p - c) parallel to n
    ell = Ellipsoid()
    u = unit(rng.normal(size=(2000, 3)))
    q = u * ell.axes
    p = q + rng.normal(size=q.shape) * 0.15
    c = ell.closest(p)
    n = unit(ell.grad(c))
    d3 = np.abs(np.cross(p - c, n)).max()
    # brute force: no random surface point may be closer than the foot point
    th = rng.uniform(0, 2 * np.pi, 200000)
    cz = rng.uniform(-1, 1, 200000)
    sz = np.sqrt(1 - cz * cz)
    cloud = np.stack([sz * np.cos(th), sz * np.sin(th), cz], -1) * ell.axes
    dist_foot = np.linalg.norm(p - c, axis=-1)
    dist_cloud = np.array([np.linalg.norm(cloud - pi, axis=-1).min() for pi in p[:300]])
    worst = (dist_foot[:300] - dist_cloud).max()
    print(f"ellipsoid closest point: |F| {np.abs(ell.f(c)).max():.1e}, "
          f"tangential residual {d3:.1e}, foot minus best-of-cloud distance "
          f"max {worst:.1e} (must be <= 0)")
    # the fit recovers a synthetic quadratic exactly (points, and with normals)
    h = 0.05
    x0 = np.zeros((1, 3))
    n0 = np.array([[0.0, 0.0, 1.0]])
    a = np.array([0.3e-3, 0.05, -0.02, 1.2, -0.4, 0.7])  # a0 .. a5
    xy = rng.uniform(-1.5 * h, 1.5 * h, size=(1, 12, 2))
    z = (a[0] + a[1] * xy[..., 0] + a[2] * xy[..., 1] + a[3] * xy[..., 0] ** 2
         + a[4] * xy[..., 0] * xy[..., 1] + a[5] * xy[..., 1] ** 2)
    pts = np.concatenate([xy, z[..., None]], -1)[0]
    gx = a[1] + 2 * a[3] * xy[0, :, 0] + a[4] * xy[0, :, 1]
    gy = a[2] + a[4] * xy[0, :, 0] + 2 * a[5] * xy[0, :, 1]
    nrm = unit(np.stack([-gx, -gy, np.ones_like(gx)], -1))
    nbidx = np.arange(12)[None]
    mask = np.ones((1, 12), bool)
    for use_n in (False, True):
        pfit, kfit, ok, n_s = fit_patches(pts, nrm, h, x0, n0, nbidx, mask, use_n)
        p1, q1 = a[1], a[2]
        kex = -((1 + q1 * q1) * 2 * a[3] - 2 * p1 * q1 * a[4] + (1 + p1 * p1) * 2 * a[5]) \
            / (1 + p1 * p1 + q1 * q1) ** 1.5
        print(f"synthetic quadratic (normals={use_n}): a0 error {abs(pfit[0, 2] - a[0]):.1e}, "
              f"kappa error {abs(kfit[0] - kex):.1e}, ok={ok[0]}")
    # sphere, exact symmetric stencil: kappa -> 2/R at second order
    for hh in (0.1, 0.05, 0.025):
        sph = Ellipsoid(1.0, 1.0, 1.0)
        g = np.array(list(itertools.product((-1, 0, 1), repeat=2)), float) * hh
        pts = np.stack([g[:, 0], g[:, 1], np.sqrt(1 - g[:, 0] ** 2 - g[:, 1] ** 2)], -1)
        nrm = unit(pts)
        pfit, kfit, ok, n_s = fit_patches(pts, nrm, hh, np.array([[0, 0, 1.0]]),
                                          np.array([[0, 0, 1.0]]), np.arange(9)[None],
                                          np.ones((1, 9), bool), False)
        print(f"sphere symmetric 3x3, h={hh}: kappa_fit - 2 = {kfit[0] - 2:.3e}, "
              f"|a0| = {abs(pfit[0, 2] - 1):.1e}")


# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--surfaces", nargs="+", default=["torus", "ellipsoid"],
                    choices=list(SURFACE_IDS))
    ap.add_argument("--levels", nargs="+", type=int, default=[20, 40, 80, 160, 320])
    ap.add_argument("--rings", nargs="+", type=int, default=[1, 2])
    ap.add_argument("--models", nargs="+", default=["E0", "E1", "E2", "E3"],
                    choices=list(MODEL_IDS))
    ap.add_argument("--realisations", type=int, default=20)
    ap.add_argument("--centres", type=int, default=4000,
                    help="patch centres per (surface, level, realisation)")
    ap.add_argument("--all-centres", action="store_true")
    ap.add_argument("--normals", action="store_true",
                    help="add the tangential components of the exact sample normals")
    ap.add_argument("--symmetric", action="store_true",
                    help="sample the undisplaced cell centres (demonstration)")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--tag", default="table3")
    ap.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--replot", action="store_true",
                    help="redraw <tag>_convergence.png from <tag>.csv without running")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)

    if args.selftest:
        selftest()
        return 0

    if args.replot:
        summary = load_long_csv(os.path.join(args.outdir, f"{args.tag}.csv"))
        surfaces = [s for s in SURFACE_IDS if any(k[0] == s for k in summary)]
        levels = sorted({k[1] for k in summary})
        rings_list = sorted({k[2] for k in summary})
        models = [m for m in MODEL_IDS if any(k[3] == m for k in summary)]
        normals = any(summary[k]["normals"] for k in summary)
        symmetric = any(summary[k]["symmetric"] for k in summary)
        tag = ("points and normals" if normals else "points only") + \
              (", symmetric samples" if symmetric else ", displaced samples")
        p = os.path.join(args.outdir, f"{args.tag}_convergence.png")
        write_plot(p, summary, surfaces, levels, rings_list, models, tag)
        print("redrawn " + p)
        return 0

    print(f"qsvd_onering_convergence: surfaces={args.surfaces} levels={args.levels} "
          f"rings={args.rings} models={args.models} realisations={args.realisations} "
          f"centres={'all' if args.all_centres else args.centres} normals={args.normals} "
          f"symmetric={args.symmetric} seed={args.seed}", flush=True)
    t0 = time.time()
    results = run_study(args)
    summary = summarise(results)
    rows = orders(summary, args.surfaces, args.levels, args.rings, args.models)
    os.makedirs(args.outdir, exist_ok=True)
    written = []
    p = os.path.join(args.outdir, f"{args.tag}.csv")
    write_long_csv(p, summary, args)
    written.append(p)
    p = os.path.join(args.outdir, f"{args.tag}_orders.csv")
    write_orders_csv(p, rows, args.levels, args)
    written.append(p)
    if not args.normals and not args.symmetric and 1 in args.rings:
        p = os.path.join(args.outdir, f"{args.tag}.tex")
        if write_tex(p, summary, rows, args.levels):
            written.append(p)
    if not args.no_plot:
        p = os.path.join(args.outdir, f"{args.tag}_convergence.png")
        tag = ("points only" if not args.normals else "points and normals") + \
              (", symmetric samples" if args.symmetric else ", displaced samples")
        write_plot(p, summary, args.surfaces, args.levels, args.rings, args.models, tag)
        written.append(p)
    runtime = time.time() - t0
    print("\nMeasured orders (RMSE; three finest levels, pairwise rates coarse to fine):")
    print(markdown_orders(rows, summary, args.levels, args.rings))
    print("\nSample statistics per (surface, N, rings):")
    for (s, N, rg, m), v in sorted(summary.items()):
        if m != args.models[0]:
            continue
        print(f"  {s:9s} N={N:4d} rings={rg}: cut cells {v['n_cut_mean']:.0f}, patches "
              f"{v['n_patches']}, samples/patch mean {v['n_samples_mean']:.2f} "
              f"min {v['n_samples_min']} max {v['n_samples_max']}, excluded "
              f"{v['frac_excluded']:.2e}, RMS kappa {v['kappa_rms']:.3f}")
    print(f"\nRuntime {runtime:.1f} s. Files written:")
    for p in written:
        print("  " + p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
