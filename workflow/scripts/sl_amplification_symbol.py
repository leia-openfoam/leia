#!/usr/bin/env python3
"""The Fourier amplification symbol of the semi-Lagrangian quadratic fit. Research instrument.

WHY THIS EXISTS. The frozen-velocity power iteration measures rho(B) > 1 on EVERY mesh,
production hexahedra included: 1.00441 at Co = 0.0164, flat to 0.5 % across n = 32/48/64
(METHOD.md 8.2). That h-invariance is the signature of a CONSTANT-COEFFICIENT operator. On a
uniform hexahedral mesh away from boundaries the semi-Lagrangian update is exactly that, so
its symbol is available in closed form and the measured number must be reproducible with
algebra rather than with 2000 power iterations on 524 288 cells.

WHAT THIS MEASURES. The update of the anchored quadratic weighted least-squares
reconstruction (uncachedQuadraticWeightedLeastSquaresReconstruction) is

    psi^{n+1}_c = (1 - sum_j g_j) psi_c + sum_j g_j psi_j,
    g_j = w_j^2 b(d)^T M^{-1} b(d_j),   M = sum_j w_j^2 b(d_j) b(d_j)^T,

with the constant-free basis b(d) = ( d_a ; 1/2 d_a^2 ; d_a d_b, a<b ) of the code
(uncachedQuadraticWeightedLeastSquaresReconstruction.C:242-266), inverse-distance weights
w_i = 1/|x_i - x_c| of exponent 1 (.C:708-710), and the one-ring cell-point-cell stencil
(slReconstruction.C:92-99): 2 neighbours in 1D, 8 in 2D, 26 in 3D.

Writing d_j = h s_j with integer offsets s_j and d = h c, every h cancels:

    M = h^-2 S_h M_1 S_h,  S_h = diag(h ... h, h^2 ... h^2),
    g_j(c) = b(c)^T M_1^{-1} b(s_j) / |s_j|^2,   M_1 = sum_j b(s_j) b(s_j)^T / |s_j|^2.

So g_j depends ONLY on the normalised displacement c = d/h and the stencil offsets. That
DERIVES the measured h-invariance rather than observing it. Inserting psi_i = exp(i theta.i),

    G(theta) = (1 - sum_j g_j) + sum_j g_j exp(i theta.s_j),   rho = max_theta |G(theta)|.

With the cell-centre weight anchor M_1 does not depend on c, so it is inverted ONCE and every
g_j is an exact polynomial in c. With the foot-point anchor (--weight-anchor footPoint,
w_i = 1/|x_i - x_d|) the matrix does depend on c and is solved per displacement.

THE DIMENSION LADDER IS THE POINT. Coefficients: 2 in 1D, 5 in 2D, 9 in 3D, against 2 / 8 / 26
neighbours. In 1D the fit is EXACTLY DETERMINED, hence an interpolant, hence the weights
cancel -- so 1D cannot show the least-squares effect at all. The reassurance of METHOD.md 8.2
item 3 (the anchored 1D quadratic is Lax-Wendroff, Lambda = 1 + C - C^2 > 1 while
|G|^2 = 1 - 4C^2(1-C^2)sin^4(theta/2) <= 1) is therefore a 1D statement that does NOT transfer
to the over-determined 3D operator.

VALIDITY GATE. `--validate-1d` must reproduce the Lax-Wendroff weights, Lambda and |G|^2 of
METHOD.md 8.2 item 3 SYMBOLICALLY. If it does not, this instrument is wrong and no 2D or 3D
number may be read from it. The command exits non-zero on failure.

Lambda here is the Lebesgue constant |1 - sum g_j| + sum |g_j| of the same weights, reported
next to rho because METHOD.md 8.2 records that the two DISAGREE (refining pMesh raises Lambda
4.3 % while rho-1 falls 7 %). Only rho decides.

Usage:
  sl_amplification_symbol.py validate-1d
  sl_amplification_symbol.py rho --dim 3 --courant 0.0164 [--axis x|xy|xyz]
                                 [--weight-anchor cellCentre|footPoint] [--exponent 1]
  sl_amplification_symbol.py sweep  --dim 2 --dim 3 [--outdir DIR]   # the curated table + figures
  sl_amplification_symbol.py weights --dim 3 --courant 0.0164        # the (s_j, g_j) pairs
"""
import argparse
import itertools
import os
import sys

import sympy as sp

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEFAULT_OUTDIR = os.path.join(REPO, "docs", "analysis-amplification-sl", "data")

# The production measurement this instrument must reproduce (METHOD.md 8.2, the
# blockMesh hex rows). Co is matched to 3 % across the three resolutions.
MEASURED_CO = sp.Rational(164, 10000)
MEASURED_RHO = {32: 1.00439, 48: 1.00441, 64: 1.00441}
MEASURED_LAMBDA_HEX = 1.0527

AXES = {"x": (1, 0, 0), "xy": (1, 1, 0), "xyz": (1, 1, 1)}


# ---------------------------------------------------------------------------
# the operator
# ---------------------------------------------------------------------------
def basis(d, dim):
    """The constant-free quadratic basis of the code, in the code's order.

    b(d) = ( d_a ; 1/2 d_a^2 ; d_a d_b, a<b ), length 2*dim + dim*(dim-1)/2
    (uncachedQuadraticWeightedLeastSquaresReconstruction.C:242-266).
    """
    lin = [d[a] for a in range(dim)]
    sq = [sp.Rational(1, 2) * d[a] ** 2 for a in range(dim)]
    cross = [d[a] * d[b] for a in range(dim) for b in range(a + 1, dim)]
    return sp.Matrix(lin + sq + cross)


def lattice(name):
    """A stencil as (physical offsets, integer lattice offsets, dimension).

    The wave argument of Section "from weights to a growth factor" needs the SAME rule to
    repeat from cell to cell. That requires a repeating lattice; it does NOT require a
    Cartesian one. So two different sets of offsets are needed and they are not the same
    object on a non-Cartesian mesh:

      physical  -- real positions, in units of the centre-to-centre spacing. These build
                   the basis and the inverse-distance weights, i.e. the FIT.
      integer   -- coordinates in the lattice's own primitive vectors. These are what the
                   Fourier phase exp(i theta . m) is written in, i.e. the SYMBOL.

    On a Cartesian mesh the two coincide, which is why the distinction can be missed.

      cartesian1d/2d/3d  the full 3^D block minus the owner: 2, 8, 26 neighbours
      hex2d              a regular hexagonal tiling. Cell centres form a TRIANGULAR
                         lattice with primitive vectors a1 = (1,0), a2 = (1/2, sqrt3/2).
                         Each vertex of such a tiling is shared by exactly THREE cells, so
                         the vertex neighbours ARE the six edge neighbours and nothing
                         more: cell-point-cell and cell-face-cell coincide. All six sit at
                         distance one, so the inverse-distance weights are all equal and
                         the weighting has no effect at all.
    """
    if name.startswith("cartesian"):
        dim = int(name[-2])
        offs = stencil_offsets(dim)
        return [tuple(float(x) for x in o) for o in offs], offs, dim
    if name == "hex2d":
        import math
        a1 = (1.0, 0.0)
        a2 = (0.5, math.sqrt(3.0)/2.0)
        ints = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
        phys = [(m1*a1[0] + m2*a2[0], m1*a1[1] + m2*a2[1]) for m1, m2 in ints]
        return phys, ints, 2
    raise ValueError(f"unknown lattice {name!r}")


def stencil_offsets(dim):
    """One ring of the cell-point-cell stencil on a uniform hex mesh, h = 1.

    2 offsets in 1D, 8 in 2D, 26 in 3D -- the full 3^dim cube minus the owner.
    """
    return [s for s in itertools.product((-1, 0, 1), repeat=dim) if any(s)]


def normal_matrix(offsets, dim, c=None, exponent=1):
    """M_1 = sum_j w_j^2 b(s_j) b(s_j)^T, with w_j^2 the squared inverse distance.

    c is None  -> the PRODUCTION cell-centre anchor, w_j = 1/|s_j|; M_1 is then independent
                  of the displacement and is inverted once.
    c is given -> the foot-point anchor, w_j = 1/|s_j - c|.
    """
    ncoeff = 2 * dim + dim * (dim - 1) // 2
    M = sp.zeros(ncoeff, ncoeff)
    for s in offsets:
        bs = basis(sp.Matrix(s), dim)
        M += w2_of(s, c, exponent) * (bs * bs.T)
    return M


def w2_of(s, c, exponent):
    """The squared weight of stencil offset s: |r|^(-2*exponent), r measured from the anchor."""
    if c is None:
        r2 = sum(sp.Integer(si) ** 2 for si in s)
    else:
        r2 = sum((sp.Integer(si) - c[i]) ** 2 for i, si in enumerate(s))
    return sp.Pow(r2, -exponent)


def fit_weights_numeric(c, dim, anchor="cellCentre", exponent=1, phys=None):
    """g_j in float64. The fast path, and the only practical one off-axis.

    A diagonal displacement carries 1/sqrt(n), so with the foot-point anchor the entries
    of M are ALGEBRAIC rather than rational and an exact 9x9 LU solve becomes intractable.
    The symbolic path exists for the 1D validation gate, where it is cheap and where
    exactness is the point; every reported number is evaluated numerically anyway.
    """
    import numpy as np

    offsets = stencil_offsets(dim) if phys is None else phys
    cv = np.array([float(x) for x in c])
    sv = np.array(offsets, dtype=float)
    ncoeff = 2 * dim + dim * (dim - 1) // 2

    def basis_rows(pts):
        lin = [pts[:, a] for a in range(dim)]
        sq = [0.5 * pts[:, a] ** 2 for a in range(dim)]
        cross = [pts[:, a] * pts[:, b] for a in range(dim) for b in range(a + 1, dim)]
        return np.stack(lin + sq + cross, axis=1)

    B = basis_rows(sv)                                   # (nslots, ncoeff)
    ref = np.zeros(dim) if anchor == "cellCentre" else cv
    w2 = np.sum((sv - ref) ** 2, axis=1) ** (-exponent)
    Mn = (B * w2[:, None]).T @ B
    bc = basis_rows(cv.reshape(1, -1))[0]
    z = np.linalg.solve(Mn, bc)
    return offsets, list(w2 * (B @ z))


def fit_weights(c, dim, anchor="cellCentre", exponent=1, simplify=True):
    """The per-neighbour weights g_j of the update, exact.

    Returns (offsets, [g_j]). c is the normalised displacement d/h as a sympy Matrix;
    the departure foot is x_c + h c, so for a velocity U and step dt, c = -U dt / h.

    simplify=False skips the symbolic tidy-up. The VALUES are identical -- it only
    changes how the exact expression is written -- so it is used wherever the result is
    evaluated numerically. With a symbolic Courant number sp.simplify on 26 expressions
    from a 9x9 LU solve dominates the runtime.
    """
    offsets = stencil_offsets(dim)
    cw = None if anchor == "cellCentre" else c
    M = normal_matrix(offsets, dim, c=cw, exponent=exponent)
    bc = basis(c, dim)
    # z = M^{-1} b(c), solved rather than inverted: exact and much faster at 9x9.
    z = M.LUsolve(bc)
    g = [w2_of(s, cw, exponent) * (z.T * basis(sp.Matrix(s), dim))[0] for s in offsets]
    if simplify:
        g = [sp.simplify(x) for x in g]
    return offsets, g


def g_floats(g):
    """Evaluate the exact g_j to float ONCE.

    sp.N on the unsimplified sum of 26 exact rationals from a 9x9 LU solve builds a
    very large expression tree and is the dominant cost of a sweep. Every consumer
    below needs floats only, so the conversion happens here, per term.
    """
    return [float(sp.re(sp.N(x))) for x in g]


def lebesgue_f(gf):
    """Lambda from already-evaluated floats."""
    return abs(1.0 - sum(gf)) + sum(abs(x) for x in gf)


def symbol(offsets, g, theta):
    """G(theta) = (1 - sum_j g_j) + sum_j g_j exp(i theta . s_j)."""
    total = sum(g)
    return (1 - total) + sum(
        gj * sp.exp(sp.I * sum(theta[a] * s[a] for a in range(len(s))))
        for s, gj in zip(offsets, g)
    )


def lebesgue(g):
    """Lambda = |1 - sum_j g_j| + sum_j |g_j|, the row-norm bound the code writes."""
    return sp.Abs(1 - sum(g)) + sum(sp.Abs(gj) for gj in g)


# ---------------------------------------------------------------------------
# rho = max_theta |G|
# ---------------------------------------------------------------------------
def rho_numeric(offsets, g, dim, ngrid=64, refine=True):
    """max over theta in [-pi,pi]^dim of |G(theta)|, by grid then local refinement."""
    import numpy as np

    gv = np.asarray(g, dtype=float) if not isinstance(g[0], sp.Basic) \
        else np.array(g_floats(g))
    sv = np.array(offsets, dtype=float)
    const = 1.0 - gv.sum()

    axes = [np.linspace(-np.pi, np.pi, ngrid, endpoint=False) for _ in range(dim)]
    grids = np.meshgrid(*axes, indexing="ij")
    phase = sum(grids[a][..., None] * sv[None, ..., a] for a in range(dim))
    G = const + (gv * np.exp(1j * phase)).sum(axis=-1)
    mag = np.abs(G)
    flat = int(mag.argmax())
    idx = np.unravel_index(flat, mag.shape)
    best_theta = np.array([grids[a][idx] for a in range(dim)])
    best = float(mag[idx])

    if refine:
        from scipy.optimize import minimize

        def negmag(t):
            return -abs(const + (gv * np.exp(1j * (sv * t).sum(axis=1))).sum())

        res = minimize(negmag, best_theta, method="Nelder-Mead",
                       options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 20000})
        if -res.fun > best:
            best, best_theta = float(-res.fun), res.x
    return best, best_theta


def displacement(courant, dim, axis):
    """c = -Co * e, the normalised departure displacement for Courant number Co.

    The foot is x_d = x_c - U dt, so c = d/h points AGAINST the flow. For a diagonal
    direction the per-axis Courant number is Co/sqrt(naxes), so that |U| dt / h = Co.
    """
    e = AXES[axis][:dim]
    n = sum(1 for v in e if v)
    scale = sp.Rational(1, 1) if n == 1 else 1 / sp.sqrt(n)
    return sp.Matrix([-sp.nsimplify(courant) * v * scale for v in e])


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------
def cmd_validate_1d(_args):
    """The pipeline gate: 1D must be Lax-Wendroff, symbolically, with no numerics."""
    C = sp.symbols("C", positive=True)
    c = sp.Matrix([-C])                      # foot at x_c - C h
    offsets, g = fit_weights(c, 1)
    gmap = {o[0]: sp.simplify(gj) for o, gj in zip(offsets, g)}
    centre = sp.simplify(1 - sum(g))

    want = {+1: sp.simplify((C**2 - C) / 2), -1: sp.simplify((C**2 + C) / 2)}
    want_centre = sp.simplify(1 - C**2)

    ok = True
    print("1D, 3-point stencil: the update weights")
    for s in (-1, +1):
        got = sp.simplify(gmap[s] - want[s])
        good = got == 0
        ok &= good
        print(f"  g[{s:+d}] = {sp.simplify(gmap[s])}"
              f"   expected {want[s]}   {'OK' if good else 'MISMATCH'}")
    good = sp.simplify(centre - want_centre) == 0
    ok &= good
    print(f"  centre = {sp.simplify(centre)}   expected {want_centre}"
          f"   {'OK' if good else 'MISMATCH'}")

    # Lambda = 1 + C - C^2 on 0 < C < 1. The absolute values are resolved by hand there,
    # because on that interval the signs are known: the centre weight and g[-1] are
    # positive and g[+1] is negative, so |1-sum g| = 1-C^2, |g[-1]| = (C^2+C)/2 and
    # |g[+1]| = (C-C^2)/2. sympy cannot drop Abs without that assumption.
    lam = sp.simplify((1 - C**2) + (C - C**2) / 2 + (C**2 + C) / 2)
    good = sp.simplify(lam - (1 + C - C**2)) == 0
    ok &= good
    print(f"\n  Lambda = {lam}   expected {1 + C - C**2}   {'OK' if good else 'MISMATCH'}")

    # |G|^2 = 1 - 4 C^2 (1 - C^2) sin^4(theta/2)
    th = sp.symbols("theta", real=True)
    G = symbol(offsets, g, [th])
    mag2 = sp.simplify(sp.expand(sp.re(sp.expand_complex(G))**2
                                 + sp.im(sp.expand_complex(G))**2))
    want_mag2 = 1 - 4 * C**2 * (1 - C**2) * sp.sin(th / 2)**4
    diff = sp.simplify(sp.expand_trig(sp.expand(mag2 - want_mag2)))
    good = sp.simplify(diff) == 0
    ok &= good
    print(f"  |G|^2  = {sp.simplify(mag2)}")
    print(f"  expected {want_mag2}   {'OK' if good else 'MISMATCH'}")

    print("\nVALIDATE-1D: " + ("PASS" if ok else "FAIL"))
    print("1D is exactly determined (2 neighbours, 2 coefficients), so the fit is an"
          "\ninterpolant and the weights cancel. 1D cannot show the least-squares effect.")
    return 0 if ok else 1


def cmd_weights(args):
    """Print the (s_j, g_j) pairs, for comparison with the code's own g_j."""
    c = displacement(args.courant, args.dim, args.axis)
    offsets, g = fit_weights(c, args.dim, args.weight_anchor, args.exponent)
    print(f"dim={args.dim} Co={args.courant} axis={args.axis} "
          f"anchor={args.weight_anchor} exponent={args.exponent}")
    print(f"c = {[sp.nsimplify(x) for x in c]}")
    print(f"{'s_j':>14}  {'g_j (exact)':>34}  {'g_j':>14}")
    for s, gj in zip(offsets, g):
        print(f"{str(s):>14}  {str(sp.nsimplify(gj)):>34}  {float(sp.N(gj)):>14.9e}")
    print(f"\n  sum g_j      = {float(sp.N(sum(g))):.12e}")
    print(f"  1 - sum g_j  = {float(sp.N(1 - sum(g))):.12e}")
    print(f"  Lambda       = {float(sp.N(lebesgue(g))):.9f}")
    return 0


def cmd_rho(args):
    c = displacement(args.courant, args.dim, args.axis)
    offsets, gf = fit_weights_numeric(c, args.dim, args.weight_anchor, args.exponent)
    lam = lebesgue_f(gf)
    r, th = rho_numeric(offsets, gf, args.dim, ngrid=args.ngrid)
    print(f"dim={args.dim} Co={args.courant} axis={args.axis} "
          f"anchor={args.weight_anchor} exponent={args.exponent}")
    print(f"  neighbours   = {len(offsets)}")
    print(f"  coefficients = {2*args.dim + args.dim*(args.dim-1)//2}")
    print(f"  Lambda       = {lam:.6f}")
    print(f"  rho          = {r:.8f}   (rho-1 = {r-1:.3e})")
    print(f"  argmax theta = {[f'{t:+.4f}' for t in th]}  "
          f"(in units of pi: {[f'{t/3.141592653589793:+.3f}' for t in th]})")
    if args.dim == 3 and args.axis == "x" and args.weight_anchor == "cellCentre" \
            and abs(float(args.courant) - float(MEASURED_CO)) < 1e-9:
        print(f"\n  MEASURED (METHOD.md 8.2, blockMesh hex): rho = {MEASURED_RHO[64]} "
              f"at n = 64, Lambda_max = {MEASURED_LAMBDA_HEX}")
        print(f"  symbol - measured: rho {r - MEASURED_RHO[64]:+.3e}, "
              f"Lambda {lam - MEASURED_LAMBDA_HEX:+.4f}")
    return 0


def cmd_sweep(args):
    """The curated table: rho and Lambda against Courant number, per dimension and anchor."""
    import csv

    outdir = args.outdir or DEFAULT_OUTDIR
    tables = os.path.join(outdir, "tables")
    os.makedirs(tables, exist_ok=True)

    cos = [sp.Rational(k, 1000) for k in args.courants] if args.courants else \
        [MEASURED_CO] + [sp.Rational(k, 100) for k in (5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95)]
    rows = []
    for dim in args.dims:
        for anchor in args.anchors:
            for axis in args.axes:
                if dim == 1 and axis != "x":
                    continue
                if dim == 2 and axis == "xyz":
                    continue
                for co in cos:
                    c = displacement(co, dim, axis)
                    offsets, gf = fit_weights_numeric(c, dim, anchor, args.exponent)
                    lam = lebesgue_f(gf)
                    r, _ = rho_numeric(offsets, gf, dim, ngrid=args.ngrid)
                    rows.append(dict(dim=dim, anchor=anchor, axis=axis,
                                     courant=float(co), nNeighbours=len(offsets),
                                     nCoeff=2*dim + dim*(dim-1)//2,
                                     Lambda=lam, rho=r, rho_minus_1=r - 1.0))
                    print(f"dim {dim} {anchor:>10} {axis:>3} Co {float(co):.4f}  "
                          f"Lambda {lam:.6f}  rho {r:.8f}")

    path = os.path.join(tables, "sl_amplification_symbol.csv")
    with open(path, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    print(f"\nwrote {path}  ({len(rows)} rows)")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("validate-1d", help="the pipeline gate; exits non-zero on failure")

    def common(sp_):
        sp_.add_argument("--dim", type=int, default=3, choices=(1, 2, 3))
        sp_.add_argument("--courant", type=str, default=str(float(MEASURED_CO)))
        sp_.add_argument("--axis", default="x", choices=sorted(AXES))
        sp_.add_argument("--weight-anchor", default="cellCentre",
                         choices=("cellCentre", "footPoint"))
        sp_.add_argument("--exponent", type=int, default=1)
        sp_.add_argument("--ngrid", type=int, default=64)

    common(sub.add_parser("rho", help="the spectral radius at one displacement"))
    common(sub.add_parser("weights", help="the (s_j, g_j) pairs"))

    s = sub.add_parser("sweep", help="the curated table over Courant numbers")
    s.add_argument("--dims", type=int, nargs="+", default=[1, 2, 3])
    s.add_argument("--anchors", nargs="+", default=["cellCentre", "footPoint"])
    s.add_argument("--axes", nargs="+", default=["x", "xyz"])
    s.add_argument("--courants", type=int, nargs="*", default=None,
                   help="Courant numbers in thousandths, e.g. 164 250 500")
    s.add_argument("--exponent", type=int, default=1)
    s.add_argument("--ngrid", type=int, default=48)
    s.add_argument("--outdir", default=None)

    args = p.parse_args(argv)
    if args.cmd == "validate-1d":
        return cmd_validate_1d(args)
    if args.cmd == "rho":
        args.courant = sp.nsimplify(args.courant, rational=True)
        return cmd_rho(args)
    if args.cmd == "weights":
        args.courant = sp.nsimplify(args.courant, rational=True)
        return cmd_weights(args)
    if args.cmd == "sweep":
        return cmd_sweep(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())
