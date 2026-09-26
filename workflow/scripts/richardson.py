#!/usr/bin/env python3
"""Observed orders, Richardson extrapolation and the Grid Convergence Index (GCI) for
the method gates (CLAUDE.md "Method gates").

Two kinds of quantity:

* An ERROR with an exact reference (shape error, gradient band error, volume error):
  the pairwise orders p_ij = ln(e_i/e_j)/ln(h_i/h_j) between consecutive rungs and the
  least-squares slope of ln e against ln h over all rungs. The error itself is the
  measurement; no extrapolation.

* A QUANTITY without an exact value (oscillation period, damping rate): the procedure of
  Celik et al., "Procedure for estimation and reporting of uncertainty due to
  discretization in CFD applications", J. Fluids Eng. 130, 078001 (2008), with the
  Roache (1994) safety factor Fs = 1.25. Three rungs h1 < h2 < h3 (fine to coarse),
  r21 = h2/h1, r32 = h3/h2, eps21 = f2 - f1, eps32 = f3 - f2, s = sign(eps32/eps21):

      p   = |ln|eps32/eps21| + q(p)| / ln r21,  q(p) = ln((r21^p - s)/(r32^p - s))
      f_ext = (r21^p f1 - f2)/(r21^p - 1)
      GCI_fine = Fs |(f1 - f2)/f1| / (r21^p - 1)

  solved by fixed-point iteration, so a non-constant refinement ratio (integer N on a
  sqrt(2) or 1.3 ladder) is handled exactly. The convergence type comes from
  R = eps21/eps32: monotone for 0 < R < 1, oscillatory for -1 < R < 0, divergent for
  |R| > 1. The asymptotic-range ratio GCI_32/(r21^p GCI_21) is about 1 inside the
  asymptotic range.

`--self-test` checks every formula on synthetic data and exits nonzero on a failure.
"""
import math
import sys

FS = 1.25


def _sorted_fine_to_coarse(h, f):
    pairs = sorted(zip(h, f), key=lambda x: x[0])
    return [p[0] for p in pairs], [p[1] for p in pairs]


def pairwise_orders(h, e):
    """Orders between consecutive rungs, coarse to fine: [(h_coarse, h_fine, p), ...].
    p is None where an error is not positive."""
    hs, es = _sorted_fine_to_coarse(h, e)
    hs, es = hs[::-1], es[::-1]          # coarse to fine
    out = []
    for i in range(1, len(hs)):
        a, b = es[i - 1], es[i]
        p = (math.log(a / b) / math.log(hs[i - 1] / hs[i])
             if (a is not None and b is not None and a > 0 and b > 0) else None)
        out.append((hs[i - 1], hs[i], p))
    return out


def lsq_order(h, e):
    """Least-squares slope of ln e against ln h; None with fewer than two positive errors."""
    pts = [(math.log(x), math.log(y)) for x, y in zip(h, e)
           if x and x > 0 and y is not None and y > 0]
    if len(pts) < 2:
        return None
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    if sxx == 0:
        return None
    return sum((p[0] - mx) * (p[1] - my) for p in pts) / sxx


def celik(h, f, max_iter=200, tol=1e-13):
    """Celik et al. (2008) three-rung procedure. Returns a dict with p, f_ext, gci_fine,
    gci_coarse, asymptotic_ratio, conv_type, R. Needs exactly three rungs."""
    if len(h) != 3 or len(f) != 3 or any(x is None for x in f):
        return {"conv_type": "insufficient"}
    (h1, h2, h3), (f1, f2, f3) = _sorted_fine_to_coarse(h, f)
    r21, r32 = h2 / h1, h3 / h2
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 == 0:
        return {"conv_type": "degenerate", "R": None}
    R = e21 / e32
    conv = ("monotone" if 0 < R < 1 else
            "oscillatory" if -1 < R < 0 else "divergent")
    s = 1.0 if e32 / e21 > 0 else -1.0
    p = abs(math.log(abs(e32 / e21))) / math.log(r21)
    for _ in range(max_iter):
        try:
            q = math.log((r21 ** p - s) / (r32 ** p - s))
        except (ValueError, ZeroDivisionError):
            return {"conv_type": conv, "R": R, "p": None}
        p_new = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        if abs(p_new - p) < tol:
            p = p_new
            break
        p = p_new
    rp21, rp32 = r21 ** p, r32 ** p
    f_ext = (rp21 * f1 - f2) / (rp21 - 1.0)
    gci_fine = FS * abs((f1 - f2) / f1) / (rp21 - 1.0) if f1 else None
    gci_coarse = FS * abs((f2 - f3) / f2) / (rp32 - 1.0) if f2 else None
    ratio = (gci_coarse / (rp21 * gci_fine)) if (gci_fine and gci_coarse) else None
    return {"p": p, "f_ext": f_ext, "gci_fine": gci_fine, "gci_coarse": gci_coarse,
            "asymptotic_ratio": ratio, "conv_type": conv, "R": R,
            "r21": r21, "r32": r32}


def _self_test():
    fails = []

    def check(name, got, want, rtol):
        ok = got is not None and abs(got - want) <= rtol * max(abs(want), 1e-300)
        print(f"  {'PASS' if ok else 'FAIL'} {name}: got {got!r}, want {want!r}")
        if not ok:
            fails.append(name)

    # Integer ladders with non-constant ratios, as the gates use them.
    ladders = {"2D droplets": [100, 142, 200], "2D shear": [68, 96, 136],
               "3D droplets": [60, 78, 102], "3D shear": [68, 90, 118]}
    for lname, Ns in ladders.items():
        h = [1.0 / n for n in Ns]
        for p_true in (1.0, 2.0, 3.0):
            f0, C = 0.3, 5.0
            f = [f0 + C * x ** p_true for x in h]
            res = celik(h, f)
            check(f"{lname} celik p (p={p_true:g})", res["p"], p_true, 1e-8)
            check(f"{lname} celik f_ext (p={p_true:g})", res["f_ext"], f0, 1e-10)
            if res["conv_type"] != "monotone":
                print(f"  FAIL {lname} conv_type (p={p_true:g}): {res['conv_type']}")
                fails.append(f"{lname} conv_type {res['conv_type']}")
            e = [C * x ** p_true for x in h]
            check(f"{lname} lsq order (p={p_true:g})", lsq_order(h, e), p_true, 1e-10)
            for hc, hf, po in pairwise_orders(h, e):
                check(f"{lname} pairwise order {hc:.4g}->{hf:.4g} (p={p_true:g})", po, p_true, 1e-10)
    # Oscillatory and divergent sequences must be classified, not extrapolated blindly.
    # h and f are listed as PAIRS, fine to coarse: R = eps21/eps32.
    h = [1.0 / 200, 1.0 / 142, 1.0 / 100]
    for fvals, want in (([1.0, 1.1, 0.95], "oscillatory"),   # R = 0.1/-0.15
                        ([1.0, 1.2, 1.3], "divergent"),      # R = 0.2/0.1
                        ([1.0, 1.1, 1.3], "monotone")):      # R = 0.1/0.2
        got = celik(h, fvals)["conv_type"]
        print(f"  {'PASS' if got == want else 'FAIL'} classification {fvals}: got {got}, want {want}")
        if got != want:
            fails.append(f"classification {fvals}")
    # Asymptotic-range ratio: GCI_32/(r21^p GCI_21) = f1/f2 exactly for f = f0 + C h^p, so it
    # tends to 1 only as the relative discretization error becomes small.
    Ns = (68, 96, 136)
    f = [0.3 + 5 * (1 / n) ** 2 for n in Ns]
    res = celik([1 / n for n in Ns], f)
    check("asymptotic ratio = f1/f2", res["asymptotic_ratio"], f[2] / f[1], 1e-8)   # f1 = finest (N = 136)
    print(f"[richardson] self-test: {len(fails)} failure(s)" + (f": {fails}" if fails else ""))
    return 1 if fails else 0


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        sys.exit(_self_test())
    print(__doc__)
