#!/usr/bin/env python3
"""Lambda and max|G| for BOUNDARY cells, symbolically, with the code's rank ladder
(quadratic -> linear -> constant). Tests whether the measured Lambda_max = 1.0527 and
rho = 1.00441 come from the boundary rather than the interior."""
import os, sys, itertools, numpy as np, sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "workflow", "scripts")); import sl_amplification_symbol as M

Co = sp.Rational(41,2500)

def basis_n(d, ncoeff):
    """The code's basis truncated to the fitted rank: 9 = quadratic, 3 = linear."""
    full = M.basis(d,3)
    return sp.Matrix(full[:3]) if ncoeff == 3 else full

def run(label, offsets, c, faces=()):
    slots = [sp.Matrix(o) for o in offsets] + [sp.Matrix(f) for f in faces]
    for ncoeff in (9, 3):
        Mn = sp.zeros(ncoeff, ncoeff)
        for s in slots:
            b = basis_n(s, ncoeff); Mn += (b*b.T)/sum(x**2 for x in s)
        if Mn.det() != 0:
            break
    else:
        print(f"{label:34s} SINGULAR at every rank"); return
    z = Mn.LUsolve(basis_n(c, ncoeff))
    g = [sp.simplify((z.T*basis_n(s,ncoeff))[0]/sum(x**2 for x in s)) for s in slots]
    ncell = len(offsets)
    gcell, gface = g[:ncell], g[ncell:]
    centre = 1 - sum(g) + sum(gface)      # zeroGradient faces carry psi_c
    lam = float(sp.N(sp.Abs(centre) + sum(sp.Abs(x) for x in gcell)))
    gv = np.array([float(sp.N(x)) for x in gcell]); sv = np.array(offsets, float)
    cst = float(sp.N(centre))
    n=64; ax=[np.linspace(-np.pi,np.pi,n,endpoint=False)]*3
    G=np.meshgrid(*ax,indexing="ij"); ph=sum(G[a][...,None]*sv[None,...,a] for a in range(3))
    mg=np.abs(cst+(gv*np.exp(1j*ph)).sum(axis=-1))
    fit = "quadratic" if ncoeff==9 else "LINEAR"
    print(f"{label:34s} {len(offsets):2d}+{len(faces)} slots  {fit:9s} "
          f"Lambda={lam:7.4f}  max|G|={mg.max():.6f}")

cube  = [s for s in itertools.product((-1,0,1),repeat=3) if any(s)]
face1 = [s for s in cube if s[0] >= 0]
edge  = [s for s in cube if s[0] >= 0 and s[1] >= 0]
corner= [s for s in cube if s[0] >= 0 and s[1] >= 0 and s[2] >= 0]
H = sp.Rational(1,2)

for axis in ("x","xyz"):
    c = M.displacement(Co,3,axis)
    print(f"=== displacement along {axis}, Co={float(Co)}  (MEASURED: Lambda_max=1.0527, rho=1.00441) ===")
    run("interior (27-pt cube)", cube, c)
    run("face cell, stencilBoundaryFaces=exclude", face1, c)
    run("face cell + zeroGradient face", face1, c, faces=[(-H,0,0)])
    run("edge cell + 2 faces", edge, c, faces=[(-H,0,0),(0,-H,0)])
    run("corner cell + 3 faces", corner, c, faces=[(-H,0,0),(0,-H,0),(0,0,-H)])
    print()
