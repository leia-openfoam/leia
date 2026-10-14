#!/usr/bin/env python3
"""Four ways to make a boundary cell's quadratic fit well posed, compared on one footing."""
import sys, itertools, numpy as np
sys.path.insert(0,"workflow/scripts"); import sl_amplification_symbol as M

def basis_rows(p):
    return np.stack([p[:,0],p[:,1],p[:,2],
                     .5*p[:,0]**2,.5*p[:,1]**2,.5*p[:,2]**2,
                     p[:,0]*p[:,1],p[:,0]*p[:,2],p[:,1]*p[:,2]],axis=1)

def arm(label, cells, faces, Co):
    c=np.array([-Co,0,0]); slots=np.array(list(cells)+list(faces),float); ncell=len(cells)
    w2=np.sum((slots)**2,axis=1)**-1.0
    for nc in (9,3):
        B=basis_rows(slots)[:,:nc]; Mn=(B*w2[:,None]).T@B
        if np.linalg.cond(Mn)<1e12: break
    bc=basis_rows(c.reshape(1,-1))[0][:nc]
    z=np.linalg.solve(Mn,bc); g=w2*(B@z)
    gc,gf=g[:ncell],g[ncell:]
    centre=1-g.sum()+gf.sum()        # zeroGradient faces carry psi_c
    lam=abs(centre)+np.abs(gc).sum()
    sv=np.array(list(cells),float); n=56
    ax=[np.linspace(-np.pi,np.pi,n,endpoint=False)]*3
    G=np.meshgrid(*ax,indexing="ij"); ph=sum(G[a][...,None]*sv[None,...,a] for a in range(3))
    mg=np.abs(centre+(gc*np.exp(1j*ph)).sum(axis=-1))
    print(f"  {label:52s} {len(cells):2d}c+{len(faces)}f  {'quad' if nc==9 else 'LIN '}"
          f"  Lambda={lam:7.4f}  max|G|={mg.max():.6f}  cond={np.linalg.cond(Mn):.1e}")

cube=[s for s in itertools.product((-1,0,1),repeat=3) if any(s)]
face1=[s for s in cube if s[0]>=0]                      # wall at -x: 17 interior neighbours
layer2=[(2,j,k) for j in (-1,0,1) for k in (-1,0,1)]    # the SECOND interior layer, 9 cells
H=0.5

for Co in (0.0164, 0.25):
    print(f"=== boundary (face) cell, wall at -x, Co = {Co} ===")
    arm("interior cell, for reference", cube, [], Co)
    arm("A exclude faces (rank-deficient -> linear)", face1, [], Co)
    arm("B include the zeroGradient face (PRODUCTION)", face1, [(-H,0,0)], Co)
    arm("C one-sided: add the 2nd interior layer, no face", face1+layer2, [], Co)
    arm("D both: 2nd layer AND the face", face1+layer2, [(-H,0,0)], Co)
    print()
