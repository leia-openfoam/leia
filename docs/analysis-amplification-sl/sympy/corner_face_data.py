#!/usr/bin/env python3
"""Reconstruct from a cell's own corner and face values instead of neighbour centres.

Corner and face values are not stored; they are interpolated from cell centres, so every
datum is a LINEAR COMBINATION of the 27 cell values of the 3x3x3 block. The whole scheme is
then still one weight per cell, and its growth factor follows exactly as before.
"""
import itertools, sys
import numpy as np

OFF = [s for s in itertools.product((-1,0,1), repeat=3)]
IDX = {s:i for i,s in enumerate(OFF)}
E = lambda s: np.eye(27)[IDX[tuple(s)]]
C0 = E((0,0,0))

def basis(d):
    d = np.atleast_2d(np.asarray(d, float))
    return np.stack([d[:,0],d[:,1],d[:,2], .5*d[:,0]**2,.5*d[:,1]**2,.5*d[:,2]**2,
                     d[:,0]*d[:,1], d[:,0]*d[:,2], d[:,1]*d[:,2]], axis=1)

def centres(cells):            # (position, value-weights) of neighbour centres
    return [(np.array(s,float), E(s)) for s in cells if any(s)]

def faces(avail):
    out=[]
    for a in range(3):
        for sg in (-1,1):
            d=np.zeros(3); d[a]=sg/2; nb=[0,0,0]; nb[a]=sg
            v = 0.5*(C0+E(nb)) if tuple(nb) in avail else C0      # zeroGradient: own value
            out.append((d,v))
    return out

def corners(avail):
    out=[]
    for sg in itertools.product((-1,1),repeat=3):
        cells=[c for c in itertools.product(*[(0,s) for s in sg]) if c in avail]
        out.append((np.array(sg,float)/2, sum(E(c) for c in cells)/len(cells)))
    return out

def scheme(data, c):
    """Weights over the 27 cells of the anchored IDW quadratic fit through `data`, read at c."""
    P = np.array([d for d,_ in data]); V = np.array([v for _,v in data])
    w2 = np.sum(P**2,axis=1)**-1.0
    B = basis(P); M = (B*w2[:,None]).T@B
    rhs = (B*w2[:,None]).T@(V - C0)              # 9 x 27
    return C0 + basis(c)[0]@np.linalg.solve(M, rhs)

def gain(w, n=40):
    S = np.array(OFF,float)
    ax=[np.linspace(-np.pi,np.pi,n,endpoint=False)]*3
    G=np.meshgrid(*ax,indexing="ij"); ph=sum(G[a][...,None]*S[None,...,a] for a in range(3))
    m = np.abs((w*np.exp(1j*ph)).sum(axis=-1))
    return m.max(), m.min()

def quad_err(w, c):
    q = lambda x: .37+.7*x[:,0]-.4*x[:,1]+.65*x[:,2]+.5*(1.3*x[:,0]**2-.9*x[:,1]**2+.4*x[:,2]**2)\
                  +.5*x[:,0]*x[:,1]-.3*x[:,1]*x[:,2]
    S = np.array(OFF,float)
    return abs(w@q(S) - q(np.array([c]))[0])

def threshold(build, direction):
    d=np.array(direction,float); d/=np.linalg.norm(d); lo,hi=0.0,1.5
    for _ in range(36):
        mid=(lo+hi)/2
        if gain(build(-mid*d), n=28)[0] > 1+1e-9: hi=mid
        else: lo=mid
    return lo

ALL = set(OFF)
base = lambda c: scheme(centres(OFF), c)
armA = lambda c: scheme(faces(ALL)+corners(ALL), c)
armB = lambda c: scheme(centres(OFF)+faces(ALL)+corners(ALL), c)
def armC(c):   # control: corner/face values SAMPLED from the baseline fit at those points
    data=[(d, scheme(centres(OFF), d)) for d,_ in faces(ALL)+corners(ALL)]
    return scheme(data, c)

def main():
    print("INTERIOR CELL, perfect cube lattice\n")
    arms=[("baseline: 26 neighbour centres", base),
          ("A: own 6 faces + 8 corners (averaged)", armA),
          ("B: 26 centres + own faces + corners", armB),
          ("C: faces+corners SAMPLED from baseline", armC)]
    print(f"{'arm':42s} {'Co':>6} {'max|G|':>10} {'min|G|':>10} {'quad err':>10}")
    for lab,f in arms:
        for Co in (0.0164,0.25,0.5):
            c=[-Co,0,0]; w=f(c); mx,mn=gain(w)
            print(f"{lab:42s} {Co:6.4f} {mx:10.6f} {mn:10.6f} {quad_err(w,c):10.2e}")
        lab=""
    print("\ncontrol C equals baseline to round-off:",
          max(np.abs(armC([-Co,0,0])-base([-Co,0,0])).max() for Co in (0.0164,0.25,0.5)))
    print("\nCourant limit (first Co with max|G| > 1):")
    for lab,f in arms[:3]:
        print(f"  {lab:42s} axis {threshold(f,(1,0,0)):.4f}   diagonal {threshold(f,(1,1,1)):.4f}")

    print("\nBOUNDARY CELL, wall on the -x side, departure point towards the wall")
    av = {s for s in OFF if s[0] >= 0}
    bnd_base = lambda c: scheme(centres([s for s in OFF if s[0]>=0])
                                + [(np.array([-.5,0,0]), C0)], c)   # production: 17 + wall face
    bnd_A = lambda c: scheme(faces(av)+corners(av), c)
    bnd_B = lambda c: scheme(centres([s for s in OFF if s[0]>=0])+faces(av)+corners(av), c)
    S=np.array(OFF,float); ok=S[:,0]>=0
    print(f"{'arm':42s} {'Co':>6} {'max|G|':>10} {'sum|g| nbrs':>12}")
    for lab,f in (("production: 17 centres + wall face", bnd_base),
                  ("A: own faces + corners (wall: one-sided)", bnd_A),
                  ("B: 17 centres + own faces + corners", bnd_B)):
        for Co in (0.0164,0.25,0.5):
            w=f([-Co,0,0]); mx,_=gain(w)
            nb=np.abs(np.delete(w, IDX[(0,0,0)])).sum()
            print(f"{lab:42s} {Co:6.4f} {mx:10.6f} {nb:12.6f}")
        lab=""
    return 0

if __name__=="__main__":
    sys.exit(main())
