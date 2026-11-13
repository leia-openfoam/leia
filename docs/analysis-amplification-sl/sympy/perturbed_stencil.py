#!/usr/bin/env python3
"""Does a PERTURBED hexahedral stencil amplify, and can the weights repair it?

WHY. rho = 1 in the interior is a property of the PERFECT lattice. A mesher does not make
one: an unstructured hexahedral mesh has cell centres displaced and faces tilted. This asks
what that costs, and whether the weights -- the one thing we are free to choose -- can buy it
back. Section "Can the weights be optimised instead" of ms.tex.
"""
import itertools, sys
import numpy as np
from scipy.optimize import minimize

def basis(p):
    return np.stack([p[:,0],p[:,1],p[:,2],.5*p[:,0]**2,.5*p[:,1]**2,.5*p[:,2]**2,
                     p[:,0]*p[:,1],p[:,0]*p[:,2],p[:,1]*p[:,2]],axis=1)

CUBE = np.array([s for s in itertools.product((-1,0,1),repeat=3) if any(s)], float)

def max_gain(phys, ints, c, w2, n=40):
    B = basis(phys); M = (B*w2[:,None]).T@B
    g = w2*(B@np.linalg.solve(M, basis(np.array([c]))[0]))
    const = 1-g.sum(); sv = np.array(ints,float)
    ax=[np.linspace(-np.pi,np.pi,n,endpoint=False)]*3
    G=np.meshgrid(*ax,indexing="ij")
    ph=sum(G[a][...,None]*sv[None,...,a] for a in range(3))
    return float(np.abs(const+(g*np.exp(1j*ph)).sum(axis=-1)).max())

def main():
    rng = np.random.default_rng(7)
    idw = lambda P: np.sum(P**2,axis=1)**-1.0
    print("1. A perturbed stencil amplifies; a perfect one does not.\n")
    print(f"{'displacement':>16} {'Co=0.0164':>12} {'Co=0.3':>10}")
    for amp in (0.0,0.01,0.02,0.05,0.10,0.20):
        a=b=0.0
        for _ in range(12 if amp else 1):
            P = CUBE + (rng.random(CUBE.shape)-0.5)*2*amp if amp else CUBE.copy()
            a=max(a,max_gain(P,CUBE,[-0.0164,0,0],idw(P)))
            b=max(b,max_gain(P,CUBE,[-0.3,0,0],idw(P)))
        print(f"{amp*100:13.0f} %  {a:12.6f} {b:10.6f}")

    rng = np.random.default_rng(7)
    P = CUBE + (rng.random(CUBE.shape)-0.5)*2*0.10
    w0 = idw(P); x0 = np.log(w0)
    r = minimize(lambda lw: max_gain(P,CUBE,[-0.3,0,0],np.exp(lw),n=24), x0,
                 method="Nelder-Mead", options={"maxiter":4000,"fatol":1e-10,"xatol":1e-8})
    w1 = np.exp(r.x)
    print("\n2. Optimising all 26 weights at ONE displacement (10 % perturbed, Co = 0.3):")
    print(f"   inverse distance  max|G| = {max_gain(P,CUBE,[-0.3,0,0],w0):.6f}")
    print(f"   optimised         max|G| = {max_gain(P,CUBE,[-0.3,0,0],w1):.6f}  "
          f"({r.nit} iterations)")
    print("\n3. ...and it does not generalise:\n")
    print(f"{'displacement':>26} {'inverse distance':>18} {'optimised':>12}")
    for lab,c in (("Co 0.3  along x",[-0.3,0,0]), ("Co 0.3  along y",[0,-0.3,0]),
                  ("Co 0.3  diagonal",[-0.1732]*3), ("Co 0.5  along x",[-0.5,0,0]),
                  ("Co 0.1  along x",[-0.1,0,0])):
        print(f"{lab:>26} {max_gain(P,CUBE,c,w0):18.6f} {max_gain(P,CUBE,c,w1):12.6f}")
    print(f"\n   weight spread: inverse distance {w0.min():.2e}..{w0.max():.2e}, "
          f"optimised {w1.min():.2e}..{w1.max():.2e}")
    # quadratic reproduction is weight-independent; confirm it survives
    q = lambda p: .37+.7*p[:,0]-.4*p[:,1]+.65*p[:,2]+.5*(1.3*p[:,0]**2-.9*p[:,1]**2)+.5*p[:,0]*p[:,1]
    for nm,w2 in (("inverse distance",w0),("optimised",w1)):
        B=basis(P); M=(B*w2[:,None]).T@B
        coef=np.linalg.solve(M,(B*w2[:,None]).T@(q(P)-q(np.zeros((1,3)))))
        c=np.array([[-0.3,0,0]])
        print(f"   quadratic reproduced, {nm:17s}: error "
              f"{abs(q(np.zeros((1,3)))[0]+basis(c)[0]@coef-q(c)[0]):.2e}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
