#!/usr/bin/env python3
"""On an IRREGULAR (perturbed) hex stencil: does raising the polynomial degree, or moving
the weights to the foot, change the amplification? Worst case over flow directions."""
import itertools, sys
import numpy as np

CUBE = np.array([s for s in itertools.product((-1,0,1),repeat=3) if any(s)], float)

def monomials(deg):
    return [e for e in itertools.product(range(deg+1),repeat=3) if 0 < sum(e) <= deg]

def basis(d, deg):
    d=np.atleast_2d(d)
    return np.stack([np.prod(d**np.array(e),axis=1) for e in monomials(deg)],axis=1)

def weights(P, c, deg, anchor):
    ref = np.zeros(3) if anchor=="centre" else np.array(c)
    w2 = np.sum((P-ref)**2,axis=1)**-1.0
    B = basis(P,deg); M=(B*w2[:,None]).T@B
    return w2*(B@np.linalg.solve(M, basis(np.array(c),deg)[0]))

def max_gain(g, n=32):
    const=1-g.sum()
    ax=[np.linspace(-np.pi,np.pi,n,endpoint=False)]*3
    G=np.meshgrid(*ax,indexing="ij"); ph=sum(G[a][...,None]*CUBE[None,...,a] for a in range(3))
    return float(np.abs(const+(g*np.exp(1j*ph)).sum(axis=-1)).max())

def exact_err(P, c, deg, anchor, fdeg):
    rng=np.random.default_rng(3); coef={e:rng.normal() for e in monomials(fdeg)}
    f=lambda X: sum(v*np.prod(np.atleast_2d(X)**np.array(e),axis=1) for e,v in coef.items())
    g=weights(P,c,deg,anchor)
    return abs(g@(f(P)-0.0) - f(np.array(c))[0])        # f(0)=0, anchored

def main():
    rng=np.random.default_rng(7)
    P = CUBE + (rng.random(CUBE.shape)-0.5)*2*0.10       # 10 % perturbed, as in pass 3
    dirs = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1),(1,-1,0),(-1,1,1)]
    arms = [("linear,    centre weights",1,"centre"),
            ("quadratic, centre weights (production)",2,"centre"),
            ("quadratic, foot weights",2,"foot"),
            ("cubic,     centre weights",3,"centre"),
            ("cubic,     foot weights",3,"foot")]
    print("10 % perturbed hex stencil, worst max|G| over 9 flow directions\n")
    print(f"{'reconstruction':40s} {'coeffs':>6} {'Co 0.1':>9} {'Co 0.3':>9} {'Co 0.5':>9}  exact for")
    for lab,deg,anc in arms:
        row=[]
        for Co in (0.1,0.3,0.5):
            worst=0.0
            for d in dirs:
                v=np.array(d,float); v/=np.linalg.norm(v)
                worst=max(worst, max_gain(weights(P,(-Co*v).tolist(),deg,anc)))
            row.append(worst)
        ex = [k for k in (1,2,3) if exact_err(P,[-0.3,0,0],deg,anc,k) < 1e-10]
        print(f"{lab:40s} {len(monomials(deg)):6d} {row[0]:9.6f} {row[1]:9.6f} {row[2]:9.6f}  "
              f"degree <= {max(ex) if ex else 0}")
    print("\nperfect lattice, same arms, Co 0.5 worst over directions:")
    for lab,deg,anc in arms:
        try:
            worst=max(max_gain(weights(CUBE,(-0.5*np.array(d,float)/np.linalg.norm(d)).tolist(),deg,anc))
                      for d in dirs)
            print(f"   {lab:40s} {worst:9.6f}")
        except np.linalg.LinAlgError:
            B=basis(CUBE,deg); print(f"   {lab:40s}  SINGULAR: rank {np.linalg.matrix_rank(B)} of "
                                     f"{B.shape[1]} (on {{-1,0,1}} x^3 = x)")
    # condition numbers explain the perturbed cubic
    rng=np.random.default_rng(7); P=CUBE+(rng.random(CUBE.shape)-0.5)*2*0.10
    for deg in (1,2,3):
        B=basis(P,deg); w2=np.sum(P**2,axis=1)**-1.0
        print(f"   condition of the normal matrix, degree {deg}, perturbed: "
              f"{np.linalg.cond((B*w2[:,None]).T@B):.2e}")
    return 0

if __name__=="__main__":
    sys.exit(main())
