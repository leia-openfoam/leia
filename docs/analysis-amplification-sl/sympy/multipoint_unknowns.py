#!/usr/bin/env python3
"""Semi-Lagrangian transport of centre, face-centre and vertex values as INDEPENDENT unknowns.

Each unknown is moved along its own trajectory. Its foot is located in the host cell that
contains it, and the host's own values are interpolated there. Several unknowns per cell
make the update of one Fourier wave a small matrix; the growth factor is the largest
eigenvalue magnitude over all waves.
"""
import itertools, sys
import numpy as np

INTERP15=False
def setup(d, types):
    """types: subset of {0,1/2}^d positions kept as unknowns. Returns host-fit weights fn."""
    Q = [q for q in itertools.product((-0.5,0.0,0.5), repeat=d)
         if tuple(abs(x) for x in q) in types]                  # host-local points of kept types
    full_tensor = len(Q) == 3**d
    if full_tensor:
        def L1(x, node):
            others=[o for o in (-0.5,0.0,0.5) if o!=node]
            return np.prod([(x-o)/(node-o) for o in others])
        def wts(xi): return np.array([np.prod([L1(xi[i],q[i]) for i in range(d)]) for q in Q])
    else:
        mons=[e for e in itertools.product(range(3),repeat=d) if sum(e)<=2]   # complete quadratic
        if INTERP15 and d==3 and len(Q)==15:
            mons=[e for e in itertools.product(range(2),repeat=3)]+[(2,0,0),(0,2,0),(0,0,2),
                  (1,2,2),(2,1,2),(2,2,1),(2,2,2)]
        B=np.array([[np.prod(np.array(q)**np.array(e)) for e in mons] for q in Q])
        Bp=np.linalg.pinv(B)
        def wts(xi): return np.array([np.prod(np.array(xi)**np.array(e)) for e in mons])@Bp
    return Q, wts, full_tensor

def matrix_symbol(d, types, c):
    T=sorted(types); tix={t:i for i,t in enumerate(T)}
    Q,wts,_=setup(d,set(T))
    rows=[]                                   # per unknown: list of (target type, cell offset, weight)
    for t in T:
        f=np.array(t)-c; k=np.floor(f+0.5); xi=f-k
        w=wts(xi); ent=[]
        for q,wq in zip(Q,w):
            p=tuple(abs(x) for x in q); off=k+np.array([-1 if x==-0.5 else 0 for x in q])
            ent.append((tix[p],off,wq))
        rows.append(ent)
    return T,rows

def rho(d, types, c, n):
    T,rows=matrix_symbol(d,types,np.array(c,float)); nt=len(T)
    ax=np.linspace(-np.pi,np.pi,n,endpoint=False)
    TH=np.stack(np.meshgrid(*([ax]*d),indexing="ij"),axis=-1).reshape(-1,d)
    A=np.zeros((len(TH),nt,nt),complex)
    for i,ent in enumerate(rows):
        for j,off,w in ent:
            A[:,i,j]+=w*np.exp(1j*(TH@off))
    ev=np.linalg.eigvals(A)
    mags=np.abs(ev)
    # physical eigenvalue error along the flow, small theta
    return mags.max()

def phys_err(d,types,c,theta_dir,eps):
    T,rows=matrix_symbol(d,types,np.array(c,float)); nt=len(T)
    th=np.array(theta_dir,float)*eps
    A=np.zeros((nt,nt),complex)
    for i,ent in enumerate(rows):
        for j,off,w in ent: A[i,j]+=w*np.exp(1j*(th@off))
    ev=np.linalg.eigvals(A); ex=np.exp(-1j*(th@np.array(c)))
    return abs(ev[np.argmin(abs(ev-ex))]-ex)

def cc_scheme_err(d,c,theta_dir,eps):          # current cell-centre quadratic LSQ, scalar symbol
    S=np.array([s for s in itertools.product((-1,0,1),repeat=d) if any(s)],float)
    mons=[e for e in itertools.product(range(3),repeat=d) if 0<sum(e)<=2]
    b=lambda X: np.stack([np.prod(np.atleast_2d(X)**np.array(e),axis=1)*(0.5 if max(e)==2 and sum(e)==2 else 1) for e in mons],axis=1)
    w2=np.sum(S**2,axis=1)**-1.0; B=b(S); M=(B*w2[:,None]).T@B
    g=w2*(B@np.linalg.solve(M,b(-np.array(c))[0]))
    th=np.array(theta_dir,float)*eps
    G=(1-g.sum())+(g*np.exp(1j*(S@th))).sum()
    return abs(G-np.exp(-1j*(th@np.array(c))))

def cc_rho(d,c,n):
    S=np.array([s for s in itertools.product((-1,0,1),repeat=d) if any(s)],float)
    mons=[e for e in itertools.product(range(3),repeat=d) if 0<sum(e)<=2]
    b=lambda X: np.stack([np.prod(np.atleast_2d(X)**np.array(e),axis=1) for e in mons],axis=1)
    w2=np.sum(S**2,axis=1)**-1.0; B=b(S); M=(B*w2[:,None]).T@B
    g=w2*(B@np.linalg.solve(M,b(np.array(c))[0]))
    ax=np.linspace(-np.pi,np.pi,n,endpoint=False)
    TH=np.stack(np.meshgrid(*([ax]*d),indexing="ij"),axis=-1).reshape(-1,d)
    return np.abs((1-g.sum())+(np.exp(1j*TH@S.T)*g).sum(axis=1)).max()

def main():
    h=0.5
    global INTERP15
    cases=[("2D: centre + faces + vertices (3x3 tensor, biquadratic)",2,
            {(0,0),(h,0),(0,h),(h,h)},40),
           ("3D: centre + 6 faces + 8 vertices (15 pts, LSQ quadratic)",3,
            {(0,0,0),(h,0,0),(0,h,0),(0,0,h),(h,h,h)},20),
           ("3D: + 12 edge midpoints (3x3x3 tensor, triquadratic)",3,
            {t for t in itertools.product((0,h),repeat=3)},20)]
    dirs=lambda d: [v for v in ([(1,0),(1,1),(2,1)] if d==2 else [(1,0,0),(1,1,0),(1,1,1),(2,1,0)])]
    cases.insert(2,("3D: centre + faces + vertices, INTERPOLATORY 15-term space",3,
            {(0,0,0),(h,0,0),(0,h,0),(0,0,h),(h,h,h)},20))
    for lab,d,types,n in cases:
        INTERP15 = "INTERPOLATORY" in lab
        if INTERP15:
            Q,_,_=setup(3,types); 
            mons=[e for e in itertools.product(range(2),repeat=3)]+[(2,0,0),(0,2,0),(0,0,2),(1,2,2),(2,1,2),(2,2,1),(2,2,2)]
            V=np.array([[np.prod(np.array(q)**np.array(e)) for e in mons] for q in Q])
            print(f"\n   [15x15 nodal matrix: rank {np.linalg.matrix_rank(V)}, condition {np.linalg.cond(V):.1f}]",end="")
        print(f"\n{lab}   unknowns per cell: {len(types)}")
        print(f"   {'Co':>5} " + " ".join(f"{str(v):>11}" for v in dirs(d)) + f"   {'current':>9}")
        for Co in (0.0164,0.25,0.5,0.75,1.0,1.25):
            r=[]
            for v in dirs(d):
                u=np.array(v,float); u/=np.linalg.norm(u)
                r.append(rho(d,types,(Co*u).tolist(),n))
            u=np.array(dirs(d)[0],float)
            cur=max(cc_rho(d,(Co*np.array(v,float)/np.linalg.norm(v)).tolist(),n) for v in dirs(d))
            print(f"   {Co:5.3f} " + " ".join(f"{x:11.6f}" for x in r) + f"   {cur:9.6f}")
        # accuracy: error of the physical eigenvalue for a long wave along x, Co 0.3
        c=[0.3]+[0]*(d-1); td=[1]+[0]*(d-1)
        e1=[phys_err(d,types,c,td,eps) for eps in (0.2,0.1,0.05)]
        e0=[cc_scheme_err(d,c,td,eps) for eps in (0.2,0.1,0.05)]
        print(f"   one-step error, long wave, Co 0.3:  new {e1[0]:.2e} {e1[1]:.2e} {e1[2]:.2e}"
              f"  (x{e1[0]/e1[1]:.1f}, x{e1[1]/e1[2]:.1f} per halving)")
        print(f"                                   current {e0[0]:.2e} {e0[1]:.2e} {e0[2]:.2e}"
              f"  (x{e0[0]/e0[1]:.1f}, x{e0[1]/e0[2]:.1f})")
    return 0

if __name__=="__main__":
    sys.exit(main())
