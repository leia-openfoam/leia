#!/usr/bin/env python3
"""HYPOTHESIS: a boundary cell amplifies when the departure point lies OUTSIDE the data,
i.e. the fit EXTRAPOLATES. If so, the sign of the flow decides, and the boundary face slot
helps only because it sits on the far side of the foot and turns extrapolation into
interpolation -- not because it carries information (it carries none: it is a copy of psi_c)."""
import sys, itertools, numpy as np
sys.path.insert(0,"workflow/scripts")

def basis_rows(p):
    return np.stack([p[:,0],p[:,1],p[:,2],.5*p[:,0]**2,.5*p[:,1]**2,.5*p[:,2]**2,
                     p[:,0]*p[:,1],p[:,0]*p[:,2],p[:,1]*p[:,2]],axis=1)

def mg(cells, faces, c, wface=1.0):
    slots=np.array(list(cells)+list(faces),float); ncell=len(cells)
    w2=np.sum(slots**2,axis=1)**-1.0
    w2[ncell:] *= wface                      # optional cap on the face weight
    for nc in (9,3):
        B=basis_rows(slots)[:,:nc]; Mn=(B*w2[:,None]).T@B
        if np.linalg.cond(Mn)<1e12: break
    z=np.linalg.solve(Mn, basis_rows(np.array(c,float).reshape(1,-1))[0][:nc])
    g=w2*(B@z); gc,gf=g[:ncell],g[ncell:]; centre=1-g.sum()+gf.sum()
    sv=np.array(list(cells),float); n=56
    ax=[np.linspace(-np.pi,np.pi,n,endpoint=False)]*3
    G=np.meshgrid(*ax,indexing="ij"); ph=sum(G[a][...,None]*sv[None,...,a] for a in range(3))
    return np.abs(centre+(gc*np.exp(1j*ph)).sum(axis=-1)).max()

cube=[s for s in itertools.product((-1,0,1),repeat=3) if any(s)]
face1=[s for s in cube if s[0]>=0]        # wall at -x
H=0.5
print("Wall at -x. The foot is at x = c_x in cell widths.")
print("c_x < 0 : flow INTO the wall, foot on the WALL side  -> outside the cell data")
print("c_x > 0 : flow AWAY from the wall, foot INSIDE the cell data\n")
print(f"{'Co':>7} {'c_x':>7} | {'no face (A)':>12} {'with face (B)':>14} | interior")
for Co in (0.0164,0.1,0.25,0.5):
    for sgn,lab in ((-1,"toward wall"),(+1,"into domain")):
        c=(sgn*Co,0,0)
        a=mg(face1,[],c); b=mg(face1,[(-H,0,0)],c); i=mg(cube,[],c)
        print(f"{Co:7.4f} {sgn*Co:+7.4f} | {a:12.6f} {b:14.6f} | {i:.6f}   ({lab})")
print("\nFace-weight sweep at Co = 0.25, foot toward the wall (the extrapolating case):")
for wf in (1.0,0.5,0.25,0.1,0.01):
    print(f"   face weight x {wf:<5} -> max|G| = {mg(face1,[(-H,0,0)],(-0.25,0,0),wf):.6f}")
