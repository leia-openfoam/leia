#!/usr/bin/env python3
"""Split the translating-droplet disturbance U - U0 into a coherent DRIFT and a FLUCTUATION.

A zero-mean multipolar spurious current sitting on the interface cannot move the droplet;
a BIASED one propels it. This separates the two, per phase, and profiles the field against
each boundary so an apparent "reflection" in a log-scaled image can be checked as a number.

Measured on config/bestConfigAnim2D (rhoLENT + upwind, ratio 838.8, N=128): the droplet
drift starts at exactly 0 with fluctuation 0.0030*U0, and ends at drift 0.1289*U0 against
fluctuation 0.0151*U0 -- the ordering inverts, and the drift becomes 8.5x the fluctuation.
Droplet drift / ambient drift = 28 against the momentum-conserving prediction
rho_d V_d / (rho_a V_a) = 27: the capillary force is internal, so it cannot change total
momentum, but it drives the phases apart.

The boundary profile at t = 40.5 ms rises monotonically AWAY from every boundary toward
the droplet (inlet strip 0.00033 -> 0.00685 at 32 cells), so there is no boundary-localised
structure. There cannot be: the pressure equation is elliptic, so nothing propagates and
nothing reflects; apparent arcs are a smooth decay stretched over four decades of log colour.

Usage:  python3 workflow/scripts/translating_drift_decomposition.py [--case DIR] [--u0 0.05]
"""
import argparse
import numpy as np, re, os, sys
sys.path.insert(0,'workflow/scripts')
from make_parasitic_animation import read_field
_ap=argparse.ArgumentParser()
_ap.add_argument('--case', default='studies/bestConfigAnim2D/translatingDroplet2D_00000')
_ap.add_argument('--u0', type=float, default=0.05)
_a=_ap.parse_args()
D=_a.case
ts=sorted([d for d in os.listdir(D) if re.fullmatch(r'0\.\d+',d)],key=float)
n=128; L=0.01; h=L/n; U0=np.array([_a.u0,0,0])
xs=(np.arange(n)+0.5)*h; X,Y=np.meshgrid(xs,xs)
print(f"{'t[ms]':>7}{'max|dU|/U0':>12}{'drop drift/U0':>15}{'drop fluct/U0':>15}{'amb drift/U0':>14}{'amb fluct/U0':>14}")
sel=[0,8,24,40,60,80,100,119]
prof={}
for k in sel:
    t=ts[k]
    U=read_field(f'{D}/{t}/U',n*n).reshape(n,n,3)-U0
    al=read_field(f'{D}/{t}/alpha.water',n*n).reshape(n,n)
    m=np.hypot(U[...,0],U[...,1])
    din=al>0.5; amb=al<0.01
    dd=U[din].mean(axis=0)[:2]; df=(U[din][:,:2]-U[din].mean(axis=0)[:2])
    ad=U[amb].mean(axis=0)[:2]; af=(U[amb][:,:2]-U[amb].mean(axis=0)[:2])
    print(f"{float(t)*1e3:>7.2f}{m.max()/_a.u0:>12.4f}"
          f"{np.hypot(*dd)/_a.u0:>15.4f}{np.hypot(df[:,0],df[:,1]).mean()/_a.u0:>15.4f}"
          f"{np.hypot(*ad)/_a.u0:>14.4f}{np.hypot(af[:,0],af[:,1]).mean()/_a.u0:>14.4f}")
    prof[t]=(m,al)
# boundary structure: mean |dU| in columns/rows at given distance from each boundary
t=ts[80]; m,al=prof[t]
print(f"\nBOUNDARY PROFILE at t={float(t)*1e3:.1f} ms  (mean |U-U0|/U0 over the strip)")
print(f"{'cells from bdry':>16}{'inlet(x=0)':>12}{'outlet(x=L)':>13}{'wall y=0':>11}{'wall y=L':>11}")
for d in [0,1,2,4,8,16,32]:
    print(f"{d:>16}{m[:,d].mean()/_a.u0:>12.5f}{m[:,n-1-d].mean()/_a.u0:>13.5f}"
          f"{m[d,:].mean()/_a.u0:>11.5f}{m[n-1-d,:].mean()/_a.u0:>11.5f}")
ay=np.where(al.sum(axis=1)>0)[0]; ax=np.where(al.sum(axis=0)>0)[0]
print(f"droplet occupies columns {ax.min()}..{ax.max()}, rows {ay.min()}..{ay.max()}")
