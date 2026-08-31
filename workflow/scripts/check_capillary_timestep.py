#!/usr/bin/env python3
"""Verify that every hydrodynamic study's time step obeys the CAPILLARY limit.

THE CRITERION (Brackbill, Kothe & Zemach 1992, the standard explicit-CSF limit):
the step must resolve the fastest capillary wave the grid can carry, whose
wavelength is 2h:

    dt_cap = sqrt( (rho_1 + rho_2) * h^3 / (2*pi*sigma) )

WHAT THE REPOSITORY ACTUALLY DOES. materialize derives
MAX_DELTA_T = CAPILLARY_DT_COEFF / nRef^1.5 with nRef = CAPILLARY_REF_LENGTH*N/L,
which is a coefficient times h^1.5 -- the SAME scaling as the criterion, since
h = L/N. This script makes that identity explicit and reports the safety factor
dt/dt_cap for each rung, so the number is checked rather than trusted.

A study is FLAGGED if dt > dt_cap (unstable by the criterion) and reported with
its margin otherwise. Nothing here changes a step; it only measures.
"""
import argparse, math, sys

def dt_capillary(h, rho1, rho2, sigma):
    return math.sqrt((rho1 + rho2) * h**3 / (2.0 * math.pi * sigma))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rho1", type=float, default=998.2)
    ap.add_argument("--rho2", type=float, default=1.19)
    ap.add_argument("--sigma", type=float, default=0.07274)
    ap.add_argument("--coeff", type=float, default=0.010861,
                    help="CAPILLARY_DT_COEFF")
    ap.add_argument("--ref-length", type=float, default=0.01,
                    help="CAPILLARY_REF_LENGTH")
    a = ap.parse_args()

    # The coefficient expressed as dt = C_repo * h^1.5, independent of the box:
    #   nRef = Lref*N/L and h = L/N  =>  nRef = Lref/h
    #   dt   = coeff / nRef^1.5 = coeff * h^1.5 / Lref^1.5
    C_repo = a.coeff / a.ref_length**1.5
    C_crit = math.sqrt((a.rho1 + a.rho2) / (2.0 * math.pi * a.sigma))
    print("rho1=%g rho2=%g sigma=%g" % (a.rho1, a.rho2, a.sigma))
    print("repo law      dt = %.4f * h^1.5   (CAPILLARY_DT_COEFF=%g, ref L=%g)"
          % (C_repo, a.coeff, a.ref_length))
    print("Brackbill     dt = %.4f * h^1.5" % C_crit)
    print("SAFETY FACTOR dt/dt_cap = %.4f  %s\n"
          % (C_repo / C_crit,
             "STABLE (below the limit)" if C_repo < C_crit else "*** ABOVE THE LIMIT ***"))

    print("%-28s %-8s %12s %12s %12s %9s" %
          ("study / box", "N", "h", "dt (repo)", "dt_cap", "dt/dt_cap"))
    rows = [
        ("stationaryDroplet2D  L=0.01", 0.01, [32, 64, 128, 256]),
        ("stationary/translating 3D L=0.006", 0.006, [60, 76, 95]),
    ]
    bad = 0
    for label, L, Ns in rows:
        for N in Ns:
            h = L / N
            dt = C_repo * h**1.5
            dc = dt_capillary(h, a.rho1, a.rho2, a.sigma)
            flag = "" if dt < dc else "   *** UNSTABLE ***"
            if dt >= dc:
                bad += 1
            print("%-28s %-8d %12.6e %12.6e %12.6e %9.4f%s"
                  % (label.split()[0] + " " + label.split()[-1], N, h, dt, dc,
                     dt / dc, flag))
    print()
    if bad:
        print("%d rung(s) violate the capillary limit" % bad)
        return 1
    print("All rungs are below the capillary limit; the repo law is the "
          "criterion\nwith a constant safety factor, and the h^1.5 scaling is "
          "exact.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
