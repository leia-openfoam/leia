#!/usr/bin/env python3
"""The parallel identity of the port (docs/davof/PORT-DESIGN.md 4.2, 2026-11-06): the
studies port_parN32np<k>_{fork,leia}, k in 1, 2, 4, 8, 16 (sphereRotation3D, rotation,
N 32, the reference setting + the four round-4 keys, END_TIME 0.05, writePrecision 17).

  (a) leia at np k against the fork at np k: port_identity_check.py on the two studies
      (the reconstructed fields of the last time, the function-object files, the CSV
      numbers, the logs);
  (b) leia at np k against leia at np 1: the per-step L1(alphaf) of advectorVoF.csv,
      max over the steps of |L1(np k) - L1(np 1)|, the criterion 1e-12 (the fork's own
      parallel scatter, REPORT-swept.md 2026-11-02); the same for the fork's np k against
      its np 1, for the record.

Usage: port_parallel_check.py <studies_dir>   (e.g. studies/davof)
"""
import csv, glob, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
NPS = (1, 2, 4, 8, 16)

def l1_series(study):
    cases = sorted(c for c in glob.glob(os.path.join(study, "sphereRotation3D_0*")) if "." not in os.path.basename(c))
    if len(cases) != 1:
        return None
    rows = list(csv.DictReader(open(os.path.join(cases[0], "advectorVoF.csv"))))
    col = next((c for c in rows[0] if c in ("L1_alphaf", "L1af", "alphafL1", "L1alphaf")), None)
    if col is None:
        col = next(c for c in rows[0] if "L1" in c and "alphaf" in c.lower())
    return [(float(r.get("time", r.get("t", i))), float(r[col])) for i, r in enumerate(rows) if r[col] not in ("", None)], col

def main():
    sd = sys.argv[1]
    print("(a) leia np k against fork np k (port_identity_check.py)")
    verdict = {}
    for k in NPS:
        f, l = os.path.join(sd, f"port_parN32np{k}_fork"), os.path.join(sd, f"port_parN32np{k}_leia")
        if not (os.path.isdir(f) and os.path.isdir(l)):
            print(f"  np {k}: missing"); continue
        r = subprocess.run([sys.executable, os.path.join(HERE, "port_identity_check.py"), f, l], capture_output=True, text=True)
        last = r.stdout.strip().split("\n")[-1]
        lines = [x for x in r.stdout.split("\n") if x.startswith("  ")]
        verdict[k] = last
        print(f"  np {k}: {last}")
        for x in lines: print("     " + x[:200])
    print("(b) np k against np 1: max over the steps of |L1(alphaf)(np k) - L1(alphaf)(np 1)|")
    for route in ("leia", "fork"):
        base = l1_series(os.path.join(sd, f"port_parN32np1_{route}"))
        if base is None:
            print(f"  {route}: np 1 missing"); continue
        b, col = base
        for k in NPS[1:]:
            s = l1_series(os.path.join(sd, f"port_parN32np{k}_{route}"))
            if s is None:
                print(f"  {route} np {k}: missing"); continue
            s = s[0]
            n = min(len(b), len(s))
            d = max(abs(b[i][1] - s[i][1]) for i in range(n))
            rel = d / max(abs(b[n-1][1]), 1e-300)
            print(f"  {route} np {k} vs np 1: max |dL1({col})| = {d:.3e} over {n} steps (relative to the end value {rel:.3e}) -> {'within 1e-12' if d <= 1e-12 else 'ABOVE 1e-12'}")

if __name__ == "__main__":
    main()
