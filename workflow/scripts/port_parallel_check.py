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
    """{case index: ([(t, L1), ...], column)} for every case of the study (both rows)"""
    out = {}
    for c in sorted(glob.glob(os.path.join(study, "sphereRotation3D_0*"))):
        if "." in os.path.basename(c) or not os.path.isfile(os.path.join(c, "advectorVoF.csv")): continue
        rows = list(csv.DictReader(open(os.path.join(c, "advectorVoF.csv"))))
        if not rows: continue
        col = next((k for k in rows[0] if k in ("L1_alphaf", "L1af", "alphafL1", "L1alphaf")), None)
        if col is None:
            col = next(k for k in rows[0] if "L1" in k and "alphaf" in k.lower())
        out[os.path.basename(c)[-5:]] = ([(i, float(r[col])) for i, r in enumerate(rows) if r[col] not in ("", None)], col)
    return out or None


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
    print("(b) np k against np 1: max over the steps of |L1(alphaf)(np k) - L1(alphaf)(np 1)|, per case (row)")
    for route in ("leia", "fork"):
        base = l1_series(os.path.join(sd, f"port_parN32np1_{route}"))
        if base is None:
            print(f"  {route}: np 1 missing"); continue
        for k in NPS[1:]:
            s = l1_series(os.path.join(sd, f"port_parN32np{k}_{route}"))
            if s is None:
                print(f"  {route} np {k}: missing"); continue
            for idx in sorted(base):
                if idx not in s: print(f"  {route} np {k} case {idx}: missing"); continue
                b, col = base[idx]; q = s[idx][0]
                n = min(len(b), len(q))
                d = max(abs(b[i][1] - q[i][1]) for i in range(n))
                print(f"  {route} np {k} vs np 1, case {idx}: max |dL1({col})| = {d:.3e} over {n} steps, end value {b[n-1][1]:.6e} -> {'within 1e-12' if d <= 1e-12 else 'ABOVE 1e-12'}")


if __name__ == "__main__":
    main()
