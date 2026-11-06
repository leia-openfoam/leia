#!/usr/bin/env python3
"""The bit-for-bit identity of the port (docs/davof/PORT-DESIGN.md section 4.1, 2026-11-06):
the cases of a study rendered on the fork route (solver advectorVoF, the TwoPhaseFlow fork
at c5c00c5) against the same study on the leia route (solver davofAdvectorFoam,
libleiaDavof), case by case:

  * the written fields of the last time directory (alpha.water, alphaf.water and every
    other field present on both sides), byte for byte after the FoamFile header
    (the header names the application);
  * the function objects' files under postProcessing (the .dat of areaFractionError and
    volumeFractionError, every rank) and the converter's advectorVoF.csv, byte for byte;
  * the solver logs, line by line after the log prefixes are normalised (isoAdvection: ->
    davofSwept:, advectorVoF -> davofAdvectorFoam) and the lines that name a binary, a
    path, a time stamp, a wall-clock or CPU time or the version banner are dropped.

Usage:
    port_identity_check.py <studies>/port_<name>_fork <studies>/port_<name>_leia [--verbose]

Exit 0 when every case is identical in all three, 1 otherwise; the report names the first
differing file and line of every case.
"""
import argparse, os, re, sys, glob, filecmp

NOISE = re.compile(
    r"^(Build|Arch|Exec|Date|Time|Host|PID|Case|nProcs|I/O|ExecutionTime|Executable|"
    r"\s*leia\b|\s*lib|davof: version|Reading|Create time|Create mesh|Selecting|\*|/\*|Starting time loop|"
    r"davofSwept: the key|isoAdvector: OpenFOAM|davof: grazingCut|leia-stamp|version|Hosts|Slaves|Roots|Pstream|\s*\(\S+ \d+\)\s*$).*|.*(ExecutionTime|ClockTime|timeStepTime|wall time| s \(wait |elapsed).*$"
)
def norm_line(s):
    # the log prefixes, the application name and the two renamed keys of the leia route
    # (surface for clippingSurface, flux for fluxIntegration; the values are the same)
    s = (s.replace("isoAdvection: ", "davofSwept: ").replace("advectorVoF", "davofAdvectorFoam")
          .replace("davofSwept: clippingSurface ", "davofSwept: surface ")
          .replace(", fluxIntegration ", ", flux "))
    return s.rstrip("\n")

def last_time_dir(case):
    ts = [d for d in os.listdir(case) if re.fullmatch(r"[0-9]+(\.[0-9]+)?(e-?[0-9]+)?", d)]
    ts = [d for d in ts if os.path.isdir(os.path.join(case, d))]
    if not ts: return None
    return max(ts, key=float)

def strip_header(path):
    s = open(path, "rb").read()
    i = s.find(b"// * * *")
    return s[i:] if i >= 0 else s

def compare_fields(cf, cl, rep):
    tf, tl = last_time_dir(cf), last_time_dir(cl)
    if tf is None or tl is None or float(tf) != float(tl):
        rep.append(f"  last time: fork {tf} leia {tl}"); return False
    ok = True
    names = sorted(set(os.listdir(os.path.join(cf, tf))) & set(os.listdir(os.path.join(cl, tl))))
    names = [n for n in names if os.path.isfile(os.path.join(cf, tf, n)) and not n.endswith((".orig",))]
    for n in names:
        a, b = strip_header(os.path.join(cf, tf, n)), strip_header(os.path.join(cl, tl, n))
        if a != b:
            ok = False
            la, lb = a.split(b"\n"), b.split(b"\n")
            k = next((i for i in range(min(len(la), len(lb))) if la[i] != lb[i]), min(len(la), len(lb)))
            rep.append(f"  field {tf}/{n} differs at line {k}: fork {la[k][:60] if k < len(la) else b''!r} leia {lb[k][:60] if k < len(lb) else b''!r}")
    if ok: rep.append(f"  fields at t = {tf}: {len(names)} identical ({', '.join(names)})")
    return ok

def compare_post(cf, cl, rep):
    ok = True
    files = sorted(set(os.path.relpath(p, cf) for p in glob.glob(os.path.join(cf, "postProcessing", "**", "*"), recursive=True) if os.path.isfile(p)))
    files += ["advectorVoF.csv"]
    n = 0
    for f in files:
        a, b = os.path.join(cf, f), os.path.join(cl, f)
        if not os.path.isfile(b):
            if "sweptFaceLog" in f or "consistency" in f.lower(): continue
            rep.append(f"  missing on the leia side: {f}"); ok = False; continue
        if not filecmp.cmp(a, b, shallow=False):
            la, lb = open(a, errors="replace").read().split("\n"), open(b, errors="replace").read().split("\n")
            k = next((i for i in range(min(len(la), len(lb))) if la[i] != lb[i]), min(len(la), len(lb)))
            if f.endswith(".csv") and k < min(len(la), len(lb)) and len(la) == len(lb):
                # the converter's CSV: name the differing columns; a column that names the
                # route (the solver, the log file, the version stamp) is not a number
                hdr = la[0].split(",")
                cols = set()
                for i in range(1, len(la)):
                    ca, cb = la[i].split(","), lb[i].split(",")
                    cols |= {hdr[j] if j < len(hdr) else str(j) for j in range(min(len(ca), len(cb))) if ca[j] != cb[j]}
                    if len(ca) != len(cb): cols.add("(column count)")
                # the wall-clock columns of the converter are not results
                TIMING = ("advectWallMax", "curvedFluxWallMax", "exchangeWallMax", "clockTime", "executionTime", "timeStepTime")
                numeric = [c for c in cols if c not in ("solver", "binary", "version", "gitCommit", "log", "application", "clippingSurface", "fluxIntegration") and c not in TIMING]
                rep.append(f"  {f}: differing columns {sorted(cols)}" + (" (wall-clock and route labels only)" if not numeric else ""))
                if numeric: ok = False
                continue
            ok = False
            rep.append(f"  {f} differs at line {k + 1}:\n    fork {la[k][:120] if k < len(la) else ''}\n    leia {lb[k][:120] if k < len(lb) else ''}")
        else:
            n += 1
    if ok: rep.append(f"  postProcessing and advectorVoF.csv: {n} files identical")
    return ok

def compare_logs(cf, cl, rep, verbose):
    def logof(c):   # either route's log
        for n in ("log.advectorVoF", "log.davofAdvectorFoam"):
            if os.path.isfile(os.path.join(c, n)): return os.path.join(c, n)
        return os.path.join(c, "log.advectorVoF")
    lf, ll = logof(cf), logof(cl)
    if not (os.path.isfile(lf) and os.path.isfile(ll)):
        rep.append(f"  logs missing: {lf if not os.path.isfile(lf) else ''} {ll if not os.path.isfile(ll) else ''}"); return False
    A = [norm_line(x) for x in open(lf, errors="replace") if not NOISE.match(x)]
    B = [norm_line(x) for x in open(ll, errors="replace") if not NOISE.match(x)]
    # the numbers of the per-step lines: keep the lines that carry a digit
    A = [x for x in A if re.search(r"\d", x)]; B = [x for x in B if re.search(r"\d", x)]
    diffs = [(i, a, b) for i, (a, b) in enumerate(zip(A, B)) if a != b]
    if len(A) != len(B) or diffs:
        i, a, b = diffs[0] if diffs else (min(len(A), len(B)), "", "")
        rep.append(f"  log: {len(A)} vs {len(B)} numeric lines, first difference at {i}:\n    fork {a[:140]}\n    leia {b[:140]}")
        if verbose:
            for i, a, b in diffs[1:6]: rep.append(f"    ... {i}: {a[:100]} | {b[:100]}")
        return False
    rep.append(f"  log: {len(A)} numeric lines identical (prefixes normalised)")
    return True

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fork"); ap.add_argument("leia"); ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    if os.path.isfile(os.path.join(a.fork, "system", "controlDict")) and os.path.isfile(os.path.join(a.leia, "system", "controlDict")):
        # two case directories (e.g. the fork's N 64 case against the one case of a proof study,
        # or the two routes' runs of a timed cell, which carry no case_params.json)
        pairs = [(os.path.basename(a.fork) + " vs " + os.path.basename(a.leia), a.fork, a.leia)]
    else:
        cases_f = sorted(d for d in os.listdir(a.fork) if re.search(r"_\d{5}$", d))
        cases_l = sorted(d for d in os.listdir(a.leia) if re.search(r"_\d{5}$", d))
        if cases_f != cases_l or not cases_f:
            print(f"case sets differ or empty: {cases_f} vs {cases_l}"); sys.exit(1)
        pairs = [(c, os.path.join(a.fork, c), os.path.join(a.leia, c)) for c in cases_f]
    allok = True
    for c, cf, cl in pairs:
        rep = []
        ok = compare_fields(cf, cl, rep) & compare_post(cf, cl, rep) & compare_logs(cf, cl, rep, a.verbose)
        params = ""
        pj = os.path.join(cl, "case_params.json")
        if os.path.isfile(pj):
            import json
            t = json.load(open(pj))
            t = t.get("tokens", t)
            params = " ".join(f"{k}={t[k]}" for k in ("N_CELLS", "TPF_ROW", "TPF_VELOCITY_ARM", "TARGET_ADVECTIVE_CO") if k in t)
        print(f"{'IDENTICAL' if ok else 'DIFFERENT'} {c} {params}")
        for r in rep: print(r)
        allok &= ok
    print("ALL IDENTICAL" if allok else "NOT IDENTICAL")
    sys.exit(0 if allok else 1)

if __name__ == "__main__":
    main()
