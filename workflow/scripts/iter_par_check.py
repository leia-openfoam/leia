#!/usr/bin/env python3
"""The parallel identity of the face-fraction corrector (leia davof-iter, 2026-11-09;
REPORT-swept.md "The iterative face-fraction corrector", P1): the fields alpha.water and
alphaf.water of the last written time of two cases (np 1 and np 8, the parallel one
reconstructed by the workflow) compared value by value; the max and the L1 of the
difference over the internal cells / faces and over the patch faces, and the number of
values differing by more than a tolerance.

    python3 workflow/scripts/iter_par_check.py <case_np1> <case_np8> [--tol 1e-12]

Both cases render from the same template with the same tokens except np, so the cell
and face order of the reconstructed mesh is the serial mesh's (reconstructPar maps the
processor fields back through the decomposition's addressing).
"""
import argparse
import glob
import os
import re
import sys


def last_time(case):
    ts = []
    for d in os.listdir(case):
        try:
            v = float(d)
        except ValueError:
            continue
        if os.path.isdir(os.path.join(case, d)) and v > 0:
            ts.append((v, d))
    if not ts:
        raise SystemExit("no written time in " + case)
    return max(ts)[1]


def read_field(path):
    """{'internal': [floats], 'patches': {name: [floats] or ('uniform', v)}} of an ascii
    OpenFOAM scalar field"""
    text = open(path, errors="replace").read()
    out = {"internal": None, "patches": {}}

    def parse_list(s, pos):
        m = re.compile(r"\s*(\d+)\s*\(", re.S).match(s, pos)
        if not m:
            return None, pos
        n = int(m.group(1))
        start = m.end()
        end = s.index(")", start)
        vals = [float(x) for x in s[start:end].split()]
        assert len(vals) == n, (len(vals), n, path)
        return vals, end + 1

    m = re.search(r"internalField\s+nonuniform\s+List<scalar>", text)
    if m:
        out["internal"], _ = parse_list(text, m.end())
    else:
        m = re.search(r"internalField\s+uniform\s+([^;]+);", text)
        out["internal"] = ("uniform", float(m.group(1)))
    bpos = text.find("boundaryField")
    for pm in re.finditer(r"\n\s{4}(\w+)\s*\{", text[bpos:]):
        name = pm.group(1)
        blk_start = bpos + pm.end()
        depth = 1
        i = blk_start
        while depth and i < len(text):
            if text[i] == "{":
                depth += 1
            elif text[i] == "}":
                depth -= 1
            i += 1
        blk = text[blk_start:i]
        vm = re.search(r"value\s+nonuniform\s+List<scalar>", blk)
        if vm:
            vals, _ = parse_list(blk, vm.end())
            out["patches"][name] = vals
        else:
            vm = re.search(r"value\s+uniform\s+([^;]+);", blk)
            if vm:
                out["patches"][name] = ("uniform", float(vm.group(1)))
    return out


def compare(a, b, tol):
    if isinstance(a, tuple) or isinstance(b, tuple):
        if isinstance(a, tuple) and isinstance(b, tuple):
            d = abs(a[1] - b[1])
            return d, d, int(d > tol), 1
        lst = b if isinstance(a, tuple) else a
        u = a[1] if isinstance(a, tuple) else b[1]
        ds = [abs(x - u) for x in lst]
    else:
        assert len(a) == len(b), (len(a), len(b))
        ds = [abs(x - y) for x, y in zip(a, b)]
    return max(ds), sum(ds) / len(ds), sum(1 for d in ds if d > tol), len(ds)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case1")
    ap.add_argument("case8")
    ap.add_argument("--tol", type=float, default=1e-12)
    a = ap.parse_args()
    t1, t8 = last_time(a.case1), last_time(a.case8)
    print("last times: %s (np1) %s (np8)" % (t1, t8))
    ok = True
    for name in ("alpha.water", "alphaf.water"):
        f1 = read_field(os.path.join(a.case1, t1, name))
        f8 = read_field(os.path.join(a.case8, t8, name))
        mx, l1, nbad, n = compare(f1["internal"], f8["internal"], a.tol)
        print("%-13s internal: n %d, max |diff| %.3e, L1 %.3e, above %g: %d" % (name, n, mx, l1, a.tol, nbad))
        ok = ok and nbad == 0
        for p in sorted(f1["patches"]):
            if p not in f8["patches"]:
                print("  patch %s missing in np8" % p)
                ok = False
                continue
            mx, l1, nbad, n = compare(f1["patches"][p], f8["patches"][p], a.tol)
            print("  patch %-16s n %d, max |diff| %.3e, above %g: %d" % (p, n, mx, a.tol, nbad))
            ok = ok and nbad == 0
    print("PARALLEL IDENTITY:", "held (every value within %g)" % a.tol if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
