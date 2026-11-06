#!/usr/bin/env python3
"""The field comparison of two runs of the matrix (2026-11-07): the largest |dalpha| and
|dalphaf| over the cells and the internal faces of the last written time (both serial, or
both decomposed on the same ranks: processor by processor), the number of differing values,
and whether the two are bit for bit (every value equal). Prints one CSV line
    ref,case,time,maxDalpha,nDalpha,maxDalphaf,nDalphaf,identical
    python3 enright_compare.py <ref case dir> <case dir> [field name alpha.water] [alphaf.water]
"""
import glob, os, re, sys


def field_values(fn):
    txt = open(fn).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", txt)
    if m:
        start = m.end(); end = txt.index("\n)", start)
        return [float(x) for x in txt[start:end].split()]
    m = re.search(r"internalField\s+uniform\s+(\S+);", txt)
    if m:
        return [float(m.group(1))]
    return None


def last_time(case):
    ds = glob.glob(case + "/processor0/[0-9]*") or glob.glob(case + "/[0-9]*")
    ts = []
    for d in ds:
        try: ts.append((float(os.path.basename(d)), os.path.basename(d)))
        except ValueError: pass
    ts.sort()
    return ts[-1][1] if ts else None


def files(case, t, name):
    pr = sorted(glob.glob(case + "/processor*"), key=lambda p: int(p.split("processor")[-1]))
    if pr:
        return [p + "/" + t + "/" + name for p in pr]
    return [case + "/" + t + "/" + name]


def main():
    ref, case = sys.argv[1], sys.argv[2]
    names = sys.argv[3:] if len(sys.argv) > 3 else ["alpha.water", "alphaf.water"]
    tr, tc = last_time(ref), last_time(case)
    out = [os.path.basename(ref.rstrip("/")), os.path.basename(case.rstrip("/")), "%s/%s" % (tr, tc)]
    identical = True
    for name in names:
        mx = 0.0; nd = 0; ok = True
        fr, fc = files(ref, tr, name), files(case, tc, name)
        if len(fr) != len(fc) or not all(os.path.isfile(f) for f in fr + fc):
            out += ["nan", "nan"]; identical = False; continue
        for a, b in zip(fr, fc):
            va, vb = field_values(a), field_values(b)
            if va is None or vb is None or len(va) != len(vb):
                ok = False; break
            for x, y in zip(va, vb):
                if x != y:
                    nd += 1; mx = max(mx, abs(x - y))
        if not ok:
            out += ["nan", "nan"]; identical = False; continue
        out += ["%.3e" % mx, str(nd)]
        if nd: identical = False
    out.append("1" if identical else "0")
    print(",".join(out))


if __name__ == "__main__":
    main()
