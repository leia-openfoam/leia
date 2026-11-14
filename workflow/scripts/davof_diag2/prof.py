#!/usr/bin/env python3
"""OpenFOAM profiling (uniform/profiling of the last write time) of a parallel case: the trigger
tree per rank, summarised as mean and max over the ranks of totalTime (seconds) and the calls."""
import glob, os, re, sys

rx_block = re.compile(r"trigger(\d+)\s*\{(.*?)\}", re.S)

def short(desc):
    d = desc.strip('"')
    d = re.sub(r"\[with .*", "", d)
    d = re.sub(r"\(const .*", "", d)
    d = re.sub(r"^(void|bool|Foam::tmp<[^>]*>|Foam::label|Foam::scalar)\s+", "", d)
    return d.replace("Foam::", "").strip()

def read(path):
    s = open(path).read()
    out = {}
    for m in rx_block.finditer(s):
        body = m.group(2)
        g = lambda k: (re.search(k + r"\s+([^;]+);", body) or [None, None])[1]
        out[int(m.group(1))] = dict(id=int(m.group(1)), parent=int(g("parentId")) if g("parentId") is not None else -1,
                                    desc=short(g("description") or ""), calls=int(g("calls")), total=float(g("totalTime")),
                                    child=float(g("childTime")))
    return out

def case_profile(case):
    """-> {desc: (calls, mean total over ranks, max total over ranks, parent desc)}"""
    procs = sorted(glob.glob(case + "/processor*")) or [case]
    per = []
    for p in procs:
        times = [d for d in glob.glob(p + "/*") if os.path.exists(d + "/uniform/profiling")]
        if not times:
            continue
        last = max(times, key=lambda d: float(os.path.basename(d)))
        per.append(read(last + "/uniform/profiling"))
    if not per:
        return {}
    res = {}
    keys = {}
    for tr in per[0].values():
        parent = per[0][tr["parent"]]["desc"] if tr["parent"] in per[0] else "-"
        keys[tr["id"]] = (tr["desc"], parent)
    for i, (desc, parent) in keys.items():
        tot = [r[i]["total"] for r in per if i in r]
        calls = [r[i]["calls"] for r in per if i in r]
        # the rank name is part of the time.run() description: strip it
        d = re.sub(r"processor\d+$", "processorN", desc)
        res[d] = (max(calls), sum(tot)/len(tot), max(tot), parent)
    return res

if __name__ == "__main__":
    for c in sys.argv[1:]:
        print("##", c)
        for d, (n, mean, mx, par) in case_profile(c).items():
            print("%-95s calls %7d mean %10.3f max %10.3f  <- %s" % (d[:95], n, mean, mx, par[:50]))
