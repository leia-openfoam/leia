import re, glob, os
S = "/work/scratch/tm83tomy/davof/agents/pointiso/leia/studies/davof/"
def lst(path):
    t = open(path).read()
    mu = re.search(r"internalField\s+uniform\s+(\S+);", t)
    if mu:
        return float(mu.group(1))
    i = t.index("internalField")
    m = re.compile(r"\s(\d+)\s*\(").search(t, i)
    n = int(m.group(1)); body = t[m.end():]
    body = body[:body.index(")")]
    return [float(x) for x in body.split()][:n]
def labels(path):
    t = open(path).read()
    i = t.index("object")
    m = re.compile(r"\s(\d+)\s*\(").search(t, i)
    n = int(m.group(1)); body = t[m.end():]
    body = body[:body.index(")")]
    return [int(x) for x in body.split()][:n]
ok = True
c1 = S + "contactLineD_N128/contactLine2D_00003"
c8 = S + "cl_gr2_np8/contactLine2D_00000"
import json
for c in (c1, c8):
    tk = json.load(open(c + "/case_params.json"))["tokens"]
    print(c, tk["TPF_ROW"], tk["TPF_VELOCITY_ARM"], tk["TARGET_ADVECTIVE_CO"], tk["N_CELLS"])
times = sorted([d for d in os.listdir(c1) if re.match(r"^[0-9.eE+-]+$", d) and d != "0"], key=float)
for t in times:
    a1 = lst(c1 + "/%s/alpha.water" % t)
    a8 = [None] * len(a1)
    for p in sorted(glob.glob(c8 + "/processor*")):
        al = lst(p + "/%s/alpha.water" % t)
        ad = labels(p + "/constant/polyMesh/cellProcAddressing")
        if isinstance(al, float):
            al = [al] * len(ad)
        for v, g in zip(al, ad):
            a8[g] = v
    assert None not in a8
    md = max(abs(x - y) for x, y in zip(a1, a8))
    print("time %s, %d cells, max |alpha(np1) - alpha(np8)| = %.3e, interface cells %d" % (t, len(a1), md, sum(1 for x in a1 if 1e-8 < x < 1 - 1e-8)))
    ok = ok and md <= 1e-12
print("Gr2", "PASS" if ok else "FAIL")