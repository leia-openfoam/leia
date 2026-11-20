import re, glob, sys
C = sys.argv[1]
def patches(path):
    t = open(path).read()
    out = {}
    for m in re.finditer(r"(procBoundary\d+to\d+)\s*\{(.*?)\n    \}", t, re.S):
        body = m.group(2)
        v = re.search(r"value\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(([^)]*)\)", body, re.S)
        out[m.group(1)] = [float(x) for x in v.group(2).split()] if v else None
        out[m.group(1) + "_type"] = re.search(r"type\s+(\w+)", body).group(1)
    return out
P = {i: patches(C + "processor%d/0/phi" % i) for i in range(8)}
worst = 0
for a in range(8):
    for b in range(8):
        k = "procBoundary%dto%d" % (a, b)
        k2 = "procBoundary%dto%d" % (b, a)
        if k in P[a] and k2 in P[b] and P[a][k] and P[b][k2]:
            s = [x + y for x, y in zip(P[a][k], P[b][k2])]
            m = max(abs(x) for x in s)
            print(k, "faces", len(s), "max |phi_a + phi_b| =", m, "max |phi| =", max(abs(x) for x in P[a][k]), "type", P[a][k + "_type"])
            worst = max(worst, m)
print("worst", worst)