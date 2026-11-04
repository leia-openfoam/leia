"""per-rank cell and partial-cell counts of a decomposed case (ascii or binary alpha fields; no numpy)"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamfield import read_scalars, read_labels
case = sys.argv[1]; t = sys.argv[2] if len(sys.argv) > 2 else "0"
procs = sorted(glob.glob(case + "/processor*"), key=lambda p: int(p.split("processor")[-1]))
ni = []; nc = []
for p in procs:
    a, u = read_scalars(p + "/" + t + "/alpha.water")
    if a is None:
        own = read_labels(p + "/constant/polyMesh/owner"); n = max(own) + 1 if own else 0
        nc.append(n); ni.append(0)
    else:
        nc.append(len(a)); ni.append(sum(1 for x in a if 1e-8 < x < 1 - 1e-8))
mean = sum(ni)/len(ni)
print("%s t=%s: ranks %d, cells/rank %d..%d, partial cells (1e-8<a<1-1e-8) per rank: total %d, mean %.1f, max %d, min %d, ranks with 0: %d, max/mean %.2f" % (case.split("/")[-2] + "/" + case.split("/")[-1], t, len(ni), min(nc), max(nc), sum(ni), mean, max(ni), min(ni), sum(1 for x in ni if x == 0), max(ni)/max(mean, 1e-9)))
print("  per rank:", " ".join(str(x) for x in ni))
