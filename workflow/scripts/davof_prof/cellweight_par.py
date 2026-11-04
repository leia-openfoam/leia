"""cellweight_par.py <case> <time> <W> [<layers>]: the cell weights (cellweight.py) of a DECOMPOSED
case, written per processor at <time> from each processor's alpha.water and polyMesh (the layers
grown inside each processor; a band cell's layer across a processor boundary is missed, at most
LAYERS cells deep). No numpy. Claude, 2026-11-04."""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamfield import read_labels, read_scalars, band_weights, write_weights
case, t, W = sys.argv[1], sys.argv[2], float(sys.argv[3])
layers = int(sys.argv[4]) if len(sys.argv) > 4 else 2
tot = totb = totw = 0; per = []
for p in sorted(glob.glob(case + "/processor*"), key=lambda p: int(p.split("processor")[-1])):
    own = read_labels(p + "/constant/polyMesh/owner"); nei = read_labels(p + "/constant/polyMesh/neighbour")
    nc = max(own) + 1 if own else 0
    a, u = read_scalars(p + "/" + t + "/alpha.water")
    if a is None: a = [u]*nc
    w, nb, nw = band_weights(own, nei, a, W, layers)
    write_weights(p + "/" + t + "/cellWeight", t, w)
    tot += len(w); totb += nb; totw += nw; per.append(nb)
print("cellWeight per processor at t=%s: %d cells, %d partial, %d weighted (W %g, %d layers); partial per rank max %d mean %.1f" % (t, tot, totb, totw, W, layers, max(per), sum(per)/len(per)))
