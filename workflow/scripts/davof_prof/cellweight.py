"""cellweight.py <case> <time> <W> [<layers>]: the cell weights of an interface-weighted
decomposition (Claude, 2026-11-04): a cell within <layers> face-neighbour layers of a partial
cell (1e-8 < alpha < 1 - 1e-8) at <time> gets the weight W, every other cell 1; written as the
volScalarField cellWeight in <case>/<time>/ (serial, undecomposed case; no numpy)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamfield import read_labels, read_scalars, band_weights, write_weights
case, t, W = sys.argv[1], sys.argv[2], float(sys.argv[3])
layers = int(sys.argv[4]) if len(sys.argv) > 4 else 2
own = read_labels(os.path.join(case, "constant/polyMesh/owner"))
nei = read_labels(os.path.join(case, "constant/polyMesh/neighbour"))
a, u = read_scalars(os.path.join(case, t, "alpha.water"))
if a is None: a = [u]*(max(own) + 1)
w, nb, nw = band_weights(own, nei, a, W, layers)
write_weights(os.path.join(case, t, "cellWeight"), t, w)
print("cellWeight: %d cells, %d partial, %d weighted (W = %g, %d layers), total weight %g, interface share of the weight %.3f" % (len(w), nb, nw, W, layers, sum(w), (sum(w) - (len(w) - nw))/sum(w)))
