"""the in-code timers of the swept advection from a comparison log (max over ranks per step)"""
import sys, re, statistics as st
for path in sys.argv[1:]:
    log = open(path + "/log.advectorVoF").read()
    steps = [float(x) for x in re.findall(r"timeStepTime = (\S+) s", log)]
    rows = re.findall(r"wall time max over (\d+) ranks: advect (\S+) s, curvedFlux (\S+) s, exchange (\S+) s \(wait (\S+), stencil (\S+), map (\S+), pack (\S+), transfer (\S+), unpack (\S+)\)", log)
    ic = [int(x) for x in re.findall(r"davof: (\d+) interface cells", log)]
    pieces = [int(x) for x in re.findall(r"cells: pieces (\d+)", log)]
    partial = [int(x) for x in re.findall(r"partial cells (\d+)", log)]
    wispno = [int(x) for x in re.findall(r"wisps without a neighbour surface (\d+)", log)]
    wclip = [int(x) for x in re.findall(r"pieces clipped by a neighbour's plane (\d+)", log)]
    faces = [int(x) for x in re.findall(r"swept faces (\d+)", log)]
    tris = [int(x) for x in re.findall(r"clipped triangles (\d+)", log)]
    cand = [int(x) for x in re.findall(r"candidates added (\d+)", log)]
    second = [int(x) for x in re.findall(r"\((\d+) second ring", log)]
    remote = [int(x) for x in re.findall(r"remote cells (\d+)", log)]
    print(path.split("/")[-2] + "/" + path.split("/")[-1], "steps", len(steps), "mean step %.3f" % st.mean(steps))
    if rows:
        cols = ["advect", "curvedFlux", "exchange", "wait", "stencil", "map", "pack", "transfer", "unpack"]
        vals = [[float(r[i+1]) for r in rows] for i in range(9)]
        print("  ranks", rows[0][0], " ".join("%s %.4f" % (c, st.mean(v)) for c, v in zip(cols, vals)), " (means over the turn, s/step)")
        print("  advect/step %.3f of %.3f -> %.0f%%; exchange %.1f%% of advect; wait %.1f%%" % (st.mean(vals[0]), st.mean(steps), 100*st.mean(vals[0])/st.mean(steps), 100*st.mean(vals[2])/st.mean(vals[0]), 100*st.mean(vals[3])/st.mean(vals[0])))
    def q(a): return "%d / %d / %d (min/mean/max)" % (min(a), st.mean(a), max(a)) if a else "-"
    print("  interface cells", q(ic), "| swept faces", q(faces), "| second ring", q(second), "| candidates added", q(cand))
    print("  pieces", q(pieces), "| clipped triangles", q(tris), "| partial cells (wisps met)", q(partial), "| wisps w/o neighbour surface", q(wispno), "| pieces clipped by neighbour plane", q(wclip), "| remote cells", q(remote))
    # correlation of the step time with the counters
    import math
    def corr(a, b):
        n = min(len(a), len(b)); a = a[:n]; b = b[:n]
        ma, mb = st.mean(a), st.mean(b)
        sa = math.sqrt(sum((x-ma)**2 for x in a)); sb = math.sqrt(sum((x-mb)**2 for x in b))
        return sum((x-ma)*(y-mb) for x, y in zip(a, b))/(sa*sb) if sa and sb else float('nan')
    print("  corr(step time, .): partial %.2f, wispsNoSurf %.2f, pieces %.2f, tris %.2f, faces %.2f, interface %.2f" % (corr(steps, partial), corr(steps, wispno), corr(steps, pieces), corr(steps, tris), corr(steps, faces), corr(steps, ic)))
