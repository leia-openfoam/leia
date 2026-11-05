"""The per-face-class statistics of a log (key sweptClassStats): the means per step over the steps
(the first and the last step excluded), printed as a markdown table.
    python3 r3_cls.py <log.advectorVoF> [label]"""
import re, sys, statistics as st
log = open(sys.argv[1]).read()
lab = sys.argv[2] if len(sys.argv) > 2 else sys.argv[1]
rx = re.compile(r"class (no-sweep|donor-only|multi-cell|boundary): faces (\d+), touched cells (\d+), pieces (\d+), volume clips (\d+), area clips (\d+), triangles (\d+), time volume ([0-9.e+-]+) area ([0-9.e+-]+) trace ([0-9.e+-]+) pre-test ([0-9.e+-]+)")
rows = {}
for m in rx.finditer(log):
    rows.setdefault(m.group(1), []).append([float(x) for x in m.groups()[1:]])
steps = [float(x) for x in re.findall(r"timeStepTime = (\S+) s", log)]
print("%s: %d steps, step mean %.4f median %.4f s" % (lab, len(steps), st.mean(steps[1:-1]) if len(steps) > 2 else float('nan'), st.median(steps) if steps else float('nan')))
print("| class | faces/step | touched cells/face | pieces/face | volume clips/face | area clips/face | triangles/face | t volume s/step | t area | t trace | t pre-test | sum | share of the four |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
tot = 0.0; per = {}
for k in ("no-sweep", "donor-only", "multi-cell", "boundary"):
    r = rows.get(k, [])
    if len(r) > 2: r = r[1:-1]
    if not r: continue
    m = [st.mean(c) for c in zip(*r)]
    s = m[6] + m[7] + m[8] + m[9]; per[k] = (m, s); tot += s
for k, (m, s) in per.items():
    f = max(m[0], 1e-300)
    print("| %s | %.0f | %.2f | %.2f | %.1f | %.1f | %.1f | %.4f | %.4f | %.4f | %.4f | %.4f | %.1f %% |" % (k, m[0], m[1]/f, m[2]/f, m[3]/f, m[4]/f, m[5]/f, m[6], m[7], m[8], m[9], s, 100*s/tot if tot else 0))
print("the four timed parts: %.4f s/step in total" % tot)
