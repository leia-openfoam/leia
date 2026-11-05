"""The two HPCToolkit profiles side by side (Claude, 2026-11-04): the first profile (fork 720cfc6,
prof/an/<run>.*) against the reference setting (fork 407f969, prof/an2/<run>.*), the same runs
(serial ring32_full, ring64_full, ring64_20, ring128_20; par8 ring64_full; par64 ring128_full).
    python3 hp_side.py <prof dir> <csv out dir>
writes hp2_categories.csv (run x group: s/step/rank old and new), hp2_ranks.csv, hp2_hotspots.csv,
hp2_timing.csv, and prints the markdown tables (a) to (d)."""
import sys, os, glob, csv, re, statistics
P, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
GROUPS = [
    ("pieces: the clip by the cell planes and the PLIC plane", ["pieces: clip"]),
    ("pieces: the cap closure and the volume", ["pieces: cap", "pieces: sweptPieceVolume other"]),
    ("face update: the pre-image fan clipped (sweptPieceArea)", ["face update: pre-image"]),
    ("face update: the donor plane", ["face update: donor"]),
    ("swept cells: the ring search", ["swept cells: ring"]),
    ("swept cells: the cell planes, state, wisp rule", ["swept cells: cell planes"]),
    ("swept cells: the pre-test", ["swept cells: the pre-test"]),
    ("swept cells: the shortcut and the check (off)", ["swept cells: the state shortcut", "pre-test check"]),
    ("swept cells: the pre-image fan and the rest", ["swept cells: pre-image fan", "swept cells: other", "one-donor"]),
    ("trace: the velocity at the traced points (findTetrahedron)", ["trace: velocity interpolate"]),
    ("trace: the cell search (ring, octree)", ["trace: cell search"]),
    ("trace: the midpoint trace, the triangle list", ["trace: midpoint", "trace: triangle"]),
    ("the point velocity (volPointInterpolation of U)", ["curvedFlux: point velocity"]),
    ("curvedFlux: the band, the wisp surfaces, the exchanges, the face loop", ["curvedFlux:"]),
    ("the remote exchange (serial: the stencil setup)", ["remote exchange"]),
    ("the surface build (fit, anchor, volume rule)", ["surface build"]),
    ("the reconstruction (davof, ring normal, positioning)", ["reconstruction"]),
    ("the bounding", ["bounding"]),
    ("isoAdvector: timeIntegratedFlux", ["timeIntegratedFlux"]),
    ("isoAdvector: advect, the rest (alpha update, counters)", ["advect other"]),
    ("setFlow", ["flow setup"]),
    ("the Courant number, other function objects", ["Courant"]),
    ("I/O and start-up", ["I/O"]),
    ("the solver loop (main)", ["main:"]),
    ("MPI library", ["MPI library"]),
    ("Pstream wrappers", ["Pstream"]),
    ("(uncategorised)", ["(uncategorised)"]),
]
def group(cat):
    for g, subs in GROUPS:
        if any(s in cat for s in subs): return g
    return "(other: %s)" % cat
def cats(path):
    out = {}; tot = None
    for r in csv.DictReader(open(path)):
        if r["category"] == "(total)": tot = r; continue
        g = group(r["category"])
        e = out.setdefault(g, [0.0, 0.0, 0.0])
        e[0] += float(r["s_per_step_per_rank"]); e[1] += float(r["pct"]); e[2] += float(r["cputime_s"])
    return out, tot
runs = ["serial_ring32_full", "serial_ring64_full", "serial_ring64_20", "serial_ring128_20", "par8_ring64_full", "par64_ring128_full"]
table = {}
for run in runs:
    for lab, d in (("old", "an"), ("new", "an2")):
        f = "%s/%s/%s.cat.csv" % (P, d, run)
        if os.path.isfile(f): table[(run, lab)] = cats(f)
with open(OUT + "/hp2_categories.csv", "w") as f:
    f.write("group," + ",".join("%s_%s_s_per_step_per_rank,%s_%s_pct" % (r, l, r, l) for r in runs for l in ("old", "new")) + "\n")
    for g, _ in GROUPS:
        f.write('"%s"' % g)
        for r in runs:
            for l in ("old", "new"):
                x = table.get((r, l)); v = x[0].get(g) if x else None
                f.write(",%.4f,%.2f" % (v[0], v[1]) if v else ",,")
        f.write("\n")
    f.write("(total)")
    for r in runs:
        for l in ("old", "new"):
            x = table.get((r, l)); f.write(",%s,100" % x[1]["s_per_step_per_rank"] if x else ",,")
    f.write("\n")
def md_a(run, title):
    o, n = table.get((run, "old")), table.get((run, "new"))
    if not o or not n: print("(%s: missing %s)" % (run, "old" if not o else "new")); return
    print("\n%s\n" % title)
    print("| category | before (720cfc6) s/step | % | after (407f969) s/step | % | ratio |")
    print("|---|---|---|---|---|---|")
    for g, _ in GROUPS:
        a, b = o[0].get(g), n[0].get(g)
        if not a and not b: continue
        av, bv = (a[0] if a else 0.0), (b[0] if b else 0.0)
        ap, bp = (a[1] if a else 0.0), (b[1] if b else 0.0)
        print("| %s | %.4f | %.1f | %.4f | %.1f | %s |" % (g, av, ap, bv, bp, ("%.2f" % (bv/av)) if av > 0 else "-"))
    print("| total | %s | 100 | %s | 100 | %.2f |" % (o[1]["s_per_step_per_rank"], n[1]["s_per_step_per_rank"], float(n[1]["s_per_step_per_rank"])/float(o[1]["s_per_step_per_rank"])))
md_a("serial_ring64_full", "### (a) The serial breakdown per step, N = 64, one turn, before and after")
md_a("serial_ring32_full", "### (a') N = 32, one turn")
md_a("par8_ring64_full", "### (c') N = 64 on 8 ranks, one turn (per rank)")
md_a("par64_ring128_full", "### (c'') N = 128 on 64 ranks, one turn (per rank)")
# (b) the scaling with N: new only, the swept groups
print("\n### (b) The scaling with N (after); the ratio to before in brackets\n")
print("| group | N 32 full | N 64 full | x | N 64, 20 | N 128, 20 | x |")
print("|---|---|---|---|---|---|---|")
def v(run, g, lab="new"):
    x = table.get((run, lab)); return x[0].get(g, [0.0])[0] if x else float("nan")
for g, _ in GROUPS[:20]:
    a, b, c, d = v("serial_ring32_full", g), v("serial_ring64_full", g), v("serial_ring64_20", g), v("serial_ring128_20", g)
    if a == 0 and b == 0: continue
    print("| %s | %.4f | %.4f | %.2f | %.4f | %.4f | %.2f |" % (g, a, b, b/a if a else float("nan"), c, d, d/c if c else float("nan")))
# (c) ranks
with open(OUT + "/hp2_ranks.csv", "w") as f:
    first = True
    for lab, d in (("old", "an"), ("new", "an2")):
        for g in sorted(glob.glob("%s/%s/par*.ranks.csv" % (P, d))):
            run = os.path.basename(g)[:-10]; lines = open(g).read().splitlines()
            if first: f.write("profile,run," + lines[0] + "\n"); first = False
            for l in lines[1:]: f.write("%s,%s,%s\n" % (lab, run, l))
print("\n### (c) Parallel: MPI against compute, the imbalance (CPUTIME over the turn per rank)\n")
print("| | N 64, 8 ranks before | after | N 128, 64 ranks before | after |")
print("|---|---|---|---|---|")
def ranks(run, d):
    f = "%s/%s/%s.ranks.csv" % (P, d, run)
    if not os.path.isfile(f): return None
    rows = list(csv.DictReader(open(f)))
    tot = [float(r["cputime_s"]) for r in rows]; mpi = [float(r["mpi_s"]) for r in rows]; comp = [float(r["compute_s"]) for r in rows]
    steps = 697 if "64" in run.split("_")[1] else 1401
    return dict(step=sum(tot)/len(tot)/steps, mpi_share=100*sum(mpi)/sum(tot), mpi="%.0f / %.0f / %.0f" % (statistics.mean(mpi), max(mpi), min(mpi)),
                comp="%.0f / %.0f (%.2f)" % (max(comp), statistics.mean(comp), max(comp)/statistics.mean(comp)), idle=sum(1 for c in comp if c < 0.1*max(comp)), n=len(comp),
                balanced=sum(comp)/len(comp)/steps, busiest=max(range(len(comp)), key=lambda i: comp[i]), comps=comp, steps=steps)
R = {(run, lab): ranks(run, d) for run in ("par8_ring64_full", "par64_ring128_full") for lab, d in (("old", "an"), ("new", "an2"))}
def cell(run, lab, k, fmt="%s"):
    x = R.get((run, lab)); return (fmt % x[k]) if x else "-"
for key, name, fmt in (("step", "step (CPUTIME per rank / steps), s", "%.3f"), ("mpi_share", "MPI library, % of all CPUTIME", "%.1f"), ("mpi", "MPI per rank mean / max / min, s", "%s"),
                       ("comp", "compute per rank max / mean (max/mean), s", "%s"), ("idle", "ranks with compute < 10 % of the max", "%d"), ("balanced", "the balanced step (sum of compute / ranks / steps), s", "%.3f")):
    print("| %s | %s | %s | %s | %s |" % (name, cell("par8_ring64_full", "old", key, fmt), cell("par8_ring64_full", "new", key, fmt), cell("par64_ring128_full", "old", key, fmt), cell("par64_ring128_full", "new", key, fmt)))
for run in ("par8_ring64_full", "par64_ring128_full"):
    x = R.get((run, "new"))
    if x: print("- %s after: compute per rank s: %s" % (run, " ".join("%.0f" % c for c in x["comps"])))
# (d) hotspots
with open(OUT + "/hp2_hotspots.csv", "w") as f:
    f.write("profile,run,rank_in_run,function,module,file,line,excl_s,excl_pct,incl_s,incl_pct,llc_misses\n")
    for lab, d in (("old", "an"), ("new", "an2")):
        for g in sorted(glob.glob("%s/%s/*.flat.csv" % (P, d))):
            run = os.path.basename(g)[:-9]
            for i, r in enumerate(csv.DictReader(open(g))):
                if i >= 25: break
                f.write('%s,%s,%d,"%s",%s,%s,%s,%s,%s,%s,%s,%s\n' % (lab, run, i + 1, r["function"].replace('"', "'")[:160], r["module"], os.path.basename(r["file"]), r["line"], r["excl_s"], r["excl_pct"], r["incl_s"], r["incl_pct"], r["llc_miss_samples"]))
print("\n### (d) The top ten exclusive hotspots after (serial, N = 64, one turn)\n")
f = "%s/an2/serial_ring64_full.flat.csv" % P
if os.path.isfile(f):
    rows = list(csv.DictReader(open(f)))
    T = sum(float(r["excl_s"]) for r in rows)
    print("| # | %% | s | function (file:line) | incl s |")
    print("|---|---|---|---|---|")
    for i, r in enumerate(rows[:14]):
        print("| %d | %s | %s | %s (%s:%s) | %s |" % (i + 1, r["excl_pct"], r["excl_s"], r["function"][:110], os.path.basename(r["file"]), r["line"], r["incl_s"]))
    print("total exclusive %.1f s" % T)
    # the hot lines of the fork
    lf = "%s/an2/serial_ring64_full.lines.csv" % P
    print("\nthe hot source lines (fork + OpenFOAM):")
    for i, r in enumerate(csv.DictReader(open(lf))):
        if i >= 20: break
        print("  %s:%s %s s %s %% %s" % (os.path.basename(r["file"]), r["line"], r["excl_s"], r["excl_pct"], r["function"][:80]))
with open(OUT + "/hp2_timing.csv", "w") as f:
    first = True
    for g in sorted(glob.glob(P + "/runs2/*/timing.csv")):
        lines = open(g).read().splitlines()
        if first: f.write(lines[0] + "\n"); first = False
        for l in lines[1:]: f.write(l + "\n")
print("\nwrote", OUT, "hp2_categories.csv hp2_ranks.csv hp2_hotspots.csv hp2_timing.csv")
