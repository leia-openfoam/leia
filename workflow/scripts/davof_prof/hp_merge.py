"""merge the per-run analysis outputs of prof/an into the CSVs of the report (Claude, 2026-11-04)
    python3 hp_merge.py <prof dir> <csv out dir>
writes hp_categories.csv (run x category: s/step/rank, pct), hp_ranks.csv (per rank of the parallel runs),
hp_hotspots.csv (top 25 exclusive functions of every run), hp_timing.csv (the jobs' timing.csv rows)
"""
import sys, os, glob, csv, re
P, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
runs = sorted(glob.glob(P + "/an/*.cat.csv"))
cats = []
table = {}
tot = {}
for f in runs:
    run = os.path.basename(f)[:-8]
    rows = list(csv.DictReader(open(f)))
    for r in rows:
        c = r["category"]
        if c not in cats: cats.append(c)
        table[(run, c)] = r
runnames = [os.path.basename(f)[:-8] for f in runs]
with open(OUT + "/hp_categories.csv", "w") as f:
    f.write("category," + ",".join("%s_s_per_step_per_rank,%s_pct,%s_max_rank_s,%s_mpi_at_site_s" % (r, r, r, r) for r in runnames) + "\n")
    for c in cats:
        f.write('"%s"' % c)
        for r in runnames:
            x = table.get((r, c))
            f.write(",%s,%s,%s,%s" % ((x["s_per_step_per_rank"], x["pct"], x.get("max_rank_s", ""), x.get("mpi_at_site_s", "")) if x else ("", "", "", "")))
        f.write("\n")
with open(OUT + "/hp_ranks.csv", "w") as f:
    first = True
    for g in sorted(glob.glob(P + "/an/par*.ranks.csv")):
        run = os.path.basename(g)[:-10]
        lines = open(g).read().splitlines()
        if first: f.write("run," + lines[0] + "\n"); first = False
        for l in lines[1:]: f.write(run + "," + l + "\n")
with open(OUT + "/hp_hotspots.csv", "w") as f:
    f.write("run,rank_in_run,function,module,file,line,excl_s,excl_pct,incl_s,incl_pct,llc_misses\n")
    for g in sorted(glob.glob(P + "/an/*.flat.csv")):
        run = os.path.basename(g)[:-9]
        for i, r in enumerate(csv.DictReader(open(g))):
            if i >= 25: break
            f.write('%s,%d,"%s",%s,%s,%s,%s,%s,%s,%s,%s\n' % (run, i + 1, r["function"].replace('"', "'")[:160], r["module"], os.path.basename(r["file"]), r["line"], r["excl_s"], r["excl_pct"], r["incl_s"], r["incl_pct"], r["llc_miss_samples"]))
with open(OUT + "/hp_timing.csv", "w") as f:
    first = True
    for g in sorted(glob.glob(P + "/runs/*/timing.csv")):
        lines = open(g).read().splitlines()
        if first: f.write(lines[0] + "\n"); first = False
        for l in lines[1:]: f.write(l + "\n")
print("wrote", OUT, runnames)
