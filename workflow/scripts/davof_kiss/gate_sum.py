import csv, sys, collections
S = "/work/scratch/tm83tomy/davof/agents/diag1/leia/studies/davof/sw_tables/"
for v in ("S1S2low", "S1S2all"):
    rows = list(csv.DictReader(open(S + f"kiss_swg_{v}_table.csv")))
    print(f"\n### {v}: per arm and Co (np1 and np8 pooled), cases / P1 ok / L1(alpha_f) min..max / Linf max / min alpha min / limiter faces max")
    print("| arm | Co | cases | P1 ok | L1 min..max | Linf max | min alpha | limiter faces max |")
    print("|---|---|---|---|---|---|---|---|")
    g = collections.defaultdict(list)
    for r in rows:
        g[(r["arm"], float(r["Co"]))].append(r)
    def fl(x):
        try: return float(x)
        except: return None
    for k in sorted(g):
        rs = g[k]
        l1 = [fl(r["L1af"]) for r in rs if fl(r["L1af"]) is not None]
        li = [fl(r["Linfaf"]) for r in rs if fl(r["Linfaf"]) is not None]
        ma = [fl(r["minAlpha"]) for r in rs if fl(r["minAlpha"]) is not None]
        lf = [fl(r["limFaces"]) for r in rs if fl(r["limFaces"]) is not None]
        ok = sum(1 for r in rs if "P1 ok" in r["verdict"] and not r["verdict"].startswith("OPEN"))
        op = any(r["verdict"].startswith("OPEN") for r in rs)
        print(f"| {k[0]} | {k[1]}{' (open cell)' if op else ''} | {len(rs)} | {ok} | {min(l1):.2g}..{max(l1):.2g} | {max(li):.2g} | {min(ma):.2g} | {max(lf) if lf else 'n/a'} |")
    nc = [r for r in rows if r["status"] != "completed"]
    print("not completed:", [(r["study"].split("_plane_")[1], r["N"], r["Co"], r["status"]) for r in nc])
    print("P1 ok cases:", [(r["study"].split("_plane_")[1], r["N"], r["Co"]) for r in rows if r["verdict"].startswith("P1 ok")])
