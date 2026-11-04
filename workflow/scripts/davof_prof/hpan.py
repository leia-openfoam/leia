"""hpan.py -- the analysis of one HPCToolkit database (Claude, 2026-11-04):
    python3 hpan.py --db DB --out PREFIX [--steps N] [--rules rules.py] [--src REGEX] [--tree-min 0.3]
writes PREFIX.summary.txt, PREFIX.flat.csv, PREFIX.lines.csv, PREFIX.tree.txt,
PREFIX.ranks.csv and, with --rules, PREFIX.cat.csv and PREFIX.catranks.csv.
Exclusive time of a function = the samples in its own lexical scope, not in a
called or an inlined function (the latter counted for the inlined function).
"""
import argparse, os, re, sys, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hpdb

MPI_MOD = re.compile(r"lib(mpi|open-pal|open-rte|pmix|ucx|ucp|uct|ucs|ucm|hwloc|event|psm|fabric|ibverbs|rdmacm|mlx|hcoll|pthread)")
MPI_FN = re.compile(r"^(MPI_|PMPI_|ompi_|opal_|mca_|pmix|ucp_|uct_|ucs_|mlx5_|ibv_|__GI___poll|poll$|epoll)")
PSTREAM_FN = re.compile(r"Foam::(UPstream|UIPstream|UOPstream|Pstream|IPstream|OPstream|PstreamBuffers|processorFvPatchField|processorPolyPatch|processorFvPatch|processorLduInterface)")


def short(name, n=110):
    name = re.sub(r"\(anonymous namespace\)::", "", name)
    return name if len(name) <= n else name[:n - 3] + "..."


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--steps", type=float, default=1.0)
    ap.add_argument("--rules", default=None)
    ap.add_argument("--src", default=r"TwoPhaseFlow")
    ap.add_argument("--tree-min", type=float, default=0.3)
    ap.add_argument("--tree-depth", type=int, default=16)
    a = ap.parse_args()
    db = hpdb.Database(a.db)
    profs = db.profiles()
    met = {n: d for n, d in db.metrics.items()}
    cpu = [n for n in met if n.startswith("CPUTIME")][0]
    llc = [n for n in met if "LLC" in n]
    llc = llc[0] if llc else None
    N = db.nCtx
    for p in profs:
        if p.nCtxs: N = max(N, int(p._load()[1]["c"].max()) + 1)
    # ---- static context attributes
    ids = np.arange(N)
    parent = np.full(N, -1, dtype=np.int64)
    isfn = np.zeros(N, dtype=bool)          # a function node: lexicalType 0 (entry, call or inlined call)
    fnnode = np.full(N, -1, dtype=np.int64)  # the nearest function node at or above (own lexical scope)
    fname = {}
    srcloc = np.full(N, -1, dtype=np.int64)  # ctx id of the nearest context with a source line at or above
    depth = np.zeros(N, dtype=np.int32)
    order = []
    stack = list(db.roots)
    while stack:
        c = stack.pop()
        order.append(c.id)
        par = c.parent.id if c.parent is not None else -1
        parent[c.id] = par
        depth[c.id] = depth[par] + 1 if par >= 0 else 0
        if c.lexicalType == 0 or c.relation != 0 or par < 0:
            # a function node: an entry, a call target or an inlined call target
            isfn[c.id] = True; fnnode[c.id] = c.id
            if c.function: fname[c.id] = c.function.name
            elif c.module: fname[c.id] = "<unknown procedure> 0x%x [%s]" % (c.offset, os.path.basename(c.module))
            elif c.file: fname[c.id] = "<unknown procedure> %s:%d" % (os.path.basename(c.file), c.line)
            else: fname[c.id] = "<unknown procedure>"
        else:
            fnnode[c.id] = fnnode[par] if par >= 0 else -1
        if c.file is not None: srcloc[c.id] = c.id
        else: srcloc[c.id] = srcloc[par] if par >= 0 else -1
        stack.extend(c.children)
    # outermost function nodes per function name (recursion-safe inclusive)
    outer = np.zeros(N, dtype=bool)
    for cid in order:
        if not isfn[cid]: continue
        nm = fname[cid]; p = parent[cid]; seen = False
        while p >= 0:
            if isfn[p] and fname[p] == nm: seen = True; break
            p = parent[p]
        outer[cid] = not seen
    # module of a function node
    def modof(cid):
        c = db.ctx.get(cid)
        if c is None: return ""
        if c.function and c.function.module: return os.path.basename(c.function.module)
        if c.module: return os.path.basename(c.module)
        return ""
    # MPI classification of function nodes: nearest ancestor function node (incl. itself) that is MPI or Pstream
    cls = {}  # fn node -> "mpi" | "pstream" | ""
    for cid in order:
        if not isfn[cid]: continue
        nm = fname[cid]; m = modof(cid)
        if MPI_MOD.search(m) or MPI_FN.search(nm): cls[cid] = "mpi"
        elif PSTREAM_FN.search(nm): cls[cid] = "pstream"
        else:
            p = fnnode[parent[cid]] if parent[cid] >= 0 else -1
            cls[cid] = cls.get(p, "") if p >= 0 else ""
    # category rules
    rules = []
    alloc_rx = None
    if a.rules:
        ns = {}
        exec(open(a.rules).read(), ns)
        rules = ns["RULES"]  # list of (category, kind, pattern[, (lo, hi)]) kind = fn | line | mod
        alloc_rx = re.compile(ns["ALLOC"]) if "ALLOC" in ns else None
    # allocation: a sample in an allocator function (or its callees)
    isalloc = np.zeros(N, dtype=bool)
    ismpi = np.zeros(N, dtype=bool)
    for cid in order:
        p = parent[cid]
        if isfn[cid]:
            isalloc[cid] = bool(alloc_rx and alloc_rx.search(fname[cid])) or (p >= 0 and isalloc[p])
            ismpi[cid] = cls.get(cid) == "mpi" or (p >= 0 and ismpi[p])
        else:
            isalloc[cid] = isalloc[p] if p >= 0 else False
            ismpi[cid] = ismpi[p] if p >= 0 else False
    cat_of_ctx = np.full(N, -1, dtype=np.int64)
    catnames = [r[0] for r in rules]
    catindex = {}
    for i, r in enumerate(rules): catindex.setdefault(r[0], i)
    cat_site = np.full(N, -1, dtype=np.int64)   # the category with the MPI/Pstream rules dropped (the call site of MPI time)
    if rules:
        excl = set(ns.get("MPI_SITE_EXCLUDE", ()))
        comp = []
        for cat, kind, pat, *rest in rules:
            comp.append((catindex[cat], kind, re.compile(pat), rest[0] if rest else None, cat in excl))
        def match(c, site):
            for ci, kind, rx, rng, isx in comp:
                if site and isx: continue
                if kind == "fn" and c.lexicalType == 0 and c.function and rx.search(c.function.name): return ci
                if kind == "line" and c.file is not None and rx.search(c.file) and rng[0] <= c.line <= rng[1]: return ci
                if kind == "mod" and c.lexicalType == 0 and c.function and c.function.module and rx.search(os.path.basename(c.function.module)): return ci
            return -1
        for cid in order:  # parents before children
            c = db.ctx[cid]
            m = match(c, False)
            cat_of_ctx[cid] = m if m >= 0 else (cat_of_ctx[parent[cid]] if parent[cid] >= 0 else -1)
            m = match(c, True)
            cat_site[cid] = m if m >= 0 else (cat_site[parent[cid]] if parent[cid] >= 0 else -1)
    # ---- per profile values
    mids = met[cpu]
    exec_id, fs_id = mids["execution"], mids["function"]
    lmids = met[llc] if llc else None
    nP = len(profs)
    # own samples of a context = its function-scope value (the lexical subtree,
    # not crossing a call or an inlined call) minus its lexical children's:
    # hpcprof elides the instruction contexts from meta.db, their samples sit here
    lexmask = np.zeros(N, dtype=bool)
    for cid in order:
        if parent[cid] >= 0 and db.ctx[cid].relation == 0: lexmask[cid] = True
    lexpar = parent[lexmask]
    intree = np.zeros(N, dtype=bool); intree[list(db.ctx)] = True
    def own_of(fs):
        own = fs.copy()
        np.subtract.at(own, lexpar, fs[lexmask])
        own[~intree] = 0.0   # the elided leaves: already in their parents' function scope
        own[0] = 0.0
        return own
    tot = np.zeros(nP); pt_sum = np.zeros(N); ex_sum = np.zeros(N); ll_sum = np.zeros(N); own_check = 0.0
    rank_rows = []
    fn_ex_rank = {}   # per function name: array over profiles of exclusive
    cat_rank = np.zeros((nP, max(1, len(catindex))))
    cat_alloc = np.zeros(max(1, len(catindex))); cat_mpi = np.zeros(max(1, len(catindex)))
    cat_mpisite = np.zeros(max(1, len(catindex)))   # MPI samples by the call site's category
    alloc_total = 0.0
    # callers: exclusive time of a function node credited to its caller function
    caller_agg = collections.defaultdict(float)
    names_sorted = None
    name_id = {}
    for cid in order:
        if isfn[cid]: name_id.setdefault(fname[cid], len(name_id))
    fnid_of_ctx = np.array([name_id[fname[fnnode[i]]] if fnnode[i] >= 0 and fnnode[i] in fname else -1 for i in range(N)], dtype=np.int64)
    nF = len(name_id)
    fn_ex = np.zeros((nP, nF))
    for k, p in enumerate(profs):
        fs = p.dense(fs_id, N); ex = p.dense(exec_id, N)
        pt = own_of(fs)
        pt_sum += pt; ex_sum += ex
        tot[k] = ex[0] if ex[0] > 0 else pt.sum()
        own_check += pt.sum()
        ok = fnid_of_ctx >= 0
        fn_ex[k] = np.bincount(fnid_of_ctx[ok], weights=pt[ok], minlength=nF)
        alloc_total += pt[isalloc].sum()
        if rules:
            okc = cat_of_ctx >= 0
            cat_rank[k, :] = np.bincount(cat_of_ctx[okc], weights=pt[okc], minlength=len(catindex))
            cat_alloc += np.bincount(cat_of_ctx[okc & isalloc], weights=pt[okc & isalloc], minlength=len(catindex))
            cat_mpi += np.bincount(cat_of_ctx[okc & ismpi], weights=pt[okc & ismpi], minlength=len(catindex))
            oks = (cat_site >= 0) & ismpi
            cat_mpisite += np.bincount(cat_site[oks], weights=pt[oks], minlength=len(catindex))
        mpi = float(pt[ismpi].sum()) if nP > 1 else 0.0
        pst = sum(pt[i] for i in range(N) if fnnode[i] >= 0 and cls.get(fnnode[i]) == "pstream") if nP > 1 else 0.0
        rank_rows.append((p.rank(), p.label(), tot[k], mpi, pst))
        if llc:
            ll_sum += own_of(p.dense(lmids["function"], N))
    T = tot.sum()
    S = a.steps
    # ---- flat profile by function
    inv = {v: k for k, v in name_id.items()}
    ex_by_name = fn_ex.sum(axis=0)
    inc_by_name = np.zeros(nF); ll_by_name = np.zeros(nF)
    fn_loc = {}
    for cid in order:
        if not isfn[cid]: continue
        i = name_id[fname[cid]]
        if outer[cid]: inc_by_name[i] += ex_sum[cid]
        c = db.ctx[cid]
        if i not in fn_loc and c.function: fn_loc[i] = (modof(cid), c.function.file or "", c.function.line)
    okf = fnid_of_ctx >= 0
    ll_by_name = np.bincount(fnid_of_ctx[okf], weights=ll_sum[okf], minlength=nF)
    with open(a.out + ".flat.csv", "w") as f:
        f.write("function,module,file,line,excl_s,excl_pct,incl_s,incl_pct,excl_s_per_step,llc_miss_samples,excl_max_rank_s,excl_mean_rank_s\n")
        for i in np.argsort(-ex_by_name):
            if ex_by_name[i] <= 0 and inc_by_name[i] <= 0: continue
            m, fl, ln = fn_loc.get(i, ("", "", 0))
            f.write('"%s",%s,%s,%d,%.3f,%.2f,%.3f,%.2f,%.4f,%.0f,%.3f,%.3f\n' % (inv[i].replace('"', "'"), m, fl, ln, ex_by_name[i], 100*ex_by_name[i]/T, inc_by_name[i], 100*inc_by_name[i]/T, ex_by_name[i]/S/nP if nP > 1 else ex_by_name[i]/S, ll_by_name[i], fn_ex[:, i].max(), fn_ex[:, i].mean()))
    # ---- per line (files matching --src), exclusive = point samples whose nearest source line is this one
    srx = re.compile(a.src)
    line_agg = collections.defaultdict(float); line_fn = {}
    for i in range(N):
        s = srcloc[i]
        if s < 0 or pt_sum[i] == 0: continue
        c = db.ctx[s]
        if not srx.search(c.file or ""): continue
        key = (c.file, c.line)
        line_agg[key] += pt_sum[i]
        line_fn.setdefault(key, fname.get(fnnode[i], ""))
    with open(a.out + ".lines.csv", "w") as f:
        f.write("file,line,function,excl_s,excl_pct\n")
        for (fl, ln), v in sorted(line_agg.items(), key=lambda x: -x[1]):
            f.write('%s,%d,"%s",%.3f,%.3f\n' % (fl, ln, line_fn[(fl, ln)].replace('"', "'"), v, 100*v/T))
    # ---- the call tree (function nodes, children merged by name)
    kids = collections.defaultdict(list)
    for cid in order:
        if isfn[cid] and parent[cid] >= 0:
            p = fnnode[parent[cid]]
            if p >= 0: kids[p].append(cid)
    lines = []
    def walk(nodes, d):
        inc = sum(ex_sum[n] for n in nodes)
        if 100*inc/T < a.tree_min or d > a.tree_depth: return
        ex = sum(fn_ex_node[n] for n in nodes)
        lines.append("%s%6.2f%% %9.3f s  (excl %5.2f%%)  %s" % ("  "*d, 100*inc/T, inc, 100*ex/T, short(fname[nodes[0]])))
        groups = collections.defaultdict(list)
        for n in nodes:
            for k in kids[n]: groups[fname[k]].append(k)
        for nm, g in sorted(groups.items(), key=lambda x: -sum(ex_sum[n] for n in x[1])):
            walk(g, d + 1)
    fn_ex_node = np.bincount(fnnode[fnnode >= 0], weights=pt_sum[fnnode >= 0], minlength=N)
    for r in db.roots: walk([r.id], 0)
    with open(a.out + ".tree.txt", "w") as f:
        f.write("inclusive %% of %.3f s (sum over %d profiles), nodes >= %.2f%%\n" % (T, nP, a.tree_min))
        f.write("\n".join(lines) + "\n")
    # ---- ranks (the profiles of one MPI rank, its threads, summed)
    byrank = collections.OrderedDict()
    for k, (r, lab, t, mpi, pst) in enumerate(rank_rows):
        e = byrank.setdefault(r, [0.0, 0.0, 0.0, np.zeros(max(1, len(catindex))), 0])
        e[0] += t; e[1] += mpi; e[2] += pst; e[3] += cat_rank[k]; e[4] += 1
    rk = sorted(byrank)
    r_tot = np.array([byrank[r][0] for r in rk]); r_mpi = np.array([byrank[r][1] for r in rk]); r_pst = np.array([byrank[r][2] for r in rk])
    r_comp = r_tot - r_mpi - r_pst
    nR = len(rk)
    with open(a.out + ".ranks.csv", "w") as f:
        f.write("rank,threads,cputime_s,mpi_s,pstream_s,compute_s,compute_s_per_step" + "".join(',"%s"' % c for c in catindex) + "\n")
        for j, r in enumerate(rk):
            f.write("%d,%d,%.3f,%.3f,%.3f,%.3f,%.4f" % (r, byrank[r][4], r_tot[j], r_mpi[j], r_pst[j], r_comp[j], r_comp[j]/S) + "".join(",%.3f" % byrank[r][3][i] for i in range(len(catindex))) + "\n")
    # callers of the hot functions
    hot = set(int(i) for i in np.argsort(-ex_by_name)[:40])
    for cid in order:
        if not isfn[cid]: continue
        i = name_id[fname[cid]]
        if i not in hot: continue
        p = fnnode[parent[cid]] if parent[cid] >= 0 else -1
        caller_agg[(fname[cid], fname[p] if p >= 0 else "<root>")] += fn_ex_node[cid]
    with open(a.out + ".callers.txt", "w") as f:
        cur = None
        for (fn, cl), v in sorted(caller_agg.items(), key=lambda x: (-ex_by_name[name_id[x[0][0]]], -x[1])):
            if v < 0.002*T: continue
            if fn != cur:
                cur = fn; f.write("%s  (excl %.3f s, %.2f%%)\n" % (short(fn, 100), ex_by_name[name_id[fn]], 100*ex_by_name[name_id[fn]]/T))
            f.write("    %8.3f s %5.2f%%  <- %s\n" % (v, 100*v/T, short(cl, 100)))
    # ---- categories
    if rules:
        with open(a.out + ".cat.csv", "w") as f:
            f.write("category,cputime_s,pct,s_per_step_per_rank,max_rank_s,mean_rank_s,alloc_s,mpi_s,mpi_at_site_s,llc_misses\n")
            okc = cat_of_ctx >= 0
            ll_cat = np.bincount(cat_of_ctx[okc], weights=ll_sum[okc], minlength=len(catindex))
            rest = T - cat_rank.sum()
            for c, i in catindex.items():
                v = cat_rank[:, i].sum()
                rv = np.array([byrank[r][3][i] for r in rk])
                f.write('"%s",%.3f,%.2f,%.4f,%.3f,%.3f,%.3f,%.3f,%.3f,%.0f\n' % (c, v, 100*v/T, v/S/nR, rv.max(), rv.mean(), cat_alloc[i], cat_mpi[i], cat_mpisite[i], ll_cat[i]))
            f.write("(uncategorised),%.3f,%.2f,%.4f,,,,,,\n" % (rest, 100*rest/T, rest/S/nR))
            f.write("(total),%.3f,100,%.4f,%.3f,%.3f,%.3f,%.3f,,%.0f\n" % (T, T/S/nR, r_tot.max(), r_tot.mean(), alloc_total, r_mpi.sum(), ll_sum.sum()))
        with open(a.out + ".catranks.csv", "w") as f:
            f.write("rank," + ",".join('"%s"' % c for c in catindex) + "\n")
            for k, (r, lab, t, mpi, pst) in enumerate(rank_rows):
                f.write("%d," % r + ",".join("%.3f" % cat_rank[k, i] for i in range(len(catindex))) + "\n")
        # the uncategorised functions, for completing the rules
        unc = collections.defaultdict(float)
        for i in range(N):
            if cat_of_ctx[i] < 0 and fnnode[i] >= 0: unc[fname[fnnode[i]]] += pt_sum[i]
        with open(a.out + ".uncat.txt", "w") as f:
            for n, v in sorted(unc.items(), key=lambda x: -x[1])[:60]: f.write("%8.3f %5.2f%% %s\n" % (v, 100*v/T, short(n, 120)))
    # ---- summary
    with open(a.out + ".summary.txt", "w") as f:
        f.write("database %s\nmetrics %s\nprofiles %d (ranks)\ncontexts %d, function nodes %d, functions %d\n" % (a.db, list(met), nP, N, int(isfn.sum()), nF))
        f.write("CPUTIME total %.3f s (own-sample check %.3f s) over %d ranks (%d profiles incl. helper threads): per rank mean %.3f max %.3f min %.3f; steps %g -> %.4f s/step/rank mean, %.4f max\n" % (T, own_check, nR, nP, r_tot.mean(), r_tot.max(), r_tot.min(), S, r_tot.mean()/S, r_tot.max()/S))
        if nR > 1:
            f.write("MPI (libmpi & co) per rank: mean %.3f max %.3f min %.3f s = %.1f%% of CPUTIME; Pstream wrappers mean %.3f s; compute per rank mean %.3f max %.3f min %.3f -> imbalance max/mean %.3f; ranks with compute < 10%% of max: %d\n" % (r_mpi.mean(), r_mpi.max(), r_mpi.min(), 100*r_mpi.sum()/T, r_pst.mean(), r_comp.mean(), r_comp.max(), r_comp.min(), r_comp.max()/r_comp.mean(), int((r_comp < 0.1*r_comp.max()).sum())))
            f.write("per-rank compute s: %s\n" % " ".join("%.1f" % x for x in r_comp))
        f.write("top 25 exclusive:\n")
        for i in np.argsort(-ex_by_name)[:25]:
            m, fl, ln = fn_loc.get(i, ("", "", 0))
            f.write("  %6.2f%% %9.3f s  %s  [%s %s:%d]\n" % (100*ex_by_name[i]/T, ex_by_name[i], short(inv[i], 90), m, os.path.basename(fl), ln))
    print(open(a.out + ".summary.txt").read())


if __name__ == "__main__":
    main()
