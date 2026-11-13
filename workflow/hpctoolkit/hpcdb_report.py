#!/usr/bin/env python3
"""Text report of an HPCToolkit database (meta.db and profile.db, format version 4.0).

HPCToolkit 2025 has no text tool for its database (hpcviewer is a GUI). This script reads
the two files as share/doc/hpctoolkit/FORMATS.md of the installation documents them, and
prints, for the time metric:

  1. the profiles (one per thread) and their total time,
  2. the functions by exclusive time and by inclusive time,
  3. the calling-context tree, reduced to functions, cut at a share of the total,
  4. the time below the MPI library, by MPI function and by its callers.

Exclusive time of a function: the samples whose innermost function context is that
function. Inclusive time: the HPCToolkit `execution` value of the outermost context of the
function on each call path, so recursion is not counted twice. A function that the compiler
inlined has no context of its own in a build without -g; its time belongs to its caller.

    hpcdb_report.py <database dir> [--metric REALTIME] [--tree-cut 0.5] [--top 40]
                    [--focus <substring of a function name>] [--depth 60]
"""
import argparse
import os
import struct
import sys
from collections import defaultdict

import numpy as np


def u(buf, off, fmt):
    return struct.unpack_from("<" + fmt, buf, off)[0]


def cstr(buf, p):
    if p == 0:
        return ""
    e = buf.index(b"\0", p)
    return buf[p:e].decode("utf-8", "replace")


class Meta:
    def __init__(self, path):
        b = self.b = open(path, "rb").read()
        assert b[:10] == b"HPCTOOLKIT" and b[10:14] == b"meta", "not a meta.db"
        sec = [(u(b, 0x10 + 16 * i, "Q"), u(b, 0x18 + 16 * i, "Q")) for i in range(8)]
        (_, pGen), (_, pId), (_, pMet), (_, pCtx), _, (_, pMod), (_, pFil), (_, pFun) = sec
        self.title = cstr(b, u(b, pGen, "Q"))
        # identifier kinds
        pp, nk = u(b, pId, "Q"), u(b, pId + 8, "B")
        self.kinds = [cstr(b, u(b, pp + 8 * i, "Q")) for i in range(nk)]
        # metrics: (metric name, scope name) -> propMetricId
        pM, nM, szM, szSI = u(b, pMet, "Q"), u(b, pMet + 8, "I"), u(b, pMet + 12, "B"), u(b, pMet + 13, "B")
        self.metrics = {}
        for i in range(nM):
            m = pM + i * szM
            name = cstr(b, u(b, m, "Q"))
            pSI, nSI = u(b, m + 8, "Q"), u(b, m + 0x18, "H")
            for j in range(nSI):
                s = pSI + j * szSI
                scope = cstr(b, u(b, u(b, s, "Q"), "Q"))
                self.metrics[(name, scope)] = u(b, s + 8, "H")
        # context tree: arrays indexed by ctxId
        self.parent, self.kind, self.name, self.module, self.rel = {}, {}, {}, {}, {}
        pE, nE, szE = u(b, pCtx, "Q"), u(b, pCtx + 8, "H"), u(b, pCtx + 10, "B")
        stack = []
        for i in range(nE):
            e = pE + i * szE
            cid = u(b, e + 0x10, "I")
            self.parent[cid] = 0
            self.kind[cid] = "entry"
            self.name[cid] = "<" + cstr(b, u(b, e + 0x18, "Q")) + ">"
            self.module[cid] = ""
            self.rel[cid] = 0
            stack.append((u(b, e + 8, "Q"), u(b, e, "Q"), cid))
        lex = {0: "function", 1: "loop", 2: "line", 3: "instruction"}
        while stack:
            p, sz, par = stack.pop()
            end = p + sz
            while p and p < end:
                szc, pc, cid = u(b, p, "Q"), u(b, p + 8, "Q"), u(b, p + 0x10, "I")
                flags, rel, lt, nfw = struct.unpack_from("<BBBB", b, p + 0x14)
                f = p + 0x20
                name, mod = "", ""
                if flags & 1:
                    pf = u(b, f, "Q")
                    f += 8
                    name = cstr(b, u(b, pf, "Q"))
                    pm = u(b, pf + 8, "Q")
                    mod = cstr(b, u(b, pm + 8, "Q")) if pm else ""
                    if not name:
                        name = "<unknown function at 0x%x>" % u(b, pf + 0x10, "Q")
                if flags & 2:
                    f += 16
                if flags & 4:
                    pm = u(b, f, "Q")
                    if pm and not mod:
                        mod = cstr(b, u(b, pm + 8, "Q"))
                self.parent[cid] = par
                self.kind[cid] = lex.get(lt, "?")
                self.name[cid] = name if lt == 0 else ""
                if lt == 0 and not name:
                    self.name[cid] = "<unknown function>"
                self.module[cid] = os.path.basename(mod)
                self.rel[cid] = rel
                if pc:
                    stack.append((pc, szc, cid))
                p += 0x20 + 8 * nfw
        # nearest function (or entry) context at or above each context
        self.fn = {}
        for c in self.parent:
            x = c
            while self.kind[x] not in ("function", "entry"):
                x = self.parent[x]
            self.fn[c] = x

    def key(self, c):
        return (self.name[c], self.module[c])

    def fn_parent(self, c):
        """The function context that encloses the function context c (0 at the top)."""
        p = self.parent[c]
        return self.fn[p] if p else 0


def read_profiles(path, meta):
    b = open(path, "rb").read()
    assert b[:10] == b"HPCTOOLKIT" and b[10:14] == b"prof", "not a profile.db"
    pPI = u(b, 0x18, "Q")
    pP, nP, szP = u(b, pPI, "Q"), u(b, pPI + 8, "I"), u(b, pPI + 12, "B")
    vdt = np.dtype([("m", "<u2"), ("v", "<f8")])
    idt = np.dtype([("c", "<u4"), ("s", "<u8")])
    out = []
    for i in range(nP):
        p = pP + i * szP
        nV, pV, nC, pC = u(b, p, "Q"), u(b, p + 8, "Q"), u(b, p + 0x10, "I"), u(b, p + 0x18, "Q")
        pT, fl = u(b, p + 0x20, "Q"), u(b, p + 0x28, "I")
        if fl & 1 or pT == 0:
            continue
        nI = u(b, pT, "H")
        ident = {}
        for j in range(nI):
            q = pT + 8 + 16 * j
            k = u(b, q, "B")
            ident[meta.kinds[k] if k < len(meta.kinds) else str(k)] = (u(b, q + 4, "I"), u(b, q + 8, "Q"))
        vals = np.frombuffer(b, dtype=vdt, count=nV, offset=pV)
        idx = np.frombuffer(b, dtype=idt, count=nC, offset=pC)
        out.append((ident, vals, idx))
    return out


def per_ctx(vals, idx, mid):
    """ctxId -> value of metric mid."""
    sel = np.nonzero(vals["m"] == mid)[0]
    if sel.size == 0:
        return {}
    starts = idx["s"]
    pos = np.searchsorted(starts, sel, side="right") - 1
    return dict(zip(idx["c"][pos].tolist(), vals["v"][sel].tolist()))


MPI_MODS = ("libmpi", "libopen-pal", "libopen-rte", "mca_", "libpmi", "libucx", "libucp", "libuct",
            "libucs", "libevent", "libhwloc", "libfabric", "libpsm")


def is_mpi(mod):
    return mod.startswith(MPI_MODS)


def short(name, n):
    return name if len(name) <= n else name[: n - 3] + "..."


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("db")
    ap.add_argument("--metric", default="REALTIME")
    ap.add_argument("--tree-cut", type=float, default=0.5, help="percent of the total")
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--depth", type=int, default=60)
    ap.add_argument("--width", type=int, default=150)
    ap.add_argument("--focus", default="", help="print the tree only below functions with this substring")
    a = ap.parse_args()

    meta = Meta(os.path.join(a.db, "meta.db"))
    names = sorted({m for m, _ in meta.metrics})
    mname = next((m for m in names if m.startswith(a.metric)), None)
    if mname is None:
        sys.exit("metric %s not in %s" % (a.metric, names))
    m_point, m_exec = meta.metrics[(mname, "point")], meta.metrics[(mname, "execution")]
    profs = read_profiles(os.path.join(a.db, "profile.db"), meta)

    def label(ident):
        return " ".join("%s=%d" % (k, v[0]) for k, v in ident.items())

    def rank_of(ident):
        return ident.get("RANK", (0, 0))[0]

    print("database: %s" % os.path.abspath(a.db))
    print("metric: %s   metrics in the database: %s" % (mname, ", ".join(names)))
    print("\n== 1. Profiles (one per thread) ==")
    data = []
    for ident, vals, idx in sorted(profs, key=lambda p: (rank_of(p[0]), p[0].get("THREAD", (0, 0))[0])):
        pt, ex = per_ctx(vals, idx, m_point), per_ctx(vals, idx, m_exec)
        tot = sum(pt.values())
        main_thread = ident.get("THREAD", (0, 0))[0] == 0
        print("  %-40s total %10.2f s  %s" % (label(ident), tot, "" if main_thread else "(not a main thread, excluded)"))
        if main_thread:
            data.append((rank_of(ident), pt, ex, tot))
    nR = len(data)
    T = sum(d[3] for d in data)
    print("  main threads: %d   sum of their time: %.2f s   mean: %.2f s" % (nR, T, T / max(nR, 1)))
    if nR > 1:
        tt = [d[3] for d in data]
        print("  per rank: " + "  ".join("r%d %.2f" % (d[0], d[3]) for d in data))
        print("  max/mean = %.3f" % (max(tt) / (sum(tt) / nR)))

    # ---- flat views -------------------------------------------------------------------
    excl = defaultdict(lambda: [0.0] * nR)
    incl = defaultdict(lambda: [0.0] * nR)
    for r, (_, pt, ex, _) in enumerate(data):
        for c, v in pt.items():
            excl[meta.key(meta.fn[c])][r] += v
        for c, v in ex.items():
            if meta.kind.get(c) != "function":
                continue
            k = meta.key(c)
            x = meta.fn_parent(c)
            rec = False
            while x:
                if meta.key(x) == k:
                    rec = True
                    break
                x = meta.fn_parent(x)
            if not rec:
                incl[k][r] += v

    def cols(v):
        if nR == 1:
            return ""
        return "  [" + " ".join("%7.2f" % x for x in v) + "]"

    def flat(title, d):
        print("\n== %s ==" % title)
        print("  %9s %6s  %s%s" % ("time [s]", "%", "function  {module}", "   [per rank, s]" if nR > 1 else ""))
        for k, v in sorted(d.items(), key=lambda kv: -sum(kv[1]))[: a.top]:
            s = sum(v)
            print("  %9.2f %6.2f  %s  {%s}%s" % (s, 100 * s / T, short(k[0], a.width), k[1], cols(v)))

    flat("2a. Functions by EXCLUSIVE time (sum over the main threads)", excl)
    flat("2b. Functions by INCLUSIVE time (sum over the main threads)", incl)

    # exclusive time by load module
    bymod = defaultdict(lambda: [0.0] * nR)
    for k, v in excl.items():
        for r in range(nR):
            bymod[k[1]][r] += v[r]
    print("\n== 2c. EXCLUSIVE time by load module ==")
    for k, v in sorted(bymod.items(), key=lambda kv: -sum(kv[1]))[:25]:
        print("  %9.2f %6.2f  %s%s" % (sum(v), 100 * sum(v) / T, k or "<none>", cols(v)))

    # ---- merged function tree ---------------------------------------------------------
    # node = tuple of function keys from the top; siblings with the same key are merged
    path_cache = {0: ()}

    def path(c):
        """Path of function keys of the function context c."""
        if c in path_cache:
            return path_cache[c]
        p = path(meta.fn_parent(c)) + (meta.key(c),)
        path_cache[c] = p
        return p

    t_incl = defaultdict(lambda: [0.0] * nR)
    t_excl = defaultdict(lambda: [0.0] * nR)
    for r, (_, pt, ex, _) in enumerate(data):
        for c, v in ex.items():
            if meta.kind.get(c) in ("function", "entry"):
                t_incl[path(c)][r] += v
        for c, v in pt.items():
            t_excl[path(meta.fn[c])][r] += v
    children = defaultdict(list)
    for p in t_incl:
        if p:
            children[p[:-1]].append(p)

    cut = a.tree_cut / 100 * T

    def show(p, depth):
        s = sum(t_incl[p])
        if s < cut or depth > a.depth:
            return
        e = sum(t_excl.get(p, [0.0]))
        print("  %8.2f %6.2f %8.2f  %s%s  {%s}%s" % (
            s, 100 * s / T, e, "| " * depth, short(p[-1][0], a.width), p[-1][1], cols(t_incl[p])))
        for ch in sorted(children.get(p, []), key=lambda q: -sum(t_incl[q])):
            show(ch, depth + 1)

    print("\n== 3. Calling-context tree, functions only, cut at %.2f %% of the total ==" % a.tree_cut)
    print("  %8s %6s %8s  %s" % ("incl [s]", "%", "excl [s]", "function  {module}"))
    if a.focus:
        seen = set()
        for p in sorted(t_incl, key=len):
            if a.focus in p[-1][0] and not any(p[:i] in seen for i in range(1, len(p))):
                seen.add(p)
        agg = defaultdict(float)
        for p in seen:
            agg[p[-1]] += sum(t_incl[p])
        for p in sorted(seen, key=lambda q: -sum(t_incl[q])):
            print("  -- called from: %s" % " <- ".join(short(k[0], 60) for k in reversed(p[-4:-1])))
            show(p, 0)
    else:
        for p in sorted(children.get((), []), key=lambda q: -sum(t_incl[q])):
            show(p, 0)

    # ---- MPI ---------------------------------------------------------------------------
    print("\n== 4. Time below the MPI library (outermost MPI function on each call path) ==")
    mpi_tot = [0.0] * nR
    by_fn = defaultdict(lambda: [0.0] * nR)
    by_site = defaultdict(lambda: [0.0] * nR)
    by_leia = defaultdict(lambda: [0.0] * nR)
    for p, v in t_incl.items():
        if not p or not is_mpi(p[-1][1]) or any(is_mpi(k[1]) for k in p[:-1]):
            continue
        for r in range(nR):
            mpi_tot[r] += v[r]
            by_fn[p[-1][0]][r] += v[r]
        callers = [k for k in reversed(p[:-1])]
        site = " <- ".join(short(k[0], 70) for k in callers[:3])
        leia = next((k[0] for k in callers if k[1].startswith("libleia") or k[1].startswith("leia")), "<none>")
        for r in range(nR):
            by_site[(p[-1][0], site)][r] += v[r]
            by_leia[leia][r] += v[r]
    M = sum(mpi_tot)
    print("  total: %.2f s = %.2f %% of the time of the main threads%s" % (M, 100 * M / T, cols(mpi_tot)))
    if nR > 1:
        print("  share per rank [%]: " + " ".join("%6.2f" % (100 * mpi_tot[r] / data[r][3]) for r in range(nR)))
    print("  -- by MPI function")
    for k, v in sorted(by_fn.items(), key=lambda kv: -sum(kv[1]))[:20]:
        print("  %9.2f %6.2f  %s%s" % (sum(v), 100 * sum(v) / T, k, cols(v)))
    print("  -- by the first leia function above the MPI call")
    for k, v in sorted(by_leia.items(), key=lambda kv: -sum(kv[1]))[:30]:
        print("  %9.2f %6.2f  %s%s" % (sum(v), 100 * sum(v) / T, short(k, a.width), cols(v)))
    print("  -- by MPI function and its three nearest callers")
    for k, v in sorted(by_site.items(), key=lambda kv: -sum(kv[1]))[:40]:
        print("  %9.2f %6.2f  %s <- %s%s" % (sum(v), 100 * sum(v) / T, k[0], k[1], cols(v)))


if __name__ == "__main__":
    main()
