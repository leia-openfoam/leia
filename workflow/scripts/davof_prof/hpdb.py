"""hpdb.py -- a reader for HPCToolkit's database format v4 (meta.db and
profile.db, as written by hpcprof 2025.0.1; the layout follows the installed
share/doc/hpctoolkit/FORMATS.md). Claude, 2026-11-04.

    db = Database("hpctoolkit-database")
    db.metrics            # name -> {scope name -> propMetricId}
    db.ctx[ctxId]         # Ctx: parent, relation, lexicalType, function, file, line, module, offset
    for prof in db.profiles(): prof.ids, prof.values(metricId) -> dict ctxId -> value
"""
import struct, os, mmap
import numpy as np


class Ctx:
    __slots__ = ("id", "parent", "relation", "lexicalType", "function", "file", "line", "module", "offset", "propagation", "children")

    def __init__(self):
        self.children = []


class Function:
    __slots__ = ("name", "module", "offset", "file", "line")


class _File:
    def __init__(self, path):
        self.f = open(path, "rb")
        self.m = mmap.mmap(self.f.fileno(), 0, access=mmap.ACCESS_READ)
        magic = self.m[:10]
        assert magic == b"HPCTOOLKIT", (path, magic)
        self.format = self.m[10:14].decode()
        self.major, self.minor = self.m[14], self.m[15]
        assert self.major == 4, self.major

    def u8(self, o): return self.m[o]
    def u16(self, o): return struct.unpack_from("<H", self.m, o)[0]
    def u32(self, o): return struct.unpack_from("<I", self.m, o)[0]
    def u64(self, o): return struct.unpack_from("<Q", self.m, o)[0]
    def cstr(self, o):
        if o == 0: return None
        e = self.m.find(b"\0", o)
        return self.m[o:e].decode(errors="replace")
    def section(self, i):  # i = 0 for the first section at 0x10
        return self.u64(0x10 + 16*i), self.u64(0x18 + 16*i)


class Database:
    def __init__(self, path):
        self.path = path
        self.meta = _File(os.path.join(path, "meta.db"))
        assert self.meta.format == "meta", self.meta.format
        self._read_strings_tables()
        self._read_metrics()
        self._read_contexts()
        self.prof = _File(os.path.join(path, "profile.db"))
        assert self.prof.format == "prof", self.prof.format
        self._read_profiles()

    # ---- meta.db
    def _read_strings_tables(self):
        f = self.meta
        # identifier kind names
        sz, p = f.section(1)
        pp, nk = f.u64(p), f.u8(p + 8)
        self.kindNames = [f.cstr(f.u64(pp + 8*i)) for i in range(nk)]
        # modules
        sz, p = f.section(5)
        pm, n, szm = f.u64(p), f.u32(p + 8), f.u16(p + 12)
        self.modules = {}
        for i in range(n):
            o = pm + i*szm
            self.modules[o] = f.cstr(f.u64(o + 8))
        # files
        sz, p = f.section(6)
        pf, n, szf = f.u64(p), f.u32(p + 8), f.u16(p + 12)
        self.files = {}
        for i in range(n):
            o = pf + i*szf
            self.files[o] = f.cstr(f.u64(o + 8))
        # functions
        sz, p = f.section(7)
        pfn, n, szfn = f.u64(p), f.u32(p + 8), f.u16(p + 12)
        self.functions = {}
        for i in range(n):
            o = pfn + i*szfn
            fn = Function()
            fn.name = f.cstr(f.u64(o))
            fn.module = self.modules.get(f.u64(o + 8))
            fn.offset = f.u64(o + 16)
            fn.file = self.files.get(f.u64(o + 24))
            fn.line = f.u32(o + 32)
            if fn.name is None:
                fn.name = "<unknown %s+0x%x>" % (os.path.basename(fn.module or "?"), fn.offset)
            self.functions[o] = fn

    def _read_metrics(self):
        f = self.meta
        sz, p = f.section(2)
        pM, nM, szM, szPSI, szSS = f.u64(p), f.u32(p + 8), f.u8(p + 12), f.u8(p + 13), f.u8(p + 14)
        pS, nS, szS = f.u64(p + 16), f.u16(p + 24), f.u8(p + 26)
        self.scopes = {}
        for i in range(nS):
            o = pS + i*szS
            self.scopes[o] = (f.cstr(f.u64(o)), f.u8(o + 8), f.u8(o + 9))  # name, type, propagationIndex
        self.metrics = {}        # name -> {scope: propMetricId}
        self.metricStats = {}    # name -> {(scope, formula, combine): statMetricId}
        for i in range(nM):
            o = pM + i*szM
            name = f.cstr(f.u64(o))
            pPSI, pSS, nPSI, nSS = f.u64(o + 8), f.u64(o + 16), f.u16(o + 24), f.u16(o + 26)
            d = {}
            for j in range(nPSI):
                q = pPSI + j*szPSI
                d[self.scopes[f.u64(q)][0]] = f.u16(q + 8)
            self.metrics[name] = d
            s = {}
            for j in range(nSS):
                q = pSS + j*szSS
                s[(self.scopes[f.u64(q)][0], f.cstr(f.u64(q + 8)), f.u8(q + 16))] = f.u16(q + 18)
            self.metricStats[name] = s

    def _read_contexts(self):
        f = self.meta
        sz, p = f.section(3)
        pE, nE, szE = f.u64(p), f.u16(p + 8), f.u8(p + 10)
        self.ctx = {}
        self.roots = []
        stack = []
        for i in range(nE):
            o = pE + i*szE
            c = Ctx()
            c.id = f.u32(o + 16); c.parent = None; c.relation = 0; c.lexicalType = 0
            c.function = None; c.file = None; c.line = 0; c.module = None; c.offset = 0; c.propagation = 0
            name = f.cstr(f.u64(o + 24))
            fn = Function(); fn.name = "<entry: %s>" % name; fn.module = None; fn.offset = 0; fn.file = None; fn.line = 0
            c.function = fn
            self.ctx[c.id] = c
            self.roots.append(c)
            stack.append((c, f.u64(o + 8), f.u64(o)))
        while stack:
            parent, pC, szC = stack.pop()
            o, end = pC, pC + szC
            while o < end:
                c = Ctx()
                szCh, pCh = f.u64(o), f.u64(o + 8)
                c.id = f.u32(o + 16)
                flags, c.relation, c.lexicalType, nFlex = f.u8(o + 20), f.u8(o + 21), f.u8(o + 22), f.u8(o + 23)
                c.propagation = f.u16(o + 24)
                c.parent = parent
                c.function = None; c.file = None; c.line = 0; c.module = None; c.offset = 0
                cur = o + 32
                if flags & 1:
                    c.function = self.functions.get(f.u64(cur)); cur += 8
                if flags & 2:
                    c.file = self.files.get(f.u64(cur)); cur += 8
                    c.line = f.u32(cur); cur += 4
                if flags & 4:
                    cur = (cur + 7) & ~7
                    c.module = self.modules.get(f.u64(cur)); cur += 8
                    c.offset = f.u64(cur); cur += 8
                self.ctx[c.id] = c
                parent.children.append(c)
                if pCh and szCh: stack.append((c, pCh, szCh))
                o += 32 + 8*nFlex
        self.nCtx = max(self.ctx) + 1

    # ---- profile.db
    def _read_profiles(self):
        f = self.prof
        sz, p = f.section(0)
        pP, nP, szP = f.u64(p), f.u32(p + 8), f.u8(p + 12)
        self._profiles = []
        for i in range(nP):
            o = pP + i*szP
            nValues, pValues, nCtxs, pCtxIdx = f.u64(o), f.u64(o + 8), f.u32(o + 16), f.u64(o + 24)
            pTuple, flags = f.u64(o + 32), f.u32(o + 40)
            ids = []
            if pTuple:
                nIds = f.u16(pTuple)
                for j in range(nIds):
                    q = pTuple + 8 + 16*j
                    kind, fl, lid, pid = f.u8(q), f.u16(q + 2), f.u32(q + 4), f.u64(q + 8)
                    ids.append((self.kindNames[kind] if kind < len(self.kindNames) else str(kind), lid, pid, fl & 1))
            self._profiles.append(Profile(self, i, nValues, pValues, nCtxs, pCtxIdx, ids, flags & 1))

    def profiles(self, summary=False):
        return [p for p in self._profiles if p.isSummary == summary]


_valdt = np.dtype({"names": ["m", "v"], "formats": ["<u2", "<f8"], "offsets": [0, 2], "itemsize": 10})
_idxdt = np.dtype({"names": ["c", "s"], "formats": ["<u4", "<u8"], "offsets": [0, 4], "itemsize": 12})


class Profile:
    def __init__(self, db, index, nValues, pValues, nCtxs, pCtxIdx, ids, isSummary):
        self.db, self.index, self.ids, self.isSummary = db, index, ids, isSummary
        self.nValues, self.pValues, self.nCtxs, self.pCtxIdx = nValues, pValues, nCtxs, pCtxIdx
        self._loaded = None

    def rank(self):
        for k, lid, pid, phys in self.ids:
            if k == "RANK": return lid
        return 0

    def label(self):
        return " ".join("%s %d" % (k, pid if phys else lid) for k, lid, pid, phys in self.ids)

    def _load(self):
        if self._loaded is None:
            m = self.db.prof.m
            vals = np.frombuffer(m, dtype=_valdt, count=self.nValues, offset=self.pValues)
            idx = np.frombuffer(m, dtype=_idxdt, count=self.nCtxs, offset=self.pCtxIdx)
            self._loaded = (vals, idx)
        return self._loaded

    def dense(self, metricId, n=None):
        """values of metricId as a dense array over ctxId (0 = the global context)"""
        vals, idx = self._load()
        n = n or self.db.nCtx
        if self.nCtxs: n = max(n, int(idx["c"].max()) + 1)
        out = np.zeros(n)
        if self.nValues == 0: return out
        starts = idx["s"].astype(np.int64)
        ends = np.append(starts[1:], self.nValues)
        sel = vals["m"] == metricId
        # which ctx does each value belong to: searchsorted over starts
        pos = np.nonzero(sel)[0]
        if pos.size == 0: return out
        which = np.searchsorted(starts, pos, side="right") - 1
        out[idx["c"][which]] = vals["v"][pos]
        return out
