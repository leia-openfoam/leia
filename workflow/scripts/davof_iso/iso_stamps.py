"""The stamps of the timed verification: per timed cell (iso/time/<flow>_[<class>_]N<N>_co<CO>)
the fork commit and build id of the binary its runs used (the log's Exec line -> the
worktree; the sha256 of advectorVoF and libVoF.so there), the leia commit that last touched
the case template, the study configuration and the runner, and the sha256 over the rendered
dictionaries of the prepared cases of the study (system/*Dict, fvSolution, fvSchemes,
controlDict, constant/transportProperties, 0/U, 0/alpha.water headers excluded? no: whole
files). Writes iso/csv/iso_stamps.csv. Claude, 2026-11-07."""
import csv, glob, hashlib, os, re, subprocess
W = "/work/scratch/tm83tomy/davof/agents/swept"
I = W + "/iso"; S = I + "/leia/studies/davof"
def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()[:16]
def git(args, cwd):
    return subprocess.run(["git"] + args, cwd=cwd, capture_output=True, text=True).stdout.strip()
build_cache = {}
def build_id(wt):
    if wt not in build_cache:
        b = glob.glob(wt + "/platforms/*/bin/advectorVoF"); l = glob.glob(wt + "/platforms/*/lib/libVoF.so")
        build_cache[wt] = (git(["rev-parse", "--short", "HEAD"], wt), sha(b[0]) if b else "-", sha(l[0]) if l else "-")
    return build_cache[wt]
dict_cache = {}
def dict_hash(study):
    if study not in dict_cache:
        hs = []
        for c in sorted(glob.glob(S + "/" + study + "/*_0*")):
            if not os.path.isdir(c) or "." in os.path.basename(c): continue
            for f in ("system/controlDict", "system/fvSolution", "system/fvSchemes", "system/blockMeshDict", "system/initAlphaFieldDict", "system/meshDict", "system/decomposeParDict", "constant/transportProperties"):
                if os.path.exists(c + "/" + f): hs.append(os.path.basename(c) + "/" + f + ":" + sha(c + "/" + f))
        dict_cache[study] = (hashlib.sha256("\n".join(hs).encode()).hexdigest()[:16], len(hs))
    return dict_cache[study]
rows = []
for d in sorted(glob.glob(I + "/time/*_N*_co*")):
    b = os.path.basename(d)
    m = re.match(r"(\w+?)_(hex|perturbed|poly|tet)_N(\d+)_co(\d)$", b) or re.match(r"(\w+)_N(\d+)_co(\d)$", b)
    if not m: continue
    if len(m.groups()) == 4: flow, mc, N, co = m.groups()
    else: flow, N, co = m.groups(); mc = "hex"
    case = "iso%s3D" % flow if not flow.endswith("2D") else "iso" + flow
    study = "%s_%sN%s" % (case, mc, N)
    for run in sorted(glob.glob(d + "/*_r*")):
        log = run + "/log.advectorVoF"
        if not os.path.exists(log): continue
        wt = "-"
        with open(log, errors="replace") as f:
            for line in f:
                if line.startswith("Exec"):
                    mm = re.search(r"(/work/\S+?/TwoPhaseFlow[^/]*)/platforms", line); wt = mm.group(1) if mm else line.split(":", 1)[1].strip(); break
                if line.startswith("Create time"): break
        fc, fb, fl = build_id(wt) if os.path.isdir(wt) else ("-", "-", "-")
        lc = git(["log", "-1", "--format=%h", "--", "cases/davof/" + case, "cases/davof/" + case + ".parameter", "cases/davof/" + case + "_poly.parameter", "config/davof/" + study + ".yaml", "workflow/scripts/davof_iso/iso-time.sbatch", "workflow/scripts/tpf_rows.py", "workflow/scripts/materialize.py", "workflow/scripts/foam_param.py"], I + "/leia")
        dh, nf = dict_hash(study)
        method, rep = os.path.basename(run).rsplit("_r", 1)
        rows.append(dict(flow=flow, mesh=mc, N=N, co_index=co, method=method, rep=rep, fork_worktree=os.path.basename(wt), fork_commit=fc, advectorVoF_sha256_16=fb, libVoF_sha256_16=fl, leia_commit_inputs=lc, leia_head=git(["rev-parse", "--short", "HEAD"], I + "/leia"), study=study, dict_hash_16=dh, dict_files=nf))
os.makedirs(I + "/csv", exist_ok=True)
with open(I + "/csv/iso_stamps.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("wrote", len(rows), "rows")
seen = set()
for r in rows:
    k = (r["flow"], r["mesh"], r["N"], r["fork_commit"], r["advectorVoF_sha256_16"], r["leia_commit_inputs"], r["dict_hash_16"])
    if k in seen: continue
    seen.add(k); print(k, r["study"], r["dict_files"])
