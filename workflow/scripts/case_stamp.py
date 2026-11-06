#!/usr/bin/env python3
"""The case stamp: a study recomputes a case only when something that decides its result
changed (Tomislav, 2026-11-07: "only re-calculate studies if we change something in the
methods"; docs/davof/PORT-DESIGN.md section 5, workflow/README.md "Recompute only on
change").

Every case directory carries three small JSON files:

  .leia_render_stamp.json   the dictionaries' part of the stamp: the hash of the case's
                            rendered dictionaries (system/, constant/ without polyMesh,
                            0.org/, the generated decomposeParDict), of which the mesh
                            inputs (blockMeshDict, meshDict, triSurface) are named
                            separately, plus np and the mode. Written by the Snakefile at
                            parse time from the templates and the variation's tokens,
                            rendered in memory exactly as materialize.py renders them to
                            disk; input of the generate_case rule, so a change re-materialises
                            the case (and re-meshes, re-initialises and re-solves it).
  .leia_stamp_current.json  the whole stamp: the render stamp, the binaries (the solver and
                            the field initialiser, their sha256 and that of every library
                            they load from this clone's platforms/, from the TwoPhaseFlow
                            fork's and from OpenFOAM's, by the ELF DT_NEEDED closure), the
                            OpenFOAM version, the solve command, the post-solve command and
                            the scripts it names, the mesh command. Its id is the sha256 of
                            these trigger fields; the record fields (leia's git describe,
                            the fork's, the tokens, the time, the grazingCut key) are kept
                            for the tables and do not enter the id. Written at parse time,
                            only when its content changed (so its mtime moves only then);
                            input of the solve rule, so Snakemake (--rerun-triggers mtime)
                            schedules a solve exactly when the stamp changed.
  .leia_stamp.json          the stamp a completed solve ran under: written at the end of the
                            solve rule from the current stamp, with the binaries and the
                            dictionaries hashed again at run time (a rebuild between parse
                            and run is recorded, not hidden), the exit code and the job id.
                            The solve rule's first step compares it with the current stamp
                            and skips the solver when the ids are equal and a result is
                            there (`check`), so a scheduled solve over an unchanged case
                            (an mtime sweep, a touched file) costs a few seconds and no
                            result. Every previous stamp is appended to
                            .leia_stamp_history.jsonl, and the previous run's log, CSV,
                            launch record and postProcessing are moved to
                            previous/<old id>/ before a recompute (`archive`).

A finished case without a stored stamp (run before 2026-11-07) is ADOPTED at parse time:
the current stamp is stored as its stamp with "adopted": true, and the stamp files are
back-dated to the case's result so that nothing is rescheduled; the table shows the
adoption. aggregate.py carries stampId, stampBinary, stampDicts, stampLeia, stampTpf,
stampGrazingCut, stampAdopted and stampCompletedAt per case.

Command line:
    case_stamp.py check <case> --csv <solver csv>        exit 0: up to date (skip the solver)
    case_stamp.py write <case> --apps <solver> [<app>..] --rc <rc>   the stored stamp
    case_stamp.py archive <case> --solver <solver> --csv <csv>       previous/<old id>/
    case_stamp.py report <study dir> [--verbose]         what would be recomputed, and why
The Snakefile imports refresh_study() for the parse-time part.
"""
import datetime, functools, glob, hashlib, json, os, re, shutil, subprocess, sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
RENDER_FILE = ".leia_render_stamp.json"
CURRENT_FILE = ".leia_stamp_current.json"
STAMP_FILE = ".leia_stamp.json"
HISTORY_FILE = ".leia_stamp_history.jsonl"
VERSION = 1
# the dictionary files of a case: these directories, every file (polyMesh excluded: it is
# the mesh rule's output, not an input)
DICT_DIRS = ("system", "constant", "0.org")
MESH_INPUTS = re.compile(r"^(system/(blockMeshDict|meshDict|snappyHexMeshDict|topoSetDict.*|extrudeMeshDict|refineMeshDict)|constant/triSurface/.*)$")
_TOKEN = re.compile(r"@!([A-Z0-9_]+)!@")


# ----------------------------------------------------------------------------- hashing
def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


@functools.lru_cache(maxsize=None)
def _sha_file_cached(path, mtime_ns, size):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha_file(path):
    st = os.stat(path)
    return _sha_file_cached(os.path.realpath(path), st.st_mtime_ns, st.st_size)


def files_hash(files):
    """files: {relative path: bytes}; the hash of the sorted names and contents."""
    h = hashlib.sha256()
    for rel in sorted(files):
        h.update(rel.encode()); h.update(b"\0"); h.update(str(len(files[rel])).encode()); h.update(b"\0")
        h.update(files[rel]); h.update(b"\0")
    return h.hexdigest()


# ----------------------------------------------------------------------------- the dictionaries
def _decompose_par_dict_text(np_):
    import materialize   # the same text materialize.py writes
    return materialize.decompose_par_dict_text(np_)


def rendered_files(base_case, tokens, np_):
    """The case's dictionary files as materialize.py renders them, in memory: every file
    under DICT_DIRS of the template case, the *.template files with the tokens substituted
    and the suffix dropped, plus the generated system/decomposeParDict."""
    import materialize
    tokens = materialize.derived_tokens(tokens)
    out = {}
    for d in DICT_DIRS:
        top = os.path.join(base_case, d)
        if not os.path.isdir(top):
            continue
        for root, dirs, files in os.walk(top):
            dirs[:] = [x for x in dirs if not (d == "constant" and root == top and x == "polyMesh")]
            for f in files:
                p = os.path.join(root, f)
                rel = os.path.relpath(p, base_case).replace(os.sep, "/")
                with open(p, "rb") as fh:
                    b = fh.read()
                if f.endswith(".template"):
                    text = b.decode("utf-8", errors="surrogateescape")
                    for key, val in tokens.items():
                        text = text.replace(f"@!{key}!@", str(val))
                    b = text.encode("utf-8", errors="surrogateescape")
                    rel = rel[:-len(".template")]
                out[rel] = b
    out["system/decomposeParDict"] = _decompose_par_dict_text(np_).encode()
    return out


def case_files(case_dir):
    """The same selection read from a materialised case directory."""
    out = {}
    for d in DICT_DIRS:
        top = os.path.join(case_dir, d)
        if not os.path.isdir(top):
            continue
        for root, dirs, files in os.walk(top):
            dirs[:] = [x for x in dirs if not (d == "constant" and root == top and x == "polyMesh")]
            for f in files:
                p = os.path.join(root, f)
                with open(p, "rb") as fh:
                    out[os.path.relpath(p, case_dir).replace(os.sep, "/")] = fh.read()
    return out


def render_stamp_of(files, np_, mode, mesh):
    mesh_files = {k: v for k, v in files.items() if MESH_INPUTS.match(k)}
    gc = ""
    fvs = files.get("system/fvSolution")
    if fvs:
        m = re.search(rb"^\s*grazingCut\s+(\w+)\s*;", fvs, re.M)
        gc = m.group(1).decode() if m else ""
    return {"dicts": files_hash(files), "mesh": files_hash(mesh_files), "np": int(np_), "mode": mode,
            "meshKind": mesh, "nFiles": len(files), "grazingCut": gc}


# ----------------------------------------------------------------------------- the binaries
def _wm_options():
    opt = os.environ.get("WM_OPTIONS")
    if opt:
        return opt
    ds = sorted(os.path.basename(d) for d in glob.glob(os.path.join(REPO, "platforms", "*")) if os.path.isdir(os.path.join(d, "bin")))
    return ds[0] if ds else "linux64GccDPInt32Opt"


def tpf_dir():
    return os.environ.get("LEIA_TPF_DIR") or os.path.join(os.path.dirname(REPO), "TwoPhaseFlow")


def openfoam_dir():
    d = os.environ.get("WM_PROJECT_DIR")
    if d and os.path.isdir(d):
        return d
    ver = ""
    vf = os.path.join(REPO, "platforms", _wm_options(), ".openfoam-version")
    if os.path.isfile(vf):
        ver = open(vf).read().strip()
    cands = [os.path.join(os.path.expanduser("~"), "OpenFOAM", f"OpenFOAM-{ver}")] if ver else []
    cands += sorted(glob.glob(os.path.join(os.path.expanduser("~"), "OpenFOAM", "OpenFOAM-v*")), reverse=True)
    for c in cands:
        if os.path.isdir(os.path.join(c, "platforms")):
            return c
    return ""


def openfoam_version():
    v = os.environ.get("WM_PROJECT_VERSION")
    if v:
        return v
    vf = os.path.join(REPO, "platforms", _wm_options(), ".openfoam-version")
    return open(vf).read().strip() if os.path.isfile(vf) else os.path.basename(openfoam_dir()).replace("OpenFOAM-", "")


def bin_dirs():
    opt = _wm_options()
    return [os.path.join(REPO, "platforms", opt, "bin"), os.path.join(tpf_dir(), "platforms", opt, "bin")]


def lib_dirs():
    opt = _wm_options()
    dirs = [os.path.join(REPO, "platforms", opt, "lib"), os.path.join(tpf_dir(), "platforms", opt, "lib")]
    of = openfoam_dir()
    if of:
        dirs += [os.path.join(of, "platforms", opt, "lib"), os.path.join(of, "platforms", opt, "lib", "sys-openmpi"),
                 os.path.join(of, "platforms", opt, "lib", "dummy")]
        tp = os.environ.get("WM_THIRD_PARTY_DIR") or os.path.join(os.path.dirname(of), "ThirdParty-" + os.path.basename(of).replace("OpenFOAM-", ""))
        dirs += glob.glob(os.path.join(tp, "platforms", opt.replace("Opt", "") + "*", "lib"))
    return [d for d in dirs if os.path.isdir(d)]


def own_lib_dirs():
    """The directories whose libraries are followed transitively: this clone's and the fork's."""
    opt = _wm_options()
    return [os.path.realpath(d) for d in (os.path.join(REPO, "platforms", opt, "lib"), os.path.join(tpf_dir(), "platforms", opt, "lib")) if os.path.isdir(d)]


def locate(name):
    """An application: on PATH when the environment is sourced (a job), else in this clone's
    or the fork's platforms/ bin (the parse-time driver without the OpenFOAM environment)."""
    if os.sep in name:
        return name if os.path.isfile(name) else ""
    p = shutil.which(name)
    if p:
        return os.path.realpath(p)
    for d in bin_dirs():
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return os.path.realpath(p)
    return ""


@functools.lru_cache(maxsize=None)
def needed(path):
    try:
        out = subprocess.run(["readelf", "-d", path], capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return ()
    return tuple(re.findall(r"\(NEEDED\)\s+Shared library: \[([^\]]+)\]", out))


def _resolve(name, dirs):
    for d in dirs:
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return os.path.realpath(p)
    return ""


def library_closure(app_path):
    """{library name: sha256} of every library the application loads from this clone's
    platforms/, the fork's (followed transitively) and OpenFOAM's (as needed by the
    application or by those libraries); the system's libraries are not recorded."""
    dirs, own = lib_dirs(), own_lib_dirs()
    libs, queue, seen = {}, [app_path], set()
    while queue:
        p = queue.pop(0)
        for n in needed(p):
            if n in libs:
                continue
            r = _resolve(n, dirs)
            if not r:
                continue
            libs[n] = sha_file(r)
            if os.path.dirname(r) in own and r not in seen:
                seen.add(r); queue.append(r)
    return libs


def binary_stamp(apps):
    """{app: {path, sha256, libs}} for the named applications (the solver, the initialiser)."""
    out = {}
    for a in apps:
        if not a:
            continue
        p = locate(a)
        if not p:
            out[a] = {"path": "", "sha256": "", "libs": {}, "missing": True}
            continue
        out[a] = {"path": p, "sha256": sha_file(p), "libs": library_closure(p)}
    return out


# ----------------------------------------------------------------------------- the record
def _git(cwd, *args):
    try:
        return subprocess.run(("git",) + args, cwd=cwd, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


@functools.lru_cache(maxsize=None)
def leia_describe():
    d = _git(REPO, "describe", "--always")
    dirty = _git(REPO, "status", "--porcelain", "--", "src", "applications", "workflow", "cases", "config")
    return (d + ("-dirty" if dirty else "")) if d else ""


@functools.lru_cache(maxsize=None)
def tpf_describe():
    t = tpf_dir()
    if not os.path.isdir(os.path.join(t, ".git")) and not os.path.isfile(os.path.join(t, ".git")):
        return ""
    d = _git(t, "describe", "--always")
    dirty = _git(t, "status", "--porcelain", "--", "src", "applications", "apps")
    return (d + ("-dirty" if dirty else "")) if d else ""


def _script_hashes(*commands):
    """The scripts a command names (paths ending in .py or .sh that exist), hashed."""
    out = {}
    for c in commands:
        for tok in re.findall(r"(\S+\.(?:py|sh))\b", c or ""):
            p = tok if os.path.isabs(tok) else os.path.join(REPO, tok)
            if os.path.isfile(p):
                out[os.path.relpath(p, REPO)] = sha_file(p)
    return out


def stamp_id(trigger):
    return sha_bytes(json.dumps(trigger, sort_keys=True).encode())[:12]


def current_stamp(render, spec, apps=None):
    """The whole stamp of a case from its render stamp and the study's solve specification
    (solver, setfields_app, setfields_args, solver_cmd, post_solve, mesh_cmd)."""
    apps = apps or [spec.get("solver", ""), spec.get("setfields_app", "")]
    binaries = binary_stamp(apps)
    trigger = {
        "render": {k: render[k] for k in ("dicts", "mesh", "np", "mode", "meshKind")},
        "binaries": {a: {"sha256": b["sha256"], "libs": b["libs"]} for a, b in binaries.items()},
        "openfoam": openfoam_version(),
        "solve": {"solverCmd": spec.get("solver_cmd", ""), "postSolve": spec.get("post_solve", ""),
                  "meshCmd": spec.get("mesh_cmd", ""), "setfieldsArgs": spec.get("setfields_args", ""),
                  "scripts": _script_hashes(spec.get("post_solve", ""), spec.get("mesh_cmd", ""))},
    }
    return {
        "version": VERSION, "id": stamp_id(trigger), "trigger": trigger,
        "record": {"leia": leia_describe(), "tpf": tpf_describe() if "advectorVoF" in apps or "initDavofFields" in apps else "",
                   "binaryPaths": {a: b["path"] for a, b in binaries.items()},
                   "grazingCut": render.get("grazingCut", ""), "nFiles": render.get("nFiles"),
                   "computedAt": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
                   "host": os.uname().nodename},
    }


# ----------------------------------------------------------------------------- files
def load(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def write_if_changed(path, obj, backdate_to=None, compare_key=None):
    """Write the JSON only when its content differs (the mtime is Snakemake's trigger).
    compare_key: compare this key of the old and the new object instead of the whole
    text (the current stamp is compared on its trigger: its record carries the time of
    the computation, which must not move the file's mtime). backdate_to: a file whose
    mtime minus one second the new file gets, when it was written (adoption: a stamp
    for a result that exists must not reschedule it)."""
    text = json.dumps(obj, indent=1, sort_keys=True)
    old = None
    try:
        with open(path) as fh:
            old = fh.read()
    except OSError:
        pass
    if old == text:
        return False
    if compare_key and old is not None:
        try:
            if json.loads(old).get(compare_key) == obj.get(compare_key):
                return False
        except ValueError:
            pass
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)
    if backdate_to and os.path.isfile(backdate_to):
        t = os.stat(backdate_to).st_mtime - 1
        os.utime(path, (t, t))
    return True


def result_present(case_dir, csv_name):
    """A result: the solver's CSV exists and the launch record carries an exit code."""
    launch = os.path.join(case_dir, ".leia_launch")
    if not os.path.isfile(os.path.join(case_dir, csv_name)) or not os.path.isfile(launch):
        return False
    return any(l.startswith("rc=") for l in open(launch))


def diff_reasons(stored, current):
    """Why two stamps differ, as short phrases."""
    if not stored:
        return ["no stored stamp"]
    if stored.get("id") == current.get("id"):
        return []
    a, b = stored.get("trigger", {}), current.get("trigger", {})
    out = []
    ra, rb = a.get("render", {}), b.get("render", {})
    if ra.get("dicts") != rb.get("dicts"):
        out.append("dictionaries changed" + (" (the mesh inputs among them)" if ra.get("mesh") != rb.get("mesh") else ""))
    for k in ("np", "mode", "meshKind"):
        if ra.get(k) != rb.get(k):
            out.append(f"{k} {ra.get(k)} -> {rb.get(k)}")
    ba, bb = a.get("binaries", {}), b.get("binaries", {})
    for app in sorted(set(ba) | set(bb)):
        x, y = ba.get(app, {}), bb.get(app, {})
        if x.get("sha256") != y.get("sha256"):
            out.append(f"binary {app} rebuilt")
        la, lb = x.get("libs", {}), y.get("libs", {})
        changed = sorted(n for n in set(la) | set(lb) if la.get(n) != lb.get(n))
        if changed:
            out.append(f"libraries of {app} rebuilt: " + ", ".join(changed[:6]) + (" ..." if len(changed) > 6 else ""))
    if a.get("openfoam") != b.get("openfoam"):
        out.append(f"OpenFOAM {a.get('openfoam')} -> {b.get('openfoam')}")
    sa, sb = a.get("solve", {}), b.get("solve", {})
    for k in ("solverCmd", "postSolve", "meshCmd", "setfieldsArgs"):
        if sa.get(k) != sb.get(k):
            out.append(f"{k} changed")
    if sa.get("scripts") != sb.get("scripts"):
        out.append("post-solve scripts changed: " + ", ".join(sorted(set(sa.get("scripts", {})) | set(sb.get("scripts", {}))))[:120])
    return out or ["stamps differ"]


# ----------------------------------------------------------------------------- parse time
def refresh_study(study_dir, case_name, indices, base_case, variations, spec, csv_name):
    """For every case of a study: write the render stamp and the current stamp (content
    change only), adopt a finished case without a stored stamp, and return
    {idx: (status, reasons)} with status in new | same | changed | adopted."""
    status = {}
    for idx, tokens in zip(indices, variations):
        d = os.path.join(study_dir, f"{case_name}_{idx}")
        render = render_stamp_of(rendered_files(base_case, tokens, spec["np"]), spec["np"], spec["mode"], spec["mesh"])
        cur = current_stamp(render, spec)
        stored = load(os.path.join(d, STAMP_FILE))
        finished = result_present(d, csv_name)
        params = os.path.join(d, "case_params.json")
        csv = os.path.join(d, csv_name)
        if finished and stored is None:
            # adoption: the result stands for the current stamp; nothing is rescheduled
            adopted = dict(cur, adopted=True, completedAt="", rc="", jobId="")
            # only when the materialised dictionaries are the ones the templates render now
            on_disk = files_hash(case_files(d)) if os.path.isdir(os.path.join(d, "system")) else ""
            if on_disk == render["dicts"]:
                write_if_changed(os.path.join(d, STAMP_FILE), adopted, backdate_to=csv)
                write_if_changed(os.path.join(d, RENDER_FILE), render, backdate_to=params)
                write_if_changed(os.path.join(d, CURRENT_FILE), cur, backdate_to=csv)
                status[idx] = ("adopted", [])
                continue
            stored = None
        os.makedirs(d, exist_ok=True)
        write_if_changed(os.path.join(d, RENDER_FILE), render)
        write_if_changed(os.path.join(d, CURRENT_FILE), cur, compare_key="trigger")
        if stored is None:
            status[idx] = ("new", [])
        elif stored.get("id") == cur["id"]:
            status[idx] = ("same", [])
        else:
            status[idx] = ("changed", diff_reasons(stored, cur))
    return status


def summary_line(status):
    n = {k: sum(1 for s, _ in status.values() if s == k) for k in ("same", "adopted", "changed", "new")}
    return (f"[leia] stamps: {len(status)} cases: {n['same']} up to date, {n['adopted']} adopted, "
            f"{n['changed']} to recompute (stamp changed), {n['new']} without a result")


# ----------------------------------------------------------------------------- run time
def cmd_check(case_dir, csv_name):
    stored, cur = load(os.path.join(case_dir, STAMP_FILE)), load(os.path.join(case_dir, CURRENT_FILE))
    if cur is None:
        print("[leia] stamp: no current stamp (the Snakefile writes it at parse time); the solver runs"); return 1
    if stored is None:
        print(f"[leia] stamp: no stored stamp; the solver runs under {cur['id']}"); return 1
    if not result_present(case_dir, csv_name):
        print(f"[leia] stamp: no result in the case; the solver runs under {cur['id']}"); return 1
    r = diff_reasons(stored, cur)
    if r:
        print(f"[leia] stamp changed {stored['id']} -> {cur['id']}: " + "; ".join(r) + "; the solver runs"); return 1
    print(f"[leia] stamp {cur['id']} unchanged and a result is present: the solver is NOT run (up to date)")
    return 0


def cmd_archive(case_dir, solver, csv_name):
    stored = load(os.path.join(case_dir, STAMP_FILE))
    if stored is None or not result_present(case_dir, csv_name):
        return 0
    dest = os.path.join(case_dir, "previous", stored.get("id", "unstamped"))
    os.makedirs(dest, exist_ok=True)
    moved = []
    for f in (f"log.{solver}", csv_name, ".leia_launch", ".leia_solver_failed", "log.postSolve", "leia.version", "postProcessing", STAMP_FILE):
        p = os.path.join(case_dir, f)
        if os.path.exists(p):
            shutil.move(p, os.path.join(dest, f)); moved.append(f)
    print(f"[leia] stamp: the previous result ({stored.get('id')}) moved to previous/{stored.get('id')}/: " + ", ".join(moved))
    return 0


def cmd_write(case_dir, apps, rc):
    cur = load(os.path.join(case_dir, CURRENT_FILE))
    if cur is None:
        print("[leia] stamp: no current stamp to store", file=sys.stderr); return 1
    # the binaries and the dictionaries as they are NOW (a rebuild or an edit between the
    # parse and the run is recorded, and the id follows the run, not the plan)
    binaries = binary_stamp(apps)
    trigger = json.loads(json.dumps(cur["trigger"]))
    trigger["binaries"] = {a: {"sha256": b["sha256"], "libs": b["libs"]} for a, b in binaries.items()}
    trigger["openfoam"] = openfoam_version()
    on_disk = files_hash(case_files(case_dir))
    note = {}
    if on_disk != trigger["render"]["dicts"]:
        note["dictsAtParse"] = trigger["render"]["dicts"]
        trigger["render"]["dicts"] = on_disk
        mesh_files = {k: v for k, v in case_files(case_dir).items() if MESH_INPUTS.match(k)}
        trigger["render"]["mesh"] = files_hash(mesh_files)
    if trigger["binaries"] != cur["trigger"].get("binaries"):
        note["binariesAtParse"] = cur["trigger"].get("binaries")
    stored = dict(cur, trigger=trigger, id=stamp_id(trigger))
    stored["record"] = dict(cur.get("record", {}), binaryPaths={a: b["path"] for a, b in binaries.items()},
                            leiaAtRun=leia_describe(), tpfAtRun=tpf_describe() if any(a in ("advectorVoF", "initDavofFields") for a in apps) else "")
    stored.update(completedAt=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
                  rc=rc, jobId=os.environ.get("SLURM_JOB_ID", ""), adopted=False, **note)
    old = load(os.path.join(case_dir, STAMP_FILE))
    if old is not None:
        with open(os.path.join(case_dir, HISTORY_FILE), "a") as fh:
            fh.write(json.dumps(old, sort_keys=True) + "\n")
    with open(os.path.join(case_dir, STAMP_FILE), "w") as fh:
        json.dump(stored, fh, indent=1, sort_keys=True)
    if stored["id"] != cur["id"]:
        # the plan's stamp follows the run's, so the next parse finds them equal
        write_if_changed(os.path.join(case_dir, CURRENT_FILE), dict(cur, trigger=trigger, id=stored["id"]), compare_key="trigger")
        print(f"[leia] stamp: stored {stored['id']} (the plan's {cur['id']} differed: {', '.join(note) or 'binaries'})")
    else:
        print(f"[leia] stamp: stored {stored['id']}")
    return 0


def cmd_report(study_dir, verbose):
    cases = sorted(d for d in glob.glob(os.path.join(study_dir, "*_[0-9][0-9][0-9][0-9][0-9]")) if os.path.isdir(d))
    status = {}
    for d in cases:
        stored, cur = load(os.path.join(d, STAMP_FILE)), load(os.path.join(d, CURRENT_FILE))
        if cur is None:
            status[d] = ("new", ["no current stamp: run `snakemake -n` on the study first"])
        elif stored is None:
            status[d] = ("new", [])
        elif stored.get("adopted") and stored.get("id") == cur.get("id"):
            status[d] = ("adopted", [])
        elif stored.get("id") == cur.get("id"):
            status[d] = ("same", [])
        else:
            status[d] = ("changed", diff_reasons(stored, cur))
    print(summary_line(status))
    for d, (s, r) in status.items():
        if verbose or s in ("changed", "new"):
            st = load(os.path.join(d, STAMP_FILE)) or {}
            cu = load(os.path.join(d, CURRENT_FILE)) or {}
            print(f"  {os.path.basename(d)}: {s} {st.get('id', '-')} -> {cu.get('id', '-')}" + ("; " + "; ".join(r) if r else ""))
    return 0


def main(argv):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check"); c.add_argument("case"); c.add_argument("--csv", required=True)
    w = sub.add_parser("write"); w.add_argument("case"); w.add_argument("--apps", nargs="+", required=True); w.add_argument("--rc", type=int, default=0)
    a = sub.add_parser("archive"); a.add_argument("case"); a.add_argument("--solver", required=True); a.add_argument("--csv", required=True)
    r = sub.add_parser("report"); r.add_argument("study"); r.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    if args.cmd == "check":
        return cmd_check(args.case, args.csv)
    if args.cmd == "write":
        return cmd_write(args.case, args.apps, args.rc)
    if args.cmd == "archive":
        return cmd_archive(args.case, args.solver, args.csv)
    return cmd_report(args.study, args.verbose)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
