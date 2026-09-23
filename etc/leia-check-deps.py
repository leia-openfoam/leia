#!/usr/bin/env python3
"""Gate: the header includes between the parts of src/leiaLevelSet must follow the
library link graph of docs/plan-library-split-and-build-policy.md (WP3).

A part may include headers of its own library, of libleiaCore, and of the libraries
its Make/options links. Anything else is an undeclared dependency: it would compile
(the root lnInclude holds every header) and fail only at link time (--no-undefined)
or, for a header-only use, never. Allwmake runs this before the first wmake and
stops on a violation. Prints the measured table; exit 1 on a violation or on a
part directory this table does not know."""
import os, re, sys, collections

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "leiaLevelSet"))
LIB = {   # part directory -> library
    "(root)": "leiaCore", "profile": "leiaCore", "narrowBand": "leiaCore",
    "phaseIndicator": "leiaCore", "velocityModel": "leiaCore", "schemes": "leiaCore",
    "sdplsSource": "leiaSdplsSource", "semiLagrangian": "leiaSemiLagrangian",
    "velocityExtension": "leiaVelocityExtension", "redistancer": "leiaRedistancer",
    "volumeCorrection": "leiaVolumeCorrection", "surfaceTensionForce": "leiaSurfaceTension",
    "advection": "leiaAdvection",
}
LINKS = {  # library -> the leia libraries in its LIB_LIBS (libleiaCore is implicit)
    "leiaCore": set(), "leiaSdplsSource": set(), "leiaSemiLagrangian": set(),
    "leiaRedistancer": set(), "leiaVolumeCorrection": set(), "leiaSurfaceTension": set(),
    "leiaVelocityExtension": {"leiaSemiLagrangian"},
    "leiaAdvection": {"leiaSdplsSource", "leiaSemiLagrangian", "leiaVelocityExtension"},
}

def part_of(path):
    rel = os.path.relpath(path, ROOT)
    return "(root)" if os.sep not in rel else rel.split(os.sep)[0]

def sources():
    for dp, dn, fn in os.walk(ROOT):
        dn[:] = [d for d in dn if d not in ("Make", "lnInclude")]
        for f in fn:
            if f.endswith((".C", ".H")) and not f.startswith("leiaStamp_"):
                yield os.path.join(dp, f)

header_part = {os.path.basename(p): part_of(p) for p in sources() if p.endswith(".H")}
inc = re.compile(r'#include\s+"([^"]+)"')
uses = collections.defaultdict(lambda: collections.defaultdict(set))
unknown = set()
for p in sources():
    src = part_of(p)
    if src not in LIB:
        unknown.add(src); continue
    for m in inc.finditer(open(p, errors="ignore").read()):
        dst = header_part.get(os.path.basename(m.group(1)))
        if dst is not None and dst != src:
            uses[src][dst].add(os.path.basename(p))

bad = []
print("leia-check-deps: header includes between the parts of src/leiaLevelSet")
for src in sorted(LIB):
    row = uses.get(src, {})
    allowed = LINKS[LIB[src]] | {LIB[src], "leiaCore"}
    for dst in sorted(row):
        ok = LIB[dst] in allowed
        print(f"  {src:20s} -> {dst:20s} ({LIB[dst]:22s}) {'ok ' if ok else 'BAD'} {', '.join(sorted(row[dst]))}")
        if not ok:
            bad.append((src, dst))
for u in sorted(unknown):
    print(f"  {u:20s} : part directory not in the table of etc/leia-check-deps.py")
if bad or unknown:
    print("leia-check-deps: FAIL -- add the library to LIB_LIBS and to LINKS here, or move the include")
    sys.exit(1)
print("leia-check-deps: PASS")
