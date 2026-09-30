#!/usr/bin/env python3
"""Check that a leia build produced every target that a Make/files declares.

`wmake all applications` does not stop, and does not always exit non-zero, when one
application fails to compile, so the exit code of ./Allwmake is not evidence that every
solver and every library was built. This script reads every src/**/Make/files and
applications/**/Make/files, takes its `EXE = $(FOAM_USER_APPBIN)/<name>` or
`LIB = $(FOAM_USER_LIBBIN)/<name>` line, and checks that the file exists in the clone's own
install (etc/leia-env.sh sets FOAM_USER_APPBIN and FOAM_USER_LIBBIN to <clone>/platforms/).

Usage (after ./Allwmake, in a shell that sourced OpenFOAM's etc/bashrc and etc/leia-env.sh):
    python3 etc/leia-check-build.py
Exit 0 when every target exists, 1 when one is missing, 2 when the environment is not set.
Written 2026-09-30 for the CI (.github/scripts/ci-build-and-smoke.sh).
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TARGET = re.compile(r"^\s*(EXE|LIB)\s*=\s*\$\((FOAM_USER_APPBIN|FOAM_USER_LIBBIN)\)/(\S+)\s*$")


def main():
    appbin, libbin = os.environ.get("FOAM_USER_APPBIN"), os.environ.get("FOAM_USER_LIBBIN")
    if not appbin or not libbin:
        print("leia-check-build: FOAM_USER_APPBIN/FOAM_USER_LIBBIN not set; source etc/bashrc and etc/leia-env.sh")
        return 2
    dirs = {"FOAM_USER_APPBIN": appbin, "FOAM_USER_LIBBIN": libbin}
    found, missing, unparsed = [], [], []
    for top in ("src", "applications"):
        for d, _, files in os.walk(os.path.join(ROOT, top)):
            if os.path.basename(d) != "Make" or "files" not in files or "lnInclude" in d:
                continue
            rel = os.path.relpath(os.path.dirname(d), ROOT)
            targets = [TARGET.match(l) for l in open(os.path.join(d, "files"))]
            targets = [m for m in targets if m]
            if not targets:
                unparsed.append(rel)
                continue
            kind, var, name = targets[-1].groups()
            path = os.path.join(dirs[var], name + (".so" if kind == "LIB" else ""))
            (found if os.path.exists(path) else missing).append((kind, name, rel))
    for kind, name, rel in sorted(missing):
        print(f"leia-check-build: MISSING {kind} {name}  ({rel}/Make/files)")
    for rel in unparsed:
        print(f"leia-check-build: no EXE/LIB line in {rel}/Make/files")
    n_exe = sum(1 for k, _, _ in found if k == "EXE")
    n_lib = len(found) - n_exe
    print(f"leia-check-build: {n_exe} executables and {n_lib} libraries present, "
          f"{len(missing)} missing, {len(unparsed)} unreadable Make/files")
    return 1 if missing or unparsed else 0


if __name__ == "__main__":
    sys.exit(main())
