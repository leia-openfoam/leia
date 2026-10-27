#!/usr/bin/env python3
"""Restructure the transport part of ms.tex (2026-10-27).

Reads the backup _aside/ms-2026-10-27-before-restructure.tex, replaces the old transport section
(sec:tpf) by sec_tpf_main.tex (short result statements, tables and figures), and moves the old
section's prose paragraphs verbatim into a new appendix "Implementation record of the transport
studies" so that no measurement, number or qualification is lost. Everything else is unchanged.
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "_aside" / "ms-2026-10-27-before-restructure.tex"
MAIN = HERE / "sec_tpf_main.tex"
OUT = HERE / "ms.tex"

if OUT.exists() and "app:transport" in OUT.read_text(encoding="utf-8"):
    sys.exit("ms.tex already restructured (one-time script, 2026-10-27); later edits live in ms.tex itself")
text = SRC.read_text(encoding="utf-8")
lines = text.split("\n")

# locate the old section boundaries
i_tpf = next(i for i, l in enumerate(lines) if l.startswith(r"\section{Transport: the plane-clipped loop"))
i_tr = next(i for i, l in enumerate(lines) if l.startswith(r"\section{The paraboloid update"))
i_bib = next(i for i, l in enumerate(lines) if l.startswith(r"\bibliographystyle"))

old = "\n".join(lines[i_tpf:i_tr])
paras = [p for p in re.split(r"\n\s*\n", old) if p.strip()]

expected = [
    r"\section{Transport",                    # 0
    "The state has been transported",          # 1
    r"\emph{The normal in the loop.}",         # 2
    r"\begin{table}",                          # 3 plane table
    r"\begin{table}",                          # 4 sphere table
    r"\begin{figure}",                         # 5 history figure
    r"\emph{Plane} (\cref{tab:tpf-plane})",    # 6
    r"\emph{Sphere}",                          # 7
    r"\emph{Perturbed hexahedra.}",            # 8
    r"\emph{A second-order correction",        # 9
    "The per-face fields",                     # 10
    r"\emph{The dual advection with one surface",  # 11
    "The exact sphere as every cell",          # 12
    "The six-coefficient surface",             # 13
    "The surface that is stable",              # 14
    r"\emph{The volume rules.}",               # 15
    r"\emph{Parallel.}",                       # 16
    r"\emph{The control over the ladder}",     # 17
    r"\begin{table}",                          # 18 curved table
    r"\begin{table}",                          # 19 curved-plane table
    r"\emph{The ladders}",                     # 20
    "With the ring normal in the frame",       # 21
    r"\emph{The flux as a quadrature in time.}",  # 22
    "The order of the closed loop is thus",    # 23
    r"\emph{Plane} (\cref{tab:tpf-curved-plane})",  # 24
    "With the cell's own normal in the frame", # 25
    r"\emph{What this says.}",                 # 26
]
if len(paras) != len(expected):
    sys.exit(f"paragraph count {len(paras)} != {len(expected)}")
for k, (p, e) in enumerate(zip(paras, expected)):
    if not p.lstrip().startswith(e):
        sys.exit(f"paragraph {k} does not start with {e!r}: {p[:60]!r}")


def strip_lead(p):
    """Drop a run-in \\emph{...} heading at the start of a moved paragraph (the appendix has its own headings)."""
    return re.sub(r"^\\emph\{[^}]*\}\s*(\([^)]*\)\s*)?\.?\s*", "", p.lstrip(), count=1)


appendix = [
    r"\section{Implementation record of the transport studies}",
    r"\label{app:transport}",
    "",
    r"This appendix is the record of the measurements of \cref{sec:tpf}, in the order of that section: what was built, how each quantity was computed, which defects were found on the way, and every number behind the tables there. It is written for whoever re-runs or extends the studies.",
    "",
    r"\subsection{The plane-clipped loop in TwoPhaseFlow}",
    r"\label{app:tpf-plane-loop}",
    "",
    paras[1].strip(), "",
    strip_lead(paras[2]), "",
    r"\emph{The translated plane.} " + strip_lead(paras[6]), "",
    r"\emph{The translated sphere.} " + strip_lead(paras[7]), "",
    r"\subsection{Perturbed hexahedra and the wetted area vector per face}",
    r"\label{app:tpf-warped}",
    "",
    strip_lead(paras[8]), "",
    r"\subsection{The sagitta correction within the plane loop}",
    r"\label{app:tpf-sagitta}",
    "",
    strip_lead(paras[9]), "",
    paras[10].strip(), "",
    r"\subsection{The dual advection: quadrature, the three surfaces, the volume anchor}",
    r"\label{app:tpf-dual}",
    "",
    strip_lead(paras[11]), "",
    paras[12].strip(), "",
    paras[13].strip(), "",
    paras[14].strip(), "",
    r"\subsection{The volume rules}",
    r"\label{app:tpf-rules}",
    "",
    strip_lead(paras[15]), "",
    r"\subsection{Parallel exchange}",
    r"\label{app:tpf-parallel}",
    "",
    strip_lead(paras[16]), "",
    r"\subsection{The exact-sphere control and the reference point of the normal error}",
    r"\label{app:tpf-control}",
    "",
    strip_lead(paras[17]), "",
    r"\subsection{The ladders}",
    r"\label{app:tpf-ladders}",
    "",
    strip_lead(paras[20]), "",
    paras[21].strip(), "",
    paras[23].strip(), "",
    r"\emph{The translated plane, and the defects found by it.} " + strip_lead(paras[24]), "",
    paras[25].strip(), "",
    r"\subsection{The flux as a quadrature in time}",
    r"\label{app:tpf-flux}",
    "",
    strip_lead(paras[22]), "",
]

main = MAIN.read_text(encoding="utf-8").rstrip("\n")

new = lines[:i_tpf] + main.split("\n") + [""] + lines[i_tr:i_bib] + appendix + lines[i_bib:]
OUT.write_text("\n".join(new), encoding="utf-8")
print(f"wrote {OUT}: {len(new)} lines; moved {len(appendix)} appendix lines; main section {len(main.split())} words")
