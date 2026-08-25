# Agent guide

This repo's agent instructions live in **[CLAUDE.md](CLAUDE.md)** (execution
environment, git discipline, running studies) and the cluster workflow in
**[CLUSTER.md](CLUSTER.md)** (Lichtenberg / TU Darmstadt: SSH, SLURM account
`special00004`, OpenFOAM-v2512, rsync helpers).

Read both before making changes or running simulations. Key rules:

- Run OpenFOAM/Python/Snakemake in **WSL**, not Windows (`wsl bash -lc '...'`).
- Code moves by **git** (hub: GitHub `leia-openfoam/leia`); raw simulation output
  (`studies/`, `runs/`) is git-ignored and moves by **rsync** — never `git add` it.
- Build against `$HOME/OpenFOAM/OpenFOAM-v2512`; run studies with
  `make studies PROFILE=profiles/{local,slurm}`.

## Method constraint: unstructured FVM only

**Never propose or pursue a method that depends on a structured or Cartesian
discretisation, on a named mode, or on a specific interface shape.** Everything
must work for the **unstructured Finite Volume Method as OpenFOAM implements it**
— or at minimum for unstructured methods that admit OpenFOAM meshes, use
**compact stencils**, and parallelise efficiently under **MPI domain
decomposition**.

Ruled out, non-exhaustively:

- filtering or projecting out named azimuthal/Fourier modes (e.g. "remove the
  m = 4, 8, 16 mesh-locked curvature content") — those modes exist only because
  the mesh is Cartesian and the interface is a circle or sphere;
- height-function constructions needing mesh-aligned columns;
- anything requiring a structured index space (i, j, k), a global transform, a
  wide or unbounded stencil, or a non-local solve that will not decompose cleanly
  across ranks;
- corrections tuned to one benchmark geometry.

**Why.** A mode- or mesh-specific correction can look excellent on the
stationary-droplet test and be worthless in general — and worse, it hides the
underlying defect behind a tuning parameter nobody can justify later. The target
is a method for arbitrary unstructured meshes and arbitrary interface topology; a
fix that does not survive that transition is not a fix.

**No partial solutions.** Do not bank a partial gain while the underlying defect
remains — e.g. do not raise the capillary time-step coefficient for a 2x speedup
when the growth it papers over is still there. That bakes a tuned number into
every case file for a problem that is about to be removed.

**How to apply.** When proposing a fix, state explicitly how it behaves on a
polyhedral unstructured mesh under MPI decomposition *before* proposing it.
Prefer formulations stated per cell or per face and summed — those are
mesh-agnostic by construction. Worked example that passes the test: the
variational capillary force, `f_c = -sigma * dA_h/dpsi_c` derived from a discrete
surface-energy functional `A_h = sum_c V_c |grad_h alpha_c|` — a sum over cells,
compact stencil, no structure assumed, no modes named, and zero net work around
any closed cycle by construction.
