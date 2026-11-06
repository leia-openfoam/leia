#!/usr/bin/env python3
"""Seed the face-fraction state with a random perturbation (2026-11-08, the
stability test of the ring normal's stencil variants on the translated plane,
config/davof/stencil_plane_growth.yaml): every internal face with
1e-8 < alphaf < 1 - 1e-8 (the partial faces of the interface cells) gets
alphaf + eps u, u uniform in [-1, 1] from a fixed seed (one draw per partial
face in face order; the same seed gives the same perturbation on every
variant and resolution), clipped to [0, 1]; the full and empty faces and the
cell fractions are left exact. The perturbation enters the loop only through
the normal m_k = -sum_f alpha_f S_f and the planes positioned along it: the
growth of L1(alphaf - alphaf_exact) per step then measures the contraction
of the normal's feedback (the tilt diagnostic of 2026-10-26 measured the own
normal's growth from round-off; a seeded 1e-8 amplitude lifts the mode above
round-off from step 1 so that 50 steps give the factor). The field is
rewritten with 17 significant digits. Runs in the case directory after the
field initialisation and before decomposePar (a study's setfields_args,
e.g. "; python3 {REPO}/workflow/scripts/davof_perturb_alphaf.py
0/alphaf.water 1e-8 20261108").

    python3 davof_perturb_alphaf.py <field file> <eps> [<seed>]
"""
import random
import re
import sys


def main(argv):
    path, eps = argv[1], float(argv[2])
    seed = int(argv[3]) if len(argv) > 3 else 20261108
    text = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?\s*(\d+)\s*\n?\(", text)
    if not m:
        sys.exit(f"davof_perturb_alphaf: no nonuniform internalField in {path}")
    n = int(m.group(1))
    start = m.end()
    end = text.index(")", start)
    vals = text[start:end].split()
    if len(vals) != n:
        sys.exit(f"davof_perturb_alphaf: {len(vals)} values, header says {n}")
    rng = random.Random(seed)
    out, k, amax = [], 0, 0.0
    for v in vals:
        a = float(v)
        if 1e-8 < a < 1 - 1e-8:
            da = eps * rng.uniform(-1.0, 1.0)
            a = min(1.0, max(0.0, a + da))
            amax = max(amax, abs(da))
            k += 1
        out.append(repr(a))
    text = text[:start] + "\n" + "\n".join(out) + "\n" + text[end:]
    open(path, "w").write(text)
    print(f"davof_perturb_alphaf: {k} of {n} faces of {path} perturbed, "
          f"amplitude {eps} (max |da| {amax:.3e}), seed {seed}")


if __name__ == "__main__":
    main(sys.argv)
