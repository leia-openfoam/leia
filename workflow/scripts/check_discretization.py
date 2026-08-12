#!/usr/bin/env python3
"""Assert that every case of a study group ran ONE discretization.

This is the check that would have caught the confound that cost a study round:
the 2D SDPLS arms ran `div(phi,psi) Gauss linearUpwind grad(psi)` while the new
3D configs pinned `Gauss linear`, so the 2D and 3D numbers were not comparable
-- and nothing said so. Comparing arms across different discretizations is not a
method comparison.

It reads the RENDERED `system/fvSchemes` / `system/controlDict` of each case
(via fvschemes.py), never the study tokens, because the tokens are not a
faithful record: 2Dvortex hardcodes its div scheme and has no DIV token at all.

Usage (from the repo root):
    python3 workflow/scripts/check_discretization.py                # all studies
    python3 workflow/scripts/check_discretization.py sdplsStability benchVortexEulerT2
    make check-discretization

Exit status is 1 if any case deviates or if a group is internally inconsistent,
so it is usable as a gate.
"""
import collections
import glob
import json
import os
import sys

import fvschemes

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
STUDIES = os.path.join(REPO, "studies")

# What every finite-volume div(phi,psi) study is expected to have run after the
# 2026-08 uniformisation. A study whose cases disagree with this is reported;
# a study whose cases disagree with EACH OTHER is reported even more loudly.
EXPECTED = {
    "divPsi": "Gauss linearUpwind grad(psi)",
    "divPsiGradScheme": "cellLimited leastSquares 1",
}
EXPECTED_MIN_NDEFCORR = 3

# Studies that vary the discretization ON PURPOSE. For these, "the cases used N
# different discretizations" is the experiment, not the bug this gate exists to
# catch -- but they are still checked for everything else, and the deviation is
# printed so an ablation can never be mistaken for a uniform study.
ABLATION_STUDIES = {"sdplsOrderAblation"}


def _is_semi_lagrangian(case_dir):
    """A semi-Lagrangian case never assembles div(phi,psi) -- its update is a
    pointwise assignment -- yet the shared case templates still DEFINE the key,
    so reading fvSchemes alone would flag SL studies as deviations."""
    meta = os.path.join(case_dir, "case_params.json")
    if os.path.isfile(meta):
        try:
            tokens = json.load(open(meta)).get("tokens", {})
        except (OSError, ValueError):
            tokens = {}
        if tokens.get("ADVECTION") == "semiLagrangian":
            return True
    return any(os.path.isfile(os.path.join(case_dir, n)) for n in
               ("leiaSemiLagrangeLevelSetFoam.csv",
                "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"))


def _cases(study_dir):
    return sorted(d for d in glob.glob(os.path.join(study_dir, "*"))
                  if os.path.isdir(d) and os.path.isfile(
                      os.path.join(d, "system", "fvSchemes")))


def check_study(study_dir):
    """(n_fv_cases_checked, list of complaint strings) for one study."""
    name = os.path.basename(study_dir)
    cases = _cases(study_dir)
    if not cases:
        return 0, []

    seen = collections.defaultdict(list)   # discretization tuple -> case names
    complaints = []
    for case in cases:
        if _is_semi_lagrangian(case):
            continue
        d = fvschemes.read_discretization(case)
        # A case with no div(phi,psi) key at all is likewise not FV advection.
        if not d["divPsi"]:
            continue
        key = (d["divPsi"], d["divPsiGradScheme"], d["nDefCorr"])
        seen[key].append(os.path.basename(case))

    if not seen:
        return 0, []      # no FV-advection case here (e.g. a semi-Lagrangian study)

    ablation = name in ABLATION_STUDIES
    if len(seen) > 1 and not ablation:
        complaints.append(
            f"{name}: cases used {len(seen)} DIFFERENT discretizations -- "
            "arms are not comparable:")
        for key, names in sorted(seen.items(), key=lambda kv: -len(kv[1])):
            complaints.append(
                f"    {len(names):3d} case(s)  div={key[0]!r} "
                f"grad={key[1]!r} nDefCorr={key[2]!r}"
                f"   e.g. {names[0]}")
    elif ablation:
        # Say what it varied. Silence here would make an ablation study
        # indistinguishable from a uniform one in the gate's output.
        print(f"[check-disc] {name}: ABLATION study, {len(seen)} discretization(s) "
              f"by design:")
        for key, names in sorted(seen.items(), key=lambda kv: -len(kv[1])):
            print(f"    {len(names):3d} case(s)  grad={key[1]!r}")

    for (div, grad, ndc), names in seen.items():
        if ablation and div == EXPECTED["divPsi"]:
            # The div scheme is still pinned; only the reconstruction gradient
            # is swept, so check everything except that one field.
            try:
                n = int(str(ndc).split()[0])
            except (ValueError, IndexError):
                n = -1
            if n < EXPECTED_MIN_NDEFCORR:
                complaints.append(
                    f"{name}: nDefCorr = {ndc!r}, expected >= "
                    f"{EXPECTED_MIN_NDEFCORR}  ({len(names)} case(s))")
            continue
        if div != EXPECTED["divPsi"]:
            complaints.append(
                f"{name}: div(phi,psi) = {div!r}, expected "
                f"{EXPECTED['divPsi']!r}  ({len(names)} case(s))")
        if grad != EXPECTED["divPsiGradScheme"]:
            complaints.append(
                f"{name}: the gradient named by div(phi,psi) resolves to "
                f"{grad!r}, expected {EXPECTED['divPsiGradScheme']!r}  "
                f"({len(names)} case(s))")
        try:
            n = int(str(ndc).split()[0])
        except (ValueError, IndexError):
            n = -1
        if n < EXPECTED_MIN_NDEFCORR:
            complaints.append(
                f"{name}: nDefCorr = {ndc!r}, expected >= "
                f"{EXPECTED_MIN_NDEFCORR} (linearUpwind's second-order part is "
                f"an explicit source)  ({len(names)} case(s))")
    return sum(len(v) for v in seen.values()), complaints


def main(argv):
    names = argv or sorted(
        os.path.basename(d) for d in glob.glob(os.path.join(STUDIES, "*"))
        if os.path.isdir(d))
    if not os.path.isdir(STUDIES):
        print(f"[check-disc] no {STUDIES} -- nothing to check")
        return 0

    all_complaints, n_fv, n_studies = [], 0, 0
    for name in names:
        d = os.path.join(STUDIES, name)
        if not os.path.isdir(d):
            print(f"[check-disc] no study directory {name}; skip")
            continue
        n, complaints = check_study(d)
        n_fv += n
        n_studies += 1 if n else 0
        if complaints:
            all_complaints += complaints
    if all_complaints:
        print("[check-disc] DISCRETIZATION MISMATCH")
        for c in all_complaints:
            print("  " + c)
        return 1
    if not n_fv:
        print("[check-disc] no finite-volume div(phi,psi) cases found "
              "(nothing to check -- semi-Lagrangian studies are skipped)")
        return 0
    print(f"[check-disc] OK: {n_fv} FV cases across {n_studies} studies, "
          f"one discretization "
          f"({EXPECTED['divPsi']}, grad {EXPECTED['divPsiGradScheme']}, "
          f"nDefCorr >= {EXPECTED_MIN_NDEFCORR})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
