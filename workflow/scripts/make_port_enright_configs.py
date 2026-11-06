#!/usr/bin/env python3
"""The identity acceptance of the Enright agent's decisions 1-3 in leia (2026-11-07): the
configs config/davof/port_enrightKeys_{rotation,diagonal,plane}_{fork,leia}.yaml, the
eight-case identity set of PORT-DESIGN 4.1 (the rotation N 32 and 64, the diagonal, the
plane N 20 and 40 on both arms and both rows) with every setting of
tpf_rows.ENRIGHT_SETTINGS as an axis, on the fork route (advectorVoF built from
davof-swept-enright 25e6ebc, LEIA_TPF_DIR) and on the leia route (davofAdvectorFoam); and
port_enrightKeys_par{1,2,4,8}_{fork,leia}.yaml, the parallel identity of the recommended
setting (the rotation N 32, 20 steps, both rows). Written from the port_* configs with
the study names and the axis changed, so that everything else is the acceptance's.

    python3 workflow/scripts/make_port_enright_configs.py
"""
import os, re, sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
C = os.path.join(REPO, "config", "davof")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tpf_rows

HEAD = ("# The identity acceptance of the Enright agent's decisions 1-3 (fork davof-swept-enright\n"
        "# 25e6ebc, leia 2026-11-07; workflow/scripts/make_port_enright_configs.py writes this file):\n"
        "# the port's identity set with the setting TPF_ENRIGHT_SETTING as an axis; the fork route\n"
        "# runs advectorVoF built from 25e6ebc (LEIA_TPF_DIR), the leia route davofAdvectorFoam;\n"
        "# the two must be bit for bit for every setting (port_identity_check.py).\n")


def derive(src, dst, study, axis):
    s = open(os.path.join(C, src)).read()
    s = re.sub(r"^(#.*\n)+", HEAD, s, count=1)
    s = re.sub(r"^study_name: .*$", f"study_name: {study}", s, flags=re.M)
    # grazingCut fork: the identity against the fork (PORT-DESIGN 1.2 item 1; the fork's
    # binary ignores the key); without it the rotation N 64 separates at a grazing cut
    # (measured 2026-11-06: 7 of 20 rotation cases differed, all at N 64, the diagonal, the
    # plane and the parallel studies identical)
    s = s.replace("axes_override:\n", "axes_override:\n  TPF_ENRIGHT_SETTING: [" + ", ".join(axis) + "]\n  TPF_GRAZING_CUT: [fork]\n", 1)
    assert "TPF_ENRIGHT_SETTING" in s
    with open(os.path.join(C, dst), "w") as fh:
        fh.write(s)
    print("written:", dst)


def main():
    settings = list(tpf_rows.ENRIGHT_SETTINGS)
    for fam in ("rotation", "diagonal", "plane"):
        for route in ("fork", "leia"):
            derive(f"port_{fam}_{route}.yaml", f"port_enrightKeys_{fam}_{route}.yaml", f"port_enrightKeys_{fam}_{route}", settings)
    for np_ in (1, 2, 4, 8):
        for route in ("fork", "leia"):
            derive(f"port_parN32np{np_}_{route}.yaml", f"port_enrightKeys_par{np_}_{route}.yaml", f"port_enrightKeys_par{np_}_{route}", ["recommended"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
