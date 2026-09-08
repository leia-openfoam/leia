#!/usr/bin/env python3
"""Popinet (JCP 228, 2009, Sec. 6.2.2) translating droplet: the SI parameter set.

The benchmark is defined by dimensionless groups only: We = rho U^2 D / sigma = 0.4,
La = sigma D / (rho nu^2) in {120, 1200, 12000, inf}, density ratio 1, viscosity ratio 1,
D = 0.4 of the box height, one diameter of travel (T_U = D / U). Popinet runs it with
rho = sigma = U = 1 and D = 0.4. This repository runs a DIMENSIONAL twin in SI units, so
that every dictionary carries a real droplet and a real velocity. The set is fixed here,
in one place, and every Popinet config is checked against it (`check`).

The three free choices (five dimensional parameters rho, mu, sigma, D, U minus two groups
La, We leave three), taken 2026-09-08:

    D    = 1.0e-3 m       droplet diameter (R = 0.5 mm)
    rho  = 1000 kg/m^3    both phases
    nu   = 1.0e-6 m^2/s   both phases (mu = 1.0e-3 Pa s, water)

Everything else follows from La and We:

    sigma = La rho nu^2 / D          = 0.012 N/m      (an oil/water interfacial tension)
    U     = sqrt(We sigma / (rho D)) = 0.069282 m/s
    Re    = rho U D / mu = sqrt(La We) = 69.28
    H     = D / 0.4 = 2.5 mm (box height and depth), L = 2 H = 5 mm (box length)
    T_U   = D / U = 14.434 ms (horizon: one diameter of travel)

Time step. The repository law is dt = CAPILLARY_DT_COEFF / nRef^1.5 with
nRef = CAPILLARY_REF_LENGTH * N_CELLS / DOMAIN_LENGTH (workflow/scripts/materialize.py).
With CAPILLARY_REF_LENGTH = DOMAIN_LENGTH = H, nRef = N_CELLS and

    CAPILLARY_DT_COEFF = c * sqrt((rho1 + rho2) H^3 / (2 pi sigma)),   c = 0.2323,

the same safety factor as every production run (0.010861 on the 0.01 m water/air box):
dt is 0.2323 of sqrt((rho1 + rho2) h^3 / (2 pi sigma)) at the interface cell h = H / N.
For a polyhedral mesh N_CELLS is pinned to H / h_interface by the band check
(workflow/scripts/check_refined_band.py --mode poly), so the same law holds there.

Modes
    popinet_si.py print                      the set, the scale factors, dt per N
    popinet_si.py yaml                       axes_override lines with the SI values
    popinet_si.py check config/a.yaml ...    verify the dimensional axes of each config
                                             against the set; exit 1 on any mismatch
    popinet_si.py rewrite config/a.yaml ...  convert a config from the OLD dimensionless
                                             set (rho = sigma = U = 1, D = 0.4) to SI, in
                                             place; a config already in SI is left alone

Options --D --rho --nu --La --We --safety --aspect change the set for print/yaml/check.
Stdlib only.
"""
import argparse
import math
import re
import sys

# The dimensionless reference Popinet used (and this repository used until 2026-09-08).
ND = dict(D=0.4, rho=1.0, sigma=1.0, U=1.0, H=1.0)

SURFACE_RENAME = {
    "box2x1x1.fms": "box5x2p5x2p5mm.fms",
    "box2x1x1.stl": "box5x2p5x2p5mm.stl",
}


def fnum(x):
    """Number for a YAML list or an OpenFOAM dictionary: 10 significant digits, and a
    decimal point in the mantissa so that PyYAML (YAML 1.1) reads it as a float."""
    if x == 0:
        return "0"
    s = "{:.10g}".format(x)
    if "e" in s:
        m, e = s.split("e")
        if "." not in m:
            m += ".0"
        s = m + "e" + e
    return s


class SISet:
    def __init__(self, D=1.0e-3, rho=1000.0, nu=1.0e-6, La=12000.0, We=0.4,
                 safety=0.2323, aspect=2.0):
        self.D, self.rho, self.nu, self.La, self.We = D, rho, nu, La, We
        self.safety, self.aspect = safety, aspect
        self.mu = rho * nu
        self.sigma = La * rho * nu ** 2 / D
        self.U = math.sqrt(We * self.sigma / (rho * D))
        self.Re = rho * self.U * D / self.mu
        self.Oh = 1.0 / math.sqrt(La)
        self.H = D / 0.4                      # Popinet: D = 0.4 of the unit square
        self.L = aspect * self.H
        self.TU = D / self.U
        self.R = 0.5 * D
        # CAPILLARY_DT_COEFF for CAPILLARY_REF_LENGTH = DOMAIN_LENGTH = H (nRef = N)
        self.coeff = safety * math.sqrt(2.0 * rho * self.H ** 3 / (2.0 * math.pi * self.sigma))
        # scale factors from the dimensionless reference
        self.lamL = self.H / ND["H"]
        self.lamU = self.U / ND["U"]
        self.lamT = self.lamL / self.lamU
        self.lamNu = self.lamL * self.lamU
        self.lamSigma = self.sigma / ND["sigma"]
        self.lamRho = rho / ND["rho"]

    def nu_for_La(self, La):
        """Viscosity of another Laplace number at the same rho, sigma, D (the La sweep)."""
        return 0.0 if math.isinf(La) else math.sqrt(self.sigma * self.D / (self.rho * La))

    def dt(self, N):
        return self.coeff / N ** 1.5

    def brackbill(self, h):
        return math.sqrt(2.0 * self.rho * h ** 3 / (2.0 * math.pi * self.sigma))

    def print(self, out=sys.stdout):
        s = self
        w = out.write
        w("Popinet translating droplet, SI set\n")
        w(f"  D      = {fnum(s.D)} m        R = {fnum(s.R)} m\n")
        w(f"  rho    = {fnum(s.rho)} kg/m^3   (both phases)\n")
        w(f"  nu     = {fnum(s.nu)} m^2/s    mu = {fnum(s.mu)} Pa s (both phases)\n")
        w(f"  sigma  = {fnum(s.sigma)} N/m      = La rho nu^2 / D\n")
        w(f"  U      = {fnum(s.U)} m/s      = sqrt(We sigma / (rho D))\n")
        w(f"  La = {fnum(s.La)}  We = {fnum(s.We)}  Re = {s.Re:.4f} (= sqrt(La We) = {math.sqrt(s.La * s.We):.4f})  Oh = {s.Oh:.4e}\n")
        w(f"  box    = {fnum(s.L)} x {fnum(s.H)} x {fnum(s.H)} m  (DOMAIN_LENGTH = H = {fnum(s.H)}, POPINET_XLEN = {fnum(s.aspect)} heights)\n")
        w(f"  T_U    = D / U = {fnum(s.TU)} s   (END_TIME for one diameter of travel; T_U/2 = {fnum(s.TU / 2)})\n")
        w(f"  CAPILLARY_DT_COEFF = {fnum(s.coeff)}  = {s.safety} sqrt((rho1+rho2) H^3 / (2 pi sigma))\n")
        w("  scale factors from the dimensionless twin (rho = sigma = U = 1, D = 0.4):\n")
        w(f"    length x {fnum(s.lamL)}   velocity x {fnum(s.lamU)}   time x {fnum(s.lamT)}   nu x {fnum(s.lamNu)}   pressure x {fnum(s.rho * s.U ** 2)}\n")
        w("  La sweep at fixed rho, sigma, D:\n")
        for La in (120.0, 1200.0, 12000.0, float("inf")):
            w(f"    La = {La:<8g} nu = {fnum(s.nu_for_La(La))}\n")
        w("  per rung (h = H / N at the interface; dt = coeff / N^1.5):\n")
        for N in (64, 96, 128, 256):
            h, dt = s.H / N, s.dt(N)
            w(f"    N = {N:>3}: h = {h:.4e} m  dt = {fnum(dt)} s  steps to T_U = {s.TU / dt:7.1f}"
              f"  dt/Brackbill = {dt / s.brackbill(h):.4f}  Co = {s.U * dt / h:.4f}"
              f"  poly MAX_CELL_SIZE = h / 2^(-1/3) = {fnum(h / 2 ** (-1.0 / 3.0))}\n")

    def axes(self):
        """The SI values of the dimensional axes, in axes_override form."""
        s = self
        return [
            ("DOMAIN_LENGTH", [fnum(s.H)]),
            ("CAPILLARY_REF_LENGTH", [fnum(s.H)]),
            ("CAPILLARY_DT_COEFF", [fnum(s.coeff)]),
            ("DROPLET_RADIUS", [fnum(s.R)]),
            ("DROPLET_OFFSET_X", ["0"]),
            ("TRANSLATION_SPEED", [fnum(s.U)]),
            ("SIGMA", [fnum(s.sigma)]),
            ("POPINET_RHO", [fnum(s.rho)]),
            ("POPINET_NU", [fnum(s.nu)]),
            ("POPINET_XLEN", [fnum(s.aspect)]),
            ("END_TIME", [fnum(s.TU)]),
            ("WRITE_INTERVAL", [fnum(s.TU / 2)]),
        ]


# ----------------------------------------------------------------------------------
# config files: a line-based reader/writer for the axes_override block (stdlib only)
# ----------------------------------------------------------------------------------
_AXIS = re.compile(r"^(\s+)([A-Z0-9_]+):\s*\[(.*)\](\s*#.*)?$")


def _split_list(body):
    return [v.strip() for v in body.split(",") if v.strip()]


def _num(v):
    try:
        return float(v.strip('"').strip("'"))
    except ValueError:
        return None


def read_axes(path):
    """{KEY: (lineno, [values])} for every list-valued key of axes_override."""
    axes = {}
    inside = False
    for i, line in enumerate(open(path)):
        if line.startswith("axes_override:"):
            inside = True
            continue
        if inside and line and not line[0].isspace() and line.strip():
            inside = False
        m = _AXIS.match(line) if inside else None
        if m:
            axes[m.group(2)] = (i, _split_list(m.group(3)))
    return axes


def _close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(abs(a), abs(b), 1e-300)


def check(path, s, out=sys.stdout):
    """Verify the dimensional axes of one config against the set. Returns [problems]."""
    axes = read_axes(path)
    probs = []

    def one(key, expect, allow=()):
        if key not in axes:
            return
        for v in axes[key][1]:
            x = _num(v)
            if x is None:
                probs.append(f"{key}: non-numeric {v}")
            elif not (_close(x, expect) or any(_close(x, a) for a in allow)):
                probs.append(f"{key}: {v} (SI set: {fnum(expect)})")

    one("DOMAIN_LENGTH", s.H)
    one("CAPILLARY_REF_LENGTH", s.H)
    one("CAPILLARY_DT_COEFF", s.coeff)
    one("DROPLET_RADIUS", s.R)
    one("TRANSLATION_SPEED", s.U)
    one("SIGMA", s.sigma, allow=(0.0,))          # SIGMA 0 = the passive-transport control
    one("POPINET_RHO", s.rho)
    one("POPINET_XLEN", s.aspect)
    if "POPINET_NU" in axes:
        for v in axes["POPINET_NU"][1]:
            x = _num(v)
            if x is None:
                probs.append(f"POPINET_NU: non-numeric {v}")
            elif x == 0.0:
                continue                           # inviscid arm of the La sweep
            else:
                La = s.sigma * s.D / (s.rho * x ** 2)
                if not any(_close(La, ref, 1e-4) for ref in (120.0, 1200.0, 12000.0)):
                    probs.append(f"POPINET_NU: {v} gives La = {La:.4g}, not one of Popinet's (120, 1200, 12000)")
    if "END_TIME" in axes:
        for v in axes["END_TIME"][1]:
            x = _num(v)
            if x is None or x > s.TU * (1 + 1e-6):
                probs.append(f"END_TIME: {v} exceeds T_U = {fnum(s.TU)} (one diameter of travel)")
            elif x < 1e-3 * s.TU:
                probs.append(f"END_TIME: {v} is below 1e-3 T_U -- a dimensionless left-over?")
    if "MAX_CELL_SIZE" in axes and "N_CELLS" in axes:
        n = _num(axes["N_CELLS"][1][0])
        for v in axes["MAX_CELL_SIZE"][1]:
            x = _num(v)
            if x is None or n is None or n <= 0:
                continue
            h_if = x * 2 ** (-1.0 / 3.0)           # cfMesh interface cell = 2^(-1/3) maxCellSize
            if abs(h_if - s.H / n) / (s.H / n) > 0.05:
                probs.append(f"MAX_CELL_SIZE: {v} gives an interface cell {h_if:.4e} m against "
                             f"H/N_CELLS = {s.H / n:.4e} m (> 5 % off: the capillary dt pin is wrong)")
    if "SURFACE_FILE" in axes:
        for v in axes["SURFACE_FILE"][1]:
            if v in SURFACE_RENAME:
                probs.append(f"SURFACE_FILE: {v} is the dimensionless box; use {SURFACE_RENAME[v]}")
    out.write(f"{'OK  ' if not probs else 'FAIL'} {path}\n")
    for p in probs:
        out.write(f"      {p}\n")
    return probs


def rewrite(path, s, out=sys.stdout):
    """Convert one config from the dimensionless set to SI, in place."""
    lines = open(path).read().split("\n")
    axes = read_axes(path)
    if "DOMAIN_LENGTH" not in axes:
        out.write(f"skip {path}: no DOMAIN_LENGTH axis\n")
        return False
    dl = _num(axes["DOMAIN_LENGTH"][1][0])
    if _close(dl, s.H):
        out.write(f"skip {path}: already SI (DOMAIN_LENGTH = {fnum(s.H)})\n")
        return False
    if not _close(dl, ND["H"]):
        raise SystemExit(f"{path}: DOMAIN_LENGTH = {dl} is neither the dimensionless 1 nor the SI {fnum(s.H)}")

    def scale(factor):
        return lambda v: fnum(_num(v) * factor)

    def coeff(v):
        x = _num(v)
        if not _close(x, 0.13107, 1e-3):
            out.write(f"  note {path}: CAPILLARY_DT_COEFF was {v}, not 0.13107; scaled by time\n")
            return fnum(x * s.lamT)
        return fnum(s.coeff)

    def sigma(v):
        x = _num(v)
        return "0" if x == 0 else fnum(x * s.lamSigma)

    def surface(v):
        return SURFACE_RENAME.get(v, v)

    rules = {
        "DOMAIN_LENGTH": scale(s.lamL), "CAPILLARY_REF_LENGTH": scale(s.lamL),
        "DROPLET_RADIUS": scale(s.lamL), "DROPLET_OFFSET_X": scale(s.lamL),
        "MAX_CELL_SIZE": scale(s.lamL),
        "TRANSLATION_SPEED": scale(s.lamU),
        "POPINET_NU": scale(s.lamNu), "POPINET_RHO": scale(s.lamRho),
        "SIGMA": sigma, "CAPILLARY_DT_COEFF": coeff,
        "END_TIME": scale(s.lamT), "WRITE_INTERVAL": scale(s.lamT), "FIXED_DELTA_T": scale(s.lamT),
        "SURFACE_FILE": surface,
    }
    changed = 0
    for key, (i, vals) in axes.items():
        if key not in rules:
            continue
        m = _AXIS.match(lines[i])
        new = [rules[key](v) for v in vals]
        if new != vals:
            lines[i] = f"{m.group(1)}{key}: [{', '.join(new)}]{m.group(4) or ''}"
            changed += 1
    open(path, "w").write("\n".join(lines))
    out.write(f"rewrote {path}: {changed} axes\n")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=("print", "yaml", "check", "rewrite"))
    ap.add_argument("configs", nargs="*")
    ap.add_argument("--D", type=float, default=1.0e-3)
    ap.add_argument("--rho", type=float, default=1000.0)
    ap.add_argument("--nu", type=float, default=1.0e-6)
    ap.add_argument("--La", type=float, default=12000.0)
    ap.add_argument("--We", type=float, default=0.4)
    ap.add_argument("--safety", type=float, default=0.2323)
    ap.add_argument("--aspect", type=float, default=2.0)
    a = ap.parse_args()
    s = SISet(a.D, a.rho, a.nu, a.La, a.We, a.safety, a.aspect)
    if a.mode == "print":
        s.print()
    elif a.mode == "yaml":
        for k, v in s.axes():
            print(f"  {k}: [{', '.join(v)}]")
    elif a.mode == "check":
        if not a.configs:
            sys.exit("check: give at least one config")
        bad = sum(1 for c in a.configs if check(c, s))
        sys.exit(1 if bad else 0)
    elif a.mode == "rewrite":
        if not a.configs:
            sys.exit("rewrite: give at least one config")
        for c in a.configs:
            rewrite(c, s)


if __name__ == "__main__":
    main()
