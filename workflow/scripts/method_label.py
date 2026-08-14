#!/usr/bin/env python3
"""The single definition of a method label, in its three renderings.

There used to be two divergent producers: the ``method_of()`` closure inside
``aggregate._write_error_table`` (curated CSV, figure legends, LaTeX tables) and
``make_bench_fields_fig.method_of()`` (PNG filenames). They disagreed about both
separators and content -- the filename one ignores ``SOURCE_SCHEME``, so the two
SDPLS linearization arms of a study silently overwrote each other's montage.

Three renderings, one source of truth:

``method_label(rec)``
    FULLY PRECISE. Identifies the method *and* the discretization it ran with:
    the level-set div scheme, the gradient that scheme reconstructs with, and
    the defect-correction count. This is the ``method`` column of the curated
    CSV and the cell text of the generated LaTeX tables. Restricted to
    LaTeX-safe characters -- ``A-Za-z0-9 : / + ( ) - .`` and never any of
    ``_ ^ & % # $ { } \\ ~ |`` -- because several emitters write it into a
    tabular cell unescaped.

``method_slug(rec)``
    FILENAME-SAFE and deliberately compact: ``[A-Za-z0-9_]`` only. Used for
    per-method figure filenames, which are referenced by hand from the article
    and the deck, so it stays as close to the historical form as possible while
    still disambiguating arms that would otherwise collide.

``latex_escape(s)``
    For any other string headed into a tabular cell.

Why the discretization belongs *in* the label: a 3D SDPLS study run with
``div(phi,psi) Gauss linear`` and a 2D one run with ``Gauss linearUpwind
grad(psi)`` are not the same method, and comparing them cost a full study round
because nothing in the curated CSV recorded the difference.
"""
import re

__all__ = ["method_label", "method_slug", "latex_escape"]

_LIN_TAG = {
    "explicit": "expl",
    "simpleLinearImplicit": "simpleImp",
    "strictNegativeSpLinearImplicit": "strictNegSp",
}

# LaTeX specials that break a tabular cell.
_LATEX = {"&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
          "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
          "^": r"\textasciicircum{}", "\\": r"\textbackslash{}"}


def latex_escape(s):
    return "".join(_LATEX.get(c, c) for c in str(s))


def _is_semi_lagrangian(rec):
    return rec.get("ADVECTION", "eulerian") == "semiLagrangian"


def _sl_parts(rec):
    """Semi-Lagrangian label. UNCHANGED from the historical form: the SL update
    is a pointwise assignment with no div(phi,psi) term at all, so none of the
    discretization components below apply, and the cross-method comparison
    figures keep matching these strings."""
    label = f"SL:{rec.get('SL_RECONSTRUCTION', '')}"
    scheme = rec.get("SL_SCHEME", "pointValue")
    if scheme == "fluxForm":
        label += "+flux"
    elif scheme == "normalProjected":
        label += "+nSL"
    lim = rec.get("SL_LIMITER", "none")
    if lim not in ("", "none"):
        label += "+lim:" + ("venk" if lim == "venkatakrishnan" else lim)
    if rec.get("SL_FIT", "normalEquations") == "householderQR":
        label += "+qr"
    if rec.get("mesh") == "perturbed":
        label += "/pert"
    return label


# The historical reconstruction gradient for div(phi,psi). Anything else is an
# ABLATION arm and must be visible in the label, or the variants collapse to one
# `method` string at each h and the order fit runs straight through them -- the
# same defect the beta sweep had before _beta() existed.
_DEFAULT_DIV_GRAD_SCHEME = "cellLimited leastSquares 1"
_GRAD_TAG = {
    "leastSquares":               "lsq",
    "cellLimited Gauss linear 1": "clGauss",
    "faceLimited leastSquares 1": "flLsq",
    "Gauss linear":               "gauss",
}


def _compact_div(rec):
    """`Gauss linearUpwind grad(psi)` -> `linearUpwind/grad(psi)`.

    The named gradient is load-bearing, not decoration: `linearUpwind grad(psi)`
    means something different depending on what `grad(psi)` resolves to, and in
    a case with no `grad(psi)` key it silently falls back to gradSchemes
    `default`. When that resolution is not the historical default, the RESOLVED
    scheme is appended too -- omitted otherwise, so every pre-ablation label
    stays byte-identical."""
    div = (rec.get("divPsi") or "").strip()
    if not div:
        return ""
    toks = div.split()
    if toks and toks[0] == "Gauss":
        toks = toks[1:]
    if not toks:
        return ""
    scheme = toks[0]
    named = (rec.get("divPsiGrad") or "").strip()
    out = f"{scheme}/{named}" if named else scheme
    resolved = " ".join((rec.get("divPsiGradScheme") or "").split())
    if resolved and resolved != _DEFAULT_DIV_GRAD_SCHEME:
        out += ":" + _GRAD_TAG.get(resolved, resolved.replace(" ", "-"))
    return out


def _dc(rec):
    """Defect-correction passes; `3 (solver default)` -> `3`."""
    m = re.match(r"\s*(\d+)", str(rec.get("nDefCorr", "")))
    return m.group(1) if m else ""


def _beta(rec):
    """The sdplsBeta target ``beta``, rendered ONLY when it differs from the
    default 1.0.

    Same idiom as every other component in this module -- ``VELOCITY_EXTENSION
    none`` and ``REDISTANCER noRedistancing`` are likewise omitted -- so every
    label produced before a beta sweep existed stays byte-identical, and only a
    genuinely off-default arm introduces a new string.

    Without this the arms of a beta sweep collapse to ONE ``method`` at one
    ``h``, and ``make_convergence_table.py`` fits a single regression straight
    through them and reports the slope as a convergence order.
    """
    if rec.get("SDPLS_SOURCE", "noSource") != "beta":
        return ""
    try:
        v = float(rec.get("SDPLS_BETA", 1.0))
    except (TypeError, ValueError):
        return ""
    return "" if v == 1.0 else f"{v:g}"


def _mollifier(rec):
    """The SDPLS cut-off, rendered ONLY when one is active.

    ``mollifier1`` (dictionary name ``m1``) is eq. (24) of arXiv:2208.01269: it
    forces the source to zero away from the interface. Turning it on changes the
    equation being solved, so a mollified arm is a DIFFERENT method and must not
    share a label with an unmollified one -- otherwise a study that sweeps it
    collapses both arms into one series per h and the fitted slope is a
    regression through two different methods. Exactly the failure `_beta` and
    the reconstruction-gradient tag already guard against.

    Widths are part of the identity: they are physical lengths that set where
    the source stops acting, not a solver tolerance.
    """
    m = (rec.get("MOLLIFIER") or "none").strip()
    if not m or m == "none":
        return ""
    # The `band` cut-off is parameterised by a CELL COUNT, not by lengths: its
    # w1/w2 are inert, and rendering them would make every nLayers arm look
    # identical -- the collision this function exists to prevent.
    if m == "band":
        n = (rec.get("SDPLS_NLAYERS") or "").strip()
        return f"{m}({n})" if n else m
    w1 = (rec.get("SDPLS_W1") or "").strip()
    w2 = (rec.get("SDPLS_W2") or "").strip()
    return f"{m}({w1},{w2})" if (w1 and w2) else m


def _euler_parts(rec):
    parts = ["euler"]
    ve = rec.get("VELOCITY_EXTENSION", "none")
    if ve and ve != "none":
        parts.append(f"VE:{ve}")
    ss = rec.get("SDPLS_SOURCE", "noSource")
    if ss and ss != "noSource":
        # The linearization is part of the method: the three discretizations of
        # the same continuum source are separate arms, and they disagreed badly
        # until the 2026-08 sign fix.
        raw = rec.get("SOURCE_SCHEME", "")
        tag = _LIN_TAG.get(raw, raw)
        b = _beta(rec)
        name = f"{ss}({b})" if b else ss
        parts.append(f"SDPLS:{name}/{tag}" if tag else f"SDPLS:{name}")
        moll = _mollifier(rec)
        if moll:
            parts.append(f"moll:{moll}")
    rd = rec.get("REDISTANCER", "noRedistancing")
    if rd and rd != "noRedistancing":
        if rd == "PDE" and rec.get("REDIST_FREEZE", "false") == "true":
            rd = "PDEfrozen"   # bulk-only fill, frozen Dirichlet anchors
        parts.append(f"RD:{rd}")
    return parts


def method_label(rec):
    """Fully precise, LaTeX-safe method + discretization identifier."""
    if _is_semi_lagrangian(rec):
        return _sl_parts(rec)
    parts = _euler_parts(rec)
    div = _compact_div(rec)
    if div:
        parts.append(f"div:{div}")
    dc = _dc(rec)
    if dc:
        parts.append(f"dc:{dc}")
    return "+".join(parts)


def method_slug(rec):
    """Filename-safe, compact. Historical form plus the SDPLS linearization.

    Deliberately does NOT carry the div/dc components: these names are
    referenced by hand from methodComparison.tex and the comparison deck, and
    the full provenance is one column away in the curated CSV. It DOES carry the
    linearization, because without it the two SDPLS arms of a study write to the
    same PNG and the second silently overwrites the first.
    """
    if _is_semi_lagrangian(rec):
        base = f"SL_{rec.get('SL_RECONSTRUCTION', '')}"
        scheme = rec.get("SL_SCHEME", "pointValue")
        if scheme == "fluxForm":
            base += "_flux"
        elif scheme == "normalProjected":
            base += "_nSL"
        return re.sub(r"[^A-Za-z0-9_]+", "_", base).strip("_")
    parts = ["euler"]
    ve = rec.get("VELOCITY_EXTENSION", "none")
    if ve and ve != "none":
        parts.append(f"VE_{ve}")
    ss = rec.get("SDPLS_SOURCE", "noSource")
    if ss and ss != "noSource":
        raw = rec.get("SOURCE_SCHEME", "")
        tag = _LIN_TAG.get(raw, raw)
        b = _beta(rec)
        name = f"{ss}_{b}" if b else ss
        parts.append(f"SDPLS_{name}_{tag}" if tag else f"SDPLS_{name}")
    rd = rec.get("REDISTANCER", "noRedistancing")
    if rd and rd != "noRedistancing":
        if rd == "PDE" and rec.get("REDIST_FREEZE", "false") == "true":
            rd = "PDEfrozen"
        parts.append(f"RD_{rd}")
    # OFF-DEFAULT reconstruction gradient only, exactly as _beta() renders only
    # an off-default beta. sdplsOrderAblation sweeps this, and without it all
    # four variants of an arm collapse to ONE filename -- make_bench_fields_fig
    # then picks among genuinely different discretizations by glob order and
    # labels the atlas for a configuration it is not. Suppressing the default
    # keeps every pre-ablation filename byte-identical, which matters because
    # these names are referenced by hand from methodComparison.tex and the decks.
    resolved = " ".join((rec.get("divPsiGradScheme") or "").split())
    if resolved and resolved != _DEFAULT_DIV_GRAD_SCHEME:
        parts.append(_GRAD_TAG.get(resolved, resolved.replace(" ", "-")))
    # Same reason as in method_label: a mollified arm is a different method, and
    # two arms sharing a filename means one montage silently overwrites the other.
    moll = _mollifier(rec)
    if moll:
        parts.append("moll" + moll.split("(")[0])
    return re.sub(r"[^A-Za-z0-9_]+", "_", "_".join(parts)).strip("_")
