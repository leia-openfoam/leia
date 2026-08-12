#!/usr/bin/env python3
"""Read the discretization a case ACTUALLY ran with, from its rendered dicts.

Why this exists
---------------
Study tokens are not a reliable record of the discretization:

* ``cases/2Dvortex`` hardcodes ``div(phi,psi) Gauss linearUpwind grad(psi)`` and
  has no ``@!DIV!@`` token at all, so ``DIV`` never reaches ``case_params.json``
  for any ``benchVortex*`` study -- the value in ``cases/2Dvortex.parameter`` is
  dropped as an unreferenced token.
* Where a ``DIV`` token does exist it may hold an unexpanded OpenFOAM alias
  (``$LinearUpwind``), whose meaning lives in the case's ``aliasDict`` and whose
  inner gradient name cannot be read off the token.
* ``gradPsiSdpls`` / ``gradUSdpls`` -- the SDPLS source's own gradient schemes --
  are hardcoded in every template and are not tokens at all.

This cost a whole 3D study round: the 2D and 3D SDPLS arms silently differed in
``div(phi,psi)`` (``linearUpwind grad(psi)`` vs ``linear``), which is a different
method, and nothing in the curated CSV recorded it.

``materialize.py`` renders every ``*.template`` into the case directory and then
deletes the template, so ``<case>/system/fvSchemes`` and ``<case>/system/controlDict``
are the ground truth for what the solver read. Parse those.
"""
import os
import re

__all__ = ["read_discretization", "COLUMNS"]

#: Column order for the curated error CSV.
COLUMNS = ["divPsi", "divPsiGrad", "divPsiGradScheme",
           "gradPsiSdpls", "gradUSdpls", "nDefCorr",
           # Whether the velocity ACTUALLY reversed -- see read_discretization.
           "oscillationRendered"]

_COMMENT_BLOCK = re.compile(r"/\*.*?\*/", re.S)
_COMMENT_LINE = re.compile(r"//[^\n]*")


def _strip_comments(text):
    return _COMMENT_LINE.sub("", _COMMENT_BLOCK.sub("", text))


def _block(text, name):
    """Body of a top-level ``name { ... }`` block, or '' if absent."""
    m = re.search(r"^\s*" + re.escape(name) + r"\s*\{", text, re.M)
    if not m:
        return ""
    depth, start = 1, m.end()
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start:i]
    return text[start:]


def _entries(body):
    """``{key: value}`` for ``key  value;`` lines. Keys may contain parens,
    e.g. ``div(phi,psi)`` or ``grad(psi)``."""
    out = {}
    for m in re.finditer(r"([A-Za-z_][\w()\[\],.]*)\s+([^;{}]+?)\s*;", body):
        out[m.group(1)] = " ".join(m.group(2).split())
    return out


def _aliases(case_dir):
    """``aliasDict`` entries. OpenFOAM expands ``$Alias`` at run time via the
    ``#include "aliasDict"`` in fvSchemes; the workflow keeps them verbatim, so
    resolve here to record what the value actually means."""
    path = os.path.join(case_dir, "system", "aliasDict")
    if not os.path.isfile(path):
        return {}
    try:
        return _entries(_strip_comments(open(path).read()))
    except OSError:
        return {}


def _resolve(value, aliases):
    """Expand a leading ``$Alias`` token (possibly repeatedly, bounded)."""
    for _ in range(4):
        m = re.match(r"^\$(\w+)$", value.strip())
        if not m or m.group(1) not in aliases:
            return value
        value = aliases[m.group(1)]
    return value


def read_discretization(case_dir):
    """Discretization provenance for one materialized case directory.

    Never raises: a missing or unparsable dictionary yields blanks, so a partial
    study still aggregates. Blank means "not recorded", never "zero"/"none".
    """
    out = {k: "" for k in COLUMNS}

    fvs = os.path.join(case_dir, "system", "fvSchemes")
    if os.path.isfile(fvs):
        try:
            text = _strip_comments(open(fvs).read())
        except OSError:
            text = ""
        if text:
            aliases = _aliases(case_dir)
            grads = _entries(_block(text, "gradSchemes"))
            divs = _entries(_block(text, "divSchemes"))

            out["gradPsiSdpls"] = _resolve(grads.get("gradPsiSdpls", ""), aliases)
            out["gradUSdpls"] = _resolve(grads.get("gradUSdpls", ""), aliases)

            div_psi = _resolve(divs.get("div(phi,psi)", ""), aliases)
            out["divPsi"] = div_psi

            # A deferred scheme names the gradient it reconstructs with, e.g.
            # `Gauss linearUpwind grad(psi)`. Record BOTH the name and the
            # scheme that name resolves to -- `linearUpwind grad(psi)` means
            # something different depending on what `grad(psi)` is set to, and
            # in a case with no `grad(psi)` key it silently falls back to
            # gradSchemes `default`.
            tokens = div_psi.split()
            named = next((t for t in tokens[2:] if t in grads
                          or re.match(r"^\$?\w*[Gg]rad", t)), "")
            if named:
                out["divPsiGrad"] = named
                scheme = grads.get(named)
                if scheme is None and "default" in grads:
                    scheme = grads["default"] + "   [via gradSchemes default]"
                out["divPsiGradScheme"] = _resolve(scheme or "", aliases)

    ctrl = os.path.join(case_dir, "system", "controlDict")
    if os.path.isfile(ctrl):
        try:
            ent = _entries(_strip_comments(open(ctrl).read()))
        except OSError:
            ent = {}
        # Absent key => the C++ default in eulerianAdvection.C /
        # leiaRedistancedLevelSetFoam.C, which is 3.
        out["nDefCorr"] = ent.get("nDefCorr", "3 (solver default)")

    # DID THE FLOW ACTUALLY REVERSE? Same reason as everything above: the token
    # is not the truth. cases/3Dshear and cases/3Ddeformation HARDCODE
    # `oscillation on` in their fvSolution.template and have no @!OSCILLATION!@
    # placeholder, so setting OSCILLATION in those studies' configs is a silent
    # no-op -- materialize drops it as an unreferenced token, it never reaches
    # case_params.json, and the curated `oscillation` column comes out blank.
    #
    # That blank is not harmless. The reversed t=T gradient error is suppressed
    # on rows recording oscillation=on (make_convergence_table.REVERSED_INVALID);
    # a blank sails straight through and gets published as a convergence order.
    # Read the rendered dictionary instead, which is what the solver parsed.
    fvsol = os.path.join(case_dir, "system", "fvSolution")
    if os.path.isfile(fvsol):
        try:
            body = _strip_comments(open(fvsol).read())
        except OSError:
            body = ""
        m = re.search(r"^\s*oscillation\s+([^;]+);", body, re.M)
        if m:
            out["oscillationRendered"] = m.group(1).strip()

    return out


if __name__ == "__main__":
    import sys
    for d in sys.argv[1:]:
        print(d)
        for k, v in read_discretization(d).items():
            print(f"    {k:20s} {v}")
