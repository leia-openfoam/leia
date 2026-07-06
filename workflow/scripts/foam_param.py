#!/usr/bin/env python3
"""Parse leia OpenFOAM ``.parameter`` study files and build the token grid.

A ``.parameter`` file is a small OpenFOAM dictionary.  Two shapes occur:

  * per-case files (e.g. ``cases/3Drotation.parameter``) wrap the sweep in a
    ``values { ... }`` block, each entry a list:  ``N_CELLS ( 32 64 128 );``
  * ``cases/default.parameter`` is flat:          ``CFL  0.3;``

Values may carry an alias token such as ``$LinearUpwind``.  Aliases are kept
*verbatim* (leading ``$`` included): OpenFOAM expands them at run time via the
``#include "aliasDict"`` in ``fvSchemes``, so this module never expands them.

The grid that actually matters for a study is the set of ``@!TOKEN!@``
placeholders that appear in the case's ``*.template`` files.  Parameters that
are not referenced by any template are ignored (sweeping them would only
produce identical duplicate cases), and referenced tokens that the case does
not vary fall back to ``cases/default.parameter``.

Stand-alone use (validation)::

    python3 foam_param.py cases/3Drotation.parameter
"""
import itertools
import os
import re
import sys
from collections import OrderedDict

_COMMENT_BLOCK = re.compile(r"/\*.*?\*/", re.DOTALL)
_COMMENT_LINE = re.compile(r"//[^\n]*")
_LIST_ENTRY = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)\s*\((.*?)\)\s*;", re.DOTALL)
_SCALAR_ENTRY = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)\s+([^();{}]+?)\s*;")
_TOKEN = re.compile(r"@!([A-Z0-9_]+)!@")

# Tokens computed by materialize._with_derived_tokens from another token (so they
# must NOT be a sweep axis or require a .parameter value); skip them in the grid.
_DERIVED_TOKENS = {"HALF_END_TIME"}


def _strip_comments(text):
    text = _COMMENT_BLOCK.sub(" ", text)
    text = _COMMENT_LINE.sub("", text)
    return text


def _values_block(text):
    """Return the body of the outermost ``values { ... }`` block, or the whole
    text when there is no such block (flat ``default.parameter``)."""
    m = re.search(r"\bvalues\b\s*\{", text)
    if not m:
        return text
    depth, i = 1, m.end()
    while i < len(text) and depth:
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
        i += 1
    return text[m.end():i - 1]


def parse_parameter_file(path):
    """Return an ``OrderedDict {KEY: [value, ...]}`` preserving file order."""
    with open(path) as fh:
        body = _values_block(_strip_comments(fh.read()))
    out = OrderedDict()
    spans = []
    for m in _LIST_ENTRY.finditer(body):
        out[m.group(1)] = m.group(2).split()
        spans.append((m.start(), m.end()))
    for m in _SCALAR_ENTRY.finditer(body):
        if any(s <= m.start() < e for s, e in spans):
            continue  # tail of a list entry
        key, val = m.group(1), m.group(2).strip()
        if val and key not in out:
            out[key] = [val]
    return out


def referenced_tokens(case_dir):
    """Set of ``@!TOKEN!@`` names appearing in any ``*.template`` under case_dir."""
    toks = set()
    for root, _dirs, files in os.walk(case_dir):
        for f in files:
            if f.endswith(".template"):
                with open(os.path.join(root, f)) as fh:
                    toks.update(_TOKEN.findall(fh.read()))
    return toks


def build_token_grid(param_file, default_file, case_dir,
                     axes_override=None, collapse_other_axes=False):
    """Resolve every template token to a value list.

    Returns ``(axes, constants, ignored)`` where *axes* is an OrderedDict of
    tokens with more than one value (sorted by name for deterministic indexing),
    *constants* maps single-valued tokens to their value, and *ignored* lists
    parameters present in the ``.parameter`` file that no template references.
    """
    params = parse_parameter_file(param_file)
    defaults = parse_parameter_file(default_file) if os.path.isfile(default_file) else {}
    referenced = referenced_tokens(case_dir)
    axes_override = {k: [str(x) for x in v] for k, v in (axes_override or {}).items()}

    axes, constants, missing = OrderedDict(), OrderedDict(), []
    for tok in sorted(referenced):
        if tok in _DERIVED_TOKENS:
            continue  # supplied per-case by materialize._with_derived_tokens
        if tok in axes_override and axes_override[tok]:
            vals = axes_override[tok]
        elif tok in params:
            vals = list(params[tok])
        elif tok in defaults:
            vals = list(defaults[tok])
        else:
            missing.append(tok)
            continue
        if len(vals) > 1 and collapse_other_axes and tok not in axes_override:
            vals = vals[:1]
        if len(vals) > 1:
            axes[tok] = vals
        else:
            constants[tok] = vals[0]
    if missing:
        raise SystemExit(
            "[foam_param] template tokens with no value in "
            f"{os.path.basename(param_file)} or default.parameter: "
            + ", ".join(missing))
    ignored = [k for k in params if k != "solver" and k not in referenced]
    return axes, constants, ignored


def expand(axes, constants):
    """Cartesian product of *axes*, merged with *constants*.

    Returns a list of OrderedDicts, one per variation, in a deterministic order.
    """
    keys = list(axes.keys())
    combos = list(itertools.product(*[axes[k] for k in keys])) if keys else [()]
    variations = []
    for combo in combos:
        d = OrderedDict(constants)
        d.update(zip(keys, combo))
        variations.append(d)
    return variations


def _main(argv):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("param_file")
    ap.add_argument("--case-dir", help="case template dir (default: inferred)")
    ap.add_argument("--default", help="default.parameter path (default: inferred)")
    args = ap.parse_args(argv)

    casedir = args.case_dir
    if not casedir:
        base = os.path.basename(args.param_file)
        base = base[:-len(".parameter")] if base.endswith(".parameter") else base
        if base.endswith("_poly"):
            base = base[:-len("_poly")]
        casedir = os.path.join(os.path.dirname(args.param_file), base)
    default_file = args.default or os.path.join(
        os.path.dirname(args.param_file), "default.parameter")

    axes, constants, ignored = build_token_grid(args.param_file, default_file, casedir)
    total = 1
    for v in axes.values():
        total *= len(v)
    print(f"param file : {args.param_file}")
    print(f"case dir   : {casedir}")
    print(f"axes ({len(axes)}):")
    for k, v in axes.items():
        print(f"  {k:14s} [{len(v):2d}] {v}")
    print("constants:")
    for k, v in constants.items():
        print(f"  {k:14s} = {v}")
    if ignored:
        print(f"ignored (not referenced by templates): {ignored}")
    print(f"TOTAL VARIATIONS = {total}")


if __name__ == "__main__":
    _main(sys.argv[1:])
