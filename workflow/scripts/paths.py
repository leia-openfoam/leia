#!/usr/bin/env python3
"""Single source of truth for the thematic docs layout.

Each research line ("theme") owns a folder under ``docs/`` holding a self-contained
article and presentation; the theme's results (figures + CSV tables) live in ONE
``<theme>-article/data/`` folder that BOTH the reveal.js deck and the Elsevier LaTeX
article consume. This module maps a theme (or a solver, via the default) to those
directories, so the Snakefile and the figure scripts never hard-code a path.

    theme  ->  docs/<theme-dir>/<slug>-article/data/{figures,tables}
           ->  docs/<theme-dir>/<slug>-presentation/<deck>.template.html
"""
import os

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# theme key -> (theme directory, subfolder slug, deck template name)
_THEMES = {
    "semi-lagrangian-level-set":          ("semi-lagrangian-level-set",          "sl-level-set",       "quadratic-semi-lagrangian-level-set.template.html"),
    "linear-semi-lagrangian-level-set":   ("linear-semi-lagrangian-level-set",   "lsl-level-set",      "linear-semi-lagrangian-level-set.template.html"),
    "velocity-extension":                 ("velocity-extension",                 "velocity-extension", "velocity-extension.template.html"),
    "geometrically-redistanced-levelset": ("geometrically-redistanced-levelset", "grl-level-set",      "geometrically-redistanced-level-set.template.html"),
    "sdpls-level-set":                    ("sdpls-level-set",                    "sdpls",              "sdpls-level-set.template.html"),
    "method-comparison":                  ("method-comparison",                  "method-comparison",  "level-set-method-comparison.template.html"),
}

# Themes whose deck reads figures/tables from ITS OWN presentation data/ folder
# (a verbatim mirror of the article data/). For these, propagate() copies the
# article data into the presentation data after every report, so a completed
# study leaves BOTH documents current. Deliberate duplication: the two folders
# are kept byte-identical, never edited by hand.
_DUAL_COPY = {"semi-lagrangian-level-set", "linear-semi-lagrangian-level-set"}

# Default theme for a solver when a study config does not set `theme:` explicitly.
_SOLVER_THEME = {
    "leiaSemiLagrangeLevelSetFoam":            "semi-lagrangian-level-set",
    "leiaSemiLagrangianLevelSetTwoPhaseFoam":  "semi-lagrangian-level-set",
    "leiaLevelSetFoam":                        "velocity-extension",
    "leiaRedistancedLevelSetFoam":             "geometrically-redistanced-levelset",
}


def theme_of(solver=None, explicit=None):
    """Resolve the theme: an explicit config value wins, else map by solver."""
    if explicit:
        if explicit not in _THEMES:
            raise KeyError(f"unknown theme {explicit!r}; known: {sorted(_THEMES)}")
        return explicit
    return _SOLVER_THEME.get(solver, "velocity-extension")


def _theme(theme):
    if theme not in _THEMES:
        raise KeyError(f"unknown theme {theme!r}; known: {sorted(_THEMES)}")
    return _THEMES[theme]


def article_dir(theme):
    tdir, slug, _ = _theme(theme)
    return os.path.join(REPO, "docs", tdir, f"{slug}-article")


def data_dir(theme):
    return os.path.join(article_dir(theme), "data")


def figs_dir(theme, make=True):
    d = os.path.join(data_dir(theme), "figures")
    if make:
        os.makedirs(d, exist_ok=True)
    return d


def tables_dir(theme, make=True):
    d = os.path.join(data_dir(theme), "tables")
    if make:
        os.makedirs(d, exist_ok=True)
    return d


def presentation_dir(theme):
    tdir, slug, _ = _theme(theme)
    return os.path.join(REPO, "docs", tdir, f"{slug}-presentation")


def deck_template(theme):
    return _theme(theme)[2]


def pres_data_dir(theme):
    return os.path.join(presentation_dir(theme), "data")


def pres_figs_dir(theme, make=True):
    d = os.path.join(pres_data_dir(theme), "figures")
    if make:
        os.makedirs(d, exist_ok=True)
    return d


def pres_tables_dir(theme, make=True):
    d = os.path.join(pres_data_dir(theme), "tables")
    if make:
        os.makedirs(d, exist_ok=True)
    return d


def propagate(theme):
    """Mirror the theme's article data/ (figures + tables) into its presentation
    data/, so the deck and the article consume identical, current results.
    No-op for themes whose deck references the article data/ directly."""
    import shutil
    if theme not in _DUAL_COPY:
        return []
    copied = []
    for src_dir, dst_dir in ((figs_dir(theme, make=False), pres_figs_dir(theme)),
                             (tables_dir(theme, make=False), pres_tables_dir(theme))):
        if not os.path.isdir(src_dir):
            continue
        for name in sorted(os.listdir(src_dir)):
            src = os.path.join(src_dir, name)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(dst_dir, name))
                copied.append(name)
    print(f"[paths] propagate {theme}: {len(copied)} file(s) -> presentation data/")
    return copied
