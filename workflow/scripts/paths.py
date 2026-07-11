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
    "semi-lagrangian-level-set": ("semi-lagrangian-level-set", "sl-level-set",       "sl.template.html"),
    "velocity-extension":        ("velocity-extension",        "velocity-extension", "index.template.html"),
}

# Default theme for a solver when a study config does not set `theme:` explicitly.
_SOLVER_THEME = {
    "leiaSemiLagrangeLevelSetFoam": "semi-lagrangian-level-set",
    "leiaLevelSetFoam":             "velocity-extension",
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
