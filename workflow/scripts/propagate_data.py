#!/usr/bin/env python3
"""Sync theme article data/ -> presentation data/ for dual-copy themes.

The semi-Lagrangian method lines keep a verbatim mirror of the article's
data/{figures,tables} inside the presentation folder, so the reveal deck and
the Elsevier article always consume identical, current results. The snakemake
report rule calls paths.propagate() after every study; this CLI covers manual
syncs and docs/build-decks.sh (self-healing before a deck build).

Usage:
    python3 workflow/scripts/propagate_data.py [theme ...]

With no arguments every dual-copy theme is synced.
"""
import sys

import paths


def main(argv):
    themes = argv or sorted(paths._DUAL_COPY)
    for theme in themes:
        paths.propagate(theme)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
