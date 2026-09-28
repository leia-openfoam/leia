#!/usr/bin/env python3
"""List the notes of the manifest that do not exist yet, grouped by owner.

Usage: python3 missing_notes.py <vault-dir>

The manifest (.quartz/manifest.md) fixes every slug and its owner (A to F, ME). A relaunch of
the note writers starts from this list, never from a re-plan; each writer reads
.quartz/writing-brief.md, writes a digest to .quartz/digests/<slug>.md before drafting, and
writes each note to disk as soon as it is done.
"""
import os
import re
import sys


def main():
    vault = sys.argv[1] if len(sys.argv) > 1 else "docs/knowledge-base"
    man = open(os.path.join(vault, ".quartz", "manifest.md"), encoding="utf-8").read()
    rows = {}
    for m in re.finditer(r"^((?:hubs|models|concepts|cases|studies|decisions|retractions|sessions)/[a-z0-9-]+)\s*\|.*\|\s*([A-Z]+)\s*$", man, re.M):
        rows[m.group(1)] = m.group(2)
    block = man.split("## concepts/ gradient control (owner ME)")[1].split("## Appendix")[0]
    for slug in re.findall(r"concepts/[a-z0-9-]+", block):
        rows.setdefault(slug, "ME")
    by_owner = {}
    for slug, owner in sorted(rows.items()):
        exists = os.path.exists(os.path.join(vault, slug + ".md"))
        by_owner.setdefault(owner, []).append((slug, exists))
    total = missing = 0
    for owner, items in sorted(by_owner.items()):
        miss = [s for s, e in items if not e]
        total += len(items); missing += len(miss)
        print(f"owner {owner}: {len(items) - len(miss)}/{len(items)} written" + (":" if miss else ""))
        for s in miss:
            print(f"  missing  {s}")
    print(f"total: {total - missing}/{total} written, {missing} missing")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
