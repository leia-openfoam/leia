#!/usr/bin/env python3
"""Write graph3d/graph.json for the three.js page from the vault's notes and wikilinks.

Usage: python3 build_graph.py <vault-dir>

Nodes: every note whose kind is not "index", plus the root index (Home). Links: every resolved
wikilink between two nodes, deduplicated, without self-links. An unresolved link is an error of
the vault: the script exits 1, so the site build fails on a broken link. Uses the parser of
check_kb.py (standard library only).
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_kb import load_notes, links_of  # noqa: E402


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    vault = sys.argv[1]
    notes = load_notes(vault)
    aliases = {}
    for rel, n in notes.items():
        for a in (n["fm"] or {}).get("aliases", []) or []:
            aliases.setdefault(a, rel)
    keep = {rel for rel, n in notes.items()
            if n["fm"] and (n["fm"].get("kind") != "index" or rel == "index")}
    nodes, links, errors = {}, set(), []
    for rel in sorted(keep):
        fm = notes[rel]["fm"]
        nodes[rel] = {
            "id": rel,
            "title": "Home" if rel == "index" else fm.get("title", rel),
            "kind": "hub" if rel == "index" else fm.get("kind"),
            "part": fm.get("part", "all"),
            "status": fm.get("status", "open"),
            "tags": fm.get("tags", []),
            "description": fm.get("description", ""),
            "url": rel,                 # Quartz URL relative to the site root
            "file": rel + ".md",        # the raw note, for the local preview
            "degree": 0,
        }
    for rel in sorted(keep):
        for target, _heading, is_embed in links_of(notes[rel]["body"]):
            if is_embed:
                continue
            t = target[:-3] if target.endswith(".md") else target
            r = t if t in notes else aliases.get(t)
            if r is None:
                errors.append(f"{rel}.md: link [[{target}]] does not resolve")
                continue
            if r in keep and r != rel:
                links.add((rel, r))
    for s, t in links:
        nodes[s]["degree"] += 1
        nodes[t]["degree"] += 1
    if errors:
        print("\n".join(errors))
        return 1
    out = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "counts": {"nodes": len(nodes), "links": len(links)},
        "parts": sorted({n["part"] for n in nodes.values()}),
        "statuses": sorted({n["status"] for n in nodes.values()}),
        "nodes": list(nodes.values()),
        "links": [{"source": s, "target": t} for s, t in sorted(links)],
    }
    os.makedirs(os.path.join(vault, "graph3d"), exist_ok=True)
    path = os.path.join(vault, "graph3d", "graph.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(f"build_graph: {len(nodes)} nodes, {len(links)} links -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
