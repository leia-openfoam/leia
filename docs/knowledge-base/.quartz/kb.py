#!/usr/bin/env python3
"""Query the leia knowledge base: the permanent memory of every agent and every session.

Per-machine memories (a Claude session's own memory, a terminal history) are not shared; this
vault is, through git. Query it before you propose a route, and write to it when you decide,
retract or fail (docs/knowledge-base/conventions.md).

Usage, from the repository root (standard library only):
  kb.py search <words...> [--kind K] [--status S] [--negative] [-n N]
        notes ranked by the words in title, aliases and description (weight 3) and in the body
        (weight 1); prints path, kind, status, outcome, date, title and description
  kb.py show <note>              the frontmatter, the verdict and the log of one note
                                 (<note> is a vault path such as decisions/psi-filter-none,
                                 a slug such as psi-filter-none, or an alias)
  kb.py negative [words...]      the negative-findings log: every route tried and failed
  kb.py decisions [words...]     the decision log
  kb.py retractions [words...]   the retraction log
  kb.py open [words...]          the notes with status open or candidate: the open questions
  kb.py since <YYYY-MM-DD>       every dated log entry since the date: the three logs and the
                                 '### YYYY-MM-DD' entries of every note's '## Log'
Words match case-insensitively; every word must occur (AND).
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VAULT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import check_kb  # noqa: E402  (the same parser as the checker and the graph)


def load():
    notes, aliases = {}, {}
    for rel, x in check_kb.load_notes(VAULT).items():
        if not x["fm"]:
            continue
        notes[rel] = x
        for a in x["fm"].get("aliases") or []:
            aliases[str(a).lower()] = rel
    return notes, aliases


def text_of(x):
    fm = x["fm"]
    head = " ".join(str(fm.get(k, "")) for k in ("title", "description"))
    head += " " + " ".join(str(a) for a in (fm.get("aliases") or []))
    return head.lower(), x["body"].lower()


def summary(rel, x):
    fm = x["fm"]
    out = fm.get("outcome", "")
    tag = f"{fm.get('kind', '')}/{fm.get('status', '')}" + (f"/{out}" if out else "")
    return f"{rel}  [{tag}, {fm.get('date', '')}]\n    {fm.get('title', '')}\n    {fm.get('description', '')}"


def resolve(notes, aliases, name):
    name = name[:-3] if name.endswith(".md") else name
    if name in notes:
        return name
    hits = [r for r in notes if r.split("/")[-1] == name]
    if len(hits) == 1:
        return hits[0]
    return aliases.get(name.lower())


def cmd_search(args):
    kind = status = None
    negative = False
    n = 15
    words = []
    it = iter(args)
    for a in it:
        if a == "--kind":
            kind = next(it)
        elif a == "--status":
            status = next(it)
        elif a == "--negative":
            negative = True
        elif a == "-n":
            n = int(next(it))
        else:
            words.append(a.lower())
    notes, _ = load()
    ranked = []
    for rel, x in notes.items():
        fm = x["fm"]
        if kind and fm.get("kind") != kind:
            continue
        if status and fm.get("status") != status:
            continue
        if negative and fm.get("outcome") != "negative":
            continue
        head, body = text_of(x)
        if words and not all(w in head or w in body for w in words):
            continue
        score = sum(3 * head.count(w) + body.count(w) for w in words) if words else 0
        ranked.append((score, rel, x))
    ranked.sort(key=lambda r: (-r[0], r[1]))
    for score, rel, x in ranked[:n]:
        print(summary(rel, x))
    print(f"-- {len(ranked)} note(s) match" + (f", the first {n} shown" if len(ranked) > n else ""))


def cmd_show(args):
    notes, aliases = load()
    rel = resolve(notes, aliases, args[0]) if args else None
    if not rel:
        print(f"no note '{args[0] if args else ''}'; try: kb.py search <words>")
        return 1
    x = notes[rel]
    print(summary(rel, x))
    for k in ("decided_by", "date_settled", "code", "sources"):
        if x["fm"].get(k):
            print(f"    {k}: {x['fm'][k]}")
    body = x["body"]
    m = re.search(r"^(>.*(?:\n>.*)*)", body, re.M)
    verdict = m.group(1) if m else next((p for p in body.split("\n\n") if p.strip() and not p.startswith("#")), "")
    print("\nVERDICT\n" + verdict.strip())
    i = body.find("\n## Log")
    if i >= 0:
        print("\nLOG" + body[i + len("\n## Log"):].rstrip())
    print(f"\nfile: {os.path.relpath(x['path'])}")
    return 0


def cmd_log(name, args):
    notes, _ = load()
    body = notes.get(name, {}).get("body", "")
    words = [a.lower() for a in args]
    lines = [l for l in body.splitlines() if l.startswith("- ")]
    lines = [l for l in lines if all(w in l.lower() for w in words)]
    for l in lines:
        print(l)
    print(f"-- {len(lines)} line(s) in {name}.md")


def cmd_open(args):
    words = [a.lower() for a in args]
    notes, _ = load()
    hits = []
    for rel, x in sorted(notes.items()):
        if x["fm"].get("status") in ("open", "candidate") and x["fm"].get("kind") != "index":
            head, body = text_of(x)
            if all(w in head or w in body for w in words):
                hits.append(rel)
                print(summary(rel, x))
    print(f"-- {len(hits)} open or candidate note(s)")


def cmd_since(args):
    if not args or not re.fullmatch(r"\d{4}-\d{2}-\d{2}", args[0]):
        print("usage: kb.py since YYYY-MM-DD")
        return 2
    since = args[0]
    notes, _ = load()
    rows = []
    for log in ("decision-log", "retraction-log", "negative-log"):
        for l in notes.get(log, {}).get("body", "").splitlines():
            m = re.match(r"- (\d{4}-\d{2}-\d{2})", l)
            if m and m.group(1) >= since:
                rows.append((m.group(1), log, l[2:]))
    for rel, x in notes.items():
        body = x["body"]
        i = body.find("\n## Log")
        if i < 0:
            continue
        for m in re.finditer(r"^### (\d{4}-\d{2}-\d{2})\s*\n(.*?)(?=^### |\Z)", body[i:], re.M | re.S):
            if m.group(1) >= since:
                text = " ".join(m.group(2).split())
                rows.append((m.group(1), rel, text[:240]))
    rows.sort()
    for d, where, text in rows:
        print(f"{d}  {where}: {text}")
    print(f"-- {len(rows)} entr(y/ies) since {since}")
    return 0


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print(__doc__)
        return 0
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "search":
        return cmd_search(args)
    if cmd == "show":
        return cmd_show(args)
    if cmd == "negative":
        return cmd_log("negative-log", args)
    if cmd == "decisions":
        return cmd_log("decision-log", args)
    if cmd == "retractions":
        return cmd_log("retraction-log", args)
    if cmd == "open":
        return cmd_open(args)
    if cmd == "since":
        return cmd_since(args)
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main() or 0)
