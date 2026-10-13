#!/usr/bin/env python3
"""What the lab notebook records and the knowledge base does not: the curation backlog.

Usage, from the repository root:  python3 docs/knowledge-base/.quartz/kb_coverage.py [--strict]

Reports, each a way the vault falls behind the record:
 1. STATUS.md sections ('### N.M title') that no note cites as 'STATUS N.M' (in the body, the
    sources or a link text). A range 'STATUS 11.24 to 11.25' cites both ends.
 2. METHOD.md section 8.1 tokens (the rows '| `TOKEN` | ...') that no decision note names in its
    title, aliases, description or body.
A report, not a gate: exit 0, or 1 with --strict when anything is uncovered. Standard library only.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VAULT = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(VAULT))
sys.path.insert(0, HERE)
import check_kb  # noqa: E402


def cited_sections(notes):
    cited = set()
    for x in notes.values():
        fm = x["fm"] or {}
        text = x["body"] + " " + " ".join(str(s) for s in (fm.get("sources") or []))
        for m in re.finditer(r"STATUS\s+(\d+(?:\.\d+)?)(?:\s*(?:to|-|–)\s*(\d+(?:\.\d+)?))?", text):
            cited.add(m.group(1))
            if m.group(2):
                a, b = m.group(1), m.group(2)
                cited.add(b)
                if "." in a and "." in b and a.split(".")[0] == b.split(".")[0]:
                    major = a.split(".")[0]
                    for k in range(int(a.split(".")[1]), int(b.split(".")[1]) + 1):
                        cited.add(f"{major}.{k}")
    return cited


def main():
    strict = "--strict" in sys.argv
    notes = {r: x for r, x in check_kb.load_notes(VAULT).items() if x["fm"]}
    problems = 0

    status = open(os.path.join(REPO, "STATUS.md"), encoding="utf-8").read().splitlines()
    sections = [(i + 1, m.group(1), m.group(2)) for i, l in enumerate(status)
                for m in [re.match(r"^### (\d+\.\d+)\s+(.*)$", l)] if m]
    cited = cited_sections(notes)
    missing = [(ln, num, title) for ln, num, title in sections if num not in cited]
    print(f"1. STATUS sections: {len(sections)}, cited by a note: {len(sections) - len(missing)}, "
          f"not cited: {len(missing)}")
    for ln, num, title in missing:
        print(f"   STATUS {num} (line {ln}): {title[:110]}")
    problems += len(missing)

    method = open(os.path.join(REPO, "METHOD.md"), encoding="utf-8").read()
    i = method.find("### 8.1")
    j = method.find("\n### ", i + 5)
    block = method[i:j if j > 0 else len(method)]
    tokens = sorted(set(re.findall(r"^\| `([A-Z][A-Z0-9_]+)`", block, re.M)))
    dec = [x for r, x in notes.items() if r.startswith("decisions/")]
    def named(tok):
        for x in dec:
            fm = x["fm"]
            hay = " ".join([str(fm.get("title", "")), str(fm.get("description", "")),
                            " ".join(str(a) for a in (fm.get("aliases") or [])), x["body"]])
            if tok in hay:
                return True
        return False
    nodec = [t for t in tokens if not named(t)]
    print(f"2. METHOD 8.1 tokens: {len(tokens)}, named by a decision note: {len(tokens) - len(nodec)}, "
          f"without one: {len(nodec)}")
    for t in nodec:
        print(f"   {t}")
    problems += len(nodec)
    return 1 if (strict and problems) else 0


if __name__ == "__main__":
    sys.exit(main())
