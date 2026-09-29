#!/usr/bin/env python3
"""Check the leia knowledge base (docs/knowledge-base) before a commit and before a site build.

Usage: python3 check_kb.py <vault-dir>

Checks, each written from a way the vault can silently break:
 1. every *.md outside templates/ has frontmatter with title, description, kind, status, part,
    tags (the kind is one of the tags), date; the vocabularies are closed;
 2. every [[target]], [[target#heading]], [[target|alias]] and ![[target]] resolves to a note
    (by its vault path, or by an alias), or to a file under assets/; a note in a folder must be
    linked with its folder, because Quartz resolves links from the vault root;
 3. every note in decisions/ has a line in decision-log.md, every note in retractions/ a line in
    retraction-log.md;
 4. every settled decision has date_settled and a non-empty decided_by;
 5. no file named log.* and no folder or file name that the repository .gitignore traps
    (a single digit, "<digit>.<digit>"); a non-Markdown file outside assets/, graph3d/,
    .obsidian/, .quartz/ and templates/ is reported (the site would publish it);
 6. "## Log" is the last "##" section of every note that is not an index.
Exit code 1 on any failure; every message is "<file>: <message>".
Standard library only (the same parser feeds build_graph.py).
"""
import os
import re
import sys

KINDS = {"hub", "concept", "model", "decision", "retraction", "case", "study", "session", "index"}
STATUSES = {"settled", "open", "retracted", "voided", "candidate"}
PARTS = {"advection", "viscosity", "surface-tension", "mass-flux", "gradient-control",
         "verification", "all"}
REQUIRED = ("title", "description", "kind", "status", "part", "tags", "date")
SKIP_DIRS = {".obsidian", ".quartz", "templates", "graph3d"}
LINK_RE = re.compile(r"(!?)\[\[([^\]\|#]+)(#[^\]\|]*)?(\|[^\]]*)?\]\]")
FENCE_RE = re.compile(r"```.*?```", re.S)
CODE_RE = re.compile(r"`[^`\n]*`")


def parse_scalar(v):
    v = v.strip()
    if len(v) >= 2 and v[0] == v[-1] and v[0] in "\"'":
        return v[1:-1]
    return v


def parse_list(v):
    v = v.strip()
    if v.startswith("[") and v.endswith("]"):
        inner = v[1:-1].strip()
        return [parse_scalar(x) for x in inner.split(",") if x.strip()] if inner else []
    return [parse_scalar(v)] if v else []


def yaml_hazards(block):
    """Unquoted scalars that a strict YAML parser rejects or misreads: ': ' or ' #' inside, or a
    leading special character. check_kb's own parser is lenient; the site's is not."""
    out = []
    for line in block.splitlines():
        m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*):\s*(.*)$", line)
        if not m:
            continue
        key, val = m.group(1), m.group(2).strip()
        if not val or (val[0] == val[-1] and val[0] in "\"'" and len(val) > 1):
            continue
        if key in ("tags", "aliases", "code", "sources", "decided_by") and val.startswith("["):
            inner = val[1:-1] if val.endswith("]") else val[1:]
            if ": " in inner or " #" in inner:
                out.append(f"frontmatter '{key}': quote the list items that contain ': ' or ' #'")
            continue
        if ": " in val or " #" in val or val[0] in "[{&*!|>%@`":
            out.append(f"frontmatter '{key}': quote this value (it contains ': ', ' #' or starts with '{val[0]}')")
    return out


def parse_frontmatter(text):
    """Minimal YAML: key: scalar | key: [a, b] | key:\n  - a. Returns (dict, body) or (None, text)."""
    if not text.startswith("---\n"):
        return None, text
    end = text.find("\n---\n", 4)
    if end < 0:
        return None, text
    block, body = text[4:end], text[end + 5:]
    fm, key = {}, None
    for line in block.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        m = re.match(r"^\s*-\s*(.*)$", line)
        if m and key is not None and isinstance(fm.get(key), list):
            fm[key].append(parse_scalar(m.group(1)))
            continue
        m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*):\s*(.*)$", line)
        if not m:
            continue
        key, val = m.group(1), m.group(2)
        if key in ("tags", "aliases", "code", "sources", "decided_by"):
            fm[key] = parse_list(val) if val.strip() else []
        else:
            fm[key] = parse_scalar(val)
    return fm, body


def note_files(vault):
    for root, dirs, files in os.walk(vault):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS and not d.startswith(".")]
        for f in files:
            yield os.path.join(root, f)


def load_notes(vault):
    """Return {vault-relative path without .md: {"fm":..., "body":..., "path":...}}."""
    notes = {}
    for p in note_files(vault):
        if not p.endswith(".md"):
            continue
        rel = os.path.relpath(p, vault)[:-3].replace(os.sep, "/")
        text = open(p, encoding="utf-8").read()
        fm, body = parse_frontmatter(text)
        notes[rel] = {"fm": fm, "body": body, "path": p, "text": text}
    return notes


def strip_code(body):
    return CODE_RE.sub("", FENCE_RE.sub("", body))


def links_of(body):
    for m in LINK_RE.finditer(strip_code(body)):
        yield m.group(2).strip(), (m.group(3) or "")[1:].strip(), bool(m.group(1))


def headings_of(body):
    return [re.sub(r"\s+", " ", h.strip()) for h in re.findall(r"^#{1,6}\s+(.+?)\s*$", body, re.M)]


def check(vault):
    errors = []
    notes = load_notes(vault)
    aliases = {}
    for rel, n in notes.items():
        fm = n["fm"]
        if fm:
            for a in fm.get("aliases", []) or []:
                aliases.setdefault(a, rel)
    rel_of_file = lambda p: os.path.relpath(p, vault).replace(os.sep, "/")

    # 5. names and stray files
    for p in note_files(vault):
        rel = rel_of_file(p)
        parts = rel.split("/")
        base = parts[-1]
        if base.startswith("log.") or base == "log":
            errors.append(f"{rel}: the name matches the repository .gitignore pattern log.*")
        for comp in parts:
            if re.fullmatch(r"[0-9]", comp) or re.match(r"^[0-9]\.[0-9]", comp):
                errors.append(f"{rel}: the path component '{comp}' is ignored by the repository .gitignore")
        if not p.endswith(".md") and parts[0] != "assets":
            errors.append(f"{rel}: a non-Markdown file outside assets/ (the site would publish it)")

    for rel, n in sorted(notes.items()):
        fm, body = n["fm"], n["body"]
        # 1. frontmatter
        if fm is None:
            errors.append(f"{rel}.md: no frontmatter")
            continue
        head = n["text"][4:n["text"].find("\n---\n", 4)]
        for h in yaml_hazards(head):
            errors.append(f"{rel}.md: {h}")
        for k in REQUIRED:
            if k not in fm or fm[k] in ("", [], None):
                errors.append(f"{rel}.md: frontmatter field '{k}' is missing or empty")
        kind, status, part = fm.get("kind"), fm.get("status"), fm.get("part")
        if kind not in KINDS:
            errors.append(f"{rel}.md: kind '{kind}' is not one of {sorted(KINDS)}")
        if status not in STATUSES:
            errors.append(f"{rel}.md: status '{status}' is not one of {sorted(STATUSES)}")
        if part not in PARTS:
            errors.append(f"{rel}.md: part '{part}' is not one of {sorted(PARTS)}")
        tags = fm.get("tags") or []
        if kind and kind not in tags:
            errors.append(f"{rel}.md: tags {tags} do not contain the kind '{kind}'")
        folder = rel.split("/")[0] if "/" in rel else ""
        expected = {"hubs": "hub", "concepts": "concept", "models": "model", "decisions": "decision",
                    "retractions": "retraction", "cases": "case", "studies": "study",
                    "sessions": "session"}.get(folder)
        if expected and kind != expected and not rel.endswith("/index"):
            errors.append(f"{rel}.md: kind '{kind}' does not match its folder '{folder}/' ({expected})")
        # 4. settled decisions
        if kind == "decision" and status == "settled":
            if not fm.get("date_settled"):
                errors.append(f"{rel}.md: a settled decision needs date_settled")
            if not fm.get("decided_by"):
                errors.append(f"{rel}.md: a settled decision needs a non-empty decided_by")
        # 6. Log last
        if kind != "index":
            h2 = re.findall(r"^##\s+(.+?)\s*$", body, re.M)
            if not h2 or h2[-1].strip() != "Log":
                errors.append(f"{rel}.md: '## Log' must be the last '##' section (last is {h2[-1] if h2 else 'none'})")
        # 2. links
        for target, heading, is_embed in links_of(body):
            t = target[:-3] if target.endswith(".md") else target
            if is_embed or "/" in t and t.split("/")[0] == "assets":
                if not os.path.isfile(os.path.join(vault, target)):
                    errors.append(f"{rel}.md: embedded file [[{target}]] does not exist")
                continue
            resolved = None
            if t in notes:
                resolved = t
            elif t in aliases:
                resolved = aliases[t]
            else:
                # a bare slug that exists in a folder: Quartz cannot resolve it from the root
                hits = [r for r in notes if r.split("/")[-1] == t]
                if hits:
                    errors.append(f"{rel}.md: link [[{target}]] must carry its folder: [[{hits[0]}]]")
                else:
                    errors.append(f"{rel}.md: link [[{target}]] does not resolve")
                continue
            if heading:
                hs = [h.lower() for h in headings_of(notes[resolved]["body"])]
                if heading.lower() not in hs:
                    errors.append(f"{rel}.md: link [[{target}#{heading}]]: heading not found in {resolved}.md")
    # 3. log coverage
    for folder, log in (("decisions", "decision-log"), ("retractions", "retraction-log")):
        logtext = notes.get(log, {}).get("body", "") if log in notes else ""
        if log not in notes:
            errors.append(f"{log}.md: missing")
            continue
        for rel in notes:
            if rel.startswith(folder + "/") and not rel.endswith("/index"):
                if f"[[{rel}" not in logtext:
                    errors.append(f"{rel}.md: has no line in {log}.md")
    return errors, notes, aliases


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    vault = sys.argv[1]
    errors, notes, _ = check(vault)
    for e in errors:
        print(e)
    n = sum(1 for r, x in notes.items() if x["fm"] and x["fm"].get("kind") != "index")
    print(f"check_kb: {len(notes)} files, {n} notes, {len(errors)} problem(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
