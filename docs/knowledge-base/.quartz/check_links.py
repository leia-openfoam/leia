#!/usr/bin/env python3
"""Check the pinned GitHub links of the knowledge base against the repository history.

Usage: python3 check_links.py <vault-dir> [--fix] [--infer] [--words] [--list=FILE] [--commits=A,B,...]

Every link of the form https://github.com/leia-openfoam/leia/blob/<commit>/<path>[#L<a>[-L<b>]]
is checked:
 1. the commit exists and holds <path> (else: FIX to the oldest listed commit that holds it);
 2. the anchored lines exist in the file at that commit;
 3. the anchor matches the text around the link. The context is the Markdown line that holds
    the link (a table row or a paragraph line); its tokens are the numbers, the backticked
    words and the long identifiers. The script counts how many of them occur in the anchored
    lines (plus two lines of margin) at the pinned commit and at every other listed commit. If
    another commit matches better by a margin of two and with a score of three or more, the
    writer read that commit's numbering: FIX the pin to it (the line numbers stay). The score
    adds one per context token (the link text and the 250 characters before the link; a table
    row as a whole) found in the anchored lines, five for a \\label named in the context, five
    when a section number of the link text ("STATUS 11.15") holds the anchored line, and three
    when an anchored heading shares two words with the link text. A heading anchor spans its
    section. If no commit matches any token, the link is
    UNVERIFIED (a section link whose text carries no token, or a wrong number); it is listed
    for a reader to check, as CHECK-same when the file is identical at every listed commit (the
    pin cannot be wrong, only the line) and CHECK-diff when it is not (the pin can be wrong).
--infer adds a per-note pass: a writer reads one tree per note, so where two or more of a
note's links into changed files vote decisively for one commit and none for another, every link of
that note into a changed file that has no evidence for its own pin is pinned to it (reported as
INFER; MIXED when the votes split).
--words adds a fallback for an unverified link into a changed file: the distinct content words
(five letters or more, no stop words) of its context are counted in the anchored lines at each
commit; a margin of three and a factor of two move the pin (FIX-words); two or more words and
twice the other commit's count verify it. --list=FILE writes the links that stay ambiguous as a
checklist for a reviewer.
A link whose anchored lines exist in its pinned version only (lines that version added) is never
moved: it points at version-specific content on purpose.
--fix rewrites the pins in place; without it the script only reports. Exit code 1 if any link
is BROKEN (a missing file or line at every listed commit).
Standard library and git only.
"""
import os
import re
import subprocess
import sys
from functools import lru_cache

LINK = re.compile(r"https://github\.com/leia-openfoam/leia/blob/([0-9a-f]{7,40})/([^\s)#\]]+)(?:#L(\d+)(?:-L(\d+))?)?")
TOK = re.compile(r"`([^`]{2,60})`|(?<![\w.])(\d+(?:\.\d+)?(?:e[-+]?\d+)?)(?![\w.])|\b([A-Za-z_][A-Za-z0-9_]*[A-Z_][A-Za-z0-9_]{3,})\b")


@lru_cache(maxsize=None)
def lines_at(commit, path):
    r = subprocess.run(["git", "show", f"{commit}:{path}"], capture_output=True)
    if r.returncode != 0:
        return None
    return r.stdout.decode("utf-8", "replace").splitlines()


HEAD_TEX = re.compile(r"^\s*\\(part|chapter|section|subsection|subsubsection|paragraph)\*?\{")
HEAD_MD = re.compile(r"^(#{1,6})\s+(.*)$")
LEVEL_TEX = {"part": 0, "chapter": 1, "section": 2, "subsection": 3, "subsubsection": 4, "paragraph": 5}
SECNUM = re.compile(r"\b(?:STATUS|METHOD|CLAUDE|section|sec\.|§)\s*(\d+(?:\.\d+){0,3})\b")


AMBIGUOUS = []
NOISE = {"MEASURED", "DERIVED", "HYPOTHESIS", "STATUS", "METHOD", "CLAUDE"}


@lru_cache(maxsize=None)
def own_lines(commit, other, path):
    """Line numbers (1-based) of the file at `commit` that do not exist at `other`: the lines
    that version added. A link into them points at version-specific content and keeps its pin."""
    import difflib
    A, B = lines_at(other, path), lines_at(commit, path)
    if A is None or B is None:
        return frozenset()
    sm = difflib.SequenceMatcher(None, A, B, autojunk=False)
    out = set()
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag in ("insert", "replace"):
            out.update(range(j1 + 1, j2 + 1))
    return frozenset(out)


def version_specific(commit, path, a, b, commits):
    return any(set(range(a, (b or a) + 1)) & own_lines(commit, k, path) for k in commits if k != commit)


def tokens(context):
    ctx = LINK.sub(" ", context)
    ctx = re.sub(r"\]\([^)\s]*\)?", " ", ctx)          # link targets, also a cut-off one
    ctx = re.sub(r"\S*(?:github\.com|/blob/|\.tex|\.md|\.H|\.C|\.yaml|\.csv)\S*", " ", ctx)
    out = set()
    for a, b, c in TOK.findall(ctx):
        t = a or b or c
        if b and "." not in b and "e" not in b and len(b) < 3:   # 1, 2, 10 ...: too common
            continue
        if b and re.fullmatch(r"20\d\d", b):                       # years
            continue
        out.add(t.strip())
    return {t for t in out if t and len(t) >= 2 and t not in NOISE and not re.fullmatch(r"[0-9a-f]{7,40}", t)}


def heading_level(line):
    m = HEAD_TEX.match(line)
    if m:
        return LEVEL_TEX[m.group(1)]
    m = HEAD_MD.match(line)
    if m:
        return len(m.group(1))
    return None


def window(L, a, b):
    """The anchored lines, two lines of margin; a heading anchor spans its section."""
    lo = max(1, a - 2)
    lev = heading_level(L[a - 1]) if b is None or b == a else None
    if lev is not None:
        hi = a
        while hi < len(L) and hi - a < 400:
            nl = heading_level(L[hi])
            if nl is not None and nl <= lev:
                break
            hi += 1
    else:
        hi = min(len(L), (b or a) + 2)
    return lo, hi


def section_of(L, a):
    """The numbered Markdown section ('11.15', '8.1', '4') that holds line a, or None."""
    for i in range(a - 1, -1, -1):
        m = HEAD_MD.match(L[i])
        if m:
            n = re.match(r"(\d+(?:\.\d+){0,3})\b", m.group(2).strip())
            if n:
                return n.group(1)
    return None


def score(commit, path, a, b, toks, link_text):
    L = lines_at(commit, path)
    if L is None or a is None or a > len(L):
        return None
    lo, hi = window(L, a, b)
    text = "\n".join(L[lo - 1:hi])
    s = sum(1 for t in toks if t in text)
    labels = re.findall(r"\b((?:sec|app|tab|fig|eq):[A-Za-z0-9_-]+)", link_text + " " + " ".join(toks))
    s += sum(5 for lab in set(labels) if ("\\label{" + lab + "}") in text)
    if path.endswith(".md"):
        want = SECNUM.findall(link_text)
        have = section_of(L, a)
        if want and have and any(have == w or have.startswith(w + ".") for w in want):
            s += 5
    lev = heading_level(L[a - 1])
    if lev is not None:
        hw = {w.lower() for w in re.findall(r"[A-Za-z][A-Za-z-]{3,}", L[a - 1])}
        cw = {w.lower() for w in re.findall(r"[A-Za-z][A-Za-z-]{3,}", link_text)}
        if len(hw & cw) >= 2:
            s += 3
    return s


STOP = set("""about above after again against along among another around because before being below
between both cannot could does doing during each every first found from further given have having
here into itself later least less made many more most much must never only other over same second
should since some such than that their them then there these they third this those three through
under until upon very were what when where which while with within without would your after
measured derived hypothesis status method claude section article table figure number value line
lines case cases study config file""".split())


def words(text):
    return {w for w in re.findall(r"[a-z][a-z-]{4,}", text.lower()) if w not in STOP}


def word_score(commit, path, a, b, ctx_words):
    L = lines_at(commit, path)
    if L is None or a > len(L):
        return None
    lo, hi = window(L, a, b)
    return len(ctx_words & words("\n".join(L[lo - 1:hi])))


def context_of(line, m):
    """The link text plus the 250 characters before the link (a table row: the whole row)."""
    if line.lstrip().startswith("|"):
        before = line
    else:
        before = line[max(0, m.start() - 250):m.start()]
    j = line.rfind("[", 0, m.start())
    link_text = line[j + 1:m.start()].rstrip("](") if j >= 0 else ""
    return before, link_text


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    fix = "--fix" in sys.argv
    commits = ["8867581", "d1e3414"]
    for x in sys.argv[1:]:
        if x.startswith("--commits="):
            commits = x.split("=", 1)[1].split(",")
    vault = args[0] if args else "docs/knowledge-base"
    stats = {"ok": 0, "fixed": 0, "unverified": 0, "broken": 0}
    report = []
    for root, dirs, files in os.walk(vault):
        dirs[:] = [d for d in dirs if d not in (".obsidian", "templates", "graph3d") and not d.startswith(".")]
        for f in files:
            if not f.endswith(".md"):
                continue
            p = os.path.join(root, f)
            src = open(p, encoding="utf-8").read().splitlines()
            changed = False
            for i, line in enumerate(src):
                new_line = line
                for m in LINK.finditer(line):
                    c, path, a, b = m.group(1), m.group(2), m.group(3), m.group(4)
                    a = int(a) if a else None
                    b = int(b) if b else None
                    rel = os.path.relpath(p, vault)
                    if lines_at(c, path) is None:
                        holder = next((k for k in commits if lines_at(k, path) is not None), None)
                        if holder:
                            stats["fixed"] += 1
                            report.append(f"FIX   {rel}:{i+1}: {path} is not in {c}; pin -> {holder}")
                            new_line = new_line.replace(m.group(0), m.group(0).replace(f"/blob/{c}/", f"/blob/{holder}/"), 1)
                        else:
                            stats["broken"] += 1
                            report.append(f"BROKEN {rel}:{i+1}: {path} exists at none of {commits + [c]}")
                        continue
                    if a is None:
                        stats["ok"] += 1
                        continue
                    before, link_text = context_of(line, m)
                    toks = tokens(before + " " + link_text)
                    here = score(c, path, a, b, toks, link_text)
                    others = {k: score(k, path, a, b, toks, link_text) for k in commits if k != c}
                    if here is None:
                        good = {k: s for k, s in others.items() if s is not None}
                        if good:
                            k = max(good, key=good.get)
                            stats["fixed"] += 1
                            report.append(f"FIX   {rel}:{i+1}: {path}#L{a} is past the end at {c}; pin -> {k}")
                            new_line = new_line.replace(m.group(0), m.group(0).replace(f"/blob/{c}/", f"/blob/{k}/"), 1)
                        else:
                            stats["broken"] += 1
                            report.append(f"BROKEN {rel}:{i+1}: {path}#L{a} exists at no listed commit")
                        continue
                    best = max(others.items(), key=lambda kv: -1 if kv[1] is None else kv[1], default=(None, None))
                    if best[1] is not None and best[1] >= here + 2 and best[1] >= 3 and not version_specific(c, path, a, b, commits):
                        stats["fixed"] += 1
                        report.append(f"FIX   {rel}:{i+1}: {path}#L{a}: {best[1]} tokens at {best[0]} against {here} at {c}; pin -> {best[0]}")
                        new_line = new_line.replace(m.group(0), m.group(0).replace(f"/blob/{c}/", f"/blob/{best[0]}/"), 1)
                    elif here == 0 and "--words" in sys.argv and not all(lines_at(k, path) == lines_at(c, path) for k in commits) \
                            and not version_specific(c, path, a, b, commits):
                        cw = words(LINK.sub(" ", before + " " + link_text))
                        wh = word_score(c, path, a, b, cw) or 0
                        wo = {k: word_score(k, path, a, b, cw) for k in commits if k != c}
                        kbest = max(wo, key=lambda k: -1 if wo[k] is None else wo[k])
                        if wo[kbest] is not None and wo[kbest] >= wh + 3 and wo[kbest] >= 2 * max(1, wh):
                            stats["fixed"] += 1
                            report.append(f"FIX-words {rel}:{i+1}: {path}#L{a}: {wo[kbest]} words at {kbest} against {wh} at {c}; pin -> {kbest}")
                            new_line = new_line.replace(m.group(0), m.group(0).replace(f"/blob/{c}/", f"/blob/{kbest}/"), 1)
                        elif wh >= 2 and wh >= 2 * (wo[kbest] or 0):
                            stats["ok"] += 1                      # verified by its content words
                        else:
                            stats["unverified"] += 1
                            report.append(f"CHECK-diff {rel}:{i+1}: {path}#L{a}{'-L'+str(b) if b else ''} at {c}: words {wh} here, {wo[kbest]} at {kbest}")
                            AMBIGUOUS.append((rel, i + 1, m.group(0), c, kbest))
                    elif here == 0:
                        stats["unverified"] += 1
                        same = all(lines_at(k, path) == lines_at(c, path) for k in commits)
                        tag = "CHECK-same" if same else "CHECK-diff"
                        report.append(f"{tag} {rel}:{i+1}: {path}#L{a}{'-L'+str(b) if b else ''} at {c}: none of {sorted(toks)[:6]} in the anchored lines")
                    else:
                        stats["ok"] += 1
                if new_line != line:
                    src[i] = new_line
                    changed = True
            if changed and fix:
                open(p, "w", encoding="utf-8").write("\n".join(src) + "\n")
    if "--infer" in sys.argv:
        infer_pass(vault, commits, fix, report, stats)
    for x in sys.argv[1:]:
        if x.startswith("--list="):
            with open(x.split("=", 1)[1], "w", encoding="utf-8") as f:
                f.write("# Ambiguous pinned links into files that differ between the commits\n"
                        "# One line per link: status | note:line | link | pinned | other. A reviewer\n"
                        "# sets status to OK, REPIN (to the other commit) or LINE <a>-<b> (a new anchor),\n"
                        "# one line at a time, saving the file after each decision.\n")
                for rel, ln, url, c, k in AMBIGUOUS:
                    f.write(f"TODO | {rel}:{ln} | {url} | {c} | {k}\n")
    for r in sorted(report):
        print(r)
    print(f"check_links: {stats['ok']} ok, {stats['fixed']} {'fixed' if fix else 'to fix'}, {stats['unverified']} unverified, {stats['broken']} broken")
    return 1 if stats["broken"] else 0


def path_of(url):
    return LINK.match(url).group(2)


def a_of(url):
    return int(LINK.match(url).group(3))


def b_of(url):
    g = LINK.match(url).group(4)
    return int(g) if g else None


def infer_pass(vault, commits, fix, report, stats):
    """A writer reads one tree per note. Where a note's links into files that differ between
    the commits give decisive evidence for one numbering (two or more decisive links, none for
    the other), every link of that note into such a file is pinned to that commit."""
    for root, dirs, files in os.walk(vault):
        dirs[:] = [d for d in dirs if d not in (".obsidian", "templates", "graph3d") and not d.startswith(".")]
        for f in files:
            if not f.endswith(".md"):
                continue
            p = os.path.join(root, f)
            rel = os.path.relpath(p, vault)
            src = open(p, encoding="utf-8").read().splitlines()
            votes = {k: 0 for k in commits}
            items = []
            for i, line in enumerate(src):
                for m in LINK.finditer(line):
                    c, path, a, b = m.group(1), m.group(2), m.group(3), m.group(4)
                    if not a or c not in commits:
                        continue
                    if all(lines_at(k, path) == lines_at(commits[0], path) for k in commits):
                        continue                      # identical file: the pin cannot be wrong
                    a = int(a); b = int(b) if b else None
                    before, link_text = context_of(line, m)
                    toks = tokens(before + " " + link_text)
                    sc = {k: score(k, path, a, b, toks, link_text) for k in commits}
                    sc = {k: (-1 if v is None else v) for k, v in sc.items()}
                    best = max(sc, key=sc.get)
                    rest = max(v for k, v in sc.items() if k != best)
                    if sc[best] >= 3 and sc[best] >= rest + 2:
                        votes[best] += 1
                    items.append((i, m.group(0), c, sc))
            winners = [k for k, v in votes.items() if v >= 2 and all(w == 0 for kk, w in votes.items() if kk != k)]
            if not winners:
                if items and any(v for v in votes.values()):
                    report.append(f"MIXED {rel}: votes {votes}; links into changed files left as they are")
                continue
            k = winners[0]
            n = 0
            for i, url, c, sc in items:
                # only a link with no evidence for its own pin moves, and never to a worse match:
                # a link to a line that exists in one version only keeps its pin
                if c != k and sc[c] == 0 and sc[k] >= sc[c] and not version_specific(c, path_of(url), a_of(url), b_of(url), commits):
                    src[i] = src[i].replace(url, url.replace(f"/blob/{c}/", f"/blob/{k}/"), 1)
                    n += 1
            if n:
                stats["fixed"] += n
                report.append(f"INFER {rel}: votes {votes}; {n} link(s) into changed files pinned -> {k}")
                if fix:
                    open(p, "w", encoding="utf-8").write("\n".join(src) + "\n")


if __name__ == "__main__":
    sys.exit(main())
