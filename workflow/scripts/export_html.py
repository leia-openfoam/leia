#!/usr/bin/env python3
"""Build the reveal.js deck into a single self-contained (standalone) HTML file.

The editable source is ``velocity-extension.template.html`` (references ``figures/*.png`` and
CDN CSS/JS). This inlines the local figures (as base64 data URIs) and the CDN
assets (reveal.js, theme, notes plugin, MathJax) into ``velocity-extension.html`` -- ONE
shareable file that opens offline, with no figures/ directory or network needed.
Stdlib only (urllib + base64) -- no browser/Playwright dependency.

Best-effort: local images are always inlined; a CDN resource that cannot be
fetched is left as its original link (still a single file that works online).
Returns the output path, or None if the source template is missing.
"""
import base64
import os
import re
import sys
import urllib.request

DEFAULT_SRC = "velocity-extension.template.html"   # editable source (external figure/CDN refs)
DEFAULT_OUT = "velocity-extension.html"            # single self-contained, shareable deck


def _fetch(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def _flatten_sections(html):
    """Linear-flow variant: remove the top-level <section> wrappers (vertical
    stacks), promoting their child slides to top level in document order. The
    wrapper's data-track attribute is stamped onto every promoted child so the
    orientation footer keeps working."""
    tag = re.compile(r"<section\b[^>]*>|</section>")
    track_re = re.compile(r'data-track="([^"]*)"')

    out = []
    pos = 0
    depth = 0
    # (is_wrapper unknown until we see a child) -> collect spans first.
    spans = []   # (open_start, open_end, close_start, close_end, track) of depth-1 wrappers
    stack = []
    for m in tag.finditer(html):
        if m.group(0).startswith("<section"):
            depth += 1
            stack.append((m.start(), m.end(), depth))
        else:
            os_, oe_, d = stack.pop()
            if d == 1:
                inner = html[oe_:m.start()]
                if "<section" in inner:   # wrapper: contains nested slides
                    tm = track_re.search(html[os_:oe_])
                    spans.append((os_, oe_, m.start(), m.end(),
                                  tm.group(1) if tm else None))
            depth -= 1

    for os_, oe_, cs_, ce_, track in spans:
        out.append(html[pos:os_])
        inner = html[oe_:cs_]
        if track:
            # stamp the track on each direct child <section> (all children are
            # depth-1 within `inner`)
            def _stamp(m2, _t=track):
                s = m2.group(0)
                if "data-track" in s:
                    return s
                return s[:-1] + f' data-track="{_t}">'
            # only stamp depth-1 opens within inner
            res, d2, p2 = [], 0, 0
            for m2 in tag.finditer(inner):
                if m2.group(0).startswith("<section"):
                    d2 += 1
                    if d2 == 1:
                        res.append(inner[p2:m2.start()])
                        res.append(_stamp(m2))
                        p2 = m2.end()
                else:
                    d2 -= 1
            res.append(inner[p2:])
            inner = "".join(res)
        out.append(inner)
        pos = ce_
    out.append(html[pos:])
    flat = "".join(out)
    flat = flat.replace("Reveal.initialize({",
                        "window.DECK_LINEAR = true;\n  Reveal.initialize({", 1)
    return flat


def export(slides_dir, src_name=DEFAULT_SRC, out_name=None, linear=False):
    """Build a standalone deck from ``src_name`` (a ``*.template.html``). The output
    name is derived from the template (``foo.template.html`` -> ``foo.html``) unless
    given explicitly, so one function builds any deck (index, sl, ...). If
    ``linear=True``, ALSO write the flat-flow ``*-linear.html`` variant (not tracked
    in git -- generated on demand). Returns the primary path or None."""
    if out_name is None:
        out_name = (src_name.replace(".template.html", ".html")
                    if ".template.html" in src_name else DEFAULT_OUT)
    src = os.path.join(slides_dir, src_name)
    if not os.path.isfile(src):
        print(f"[html] no {src}; skip standalone HTML build")
        return None
    html = open(src, encoding="utf-8").read()

    # 1. Inline LOCAL images (<img src="...">) as base64 data URIs. The path is
    #    resolved relative to the template dir, so this supports both a `figures/`
    #    dir next to the deck AND the thematic single-source layout where the deck
    #    references `../<theme>-article/data/figures/...`. http(s)/data URIs are left.
    def _img(m):
        src = m.group(1)
        p = os.path.normpath(os.path.join(slides_dir, src))
        if os.path.isfile(p):
            data = base64.b64encode(open(p, "rb").read()).decode()
            ext = (os.path.splitext(p)[1].lstrip(".") or "png").lower()
            if ext == "svg":
                ext = "svg+xml"
            return f'src="data:image/{ext};base64,{data}"'
        print(f"[html] WARNING: local image not found, left as-is: {src}")
        return m.group(0)
    html = re.sub(
        r'src="(?!https?://|data:)([^"]+\.(?:png|svg|jpe?g|gif))"', _img, html)

    # 2. Inline CDN stylesheets (<link rel="stylesheet" href="http...">).
    def _css(m):
        url = m.group(1)
        try:
            return "<style>\n" + _fetch(url) + "\n</style>"
        except Exception as exc:  # noqa: BLE001
            print(f"[html] keep CDN css (fetch failed: {exc}): {url}")
            return m.group(0)
    html = re.sub(
        r'<link\b[^>]*rel="stylesheet"[^>]*href="(https?://[^"]+)"[^>]*>', _css, html)

    # 3. Inline CDN scripts (<script ... src="http...">...</script>).
    def _js(m):
        url = m.group(1)
        try:
            js = _fetch(url).replace("</script>", "<\\/script>")  # keep it embeddable
            return "<script>\n" + js + "\n</script>"
        except Exception as exc:  # noqa: BLE001
            print(f"[html] keep CDN js (fetch failed: {exc}): {url}")
            return m.group(0)
    html = re.sub(
        r'<script\b[^>]*\bsrc="(https?://[^"]+)"[^>]*>\s*</script>', _js, html)

    out = os.path.join(slides_dir, out_name)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(html)
    remaining = len(re.findall(r'(?:href|src)="https?://', html))
    print(f"[html] wrote {out} ({len(html.encode())//1024} KiB; "
          f"{remaining} un-inlined CDN refs)")

    # Optional linear-flow variant (same inlined assets, flattened structure).
    # Not tracked in git -- generated on demand (pass linear=True / --linear).
    if linear:
        base, ext = os.path.splitext(out_name)
        lin = os.path.join(slides_dir, f"{base}-linear{ext}")
        with open(lin, "w", encoding="utf-8") as fh:
            fh.write(_flatten_sections(html))
        print(f"[html] wrote {lin} (linear flow)")
    return out


if __name__ == "__main__":
    # Usage: python3 export_html.py [slides_dir] [template_name] [--linear]
    #   e.g. python3 workflow/scripts/export_html.py doc/slides
    #        python3 workflow/scripts/export_html.py doc/slides --linear
    argv = [a for a in sys.argv[1:] if a != "--linear"]
    export(
        argv[0] if len(argv) > 0 else "doc/slides",
        argv[1] if len(argv) > 1 else DEFAULT_SRC,
        linear=("--linear" in sys.argv),
    )
