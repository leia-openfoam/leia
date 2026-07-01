#!/usr/bin/env python3
"""Export the reveal.js deck to a single self-contained (standalone) HTML file.

Inlines the local figures (as base64 data URIs) and the CDN CSS/JS (reveal.js,
theme, notes plugin, MathJax) so the result is ONE shareable .html that opens
offline, with no figures/ directory or network needed. Stdlib only (urllib +
base64) -- no browser/Playwright dependency.

Best-effort: local images are always inlined; a CDN resource that cannot be
fetched is left as its original link (still a single file that works online).
Returns the output path, or None if the deck is missing.
"""
import base64
import os
import re
import sys
import urllib.request

DEFAULT_NAME = "leia-velocity-extension.html"


def _fetch(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def export(slides_dir, html_name=DEFAULT_NAME):
    index = os.path.join(slides_dir, "index.html")
    if not os.path.isfile(index):
        print(f"[html] no {index}; skip standalone HTML export")
        return None
    html = open(index, encoding="utf-8").read()

    # 1. Inline local figures (<img src="figures/...">) as base64 data URIs.
    def _img(m):
        src = m.group(1)
        p = os.path.join(slides_dir, src)
        if os.path.isfile(p):
            data = base64.b64encode(open(p, "rb").read()).decode()
            ext = (os.path.splitext(p)[1].lstrip(".") or "png").lower()
            return f'src="data:image/{ext};base64,{data}"'
        return m.group(0)
    html = re.sub(r'src="((?:\./)?figures/[^"]+)"', _img, html)

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

    out = os.path.join(slides_dir, html_name)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(html)
    remaining = len(re.findall(r'(?:href|src)="https?://', html))
    print(f"[html] wrote {out} ({len(html.encode())//1024} KiB; "
          f"{remaining} un-inlined CDN refs)")
    return out


if __name__ == "__main__":
    export(
        sys.argv[1] if len(sys.argv) > 1 else "doc/slides",
        sys.argv[2] if len(sys.argv) > 2 else DEFAULT_NAME,
    )
