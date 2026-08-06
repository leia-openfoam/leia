#!/usr/bin/env python3
"""Export the reveal.js deck (doc/slides/velocity-extension.html) to a shareable PDF.

Uses Playwright's headless Chromium. Reveal's ``?print-pdf`` mode lays every
slide out as its own page; we wait for network idle + MathJax typesetting before
printing so equations and images are captured. Best-effort: if Playwright or its
browser is missing, the deck HTML/figures are still produced -- only the PDF is
skipped (with a warning), so the workflow never fails on account of the PDF.

The output is named after the deck (``leia-code-design.pdf`` by default), not
``index.pdf``, so it is meaningful when shared.
"""
import os
import subprocess
import sys

DEFAULT_NAME = "leia-code-design.pdf"


def _render(slides_dir, pdf_name):
    """Do the actual render with Playwright's sync API. Must run in a process
    with NO running asyncio loop (see export())."""
    from playwright.sync_api import sync_playwright

    index = os.path.join(slides_dir, "velocity-extension.html")
    out = os.path.join(slides_dir, pdf_name)
    url = "file://" + os.path.abspath(index) + "?print-pdf"
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.goto(url, wait_until="networkidle", timeout=90000)
        # MathJax v3: startup document reaches state 10 when typesetting is done.
        try:
            page.wait_for_function(
                "() => !window.MathJax || (window.MathJax.startup"
                " && MathJax.startup.document"
                " && MathJax.startup.document.state() >= 10)",
                timeout=45000,
            )
        except Exception:  # noqa: BLE001
            pass
        page.wait_for_timeout(2500)  # let reveal finish the print layout
        page.pdf(path=out, prefer_css_page_size=True, print_background=True)
        browser.close()
    print(f"[pdf] wrote {out}")


def export(slides_dir, pdf_name=DEFAULT_NAME):
    """Render the deck to a shareable PDF. Runs the Playwright work in a fresh
    subprocess: Snakemake run: blocks execute inside an asyncio loop, where the
    Playwright *sync* API refuses to run. Best-effort -- a failure only skips the
    PDF (with a warning); the figures and deck HTML are already produced."""
    index = os.path.join(slides_dir, "velocity-extension.html")
    if not os.path.isfile(index):
        print(f"[pdf] no {index}; skip PDF export")
        return None
    out = os.path.join(slides_dir, pdf_name)
    try:
        subprocess.run(
            [sys.executable, os.path.abspath(__file__), slides_dir, pdf_name],
            check=True, timeout=180,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"[pdf] export failed ({exc}); skip PDF export "
              f"(need: pip install --user playwright && python -m playwright install chromium)")
        return None
    return out


if __name__ == "__main__":
    _render(
        sys.argv[1] if len(sys.argv) > 1 else "doc/slides",
        sys.argv[2] if len(sys.argv) > 2 else DEFAULT_NAME,
    )
