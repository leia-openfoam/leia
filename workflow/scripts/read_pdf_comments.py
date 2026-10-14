#!/usr/bin/env python3
"""Read the comments out of an annotated PDF, with the text each one is attached to.

WHY THIS EXISTS. Reviewing a derivation works best when the remark sits on the sentence it
is about. PDF viewers already do that: select text, highlight it, type a note. The note and
the quoted text are stored INSIDE the PDF file, so the whole review comes back in one file
and needs no copying of passages by hand.

WHAT THIS PRINTS. One block per annotation, in reading order, each with
  - the page and the section it falls in (from the PDF's own table of contents, which
    latexmk writes from hyperref, so the section names are the document's own);
  - the text the annotation covers, pulled from under the highlight;
  - the comment typed into it, and who wrote it.

HOW TO PRODUCE THE INPUT. Any viewer that saves annotations into the file works. On this
laptop the PDF lives on the WSL side, so the simplest route is a Windows viewer reading
through \\\\wsl.localhost:

    explorer.exe "$(wslpath -w docs/analysis-amplification-sl/ms.pdf)"

Select text, highlight it, add a note, then Save (not "print to PDF", which discards
annotations). Firefox, Edge, Acrobat Reader, Okular and Xournal++ all do this.

Annotation types that carry a comment are read: Highlight, Underline, StrikeOut, Squiggly,
Text (a sticky note), FreeText, and Ink. A highlight with no typed note is still reported,
because the act of marking a passage is itself feedback.

Usage:
  read_pdf_comments.py <annotated.pdf> [--json] [--quote-chars N]
"""
import argparse
import json
import sys

try:
    import fitz                                  # PyMuPDF
except ImportError:                              # pragma: no cover
    sys.exit("[read_pdf_comments] PyMuPDF is required: pip install pymupdf")

# The annotation kinds that a reviewer can attach a note to. Numbers are PDF annotation
# subtypes as PyMuPDF reports them.
WANTED = {
    fitz.PDF_ANNOT_HIGHLIGHT:  "highlight",
    fitz.PDF_ANNOT_UNDERLINE:  "underline",
    fitz.PDF_ANNOT_STRIKE_OUT: "strikeout",
    fitz.PDF_ANNOT_SQUIGGLY:   "squiggly",
    fitz.PDF_ANNOT_TEXT:       "note",
    fitz.PDF_ANNOT_FREE_TEXT:  "freetext",
    fitz.PDF_ANNOT_INK:        "ink",
}


def section_index(doc):
    """Map each page number to the last section heading that starts at or before it."""
    toc = doc.get_toc() or []
    by_page = {}
    current = None
    for page in range(doc.page_count):
        for level, title, start in toc:
            if start - 1 <= page:
                current = title
        by_page[page] = current
    return by_page


def covered_text(page, annot, limit):
    """The document text under the annotation.

    For a text-markup annotation the vertices give the exact quads the reviewer selected, so
    the quoted text is recovered rather than guessed. For a sticky note, which has no quads,
    fall back to the line its icon sits on.
    """
    quads = None
    try:
        quads = annot.vertices
    except Exception:
        quads = None

    if quads and len(quads) >= 4 and annot.type[0] in (
        fitz.PDF_ANNOT_HIGHLIGHT, fitz.PDF_ANNOT_UNDERLINE,
        fitz.PDF_ANNOT_STRIKE_OUT, fitz.PDF_ANNOT_SQUIGGLY,
    ):
        out = []
        # vertices arrive as a flat list of corner points, four per quad
        for i in range(0, len(quads), 4):
            corners = quads[i:i + 4]
            if len(corners) < 4:
                break
            rect = fitz.Quad(corners).rect
            out.append(page.get_text("text", clip=rect).strip())
        text = " ".join(t for t in out if t)
    else:
        rect = fitz.Rect(annot.rect)
        rect.x0 = 0
        rect.x1 = page.rect.x1
        rect = rect + (0, -4, 0, 4)
        text = page.get_text("text", clip=rect).strip()

    text = " ".join(text.split())
    if limit and len(text) > limit:
        text = text[:limit].rstrip() + " ..."
    return text


def collect(path, limit):
    doc = fitz.open(path)
    sections = section_index(doc)
    rows = []
    for pno in range(doc.page_count):
        page = doc[pno]
        for annot in page.annots() or []:
            kind = WANTED.get(annot.type[0])
            if kind is None:
                continue
            info = annot.info or {}
            rows.append({
                "page": pno + 1,
                "section": sections.get(pno) or "(front matter)",
                "kind": kind,
                "author": (info.get("title") or "").strip(),
                "comment": (info.get("content") or "").strip(),
                "quoted": covered_text(page, annot, limit),
                "y": float(annot.rect.y0),
            })
    doc.close()
    rows.sort(key=lambda r: (r["page"], r["y"]))
    for i, r in enumerate(rows, 1):
        r["n"] = i
        r.pop("y", None)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pdf")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    ap.add_argument("--quote-chars", type=int, default=400,
                    help="truncate the quoted passage (0 = no limit)")
    a = ap.parse_args(argv)

    rows = collect(a.pdf, a.quote_chars)
    if a.json:
        print(json.dumps(rows, indent=2, ensure_ascii=False))
        return 0

    if not rows:
        print(f"No comments found in {a.pdf}.")
        print("If you annotated it, check the viewer SAVED into the same file -- "
              "'print to PDF' discards annotations.")
        return 0

    print(f"{len(rows)} comment(s) in {a.pdf}\n")
    for r in rows:
        head = f"[{r['n']}] p.{r['page']}  {r['kind']}"
        if r["author"]:
            head += f"  ({r['author']})"
        print(head)
        print(f"     section: {r['section']}")
        if r["quoted"]:
            print(f"     on:      “{r['quoted']}”")
        if r["comment"]:
            print(f"     says:    {r['comment']}")
        else:
            print("     says:    (marked, no note typed)")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
