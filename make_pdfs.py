#!/usr/bin/env python3
"""Render each page of the docs/ website to a PDF in docs/pdf/.

The pages are plain markdown served by Jekyll with MathJax, so the maths is
written with kramdown's ``$$...$$`` for both inline and display formulae.  This
script turns each page into a standalone print-styled HTML file (maths as
MathML, which Chrome renders natively) and prints it with headless Chrome, then
stamps page numbers on with a LaTeX-generated overlay.

Requires pandoc, google-chrome, pdflatex and pypdf.  Run from anywhere:

    python3 make_pdfs.py
"""

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DOCS = ROOT / "docs"
OUT = DOCS / "pdf"

SITE_TITLE = "Fishing for biodiversity by balanced harvesting"
SITE_DESCRIPTION = "A mizer reimplementation of Law & Plank (2023), Fish and Fisheries"
SITE_URL = "https://gustavdelius.github.io/balanced-harvest/"
PDF_URL = SITE_URL + "pdf/"

SITE_BRAND = "Balanced harvesting"


def site_pages():
    """The pages listed in docs/_data/nav.yml, in the site's own order.

    Read rather than hard-coded so that adding a page to the site's navigation
    is enough to get a PDF of it.  The file is simple enough not to need a YAML
    parser: every page is a `url:` line, "/" for index and "/name.html"
    otherwise, matching Jekyll's page.url.
    """
    nav = (DOCS / "_data" / "nav.yml").read_text(encoding="utf-8")
    names = []
    for url in re.findall(r"^\s*-?\s*url:\s*(\S+)\s*$", nav, re.M):
        stem = "index" if url == "/" else url.lstrip("/")[:-len(".html")]
        md = DOCS / (stem + ".md")
        if not md.exists():
            sys.exit(f"nav.yml lists {url}, but {md} does not exist")
        names.append(md.name)
    if not names:
        sys.exit("no pages found in docs/_data/nav.yml")
    return names


# --------------------------------------------------------------------------
# Markdown preprocessing
# --------------------------------------------------------------------------

def fix_math(text):
    """Turn kramdown's inline ``$$...$$`` into pandoc's inline ``$...$``.

    Display maths on these pages is always a ``$$`` alone on its own line, so
    those lines delimit blocks that are passed through untouched.  Everything
    else, including fenced code, is left alone apart from the delimiters.
    """
    out, in_display, in_code = [], False, False
    for line in text.split("\n"):
        stripped = line.strip()
        if stripped.startswith("```"):
            in_code = not in_code
            out.append(line)
            continue
        if not in_code and stripped == "$$":
            in_display = not in_display
            out.append(line)
            continue
        if in_code or in_display:
            out.append(line)
        else:
            out.append(line.replace("$$", "$"))
    return "\n".join(out)


def preprocess(path, pages):
    text = path.read_text(encoding="utf-8")
    text = fix_math(text)
    # Cross-page links point at the published sibling PDF.  An absolute URL,
    # because headless Chrome resolves a relative one against the local file
    # path it printed from, which would bake this machine's directories in.
    text = re.sub(r"\]\((" + "|".join(p[:-3] for p in pages) + r")\.md\)",
                  r"](" + PDF_URL + r"\1.pdf)", text)
    # A leading level-1 heading becomes the title block instead of body text.
    title = None
    lines = text.split("\n")
    if lines and lines[0].startswith("# "):
        title = lines[0][2:].strip()
        text = "\n".join(lines[1:]).lstrip("\n")
    return text, title


# --------------------------------------------------------------------------
# HTML template
# --------------------------------------------------------------------------

TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>$pagetitle$</title>
<style>
@page { size: A4; margin: 20mm 18mm 22mm 18mm; }

html { font-size: 10.5pt; }
body {
  font-family: "DejaVu Serif", Georgia, "Times New Roman", serif;
  line-height: 1.5;
  color: #1b1f23;
  margin: 0;
  -webkit-print-color-adjust: exact;
  print-color-adjust: exact;
}

/* ---- title block ---- */
header.page-head { margin: 0 0 1.6em 0; padding-bottom: 0.8em;
  border-bottom: 2px solid #159957; }
header.page-head .kicker {
  font-family: "DejaVu Sans", "Helvetica Neue", Arial, sans-serif;
  font-size: 0.72rem; letter-spacing: 0.09em; text-transform: uppercase;
  color: #159957; margin: 0 0 0.5em 0; }
header.page-head h1 {
  font-family: "DejaVu Sans", "Helvetica Neue", Arial, sans-serif;
  font-size: 1.85rem; line-height: 1.2; font-weight: 600;
  color: #155799; margin: 0; }
header.page-head .subtitle {
  font-size: 0.85rem; color: #586069; margin: 0.6em 0 0 0; font-style: italic; }

/* ---- headings ---- */
h1, h2, h3, h4, h5, h6 {
  font-family: "DejaVu Sans", "Helvetica Neue", Arial, sans-serif;
  font-weight: 600; color: #155799; line-height: 1.25;
  margin: 1.6em 0 0.5em 0; break-after: avoid; page-break-after: avoid; }
h2 { font-size: 1.32rem; padding-bottom: 0.2em; border-bottom: 1px solid #e1e4e8; }
h3 { font-size: 1.1rem; color: #1b3f66; }
h4 { font-size: 0.98rem; color: #24292e; }
h2 + p, h3 + p, h4 + p { margin-top: 0.35em; }

p { margin: 0 0 0.8em 0; orphans: 2; widows: 2; }

a { color: #155799; text-decoration: none; border-bottom: 1px solid #c8d4e0; }

/* ---- lists ---- */
ul, ol { margin: 0 0 0.9em 0; padding-left: 1.5em; }
li { margin-bottom: 0.25em; }
li > p { margin-bottom: 0.4em; }

/* ---- quotes ---- */
blockquote { margin: 1em 0; padding: 0.1em 0 0.1em 1em;
  border-left: 3px solid #159957; color: #444d56; }
blockquote p { margin-bottom: 0.4em; }

/* ---- code ---- */
code, kbd, samp {
  font-family: "DejaVu Sans Mono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.86em; background: #f3f5f7; padding: 0.12em 0.3em;
  border-radius: 3px; }
pre { background: #f6f8fa; border: 1px solid #e1e4e8; border-radius: 4px;
  padding: 0.7em 0.9em; overflow-wrap: break-word; white-space: pre-wrap;
  break-inside: avoid; page-break-inside: avoid; margin: 0 0 1em 0; }
pre code { background: none; padding: 0; font-size: 0.82em; }

/* ---- tables ---- */
table { border-collapse: collapse; width: 100%; margin: 1em 0 1.2em 0;
  font-size: 0.86em; break-inside: auto; }
th, td { border: 1px solid #dfe2e5; padding: 0.35em 0.55em;
  vertical-align: top; text-align: left; }
th { background: #f1f5f9; font-weight: 600;
  font-family: "DejaVu Sans", "Helvetica Neue", Arial, sans-serif; }
tr { break-inside: avoid; page-break-inside: avoid; }
thead { display: table-header-group; }
caption { caption-side: bottom; font-size: 0.82em; color: #586069;
  padding-top: 0.4em; text-align: left; }

/* ---- figures ---- */
img { max-width: 100%; max-height: 21cm; object-fit: contain;
  display: block; margin: 0.6em auto; }
figure { margin: 1.2em 0; break-inside: avoid; page-break-inside: avoid;
  text-align: center; }
figcaption { font-size: 0.84em; color: #586069; margin-top: 0.4em;
  font-style: italic; }

/* ---- maths ---- */
math { font-size: 1.02em; }
math[display="block"] { margin: 0.9em 0; }

hr { border: none; border-top: 1px solid #e1e4e8; margin: 1.8em 0; }

/* ---- colophon ---- */
footer.colophon { margin-top: 2.2em; padding-top: 0.7em;
  border-top: 1px solid #e1e4e8; font-size: 0.75rem; color: #6a737d;
  font-family: "DejaVu Sans", "Helvetica Neue", Arial, sans-serif; }
footer.colophon a { border: none; color: #6a737d; }
</style>
</head>
<body>
<header class="page-head">
  <p class="kicker">$kicker$</p>
  <h1>$maintitle$</h1>
  <p class="subtitle">$subtitle$</p>
</header>
$body$
<footer class="colophon">
  $sourceline$
</footer>
</body>
</html>
"""


# --------------------------------------------------------------------------
# Page-number overlay
# --------------------------------------------------------------------------

def latex_escape(s):
    for a, b in [("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"),
                 ("$", r"\$"), ("#", r"\#"), ("_", r"\_"), ("{", r"\{"),
                 ("}", r"\}"), ("~", r"\textasciitilde{}"),
                 ("^", r"\textasciicircum{}")]:
        s = s.replace(a, b)
    return s


def make_overlay(n_pages, label, workdir):
    """A LaTeX-built PDF of n_pages blank pages carrying just the footer."""
    tex = workdir / "numbers.tex"
    tex.write_text(r"""
\documentclass[a4paper,10pt]{article}
\usepackage[a4paper,margin=18mm,footskip=12mm]{geometry}
\usepackage{fancyhdr}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{xcolor}
\definecolor{rulegrey}{HTML}{6A737D}
\pagestyle{fancy}
\fancyhf{}
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}
\fancyfoot[L]{\footnotesize\sffamily\color{rulegrey}""" + label + r"""}
\fancyfoot[R]{\footnotesize\sffamily\color{rulegrey}\thepage\ of """
    + str(n_pages) + r"""}
\begin{document}
""" + "\n".join([r"\null\newpage"] * (n_pages - 1)) + r"""
\null
\end{document}
""", encoding="utf-8")
    subprocess.run(["pdflatex", "-interaction=batchmode", "-halt-on-error",
                    "numbers.tex"], cwd=workdir, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return workdir / "numbers.pdf"


def stamp_numbers(pdf_path, label, workdir):
    from pypdf import PdfReader, PdfWriter
    reader = PdfReader(str(pdf_path))
    n = len(reader.pages)
    overlay = PdfReader(str(make_overlay(n, latex_escape(label), workdir)))
    if len(overlay.pages) != n:
        print(f"    ! overlay has {len(overlay.pages)} pages, content has {n};"
              " skipping page numbers")
        return
    writer = PdfWriter()
    for page, num in zip(reader.pages, overlay.pages):
        page.merge_page(num)
        writer.add_page(page)
    writer.add_metadata({"/Title": label, "/Creator": SITE_TITLE})
    with open(pdf_path, "wb") as fh:
        writer.write(fh)


# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------

def build(md_name, pages, template_path, workdir):
    src = DOCS / md_name
    stem = src.stem
    body, heading = preprocess(src, pages)
    main_title = heading if heading is not None else SITE_TITLE

    page_url = SITE_URL + ("" if stem == "index" else stem + ".html")
    source_line = (f'Page <a href="{page_url}">{page_url}</a> of the site '
                   f'accompanying the mizer reimplementation of Law &amp; Plank '
                   f'(2023). Generated by <code>make_pdfs.py</code>.')

    md_tmp = workdir / (stem + ".md")
    md_tmp.write_text(body, encoding="utf-8")

    # The HTML is written into docs/ so that relative figure paths resolve.
    html_tmp = DOCS / ("." + stem + ".pdfsrc.html")
    subprocess.run([
        "pandoc", str(md_tmp),
        "--from", "markdown+tex_math_dollars+raw_html+pipe_tables",
        "--to", "html5", "--mathml", "--standalone",
        "--template", str(template_path),
        "--metadata", f"pagetitle={main_title}",
        "--metadata", f"maintitle={main_title}",
        "--metadata", f"kicker={SITE_BRAND}",
        "--metadata", f"subtitle={SITE_DESCRIPTION}",
        "--metadata", f"sourceline={source_line}",
        "-o", str(html_tmp),
    ], check=True)

    OUT.mkdir(parents=True, exist_ok=True)
    pdf = OUT / (stem + ".pdf")
    try:
        subprocess.run([
            "google-chrome", "--headless", "--disable-gpu", "--no-sandbox",
            "--no-pdf-header-footer", "--run-all-compositor-stages-before-draw",
            "--virtual-time-budget=20000",
            f"--print-to-pdf={pdf}", html_tmp.as_uri(),
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
           timeout=300)
    finally:
        html_tmp.unlink(missing_ok=True)

    stamp_numbers(pdf, main_title, workdir)
    from pypdf import PdfReader
    return pdf, len(PdfReader(str(pdf)).pages)


def main():
    for tool in ("pandoc", "google-chrome", "pdflatex"):
        if shutil.which(tool) is None:
            sys.exit(f"required tool not found: {tool}")

    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        template_path = workdir / "page.html"
        template_path.write_text(TEMPLATE, encoding="utf-8")
        pages = site_pages()
        for md_name in pages:
            pdf, n = build(md_name, pages, template_path, workdir)
            print(f"  {pdf.relative_to(ROOT)}  ({n} pages,"
                  f" {pdf.stat().st_size // 1024} kB)")


if __name__ == "__main__":
    main()
