"""Build the PDF of the verification note with pandoc/xelatex.

The note keeps a human-readable title block (an H1 title, a bold author line, and
two italic lines — affiliation, then version/date) followed by an Abstract section
delimited by two ``---`` rules. pandoc needs these as YAML metadata to typeset a
proper title block, so this script extracts them, writes a temporary _build.md, and
runs pandoc. The abstract is copied verbatim (list, keywords, and all), so page 1 is
a self-contained title page.

Usage:  python build.py [input.md]   (default: ExactVerification-SemilocalWeilForm-Note-v2.md)
"""
import re
import subprocess
import sys
from pathlib import Path

src = Path(sys.argv[1] if len(sys.argv) > 1 else "ExactVerification-SemilocalWeilForm-Note-v2.md")
lines = src.read_text(encoding="utf-8").splitlines()

title = re.sub(r"^#\s+", "", lines[0]).strip()
author = next(l for l in lines if l.startswith("**")).strip("* ")
italics = [l.strip() for l in lines if re.fullmatch(r"\*[^*].*\*", l.strip())]
affil = italics[0].strip("* ")          # first whole-line italic = affiliation
date = italics[1].strip("* ")           # second = version/date

# abstract: verbatim block between '## Abstract' and the following '---'
i_abs = lines.index("## Abstract")
i_end = next(i for i in range(i_abs + 1, len(lines)) if lines[i].strip() == "---")
region = lines[i_abs + 1 : i_end]
while region and not region[0].strip():
    region.pop(0)
while region and not region[-1].strip():
    region.pop()
# indent every line by 2 spaces for the YAML literal block scalar (blank -> empty)
abstract_block = "\n".join(("  " + l) if l.strip() else "" for l in region)

body = "\n".join(lines[i_end + 1 :]).lstrip("\n")

meta = f"""---
title: |
  {title}
author: |
  {author}\\
  *{affil}*
date: {date}
geometry: margin=1in
colorlinks: true
urlcolor: linkblue
linkcolor: linkblue
citecolor: linkblue
header-includes: |
  \\providecommand{{\\R}}{{{{\\mathbb{{R}}}}}}
  \\usepackage{{xcolor}}
  \\definecolor{{linkblue}}{{HTML}}{{1A5FB4}}
abstract: |
{abstract_block}
---

\\newpage

{body}
"""

build_md = src.with_name("_build.md")
build_md.write_text(meta, encoding="utf-8")
out_pdf = Path(sys.argv[2]) if len(sys.argv) > 2 else src.with_suffix(".pdf")
cmd = ["pandoc", "-s", str(build_md), "-o", str(out_pdf), "--pdf-engine=xelatex"]
print(" ".join(cmd))
rc = subprocess.run(cmd).returncode
build_md.unlink(missing_ok=True)
sys.exit(rc)
