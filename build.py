"""Build the PDF of the verification note with pandoc/xelatex.

The note's markdown keeps title/author/date/abstract as ordinary markdown
(H1, bold/italic lines, an Abstract section between two ``---`` rules).
pandoc needs them as YAML metadata to typeset a proper title block, so this
script extracts them, writes a temporary _build.md, and runs pandoc.

Usage:  python build.py [input.md]   (default: ExactVerification-SemilocalWeilForm-Note-v2.md)
"""
import re
import subprocess
import sys
from pathlib import Path

src = Path(sys.argv[1] if len(sys.argv) > 1 else "ExactVerification-SemilocalWeilForm-Note-v2.md")
text = src.read_text(encoding="utf-8")
lines = text.splitlines()

title = re.sub(r"^#\s+", "", lines[0]).strip()
author = next(l for l in lines if l.startswith("**Leonardo")).strip("* ")
affil = next(l for l in lines if l.startswith("*Independent") or l.startswith("*Avionyx")).strip("* ")
date = next(l for l in lines if re.fullmatch(r"\*\w+ \d{4}\*", l.strip())).strip("* ")
note_line = next(l for l in lines if l.startswith("**Note, version"))

# abstract: paragraph between '## Abstract' and the following '---'
i_abs = lines.index("## Abstract")
i_end = next(i for i in range(i_abs + 1, len(lines)) if lines[i].strip() == "---")
abstract = "\n".join(l for l in lines[i_abs + 1 : i_end] if l.strip())

body = "\n".join(lines[i_end + 1 :])

meta = f"""---
title: |
  {title}
author: |
  {author}\\
  *{affil}*
date: {date}
geometry: margin=1in
abstract: |
  {abstract}
header-includes: |
  \\providecommand{{\\R}}{{{{\\mathbb{{R}}}}}}
---

{note_line}

"""

build_md = src.with_name("_build.md")
build_md.write_text(meta + body, encoding="utf-8")
out_pdf = src.with_suffix(".pdf")
cmd = ["pandoc", "-s", str(build_md), "-o", str(out_pdf), "--pdf-engine=xelatex"]
print(" ".join(cmd))
rc = subprocess.run(cmd).returncode
build_md.unlink(missing_ok=True)
sys.exit(rc)
