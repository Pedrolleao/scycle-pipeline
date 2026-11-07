#!/usr/bin/env python3
"""
compare_figures.py — are two sets of figures the same picture?

PNGs are compared pixel by pixel (Pillow), SVGs as text after dropping the
<dc:date> line and the per-process clip-path ids matplotlib writes into them. Used as the acceptance test
of changes to the figure scripts that must not change what is drawn.

Usage:  python validation/compare_figures.py <before_dir> <after_dir>
Exit 1 if any figure present in <before_dir> differs or is missing after.
"""
import re, sys
from pathlib import Path
from PIL import Image, ImageChops

before, after = Path(sys.argv[1]), Path(sys.argv[2])
bad = 0
figs = sorted(p for p in before.rglob("*") if p.suffix in (".png", ".svg"))
for p in figs:
    q = after / p.relative_to(before)
    if not q.exists():
        print(f"MISSING {q}"); bad += 1; continue
    if p.suffix == ".png":
        a, b = Image.open(p).convert("RGBA"), Image.open(q).convert("RGBA")
        same = a.size == b.size and ImageChops.difference(a, b).getbbox() is None
    else:
        # matplotlib writes the creation date into the metadata and derives the
        # clip-path ids from object addresses: neither is part of the drawing.
        def strip(t):
            t = re.sub(r"\s*<dc:date>.*?</dc:date>", "", t)
            return re.sub(r"\bp[0-9a-f]{10}\b", "pID", t)
        same = strip(p.read_text()) == strip(q.read_text())
    if not same:
        print(f"DIFF {p.relative_to(before)}"); bad += 1
print(f"{len(figs) - bad} of {len(figs)} figures identical")
sys.exit(1 if bad else 0)
