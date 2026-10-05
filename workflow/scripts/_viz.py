"""Shared visual vocabulary for the pipeline figures and the HTML report.

One place for the status glyphs, the completeness glyphs and the row ordering,
so the overview grid, the per-pathway grids, the complex / process grids, the
cycle maps, the locus maps and report.html all read the same way.

Everything specific to this cycle (pathway order, palette, labels, file names)
lives in _domain.py and is re-exported here; this file is identical in the
nitrogen and sulfur sister pipelines.
"""

from __future__ import annotations

import csv
from pathlib import Path

from _domain import (CALLS_TSV, CAT_ORDER, CYCLE_LETTER, CYCLE_NAME,  # noqa: F401
                     LOCI_TSV, PATHWAY_COLOR, PATHWAY_COLOR_DARK,
                     PATHWAY_LABEL, PATHWAY_SHORT, PRETTY_REPLACE,
                     RULED_OUT_LABEL)

# Ink / chrome.
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
FAINT = "#d9d8d2"          # absent-cell ring
UNANNOTATED = "#e4e3dd"    # locus-map genes outside the target set
ACCENT = "#2a78d6"         # single hue for the completeness grids

STATUS_CODE = {
    "confirmed":     2,
    "domain-only":   1,
    "narrow-no-IPR": 1,
    "disqualified": -1,
    "absent":        0,
}
PRESENT_CODES = (1, 2)
STATUS_LABEL = {
    2:  "confirmed",
    1:  "domain-only (HMM signature, no BLAST support)",
    -1: "disqualified (failed the homology-trap gate)",
    0:  "absent",
}
STATUS_SHORT = {2: "confirmed", 1: "domain-only", -1: "disqualified", 0: "absent"}


def load_calls(path: Path) -> dict[str, dict]:
    """{target_id: row} from a per-sample calls table."""
    out: dict[str, dict] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["target_id"]] = row
    return out


def call_codes(calls: dict[str, dict]) -> dict[str, int]:
    return {tid: STATUS_CODE.get(r["status"], 0) for tid, r in calls.items()}


# ───────────────────────────── status glyphs ─────────────────────────────────

def draw_status_glyph(ax, x, y, code, hue, radius=0.36, zorder=3):
    """One presence/absence glyph. Status is carried by SHAPE, not by tint:

       2 confirmed    → solid disc in the pathway hue
       1 domain-only  → ring in the hue, left half filled
      -1 disqualified → grey ring with a cross
       0 absent       → faint empty ring
    """
    import matplotlib.patches as mp
    if code == 2:
        ax.add_patch(mp.Circle((x, y), radius, facecolor=hue, edgecolor=hue,
                               linewidth=0.6, zorder=zorder, clip_on=False))
    elif code == 1:
        ax.add_patch(mp.Circle((x, y), radius, facecolor="white", edgecolor=hue,
                               linewidth=1.1, zorder=zorder, clip_on=False))
        ax.add_patch(mp.Wedge((x, y), radius, 90, 270, facecolor=hue,
                              edgecolor="none", zorder=zorder + 0.1,
                              clip_on=False))
    elif code == -1:
        ax.add_patch(mp.Circle((x, y), radius, facecolor="white",
                               edgecolor=MUTED, linewidth=0.9, zorder=zorder,
                               clip_on=False))
        d = radius * 0.42
        for sx in (1, -1):
            ax.plot([x - d, x + d], [y - sx * d, y + sx * d], color=MUTED,
                    lw=0.9, solid_capstyle="round", zorder=zorder + 0.1,
                    clip_on=False)
    else:
        ax.add_patch(mp.Circle((x, y), radius, facecolor="none",
                               edgecolor=FAINT, linewidth=0.6,
                               zorder=zorder - 1, clip_on=False))


def draw_status_legend(ax, x, y, hue, step=1.0, radius=0.36, fontsize=8,
                       horizontal=False, char_w=0.0, labels=None):
    """Draw the four status glyphs with labels, in data coordinates, starting
    at (x, y). Vertical by default (one item per `step`); horizontal layout
    needs `char_w` (approx. data units per character) to space the items."""
    labels = labels or STATUS_LABEL
    cx, cy = x, y
    for code in (2, 1, -1, 0):
        draw_status_glyph(ax, cx, cy, code, hue, radius=radius)
        ax.text(cx + radius * 1.9, cy, labels[code], ha="left", va="center",
                fontsize=fontsize, color=INK2, clip_on=False)
        if horizontal:
            cx += radius * 1.9 + len(labels[code]) * char_w + 1.6
        else:
            cy -= step


# ─────────────────────── completeness glyphs (complex / process) ─────────────

def draw_completeness_glyph(ax, x, y, state, text="", hue=ACCENT, radius=0.36,
                            fontsize=6.5, zorder=3):
    """Complex / process cell.

      complete  → solid disc
      partial   → ring with the "n/N" fraction inside
      ruled_out → faint ring with a cross (an excluded gene is present)
      absent    → faint empty ring
    """
    import matplotlib.patches as mp
    if state == "complete":
        ax.add_patch(mp.Circle((x, y), radius, facecolor=hue, edgecolor=hue,
                               linewidth=0.6, zorder=zorder, clip_on=False))
    elif state == "partial":
        ax.add_patch(mp.Circle((x, y), radius, facecolor="white", edgecolor=hue,
                               linewidth=1.2, zorder=zorder, clip_on=False))
        if text:
            ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
                    color=INK2, zorder=zorder + 0.2, clip_on=False)
    elif state == "ruled_out":
        ax.add_patch(mp.Circle((x, y), radius, facecolor="none",
                               edgecolor=FAINT, linewidth=0.6,
                               zorder=zorder - 1, clip_on=False))
        d = radius * 0.38
        for sx in (1, -1):
            ax.plot([x - d, x + d], [y - sx * d, y + sx * d], color=MUTED,
                    lw=0.9, solid_capstyle="round", zorder=zorder,
                    clip_on=False)
    else:
        ax.add_patch(mp.Circle((x, y), radius, facecolor="none",
                               edgecolor=FAINT, linewidth=0.6,
                               zorder=zorder - 1, clip_on=False))


def draw_completeness_legend(ax, x, y, items, step=1.0, radius=0.36,
                             fontsize=8, hue=ACCENT):
    """`items` = [(state, example_text, label), ...], stacked vertically."""
    for k, (state, text, label) in enumerate(items):
        draw_completeness_glyph(ax, x, y - k * step, state, text, hue=hue,
                                radius=radius)
        ax.text(x + radius * 1.9, y - k * step, label, ha="left", va="center",
                fontsize=fontsize, color=INK2, clip_on=False)


# ───────────────────────────── row ordering ──────────────────────────────────

def order_by_profile(names: list[str], vectors: list[list[int]]) -> list[int]:
    """Order rows so genomes with similar gene content sit together (nearest-
    neighbour chain on Jaccard similarity of the presence vectors). Rows with
    nothing present go last. Deterministic; ties break on name."""
    sets = [frozenset(j for j, v in enumerate(vec) if v) for vec in vectors]
    empty = sorted((i for i, s in enumerate(sets) if not s),
                   key=lambda i: names[i])
    todo = [i for i, s in enumerate(sets) if s]
    if not todo:
        return empty

    def jac(a: frozenset, b: frozenset) -> float:
        u = len(a | b)
        return len(a & b) / u if u else 0.0

    # Start from the row whose first present column is leftmost (so the chain
    # walks the pathway blocks roughly left to right), richest first.
    cur = min(todo, key=lambda i: (min(sets[i]), -len(sets[i]), names[i]))
    order = [cur]
    todo.remove(cur)
    while todo:
        nxt = max(todo, key=lambda i: (jac(sets[cur], sets[i]),
                                       -min(sets[i]), names[i]))
        order.append(nxt)
        todo.remove(nxt)
        cur = nxt
    return order + empty
