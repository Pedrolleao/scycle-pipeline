"""Shared renderer for the complex and process-module completeness grids.

Rows = genomes, columns = complexes / process modules. Each cell is one of
four glyphs (see _viz.draw_completeness_glyph): solid = complete, ring with
"n/N" = partial, faint ring with a cross = ruled out, faint ring = absent.
The fraction is written INSIDE the partial ring, so it can never collide with
a neighbouring row.
"""

from __future__ import annotations

from pathlib import Path
from textwrap import wrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _viz import (GRID, INK, INK2, draw_completeness_glyph,
                  draw_completeness_legend)

CELL = 0.46          # inches per grid step
GUTTER = 0.7         # extra grid steps between column groups


def render_grid(samples: list[str],
                groups: list[tuple[str, list[tuple[str, str]]]],
                cells: dict[str, dict[str, tuple[str, str]]],
                title: str, subtitle: str,
                legend: list[tuple[str, str, str]],
                key: list[tuple[str, str]], key_title: str,
                out: Path) -> None:
    """groups = [(group label, [(column id, column label), …]), …]
    cells[sample][column id] = (state, text)
    legend = [(state, example text, label), …]
    key = [(column label, description), …] printed under the grid."""
    # Rows: most complete first, then most partial, then name.
    def score(s: str) -> tuple:
        states = [v[0] for v in cells.get(s, {}).values()]
        return (-states.count("complete"), -states.count("partial"), s)
    samples = sorted(samples, key=score)
    n_rows = len(samples)

    col_x: dict[str, float] = {}
    labels: list[tuple[float, str]] = []
    spans: list[tuple[str, float, float]] = []
    x = 0.0
    for gi, (glabel, cols) in enumerate(groups):
        if gi:
            x += GUTTER
        x0 = x
        for cid, clabel in cols:
            col_x[cid] = x
            labels.append((x, clabel))
            x += 1.0
        spans.append((glabel, x0, x - 1.0))
    total_w = x

    left_in = 0.3 + max(len(s) for s in samples) * 0.068
    fig_w = left_in + total_w * CELL + 3.6
    fig_h = 1.2 + n_rows * CELL + 3.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(-0.6, total_w - 0.4)
    ax.set_ylim(-0.6, n_rows - 0.4)
    ax.set_aspect("equal")

    for i in range(n_rows):
        ax.plot([-0.5, total_w - 0.5], [i, i], color=GRID, lw=0.5, zorder=0)
    for i, s in enumerate(samples):
        y = n_rows - 1 - i
        for cid, cx in col_x.items():
            state, text = cells.get(s, {}).get(cid, ("absent", ""))
            draw_completeness_glyph(ax, cx, y, state, text)

    # Group labels above the grid.
    for glabel, x0, x1 in spans:
        if not glabel:
            continue
        ax.plot([x0 - 0.4, x1 + 0.4], [n_rows - 0.25] * 2, color=INK2, lw=0.8,
                clip_on=False)
        ax.text(x0 - 0.4, n_rows - 0.05, glabel, ha="left", va="bottom",
                fontsize=8, color=INK2, clip_on=False)

    ax.set_xticks([x for x, _ in labels])
    ax.set_xticklabels([l for _, l in labels], rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=8.5, color=INK)
    ax.xaxis.set_ticks_position("bottom")
    ax.set_yticks([n_rows - 1 - i for i in range(n_rows)])
    ax.set_yticklabels(samples, fontsize=8, color=INK)
    ax.tick_params(axis="both", length=0, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(False)

    head = n_rows + (1.0 if any(g for g, _, _ in spans) else 0.3)
    ax.text(-0.5, head + 0.95, title, ha="left", va="bottom", fontsize=13,
            fontweight="bold", color=INK, clip_on=False)
    ax.text(-0.5, head + 0.35, subtitle, ha="left", va="bottom", fontsize=8,
            color=INK2, clip_on=False)

    draw_completeness_legend(ax, total_w + 0.3, n_rows - 1, legend, step=1.0,
                             fontsize=8)

    # Key under the rotated column labels.
    longest = max(len(l) for _, l in labels)
    y = -0.6 - (longest * 0.066 * 0.71 + 0.35) / CELL
    ax.text(-0.5, y, key_title, ha="left", va="top", fontsize=9,
            fontweight="bold", color=INK, clip_on=False)
    y -= 0.55
    name_w = max(len(k) for k, _ in key) * 0.135 + 0.6      # grid steps
    for name, desc in key:
        lines = wrap(desc, 96) or [""]
        ax.text(-0.5, y, name, ha="left", va="top", fontsize=8, color=INK,
                clip_on=False)
        for line in lines:
            ax.text(-0.5 + name_w, y, line, ha="left", va="top", fontsize=8,
                    color=INK2, clip_on=False)
            y -= 0.40
        y -= 0.08

    out.parent.mkdir(parents=True, exist_ok=True)
    dpi = 300 if out.suffix.lower() == ".png" else 100
    fig.savefig(out, bbox_inches="tight", pad_inches=0.25, dpi=dpi)
    plt.close(fig)
