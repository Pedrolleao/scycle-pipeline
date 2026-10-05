#!/usr/bin/env python3
"""
make_cycle_map.py — draw each genome's calls onto the element cycle.

The diagram (compounds, reaction steps, layout) is defined in _cycle_model.py;
this script is identical in the nitrogen and sulfur sister pipelines. Every
reaction arrow is coloured by its pathway
and drawn by state — complete (solid, heavy), partial (dashed), absent (thin
grey) — and the genes behind it are written next to it: bold when found,
grey when not.

Two uses:
  one genome  → --sample NAME   --out results/NAME/report/<x>cycle_map.svg [.png]
  all genomes → --samples A,B,… --out results/figures/<x>cycle_maps.svg [.png]
                (small multiples; above --max-panels genomes the grid is
                 replaced by a pointer to the per-genome maps and report.html)

--out takes one or more paths; the format follows each extension.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath

import _cycle_model as model
from _completeness import pretty
from _viz import (AXIS, CALLS_TSV, CAT_ORDER, CYCLE_LETTER, CYCLE_NAME, INK,
                  INK2, MUTED, PATHWAY_COLOR, PATHWAY_LABEL, call_codes,
                  load_calls)

ABSENT_EDGE = "#d5d3cc"
ABSENT_TEXT = "#aeaca4"
UNITS_W = model.XLIM[1] - model.XLIM[0]
UNITS_H = model.YLIM[1] - model.YLIM[0]


def complete_processes(path: Path) -> list[str]:
    """Names of the process modules scored `complete` for this genome."""
    if not path.exists():
        return []
    with open(path) as fh:
        return [pretty(r["synergy_id"])
                for r in csv.DictReader(fh, delimiter="\t")
                if r.get("status") == "complete"]


def _token_style(code: int) -> tuple[str, str, str]:
    """(suffix, colour, weight) for a gene name."""
    if code == 2:
        return "", INK, "bold"
    if code == 1:
        return "°", INK, "bold"
    if code == -1:
        return "×", MUTED, "normal"
    return "", ABSENT_TEXT, "normal"


def _draw_label(ax, renderer, x, y, ha, lines, fontsize):
    """Write the gene tokens of one step, line by line, each token styled by
    its call. Tokens are separate Text objects, so widths are measured."""
    inv = ax.transData.inverted()
    gap = 1.3
    for k, line in enumerate(lines):
        yy = y - k * model.LINE_STEP
        texts, widths = [], []
        for tok in line:
            suffix, colour, weight = _token_style(tok["code"])
            t = ax.text(0, yy, tok["label"] + suffix, fontsize=fontsize,
                        color=colour, fontweight=weight, ha="left",
                        va="center", zorder=6)
            bb = t.get_window_extent(renderer)
            (x0, _), (x1, _) = inv.transform([(bb.x0, 0), (bb.x1, 0)])
            texts.append(t)
            widths.append(x1 - x0)
        total = sum(widths) + gap * (len(widths) - 1)
        cx = x - total if ha == "right" else x - total / 2 if ha == "center" else x
        for t, w in zip(texts, widths):
            t.set_x(cx)
            cx += w + gap


def draw_cycle(ax, renderer, codes: dict[str, int], fontsize: float,
               title: str = "", subtitle: str = "",
               ctx: dict | None = None) -> None:
    ax.set_xlim(*model.XLIM)
    ax.set_ylim(*model.YLIM)
    ax.set_aspect("equal")
    ax.axis("off")
    lay = model.layout()
    ev = model.evaluate(codes, ctx)

    # Arrows — absent first so lit arrows sit on top.
    order = {"absent": 0, "partial": 1, "complete": 2}
    for st in sorted(lay["steps"], key=lambda s: order[ev[s["id"]]["state"]]):
        state = ev[st["id"]]["state"]
        hue = PATHWAY_COLOR[st["pathway"]]
        colour = ABSENT_EDGE if state == "absent" else hue
        lw = {"complete": 2.6, "partial": 1.5, "absent": 0.9}[state]
        ls = (0, (3.2, 2.2)) if state == "partial" else "solid"
        z = 2 + order[state]
        arrows = st["arrows_rev"] if ev[st["id"]]["reversed"] else st["arrows"]
        for a in arrows:
            if not a["line"]:
                continue
            verts = a["line"]
            ax.add_patch(mpatches.PathPatch(
                MPath(verts, [MPath.MOVETO] + [MPath.LINETO] * (len(verts) - 1)),
                facecolor="none", edgecolor=colour, linewidth=lw, linestyle=ls,
                capstyle="butt", joinstyle="round", zorder=z))
            ax.add_patch(mpatches.Polygon(a["head"], closed=True,
                                          facecolor=colour, edgecolor=colour,
                                          linewidth=0.4, zorder=z))

    # Compounds.
    for n in lay["nodes"]:
        ax.add_patch(mpatches.FancyBboxPatch(
            (n["x"] - n["w"] / 2 + 1.2, n["y"] - n["h"] / 2 + 1.2),
            n["w"] - 2.4, n["h"] - 2.4,
            boxstyle="round,pad=1.2,rounding_size=2.4",
            facecolor="white", edgecolor=AXIS, linewidth=0.9, zorder=5))
        ax.text(n["x"], n["y"], n["label"], ha="center", va="center",
                fontsize=fontsize + 1.5, color=INK, zorder=6)

    # Gene labels.
    for st in lay["steps"]:
        x, y, ha = st["label"]
        _draw_label(ax, renderer, x, y, ha, ev[st["id"]]["lines"], fontsize)

    if title:
        ax.text(model.XLIM[0] + 1, model.YLIM[1] + 7.5, title, ha="left",
                va="bottom", fontsize=fontsize + 3.5, fontweight="bold",
                color=INK, clip_on=False)
    if subtitle:
        ax.text(model.XLIM[0] + 1, model.YLIM[1] + 2.5, subtitle, ha="left",
                va="bottom", fontsize=fontsize, color=INK2, clip_on=False)


def draw_legend(ax, fontsize: float) -> None:
    """Legend strip in its own axes (data units = inches-ish, y up)."""
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 20)
    ax.axis("off")
    # Arrow states.
    y = 15.5
    demo = "#6b6a66"
    for k, (lab, lw, ls, colour) in enumerate([
            ("complete step — every gene of one route found", 2.6, "solid", demo),
            ("partial — some route genes found", 1.5, (0, (3.2, 2.2)), demo),
            ("absent", 0.9, "solid", ABSENT_EDGE)]):
        yy = y - k * 4.6
        ax.plot([0, 6.5], [yy, yy], color=colour, lw=lw, ls=ls,
                solid_capstyle="butt")
        ax.add_patch(mpatches.Polygon([(8.4, yy), (6.3, yy + 1.1), (6.3, yy - 1.1)],
                                      closed=True, facecolor=colour,
                                      edgecolor=colour, linewidth=0.3))
        ax.text(10, yy, lab, ha="left", va="center", fontsize=fontsize,
                color=INK2)
    # Gene-name styles.
    x0 = 56
    for k, (sample, colour, weight, lab) in enumerate([
            ("gene", INK, "bold", "confirmed"),
            ("gene°", INK, "bold", "domain-only (HMM signature, no BLAST support)"),
            ("gene×", MUTED, "normal", "disqualified (failed the homology-trap gate)"),
            ("gene", ABSENT_TEXT, "normal", "absent")]):
        yy = 17.0 - k * 4.0
        ax.text(x0, yy, sample, ha="left", va="center", fontsize=fontsize,
                color=colour, fontweight=weight)
        ax.text(x0 + 6.5, yy, lab, ha="left", va="center", fontsize=fontsize,
                color=INK2)


def draw_pathway_key(ax, fontsize: float) -> None:
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 15)
    ax.axis("off")
    per_row = 2 if max(len(PATHWAY_LABEL[c]) for c in CAT_ORDER) > 36 else 3
    for i, cat in enumerate(CAT_ORDER):
        col, row = i % per_row, i // per_row
        x, y = col * (102.0 / per_row), 12.5 - row * 4.6
        ax.plot([x, x + 3.2], [y, y], color=PATHWAY_COLOR[cat], lw=2.6,
                solid_capstyle="butt")
        ax.text(x + 4.2, y, PATHWAY_LABEL[cat], ha="left", va="center",
                fontsize=fontsize, color=INK2)


def render(samples: list[str], results_dir: Path, outs: list[Path],
           max_panels: int) -> None:
    n = len(samples)
    single = n == 1
    if not single and n > max_panels:
        fig = plt.figure(figsize=(8, 2))
        fig.text(0.5, 0.5,
                 f"{n} genomes — too many for one page of {CYCLE_LETTER}-cycle maps "
                 f"(limit {max_panels}).\nSee the per-sample report/ folders "
                 "for each genome, or open report.html.",
                 ha="center", va="center", fontsize=10, color=INK2)
        _save(fig, outs)
        return

    F = model.FIG
    ncols = 1 if single else min(F["grid_cols"], n)
    nrows = -(-n // ncols)
    panel_w = F["single_w"] if single else F["grid_w"]
    fontsize = F["single_fs"] if single else F["grid_fs"]
    unit = panel_w / UNITS_W                   # inches per data unit
    panel_h = UNITS_H * unit
    head_h = 14 * unit                         # title + subtitle band
    gap_x, gap_y = 0.25, 0.2
    margin = 0.3
    legend_h, key_h = 1.05, 0.8
    top_h = 0.0 if single else 0.55            # figure title (grid only)

    fig_w = 2 * margin + ncols * panel_w + (ncols - 1) * gap_x
    fig_h = (2 * margin + top_h + nrows * (panel_h + head_h)
             + (nrows - 1) * gap_y + legend_h + key_h + 0.2)
    fig = plt.figure(figsize=(fig_w, fig_h))
    renderer = fig.canvas.get_renderer()

    if not single:
        fig.text(margin / fig_w, 1 - (margin + 0.05) / fig_h,
                 f"{CYCLE_NAME.capitalize()}-cycle maps — {n} genomes",
                 ha="left", va="top", fontsize=14, fontweight="bold", color=INK)

    for i, s in enumerate(samples):
        col, row = i % ncols, i // ncols
        left = margin + col * (panel_w + gap_x)
        top = margin + top_h + row * (panel_h + head_h + gap_y) + head_h
        ax = fig.add_axes([left / fig_w, 1 - (top + panel_h) / fig_h,
                           panel_w / fig_w, panel_h / fig_h])
        calls = load_calls(results_dir / s / "calls" / CALLS_TSV)
        done = complete_processes(results_dir / s / "calls" / "synergy_completeness.tsv")
        ctx = model.genome_context(calls)
        sub = ("complete: " + " · ".join(done)) if done else "no complete process module"
        limit = int(98 * UNITS_W / 117)
        if len(sub) > limit:
            sub = sub[:limit - 3].rsplit(" · ", 1)[0] + " · …"
        if model.context_note(ctx):
            sub = model.context_note(ctx) + "  ·  " + sub
        draw_cycle(ax, renderer, call_codes(calls), fontsize, title=s,
                   subtitle=sub, ctx=ctx)

    leg_w = min(fig_w - 2 * margin, 7.4)
    ax_leg = fig.add_axes([margin / fig_w, (margin + key_h + 0.1) / fig_h,
                           leg_w / fig_w, legend_h / fig_h])
    draw_legend(ax_leg, 8.0)
    ax_key = fig.add_axes([margin / fig_w, margin / fig_h,
                           leg_w / fig_w, key_h / fig_h])
    draw_pathway_key(ax_key, 8.0)
    _save(fig, outs)


def _save(fig, outs: list[Path]) -> None:
    for out in outs:
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300 if out.suffix.lower() == ".png" else 100)
        print(f"[cycle_map] wrote {out}")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True, type=Path)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--sample", help="one genome → single map")
    g.add_argument("--samples", help="comma-separated → small multiples")
    ap.add_argument("--max-panels", type=int, default=48)
    ap.add_argument("--out", required=True, type=Path, nargs="+")
    args = ap.parse_args()

    samples = [args.sample] if args.sample else [s for s in args.samples.split(",") if s]
    render(samples, args.results_dir, args.out, args.max_panels)


if __name__ == "__main__":
    main()
