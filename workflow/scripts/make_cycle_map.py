#!/usr/bin/env python3
"""
make_cycle_map.py — draw each genome's calls onto the element cycle, as a
metro map.

The diagram (stations = compounds, lines = pathways, stops = genes, layout) is
defined in _cycle_model.py; this script is shared with the nitrogen and sulfur
sister pipelines once their models carry a metro layout (2026-11-08: mcycle
only). Every line is coloured by its pathway and drawn by the state
of its step — complete (solid), partial (hatched: the line at 40 % with
full-colour stripes), absent (grey) — and the genes behind it are the stops on the line: a white disc ringed
in the line colour when found (half-marked when domain-only), faint when
absent, crossed when disqualified. A chevron near the arrival station gives the
direction of travel; in a reverse-Mcr (ANME) genome the methanogenesis chain
runs the other way.

Two uses:
  one genome  → --sample NAME   --out results/NAME/report/<x>cycle_map.svg [.png]
  all genomes → --samples A,B,… --out results/figures/<x>cycle_maps.svg [.png]
                (small multiples; above --max-panels genomes the grid is
                 replaced by a pointer to the per-genome maps and the HTML report)

--out takes one or more paths; the format follows each extension.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt

import _cycle_model as model
from _domain import REPORT_HTML
from _completeness import pretty
from _viz import (CALLS_TSV, CAT_ORDER, CYCLE_LETTER, CYCLE_NAME, INK, INK2,
                  MUTED, PATHWAY_COLOR, PATHWAY_LABEL, call_codes, load_calls)

BG = "white"
ABSENT_EDGE = "#d5d3cc"
ABSENT_TEXT = "#aeaca4"
UNITS_W = model.XLIM[1] - model.XLIM[0]
UNITS_H = model.YLIM[1] - model.YLIM[0]
OFF_MARK = "∅"
HA = {"left": "left", "right": "right", "center": "center"}
VA = {"top": "top", "bottom": "bottom", "center": "center"}


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


def _halo(u2pt: float):
    return [pe.withStroke(linewidth=0.5 * u2pt, foreground=BG)]


def _line(ax, pts, colour, state, u2pt, z, stripes=None):
    """A track. A partial step is the same continuous line at 40 % opacity
    with full-colour 45° stripes across it (model.stripes)."""
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color=colour, lw=model.LINE_W * u2pt, solid_capstyle="round",
            solid_joinstyle="round", zorder=z,
            alpha=model.PARTIAL_ALPHA if state == "partial" else 1.0)
    if state == "partial":
        for q in (stripes if stripes is not None else model.stripes(pts)):
            ax.add_patch(mpatches.Polygon(q, closed=True, facecolor=colour,
                                          edgecolor="none", zorder=z + 0.05))


def _chevron(ax, c, u2pt, z):
    """A light V inside the line, pointing along the travel."""
    L = model.CHEVRON
    tx, ty = c["tx"], c["ty"]
    nx, ny = -ty, tx
    tip = (c["x"], c["y"])
    a1 = (tip[0] - (tx - nx * 0.9) * L, tip[1] - (ty - ny * 0.9) * L)
    a2 = (tip[0] - (tx + nx * 0.9) * L, tip[1] - (ty + ny * 0.9) * L)
    ax.plot([a1[0], tip[0], a2[0]], [a1[1], tip[1], a2[1]], color=BG,
            lw=0.45 * u2pt, solid_capstyle="round", solid_joinstyle="round",
            alpha=0.95, zorder=z)


def _stop(ax, sp, code, hue, u2pt, z):
    ring = hue if code in (1, 2) else MUTED if code == -1 else ABSENT_EDGE
    r = model.STOP_R * (0.8 if code == 0 else 1.0)
    ax.add_patch(mpatches.Circle((sp["x"], sp["y"]), r, facecolor=BG,
                                 edgecolor=ring,
                                 linewidth=(0.35 if code == 0 else 0.5) * u2pt,
                                 zorder=z))
    if code == 1:
        ax.add_patch(mpatches.Circle((sp["x"], sp["y"]), 0.38, facecolor=hue,
                                     edgecolor="none", zorder=z + 0.1))
    if code == -1:
        q = 0.5
        ax.plot([sp["x"] - q, sp["x"] + q], [sp["y"] - q, sp["y"] + q],
                color=MUTED, lw=0.3 * u2pt, zorder=z + 0.1)
        ax.plot([sp["x"] - q, sp["x"] + q], [sp["y"] + q, sp["y"] - q],
                color=MUTED, lw=0.3 * u2pt, zorder=z + 0.1)


def _stop_label(ax, sp, tok, fontsize, u2pt, z, off=False):
    suffix, colour, weight = _token_style(tok["code"])
    text = tok["label"] + suffix
    if off:
        text += " " + OFF_MARK
    ax.text(sp["lx"], sp["ly"], text, fontsize=fontsize, color=colour,
            fontweight=weight, ha=HA[sp["ha"]], va=VA[sp["va"]],
            rotation=sp["rot"], rotation_mode="anchor", zorder=z,
            path_effects=_halo(u2pt))


def draw_cycle(ax, codes: dict[str, int], fontsize: float, u2pt: float,
               title: str = "", subtitle: str = "",
               ctx: dict | None = None) -> None:
    """One genome's metro map. u2pt = points per data unit (line widths)."""
    ax.set_xlim(*model.XLIM)
    ax.set_ylim(*model.YLIM)
    ax.set_aspect("equal")
    ax.axis("off")
    lay = model.layout()
    ev = model.evaluate(codes, ctx)

    order = {"absent": 0, "partial": 1, "complete": 2}
    for st in sorted(lay["steps"], key=lambda s: order[ev[s["id"]]["state"]]):
        e = ev[st["id"]]
        state = e["state"]
        hue = PATHWAY_COLOR[st["pathway"]]
        colour = ABSENT_EDGE if state == "absent" else hue
        z = 1 + order[state]
        pts = list(reversed(st["path"])) if e["reversed"] else st["path"]
        _line(ax, pts, colour, state, u2pt, z, stripes=st["stripes"])
        _chevron(ax, st["chevron_rev"] if e["reversed"] else st["chevron"],
                 u2pt, z + 0.1)
        for fd in st["feeders"]:
            _line(ax, fd["path"], colour, state, u2pt, z, stripes=fd["stripes"])
            _chevron(ax, fd["chevron"], u2pt, z + 0.1)
        toks = [t for line in e["lines"] for t in line]
        for i, (sp, tok) in enumerate(zip(st["stops"], toks)):
            _stop(ax, sp, tok["code"], hue, u2pt, 5)
            _stop_label(ax, sp, tok, fontsize, u2pt, 6,
                        off=bool(e.get("off")) and i == 0)

    for n in lay["nodes"]:
        ax.add_patch(mpatches.Circle((n["x"], n["y"]), model.STATION_R,
                                     facecolor=BG, edgecolor=INK,
                                     linewidth=0.75 * u2pt, zorder=7))
        ax.text(n["x"] + n["ldx"], n["y"] + n["ldy"], n["label"],
                ha=HA[n["ha"]], va=VA[n["va"]], fontsize=fontsize + 0.5,
                fontweight="bold", color=INK, zorder=8,
                path_effects=_halo(u2pt))

    if title:
        ax.text(model.XLIM[0] + 1, model.YLIM[1] + 7.5, title, ha="left",
                va="bottom", fontsize=fontsize + 3.5, fontweight="bold",
                color=INK, clip_on=False)
    if subtitle:
        ax.text(model.XLIM[0] + 1, model.YLIM[1] + 2.5, subtitle, ha="left",
                va="bottom", fontsize=fontsize, color=INK2, clip_on=False)


# ───────────────────────────── keys ──────────────────────────────────────────

def _key_items(fontsize: float):
    """Rows of the state / glyph key: (kind, args, label)."""
    return [
        ("line", ("complete", "#6b6a66"), "complete step — every gene of one route found"),
        ("line", ("partial", "#6b6a66"), "partial — some route genes found (hatched)"),
        ("line", ("absent", ABSENT_EDGE), "absent"),
        ("stop", (2, "#6b6a66"), "gene found (confirmed)"),
        ("stop", (1, "#6b6a66"), "gene° domain-only (HMM signature, no BLAST support)"),
        ("stop", (-1, "#6b6a66"), "gene× disqualified (failed the homology-trap gate)"),
        ("stop", (0, "#6b6a66"), "gene absent"),
        ("text", (OFF_MARK,), "step off for this genome type; its genes keep their calls"),
    ]


def _key_row(ax, x, y, kind, args, label, fontsize, u2pt, scale=1.0):
    """One key row at (x, y) in data units; the sample glyph spans 8 units."""
    if kind == "line":
        state, colour = args
        _line(ax, [(x, y), (x + 8.0 * scale, y)], colour, state, u2pt * scale, 2)
        _chevron(ax, {"x": x + 7.2 * scale, "y": y, "tx": 1.0, "ty": 0.0},
                 u2pt * scale, 2.1)
    elif kind == "stop":
        code, hue = args
        _line(ax, [(x, y), (x + 8.0 * scale, y)], ABSENT_EDGE if code == 0 else hue,
              "complete", u2pt * scale, 2)
        _stop(ax, {"x": x + 4.0 * scale, "y": y}, code, hue, u2pt * scale, 5)
    else:
        ax.text(x + 4.0 * scale, y, args[0], ha="center", va="center",
                fontsize=fontsize, color=ABSENT_TEXT)
    ax.text(x + 9.5 * scale, y, label, ha="left", va="center",
            fontsize=fontsize, color=INK2)


def _key_axes(ax, aspect: float) -> float:
    """100 data units across; the height follows the axes box so circles stay
    round under an equal aspect. Returns the height in data units."""
    H = 100.0 * aspect
    ax.set_xlim(0, 100)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    return H


def draw_legend(ax, fontsize: float, u2pt: float, aspect: float) -> None:
    """Glyph key in its own axes below a grid (two columns of four rows)."""
    H = _key_axes(ax, aspect)
    step = H / 4.4
    for i, (kind, args, label) in enumerate(_key_items(fontsize)):
        col, row = i // 4, i % 4
        _key_row(ax, col * 50.0, H - step * (row + 0.6), kind, args, label,
                 fontsize, u2pt)


def draw_pathway_key(ax, fontsize: float, u2pt: float, aspect: float) -> None:
    H = _key_axes(ax, aspect)
    per_row = 2 if max(len(PATHWAY_LABEL[c]) for c in CAT_ORDER) > 30 else 3
    nrows = -(-len(CAT_ORDER) // per_row)
    step = H / (nrows + 0.4)
    for i, cat in enumerate(CAT_ORDER):
        col, row = i % per_row, i // per_row
        x, y = col * (100.0 / per_row), H - step * (row + 0.6)
        _line(ax, [(x, y), (x + 8.0, y)], PATHWAY_COLOR[cat], "complete", u2pt, 2)
        ax.text(x + 9.5, y, PATHWAY_LABEL[cat], ha="left", va="center",
                fontsize=fontsize, color=INK2)


# ───────────────────────────── figure ────────────────────────────────────────

def render(samples: list[str], results_dir: Path, outs: list[Path],
           max_panels: int) -> None:
    n = len(samples)
    single = n == 1
    if not single and n > max_panels:
        fig = plt.figure(figsize=(8, 2))
        fig.text(0.5, 0.5,
                 f"{n} genomes — too many for one page of {CYCLE_LETTER}-cycle maps "
                 f"(limit {max_panels}).\nSee the per-sample report/ folders "
                 f"for each genome, or open {REPORT_HTML}.",
                 ha="center", va="center", fontsize=10, color=INK2)
        _save(fig, outs)
        return

    F = model.FIG
    ncols = 1 if single else min(F["grid_cols"], n)
    nrows = -(-n // ncols)
    panel_w = F["single_w"] if single else F["grid_w"]
    fontsize = F["single_fs"] if single else F["grid_fs"]
    unit = panel_w / UNITS_W                   # inches per data unit
    u2pt = 72.0 * unit                         # points per data unit
    panel_h = UNITS_H * unit
    head_h = 14 * unit                         # title + subtitle band
    gap_x, gap_y = 0.25, 0.2
    margin = 0.3
    legend_h, key_h = 1.05, 1.0
    top_h = 0.0 if single else 0.55            # figure title (grid only)

    fig_w = 2 * margin + ncols * panel_w + (ncols - 1) * gap_x
    fig_h = (2 * margin + top_h + nrows * (panel_h + head_h)
             + (nrows - 1) * gap_y + legend_h + key_h + 0.2)
    fig = plt.figure(figsize=(fig_w, fig_h))

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
        draw_cycle(ax, call_codes(calls), fontsize, u2pt, title=s,
                   subtitle=sub, ctx=ctx)

    leg_w = min(fig_w - 2 * margin, 9.0)
    leg_u2pt = 72.0 * leg_w / 100.0
    ax_leg = fig.add_axes([margin / fig_w, (margin + key_h + 0.1) / fig_h,
                           leg_w / fig_w, legend_h / fig_h])
    draw_legend(ax_leg, 8.0, leg_u2pt, legend_h / leg_w)
    ax_key = fig.add_axes([margin / fig_w, margin / fig_h,
                           leg_w / fig_w, key_h / fig_h])
    draw_pathway_key(ax_key, 8.0, leg_u2pt, key_h / leg_w)
    _save(fig, outs)


def _save(fig, outs: list[Path]) -> None:
    for out in outs:
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300 if out.suffix.lower() == ".png" else 100,
                    facecolor=BG)
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
