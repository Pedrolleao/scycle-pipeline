#!/usr/bin/env python3
"""
make_locus_maps.py — draw the genomic neighbourhood of every called gene.

For one nucleotide / MAG sample: cluster the called genes that sit close
together on a contig (see _loci.py), then draw each cluster as a row of gene
arrows to a common bp scale, grouped under the pathway that contributes most
of its calls. Genes outside the target set are drawn blank (they were not
annotated). A bar marks a contig end, i.e. where an operon may simply run off
the assembly.

Outputs:
  --tsv   calls/<x>cycle_loci.tsv one row per called gene, with its locus
  --out   report/loci.svg [.png]  the figure (one or more paths)

Needs gene coordinates, so it only applies to nucleotide input (Prodigal).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import yaml

from _loci import loci_for_sample
from _cycle_model import DISPLAY
from _viz import (AXIS, CALLS_TSV, CAT_ORDER, CYCLE_LETTER, INK, INK2,
                  LOCI_TSV, MUTED, PATHWAY_COLOR, PATHWAY_LABEL, STATUS_SHORT,
                  UNANNOTATED, load_calls)

ARROW_H = 0.17          # inches
LABEL_BAND = 0.46       # room for the rotated gene names
CAPTION_H = 0.20
ROW_GAP = 0.16
SECTION_H = 0.34
LEFT, RIGHT = 0.45, 1.0
SUFFIX = {2: "", 1: "°", -1: "×"}


def write_tsv(loci: list[dict], path: Path) -> None:
    cols = ["locus_id", "pathway", "contig", "locus_start", "locus_end",
            "contig_edge", "protein_id", "gene_start", "gene_end", "strand",
            "target_id", "status", "copy", "other_calls"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for l in loci:
            edge = ("both" if l["edge_left"] and l["edge_right"] else
                    "left" if l["edge_left"] else
                    "right" if l["edge_right"] else "")
            for g in l["genes"]:
                if not g["calls"]:
                    continue
                best = g["calls"][0]
                other = ";".join(f"{c['target']}:{STATUS_SHORT[c['code']]}"
                                 for c in g["calls"][1:])
                fh.write("\t".join(str(v) for v in [
                    l["locus_id"], l["pathway"], l["contig"], l["start"],
                    l["end"], edge, g["id"], g["start"], g["end"], g["strand"],
                    best["target"], STATUS_SHORT[best["code"]],
                    "additional" if best["extra"] else "reported", other]) + "\n")


def gene_arrow(x0: float, x1: float, y: float, strand: str) -> list[tuple]:
    """Polygon for a gene arrow between x0 < x1 (inches), centred on y."""
    h = ARROW_H / 2
    head = min(0.09, (x1 - x0) * 0.6)
    if strand == "+":
        return [(x0, y - h), (x1 - head, y - h), (x1, y),
                (x1 - head, y + h), (x0, y + h)]
    return [(x1, y - h), (x0 + head, y - h), (x0, y),
            (x0 + head, y + h), (x1, y + h)]


def draw_gene(ax, x0, x1, y, strand, call) -> None:
    poly = gene_arrow(x0, x1, y, strand)
    if call is None:
        ax.add_patch(mpatches.Polygon(poly, closed=True, facecolor=UNANNOTATED,
                                      edgecolor="none", zorder=3))
        return
    hue = PATHWAY_COLOR.get(call["category"], MUTED)
    code = call["code"]
    if code == 2:
        ax.add_patch(mpatches.Polygon(poly, closed=True, facecolor=hue,
                                      edgecolor=hue, linewidth=0.5, zorder=4))
    elif code == 1:
        ax.add_patch(mpatches.Polygon(poly, closed=True, facecolor="white",
                                      edgecolor="none", zorder=4))
        ax.add_patch(mpatches.Polygon(poly, closed=True, facecolor=hue,
                                      alpha=0.30, edgecolor="none", zorder=4.1))
        ax.add_patch(mpatches.Polygon(poly, closed=True, facecolor="none",
                                      edgecolor=hue, linewidth=1.0, zorder=4.2))
    else:
        ax.add_patch(mpatches.Polygon(poly, closed=True, facecolor="white",
                                      edgecolor=MUTED, linewidth=0.9, zorder=4))


def render(sample: str, loci: list[dict], outs: list[Path], note: str = "") -> None:
    shown = [l for l in loci if l["n_present"] > 0]
    n_hidden = len(loci) - len(shown)
    sections = [(c, [l for l in shown if l["pathway"] == c]) for c in CAT_ORDER]
    sections = [(c, ls) for c, ls in sections if ls]

    max_span = max((l["end"] - l["start"] for l in shown), default=10000)
    bp_per_in = max(2500.0, max_span / 11.0)
    plot_w = max(6.2, max_span / bp_per_in)
    title = f"{sample} — {CYCLE_LETTER}-cycle gene neighbourhoods"
    fig_w = max(LEFT + plot_w + RIGHT, LEFT + len(title) * 0.108 + 0.3)
    row_h = CAPTION_H + LABEL_BAND + ARROW_H + ROW_GAP
    top_h, bottom_h = 1.1, 1.05
    body_h = sum(SECTION_H + len(ls) * row_h for _, ls in sections) or 0.6
    fig_h = top_h + body_h + bottom_h

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, fig_h)
    ax.axis("off")

    n_genes = sum(1 for l in shown for g in l["genes"]
                  if g["calls"] and g["calls"][0]["code"] in (1, 2))
    ax.text(LEFT, fig_h - 0.32, title,
            fontsize=12.5, fontweight="bold", color=INK, va="center")
    sub = (f"{len(shown)} {'locus' if len(shown) == 1 else 'loci'} · "
           f"{n_genes} called genes · "
           "blank arrows are genes outside the target set (not annotated)")
    ax.text(LEFT, fig_h - 0.58, note or sub, fontsize=7.5, color=INK2, va="center")
    if n_hidden and not note:
        ax.text(LEFT, fig_h - 0.76,
                f"{n_hidden} more {'locus has' if n_hidden == 1 else 'loci have'} "
                f"only disqualified hits and are not drawn (listed in {LOCI_TSV})",
                fontsize=7.5, color=INK2, va="center")

    y = fig_h - top_h
    if not shown and not note:
        ax.text(LEFT, y - 0.3, f"No {CYCLE_LETTER}-cycle genes were called in this sample.",
                fontsize=9, color=INK2, va="center")
    for cat, ls in sections:
        hue = PATHWAY_COLOR[cat]
        ax.add_patch(mpatches.Rectangle((LEFT, y - 0.21), 0.07, 0.17,
                                        facecolor=hue, edgecolor="none"))
        ax.text(LEFT + 0.15, y - 0.125, PATHWAY_LABEL[cat], fontsize=9.5,
                fontweight="bold", color=INK, va="center")
        ax.plot([LEFT, fig_w - 0.3], [y - 0.27, y - 0.27], color="#e1e0d9",
                lw=0.6, zorder=1)
        y -= SECTION_H
        for l in ls:
            cap = (f"{l['locus_id']} · {l['contig']} · "
                   f"{l['start']:,}–{l['end']:,} bp")
            ax.text(LEFT, y - CAPTION_H / 2, cap, fontsize=7, color=MUTED,
                    va="center")
            yc = y - CAPTION_H - LABEL_BAND - ARROW_H / 2
            x_of = lambda pos: LEFT + (pos - l["start"]) / bp_per_in   # noqa: E731
            xa, xb = x_of(l["start"]), x_of(l["end"])
            # Backbone: runs past the window unless the contig ends there.
            ext = 0.14
            if l["edge_left"]:
                xl = max(xa - 0.30, x_of(1)) if l["contig_len"] else xa - 0.1
            else:
                xl = xa - ext
            if l["edge_right"]:
                xr = (min(xb + 0.30, x_of(l["contig_len"]))
                      if l["contig_len"] else xb + 0.1)
            else:
                xr = xb + ext
            ax.plot([xl, xr], [yc, yc], color=AXIS, lw=0.9, zorder=2,
                    solid_capstyle="butt")
            for is_edge, xe in ((l["edge_left"], xl), (l["edge_right"], xr)):
                if is_edge:
                    ax.plot([xe, xe], [yc - 0.14, yc + 0.14], color=INK2,
                            lw=1.6, zorder=5, solid_capstyle="butt")
            for g in l["genes"]:
                x0, x1 = x_of(g["start"]), x_of(g["end"])
                call = g["calls"][0] if g["calls"] else None
                draw_gene(ax, x0, x1, yc, g["strand"], call)
                if call:
                    name = DISPLAY.get(call["target"], call["target"]) + SUFFIX[call["code"]]
                    ax.text((x0 + x1) / 2, yc + ARROW_H / 2 + 0.04, name,
                            rotation=40, rotation_mode="anchor", ha="left",
                            va="bottom", fontsize=7,
                            color=INK if call["code"] in (1, 2) else MUTED,
                            fontweight="bold" if call["code"] in (1, 2) else "normal",
                            zorder=6)
            y -= row_h

    # Legend + scale bar.
    yl = 0.62
    x = LEFT
    demo = "#6b6a66"
    items = [("confirmed", {"code": 2}), ("domain-only°", {"code": 1}),
             ("disqualified×", {"code": -1}), ("not annotated", None)]
    for lab, call in items:
        if call is not None:
            call = dict(call, category="__demo__")
        poly_call = call
        old = PATHWAY_COLOR.get("__demo__")
        PATHWAY_COLOR["__demo__"] = demo
        draw_gene(ax, x, x + 0.42, yl, "+", poly_call)
        if old is None:
            PATHWAY_COLOR.pop("__demo__")
        ax.text(x + 0.50, yl, lab, fontsize=7.5, color=INK2, va="center")
        x += 0.50 + len(lab) * 0.062 + 0.30
    ax.plot([x, x], [yl - 0.14, yl + 0.14], color=INK2, lw=1.6,
            solid_capstyle="butt")
    ax.text(x + 0.10, yl, "contig end", fontsize=7.5, color=INK2, va="center")
    # Scale bar: 1 kb (or 5 kb when 1 kb would be tiny).
    kb = 1000 if 1000 / bp_per_in >= 0.25 else 5000
    xs = LEFT
    ax.plot([xs, xs + kb / bp_per_in], [0.28, 0.28], color=INK2, lw=1.2,
            solid_capstyle="butt")
    ax.text(xs + kb / bp_per_in + 0.08, 0.28,
            f"{kb // 1000} kb  ·  arrow fill = pathway of the called gene; "
            "arrow direction = strand",
            fontsize=7.5, color=INK2, va="center")

    for out in outs:
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300 if out.suffix.lower() == ".png" else 100)
        print(f"[locus_maps] wrote {out}")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--results-dir", required=True, type=Path)
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--tsv", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path, nargs="+")
    ap.add_argument("--max-gap", type=int, default=5,
                    help="called genes separated by at most this many other "
                         "genes join one locus")
    ap.add_argument("--flank", type=int, default=2,
                    help="uncalled genes drawn on each side of a locus")
    args = ap.parse_args()

    targets = yaml.safe_load(open(args.targets))["targets"]
    calls = load_calls(args.results_dir / args.sample / "calls" / CALLS_TSV)
    loci = loci_for_sample(args.results_dir, args.sample, calls, targets,
                           args.max_gap, args.flank)
    has_coords = (args.results_dir / args.sample / "prodigal"
                  / f"{args.sample}.faa").exists()
    write_tsv(loci, args.tsv)
    note = "" if has_coords else (
        "No gene coordinates: locus maps need nucleotide / MAG input "
        "(this sample was a pre-called proteome).")
    render(args.sample, loci, args.out, note)
    print(f"[locus_maps] {args.sample}: {len(loci)} loci → {args.tsv}")


if __name__ == "__main__":
    main()
