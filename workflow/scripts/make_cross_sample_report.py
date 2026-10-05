#!/usr/bin/env python3
"""
make_cross_sample_report.py — assemble cross-sample matrix.tsv and heatmap.svg.

Matrix columns:
  per-target status codes:   0=absent, 1=domain-only, 2=confirmed, -1=disqualified
  per-complex completeness:  float 0–1
  per-synergy completeness:  float 0–1

Heatmap: a dot grid — rows = samples (ordered by gene-content similarity),
columns = genes blocked by pathway. Hue = pathway, glyph shape = call status.
Complexes / process modules have their own figures.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import yaml

from _viz import (CALLS_TSV, CAT_ORDER, CYCLE_LETTER, GRID, INK, INK2, MUTED, PATHWAY_COLOR,
                  PATHWAY_SHORT, PRESENT_CODES, STATUS_CODE,
                  draw_status_glyph, draw_status_legend, order_by_profile)



def load_calls(path: Path) -> dict[str, str]:
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["target_id"]] = row["status"]
    return out


def load_completeness(path: Path) -> dict[str, float]:
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = row.get("complex_id") or row.get("synergy_id")
            out[key] = float(row["completeness"])
    return out


# ───────────────────────── overview dot-grid ─────────────────────────────────
# Palette and status glyphs live in _viz.py (shared with every other figure).

def render_circle_grid(samples, ordered_targets, cat_boundaries,
                       target_codes, out) -> None:
    """Render the cross-sample gene grid as a publication-style dot grid.

    Genes run left→right grouped by pathway (colour band + label above each
    block); genomes are rows, ordered so that similar gene content sits
    together. Status is carried by glyph shape (see _viz.draw_status_glyph).
    Co-renders to whichever extension `out` carries (PNG 300 DPI, SVG vector).
    """
    order = order_by_profile(
        samples, [[1 if v in PRESENT_CODES else 0 for v in row]
                  for row in target_codes.tolist()])
    samples = [samples[i] for i in order]
    target_codes = target_codes[order, :]

    n_rows = len(samples)
    n_genes = len(ordered_targets)
    gutter = 0.6   # gap (in cell units) between pathway blocks

    # x-position per gene, with a gutter between blocks.
    gene_x: list[float] = []
    x = 0.0
    block_spans: list[tuple[str, float, float]] = []  # (cat, x_start, x_end)
    for bi, (cat, start, end) in enumerate(cat_boundaries):
        if bi > 0:
            x += gutter
        x_start = x
        for _ in range(start, end):
            gene_x.append(x)
            x += 1.0
        block_spans.append((cat, x_start, x - 1.0))
    total_w = x

    # ── layout (inches) ──────────────────────────────────────────────────
    cell = 0.30                          # inches per grid step
    head_fs, head_rot = 8.5, 32          # angled block labels
    longest = max(len(PATHWAY_SHORT[c]) for c, _, _ in cat_boundaries)
    head_rise = longest * head_fs * 0.0078 * 0.53      # ≈ text height when rotated
    head_run = len(PATHWAY_SHORT[cat_boundaries[-1][0]]) * head_fs * 0.0078 * 0.85
    left_margin = 0.35 + max(len(s) for s in samples) * 0.068
    right_margin = max(0.4, head_run - (block_spans[-1][2] - block_spans[-1][1]) * cell)
    top_margin = 0.75 + head_rise + 0.25               # title + labels + band
    bot_margin = 1.55                                  # gene labels + legend
    fig_w = left_margin + total_w * cell + right_margin
    fig_h = top_margin + n_rows * cell + bot_margin

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(
        left=left_margin / fig_w,
        right=1 - right_margin / fig_w,
        top=1 - top_margin / fig_h,
        bottom=bot_margin / fig_h,
    )

    # Hairline row guides.
    for i in range(n_rows):
        ax.plot([-0.5, total_w - 0.5], [i, i], color=GRID, lw=0.5, zorder=0)

    gene_cat: list[str] = []
    for cat, start, end in cat_boundaries:
        gene_cat.extend([cat] * (end - start))

    for j in range(n_genes):
        hue = PATHWAY_COLOR.get(gene_cat[j], MUTED)
        for i in range(n_rows):
            draw_status_glyph(ax, gene_x[j], n_rows - 1 - i,
                              int(target_codes[i, j]), hue)

    ax.set_xlim(-0.8, total_w - 0.2)
    ax.set_ylim(-0.8, n_rows - 0.2)
    ax.set_aspect("equal")

    # ── pathway headers: colour band + angled label in ink ────────────────
    band_y, band_h = n_rows - 0.30, 0.26
    for cat, x0, x1 in block_spans:
        ax.add_patch(mpatches.Rectangle(
            (x0 - 0.42, band_y), (x1 - x0) + 0.84, band_h,
            facecolor=PATHWAY_COLOR.get(cat, MUTED), edgecolor="none",
            clip_on=False, zorder=4))
        ax.text(x0 - 0.30, band_y + band_h + 0.22, PATHWAY_SHORT.get(cat, cat),
                rotation=head_rot, ha="left", va="bottom",
                rotation_mode="anchor", fontsize=head_fs, fontweight="bold",
                color=INK, clip_on=False, zorder=5)

    ax.set_yticks([n_rows - 1 - i for i in range(n_rows)])
    ax.set_yticklabels(samples, fontsize=8, color=INK)
    ax.set_xticks(gene_x)
    ax.set_xticklabels(ordered_targets, rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=7, color=INK2)
    ax.xaxis.set_ticks_position("bottom")
    ax.tick_params(axis="both", length=0, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.text(0.25 / fig_w, 1 - 0.30 / fig_h,
             f"{CYCLE_LETTER}-cycle gene presence  ·  {n_rows} genomes × {n_genes} genes",
             fontsize=13, fontweight="bold", color=INK, ha="left", va="center")
    fig.text(0.25 / fig_w, 1 - 0.55 / fig_h,
             "genomes ordered by similarity of gene content · "
             "colour = pathway, glyph = call status",
             fontsize=8, color=INK2, ha="left", va="center")

    # Status legend, one row, directly under the gene labels.
    draw_status_legend(ax, 0.0, -0.8 - 3.6, "#6b6a66", fontsize=8,
                       horizontal=True, char_w=0.205)

    out.parent.mkdir(parents=True, exist_ok=True)
    dpi = 300 if out.suffix.lower() == ".png" else 100
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ───────────────────────────── main ──────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, help="comma-separated names")
    ap.add_argument("--results-dir", required=True, type=Path)
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--matrix-out", required=True, type=Path)
    ap.add_argument("--heatmap-out", required=True, type=Path)
    args = ap.parse_args()

    samples = [s for s in args.samples.split(",") if s]
    cfg = yaml.safe_load(open(args.targets))
    targets = cfg["targets"]
    complexes = list(cfg.get("complexes", {}).keys())
    synergies = [s["name"] for s in cfg.get("synergies", [])]

    # Target ordering: group by category in CAT_ORDER.
    targets_by_cat: dict[str, list[str]] = {c: [] for c in CAT_ORDER}
    for t in targets:
        targets_by_cat.setdefault(t["category"], []).append(t["id"])
    ordered_targets: list[str] = []
    cat_boundaries: list[tuple[str, int, int]] = []   # (cat, start, end)
    for cat in CAT_ORDER:
        if not targets_by_cat[cat]:
            continue
        start = len(ordered_targets)
        ordered_targets.extend(targets_by_cat[cat])
        cat_boundaries.append((cat, start, len(ordered_targets)))

    # Build numeric matrix.
    target_codes = np.zeros((len(samples), len(ordered_targets)), dtype=int)
    complex_scores = np.zeros((len(samples), len(complexes)), dtype=float)
    synergy_scores = np.zeros((len(samples), len(synergies)), dtype=float)

    rows_tsv = []
    for i, s in enumerate(samples):
        calls_p = args.results_dir / s / "calls" / CALLS_TSV
        comp_p = args.results_dir / s / "calls" / "complex_completeness.tsv"
        syn_p = args.results_dir / s / "calls" / "synergy_completeness.tsv"
        calls = load_calls(calls_p)
        comps = load_completeness(comp_p)
        syns = load_completeness(syn_p)
        for j, tid in enumerate(ordered_targets):
            target_codes[i, j] = STATUS_CODE.get(calls.get(tid, "absent"), 0)
        for j, cid in enumerate(complexes):
            complex_scores[i, j] = comps.get(cid, 0.0)
        for j, sn in enumerate(synergies):
            synergy_scores[i, j] = syns.get(sn, 0.0)
        # TSV row.
        row = {"sample": s}
        for tid in ordered_targets:
            row[f"target__{tid}"] = STATUS_CODE.get(calls.get(tid, "absent"), 0)
        for cid in complexes:
            row[f"complex__{cid}"] = round(comps.get(cid, 0.0), 3)
        for sn in synergies:
            row[f"synergy__{sn}"] = round(syns.get(sn, 0.0), 3)
        rows_tsv.append(row)

    # Write matrix TSV.
    args.matrix_out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows_tsv[0].keys())
    with open(args.matrix_out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows_tsv:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ── Circle / dot-grid presence-absence plot ───────────────────────────
    # Publication style: rows = genomes (left y-labels); columns = the genes
    # grouped under colored functional-category HEADER bands that span each
    # block. Each cell is a CIRCLE whose fill encodes detection status; the
    # circle hue encodes the pathway (the same hue used for the block header),
    # so each block reads as a color band. Clean gene-only grid (complex /
    # synergy panels are separate figures).
    render_circle_grid(
        samples=samples,
        ordered_targets=ordered_targets,
        cat_boundaries=cat_boundaries,
        target_codes=target_codes,
        out=args.heatmap_out,
    )
    print(f"[cross_sample_report] wrote {args.matrix_out} and {args.heatmap_out}")


if __name__ == "__main__":
    main()
