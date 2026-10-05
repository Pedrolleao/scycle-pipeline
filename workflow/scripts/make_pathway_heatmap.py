#!/usr/bin/env python3
"""
make_pathway_heatmap.py — render one focused heatmap per pathway category.

Replaces the dense "everything-on-one-canvas" multisample_heatmap.svg with a
single image per pathway, using true categorical colours (not viridis) so
confirmed / domain-only / disqualified / absent are unambiguously distinct.

Usage:
  python workflow/scripts/make_pathway_heatmap.py \
      --matrix results/multisample_matrix.tsv \
      --targets config/targets.yaml \
      --pathway iron_oxidation \
      --out results/figures/pathway_iron_oxidation.svg
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import yaml


from _viz import (GRID, INK, INK2, MUTED, PATHWAY_COLOR, PATHWAY_LABEL,
                  draw_status_glyph, draw_status_legend)

COMPLEX_BAND = INK2      # neutral, so it never reads as a pathway hue


def load_matrix(path: Path) -> tuple[list[str], dict[str, dict[str, int]]]:
    """Return (sample_order, {sample: {target_id: status_code}})."""
    samples: list[str] = []
    data: dict[str, dict[str, int]] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        target_cols = [c for c in reader.fieldnames if c.startswith("target__")]
        for row in reader:
            s = row["sample"]
            samples.append(s)
            data[s] = {
                c[len("target__"):]: int(row[c])
                for c in target_cols
                if row[c] != ""
            }
    return samples, data


def pathway_targets(targets_cfg: list[dict], pathway: str) -> list[dict]:
    return [t for t in targets_cfg if t["category"] == pathway]


def render(
    samples: list[str],
    data: dict[str, dict[str, int]],
    targets: list[dict],
    pathway: str,
    complex_members: set[str],
    out: Path,
) -> None:
    target_ids = [t["id"] for t in targets]
    target_names = [t["name"] for t in targets]

    # Partition: rows where every target cell is 0 (absent) move to a footnote.
    def has_any_signal(s: str) -> bool:
        row = data.get(s, {})
        return any(row.get(tid, 0) != 0 for tid in target_ids)

    absent_samples = sorted([s for s in samples if not has_any_signal(s)])
    samples = [s for s in samples if has_any_signal(s)]

    # Sort kept samples: most confirmed → most domain-only → alpha.
    def score(s: str) -> tuple[int, int]:
        row = data.get(s, {})
        vals = [row.get(tid, 0) for tid in target_ids]
        return (sum(1 for v in vals if v == 2),
                sum(1 for v in vals if v == 1))

    samples = sorted(samples, key=score, reverse=True)
    n_rows = len(samples)
    n_cols = len(target_ids)

    if n_rows == 0:
        # No genome carries any gene of this pathway: say so instead of
        # drawing an empty grid (a zero-row axes cannot be laid out).
        fig = plt.figure(figsize=(7.5, 1.5))
        fig.text(0.04, 0.68, f"{PATHWAY_LABEL.get(pathway, pathway)} — "
                 f"{n_cols} genes · 0 genomes", fontsize=13,
                 fontweight="bold", color=INK, va="center")
        fig.text(0.04, 0.32, f"No gene of this pathway was found in any of the "
                 f"{len(absent_samples)} genomes.", fontsize=9, color=INK2,
                 va="center")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300 if out.suffix.lower() == ".png" else 100)
        plt.close(fig)
        print(f"[pathway_heatmap] wrote {out}  (0 active rows, "
              f"{len(absent_samples)} all-absent)")
        return

    # ── figure layout (true SVG, no rasterization) ───────────────────────────
    cell_w, cell_h = 0.42, 0.28          # inches per cell
    left_margin   = 2.6                  # for sample names
    right_margin  = 3.6                  # for legend
    top_margin    = 1.0                  # title + obligatory-complex band
    bot_margin    = 1.0                  # rotated gene labels (now at bottom)

    # Footnote: estimate how many wrapped lines we need for the all-absent list.
    footnote_lines = 0
    if absent_samples:
        # Rough wrap: ~70 chars/line at 8pt; +1 for "absent in this pathway:" heading.
        n_chars = sum(len(s) + 2 for s in absent_samples)
        per_line_chars = 70
        footnote_lines = 1 + max(1, -(-n_chars // per_line_chars))   # ceil-div
        # +0.7 extra to clear the rotated gene labels above the footnote.
        bot_margin += footnote_lines * 0.16 + 0.15 + 0.7

    fig_w = left_margin + n_cols * cell_w + right_margin
    fig_h = top_margin + n_rows * cell_h + bot_margin

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(
        left=left_margin / fig_w,
        right=1 - right_margin / fig_w,
        top=1 - top_margin / fig_h,
        bottom=bot_margin / fig_h,
    )

    hue = PATHWAY_COLOR.get(pathway, MUTED)

    # Hairline row guides.
    for i in range(n_rows):
        ax.plot([-0.5, n_cols - 0.5], [n_rows - 1 - i] * 2, color=GRID,
                lw=0.5, zorder=0)

    # Draw cells as status circles in the pathway hue.
    for i, s in enumerate(samples):
        row = data.get(s, {})
        for j, tid in enumerate(target_ids):
            v = row.get(tid, 0)
            draw_status_glyph(ax, j, n_rows - 1 - i, v, hue)

    # Mark obligatory-complex members above their column with a small tick band.
    for j, tid in enumerate(target_ids):
        if tid in complex_members:
            ax.add_patch(mpatches.Rectangle(
                (j - 0.4, n_rows - 0.15), 0.8, 0.18,
                facecolor=COMPLEX_BAND, edgecolor="none", clip_on=False))

    # Axes / labels.
    ax.set_xlim(-0.8, n_cols - 0.2)
    ax.set_ylim(-0.8, n_rows - 0.2)
    ax.set_aspect("equal")

    ax.set_xticks(list(range(n_cols)))
    ax.set_xticklabels(target_ids, rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=9)
    ax.xaxis.set_ticks_position("bottom")

    ax.set_yticks([n_rows - 1 - i for i in range(n_rows)])
    ax.set_yticklabels(samples, fontsize=8, color=INK)

    ax.tick_params(axis="both", length=0, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Title + subtitle.
    title = (f"{PATHWAY_LABEL.get(pathway, pathway)} — {n_cols} genes · "
             f"{n_rows} genomes")
    fig.text(left_margin / fig_w, 1 - 0.25 / fig_h, title,
             fontsize=13, fontweight="bold", va="top", color=INK)
    if complex_members:
        sub = "bar above a column = subunit of an obligatory complex"
        fig.text(left_margin / fig_w, 1 - 0.55 / fig_h, sub,
                 fontsize=8, color=INK2, va="top")

    # Legend: the four status glyphs, stacked to the right of the grid.
    draw_status_legend(ax, n_cols + 0.3, n_rows - 1, hue, step=1.0, fontsize=8)
    if complex_members:
        ly = n_rows - 1 - 4
        ax.add_patch(mpatches.Rectangle((n_cols + 0.3 - 0.36, ly - 0.09), 0.72,
                                        0.18, facecolor=COMPLEX_BAND,
                                        edgecolor="none", clip_on=False))
        ax.text(n_cols + 0.3 + 0.36 * 1.9, ly, "obligatory-complex subunit",
                ha="left", va="center", fontsize=8, color=INK2, clip_on=False)

    # Footnote: samples that had no signal in this pathway.
    if absent_samples:
        from textwrap import wrap
        n = len(absent_samples)
        # Char budget: ~8 chars per inch at 7pt
        usable_inches = fig_w - left_margin - 0.4
        wrapped = wrap(", ".join(absent_samples),
                       width=max(40, int(usable_inches * 9)))

        line_in = 0.16   # inches between lines
        # Gene labels occupy the bottom margin; start the footnote below them.
        head_in = 0.90   # gap between axes bottom and heading (clears x labels)
        axes_bottom_frac = bot_margin / fig_h

        head_y = axes_bottom_frac - head_in / fig_h
        fig.text(left_margin / fig_w, head_y,
                 f"No gene of this pathway found ({n}):",
                 fontsize=8, color=INK2, va="top")
        for k, line in enumerate(wrapped):
            y = head_y - (line_in + k * line_in) / fig_h
            fig.text(left_margin / fig_w, y, line,
                     fontsize=7, color="#666666", va="top")

    out.parent.mkdir(parents=True, exist_ok=True)
    # Format inferred from extension. SVG is vector (DPI irrelevant); PNG at 300 DPI for paper-grade rasters.
    dpi = 300 if out.suffix.lower() == ".png" else 100
    fig.savefig(out, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"[pathway_heatmap] wrote {out}  ({n_rows} active rows, "
          f"{len(absent_samples)} all-absent)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix",   required=True, type=Path)
    ap.add_argument("--targets",  required=True, type=Path)
    ap.add_argument("--pathway",  required=True)
    ap.add_argument("--out",      required=True, type=Path)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.targets))
    targets = pathway_targets(cfg["targets"], args.pathway)
    if not targets:
        raise SystemExit(f"no targets found in category '{args.pathway}'")

    # Find any obligatory complex whose membership lives entirely in this pathway.
    # Each complex member slot can be a string OR a list (any-of, P5.3.3) — flatten.
    target_id_set = {t["id"] for t in targets}
    complex_members: set[str] = set()
    for cx in cfg.get("complexes", {}).values():
        members: set[str] = set()
        for m in cx.get("members", []):
            if isinstance(m, str):
                members.add(m)
            else:
                members.update(m)
        if not cx.get("obligatory", True):
            continue
        if members.issubset(target_id_set) and members:
            complex_members |= members

    samples, data = load_matrix(args.matrix)
    render(samples, data, targets, args.pathway, complex_members, args.out)


if __name__ == "__main__":
    main()
