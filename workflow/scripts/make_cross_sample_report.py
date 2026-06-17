#!/usr/bin/env python3
"""
make_cross_sample_report.py — assemble cross-sample matrix.tsv and heatmap.svg.

Matrix columns:
  per-target status codes:   0=absent, 1=domain-only, 2=confirmed, -1=disqualified
  per-complex completeness:  float 0–1
  per-synergy completeness:  float 0–1

Heatmap: matplotlib figure, rows = samples, columns blocked into 8 pathway
categories with white gutters separating them. Cells colored by per-target
status (categorical) for the target columns, and by per-complex completeness
on a viridis scale for the complex/synergy columns. Right-side annotation
strip colors each row by its dominant application area.
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


STATUS_CODE = {
    "confirmed":     2,
    "domain-only":   1,
    "narrow-no-IPR": 1,
    "disqualified": -1,
    "absent":        0,
}

CAT_ORDER = [
    "dissimilatory_sulfate_reduction", "sulfur_oxidation",
    "assimilatory_sulfate_reduction", "thiosulfate_polysulfide",
    "organic_sulfur_dmsp", "sulfonate_taurine",
]

# Map each obligatory complex to the S-cycle process / lifestyle it confers.
APPLICATION_FOR_COMPLEX = {
    "aps_reductase":            "dissimilatory_s_reduction",
    "dsr_reductase":            "dissimilatory_s_reduction",
    "qmo_complex":              "dissimilatory_s_reduction",
    "sox_core":                 "sulfur_oxidation",
    "sox_dehydrogenase":        "sulfur_oxidation",
    "soe_sulfite_dh":           "sulfur_oxidation",
    "assim_atp_sulfurylase":    "s_assimilation",
    "assim_sulfite_reductase":  "s_assimilation",
    "tetrathionate_reductase":  "thiosulfate_polysulfide",
    "thiosulfate_reductase":    "thiosulfate_polysulfide",
    "taurine_transporter":      "organosulfonate",
}


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


def dominant_application(comp_scores: dict[str, float]) -> str:
    best, best_v = "", -1.0
    for cx, v in comp_scores.items():
        app = APPLICATION_FOR_COMPLEX.get(cx)
        if app and v > best_v:
            best, best_v = app, v
    return best if best_v >= 0.5 else "none"


# ───────────────── circle / dot-grid rendering vocabulary ───────────────────
# One distinct hue per pathway (qualitative, chroma-distinct). The same hue is
# used for that pathway's column-group header AND for the circles in its block,
# so each block reads as a color band (mirrors the reference's category coloring).
PATHWAY_COLOR = {
    "dissimilatory_sulfate_reduction": "#4477AA",   # blue
    "sulfur_oxidation":                "#EE6677",   # red/coral
    "assimilatory_sulfate_reduction":  "#228833",   # green
    "thiosulfate_polysulfide":         "#CCBB44",   # yellow/gold
    "organic_sulfur_dmsp":             "#AA3377",   # purple/magenta
    "sulfonate_taurine":               "#44AA99",   # teal
}
# Human-readable header labels for each pathway block.
PATHWAY_LABEL = {
    "dissimilatory_sulfate_reduction": "Dissimilatory SO4 reduction",
    "sulfur_oxidation":                "Sulfur oxidation",
    "assimilatory_sulfate_reduction":  "Assimilatory SO4 reduction",
    "thiosulfate_polysulfide":         "Thiosulfate / tetrathionate",
    "organic_sulfur_dmsp":             "Organic S (DMSP/DMS)",
    "sulfonate_taurine":               "Sulfonate / taurine",
}


def _draw_status_circle(ax, x, y, code, hue, radius=0.36):
    """Draw one presence/absence circle at (x, y).

    Status → fill style:
      2  confirmed             → solid filled in the pathway hue
      1  domain-only / narrow  → lighter / partial fill (45% alpha) of the hue
     -1  disqualified          → open ring, thin gray edge
      0  absent                → very faint empty circle (faint gray edge)
    """
    import matplotlib.patches as mp
    if code == 2:
        ax.add_patch(mp.Circle((x, y), radius, facecolor=hue,
                               edgecolor=hue, linewidth=0.6, zorder=3))
    elif code == 1:
        ax.add_patch(mp.Circle((x, y), radius, facecolor=hue,
                               edgecolor=hue, linewidth=0.6, alpha=0.45,
                               zorder=3))
    elif code == -1:
        ax.add_patch(mp.Circle((x, y), radius, facecolor="white",
                               edgecolor="#888888", linewidth=0.9, zorder=3))
    else:  # absent
        ax.add_patch(mp.Circle((x, y), radius, facecolor="none",
                               edgecolor="#dddddd", linewidth=0.6, zorder=2))


def render_circle_grid(samples, ordered_targets, cat_boundaries,
                       target_codes, out) -> None:
    """Render the cross-sample gene grid as a publication-style dot grid.

    Genes are laid out left→right grouped by pathway with thin gutters between
    blocks; genomes are rows. Co-renders to whichever extension `out` carries
    (PNG at 300 DPI, SVG as vector).
    """
    import matplotlib.lines as mlines

    n_rows = len(samples)
    n_genes = len(ordered_targets)
    n_blocks = len(cat_boundaries)
    gutter = 0.6   # gap (in cell units) between pathway blocks

    # Compute an x-position per gene, inserting a gutter between blocks.
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
    left_margin  = 2.9                   # genome names
    right_margin = 0.4
    top_margin   = 2.8                   # angled pathway labels + header bands + title
    bot_margin   = 3.0                   # rotated gene labels + 2-part legend
    fig_w = left_margin + total_w * cell + right_margin
    fig_h = top_margin + n_rows * cell + bot_margin

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(
        left=left_margin / fig_w,
        right=1 - right_margin / fig_w,
        top=1 - top_margin / fig_h,
        bottom=bot_margin / fig_h,
    )

    # Light gridlines: one faint horizontal line per row, faint vertical guide
    # per gene — kept very subtle for the airy reference look.
    for i in range(n_rows + 1):
        ax.axhline(i - 0.5, color="#f0f0f0", lw=0.5, zorder=0)

    # Map gene index → pathway hue.
    gene_cat: list[str] = []
    for cat, start, end in cat_boundaries:
        for _ in range(start, end):
            gene_cat.append(cat)

    # Draw circles. Rows drawn top-down (row 0 at top).
    for j in range(n_genes):
        hue = PATHWAY_COLOR.get(gene_cat[j], "#777777")
        gx = gene_x[j]
        for i in range(n_rows):
            gy = n_rows - 1 - i
            _draw_status_circle(ax, gx, gy, int(target_codes[i, j]), hue)

    ax.set_xlim(-0.8, total_w - 0.2)
    ax.set_ylim(-0.8, n_rows - 0.2)
    ax.set_aspect("equal")

    # ── pathway headers: thin color band + angled label above each block ───
    # The label is drawn at ~32° above its block so it reads regardless of how
    # narrow the block is (e.g. the 2-gene nitrite-oxidation block) — this fixes
    # white text being clipped inside thin colored bars.
    band_y = n_rows - 0.30          # just above the top row
    band_h = 0.30
    for cat, x0, x1 in block_spans:
        hue = PATHWAY_COLOR.get(cat, "#777777")
        ax.add_patch(mpatches.Rectangle(
            (x0 - 0.42, band_y), (x1 - x0) + 0.84, band_h,
            facecolor=hue, edgecolor="none", clip_on=False, zorder=4))
        ax.text(x0 - 0.42, band_y + band_h + 0.25,
                PATHWAY_LABEL.get(cat, cat), rotation=32,
                ha="left", va="bottom", rotation_mode="anchor",
                fontsize=8.5, fontweight="bold", color=hue,
                clip_on=False, zorder=5)

    # ── y labels (genome names) ────────────────────────────────────────────
    ax.set_yticks([n_rows - 1 - i for i in range(n_rows)])
    ax.set_yticklabels(samples, fontsize=8)

    # ── x labels (gene ids, ~45° rotated, bottom) ──────────────────────────
    ax.set_xticks(gene_x)
    ax.set_xticklabels(ordered_targets, rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=7)
    ax.xaxis.set_ticks_position("bottom")

    ax.tick_params(axis="both", length=0, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.suptitle(
        f"S-cycle gene presence/absence  ·  {n_rows} genomes × {n_genes} genes",
        fontsize=13, fontweight="bold", x=left_margin / fig_w, ha="left",
        y=1 - 0.35 / fig_h)

    # ── two-part legend, bottom-left ───────────────────────────────────────
    # (a) status → circle style
    grey = "#555555"
    status_handles = [
        mlines.Line2D([], [], marker="o", linestyle="none", markersize=9,
                      markerfacecolor=grey, markeredgecolor=grey,
                      label="confirmed"),
        mlines.Line2D([], [], marker="o", linestyle="none", markersize=9,
                      markerfacecolor=grey, markeredgecolor=grey, alpha=0.45,
                      label="domain-only / narrow"),
        mlines.Line2D([], [], marker="o", linestyle="none", markersize=9,
                      markerfacecolor="white", markeredgecolor="#888888",
                      markeredgewidth=1.0, label="disqualified"),
        mlines.Line2D([], [], marker="o", linestyle="none", markersize=9,
                      markerfacecolor="none", markeredgecolor="#cccccc",
                      label="absent"),
    ]
    leg1 = fig.legend(handles=status_handles, title="status",
                      loc="lower left",
                      bbox_to_anchor=(left_margin / fig_w, 0.02),
                      ncol=4, fontsize=8, title_fontsize=9,
                      frameon=False, handletextpad=0.4, columnspacing=1.4)
    leg1._legend_box.align = "left"
    fig.add_artist(leg1)

    # (b) pathway → color swatch
    pathway_handles = [
        mlines.Line2D([], [], marker="o", linestyle="none", markersize=9,
                      markerfacecolor=PATHWAY_COLOR[c],
                      markeredgecolor=PATHWAY_COLOR[c],
                      label=PATHWAY_LABEL[c].replace("\n", " "))
        for c, _, _ in cat_boundaries
    ]
    fig.legend(handles=pathway_handles, title="pathway",
               loc="lower left",
               bbox_to_anchor=(left_margin / fig_w, 0.10),
               ncol=5, fontsize=8, title_fontsize=9,
               frameon=False, handletextpad=0.4, columnspacing=1.4)

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
    dom_apps: list[str] = []

    rows_tsv = []
    for i, s in enumerate(samples):
        calls_p = args.results_dir / s / "calls" / "scycle_calls.tsv"
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
        dom_apps.append(dominant_application(comps))
        # TSV row.
        row = {"sample": s}
        for tid in ordered_targets:
            row[f"target__{tid}"] = STATUS_CODE.get(calls.get(tid, "absent"), 0)
        for cid in complexes:
            row[f"complex__{cid}"] = round(comps.get(cid, 0.0), 3)
        for sn in synergies:
            row[f"synergy__{sn}"] = round(syns.get(sn, 0.0), 3)
        row["dominant_application"] = dom_apps[-1]
        rows_tsv.append(row)

    # Write matrix TSV.
    args.matrix_out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows_tsv[0].keys())
    with open(args.matrix_out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows_tsv:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ── Circle / dot-grid presence-absence plot ───────────────────────────
    # Publication style: rows = genomes (left y-labels); columns = the 49 genes
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
