#!/usr/bin/env python3
"""
make_fig_b_comparators.py — paper figure B: pipeline F1 vs NCBIfam + KofamScan.

Two side-by-side panels. Left panel: pipeline F1 vs NCBIfam F1 across the 24
NCBIfam-eligible targets. Right panel: pipeline F1 vs KofamScan F1 across the
40 KofamScan-eligible targets. Within each panel, one row per target with
green (pipeline) and a comparator-coloured bar; significance markers
(★ Holm-α=0.05 survivor, ☆ BH-q=0.05 only) annotated at the row.

Usage:
  python workflow/scripts/figures/make_fig_b_comparators.py \
      --ncbifam validation/comparator_compare_corrected.tsv \
      --kofam   validation/kofam/kofam_compare_corrected.tsv \
      --targets-yaml config/targets.yaml \
      --out results/figures/fig_b_comparator_f1
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import yaml

PIPELINE_COLOUR = "#1b7837"   # green — matches Fig A positive seeds
NCBIFAM_COLOUR  = "#9b59b6"   # purple
KOFAM_COLOUR    = "#3498db"   # blue
GREY            = "#7f8c8d"

CATEGORY_ORDER = [
    "iron_oxidation",
    "sulfur_oxidation",
    "cyanogenesis",
    "organic_acids",
    "reductive_dissolution",
    "metal_tolerance",
]
CATEGORY_DISPLAY = {
    "iron_oxidation":         "Fe-ox",
    "sulfur_oxidation":       "S-ox",
    "cyanogenesis":           "HCN",
    "organic_acids":          "Org. acids",
    "reductive_dissolution":  "Reductive",
    "metal_tolerance":        "Metal tol.",
}


def load_category_map(config_path: Path) -> dict[str, str]:
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    return {t["id"]: t.get("category", "?") for t in cfg["targets"]}


def load_comparator(path: Path, comparator_label: str) -> list[dict]:
    """Read a *_corrected.tsv and return per-row dicts with cleaned F1s + flags.

    comparator_label is "nf" or "kf" — the column-name prefix in the file.
    """
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pipe_f1 = r["pipe_f1"]
            cmp_f1 = r[f"{comparator_label}_f1"]
            rows.append({
                "target":  r["target"],
                "pipe_f1": float(pipe_f1) if pipe_f1 not in ("", "NA") else None,
                "cmp_f1":  float(cmp_f1) if cmp_f1 not in ("", "NA") else None,
                "winner":              r["winner"],
                "winner_holm_all_05":  r.get("winner_holm_all_05", ""),
                "winner_bh_all_05":    r.get("winner_bh_all_05",   ""),
                "exp_pos": int(r["exp_pos"]),
            })
    return rows


def order_rows(rows: list[dict], cat_map: dict[str, str]) -> list[dict]:
    """Sort: category order, then pipeline F1 descending (NAs at bottom)."""
    def key(r):
        cat = cat_map.get(r["target"], "?")
        cat_idx = CATEGORY_ORDER.index(cat) if cat in CATEGORY_ORDER else len(CATEGORY_ORDER)
        f1 = r["pipe_f1"] if r["pipe_f1"] is not None else -1
        return (cat_idx, -f1, r["target"])
    return sorted(rows, key=key)


def draw_panel(ax, rows, cat_map, comparator_label, comparator_colour, title):
    """Draw one panel of the figure (one comparator). Rows already ordered."""
    n = len(rows)
    y = list(range(n))
    bar_h = 0.36

    # Build category-boundary indices for divider lines + category labels
    cat_breaks = []
    prev_cat = None
    for i, r in enumerate(rows):
        c = cat_map.get(r["target"], "?")
        if c != prev_cat:
            cat_breaks.append((i, c))
            prev_cat = c

    # Draw category background bands for visual grouping
    for k, (start, cat) in enumerate(cat_breaks):
        end = cat_breaks[k + 1][0] if k + 1 < len(cat_breaks) else n
        if k % 2 == 0:
            ax.axhspan(start - 0.5, end - 0.5, color="#000000", alpha=0.04, zorder=0)
        # Category label at the left edge of each band
        mid = (start + end - 1) / 2
        ax.text(-0.012, mid, CATEGORY_DISPLAY.get(cat, cat),
                transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=8,
                color=GREY, style="italic")

    # Bars
    for i, r in enumerate(rows):
        if r["pipe_f1"] is not None:
            ax.barh(i - bar_h / 2, r["pipe_f1"], height=bar_h,
                    color=PIPELINE_COLOUR, edgecolor="white", linewidth=0.3)
        else:
            ax.text(0.005, i - bar_h / 2, "NA", fontsize=6.5,
                    va="center", color=GREY, style="italic")
        if r["cmp_f1"] is not None:
            ax.barh(i + bar_h / 2, r["cmp_f1"], height=bar_h,
                    color=comparator_colour, edgecolor="white", linewidth=0.3)
        else:
            ax.text(0.005, i + bar_h / 2, "NA", fontsize=6.5,
                    va="center", color=GREY, style="italic")

        # Significance marker — placed at x=1.02, only when pipeline wins
        marker = ""
        if r["winner_holm_all_05"] == "pipeline":
            marker = "★"
        elif r["winner_bh_all_05"] == "pipeline":
            marker = "☆"
        if marker:
            ax.text(1.02, i, marker, fontsize=11, va="center",
                    color="#000000", weight="bold")

    ax.set_yticks(y)
    ax.set_yticklabels([r["target"] for r in rows], fontsize=8)
    ax.set_ylim(n - 0.5, -0.5)   # top to bottom
    ax.set_xlim(0, 1.08)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("F1", fontsize=9)
    ax.set_title(title, fontsize=10, weight="bold", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", labelsize=7)
    ax.tick_params(axis="y", labelsize=7.5)
    ax.grid(axis="x", color="0.85", linewidth=0.5, zorder=-1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbifam", default="validation/comparator_compare_corrected.tsv")
    ap.add_argument("--kofam",   default="validation/kofam/kofam_compare_corrected.tsv")
    ap.add_argument("--targets-yaml", default="config/targets.yaml")
    ap.add_argument("--out", default="results/figures/fig_b_comparator_f1")
    args = ap.parse_args()

    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    cat_map = load_category_map(Path(args.targets_yaml))

    nf_rows = order_rows(load_comparator(Path(args.ncbifam), "nf"), cat_map)
    kf_rows = order_rows(load_comparator(Path(args.kofam),   "kf"), cat_map)

    # Height proportional to the larger panel's row count
    h_per_row = 0.21
    fig_h = max(7.5, h_per_row * max(len(nf_rows), len(kf_rows)) + 2.5)
    fig, axes = plt.subplots(1, 2, figsize=(11, fig_h),
                              gridspec_kw={"width_ratios": [1, 1]})

    draw_panel(axes[0], nf_rows, cat_map, "nf", NCBIFAM_COLOUR,
               "vs NCBIfam (24 testable targets)")
    draw_panel(axes[1], kf_rows, cat_map, "kf", KOFAM_COLOUR,
               "vs KofamScan (40 testable targets)")

    # Legend — single shared at the top
    legend_handles = [
        mpatches.Patch(color=PIPELINE_COLOUR, label="pipeline F1"),
        mpatches.Patch(color=NCBIFAM_COLOUR,  label="NCBIfam F1"),
        mpatches.Patch(color=KOFAM_COLOUR,    label="KofamScan F1"),
        plt.Line2D([], [], color="none", marker="$\\bigstar$",
                   markersize=10, markerfacecolor="black",
                   label="Holm-α=0.05 pipeline win"),
        plt.Line2D([], [], color="none", marker="$\\star$",
                   markersize=10, markerfacecolor="black",
                   label="BH-q=0.05 only"),
    ]
    fig.legend(handles=legend_handles, loc="upper center",
               bbox_to_anchor=(0.5, 1.01), ncol=5, frameon=False, fontsize=8.5)

    fig.suptitle(
        "Fig B — pipeline vs NCBIfam and KofamScan: per-target F1 on the 34-genome panel",
        fontsize=11, y=1.035,
    )

    fig.tight_layout(rect=(0, 0, 1, 0.985))
    out_svg = out_prefix.with_suffix(".svg")
    out_png = out_prefix.with_suffix(".png")
    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    fig.savefig(out_png, format="png", dpi=200, bbox_inches="tight")
    print(f"wrote {out_svg}")
    print(f"wrote {out_png}")

    # Print survivor tallies on stderr for quick sanity check
    holm_nf = sum(1 for r in nf_rows if r["winner_holm_all_05"] == "pipeline")
    bh_nf   = sum(1 for r in nf_rows if r["winner_bh_all_05"] == "pipeline")
    holm_kf = sum(1 for r in kf_rows if r["winner_holm_all_05"] == "pipeline")
    bh_kf   = sum(1 for r in kf_rows if r["winner_bh_all_05"] == "pipeline")
    print(f"NCBIfam survivors: Holm={holm_nf}  BH={bh_nf}  (rows={len(nf_rows)})")
    print(f"KofamScan survivors: Holm={holm_kf}  BH={bh_kf}  (rows={len(kf_rows)})")


if __name__ == "__main__":
    main()
