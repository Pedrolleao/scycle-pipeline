#!/usr/bin/env python3
"""
make_fig_a_tc_calibration.py — paper figure A: per-Tier-2-HMM TC calibration.

For each of the 29 design-level Tier-2 custom HMMs (golB excluded per B10 Tier-3
demotion), reads targets/{tid}/tc_calibration.tsv (positive vs negative seed
bitscores) and manifest.yaml (calibrated TC bitscore), draws one small panel
showing the score histograms with the TC line marked. Panels grouped by
pathway category with a category-coloured title.

Usage:
  python workflow/scripts/figures/make_fig_a_tc_calibration.py \
      --targets-dir targets \
      --config config/targets.yaml \
      --out results/figures/fig_a_tc_calibration
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

CATEGORY_COLOUR = {
    "iron_oxidation":         "#c0392b",
    "sulfur_oxidation":       "#e67e22",
    "cyanogenesis":           "#8e44ad",
    "organic_acids":          "#16a085",
    "reductive_dissolution":  "#2c3e50",
    "metal_tolerance":        "#7f8c8d",
}

POS_COLOUR = "#1b7837"   # matches existing complex_heatmap "complete"
NEG_COLOUR = "#c0392b"
TC_COLOUR  = "#000000"


def load_tier2(config_path: Path) -> list[tuple[str, str]]:
    """Return [(target_id, category)] for current design-level Tier-2 HMMs.

    A target is Tier-2 here iff it has `custom_hmm: true`. golB was demoted to
    Tier-3 in B10 and no longer carries this flag, so it is naturally excluded.
    """
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    out = []
    for t in cfg["targets"]:
        if t.get("custom_hmm"):
            out.append((t["id"], t.get("category", "?")))
    return out


def load_calibration(cal_path: Path) -> tuple[list[float], list[float]]:
    pos, neg = [], []
    with open(cal_path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            score = float(row["bitscore"])
            if row["label"] == "positive":
                pos.append(score)
            elif row["label"] == "negative":
                neg.append(score)
    return pos, neg


def load_tc(manifest_path: Path) -> float | None:
    with open(manifest_path) as fh:
        man = yaml.safe_load(fh)
    return man.get("tc_bitscore")


# Display order: pathway category, then by target id within category.
CATEGORY_ORDER = [
    "iron_oxidation",
    "sulfur_oxidation",
    "cyanogenesis",
    "organic_acids",
    "reductive_dissolution",
    "metal_tolerance",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets-dir", default="targets",
                    help="dir containing per-target manifest.yaml + tc_calibration.tsv")
    ap.add_argument("--config", default="config/targets.yaml")
    ap.add_argument("--out", default="results/figures/fig_a_tc_calibration",
                    help="output prefix (writes .svg + .png)")
    args = ap.parse_args()

    targets_dir = Path(args.targets_dir)
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    tier2 = load_tier2(Path(args.config))
    # Sort: category in canonical order, then id alphabetically inside category
    tier2.sort(key=lambda x: (CATEGORY_ORDER.index(x[1]) if x[1] in CATEGORY_ORDER else 999, x[0]))

    n = len(tier2)
    # 6 cols × ceil(n/6) rows
    n_cols = 6
    n_rows = (n + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols * 2.4, n_rows * 1.8),
                             squeeze=False)
    fig.suptitle(
        "Fig A — TC bitscore calibration for the 29 Tier-2 custom HMMs\n"
        "positive seeds (green) vs calibration negatives (red); TC line in black",
        fontsize=11, y=1.005,
    )

    for idx, (tid, category) in enumerate(tier2):
        ax = axes[idx // n_cols][idx % n_cols]
        cal_path = targets_dir / tid / "tc_calibration.tsv"
        man_path = targets_dir / tid / "manifest.yaml"

        if not cal_path.exists() or not man_path.exists():
            ax.set_title(f"{tid} (missing)", fontsize=8)
            ax.axis("off")
            continue

        pos, neg = load_calibration(cal_path)
        tc = load_tc(man_path)

        # Determine plotting range
        all_scores = pos + neg + ([tc] if tc is not None else [])
        if not all_scores:
            ax.set_title(f"{tid} (no scores)", fontsize=8)
            ax.axis("off")
            continue
        xmin, xmax = min(all_scores), max(all_scores)
        # Pad ±5 % for visual breathing room
        span = max(xmax - xmin, 1.0)
        xmin -= span * 0.05
        xmax += span * 0.05

        # Bin: at most 20 bins across the visible range
        bins = max(8, min(20, len(pos) + len(neg) + 4))

        if pos:
            ax.hist(pos, bins=bins, range=(xmin, xmax), color=POS_COLOUR,
                    alpha=0.75, edgecolor="white", linewidth=0.5,
                    label=f"pos (n={len(pos)})")
        if neg:
            ax.hist(neg, bins=bins, range=(xmin, xmax), color=NEG_COLOUR,
                    alpha=0.65, edgecolor="white", linewidth=0.5,
                    label=f"neg (n={len(neg)})")
        if tc is not None:
            ax.axvline(tc, color=TC_COLOUR, linewidth=1.4, linestyle="-")
            # TC value annotation, just inside the top-right of the panel
            ax.text(0.97, 0.93, f"TC={tc:g}",
                    transform=ax.transAxes, fontsize=7,
                    ha="right", va="top",
                    bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                              edgecolor="0.6", linewidth=0.5))
        if not neg:
            # Many narrow-family HMMs have 0 calibration negatives because all
            # negatives e-filter before bitscore. Document that explicitly.
            ax.text(0.03, 0.93, "no neg\nabove e-cut",
                    transform=ax.transAxes, fontsize=6.5,
                    ha="left", va="top",
                    color="#7f8c8d", style="italic")

        cat_colour = CATEGORY_COLOUR.get(category, "#000000")
        ax.set_title(tid, fontsize=9, color=cat_colour, weight="bold")
        ax.tick_params(axis="both", labelsize=6.5, length=2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_yticks([])  # counts are not the message; keep frame clean

    # Blank out the trailing axes
    for idx in range(n, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].axis("off")

    # Single legend at the bottom: positive / negative / TC + category swatches
    handles = [
        mpatches.Patch(color=POS_COLOUR, alpha=0.75, label="positive seeds"),
        mpatches.Patch(color=NEG_COLOUR, alpha=0.65, label="calibration negatives"),
        plt.Line2D([0], [0], color=TC_COLOUR, linewidth=1.4, label="calibrated TC bitscore"),
    ]
    cat_handles = [mpatches.Patch(color=c, label=name.replace("_", " "))
                   for name, c in CATEGORY_COLOUR.items() if name in CATEGORY_ORDER]

    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.27, -0.025),
               ncol=3, frameon=False, fontsize=8)
    fig.legend(handles=cat_handles, loc="lower center", bbox_to_anchor=(0.73, -0.025),
               ncol=3, frameon=False, fontsize=7, title="title colour = category",
               title_fontsize=7)

    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    out_svg = out_prefix.with_suffix(".svg")
    out_png = out_prefix.with_suffix(".png")
    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    fig.savefig(out_png, format="png", dpi=200, bbox_inches="tight")
    print(f"wrote {out_svg}")
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
