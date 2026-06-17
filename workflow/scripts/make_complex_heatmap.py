#!/usr/bin/env python3
"""
make_complex_heatmap.py — render the per-sample obligatory-complex panel.

Differs from the pathway heatmaps: this one uses a 3-tier discrete encoding
that matches the apply_rules.py thresholds (complete ≥0.999, partial ≥0.5,
absent <0.5) and prints the n/N completeness inside each cell so the exact
value is always legible.

Usage:
  python workflow/scripts/make_complex_heatmap.py \
      --results-dir results \
      --targets config/targets.yaml \
      --samples Atferrooxidans_ATCC23270,Cmetallidurans_CH34,... \
      --out results/figures/complexes.svg
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


# Discrete status colours — match the three thresholds in
# compute_complex_completeness.py (and the 3 statuses written to
# complex_completeness.tsv).
COLOR = {
    "complete":  "#1b7837",   # solid green
    "partial":   "#fdae61",   # amber — broken obligatory complex
    "absent":    "#f2f2f2",   # near-white grey
}
LABEL = {
    "complete":  "complete (≥0.999)",
    "partial":   "partial (≥0.5, <0.999)",
    "absent":    "absent (<0.5)",
}

# Single hue used to fill the completeness circles (consistent with the
# pathway dot-grids). Completeness 0→1 is encoded as open→full fill.
FILL_HUE = "#2166ac"


def draw_completeness_circle(ax, x, y, completeness, status, radius=0.40):
    """Circle whose fill encodes 0–1 completeness (open → full).

    A faint backing ring marks the cell; an inner solid disc grows with
    completeness. Status colour tints the disc so complete/partial/absent
    stay distinguishable at a glance.
    """
    # Faint backing ring (always present, marks the slot).
    ax.add_patch(mpatches.Circle((x, y), radius, facecolor="white",
                                 edgecolor="#cccccc", linewidth=0.7, zorder=2))
    c = max(0.0, min(1.0, float(completeness)))
    if c <= 0.0:
        return
    inner = radius * (0.30 + 0.70 * c)   # min visible disc at c>0, full at c=1
    ax.add_patch(mpatches.Circle((x, y), inner, facecolor=COLOR[status],
                                 edgecolor="none", zorder=3))


def status_from_completeness(c: float) -> str:
    if c >= 0.999:
        return "complete"
    if c >= 0.5:
        return "partial"
    return "absent"


def load_sample_complexes(path: Path) -> dict[str, tuple[float, int, int, str]]:
    """Return {complex_id: (completeness, n_present, n_total, status)}."""
    out: dict[str, tuple[float, int, int, str]] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cid = r["complex_id"]
            c = float(r["completeness"])
            out[cid] = (c, int(r["n_present"]), int(r["n_total"]),
                        r.get("status") or status_from_completeness(c))
    return out


def render(
    samples: list[str],
    sample_data: dict[str, dict[str, tuple[float, int, int, str]]],
    complexes: list[tuple[str, list[str], bool, str]],   # (id, members, obligatory, application)
    out: Path,
) -> None:
    complex_ids = [c[0] for c in complexes]
    obligatory_map = {c[0]: c[2] for c in complexes}
    application_map = {c[0]: c[3] for c in complexes}

    # Sort samples: most "complete" first, then "partial", then alpha.
    def score(s: str) -> tuple[int, int]:
        d = sample_data.get(s, {})
        n_c = sum(1 for cid in complex_ids if d.get(cid, (0, 0, 0, "absent"))[3] == "complete")
        n_p = sum(1 for cid in complex_ids if d.get(cid, (0, 0, 0, "absent"))[3] == "partial")
        return (n_c, n_p)

    samples = sorted(samples, key=score, reverse=True)
    n_rows = len(samples)
    n_cols = len(complex_ids)

    # ── layout ────────────────────────────────────────────────────────────
    cell_w, cell_h = 0.85, 0.34          # wider than pathway cells (need n/N text)
    left_margin   = 2.6                  # sample names
    right_margin  = 2.2                  # legend
    top_margin    = 1.1                  # title + subtitle (labels now at bottom)
    # Reserve room for the rotated column labels + the application key below
    # (1 line per complex). +1.0 clears the bottom x labels.
    bot_margin    = 1.0 + 0.5 + n_cols * 0.20
    fig_w = left_margin + n_cols * cell_w + right_margin
    fig_h = top_margin + n_rows * cell_h + bot_margin

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(
        left=left_margin / fig_w,
        right=1 - right_margin / fig_w,
        top=1 - top_margin / fig_h,
        bottom=bot_margin / fig_h,
    )

    # Faint row gridlines.
    for i in range(n_rows):
        ax.axhline(n_rows - 1 - i, color="#f0f0f0", lw=0.5, zorder=0)

    for i, s in enumerate(samples):
        d = sample_data.get(s, {})
        for j, cid in enumerate(complex_ids):
            c, n, n_tot, status = d.get(cid, (0.0, 0, 0, "absent"))
            draw_completeness_circle(ax, j, n_rows - 1 - i, c, status)
            # n/N text below the circle so the exact value is always legible.
            ax.text(j, n_rows - 1 - i - 0.55, f"{n}/{n_tot}",
                    ha="center", va="top", fontsize=6.5, color="#555555",
                    zorder=4)

    # Mark non-obligatory complexes (currently just geobacter_eet) above the
    # column with a dashed grey band so readers know "partial here is less
    # serious than partial on an obligatory complex".
    for j, cid in enumerate(complex_ids):
        if not obligatory_map.get(cid, True):
            ax.add_patch(mpatches.Rectangle(
                (j - 0.4, n_rows - 0.15), 0.8, 0.18,
                facecolor="none", edgecolor="#888888",
                linewidth=0.8, linestyle="--", clip_on=False))

    # Axes formatting.
    ax.set_xlim(-0.8, n_cols - 0.2)
    ax.set_ylim(-0.9, n_rows - 0.2)
    ax.set_aspect("equal")

    ax.set_xticks(list(range(n_cols)))
    ax.set_xticklabels(complex_ids, rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=9)
    ax.xaxis.set_ticks_position("bottom")

    ax.set_yticks([n_rows - 1 - i for i in range(n_rows)])
    ax.set_yticklabels(samples, fontsize=8)

    ax.tick_params(axis="both", length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Title + applications strip below the column names.
    title = f"Obligatory complexes — {n_cols} complexes · {n_rows} samples"
    fig.text(left_margin / fig_w, 1 - 0.25 / fig_h, title,
             fontsize=13, fontweight="bold", va="top")
    fig.text(left_margin / fig_w, 1 - 0.55 / fig_h,
             "cell text = n_present / n_total subunits · "
             "dashed band marks non-obligatory complexes",
             fontsize=8, color="#444444", va="top", style="italic")

    # Application key below the heatmap — one line per complex.
    # Start below the rotated x labels that occupy the top of the bottom margin.
    axes_bottom_frac = bot_margin / fig_h
    head_y = axes_bottom_frac - 1.05 / fig_h
    fig.text(left_margin / fig_w, head_y,
             "Application key:",
             fontsize=9, color="#333333", fontweight="bold", va="top")
    line_in = 0.18
    for k, cid in enumerate(complex_ids):
        app = application_map.get(cid, "")
        short = app.split("(")[0].strip().rstrip(".")
        # cap length to avoid wrapping issues
        if len(short) > 92:
            short = short[:89] + "…"
        y = head_y - (line_in + k * line_in) / fig_h
        oblig = "" if obligatory_map.get(cid, True) else "  [non-obligatory]"
        fig.text(left_margin / fig_w, y,
                 f"{cid:<22s}  {short}{oblig}",
                 fontsize=8, color="#444444", va="top",
                 family="monospace")

    # Legend — circle fill grows with completeness.
    import matplotlib.lines as mlines
    handles = [
        mlines.Line2D([], [], marker="o", linestyle="none", markersize=11,
                      markerfacecolor=COLOR["complete"],
                      markeredgecolor=COLOR["complete"], label=LABEL["complete"]),
        mlines.Line2D([], [], marker="o", linestyle="none", markersize=7,
                      markerfacecolor=COLOR["partial"],
                      markeredgecolor="#cccccc", label=LABEL["partial"]),
        mlines.Line2D([], [], marker="o", linestyle="none", markersize=11,
                      markerfacecolor="white", markeredgecolor="#cccccc",
                      label=LABEL["absent"]),
        mpatches.Patch(facecolor="none", edgecolor="#888", linestyle="--",
                       label="non-obligatory complex"),
    ]
    ax.legend(
        handles=handles,
        bbox_to_anchor=(1.01, 1.0), loc="upper left",
        fontsize=8, frameon=False, borderpad=0.5,
        handlelength=1.5, handleheight=1.0,
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    # Bump DPI for raster previews; SVG is vector and ignores DPI.
    dpi = 300 if out.suffix.lower() == ".png" else 100
    fig.savefig(out, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"[complex_heatmap] wrote {out}  ({n_rows} samples, {n_cols} complexes)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True, type=Path)
    ap.add_argument("--targets",     required=True, type=Path)
    ap.add_argument("--samples",     required=True, help="comma-separated names")
    ap.add_argument("--out",         required=True, type=Path)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.targets))
    complexes_cfg = cfg.get("complexes", {})
    # Preserve declaration order from YAML.
    complexes = [
        (cid, c.get("members", []), c.get("obligatory", True),
         c.get("application", ""))
        for cid, c in complexes_cfg.items()
    ]

    samples = [s for s in args.samples.split(",") if s]
    sample_data: dict[str, dict] = {}
    for s in samples:
        p = args.results_dir / s / "calls" / "complex_completeness.tsv"
        sample_data[s] = load_sample_complexes(p)

    render(samples, sample_data, complexes, args.out)


if __name__ == "__main__":
    main()
