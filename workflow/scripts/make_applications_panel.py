#!/usr/bin/env python3
"""
make_applications_panel.py — per-genome dominant-application breakdown.

For each genome, render a compact "card" showing:
  • genome name
  • dominant_application chip (the industrial-application label assigned by
    apply_rules.py via the complex with the highest completeness ≥ 0.5)
  • the winning complex's id + n/N completeness
  • a strip of cells, one per subunit, colored by per-target status
    (same palette as the pathway heatmaps)

When dominant_application is `none` (no complex reaches 0.5), the card shows
the closest partial complex with a grey "none" chip — so the user sees what
the genome is *near*, not just "nothing".

Cards are sorted by application (groups together cohorts that share the same
industrial use) then by completeness within each application.

Usage:
  python workflow/scripts/make_applications_panel.py \
      --matrix results/multisample_matrix.tsv \
      --results-dir results \
      --targets config/targets.yaml \
      --out results/figures/applications.svg
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


# Per-subunit status colours — same palette as the pathway heatmaps.
STATUS_COLOR = {
    2:  "#1b7837",   # confirmed
    1:  "#a6d96a",   # domain-only / narrow-no-IPR
    -1: "#fdae61",   # disqualified
    0:  "#ffffff",   # absent
}

# Complex → S-cycle process / lifestyle label (mirrors make_cross_sample_report.py).
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

# Deterministic categorical colour map for process chips.
# Order matters — drives the sort.
APPLICATION_PALETTE = {
    "dissimilatory_s_reduction": "#377eb8",
    "sulfur_oxidation":          "#e41a1c",
    "s_assimilation":            "#4daf4a",
    "thiosulfate_polysulfide":   "#a65628",
    "organosulfonate":           "#984ea3",
    "none":                      "#bcbcbc",
}
APP_ORDER = list(APPLICATION_PALETTE.keys())


def load_matrix_targets(path: Path) -> tuple[list[str], dict[str, dict[str, int]]]:
    samples: list[str] = []
    data: dict[str, dict[str, int]] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        target_cols = [c for c in reader.fieldnames if c.startswith("target__")]
        for row in reader:
            s = row["sample"]
            samples.append(s)
            data[s] = {c[len("target__"):]: int(row[c]) for c in target_cols if row[c] != ""}
    return samples, data


def load_sample_complexes(path: Path) -> dict[str, tuple[float, int, int]]:
    out: dict[str, tuple[float, int, int]] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["complex_id"]] = (
                float(r["completeness"]),
                int(r["n_present"]),
                int(r["n_total"]),
            )
    return out


def pick_dominant(comp_scores: dict[str, tuple[float, int, int]]) -> tuple[str | None, str]:
    """Return (winning_complex_id, application_label).

    Logic mirrors make_cross_sample_report.py: highest-completeness complex
    whose application is in the mapping AND whose completeness ≥ 0.5.
    Fallback: closest complex (highest completeness regardless of floor) with
    application set to 'none'.
    """
    best, best_v, best_app = None, -1.0, None
    for cx, (v, *_rest) in comp_scores.items():
        app = APPLICATION_FOR_COMPLEX.get(cx)
        if app and v > best_v:
            best, best_v, best_app = cx, v, app
    if best is not None and best_v >= 0.5:
        return best, best_app

    # Fallback for "none": closest partial regardless of completeness floor.
    closest, closest_v = None, -1.0
    for cx, (v, *_rest) in comp_scores.items():
        if cx in APPLICATION_FOR_COMPLEX and v > closest_v:
            closest, closest_v = cx, v
    return closest, "none"


def render_card(
    ax: plt.Axes,
    sample: str,
    application: str,
    complex_id: str | None,
    members: list[str],
    member_status: dict[str, int],
    n_present: int,
    n_total: int,
    completeness: float,
) -> None:
    """Draw one genome card on the given Axes."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.5)
    ax.set_aspect("equal")
    ax.set_axis_off()

    # Card background.
    ax.add_patch(mpatches.Rectangle(
        (0.05, 0.05), 9.9, 6.4,
        facecolor="#fbfbfb", edgecolor="#cccccc", linewidth=0.6,
    ))

    # Sample name.
    ax.text(0.3, 5.85, sample, fontsize=9, fontweight="bold",
            color="#222222", va="top")

    # Application chip.
    chip_color = APPLICATION_PALETTE.get(application, "#cccccc")
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.3, 4.3), 4.5, 0.9,
        boxstyle="round,pad=0.05,rounding_size=0.12",
        facecolor=chip_color, edgecolor="none",
    ))
    chip_text_color = "#222222" if application in ("none",) else "#ffffff"
    ax.text(0.3 + 4.5 / 2, 4.3 + 0.45, application,
            ha="center", va="center", fontsize=9,
            fontweight="bold", color=chip_text_color)

    # Complex name + n/N.
    if complex_id:
        ax.text(0.3, 3.7,
                f"via {complex_id}   {n_present}/{n_total}  ({completeness:.2f})",
                fontsize=8, color="#444444", va="top",
                family="monospace")
    else:
        ax.text(0.3, 3.7, "no complex matches an application",
                fontsize=8, color="#888888", va="top", style="italic")

    # Subunit strip.
    if members:
        n = len(members)
        # cells span x=0.3 to x=9.7, evenly divided
        x0, x1 = 0.3, 9.7
        cell_w = (x1 - x0) / n
        for k, m in enumerate(members):
            v = member_status.get(m, 0)
            ax.add_patch(mpatches.Rectangle(
                (x0 + k * cell_w, 1.5), cell_w - 0.05, 1.4,
                facecolor=STATUS_COLOR[v],
                edgecolor="#888888", linewidth=0.5,
            ))
            ax.text(x0 + k * cell_w + (cell_w - 0.05) / 2, 1.0,
                    m, ha="center", va="top",
                    fontsize=7, color="#333333", rotation=0)


def render(
    samples: list[str],
    target_data: dict[str, dict[str, int]],
    complex_data: dict[str, dict[str, tuple[float, int, int]]],
    complex_members: dict[str, list[str]],
    out: Path,
) -> None:
    # Compute per-genome metadata.
    rows = []
    for s in samples:
        comps = complex_data.get(s, {})
        cid, app = pick_dominant(comps)
        if cid:
            comp_v, n_p, n_tot = comps.get(cid, (0.0, 0, 0))
            members = complex_members.get(cid, [])
        else:
            comp_v, n_p, n_tot = 0.0, 0, 0
            members = []
        rows.append({
            "sample": s,
            "application": app,
            "complex_id": cid,
            "members": members,
            "n_present": n_p,
            "n_total": n_tot,
            "completeness": comp_v,
        })

    # Sort: by application (per APP_ORDER) then by completeness descending.
    app_rank = {a: i for i, a in enumerate(APP_ORDER)}
    rows.sort(key=lambda r: (app_rank.get(r["application"], 999),
                              -r["completeness"], r["sample"]))

    # Grid: aim for ~5 columns
    n = len(rows)
    n_cols = 5
    n_grid_rows = -(-n // n_cols)   # ceil

    card_w, card_h = 3.5, 2.1       # inches
    pad_x, pad_y = 0.25, 0.25
    fig_w = n_cols * (card_w + pad_x) + pad_x + 0.4   # right padding for legend
    legend_w = 3.0
    fig_w += legend_w
    fig_h = n_grid_rows * (card_h + pad_y) + pad_y + 1.4   # top space for title

    fig = plt.figure(figsize=(fig_w, fig_h))

    title_y = 1 - 0.35 / fig_h
    fig.text(0.02, title_y,
             "Dominant industrial application — per-genome breakdown",
             fontsize=14, fontweight="bold", va="top")
    fig.text(0.02, title_y - 0.42 / fig_h,
             "Each card shows the application driven by the highest-completeness complex (≥0.5); "
             "cells = subunit status; grey chip = no complex reached 0.5, closest partial shown.",
             fontsize=9, color="#444444", va="top", style="italic")

    # Card axes.
    grid_top = 1 - 1.2 / fig_h
    grid_left = (pad_x) / fig_w
    card_w_frac = card_w / fig_w
    card_h_frac = card_h / fig_h
    pad_x_frac = pad_x / fig_w
    pad_y_frac = pad_y / fig_h

    for idx, r in enumerate(rows):
        gc = idx % n_cols
        gr = idx // n_cols
        left = grid_left + gc * (card_w_frac + pad_x_frac)
        bottom = grid_top - (gr + 1) * card_h_frac - gr * pad_y_frac
        ax = fig.add_axes([left, bottom, card_w_frac, card_h_frac])
        render_card(
            ax, r["sample"], r["application"], r["complex_id"],
            r["members"], target_data.get(r["sample"], {}),
            r["n_present"], r["n_total"], r["completeness"],
        )

    # Legend (status colours + application chips) — single column on the right.
    legend_left_frac = (n_cols * (card_w + pad_x) + pad_x + 0.2) / fig_w
    legend_top = title_y - 1.0 / fig_h

    fig.text(legend_left_frac, legend_top,
             "Subunit status",
             fontsize=10, fontweight="bold", va="top")
    status_items = [
        (2,  "confirmed"),
        (1,  "domain-only / narrow-no-IPR"),
        (-1, "disqualified"),
        (0,  "absent"),
    ]
    for k, (code, label) in enumerate(status_items):
        y = legend_top - (0.32 + k * 0.22) / fig_h
        # swatch
        ax = fig.add_axes([legend_left_frac, y - 0.16 / fig_h,
                           0.18 / fig_w, 0.16 / fig_h])
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()
        ax.add_patch(mpatches.Rectangle(
            (0, 0), 1, 1, facecolor=STATUS_COLOR[code],
            edgecolor="#888", linewidth=0.5))
        fig.text(legend_left_frac + 0.24 / fig_w, y - 0.08 / fig_h,
                 label, fontsize=8, color="#333333", va="center")

    fig.text(legend_left_frac, legend_top - (0.32 + 4 * 0.22 + 0.25) / fig_h,
             "Dominant application",
             fontsize=10, fontweight="bold", va="top")
    for k, app in enumerate(APP_ORDER):
        y = legend_top - (0.32 + 4 * 0.22 + 0.55 + k * 0.22) / fig_h
        ax = fig.add_axes([legend_left_frac, y - 0.16 / fig_h,
                           0.18 / fig_w, 0.16 / fig_h])
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()
        ax.add_patch(mpatches.Rectangle(
            (0, 0), 1, 1, facecolor=APPLICATION_PALETTE[app],
            edgecolor="#888", linewidth=0.5))
        fig.text(legend_left_frac + 0.24 / fig_w, y - 0.08 / fig_h,
                 app, fontsize=8, color="#333333", va="center")

    out.parent.mkdir(parents=True, exist_ok=True)
    dpi = 300 if out.suffix.lower() == ".png" else 100
    fig.savefig(out, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"[applications_panel] wrote {out}  ({n} genomes)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix",      required=True, type=Path)
    ap.add_argument("--results-dir", required=True, type=Path)
    ap.add_argument("--targets",     required=True, type=Path)
    ap.add_argument("--out",         required=True, type=Path)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.targets))
    complex_members = {
        cid: cx.get("members", [])
        for cid, cx in cfg.get("complexes", {}).items()
    }

    samples, target_data = load_matrix_targets(args.matrix)
    complex_data = {
        s: load_sample_complexes(args.results_dir / s / "calls" / "complex_completeness.tsv")
        for s in samples
    }
    render(samples, target_data, complex_data, complex_members, args.out)


if __name__ == "__main__":
    main()
