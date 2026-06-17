#!/usr/bin/env python3
"""
make_synergy_heatmap.py — render the predicted-synergy completeness panel.

Reuses the 3-tier discrete encoding from make_complex_heatmap.py but with a
different interpretation:
  complete  — all synergy members in this genome → benefit realized
  partial   — half the synergy present → co-culture candidate (find a partner
              that carries the missing members)
  absent    — synergy not realizable from this genome

Usage:
  python workflow/scripts/make_synergy_heatmap.py \
      --results-dir results \
      --targets config/targets.yaml \
      --samples Atferrooxidans_ATCC23270,Cmetallidurans_CH34,... \
      --out results/figures/synergies.svg
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


COLOR = {
    "complete":  "#1b7837",   # solid green — benefit realized
    "partial":   "#fdae61",   # amber — co-culture opportunity
    "absent":    "#f2f2f2",
}
LABEL = {
    "complete":  "complete — benefit realized",
    "partial":   "partial — co-culture candidate",
    "absent":    "absent",
}


def draw_completeness_circle(ax, x, y, completeness, status, radius=0.40):
    """Circle whose fill encodes 0–1 completeness (open → full)."""
    ax.add_patch(mpatches.Circle((x, y), radius, facecolor="white",
                                 edgecolor="#cccccc", linewidth=0.7, zorder=2))
    c = max(0.0, min(1.0, float(completeness)))
    if c <= 0.0:
        return
    inner = radius * (0.30 + 0.70 * c)
    ax.add_patch(mpatches.Circle((x, y), inner, facecolor=COLOR[status],
                                 edgecolor="none", zorder=3))


def status_from(c: float) -> str:
    if c >= 0.999:
        return "complete"
    if c >= 0.5:
        return "partial"
    return "absent"


def load_sample_synergies(path: Path) -> dict[str, tuple[float, int, int, str]]:
    out: dict[str, tuple[float, int, int, str]] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sid = r["synergy_id"]
            c = float(r["completeness"])
            n = int(r["n_present"])
            n_tot = int(r.get("n_total") or 0)
            # synergy_completeness.tsv occasionally omits n_total when 1 of 2
            # required members is the multi-target requires_missing field
            # — fall back to len(requires_present)+len(requires_missing).
            if n_tot == 0:
                pres = (r.get("requires_present") or "").split(",")
                miss = (r.get("requires_missing") or "").split(",")
                n_tot = sum(1 for x in pres + miss if x)
            out[sid] = (c, n, n_tot, r.get("status") or status_from(c))
    return out


def render(
    samples: list[str],
    sample_data: dict[str, dict[str, tuple[float, int, int, str]]],
    synergies: list[tuple[str, list[str], str]],   # (id, requires, benefit)
    out: Path,
) -> None:
    syn_ids = [s[0] for s in synergies]
    benefit_map = {s[0]: s[2] for s in synergies}

    def score(s: str) -> tuple[int, int]:
        d = sample_data.get(s, {})
        n_c = sum(1 for sid in syn_ids if d.get(sid, (0, 0, 0, "absent"))[3] == "complete")
        n_p = sum(1 for sid in syn_ids if d.get(sid, (0, 0, 0, "absent"))[3] == "partial")
        return (n_c, n_p)

    samples = sorted(samples, key=score, reverse=True)
    n_rows = len(samples)
    n_cols = len(syn_ids)

    cell_w, cell_h = 0.95, 0.34
    left_margin   = 2.6
    right_margin  = 2.4
    top_margin    = 1.1                   # title + subtitle (labels at bottom)
    bot_margin    = 1.3 + 0.5 + n_cols * 0.20   # x labels + benefit key
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
        for j, sid in enumerate(syn_ids):
            c, n, n_tot, status = d.get(sid, (0.0, 0, 0, "absent"))
            draw_completeness_circle(ax, j, n_rows - 1 - i, c, status)
            ax.text(j, n_rows - 1 - i - 0.55, f"{n}/{n_tot}",
                    ha="center", va="top", fontsize=6.5, color="#555555",
                    zorder=4)

    ax.set_xlim(-0.8, n_cols - 0.2)
    ax.set_ylim(-0.9, n_rows - 0.2)
    ax.set_aspect("equal")

    ax.set_xticks(list(range(n_cols)))
    ax.set_xticklabels(syn_ids, rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=9)
    ax.xaxis.set_ticks_position("bottom")

    ax.set_yticks([n_rows - 1 - i for i in range(n_rows)])
    ax.set_yticklabels(samples, fontsize=8)

    ax.tick_params(axis="both", length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Title.
    fig.text(left_margin / fig_w, 1 - 0.25 / fig_h,
             f"Predicted synergies — {n_cols} synergies · {n_rows} samples",
             fontsize=13, fontweight="bold", va="top")
    fig.text(left_margin / fig_w, 1 - 0.55 / fig_h,
             "cell text = n_present / n_total required members · "
             "partial = co-culture target (pair with a genome carrying the rest)",
             fontsize=8, color="#444444", va="top", style="italic")

    # Benefit key below (start under the rotated x labels).
    axes_bottom_frac = bot_margin / fig_h
    head_y = axes_bottom_frac - 1.35 / fig_h
    fig.text(left_margin / fig_w, head_y,
             "Benefit key:",
             fontsize=9, color="#333333", fontweight="bold", va="top")
    line_in = 0.18
    for k, sid in enumerate(syn_ids):
        ben = benefit_map.get(sid, "")
        short = ben.split("(")[0].strip().rstrip(".")
        if len(short) > 92:
            short = short[:89] + "…"
        y = head_y - (line_in + k * line_in) / fig_h
        fig.text(left_margin / fig_w, y,
                 f"{sid:<26s}  {short}",
                 fontsize=8, color="#444444", va="top",
                 family="monospace")

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
    ]
    ax.legend(
        handles=handles,
        bbox_to_anchor=(1.01, 1.0), loc="upper left",
        fontsize=8, frameon=False, borderpad=0.5,
        handlelength=1.5, handleheight=1.0,
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    dpi = 300 if out.suffix.lower() == ".png" else 100
    fig.savefig(out, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"[synergy_heatmap] wrote {out}  ({n_rows} samples, {n_cols} synergies)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True, type=Path)
    ap.add_argument("--targets",     required=True, type=Path)
    ap.add_argument("--samples",     required=True)
    ap.add_argument("--out",         required=True, type=Path)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.targets))
    synergies = [
        (s["name"], s.get("requires", []), s.get("benefit", ""))
        for s in cfg.get("synergies", [])
    ]

    samples = [s for s in args.samples.split(",") if s]
    sample_data: dict[str, dict] = {}
    for s in samples:
        p = args.results_dir / s / "calls" / "synergy_completeness.tsv"
        sample_data[s] = load_sample_synergies(p)

    render(samples, sample_data, synergies, args.out)


if __name__ == "__main__":
    main()
