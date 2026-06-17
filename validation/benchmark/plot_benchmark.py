#!/usr/bin/env python3
"""
plot_benchmark.py — B10 manuscript comparison figures (sister-mirrored: the identical
file runs in ncycle and scycle; it auto-detects the focal tool and column schema).

Regenerates from the final benchmark tables — run AFTER the B8 benchmark re-run:
  Fig 1  curated_panel_accuracy  — micro-F1 (ALL) + trap-precision, ±95% CI, focal vs
         comparators, significance stars (source: benchmark_results.tsv).            [B10.1]
  Fig 2  gtdb_concordance        — pairwise positive-call Jaccard heatmap (trap) +
         top per-target divergence (source: gtdb500*/CONCORDANCE.{md,tsv}).          [B10.2]
  Fig 3  direction_validation    — B9 orthogonal direction accuracy by stratum +
         confusion matrix (source: gtdb500*/DIR_ACCURACY.tsv).                       [B10.3]

PNG (300 dpi) + SVG into validation/benchmark/figures/. Each figure is skipped (with a
note) if its source table is absent, so the script degrades gracefully.

Usage:  python validation/benchmark/plot_benchmark.py
"""
from __future__ import annotations
import csv, re, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "benchmark_results.tsv"
FIGDIR = HERE / "figures"
GTDB = sorted((HERE.parents[1] / "comparators").glob("gtdb500*"))
FOCAL_COLOR, CMP_COLOR = "#1b7837", "#9e9e9e"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "svg.fonttype": "none"})


def save(fig, stem):
    FIGDIR.mkdir(exist_ok=True)
    for ext, dpi in ((".png", 300), (".svg", 100)):
        fig.savefig(FIGDIR / f"{stem}{ext}", bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"[plot] wrote figures/{stem}.png + .svg", file=sys.stderr)


def load_results():
    rows = list(csv.DictReader(open(RESULTS), delimiter="\t"))
    focal = next(c.split("delta_vs_")[1] for c in rows[0] if c.startswith("delta_vs_"))
    return rows, focal


# ── Fig 1 — curated-panel accuracy ──────────────────────────────────────────
def fig_accuracy():
    if not RESULTS.exists():
        print("[plot] no benchmark_results.tsv — skip Fig 1", file=sys.stderr); return
    rows, focal = load_results()

    def series(subset, metric):
        d = {}
        for r in rows:
            if r["resolution"] == "subunit" and r["subset"] == subset and r["metric"] == metric:
                try:
                    d[r["tool"]] = (float(r["value"]), float(r["ci_lo"]), float(r["ci_hi"]),
                                    r.get("boot_p", ""))
                except ValueError:
                    pass
        return d

    panels = [("ALL micro-F1", series("ALL", "f1")),
              ("homology-trap precision", series("trap", "precision"))]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, (title, d) in zip(axes, panels):
        if not d:
            continue
        order = [focal] + sorted((t for t in d if t != focal), key=lambda t: -d[t][0])
        vals = [d[t][0] for t in order]
        lo = [d[t][0] - d[t][1] for t in order]
        hi = [d[t][2] - d[t][0] for t in order]
        colors = [FOCAL_COLOR if t == focal else CMP_COLOR for t in order]
        x = range(len(order))
        ax.bar(x, vals, yerr=[lo, hi], color=colors, capsize=3, edgecolor="white", width=0.72)
        for i, t in enumerate(order):
            ax.text(i, vals[i] + hi[i] + 0.02, f"{vals[i]:.3f}", ha="center", va="bottom", fontsize=9)
            try:
                if t != focal and float(d[t][3]) < 0.05:
                    ax.text(i, 0.03, "*", ha="center", va="bottom", fontsize=15, color="white", weight="bold")
            except (ValueError, TypeError):
                pass
        ax.set_xticks(list(x)); ax.set_xticklabels(order, rotation=30, ha="right")
        ax.set_ylim(0, 1.12); ax.set_ylabel(title); ax.set_title(title)
        ax.axhline(d[focal][0], ls="--", lw=0.8, color=FOCAL_COLOR, alpha=0.5)
    fig.suptitle(f"Curated-panel accuracy — {focal} vs comparators "
                 "(bars = mean, whiskers = 95% bootstrap CI; * = beats focal, p<0.05)", y=1.02)
    save(fig, "curated_panel_accuracy")


# ── Fig 2 — GTDB-scale concordance ──────────────────────────────────────────
def parse_jaccard(md_path):
    """Parse the 'Pairwise agreement' table → {(a,b): jaccard_trap}, tool order."""
    jac, tools = {}, []
    for line in md_path.read_text().splitlines():
        m = re.match(r"\|\s*(\w[\w-]*)\s+vs\s+(\w[\w-]*)\s*\|.*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|$", line)
        if m:
            a, b, _jall, jtrap = m.group(1), m.group(2), m.group(3), float(m.group(4))
            jac[(a, b)] = jtrap
            for t in (a, b):
                if t not in tools:
                    tools.append(t)
    return jac, tools


def fig_concordance():
    if not GTDB:
        print("[plot] no gtdb500* dir — skip Fig 2", file=sys.stderr); return
    cdir = GTDB[0]
    md, tsv = cdir / "CONCORDANCE.md", cdir / "CONCORDANCE.tsv"
    if not md.exists() or not tsv.exists():
        print("[plot] no CONCORDANCE.{md,tsv} — skip Fig 2", file=sys.stderr); return
    jac, tools = parse_jaccard(md)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8),
                             gridspec_kw={"width_ratios": [1, 1.2], "wspace": 0.55})

    # Panel A — pairwise Jaccard (trap) heatmap
    ax = axes[0]
    n = len(tools)
    import numpy as np
    M = np.full((n, n), np.nan)
    for i, a in enumerate(tools):
        M[i, i] = 1.0
        for j, b in enumerate(tools):
            v = jac.get((a, b), jac.get((b, a)))
            if v is not None:
                M[i, j] = v
    cmap = LinearSegmentedColormap.from_list("jac", ["#b2182b", "#f7f7f7", "#2166ac"])
    im = ax.imshow(M, cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(tools, rotation=40, ha="right")
    ax.set_yticks(range(n)); ax.set_yticklabels(tools)
    for i in range(n):
        for j in range(n):
            if not np.isnan(M[i, j]):
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if M[i, j] < 0.25 or M[i, j] > 0.8 else "black")
    ax.set_title("Pairwise positive-call agreement\n(Jaccard, homology-trap loci)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Jaccard")

    # Panel B — top per-target divergence
    ax = axes[1]
    trows = list(csv.DictReader(open(tsv), delimiter="\t"))
    trows = [r for r in trows if r.get("divergence", "").strip() not in ("", "-")]
    trows.sort(key=lambda r: -int(r["divergence"]))
    top = trows[:15][::-1]
    ax.barh([r["target"] for r in top], [int(r["divergence"]) for r in top], color="#762a83")
    ax.set_xlabel("present-call spread across tools (max − min)")
    ax.set_title("Top per-target divergence at GTDB scale")
    fig.suptitle(f"GTDB-scale cross-tool concordance ({cdir.name})", y=1.03)
    save(fig, "gtdb_concordance")


# ── Fig 3 — B9 orthogonal direction validation ──────────────────────────────
def fig_direction():
    if not GTDB:
        return
    # Prefer the rigorous placement(+synteny) verdict (B9.4) over the best-hit tier (B9.3):
    # best-hit is inconclusive for close paralogs (e.g. nxrA/narG), placement is the headline.
    place = GTDB[0] / "DIR_PLACEMENT.tsv"
    bh = GTDB[0] / "DIR_ACCURACY.tsv"
    tsv, tier = (place, "phylogenetic placement" + (" + operon synteny" if (GTDB[0] / "DIR_PLACEMENT.md").exists() and "synteny" in (GTDB[0] / "DIR_PLACEMENT.md").read_text().lower() else "")) \
        if place.exists() else (bh, "sequence best-hit")
    if not tsv.exists():
        print("[plot] no DIR_PLACEMENT/DIR_ACCURACY.tsv — skip Fig 3", file=sys.stderr); return
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    cols = rows[0].keys()
    ours = next((c for c in ("scycle_dir", "ncycle_type", "ours") if c in cols), None)
    phylo = next((c for c in ("placement_dir", "verdict", "phylo_dir", "phylo_type") if c in cols), None)
    if not ours or not phylo:
        print("[plot] direction schema unrecognized — skip Fig 3", file=sys.stderr); return
    AMBIG = {"intermediate", "conflict", "subunit_conflict", "no_hit", "no_call",
             "unresolved", "none", "ambiguous", "-", ""}
    defn = [r for r in rows if r[ours] not in AMBIG and r[phylo] not in AMBIG]
    classes = sorted({r[ours] for r in defn} | {r[phylo] for r in defn})

    def acc(sub):
        a = sum(1 for r in sub if r[ours] == r[phylo])
        return a, len(sub)

    # strata over ALL rows so we can report unresolved coverage, not just agreement among definite
    def stratum(name, sub_all):
        d = [r for r in sub_all if r[ours] not in AMBIG and r[phylo] not in AMBIG]
        unres = len(sub_all) - len(d)
        return name, sub_all, d, unres
    strata = [stratum("all", rows),
              stratum("characterized", [r for r in rows if r.get("characterized") == "yes"]),
              stratum("candidate phyla", [r for r in rows if r.get("characterized") != "yes"])]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [1.1, 1]})

    ax = axes[0]
    labels, accs, anns = [], [], []
    for name, sub_all, d, unres in strata:
        if not d:
            labels.append(f"{name}\n({unres} unresolved)"); accs.append(0); anns.append("—"); continue
        a, n = acc(d)
        tail = f"\n+{unres} unresolved" if unres else ""
        labels.append(f"{name}\n(n={n} definite{tail})"); accs.append(a / n); anns.append(f"{a}/{n}")
    bars = ax.bar(labels, accs, color=["#4575b4", FOCAL_COLOR, "#fdae61"][:len(labels)],
                  edgecolor="white", width=0.6)
    for b, t in zip(bars, anns):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.02, t, ha="center", fontsize=9)
    ax.set_ylim(0, 1.12); ax.set_ylabel("agreement with orthogonal sequence call")
    ax.set_title("Direction-call accuracy by stratum\n(orthogonal, non-circular; unresolved = sequence can't place)")

    # confusion matrix
    ax = axes[1]
    import numpy as np
    idx = {c: i for i, c in enumerate(classes)}
    C = np.zeros((len(classes), len(classes)), int)
    for r in defn:
        C[idx[r[ours]], idx[r[phylo]]] += 1
    im = ax.imshow(C, cmap="Blues")
    ax.set_xticks(range(len(classes))); ax.set_xticklabels(classes, rotation=30, ha="right")
    ax.set_yticks(range(len(classes))); ax.set_yticklabels(classes)
    ax.set_xlabel("sequence-phylogeny call"); ax.set_ylabel("pipeline call")
    mx = C.max() or 1
    for i in range(len(classes)):
        for j in range(len(classes)):
            ax.text(j, i, C[i, j], ha="center", va="center",
                    color="white" if C[i, j] > mx * 0.6 else "black")
    ax.set_title("Confusion matrix (definite calls)")
    fig.suptitle(f"B9 orthogonal direction validation — {tier} ({GTDB[0].name})", y=1.03)
    save(fig, "direction_validation")


if __name__ == "__main__":
    fig_accuracy()
    fig_concordance()
    fig_direction()
    print(f"[plot] done → {FIGDIR}", file=sys.stderr)
