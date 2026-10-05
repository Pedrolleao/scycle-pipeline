#!/usr/bin/env python3
"""
make_complex_heatmap.py — render the per-sample obligatory-complex panel.

One glyph per genome × complex: solid = every subunit found, ring with "n/N"
= some subunits found, faint ring = none. Rendering is shared with the
process-module grid (_completeness_grid.py).

Usage:
  python workflow/scripts/make_complex_heatmap.py \
      --results-dir results \
      --targets config/targets.yaml \
      --samples Atferrooxidans_ATCC23270,Cmetallidurans_CH34,... \
      --out results/figures/complexes.svg
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from _completeness import load_complexes, pretty
from _completeness_grid import render_grid


def render(
    samples: list[str],
    sample_data: dict[str, dict[str, dict]],
    complexes: list[tuple[str, list[str], bool, str]],   # (id, members, obligatory, application)
    out: Path,
) -> None:
    cells = {s: {cid: (v["state"], v["text"])
                 for cid, v in sample_data.get(s, {}).items()}
             for s in samples}

    key = []
    for cid, _, obligatory, app in complexes:
        desc = app.split(";")[0].strip().rstrip(".")
        key.append((pretty(cid), desc + ("" if obligatory else "  [non-obligatory]")))

    render_grid(
        samples,
        [("", [(cid, pretty(cid)) for cid, *_ in complexes])],
        cells,
        title=f"Obligatory complexes — {len(complexes)} complexes · {len(samples)} genomes",
        subtitle="a complex counts as complete only when every subunit is found "
                 "(confirmed or domain-only)",
        legend=[("complete", "", "complete — every subunit found"),
                ("partial", "2/3", "partial — n of N subunits found"),
                ("absent", "", "absent")],
        key=key, key_title="What each complex does",
        out=out,
    )
    print(f"[complex_heatmap] wrote {out}  ({len(samples)} samples, "
          f"{len(complexes)} complexes)")


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
        sample_data[s] = load_complexes(p)

    render(samples, sample_data, complexes, args.out)


if __name__ == "__main__":
    main()
