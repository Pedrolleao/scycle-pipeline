#!/usr/bin/env python3
"""
make_synergy_heatmap.py — render the process-module ("synergy") panel.

One glyph per genome × module, taken from synergy_completeness.tsv:

  complete   solid disc — every required step found and no excluded gene
  partial    ring with "n/N" — n of the N REQUIRED steps found (excluded-gene
             rules are not counted into the fraction)
  ruled out  faint ring with a cross — an exclusion rule applies, so the
             module does not hold (nitrogen: an excluded gene is present, e.g.
             nosZ rules out the N2O-emitter flag; sulfur: dsrAB runs in the
             other direction)
  absent     faint ring

Modules are split into two column groups: plain process modules (only
"requires") and phenotype flags (requires + "forbids").

Usage:
  python workflow/scripts/make_synergy_heatmap.py \
      --results-dir results --targets config/targets.yaml \
      --samples A,B,... --out results/figures/synergies.svg
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from _completeness import load_synergies, pretty
from _completeness_grid import render_grid
from _viz import RULED_OUT_LABEL


def render(samples: list[str], sample_data: dict[str, dict],
           synergies: list[dict], out: Path) -> None:
    plain = [s for s in synergies if not s.get("forbids")]
    flags = [s for s in synergies if s.get("forbids")]
    groups = [("process modules — required steps only" if flags else "process modules",
               [(s["name"], pretty(s["name"])) for s in plain])]
    if flags:
        groups.append(("phenotype flags — some genes required, others excluded",
                       [(s["name"], pretty(s["name"])) for s in flags]))
    key = [(pretty(s["name"]),
            (s.get("benefit", "").split(". ")[0].strip().rstrip(".")))
           for s in plain + flags]
    render_grid(
        samples, groups, sample_data,
        title=f"Process modules — {len(synergies)} modules · {len(samples)} genomes",
        subtitle="does one genome carry every step of a process? "
                 "(from synergy_completeness.tsv)",
        legend=[("complete", "", "complete"),
                ("partial", "2/3", "partial — n of N required steps found"),
                ("ruled_out", "", RULED_OUT_LABEL),
                ("absent", "", "absent")],
        key=key, key_title="What each module means",
        out=out,
    )
    print(f"[synergy_heatmap] wrote {out}  ({len(samples)} samples, "
          f"{len(synergies)} synergies)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True, type=Path)
    ap.add_argument("--targets",     required=True, type=Path)
    ap.add_argument("--samples",     required=True)
    ap.add_argument("--out",         required=True, type=Path)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.targets))
    synergies = cfg.get("synergies", [])

    samples = [s for s in args.samples.split(",") if s]
    sample_data: dict[str, dict] = {}
    for s in samples:
        p = args.results_dir / s / "calls" / "synergy_completeness.tsv"
        sample_data[s] = {k: (v["state"], v["text"])
                          for k, v in load_synergies(p).items()}

    render(samples, sample_data, synergies, args.out)


if __name__ == "__main__":
    main()
