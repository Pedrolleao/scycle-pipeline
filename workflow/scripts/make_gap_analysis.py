#!/usr/bin/env python3
"""
make_gap_analysis.py — plain-text per-sample report.

Three sections:
  PRESENT (confirmed / domain-only) — what was actually found
  MISSING — CRITICAL                — obligatory complexes whose absence
                                       BLOCKS specific recovery applications
  SYNERGY GAPS                       — pairings where one half is present
                                       but the partner is absent
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from textwrap import indent

import yaml

PRESENT = {"confirmed", "domain-only"}


def load_calls(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_complexes(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fmt_present(row: dict) -> str:
    bits = [row["target_id"].ljust(14)]
    if row["pfam_hits"]:
        bits.append(row["pfam_hits"].ljust(20))
    bits.append(f"e={row['best_pfam_evalue']}" if row["best_pfam_evalue"]
                else "")
    if row["blast_acc"]:
        bits.append(f"blast {row['blast_acc']} {row['blast_pident']}%")
    bits.append(f"[{row['status']}]")
    return "  " + " ".join(b for b in bits if b)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--mode", required=True)
    ap.add_argument("--calls", required=True, type=Path)
    ap.add_argument("--complexes", required=True, type=Path)
    ap.add_argument("--synergies", required=True, type=Path)
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    targets_cfg = yaml.safe_load(open(args.targets))
    targets_by_id = {t["id"]: t for t in targets_cfg["targets"]}
    complexes_def = targets_cfg.get("complexes", {})

    calls = load_calls(args.calls)
    calls_by_id = {r["target_id"]: r for r in calls}
    complexes = load_complexes(args.complexes)
    synergies = load_complexes(args.synergies)

    lines: list[str] = []
    lines.append(f"SAMPLE: {args.sample}   mode: {args.mode}")
    lines.append("=" * 72)
    lines.append("")

    # --- PRESENT ----------------------------------------------------------
    present = [r for r in calls if r["status"] in PRESENT]
    lines.append(f"PRESENT ({len(present)} target(s)):")
    if not present:
        lines.append("  (nothing detected)")
    else:
        by_cat: dict[str, list[dict]] = {}
        for r in present:
            by_cat.setdefault(r["category"], []).append(r)
        for cat in by_cat:
            lines.append(f"  [{cat}]")
            for r in sorted(by_cat[cat], key=lambda x: x["target_id"]):
                lines.append("  " + fmt_present(r))
    lines.append("")

    # --- MISSING — CRITICAL ----------------------------------------------
    lines.append("MISSING — CRITICAL (obligatory partners absent, blocks application):")
    blocked = False
    for c in complexes:
        if c["status"] == "complete":
            continue
        if int(c["n_present"]) == 0:
            blocked = True
            lines.append(f"  {c['complex_id']}: 0/{c['n_total']} — entire complex absent")
            lines.append(f"    -> BLOCKS: {c['application']}")
            continue
        # partial
        blocked = True
        present_m = c["members_present"] or "—"
        missing_m = c["members_missing"]
        lines.append(f"  {c['complex_id']}: {c['n_present']}/{c['n_total']} "
                     f"— missing {missing_m} (present: {present_m})")
        lines.append(f"    -> BLOCKS: {c['application']}")
    if not blocked:
        lines.append("  (every obligatory complex is complete)")
    lines.append("")

    # --- SYNERGY GAPS -----------------------------------------------------
    lines.append("SYNERGY GAPS (desirable partner absent — application still possible but degraded):")
    any_synergy = False
    for s in synergies:
        if s["status"] == "complete":
            continue
        if int(s["n_present"]) == 0:
            # Skip synergies where both halves are absent — uninteresting.
            continue
        any_synergy = True
        if s.get("forbids_violated"):
            lines.append(f"  {s['synergy_id']}: ruled out by {s['forbids_violated']} "
                         "(genes present, but dsrAB runs in the other direction)")
            continue
        lines.append(f"  {s['synergy_id']}: {s['n_present']}/{s['n_total']} "
                     f"— missing {s['requires_missing']}")
        lines.append(f"    benefit lost: {s['benefit']}")
    if not any_synergy:
        lines.append("  (no partial synergies — either fully realized or both halves absent)")
    lines.append("")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    print(f"[gap_analysis] wrote {args.out}")


if __name__ == "__main__":
    main()
