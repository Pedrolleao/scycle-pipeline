#!/usr/bin/env python3
"""
compute_complex_completeness.py — given calls/scycle_calls.tsv and the
complexes/synergies blocks of targets.yaml, emit:

  calls/complex_completeness.tsv
      complex_id  completeness  status  members_present  members_missing  application

  calls/synergy_completeness.tsv
      synergy_id  completeness  status  requires_present  requires_missing  benefit
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import yaml

PRESENT_STATUSES = {"confirmed", "domain-only", "narrow-no-IPR"}


def load_calls(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["target_id"]] = row["status"]
    return out


def status_for(completeness: float) -> str:
    if completeness >= 0.999:
        return "complete"
    if completeness >= 0.5:
        return "partial"
    return "absent"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, type=Path)
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--complexes-out", required=True, type=Path)
    ap.add_argument("--synergies-out", required=True, type=Path)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.targets))
    complexes = cfg.get("complexes", {})
    synergies = cfg.get("synergies", [])
    calls = load_calls(args.calls)

    # Complexes.
    args.complexes_out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.complexes_out, "w") as fh:
        fh.write("complex_id\tcompleteness\tstatus\tmembers_present\t"
                 "members_missing\tn_present\tn_total\tapplication\n")
        for cid, body in complexes.items():
            members = body["members"]
            present = [m for m in members
                       if calls.get(m) in PRESENT_STATUSES]
            missing = [m for m in members if m not in present]
            comp = len(present) / len(members) if members else 0.0
            fh.write("\t".join([
                cid, f"{comp:.3f}", status_for(comp),
                ",".join(present), ",".join(missing),
                str(len(present)), str(len(members)),
                body.get("application", ""),
            ]) + "\n")

    # Synergies.
    with open(args.synergies_out, "w") as fh:
        fh.write("synergy_id\tcompleteness\tstatus\trequires_present\t"
                 "requires_missing\tn_present\tn_total\tbenefit\n")
        for s in synergies:
            req = s["requires"]
            present = [m for m in req if calls.get(m) in PRESENT_STATUSES]
            missing = [m for m in req if m not in present]
            comp = len(present) / len(req) if req else 0.0
            fh.write("\t".join([
                s["name"], f"{comp:.3f}", status_for(comp),
                ",".join(present), ",".join(missing),
                str(len(present)), str(len(req)),
                s.get("benefit", ""),
            ]) + "\n")


if __name__ == "__main__":
    main()
