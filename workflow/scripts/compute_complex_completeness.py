#!/usr/bin/env python3
"""
compute_complex_completeness.py — given calls/scycle_calls.tsv and the
complexes/synergies blocks of targets.yaml, emit:

  calls/complex_completeness.tsv
      complex_id  completeness  status  members_present  members_missing  application

  calls/synergy_completeness.tsv
      synergy_id  completeness  status  requires_present  requires_missing
      forbids_satisfied  forbids_violated  n_present  n_total  benefit

The two dsr modules need the same genes, so they are told apart by the dsrAB
direction that apply_rules.py tags on the dsrA / dsrB evidence_source (e.g.
`ko|dsr_oxidative`): when a direction is called, it counts as one extra slot,
and a module that needs the OPPOSITE direction is ruled out (status absent,
the direction listed under forbids_violated). No call, or an ambiguous one,
leaves both modules scored on gene presence alone, as before.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import yaml

PRESENT_STATUSES = {"confirmed", "domain-only", "narrow-no-IPR"}

# Module → the dsrAB direction it needs. (Kept here rather than in targets.yaml:
# it mirrors resolve_dsr_direction in apply_rules.py, not a marker definition.)
DSR_DIRECTION_OF = {
    "complete_sulfate_reduction":   "reductive",
    "reverse_dsr_sulfur_oxidation": "oxidative",
}


def dsr_direction(path: Path) -> str | None:
    """'reductive' | 'oxidative' when apply_rules.py made a firm direction
    call for this genome, else None (no dsrAB, or ambiguous)."""
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["target_id"] in ("dsrA", "dsrB"):
                m = re.search(r"dsr_(reductive|oxidative)\b",
                              row.get("evidence_source", ""))
                if m:
                    return m.group(1)
    return None


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
    direction = dsr_direction(args.calls)

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
                 "requires_missing\tforbids_satisfied\tforbids_violated\t"
                 "n_present\tn_total\tbenefit\n")
        for s in synergies:
            req = s["requires"]
            present = [m for m in req if calls.get(m) in PRESENT_STATUSES]
            missing = [m for m in req if m not in present]
            # Direction rule (dsr modules only, and only when a direction is called).
            satisfied, violated = [], []
            want = DSR_DIRECTION_OF.get(s["name"])
            if want and direction:
                (satisfied if direction == want else violated).append(
                    f"dsrAB={direction}")
            n_present = len(present) + len(satisfied)
            n_total = len(req) + len(satisfied) + len(violated)
            comp = n_present / n_total if n_total else 0.0
            status = "absent" if violated else status_for(comp)
            fh.write("\t".join([
                s["name"], f"{comp:.3f}", status,
                ",".join(present), ",".join(missing),
                ",".join(satisfied), ",".join(violated),
                str(n_present), str(n_total),
                s.get("benefit", ""),
            ]) + "\n")


if __name__ == "__main__":
    main()
