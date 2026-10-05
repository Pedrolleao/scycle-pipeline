"""Read complex_completeness.tsv / synergy_completeness.tsv into the display
states shared by the static grids and report.html.

States: complete · partial · ruled_out · absent  (see _viz.draw_completeness_glyph)
"""

from __future__ import annotations

import csv
from pathlib import Path

from _domain import PRETTY_REPLACE


def _split(v: str | None) -> list[str]:
    return [x for x in (v or "").split(",") if x]


def pretty(name: str) -> str:
    """Human label for a complex / module id."""
    out = name.replace("_", " ")
    for old, new in PRETTY_REPLACE:
        out = out.replace(old, new)
    return out


def load_complexes(path: Path) -> dict[str, dict]:
    """State comes from the subunit count itself (not the >=0.5 threshold in
    the TSV), so a 1/3 complex still shows its one subunit as `partial`."""
    out: dict[str, dict] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n, n_tot = int(r["n_present"]), int(r["n_total"])
            state = ("complete" if n_tot and n == n_tot
                     else "partial" if n > 0 else "absent")
            out[r["complex_id"]] = {
                "state": state, "text": f"{n}/{n_tot}",
                "present": r.get("members_present", ""),
                "missing": r.get("members_missing", ""),
            }
    return out


def load_synergies(path: Path) -> dict[str, dict]:
    """complete / partial follow the TSV status. The "n/N" text counts REQUIRED
    slots only. `ruled_out` = required genes are there but an excluded gene is
    present too, so the phenotype does not apply."""
    out: dict[str, dict] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            violated = _split(r.get("forbids_violated"))
            n_forbids = len(_split(r.get("forbids_satisfied"))) + len(violated)
            n_slots = int(r.get("n_total") or 0) - n_forbids
            filled = n_slots - len(_split(r.get("requires_missing")))
            status = r.get("status") or "absent"
            if status in ("complete", "partial"):
                state = status
            elif violated and filled > 0:
                state = "ruled_out"
            else:
                state = "absent"
            out[r["synergy_id"]] = {
                "state": state, "text": f"{filled}/{n_slots}",
                "present": r.get("requires_present", ""),
                "missing": r.get("requires_missing", ""),
                "violated": ",".join(violated),
            }
    return out
