#!/usr/bin/env python3
"""
test_regression.py — accuracy regression gate for scycle-pipeline.

Fails (exit 1) if panel accuracy drops below locked-in floors. Run AFTER the pipeline
has scored the SP1 43-genome panel into `results/multisample_matrix.tsv`
(`make regression` runs pipeline → score → gate).

Scores against the curated-FUNCTION GT by default (the SP4 authoritative reference); set env
GT_FILE=ground_truth.tsv to gate on the fully-automated KEGG-KO GT instead.

CI-lower-bound floors (harmonized with ncycle-pipeline, 2026-06-12). Switched from
point-estimate floors to genome-cluster bootstrap 95% CI lower bounds so the gate
fires on a real regression (the CI lower bound dropping below the established
baseline) rather than on point-estimate noise within the already-quantified
uncertainty. Mirrors ncycle-pipeline/validation/test_regression.py — both gates now
print the same check-row set: scored-cells, ALL-F1 CI, ALL-precision CI, ALL-FP,
trap-precision CI, hold-out-F1 CI, per-pathway-F1 CI. (Frame B trap-independence is
reported by validation/trap_independence.py; folding it into the gate is the next
harmonization step, B4.)

Baseline pinned at the 43-genome curated-FUNCTION panel (2026-06-12), with explicit
slack below each observed CI lower bound to admit normal panel evolution while
catching a real regression:

  metric                  observed CI        floor   slack   notes
  ALL F1 (whole-panel)    [0.91, 0.95]       0.85    0.06
  ALL precision           [0.93, 0.97]       0.88    0.05
  trap precision          [0.90, 0.98]       0.82    0.08    gated/trap targets stay clean
  hold-out F1             [0.91, 0.95]       0.80    0.11    20-genome generalization split
  per-pathway F1          min [0.60, 1.00]   0.45    0.15    organic_sulfur_dmsp (n=15) is the
                                                             binding small-n pathway; others ≥0.78
  ALL false-positives     28                 ≤32     —       orthogonal point check (FP-rate guard
                                                             is ALL-precision; this scales w/ panel)
  scored cells            2623               ≥1200   —       guards an empty/partial run

Usage:
  python validation/test_regression.py     # prints PASS/FAIL table, exits 0/1
  pytest validation/test_regression.py      # one test per check
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_scycle import compute_metrics  # noqa: E402

# ── CI lower-bound floors (harmonized with ncycle, 2026-06-12) ───────────────
MIN_ALL_F1_CI_LO              = 0.85
MIN_ALL_PRECISION_CI_LO       = 0.88
MIN_TRAP_PRECISION_CI_LO      = 0.82
MIN_INDEP_TRAP_PRECISION_CI_LO = 0.80   # independent-only (non-circular); observed [0.90,0.98]
MIN_HOLDOUT_F1_CI_LO          = 0.80
MIN_PATHWAY_F1_CI_LO          = 0.45

MAX_ALL_FP                = 32     # orthogonal to CI: an FP regression that leaves
                                   # F1 CI intact is still a regression. Precision is
                                   # the rate guard; this absolute count scales w/ panel.
MIN_SCORED_CELLS          = 1200   # 2623 today — guards against an empty/partial run.


def _ci_lo(metric: dict, key: str) -> float:
    """Return CI lower bound for a metric, NaN if missing."""
    ci = metric.get(f"{key}_ci")
    return ci[0] if ci is not None else float("nan")


def _ci_str(metric: dict, key: str) -> str:
    """Return human-readable point [lo, hi] string."""
    p = metric[key]
    ci = metric.get(f"{key}_ci")
    if ci is None or math.isnan(ci[0]) or math.isnan(ci[1]):
        return f"{p:.3f}"
    return f"{p:.3f} [{ci[0]:.2f}, {ci[1]:.2f}]"


def _checks(M: dict) -> list[tuple[str, str, bool]]:
    """Return [(check_name, observed_value_str, passed)]."""
    a = M["aggregates"]["ALL"]        # whole-panel (Frame A)
    trap = M["aggregates"]["trap"]
    itrap = M["indep_trap"]           # independent-only trap, seed-sourced removed (Frame B)
    hold = M["split"]["holdout"]      # 20-genome generalization split (Frame C)
    n_cells = a["TP"] + a["FP"] + a["FN"] + a["TN"]
    checks = [
        (f"scored cells >= {MIN_SCORED_CELLS}",
         str(n_cells), n_cells >= MIN_SCORED_CELLS),
        (f"ALL micro-F1 CI lo >= {MIN_ALL_F1_CI_LO}",
         _ci_str(a, "f1"), _ci_lo(a, "f1") >= MIN_ALL_F1_CI_LO),
        (f"ALL precision CI lo >= {MIN_ALL_PRECISION_CI_LO}",
         _ci_str(a, "precision"), _ci_lo(a, "precision") >= MIN_ALL_PRECISION_CI_LO),
        (f"ALL false-positives <= {MAX_ALL_FP}",
         str(a["FP"]), a["FP"] <= MAX_ALL_FP),
        (f"homology-trap precision CI lo >= {MIN_TRAP_PRECISION_CI_LO}",
         _ci_str(trap, "precision"), _ci_lo(trap, "precision") >= MIN_TRAP_PRECISION_CI_LO),
        (f"trap-independence precision CI lo >= {MIN_INDEP_TRAP_PRECISION_CI_LO}",
         _ci_str(itrap, "precision"), _ci_lo(itrap, "precision") >= MIN_INDEP_TRAP_PRECISION_CI_LO),
        (f"hold-out F1 CI lo >= {MIN_HOLDOUT_F1_CI_LO}",
         _ci_str(hold, "f1"), _ci_lo(hold, "f1") >= MIN_HOLDOUT_F1_CI_LO),
    ]
    for cat, pp in sorted(M["per_pathway"].items()):
        checks.append((f"pathway {cat} F1 CI lo >= {MIN_PATHWAY_F1_CI_LO}",
                       _ci_str(pp, "f1"), _ci_lo(pp, "f1") >= MIN_PATHWAY_F1_CI_LO))
    return checks


# ── pytest entry points (one assertion per check) ────────────────────────────

def test_regression():
    failures = [name for name, _, ok in _checks(compute_metrics()) if not ok]
    assert not failures, "accuracy regression: " + "; ".join(failures)


# ── CLI entry point (human-readable table + exit code) ───────────────────────

def main() -> int:
    results = _checks(compute_metrics())
    width = max(len(name) for name, _, _ in results)
    print("scycle-pipeline accuracy regression gate (panel = results/multisample_matrix.tsv)\n")
    n_fail = 0
    for name, value, ok in results:
        flag = "PASS" if ok else "FAIL"
        if not ok:
            n_fail += 1
        print(f"  [{flag}] {name:<{width}}  observed={value}")
    print()
    if n_fail:
        print(f"REGRESSION: {n_fail} check(s) failed — accuracy dropped below floor.")
        return 1
    print(f"OK: all {len(results)} checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
