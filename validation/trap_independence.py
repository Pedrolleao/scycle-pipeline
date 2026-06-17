#!/usr/bin/env python3
"""
trap_independence.py — audit the homology-trap claim for seed/HMM leakage.

The benchmark's confirmatory endpoint is trap PRECISION. Some panel genomes are also
curated-BLAST-seed or custom-HMM sources for the very trap targets they test, so their
positive calls are partly circular (the model was built to detect those proteins).
Independence is a per-(genome, target) property, NOT a per-genome one: e.g. Archaeoglobus
fulgidus is a dsrA/dsrB/dsrD seed source (non-independent there) but an independent test
positive for sox/phs/ttr.

This script enumerates the present trap cells, tags each independent vs seed-sourced, and
recomputes scycle's trap precision/recall on the INDEPENDENT-only subset — the load-bearing,
non-circular number. Run from the pipeline root:
    GT_FILE=curated_function_gt.tsv python validation/trap_independence.py

SEED_SRC below is derived from the actual seed organisms in resources/blast_db/*.fasta and
targets/<t>/refs.fasta (positive seeds only; negatives don't create positive-circularity),
mapped to panel genome names. NB: phsA/ttrA seeds are currently NAME-FETCHED (uncurated, the
build_blast_db.py gene-name fallback) and happen to include the panel positive — once those
lists are curated/disabled, Styphimurium becomes independent for phsA/ttrA again.
"""
from __future__ import annotations
import sys, os
from pathlib import Path

# NB: the benchmark `adapters` import is deferred into main() so that
# `from trap_independence import is_seed, SEED_SRC` is cheap and side-effect-free
# (the regression gate imports those two for its trap-independence row).

# trap target -> panel genomes that are POSITIVE curated-seed / custom-HMM sources for it
SEED_SRC: dict[str, set[str]] = {
    "dsrA": {"Afulgidus_DSM4304", "Dvulgaris_Hildenborough"},
    "dsrB": {"Afulgidus_DSM4304", "Dvulgaris_Hildenborough"},
    "dsrD": {"Afulgidus_DSM4304"},
    "soxD": {"Pdenitrificans_PD1222", "Rdenitrificans_OCh114"},
    "doxD": {"Aambivalens_LEI10", "Aferrooxidans_ATCC23270"},
    "phsA": {"Styphimurium_LT2"},            # name-fetched fallback (uncurated) — see header
    "ttrA": {"Styphimurium_LT2", "Afulgidus_DSM4304"},  # name-fetched fallback (uncurated)
    # dsrC, fccA, sdo, sreA, dmdA, otr: seeds exclude every panel positive -> independent
}


def is_seed(g: str, t: str) -> bool:
    return g in SEED_SRC.get(t, set())


def _metrics(truth, scy, cells):
    tp = fp = fn = 0
    for (g, t) in cells:
        gt = truth[(g, t)] == "present"
        pr = scy.get((g, t), False)
        tp += pr and gt
        fp += pr and not gt
        fn += (not pr) and gt
    P = tp / (tp + fp) if tp + fp else float("nan")
    R = tp / (tp + fn) if tp + fn else float("nan")
    F = 2 * P * R / (P + R) if (P == P and R == R and P + R) else float("nan")
    return tp, fp, fn, P, R, F


def main():
    sys.path.insert(0, str(Path(__file__).resolve().parent / "benchmark"))
    import adapters
    truth, genomes = adapters.load_truth()
    scy = adapters.load_scycle(genomes)
    TRAP = adapters.TRAP

    present = sorted((g, t) for (g, t), v in truth.items() if t in TRAP and v == "present")
    print(f"=== {len(present)} PRESENT trap cells (GT), tagged ===")
    for g, t in present:
        tag = "SEED-SOURCED" if is_seed(g, t) else "independent"
        pr = "present" if scy.get((g, t), False) else "ABSENT"
        print(f"  {t:6} {g:24} pred={pr:8} [{tag}]")

    trap_cells = [(g, t) for (g, t) in truth if t in TRAP]
    indep_cells = [c for c in trap_cells if not is_seed(*c)]
    print("\n=== scycle trap metrics ===")
    for label, cells in [("ALL trap", trap_cells), ("INDEPENDENT-only trap", indep_cells)]:
        tp, fp, fn, P, R, F = _metrics(truth, scy, cells)
        npos = sum(truth[c] == "present" for c in cells)
        print(f"  {label:24} cells={len(cells)} present={npos}  "
              f"TP={tp} FP={fp} FN={fn}  P={P:.3f} R={R:.3f} F1={F:.3f}")


if __name__ == "__main__":
    main()
