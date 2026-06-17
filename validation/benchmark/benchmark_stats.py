#!/usr/bin/env python3
"""
benchmark_stats.py — statistically-robust head-to-head of scycle-pipeline vs
external S-cycle tools, on the SP1 panel + curated-FUNCTION ground truth.

Methods (see prereg.md):
  - unit of independence = GENOME (cells within a genome are correlated);
  - 95% CIs via genome-level cluster bootstrap (resample genomes, B=10,000);
  - paired significance: bootstrap of Δ = metric(scycle) − metric(comparator) on the
    same resampled genomes (win iff the Δ CI excludes 0; two-sided bootstrap p);
  - secondary: exact McNemar on cell-level discordant pairs (ignores clustering);
  - Benjamini–Hochberg FDR across the pre-specified family {comparator}×{trap-precision, ALL-F1};
  - two resolutions: subunit (61 targets; trap claim) and step (12 transformations).

Ground truth defaults to the curated-FUNCTION GT (set env GT_FILE=ground_truth.tsv for the
KEGG-KO GT). Runs with whatever tools are available — scycle + kofam need no extra input; pass
--metabolic / --dram <normalized.tsv> to include those.

Usage:
  GT_FILE=curated_function_gt.tsv python validation/benchmark/benchmark_stats.py \
         [--bootstrap 10000] [--seed 1234] [--metabolic PATH] [--dram PATH]
"""
from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

import adapters as A  # same dir

REF = "scycle"                       # paired comparisons are ours-vs-each-comparator
PRIMARY = [("trap", "precision"), ("ALL", "f1")]   # pre-registered BH family (subunit res.) — FROZEN
# SECONDARY (post-hoc, 2026-06-16): TN-inclusive metrics. MCC uses all four confusion cells and is
# robust to the absent:present class imbalance; reported with CI + a SEPARATE BH correction so the
# pre-registered confirmatory family above is unaffected. Specificity is reported (CI) but not BH'd.
SECONDARY = [("ALL", "mcc"), ("trap", "mcc")]
METRICS = ("precision", "recall", "specificity", "f1", "mcc")


# ── metrics ──────────────────────────────────────────────────────────────────

def prf(tp, fp, fn):
    p = tp / (tp + fp) if (tp + fp) else float("nan")
    r = tp / (tp + fn) if (tp + fn) else float("nan")
    f = (2 * p * r / (p + r)) if (not math.isnan(p) and not math.isnan(r) and (p + r)) else float("nan")
    return p, r, f


def metric_of(tally, name):
    tp, fp, fn, tn = tally
    p, r, f = prf(tp, fp, fn)
    if name in ("precision", "recall", "f1"):
        return {"precision": p, "recall": r, "f1": f}[name]
    if name == "specificity":                       # TN / (TN + FP) — true-negative rate
        return tn / (tn + fp) if (tn + fp) else float("nan")
    if name == "mcc":                               # Matthews corr. coef. — all four cells, in [-1,1]
        denom = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
        return ((tp * tn - fp * fn) / math.sqrt(denom)) if denom else float("nan")
    raise KeyError(name)


def fmt(v):
    return "NA" if (isinstance(v, float) and math.isnan(v)) else f"{v:.3f}"


# ── build scored cells at a given resolution ─────────────────────────────────

def cells_subunit(truth):
    """{cell: expected_bool} and {genome: [cells]} at subunit resolution."""
    exp = {(g, t): (v == "present") for (g, t), v in truth.items()}
    by_g = {}
    for (g, t) in exp:
        by_g.setdefault(g, []).append((g, t))
    return exp, by_g


def collapse_to_step(d, genomes, truth, is_truth):
    """Collapse a subunit map to step resolution via STEP_DIAG (diagnostic markers).
    A step is scored for a genome iff ≥1 diagnostic marker is a ground-truth cell;
    both truth and predictions aggregate over EXACTLY those scored markers (so an
    unscored subunit prediction can't create a phantom step-level false positive).
    Step present iff any scored diagnostic marker is present."""
    out = {}
    for g in genomes:
        for step, diags in A.STEP_DIAG.items():
            scored = [t for t in diags if (g, t) in truth]
            if not scored:
                continue
            if is_truth:
                out[(g, step)] = any(truth[(g, t)] == "present" for t in scored)
            else:
                out[(g, step)] = any(d.get((g, t), False) for t in scored)
    return out


def cells_step(truth, genomes):
    truth_step = collapse_to_step(truth, genomes, truth, is_truth=True)
    by_g = {}
    for (g, s) in truth_step:
        by_g.setdefault(g, []).append((g, s))
    return truth_step, by_g


def subset_cells(by_g, exp, which):
    """Filter the genome→cells map to a subset (subunit res.: ALL/trap/nontrap)."""
    if which == "ALL":
        keep = lambda c: True
    elif which == "trap":
        keep = lambda c: c[1] in A.TRAP
    elif which == "nontrap":
        keep = lambda c: c[1] not in A.TRAP
    else:  # a specific step name
        keep = lambda c: c[1] == which
    out = {}
    for g, cells in by_g.items():
        kc = [c for c in cells if keep(c)]
        if kc:
            out[g] = kc
    return out


# ── per-genome tallies (precomputed → fast bootstrap) ────────────────────────

def genome_tallies(pred, exp, by_g):
    """g -> (TP,FP,FN,TN) over that genome's cells."""
    gt = {}
    for g, cells in by_g.items():
        TP = FP = FN = TN = 0
        for c in cells:
            p = pred.get(c, False); e = exp[c]
            if p and e: TP += 1
            elif p and not e: FP += 1
            elif (not p) and e: FN += 1
            else: TN += 1
        gt[g] = (TP, FP, FN, TN)
    return gt


def _sum(tallies):
    return tuple(sum(x) for x in zip(*tallies)) if tallies else (0, 0, 0, 0)


def boot_ci(gt, genomes, name, B, rng):
    """Genome cluster-bootstrap 95% percentile CI for one metric."""
    vals = []
    glist = list(genomes)
    for _ in range(B):
        samp = [gt[rng.choice(glist)] for _ in glist]
        v = metric_of(_sum(samp), name)
        if not math.isnan(v):
            vals.append(v)
    if not vals:
        return (float("nan"), float("nan"))
    vals.sort()
    lo = vals[int(0.025 * len(vals))]
    hi = vals[min(len(vals) - 1, int(0.975 * len(vals)))]
    return (lo, hi)


def paired_delta(gt_ref, gt_cmp, genomes, name, B, rng):
    """Bootstrap Δ = metric(ref) − metric(cmp) on the SAME resampled genomes.
    Returns (point_delta, lo, hi, two_sided_p)."""
    point = metric_of(_sum([gt_ref[g] for g in genomes]), name) - \
            metric_of(_sum([gt_cmp[g] for g in genomes]), name)
    deltas = []
    glist = list(genomes)
    for _ in range(B):
        samp = [rng.choice(glist) for _ in glist]
        a = metric_of(_sum([gt_ref[g] for g in samp]), name)
        b = metric_of(_sum([gt_cmp[g] for g in samp]), name)
        if not (math.isnan(a) or math.isnan(b)):
            deltas.append(a - b)
    if not deltas:
        return (point, float("nan"), float("nan"), float("nan"))
    deltas.sort()
    lo = deltas[int(0.025 * len(deltas))]
    hi = deltas[min(len(deltas) - 1, int(0.975 * len(deltas)))]
    n = len(deltas)
    frac_le = sum(1 for d in deltas if d <= 0) / n
    frac_ge = sum(1 for d in deltas if d >= 0) / n
    p = min(1.0, 2 * min(frac_le, frac_ge))
    return (point, lo, hi, p)


def mcnemar_exact(pred_ref, pred_cmp, exp, cells):
    """Exact (binomial) McNemar over cell-level correctness. Returns (b, c, p)."""
    b = c = 0
    for cell in cells:
        e = exp[cell]
        ok_ref = (pred_ref.get(cell, False) == e)
        ok_cmp = (pred_cmp.get(cell, False) == e)
        if ok_ref and not ok_cmp:
            b += 1
        elif ok_cmp and not ok_ref:
            c += 1
    n = b + c
    if n == 0:
        return b, c, 1.0
    k = min(b, c)
    p = min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) * (0.5 ** n))
    return b, c, p


def bh(pvals):
    """Benjamini–Hochberg adjusted q-values (preserve input order)."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        val = min(prev, pvals[i] * m / k)
        q[i] = prev = val
    return q


# ── driver ───────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--metabolic", type=Path)
    ap.add_argument("--dram", type=Path)
    ap.add_argument("--scycdb", type=Path, help="SCycDB normalized TSV (domain-DB comparator)")
    ap.add_argument("--out", type=Path,
                    default=A.ROOT / "validation" / "benchmark" / "benchmark_results.tsv")
    args = ap.parse_args()
    rng = random.Random(args.seed)

    truth, genomes = A.load_truth()
    paths = {"metabolic": args.metabolic, "dram": args.dram, "scycdb": args.scycdb}

    preds = {}
    for name, (loader, needs_path) in A.LOADERS.items():
        p = loader(paths[name], genomes) if needs_path else loader(genomes)
        if p:
            preds[name] = p
    tools = list(preds)
    comparators = [t for t in tools if t != REF]
    print(f"tools available: {tools}   (B={args.bootstrap}, n_genomes={len(genomes)})")
    if REF not in preds:
        raise SystemExit("scycle predictions missing — run the pipeline on the panel first.")

    exp_su, byg_su = cells_subunit(truth)
    exp_st, byg_st = cells_step(truth, genomes)

    out_rows = [("resolution", "tool", "subset", "metric", "value",
                 "ci_lo", "ci_hi", "delta_vs_scycle", "d_lo", "d_hi", "boot_p", "mcnemar_p")]
    family_p, family_key = [], []
    secondary_p, secondary_key = [], []

    def run_resolution(res, exp, byg, subsets, step_preds):
        print(f"\n{'='*70}\n{res.upper()} resolution\n{'='*70}")
        for which in subsets:
            sc = subset_cells(byg, exp, which)
            if not sc:
                continue
            cells = [c for cs in sc.values() for c in cs]
            gpos = sum(1 for c in cells if exp[c]);  # positives in subset
            print(f"\n[{which}]  cells={len(cells)} (present={gpos})")
            gt = {t: genome_tallies(step_preds[t], exp, sc) for t in tools}
            for t in tools:
                tot = _sum(list(gt[t].values()))
                for mname in METRICS:
                    val = metric_of(tot, mname)
                    lo, hi = boot_ci(gt[t], list(sc), mname, args.bootstrap, rng)
                    drow = ["", "", "", "", ""]
                    if t != REF:
                        d, dlo, dhi, bp = paired_delta(gt[REF], gt[t], list(sc),
                                                       mname, args.bootstrap, rng)
                        _, _, mp = mcnemar_exact(step_preds[REF], step_preds[t], exp, cells)
                        drow = [fmt(d), fmt(dlo), fmt(dhi), f"{bp:.4f}", f"{mp:.4f}"]
                        if res == "subunit" and (which, mname) in PRIMARY:
                            family_p.append(bp)
                            family_key.append((t, which, mname, d, dlo, dhi, bp))
                        if res == "subunit" and (which, mname) in SECONDARY:
                            secondary_p.append(bp)
                            secondary_key.append((t, which, mname, d, dlo, dhi, bp))
                    print(f"  {t:<10}{mname:<10} {fmt(val)}  [{fmt(lo)}, {fmt(hi)}]"
                          + (f"   Δ={drow[0]} [{drow[1]},{drow[2]}] bootP={drow[3]} mcN={drow[4]}"
                             if t != REF else ""))
                    out_rows.append((res, t, which, mname, fmt(val), fmt(lo), fmt(hi), *drow))

    run_resolution("subunit", exp_su, byg_su, ["ALL", "trap", "nontrap"], preds)
    # collapse each tool's subunit calls to step resolution before scoring at step level
    preds_step = {t: collapse_to_step(preds[t], genomes, truth, is_truth=False) for t in tools}
    step_names = sorted(A.STEP_DIAG)
    run_resolution("step", exp_st, byg_st, ["ALL", *step_names], preds_step)

    # ── BH correction over the pre-specified primary family ──
    if family_p:
        q = bh(family_p)
        print(f"\n{'='*70}\nPRIMARY FAMILY — BH-FDR corrected (subunit resolution)\n{'='*70}")
        print(f"  {'comparator':<11}{'metric':<22}{'Δ(ours−cmp)':>12}{'95% CI':>18}{'bootP':>9}{'BH q':>8}  verdict")
        for (t, which, mname, d, dlo, dhi, bp), qi in zip(family_key, q):
            sig = (qi < 0.05 and dlo > 0)
            verdict = "scycle WINS" if sig else ("n.s." if not (dlo > 0 or dhi < 0) else "see CI")
            print(f"  {t:<11}{which+' '+mname:<22}{fmt(d):>12}{('['+fmt(dlo)+','+fmt(dhi)+']'):>18}"
                  f"{bp:>9.4f}{qi:>8.4f}  {verdict}")

    # ── secondary (post-hoc) BH over the TN-inclusive MCC contrasts ──
    if secondary_p:
        q2 = bh(secondary_p)
        print(f"\n{'='*70}\nSECONDARY FAMILY — MCC, BH-FDR corrected (post-hoc; subunit res.)\n{'='*70}")
        print(f"  {'comparator':<11}{'metric':<22}{'Δ(ours−cmp)':>12}{'95% CI':>18}{'bootP':>9}{'BH q':>8}  verdict")
        for (t, which, mname, d, dlo, dhi, bp), qi in zip(secondary_key, q2):
            sig = (qi < 0.05 and dlo > 0)
            verdict = f"{REF} WINS" if sig else ("n.s." if not (dlo > 0 or dhi < 0) else "see CI")
            print(f"  {t:<11}{which+' '+mname:<22}{fmt(d):>12}{('['+fmt(dlo)+','+fmt(dhi)+']'):>18}"
                  f"{bp:>9.4f}{qi:>8.4f}  {verdict}")

    with open(args.out, "w") as fh:
        for r in out_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
