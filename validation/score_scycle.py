#!/usr/bin/env python3
"""
score_scycle.py — score pipeline calls against the curated sulfur-cycle ground truth.

Adapted from ncycle-pipeline/validation/score_ncycle.py. Scores ONLY the
(genome,target) cells present in ground_truth.tsv (curated subset). Reports a
per-target precision/recall/F1 table, a per-pathway breakdown, and a homology-trap
subset (the BLAST-gated / shared-signature sulfur targets).

Matrix status codes (scycle_results/scycle_matrix.tsv): 2 confirmed, 1 domain-only or
narrow-no-IPR (→ predicted present), 0 absent, -1 disqualified (→ absent).
"""
from __future__ import annotations
import csv, math, os, random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# Default GT = curated-FUNCTION reference (curated_function_gt.tsv), the project's authoritative
# accuracy reference since SP4 (KEGG-KO GT + UniProt/operon-verified function corrections on the
# homology-trap cells; built by build_curated_function_gt.py). Set env GT_FILE=ground_truth.tsv to
# score against the fully-automated KEGG-KO GT (the conservative contrast; see REPORT.md SP4).
GT = ROOT / "validation" / os.environ.get("GT_FILE", "curated_function_gt.tsv")
MATRIX = ROOT / "scycle_results" / "scycle_matrix.tsv"
TARGETS = ROOT / "config" / "targets.yaml"
# one pair of output tables per ground truth, so that scoring the KEGG contrast does
# not overwrite the tables of the default (curated-function) ground truth
_TAG = "" if GT.name == "curated_function_gt.tsv" else "." + GT.stem
OUT_M = ROOT / "validation" / f"scycle_metrics{_TAG}.tsv"
OUT_C = ROOT / "validation" / f"scycle_confusion{_TAG}.tsv"

# BLAST-gated / shared-signature sulfur targets (requires_blast_for_confirmation or
# broad-Pfam families) — the analogue of the N-tool's homology-trap set.
TRAP = {"dsrA", "dsrB", "soxC", "soxD", "sqr", "sdo",
        "phsA", "ttrA", "sreA", "fccA", "fccB", "tsdA", "sorA", "otr"}

# ── Frame C generalization split (harmonization, 2026-06-07) ──────────────────
# scycle's headline is the WHOLE-PANEL micro-F1 (Frame A, unchanged) plus the
# per-(genome,target) independent-subset trap precision (Frame B,
# validation/trap_independence.py). To harmonize with ncycle-pipeline we ALSO
# report an explicit train/hold-out split: the HOLDOUT genomes below are the
# truly-independent panel members — NOT a curated BLAST-seed or custom-HMM
# source for ANY target (the `INDEPENDENT` set in build_ground_truth.py). Their
# calls cannot be circular, so scoring them as a held-out set is the cleanest
# fold-based generalization estimate, analogous to ncycle's holdout_v3. This is
# a descriptive reporting cut only: it does NOT change the GT file, the
# whole-panel headline, or the pre-registered benchmark (which stays whole-panel).
HOLDOUT = {"Tdenitrificans_ATCC25259", "Sdenitrificans_DSM1251", "Dshibae_DFL12",
           "Soneidensis_MR1", "Doleivorans_Hxd3", "Dbaculatum_DSM4028",
           "Dautotrophicum_HRM2", "Carsenatis_LY1", "Pmirabilis_HI4320",
           "Dpropionicus_DSM2032", "Dacetoxidans_DSM11109", "Dsulfexigens_DSM10523",
           "Tsulfidiphilus_HLEbGr7", "Tarsenitoxydans_3As", "Bjaponicum_USDA110",
           "Aaeolicus_VF5", "Atumefaciens_H13-3", "Cnecator_H16",
           "Hydrogenobaculum_Y04AAS1", "Sazoricus_FC6"}

BOOTSTRAP_B, BOOTSTRAP_SEED = 10_000, 1234   # match benchmark_stats.py / score_ncycle.py


def boot_ci(gtallies: dict, name: str, *, B: int = BOOTSTRAP_B,
            rng: random.Random | None = None) -> tuple[float, float]:
    """Genome-cluster 95% percentile bootstrap CI (ported from score_ncycle.py).
    gtallies: dict[genome -> (TP,FP,FN,TN)]. Collapses to [point,point] when the
    scope has zero FP and zero FN (bootstrap can't manufacture variability)."""
    rng = rng or random.Random(BOOTSTRAP_SEED)
    glist = list(gtallies)
    if not glist:
        return (float("nan"), float("nan"))
    idx = {"precision": (0, 1), "recall": (0, 2), "f1": None}[name]
    vals = []
    for _ in range(B):
        TP = FP = FN = 0
        for _ in glist:
            t = gtallies[rng.choice(glist)]
            TP += t[0]; FP += t[1]; FN += t[2]
        p, r, f = prf(TP, FP, FN)
        v = {"precision": p, "recall": r, "f1": f}[name]
        if not math.isnan(v):
            vals.append(v)
    if not vals:
        return (float("nan"), float("nan"))
    vals.sort()
    return (vals[int(0.025 * len(vals))], vals[min(len(vals) - 1, int(0.975 * len(vals)))])


def load_truth():
    out = {}
    with open(GT) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[(r["genome"], r["target"])] = r["expected"]
    return out


def load_calls():
    out = {}
    with open(MATRIX) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["sample"]
            for col, val in r.items():
                if col.startswith("target__"):
                    try: out[(g, col[len("target__"):])] = int(val)
                    except (TypeError, ValueError): out[(g, col[len("target__"):])] = 0
    return out


def cat_map():
    import yaml
    cfg = yaml.safe_load(open(TARGETS))
    return {t["id"]: t["category"] for t in cfg["targets"]}


def prf(tp, fp, fn):
    p = tp/(tp+fp) if (tp+fp) else float("nan")
    r = tp/(tp+fn) if (tp+fn) else float("nan")
    f = (2*p*r/(p+r)) if (not math.isnan(p) and not math.isnan(r) and (p+r)) else float("nan")
    return p, r, f


def compute_metrics() -> dict:
    truth = load_truth(); calls = load_calls(); cats = cat_map()
    rng = random.Random(BOOTSTRAP_SEED)   # single threaded RNG across all boot_ci calls
    per = {}
    confusion = []
    for (g, t), exp in truth.items():
        pred_present = calls.get((g, t), 0) in (1, 2)
        exp_pos = (exp == "present")
        cell = ("TP" if pred_present and exp_pos else
                "FP" if pred_present and not exp_pos else
                "FN" if not pred_present and exp_pos else "TN")
        confusion.append((g, t, exp, calls.get((g, t), 0),
                          "present" if pred_present else "absent", cell))
        d = per.setdefault(t, {"TP":0, "FP":0, "FN":0, "TN":0})
        d[cell] += 1

    rows = []
    for t, d in per.items():
        p, r, f = prf(d["TP"], d["FP"], d["FN"])
        rows.append({"target": t, "category": cats.get(t, "?"),
                     "n": sum(d.values()), **d,
                     "precision": p, "recall": r, "f1": f})

    _IDX = {"TP": 0, "FP": 1, "FN": 2, "TN": 3}

    def _gtallies(tid_set: set) -> dict:
        """Per-genome (TP,FP,FN,TN) tallies over the targets in tid_set — the
        genome-cluster bootstrap unit (cells within a genome are correlated)."""
        gt: dict[str, list[int]] = {}
        for (g, t, _e, _sc, _pr, cell) in confusion:
            if t not in tid_set:
                continue
            gt.setdefault(g, [0, 0, 0, 0])[_IDX[cell]] += 1
        return {g: tuple(v) for g, v in gt.items()}

    def aggregate(tids: list[str], *, with_ci: bool = True) -> dict:
        """Micro-aggregate over `tids` with genome-cluster bootstrap 95% CIs for
        precision/recall/F1 (mirrors score_ncycle.py:aggregate). CIs collapse to
        [point,point] on a zero-FP/FN scope — see boot_ci()."""
        tid_set = set(tids)
        gt = _gtallies(tid_set)
        TP = sum(v[0] for v in gt.values()); FP = sum(v[1] for v in gt.values())
        FN = sum(v[2] for v in gt.values()); TN = sum(v[3] for v in gt.values())
        p, r, f = prf(TP, FP, FN)
        f1s = sorted(m["f1"] for m in rows if m["target"] in tid_set and not math.isnan(m["f1"]))
        med = f1s[len(f1s)//2] if f1s else float("nan")
        out = {"TP":TP, "FP":FP, "FN":FN, "TN":TN, "precision":p, "recall":r,
               "f1":f, "median_f1":med}
        if with_ci:
            out["f1_ci"]        = boot_ci(gt, "f1",        rng=rng)
            out["precision_ci"] = boot_ci(gt, "precision", rng=rng)
            out["recall_ci"]    = boot_ci(gt, "recall",    rng=rng)
        return out

    trap_tids = [t for t in per if t in TRAP]
    aggregates = {"ALL": aggregate(list(per)),
                  "trap": aggregate(trap_tids),
                  "nontrap": aggregate([t for t in per if t not in TRAP])}
    per_pathway = {}
    for cat in sorted({m["category"] for m in rows}):
        tids = [m["target"] for m in rows if m["category"] == cat]
        agg = aggregate(tids, with_ci=True)
        per_pathway[cat] = {"TP": agg["TP"], "FP": agg["FP"], "FN": agg["FN"],
                            "f1": agg["f1"], "f1_ci": agg["f1_ci"]}
    # ── Frame B: independent-only trap precision (harmonized gate row, 2026-06-12) ─
    # Trap precision over trap cells EXCLUDING seed/HMM-sourced (g,t) pairs — the
    # non-circular number. is_seed/SEED_SRC come from trap_independence.py.
    from trap_independence import is_seed
    itg: dict[str, list[int]] = {}
    for (g, t, _e, _sc, _pr, cell) in confusion:
        if t not in TRAP or is_seed(g, t):
            continue
        itg.setdefault(g, [0, 0, 0, 0])[_IDX[cell]] += 1
    itg = {g: tuple(v) for g, v in itg.items()}
    iTP = sum(v[0] for v in itg.values()); iFP = sum(v[1] for v in itg.values())
    iFN = sum(v[2] for v in itg.values()); iTN = sum(v[3] for v in itg.values())
    ip, ir, iff = prf(iTP, iFP, iFN)
    indep_trap = {"TP": iTP, "FP": iFP, "FN": iFN, "TN": iTN,
                  "precision": ip, "recall": ir, "f1": iff,
                  "precision_ci": boot_ci(itg, "precision", rng=rng),
                  "f1_ci": boot_ci(itg, "f1", rng=rng)}

    # ── Frame C: train/hold-out split (genome-cluster bootstrap, descriptive) ─
    def _split_tallies(in_holdout: bool, trap_only: bool = False) -> dict:
        gt: dict[str, list[int]] = {}
        for (g, t, _e, _sc, _pr, cell) in confusion:
            if (g in HOLDOUT) != in_holdout:
                continue
            if trap_only and t not in TRAP:
                continue
            gt.setdefault(g, [0, 0, 0, 0])[_IDX[cell]] += 1
        return {g: tuple(v) for g, v in gt.items()}

    def _agg_split(gt: dict) -> dict:
        TP = sum(v[0] for v in gt.values()); FP = sum(v[1] for v in gt.values())
        FN = sum(v[2] for v in gt.values()); TN = sum(v[3] for v in gt.values())
        p, r, f = prf(TP, FP, FN)
        return {"TP": TP, "FP": FP, "FN": FN, "TN": TN, "precision": p, "recall": r,
                "f1": f, "genomes": sorted(gt), "precision_ci": boot_ci(gt, "precision", rng=rng),
                "recall_ci": boot_ci(gt, "recall", rng=rng), "f1_ci": boot_ci(gt, "f1", rng=rng)}

    split = {"train":        _agg_split(_split_tallies(False)),
             "holdout":      _agg_split(_split_tallies(True)),
             "holdout_trap": _agg_split(_split_tallies(True, trap_only=True))}

    return {"per_target": rows, "per": per, "aggregates": aggregates,
            "per_pathway": per_pathway, "indep_trap": indep_trap,
            "split": split, "confusion": confusion}


def main():
    M = compute_metrics()
    rows = M["per_target"]; confusion = M["confusion"]

    def fmt(v): return "NA" if (isinstance(v, float) and math.isnan(v)) else f"{v:.2f}"

    def fmt_ci(ci):
        if not ci: return ""
        lo, hi = ci
        return "" if (math.isnan(lo) or math.isnan(hi)) else f" [{lo:.2f}, {hi:.2f}]"

    OUT_M.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_M, "w") as fh:
        fh.write("target\tcategory\tn\tTP\tFP\tFN\tTN\tprecision\trecall\tf1\n")
        for m in sorted(rows, key=lambda x: (x["category"], x["target"])):
            fh.write(f"{m['target']}\t{m['category']}\t{m['n']}\t{m['TP']}\t{m['FP']}\t"
                     f"{m['FN']}\t{m['TN']}\t{fmt(m['precision'])}\t{fmt(m['recall'])}\t{fmt(m['f1'])}\n")
    with open(OUT_C, "w") as fh:
        fh.write("genome\ttarget\texpected\tstatus_code\tpredicted\tcell\n")
        for c in confusion:
            fh.write("\t".join(str(x) for x in c) + "\n")

    print(f"{'target':<16}{'pathway':<32}{'TP':>3}{'FP':>3}{'FN':>3}{'TN':>3}  {'P':>4}{'R':>5}{'F1':>5}")
    for m in sorted(rows, key=lambda x: (x["category"], -(x["TP"]+x["FN"]), x["target"])):
        print(f"  {m['target']:<14}{m['category']:<32}{m['TP']:>3}{m['FP']:>3}{m['FN']:>3}{m['TN']:>3}  "
              f"{fmt(m['precision'])}{fmt(m['recall']):>5}{fmt(m['f1']):>5}")

    def agg(label, a):
        print(f"\n[{label}] cells={a['TP']+a['FP']+a['FN']+a['TN']}  "
              f"TP={a['TP']} FP={a['FP']} FN={a['FN']} TN={a['TN']}  "
              f"micro-P={fmt(a['precision'])}{fmt_ci(a.get('precision_ci'))} "
              f"micro-R={fmt(a['recall'])}{fmt_ci(a.get('recall_ci'))} "
              f"micro-F1={fmt(a['f1'])}{fmt_ci(a.get('f1_ci'))}  "
              f"median-target-F1={fmt(a['median_f1'])}")

    print("\n=== Frame A — WHOLE PANEL (headline; all panel genomes) ===")
    agg("ALL curated targets", M["aggregates"]["ALL"])
    agg("Homology-trap / BLAST-gated targets", M["aggregates"]["trap"])
    agg("Non-trap targets", M["aggregates"]["nontrap"])
    print("\nper-pathway micro-F1:")
    for cat, pp in sorted(M["per_pathway"].items()):
        print(f"  {cat:<34} F1={fmt(pp['f1'])}{fmt_ci(pp.get('f1_ci'))}  "
              f"(TP={pp['TP']} FP={pp['FP']} FN={pp['FN']})")

    # ── Frame C: train/hold-out generalization split ──────────────────────────
    def split_line(label, a):
        print(f"  {label:<34} cells={a['TP']+a['FP']+a['FN']+a['TN']}  "
              f"TP={a['TP']} FP={a['FP']} FN={a['FN']} TN={a['TN']}  "
              f"P={fmt(a['precision'])}{fmt_ci(a['precision_ci'])} "
              f"R={fmt(a['recall'])}{fmt_ci(a['recall_ci'])} "
              f"F1={fmt(a['f1'])}{fmt_ci(a['f1_ci'])}")

    s = M["split"]
    print(f"\n=== Frame C — TRAIN/HOLD-OUT split (generalization) ===")
    print(f"[TRAINING: {len(s['train']['genomes'])} genomes — seed/HMM-source panel members]")
    split_line("ALL curated targets", s["train"])
    print(f"[HOLD-OUT: {len(s['holdout']['genomes'])} genomes — independent, no seed/HMM source for any target]")
    split_line("ALL curated targets", s["holdout"])
    split_line("Homology-trap targets", s["holdout_trap"])
    print(f"  hold-out genomes: {', '.join(s['holdout']['genomes'])}")
    print("\n(Frame B — per-(genome,target) independent-subset trap precision: "
          "run validation/trap_independence.py)")
    print(f"[bootstrap] 95% CIs via genome-cluster percentile bootstrap, B={BOOTSTRAP_B}, "
          f"seed={BOOTSTRAP_SEED}; collapse to [point,point] when scope has zero FP and FN.")


if __name__ == "__main__":
    main()
