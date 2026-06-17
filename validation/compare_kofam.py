#!/usr/bin/env python3
"""
compare_kofam.py — quick per-target view of the scycle-pipeline (KO + custom HMM + BLAST
gating) vs RAW KofamScan-style KO assignment. The rigorous, statistics-backed head-to-head
lives in validation/benchmark/ (genome-cluster bootstrap + McNemar + BH-FDR); this script is
the lightweight per-target delta table. Defaults to the curated-FUNCTION GT (set env
GT_FILE=ground_truth.tsv for the KEGG-KO GT).

"Raw KofamScan" baseline: a target is called PRESENT in a genome if ANY of its KOs
(targets[].ko) has an hmmscan hit whose full-sequence score >= the stock KOfam
`ko_list` threshold — no BLAST gating, no custom HMMs, and shared KOs map to ALL
their targets (so K00370 calls BOTH nxrA and narG; comammox amoA below the KO
threshold is missed). This is exactly what mapping KofamScan output to pathways does.

The pipeline's calls come from results/multisample_matrix.tsv (codes 1/2 = present).
Both are scored against validation/ground_truth.tsv (curated cells only).

Output: stdout comparison (overall + homology-trap subset + per-target deltas).
"""
from __future__ import annotations
import csv, math, os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GT = ROOT / "validation" / os.environ.get("GT_FILE", "curated_function_gt.tsv")
MATRIX = ROOT / "results" / "multisample_matrix.tsv"
TARGETS = ROOT / "config" / "targets.yaml"
KO_LIST = ROOT / "resources" / ".cache" / "ko_list"
RESULTS = ROOT / "results"

# Homology-trap targets (must match validation/benchmark/adapters.py TRAP).
TRAP = {"dsrA","dsrB","dsrC","dsrD","soxD","fccA","sdo",
        "phsA","ttrA","sreA","dmdA","otr","doxD"}


def prf(tp, fp, fn):
    p = tp/(tp+fp) if (tp+fp) else float("nan")
    r = tp/(tp+fn) if (tp+fn) else float("nan")
    f = (2*p*r/(p+r)) if (not math.isnan(p) and not math.isnan(r) and (p+r)) else float("nan")
    return p, r, f

def fmt(v): return "NA" if (isinstance(v,float) and math.isnan(v)) else f"{v:.2f}"


def main():
    import yaml
    cfg = yaml.safe_load(open(TARGETS))
    targets = cfg["targets"]
    tgt_kos = {t["id"]: (t.get("ko") or []) for t in targets}

    # truth (skip hold-outs to compare on the same training panel)
    truth, hold = {}, set()
    with open(GT) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("source") == "holdout_v3":
                hold.add(r["genome"]); continue
            truth[(r["genome"], r["target"])] = r["expected"]
    genomes = sorted({g for g, _ in truth})

    # stock KOfam thresholds
    ko_thr = {}
    with open(KO_LIST) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[1] not in ("", "-"):
                try: ko_thr[p[0]] = float(p[1])
                except ValueError: pass

    # raw KO presence per genome: which KOs clear their stock threshold
    ko_present = {}  # genome -> set(KO)
    for g in genomes:
        tbl = RESULTS / g / "hmm" / f"{g}.hmmscan.tsv"
        present = set()
        if tbl.exists():
            with open(tbl) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip(): continue
                    c = line.split()
                    if len(c) < 8: continue
                    name = c[0]
                    if name in ko_thr:
                        try: score = float(c[7])
                        except ValueError: continue
                        if score >= ko_thr[name]:
                            present.add(name)
        ko_present[g] = present

    # pipeline calls from matrix (1/2 = present), training genomes only
    tool_present = {}
    with open(MATRIX) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["sample"]
            if g in hold: continue
            for col, val in r.items():
                if col.startswith("target__"):
                    try: tool_present[(g, col[len("target__"):])] = int(val) in (1, 2)
                    except (TypeError, ValueError): pass

    def score(pred_fn):
        agg = {"TP":0,"FP":0,"FN":0,"TN":0}
        per = {}
        for (g, t), exp in truth.items():
            pred = pred_fn(g, t)
            ep = (exp == "present")
            cell = "TP" if pred and ep else "FP" if pred and not ep else "FN" if not pred and ep else "TN"
            agg[cell]+=1
            per.setdefault(t,{"TP":0,"FP":0,"FN":0,"TN":0})[cell]+=1
        return agg, per

    raw_pred  = lambda g,t: any(k in ko_present.get(g,set()) for k in tgt_kos.get(t,[]))
    tool_pred = lambda g,t: tool_present.get((g,t), False)

    raw_agg, raw_per   = score(raw_pred)
    tool_agg, tool_per = score(tool_pred)

    def block(label, agg):
        p,r,f = prf(agg["TP"],agg["FP"],agg["FN"])
        print(f"  {label:<22} TP={agg['TP']:>3} FP={agg['FP']:>3} FN={agg['FN']:>3} TN={agg['TN']:>3}  "
              f"P={fmt(p)} R={fmt(r)} F1={fmt(f)}")

    def subset(per, keys):
        a={"TP":0,"FP":0,"FN":0,"TN":0}
        for t in keys:
            if t in per:
                for k in a: a[k]+=per[t][k]
        return a

    print(f"=== scycle-pipeline  vs  raw KofamScan (KO-only) — GT={GT.name} ===\n")
    print("ALL targets:")
    block("pipeline", tool_agg); block("raw KofamScan", raw_agg)
    print("\nHomology-trap targets (dsr/sox/fcc/phs/ttr/sre/dmdA/otr/doxD):")
    block("pipeline", subset(tool_per, TRAP)); block("raw KofamScan", subset(raw_per, TRAP))

    print("\nPer-trap-target F1 (pipeline | rawKO):")
    print(f"  {'target':<14}{'pipeline':>20}{'rawKO':>16}")
    for t in sorted(TRAP):
        if t not in tool_per: continue
        _,_,ft = prf(tool_per[t]["TP"],tool_per[t]["FP"],tool_per[t]["FN"])
        _,_,fr = prf(raw_per[t]["TP"],raw_per[t]["FP"],raw_per[t]["FN"])
        tp_t=tool_per[t]; tp_r=raw_per[t]
        print(f"  {t:<14}  F1={fmt(ft)} (FP{tp_t['FP']}/FN{tp_t['FN']})   F1={fmt(fr)} (FP{tp_r['FP']}/FN{tp_r['FN']})")


if __name__ == "__main__":
    main()
