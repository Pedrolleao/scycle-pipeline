#!/usr/bin/env python3
"""
check_report_consistency.py — do the report views agree with the calls?

For every sample in a results directory, cross-checks the calls table against the
cycle map (is every present gene drawn, with its real call?) and the gene-neighbourhood
maps (is every present gene in a drawn locus?). Also lists, for information, loci that
hold genes of several pathways and genes written next to a step of another pathway.
Exit status 1 when a present gene is missing from the cycle map or shown there as absent,
or (for samples with coordinates) missing from every drawn locus.

Usage:  python validation/check_report_consistency.py <results_dir>
"""
import sys, yaml, csv
from pathlib import Path
root = Path(__file__).resolve().parents[1]; results = Path(sys.argv[1])
sys.path.insert(0, str(root / "workflow/scripts"))
import _cycle_model as m
from _viz import load_calls, STATUS_CODE, CAT_ORDER, PRESENT_CODES
from _loci import loci_for_sample
cfg = yaml.safe_load(open(root / "config/targets.yaml")); targets = cfg["targets"]
cat = {t["id"]: t["category"] for t in targets}
calls_name = [p.name for p in (results).glob("*/calls/*_calls.tsv")][0]
tot = dict(missing_in_loci=0, muted_present=0, mixed_loci=0, cat_vs_step=0, samples=0)
for sdir in sorted(p for p in results.iterdir() if (p / "calls" / calls_name).exists()):
    s = sdir.name; tot["samples"] += 1
    calls = load_calls(sdir / "calls" / calls_name)
    codes = {t: STATUS_CODE.get(r["status"], 0) for t, r in calls.items()}
    present = {t for t, c in codes.items() if c in PRESENT_CODES}
    ctx = m.genome_context(calls); ev = m.evaluate(codes, ctx)
    # (b) present genes drawn muted/absent on the cycle map (step switched off), or not drawn at all
    drawn_codes = {}
    for sid, e in ev.items():
        for line in e["lines"]:
            for tok in line:
                drawn_codes.setdefault(tok["label"], set()).add(tok["code"])
    disp_inv = {v: k for k, v in m.DISPLAY.items()}
    muted = []
    for g in sorted(present):
        lab = m.DISPLAY.get(g, g)
        dc = drawn_codes.get(lab)
        if dc is None: muted.append(f"{g}:not-on-map")
        elif codes[g] not in dc: muted.append(f"{g}:{'/'.join(map(str,dc))}")
    # (a) present genes missing from every locus (gene map)
    loci = loci_for_sample(results, s, calls, targets)
    miss = []
    if loci or (sdir / "prodigal" / f"{s}.faa").exists():
        in_loci = {c["target"] for l in loci for g in l["genes"] for c in g["calls"] if c["code"] in PRESENT_CODES}
        shown = {c["target"] for l in loci if l["n_present"] > 0 for g in l["genes"] for c in g["calls"] if c["code"] in PRESENT_CODES}
        miss = sorted(present - shown)
    # (c) loci whose header pathway differs from some of their present genes' categories
    mixed = []
    for l in loci:
        cats = {c["category"] for g in l["genes"] for c in g["calls"] if c["code"] in PRESENT_CODES}
        if len(cats) > 1: mixed.append(f"{l['locus_id']}[{l['pathway']}]={'+'.join(sorted(cats))}")
    # (d) present genes whose active step pathway differs from their category (colour mismatch)
    cvs = []
    for st in m.STEPS:
        e = ev[st["id"]]
        if e["state"] == "absent": continue
        for line in st["lines"]:
            for slot in line:
                for g in m._slot_genes(slot):
                    if g in present and cat[g] != st["pathway"]: cvs.append(f"{g}:{cat[g]}->{st['pathway']}@{st['id']}")
    tot["missing_in_loci"] += len(miss); tot["muted_present"] += len(muted); tot["mixed_loci"] += len(mixed); tot["cat_vs_step"] += len(cvs)
    if miss or muted or mixed or cvs:
        print(f"\n## {s}  mode={ctx.get('mode','-')}  present={len(present)} loci={len(loci)}")
        if miss: print("  present but in NO drawn locus:", miss)
        if muted: print("  present but cycle map shows absent/not drawn:", muted)
        if mixed: print("  loci mixing pathways:", mixed)
        if cvs: print("  gene category != step pathway on map:", sorted(set(cvs)))
print("\nTOTAL", tot)
bad = tot["missing_in_loci"] + tot["muted_present"]
print("OK: the cycle maps and gene maps agree with the calls" if not bad else f"FAIL: {bad} present gene(s) missing or misdrawn")
sys.exit(1 if bad else 0)
