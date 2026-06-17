#!/usr/bin/env python3
"""
concordance.py — cross-tool CONCORDANCE analysis (no ground truth) for the
GTDB-scale sweep. Measures where scycle / KofamScan / SCycDB AGREE and where
they systematically DIVERGE, with a spotlight on the dsrAB reductive-vs-oxidative
DIRECTION trap (the S analogue of ncycle's nxrA↔narG mis-routing): dsrA/dsrB are
the same genes/KO in dissimilatory sulfate reducers and in reverse-Dsr sulfur
oxidizers, so a tool that only reports gene presence cannot tell the two metabolic
directions apart — scycle resolves direction via dsrD + qmo/sox companions.

We deliberately do NOT score accuracy here: arbitrary GTDB genomes have no curated
truth. Instead we quantify (a) baseline pairwise agreement, (b) the trap loci where
tools split, (c) the dsrAB-direction split, and (d) scycle's distinctive calls.

Reuses the validated loaders in adapters.py by pointing them at the GTDB results
dir — identical call semantics to the curated-panel benchmark. METABOLIC joins the
panel when its normalized TSV is passed via --metabolic (the 499-genome
comparators/gtdb500_s/metabolic.tsv: 380 backbone genomes reused from ncycle's
GTDB-500 METABOLIC run — same METABOLIC v4.0, KO output is metabolism-wide so the
sulfur KOs are already there — plus 119 S-enriched genomes run in ≤12-genome batches
to stay under METABOLIC's exit-144 scale-failure threshold; see WORKPLAN). DRAM is
excluded at scale by prereg. The panel exposes the dsrAB-direction divergence across
all 499 genomes — and METABOLIC, being direction-blind (no reductive/oxidative KO),
falls into the same trap, illustrating it rather than resolving it.

Usage:
  python validation/benchmark/concordance.py \
      --results   comparators/gtdb500_s/results \
      --scycdb    comparators/gtdb500_s/scycdb.tsv \
      --panel-map comparators/gtdb500_s/selection.tsv \
      --out       comparators/gtdb500_s/CONCORDANCE.md
"""
from __future__ import annotations
import argparse, csv, itertools, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import adapters as A

TOOLS = ["scycle", "kofam", "scycdb"]


def representable() -> dict[str, set[str]]:
    """Targets each tool can in principle call (else a 'disagreement' is just a
    coverage gap, which we report separately)."""
    tkos = A.target_kos()
    has_ko = {t for t, kos in tkos.items() if kos}
    allt = set(A.all_target_ids())
    return {
        "scycle": allt,
        "kofam": has_ko,
        "scycdb": set(A.SCYC_MAP.values()),
    }


def scycle_dsr_direction(preds, repr_, g: str) -> str | None:
    """scycle's reductive-vs-oxidative call for a dsrAB-carrying genome, recomputed
    from the matrix present-calls with the SAME presence-based rule as
    apply_rules.resolve_dsr_direction (dsrD ⇒ reductive; else Sox/sqr ⇒ oxidative;
    else qmo ⇒ reductive; else ambiguous). Returns None if dsrAB is absent."""
    P = lambda t: preds["scycle"].get((g, t), False)
    if not (P("dsrA") or P("dsrB")):
        return None
    if P("dsrD"):
        return "reductive"
    if any(P(t) for t in ("soxB", "soxA", "soxY", "soxX", "sqr", "soxC")):
        return "oxidative"
    if any(P(t) for t in ("qmoA", "qmoB", "qmoC")):
        return "reductive"
    return "ambiguous"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--scycdb", type=Path, help="SCycDB normalized TSV (optional)")
    ap.add_argument("--metabolic", type=Path, help="METABOLIC normalized TSV (optional)")
    ap.add_argument("--panel-map", type=Path, help="selection.tsv (ncbi_acc, stratum, clade)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    # point the validated loaders at the GTDB results dir
    A.RESULTS = args.results
    A.MATRIX = args.results / "multisample_matrix.tsv"

    genomes = set()
    with open(A.MATRIX) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genomes.add(r["sample"])

    # scycle + kofam always available (from the scycle run); scycdb/metabolic join
    # only when their normalized TSV exists (partial runs are fine).
    preds = {
        "scycle": A.load_scycle(genomes),
        "kofam": A.load_kofam(genomes),
    }
    if args.metabolic and args.metabolic.exists():
        preds["metabolic"] = A.load_metabolic(args.metabolic, genomes)
    if args.scycdb and args.scycdb.exists():
        preds["scycdb"] = A.load_scycdb(args.scycdb, genomes)
    global TOOLS
    TOOLS = [t for t in ["scycle", "kofam", "metabolic", "scycdb"] if t in preds]
    print(f"[concordance] tools available: {TOOLS}", file=sys.stderr)
    targets = A.all_target_ids()
    repr_ = representable()
    # tools that joined but aren't statically known to representable() (none here) → all targets
    for t in TOOLS:
        repr_.setdefault(t, set(targets))
    n_g = len(genomes)

    # ── 1. pairwise agreement ────────────────────────────────────────────
    def agreement(a, b, target_set):
        """fraction of (genome,target) cells where tool a and b make the same
        present/absent call, over targets both can represent."""
        tset = [t for t in target_set if t in repr_[a] and t in repr_[b]]
        same = tot = 0
        for g in genomes:
            for t in tset:
                tot += 1
                if preds[a].get((g, t), False) == preds[b].get((g, t), False):
                    same += 1
        return (same / tot) if tot else float("nan"), tot

    def jaccard(a, b, target_set):
        """Positive-call agreement: |present in both| / |present in either|, over
        targets both can represent. Ignores the many shared-absent cells that
        inflate raw agreement — the honest concordance metric."""
        tset = [t for t in target_set if t in repr_[a] and t in repr_[b]]
        both = either = 0
        for g in genomes:
            for t in tset:
                pa = preds[a].get((g, t), False)
                pb = preds[b].get((g, t), False)
                if pa or pb:
                    either += 1
                    if pa and pb:
                        both += 1
        return (both / either) if either else float("nan")

    trap = sorted(A.TRAP)
    pairs = list(itertools.combinations(TOOLS, 2))
    agree_all = {p: agreement(*p, targets) for p in pairs}
    agree_trap = {p: agreement(*p, trap) for p in pairs}
    jac_all = {p: jaccard(*p, targets) for p in pairs}
    jac_trap = {p: jaccard(*p, trap) for p in pairs}

    # ── 2. per-target present counts per tool + divergence ───────────────
    per_target = []
    for t in targets:
        row = {"target": t}
        counts = {}
        for tool in TOOLS:
            if t in repr_[tool]:
                counts[tool] = sum(1 for g in genomes if preds[tool].get((g, t), False))
                row[tool] = counts[tool]
            else:
                row[tool] = "-"   # tool cannot represent this target
        # divergence = spread between max and min present-count among representing tools
        row["divergence"] = (max(counts.values()) - min(counts.values())) if counts else 0
        per_target.append(row)

    # ── 3. scycle distinctiveness: cells where scycle differs from ALL others
    distinct = defaultdict(int)   # target -> n genomes where scycle != every comparator
    for g in genomes:
        for t in targets:
            sc = preds["scycle"].get((g, t), False)
            others = [preds[o].get((g, t), False) for o in TOOLS[1:] if t in repr_[o]]
            if others and all(o != sc for o in others):
                distinct[t] += 1

    # ── 4. dsrAB reductive↔oxidative direction trap (the headline locus) ─
    cmp_tools = [t for t in ("kofam", "metabolic", "scycdb") if t in preds]
    direction = {"reductive": 0, "oxidative": 0, "ambiguous": 0}
    dsrab_genomes = []
    # a comparator "reports dsrAB" if it calls dsrA or dsrB present (direction-blind)
    def cmp_reports_dsrab(g):
        for o in cmp_tools:
            if (("dsrA" in repr_[o] and preds[o].get((g, "dsrA"), False)) or
                    ("dsrB" in repr_[o] and preds[o].get((g, "dsrB"), False))):
                return True
        return False
    misattributed_oxidizers = 0   # oxidizers a direction-blind tool would read as reducers
    for g in genomes:
        d = scycle_dsr_direction(preds, repr_, g)
        if d is None:
            continue
        direction[d] += 1
        dsrab_genomes.append((g, d))
        if d == "oxidative" and cmp_reports_dsrab(g):
            misattributed_oxidizers += 1
    n_dsrab = len(dsrab_genomes)

    # ── write report ─────────────────────────────────────────────────────
    labels = {}
    if args.panel_map and args.panel_map.exists():
        for r in csv.DictReader(open(args.panel_map), delimiter="\t"):
            samp = (r.get("sample") or r.get("ncbi_acc", "")).replace(".", "_")
            labels[samp] = r
    n_enriched = sum(1 for g in genomes if labels.get(g, {}).get("stratum") == "enriched")

    L = []
    tool_label = {"scycle": "scycle", "kofam": "raw KofamScan",
                  "metabolic": "METABOLIC", "scycdb": "SCycDB"}
    tools_str = ", ".join(tool_label.get(t, t) for t in TOOLS)
    met_note = ("METABOLIC included via the 380-reused + 119-batched GTDB-500 run (see WORKPLAN); "
                if "metabolic" in preds else "METABOLIC dropped at this scale (see WORKPLAN); ")
    L.append("# scycle-pipeline — GTDB-scale Cross-Tool CONCORDANCE\n")
    L.append(f"**Genomes:** {n_g} GTDB-representative ({n_enriched} S-cycle-clade-enriched + "
             f"{n_g - n_enriched} cross-phylum backbone) · **Tools:** {tools_str} · "
             "**No ground truth** — agreement/divergence only. "
             f"({met_note}DRAM excluded at scale per prereg.)\n")
    L.append("> Concordance, not accuracy: arbitrary GTDB genomes have no curated truth. "
             "High agreement on unambiguous markers + systematic divergence at the homology/"
             "direction traps is the expected signature; the trap loci are where scycle's gating "
             "and dsrAB-direction resolution act.\n")

    L.append("\n## 1. Pairwise agreement\n")
    L.append("Over targets both tools can represent. **Raw** = identical present/absent call "
             "(inflated by shared-absent cells). **Jaccard** = positive-call agreement "
             "|both present| / |either present| (the honest concordance metric). Reported for "
             "ALL targets and the homology-trap subset.\n")
    L.append("| tool pair | raw (ALL) | raw (trap) | Jaccard (ALL) | Jaccard (trap) |")
    L.append("|---|---|---|---|---|")
    for p in pairs:
        a_all, _ = agree_all[p]
        a_tr, _ = agree_trap[p]
        L.append(f"| {p[0]} vs {p[1]} | {a_all:.3f} | {a_tr:.3f} | "
                 f"{jac_all[p]:.3f} | {jac_trap[p]:.3f} |")

    cmp_blind = {"kofam": "raw KofamScan", "metabolic": "METABOLIC", "scycdb": "SCycDB"}
    blind_str = ", ".join(cmp_blind[t] for t in cmp_tools)
    L.append("\n## 2. dsrAB reductive↔oxidative DIRECTION trap (headline)\n")
    L.append("dsrA/dsrB are the SAME genes (KO K11180/K11181) in dissimilatory sulfate reducers "
             "(REDUCTIVE Dsr, SO₃²⁻→H₂S) and in reverse-Dsr sulfur oxidizers (OXIDATIVE rDSR, "
             "H₂S/S⁰→SO₃²⁻). Gene presence alone cannot separate the two metabolic directions — "
             "the genomic companions do (dsrD ⇒ reductive; Sox/sqr with dsrD absent ⇒ oxidative). "
             "scycle resolves direction (`apply_rules.resolve_dsr_direction`); the direction-blind "
             f"comparators ({blind_str}) report dsrA/dsrB presence with **no direction**.\n")
    L.append(f"- scycle calls **dsrAB** in **{n_dsrab}** genomes, resolving direction as "
             f"**{direction['reductive']} reductive** (sulfate reducers), "
             f"**{direction['oxidative']} oxidative** (reverse-Dsr sulfur oxidizers), "
             f"**{direction['ambiguous']} ambiguous**.")
    L.append(f"- Of the **{direction['oxidative']}** oxidative (reverse-Dsr) genomes, "
             f"**{misattributed_oxidizers}** have ≥1 direction-blind comparator "
             f"({blind_str}) reporting dsrA/dsrB present — i.e. they would be **mis-attributed "
             "as sulfate reducers** by a tool that reads dsrAB presence as dissimilatory sulfate "
             "reduction. This is the dsrAB-direction analogue of ncycle's nxrA→narG mis-routing, "
             "now observed at GTDB scale.\n")

    L.append("\n## 3. Targets where scycle is DISTINCTIVE (differs from every comparator)\n")
    L.append("| target | n genomes scycle-distinct | trap? |")
    L.append("|---|---|---|")
    for t, n in sorted(distinct.items(), key=lambda x: -x[1])[:15]:
        if n:
            L.append(f"| {t} | {n} | {'trap' if t in A.TRAP else ''} |")

    L.append("\n## 4. Per-target present-call counts (by tool; '-' = tool cannot represent)\n")
    L.append("Sorted by divergence (max−min present-count among representing tools).\n")
    L.append("| target | " + " | ".join(TOOLS) + " | divergence |")
    L.append("|---|" + "|".join("---" for _ in TOOLS) + "|---|")
    for row in sorted(per_target, key=lambda r: -r["divergence"]):
        if row["divergence"] == 0 and all(row[t] in (0, "-") for t in TOOLS):
            continue
        L.append("| " + row["target"] + " | "
                 + " | ".join(str(row[t]) for t in TOOLS)
                 + f" | {row['divergence']} |")

    args.out.write_text("\n".join(L) + "\n")
    # also a machine-readable per-target TSV
    tsv = args.out.with_suffix(".tsv")
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["target", "scycle", "kofam", "metabolic",
                                           "scycdb", "divergence"], delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for row in per_target:
            w.writerow(row)
    print(f"[concordance] {n_g} genomes, {len(targets)} targets")
    print(f"[concordance] dsrAB direction: {direction}; mis-attributable oxidizers: "
          f"{misattributed_oxidizers}")
    print(f"[concordance] wrote {args.out} + {tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
