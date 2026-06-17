#!/usr/bin/env python3
"""
build_curated_function_gt.py — curated-FUNCTION ground truth for the sulfur panel
(n=43 after the SP8 expansion; inherits the genome roster from build_ground_truth.py).

WHY (SP4): the headline benchmark question is "does the pipeline beat raw KofamScan?".
Our default ground truth (`ground_truth.tsv`) is KEGG-KO-derived, so a raw-KO comparator
is near-circular (raw KO ≈ the GT) and any homology-trap advantage is invisible. This file
is the NON-circular reference: it starts from the KEGG-KO GT and applies a small set of
**function-level corrections on the homology-trap cells where KO-presence ≠ functional role**,
each verified independently (UniProt protein annotation, operon/subunit context, or
characterized biology of the reference strain) — NOT copied from either tool's output.

Curation discipline (mirrors the seed-curation rule): every correction is justified by a
verified protein FUNCTION or subunit/operon context, recorded in CORRECTIONS below with its
evidence — NOT by scycle's HMM/BLAST score. HONESTY NOTE (2026-05-27 audit): the SP4b/SP6b
corrections did cut both ways, but the SP6c specialist audit REVERTED the fccA corrections (the
only ones that penalised the pipeline), so the 3 corrections that REMAIN (doxD, sdo, dddP) are
all `absent→present` and all favour the pipeline. They are individually defensible as KEGG
under-annotation fixes (each is a characterized protein KEGG simply failed to assign a KO), but
the aggregate is one-directional — so the KEGG-KO GT (`ground_truth.tsv`) is the conservative
reference and this curated GT is best read as a sensitivity analysis, not the sole headline.

Output: validation/curated_function_gt.tsv (genome, target, expected, source, rationale)
The non-corrected cells inherit the KEGG-KO GT verbatim (source `kegg_v2`); corrected cells
carry source `curated_function` + a rationale. Run AFTER build_ground_truth.py.
"""
from __future__ import annotations
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
KEGG_GT = ROOT / "ground_truth.tsv"
OUT = ROOT / "curated_function_gt.tsv"

# ── function-level corrections vs the KEGG-KO GT ─────────────────────────────
# (genome, target): (expected, rationale).  Each verified independently of the tools.
CORRECTIONS: dict[tuple[str, str], tuple[str, str]] = {
    # doxD — thiosulfate:quinone oxidoreductase (TQO). KEGG is sparse for K16937; the
    # enzyme is biochemically characterised in Acidithiobacillus. A. ferrooxidans' own
    # DoxD is curated seed B7J3E9 (100% self-identity) and fires the DoxD custom HMM (e-186).
    ("Aferrooxidans_ATCC23270", "doxD"):
        ("present", "TQO/DoxD characterised in Acidithiobacillus; UniProt B7J3E9 is A.fer's own "
                    "DoxD (100% BLAST self-hit) + custom-HMM e-186; KEGG lacks K16937 (under-annotation)"),
    # fccA — REVERTED after the SP6c specialist audit (was 5 'present' corrections in SP4b/SP6b).
    # fccA is the diheme cyt-c subunit of flavocytochrome-c sulfide dehydrogenase (FccAB). The SP4b
    # corrections inferred fccA from fccB (K17229) presence, but the audit showed K17229 is also
    # assigned to the Sox-system soxF flavoprotein — so "fccB present" does NOT imply a cognate FccA
    # (verified: R. denitrificans' K17229 hit is soxF in a soxCDEF cluster, no FccA). FccA cannot be
    # distinguished from Sox / other c-cytochromes without operon synteny (unavailable on proteome
    # input), so fccA is left at the objective KEGG-KO baseline (absent for the panel's sox organisms)
    # rather than asserting an undetectable/uncertain positive. Genuine FccAB exists in the panel
    # (P. denitrificans, operon-verified Pden_4157/4158) but is a documented synteny-limited miss.
    # sdo — persulfide/sulfur dioxygenase. KEGG under-annotates K17725 in Acidithiobacillus, but
    # A. ferrooxidans' Sdo hits the curated Acidithiobacillus/Methylococcus "Sulfur dioxygenase"
    # seeds at 96% identity (seeds are cross-species: A. caldus / A. ferruginosus, NOT A. ferrooxidans),
    # and Acidithiobacillus is the characterized SDO genus (microbiologist-confirmed). The ~52-55%
    # "MBL fold metallo-hydrolase" hits in Nostoc/Synechocystis/T. denitrificans are left ABSENT
    # (glyoxalase-II-class paralogs; gate raised to ≥60% drops them).
    ("Aferrooxidans_ATCC23270", "sdo"):
        ("present", "Genuine persulfide dioxygenase: 96% BLAST to characterized Acidithiobacillus Sdo "
                    "seeds (cross-species, non-circular); Acidithiobacillus is the characterized SDO genus"),
    # dddP — Roseobacter denitrificans RdDddP (UniProt Q166H0) is a crystallized, experimentally
    # characterized M24B DMSP lyase, but KEGG does not assign K28073 to rde (M24 moonlighting family,
    # sparse KO). Detected by the dedicated DddP custom HMM (663 >> TC 277.5) + ≥50% DddP BLAST gate
    # (rde DddP is 77% to the non-rde seeds; the HMM/seeds exclude rde → non-circular). Function GT = present.
    ("Rdenitrificans_OCh114", "dddP"):
        ("present", "RdDddP (Q166H0) — crystallized, characterized M24B DMSP lyase; KEGG lacks K28073 "
                    "for rde (sparse KO). Detected via DddP custom HMM + BLAST gate, both excluding rde (non-circular)"),
    # dddP / Agrobacterium tumefaciens H13-3 — KEGG K28073 → present, but BBH cross-check (post-P4 audit,
    # 2026-05-28) revealed the K28073-annotated protein (WP_013637020.1) is 832 aa (vs canonical DddP
    # ~440 aa) and has ZERO significant homology to characterized DddP references (RdDddP/Q166H0,
    # Roseovarius/Rhodobacter/Cognatiyoonia/Vannielia DddPs — best hit at e<100 is 31% over 42 aa = noise).
    # The KEGG K28073 → agr assignment is a misannotation in the M24 metallopeptidase family. scycle's
    # BLAST gate correctly disqualified the hit (was a FN under KEGG-KO GT; this correction makes it TN).
    ("Atumefaciens_H13-3", "dddP"):
        ("absent", "BBH cross-check (2026-05-28): K28073-annotated WP_013637020.1 is 832 aa (canonical "
                   "DddP ~440 aa) with no significant homology to characterized DddP references; M24 "
                   "family misannotation. scycle's BLAST gate correctly disqualified; not a real DddP"),
    # ── cells DELIBERATELY left at the KEGG-KO call (documented, not corrected) ──
    #  phsA / D. vulgaris: KEGG GT already 'present'; DVU_0173 = UniProt 'Thiosulfate reductase,
    #    putative' (Mo-bis-MGD). Kept present (function GT) even though no phsBC operon partners and
    #    the pipeline gates it (a pipeline FN) — honest, medium-confidence.
    #  phsA / A. fulgidus: KEGG absent, no K08352 in KEGG, no phsBC operon → function absent (kept).
    #  fccA / 5 fccB-negative genomes (Aferrooxidans, Avinelandii, Paeruginosa, Smeliloti,
    #    Soneidensis): no catalytic flavoprotein → no FccAB → function absent = KEGG call (kept).
    #  dsrC/dsrD/soxD/otr: KEGG GT already encodes the functional call via explicit subunit/operon
    #    rules (build_ground_truth.py); function GT agrees, no correction.
}


def main() -> None:
    rows = list(csv.DictReader(open(KEGG_GT), delimiter="\t"))
    n_corr = 0
    out = [("genome", "target", "expected", "source", "rationale")]
    for r in rows:
        key = (r["genome"], r["target"])
        if key in CORRECTIONS:
            exp, why = CORRECTIONS[key]
            assert exp != r["expected"], f"no-op correction for {key} (already {exp})"
            out.append((r["genome"], r["target"], exp, "curated_function", why))
            n_corr += 1
        else:
            out.append((r["genome"], r["target"], r["expected"], r.get("source", "kegg_v2"), ""))
    with open(OUT, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(out)
    n_p = sum(1 for r in out[1:] if r[2] == "present")
    print(f"[curated_function_gt] {len(out)-1} cells ({n_p} present); "
          f"{n_corr} function-level corrections vs KEGG-KO GT -> {OUT}")
    for (g, t), (exp, _) in CORRECTIONS.items():
        print(f"    correction: {g:28} {t:6} -> {exp}")


if __name__ == "__main__":
    main()
