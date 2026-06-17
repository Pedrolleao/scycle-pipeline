#!/usr/bin/env python3
"""
select_gtdb_scyc.py — pick the reproducible 500-genome GTDB slice for the SULFUR
cross-tool CONCORDANCE study (no ground truth; we measure where tools agree/disagree).

Mirrors ncycle's comparators/gtdb_pilot/select_gtdb_pilot.py. To keep the two sister
studies directly comparable, the 380-genome cross-phylum BACKBONE is REUSED VERBATIM
from ncycle's gtdb500 selection; only the ~120 sulfur-ENRICHED genomes differ (the
clades that populate the dsrAB direction trap + Sox/thiosulfate/DMSP loci).

Input : bac120_taxonomy.tsv + ar53_taxonomy.tsv (GTDB r232; reused from ncycle's gtdb_pilot)
        + ncycle's gtdb500/selection.tsv (for the shared backbone rows)
Output: selection.tsv (accession, ncbi_acc, stratum, clade, gtdb_taxonomy)

ENRICH clades verified against the GTDB taxonomy (2026-06-12, all ≥9 species-reps):
8 reductive-Dsr reducers + 8 rDsr/Sox oxidizers + 3 DMSP (Roseobacter clade) + archaea.
"""
from __future__ import annotations
import csv, random, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
NCYCLE = Path("/home/dmin/Grants/Nitrogen_Cycle/ncycle-pipeline/comparators")
BAC = NCYCLE / "gtdb_pilot" / "bac120_taxonomy.tsv"
AR = NCYCLE / "gtdb_pilot" / "ar53_taxonomy.tsv"
N_SELECTION = NCYCLE / "gtdb500" / "selection.tsv"   # reuse its 380 backbone rows
SEED = 1234
DEF_OUT = HERE / "selection.tsv"
DEF_PER_CLADE = 5
DEF_ENRICHED_CAP = 120        # cap enriched at this many; backbone fills to TARGET_TOTAL
DEF_TARGET_TOTAL = 500

# clade label -> GTDB lineage substring (verified by the selection specialist, 2026-06-12)
ENRICH = {
    # dissimilatory sulfate / sulfur reducers (reductive Dsr)
    "Desulfovibrio": "g__Desulfovibrio", "Desulfobacter": "g__Desulfobacter",
    "Desulfobulbus": "g__Desulfobulbus", "Desulfotomaculum": "g__Desulfotomaculum",
    "Desulfomicrobium": "g__Desulfomicrobium", "Desulfuromonas": "g__Desulfuromonas",
    "Thermodesulfobacterium": "g__Thermodesulfobacterium",
    "Thermodesulfovibrio": "g__Thermodesulfovibrio",
    # sulfur/sulfide oxidizers — reverse-Dsr (rDsr)
    "Sulfurovum": "g__Sulfurovum", "Sulfurimonas": "g__Sulfurimonas",
    "Chlorobiia_GSB": "c__Chlorobiia", "Chromatiales_PSB": "o__Chromatiales",
    "Thioalkalivibrio": "g__Thioalkalivibrio", "Thioglobus_SUP05": "g__Thioglobus",
    "Beggiatoa_Thiothrix": "g__Thiothrix",
    # sulfur/sulfide oxidizers — Sox pathway
    "Thiobacillus": "g__Thiobacillus", "Acidithiobacillus": "g__Acidithiobacillus",
    "Sulfurihydrogenibium": "g__Sulfurihydrogenibium",
    "Halothiobacillus": "g__Halothiobacillus", "Thiomonas": "g__Thiomonas",
    "Sulfurospirillum": "g__Sulfurospirillum",
    # thiosulfate / DMSP cyclers (Roseobacter clade)
    "Ruegeria_DMSP": "g__Ruegeria", "Roseobacter_DMSP": "g__Roseobacter",
    "Sulfitobacter_DMSP": "g__Sulfitobacter",
    # archaea
    "Archaeoglobus": "g__Archaeoglobus", "Sulfolobaceae_Sox": "f__Sulfolobaceae",
}


def ncbi_acc(gtdb_acc: str) -> str:
    return gtdb_acc.split("_", 1)[1] if gtdb_acc[:3] in ("RS_", "GB_") else gtdb_acc


def load(path: Path) -> list[tuple[str, str]]:
    out = []
    with open(path) as fh:
        for line in fh:
            acc, lin = line.rstrip("\n").split("\t")
            out.append((acc, lin))
    return out


def species_reps(rows):
    by_sp = defaultdict(list)
    for acc, lin in rows:
        by_sp[lin.split(";")[-1]].append((acc, lin))
    reps = {}
    for sp, cands in by_sp.items():
        cands.sort(key=lambda x: (0 if x[0].startswith("RS_") else 1, x[0]))
        reps[sp] = cands[0]
    return reps


def load_backbone(path: Path) -> list[dict]:
    """Reuse ncycle's 380 backbone rows verbatim (shared cross-phylum set)."""
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("stratum") == "backbone":
                out.append(r)
    return out


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEF_OUT)
    ap.add_argument("--per-clade", type=int, default=DEF_PER_CLADE)
    ap.add_argument("--enriched-cap", type=int, default=DEF_ENRICHED_CAP)
    ap.add_argument("--target-total", type=int, default=DEF_TARGET_TOTAL)
    args = ap.parse_args()

    rng = random.Random(SEED)
    reps = species_reps(load(BAC) + load(AR))
    print(f"[select] {len(reps)} species-reps", file=sys.stderr)

    chosen: dict[str, dict] = {}
    # ---- enriched stratum (sulfur clades) ----
    for clade, pat in ENRICH.items():
        if len(chosen) >= args.enriched_cap:
            break
        hits = sorted((acc, lin) for sp, (acc, lin) in reps.items() if pat in lin)
        rng.shuffle(hits)
        for acc, lin in hits[:args.per_clade]:
            if acc not in chosen and len(chosen) < args.enriched_cap:
                chosen[acc] = {"accession": acc, "ncbi_acc": ncbi_acc(acc),
                               "stratum": "enriched", "clade": clade, "gtdb_taxonomy": lin}
    n_enriched = len(chosen)
    print(f"[select] enriched: {n_enriched} sulfur genomes across {len(ENRICH)} clades", file=sys.stderr)

    # ---- backbone stratum: REUSE ncycle's 380 backbone rows verbatim ----
    backbone = load_backbone(N_SELECTION)
    added = 0
    for r in backbone:
        if len(chosen) >= args.target_total:
            break
        if r["accession"] not in chosen:
            chosen[r["accession"]] = r
            added += 1
    print(f"[select] backbone: {added} reused from ncycle (shared cross-phylum set)", file=sys.stderr)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "ncbi_acc", "stratum", "clade", "gtdb_taxonomy"],
                           delimiter="\t")
        w.writeheader()
        for rec in chosen.values():
            w.writerow(rec)
    print(f"[select] wrote {args.out} ({len(chosen)} genomes: {n_enriched} enriched + {added} backbone)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
