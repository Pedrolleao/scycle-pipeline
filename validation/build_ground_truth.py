#!/usr/bin/env python3
"""
build_ground_truth.py — KEGG-derived sulfur-cycle ground truth (kegg_v2) for the sulfur panel.
Panel grew SP1(22) → SP5/SP6b(25) → SP8(43, 2026-05-27); the KEGG_ORG map below is the
authoritative roster (currently n=43; 20 are INDEPENDENT — not a seed/HMM source for any target).
See prereg.md CHANGELOG + EXPANSION_PLAN.md for the panel-history disclosure.

SP1 (2026-05-25): replaces the hand-curated quick-v1 (which had curation errors and a
self-tuned 10-genome panel) with an AUTHORITATIVE, reproducible ground truth derived from
KEGG per-organism KO annotation. For every KO-anchored target in config/targets.yaml, a
genome is scored `present` iff KEGG assigns any of that target's KOs to the genome (REST
`link/ko/<org>`), else `absent`. The 4 KO-less targets (dsrC, dsrD, soxD, otr — no KOfam KO)
are handled explicitly (KO_LESS below).

Honesty notes / known GT limitations:
  * KEGG KO assignment uses KOfam best-hit, so for BROAD families it can call a KO on a
    paralog (e.g. K17218/sqr on a gor/lpdA-type Pyr_redox protein). Where the pipeline's
    BLAST gate then disqualifies such a hit, the resulting "FN" may be a GT over-call rather
    than a true miss — these surface in the per-target table and are the signal that the
    seed set needs broadening (S5/SP3 in ROADMAP). We accept KEGG as the reference and report
    the disagreements rather than hand-tuning them away.
  * Independent positives (genomes NOT used as curated seed sources): Afulgidus (dsrAB/aprAB/
    sat, archaeal), Tdenitrificans + Sdenitrificans (sox/sqr/reverse-dsr), Dshibae (dmdA),
    Styphimurium (phs/ttr). These test the seeds/HMMs without circularity.

Output: validation/ground_truth.tsv  (genome, target, expected, source)
Caches KEGG KO sets in validation/.kegg_cache/ (delete to refresh).
"""
from __future__ import annotations
import sys, urllib.request
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation" / "ground_truth.tsv"
CACHE = ROOT / "validation" / ".kegg_cache"
TARGETS = ROOT / "config" / "targets.yaml"

# panel genome (proteome stem) -> KEGG organism code (all validated, KO-links > 1000)
GENOME_ORG = {
    "Dvulgaris_Hildenborough": "dvu", "Pdenitrificans_PD1222": "pde",
    "Aferrooxidans_ATCC23270": "afe", "Ecoli_K12_MG1655": "eco",
    "Synechocystis_PCC6803": "syn", "Soneidensis_MR1": "son",
    "Wsuccinogenes_DSM1740": "wsu", "Bsubtilis_168": "bsu",
    "Lacidophilus_4356": "lba", "Spneumoniae_ref": "spn",
    "Afulgidus_DSM4304": "afu", "Tdenitrificans_ATCC25259": "tbd",
    "Sdenitrificans_DSM1251": "sun", "Dshibae_DFL12": "dsh",
    "Styphimurium_LT2": "stm", "Paeruginosa_PAO1": "pae",
    "Smeliloti_1021": "sme", "Npcc7120": "ana", "Scerevisiae_S288C": "sce",
    "Hpylori_26695": "hpy", "Avinelandii_DJ": "avn", "Rpalustris_CGA009": "rpa",
    # SP5 panel expansion (2026-05-26): orphan-target fillers, microbiologist-vetted.
    "Aambivalens_LEI10": "aamb",   # Acidianus ambivalens (archaeal S-oxidizer, Sulfolobales):
                                   #   type organism for sor/K16952, doxA/K16936, sreA/K17219
                                   #   (all KO-only → clean). NOTE: also a doxD/tth custom-HMM seed
                                   #   source, so its doxD/tth calls are non-independent.
    "Snovella_DSM506": "sno",      # Starkeya novella (= Ancylobacter novellus DSM 506): textbook
                                   #   SorAB sulfite dehydrogenase → sorA/K05301 positive (KO-only,
                                   #   clean). NOTE: a soxB custom-HMM seed source → soxB non-independent.
    # SP6b: dedicated dddP-positive genome (microbiologist-vetted).
    "Rdenitrificans_OCh114": "rde",  # Roseobacter denitrificans OCh 114: RdDddP (Q166H0) is a
                                   #   crystallized, characterized M24B DMSP lyase → dddP positive via
                                   #   function-correction (KEGG lacks K28073). Excluded from DddP HMM/seeds.
    # SP8 expansion (2026-05-27): INDEPENDENT positives added to fix trap-claim circularity (audit).
    # All assembly-matched (KEGG genome == UniProt proteome) + microbiologist-vetted. See EXPANSION_PLAN.md.
    "Msedula_DSM5348": "mse",      # Metallosphaera sedula DSM 5348 (GCA_000016605.1): archaeal S-oxidizer.
                                   #   INDEPENDENT doxD + doxA (the doxD trap had ZERO independent positives).
                                   #   NB: NOT independent for tth (tth HMM/refs include Metallosphaera) — tth not in TRAP.
    "Doleivorans_Hxd3": "dol",     # Desulfosudis (Desulfococcus) oleivorans Hxd3: independent reductive SRB
                                   #   (dsrAB/dsrD/aprA/sat) + independent phsA. Not a seed source.
    "Dbaculatum_DSM4028": "dba",   # Desulfomicrobium baculatum DSM 4028: independent reductive SRB
                                   #   (dsrAB/dsrD/aprA/sat) + independent sor + phsA. Not a seed source.
    "Dautotrophicum_HRM2": "dat",  # Desulforapulum (Desulfobacterium) autotrophicum HRM2: independent
                                   #   reductive SRB (dsrAB/dsrD/aprA/sat). Not a seed source.
    "Carsenatis_LY1": "cars",      # Citrobacter arsenatis LY-1: independent enteric phsA + ttrA (non-Salmonella).
    "Pmirabilis_HI4320": "pmr",    # Proteus mirabilis HI4320: independent enteric phsA + ttrA.
    # P2 expansion (2026-05-27): independent dsr (SRB) + Sox positives to power the DRAM trap-precision
    # contrast. Microbiologist-vetted; all assembly-matched. See EXPANSION_PLAN.md (P2).
    "Dpropionicus_DSM2032": "dpr",   # Desulfobulbus propionicus: independent reductive SRB (dsrAB/dsrD/apr).
    "Dacetoxidans_DSM11109": "dao",  # Desulfobacca acetoxidans: independent reductive SRB (dsrAB/dsrD/apr).
    "Dsulfexigens_DSM10523": "dsf",  # Desulfocapsa sulfexigens: S-disproportionator, reductive Dsr (dsrAB/dsrD).
    "Tsulfidiphilus_HLEbGr7": "tgr", # Thioalkalivibrio sulfidiphilus: SOB, INDEPENDENT oxidative (reverse) dsrAB,
                                     #   dsrD-absent. NB: no SoxCD → NOT a soxD positive (vetting-confirmed).
    "Tviolascens_DSM198": "tvi",     # Thiocystis violascens: purple-S SOB, INDEPENDENT oxidative dsrAB, dsrD-absent.
                                     #   No SoxCD → not soxD. Is a soxB-expanded-set seed (soxB not a TRAP target).
    "Tarsenitoxydans_3As": "thi",    # Thiomonas arsenitoxydans 3As: Sox oxidizer carrying SoxD → independent soxD.
    # P4 expansion (2026-05-28): mixed batch — orphan reinforcements (sorA, sor, sreA, dddP), Mo-bis-MGD
    # decoy + cross-lineage Sox tests. Microbiologist-vetted (DOSSIER_P4_AUDIT 6/6 GO with named caveats).
    "Bjaponicum_USDA110": "bja",   # Bradyrhizobium diazoefficiens USDA 110 (GCA_000011365.1):
                                   #   independent SorA (K05301) reinforcement + Sox + sqr + fccB.
                                   #   Only genus-level soxC seed overlap (B. brasilense, different species).
    "Aaeolicus_VF5": "aae",        # Aquifex aeolicus VF5 (GCA_000008625.1): hyperthermophilic Aquificota SOB.
                                   #   Independent sor (K16952) reinforcement. NB: sqr seed self-hit (O67931
                                   #   is from this strain) → sqr cell must be flagged non-independent.
    "Atumefaciens_H13-3": "agr",   # Agrobacterium tumefaciens H13-3 (GCA_000192635.1): independent dddP
                                   #   reinforcement (Rhizobiaceae, not Roseobacter clade). CLEAN — no seed overlap.
    "Cnecator_H16": "reh",         # Cupriavidus necator H16 (Ralstonia eutropha; GCA_000009285.2):
                                   #   Mo-bis-MGD denitrifier decoy (narG/napA/dmsA but no phsA/ttrA/sreA).
                                   #   Also carries sorA + Sox + fccB. Audit: most likely FP route = dmsA paralog.
    "Hydrogenobaculum_Y04AAS1": "hya", # Hydrogenobaculum sp. Y04AAS1 (GCA_000020785.1): Aquificales SOB at
                                   #   near-neutral pH. ttrA (K08357) + Sox + sqr cross-lineage test. Audit
                                   #   note: ttrA call is KEGG-pipeline assignment, not biochemically characterized
                                   #   → ideal "Mo-bis-MGD stress test" candidate, do not present as validated.
    "Sazoricus_FC6": "sazo",       # Stygiolobus azoricus FC6 (GCA_009729035.1): obligate-anaerobe S⁰-respiring
                                   #   Sulfolobales archaeon. Independent sreA (K17219) reinforcement (non-Acidianus
                                   #   genus, distinct from aamb). NB: NOT independent for doxD/tth (genus-level seeds).
}

# Independent positives — for the report. NB: genome-level independence is too coarse; the
# rigorous, per-(genome,target) notion is in validation/trap_independence.py (use that for the
# trap claim). CORRECTION (2026-05-27 audit): Afulgidus is a dsrA/dsrB/dsrD seed source and was
# wrongly listed here; it is independent for sox/sulfur-oxidation but NOT for the dissimilatory
# dsr trap. These genomes are not a curated seed/HMM source for ANY target (truly independent):
INDEPENDENT = {"Tdenitrificans_ATCC25259", "Sdenitrificans_DSM1251",
               "Dshibae_DFL12", "Soneidensis_MR1",
               # SP8 independent positives (not a seed/HMM source for any target):
               "Doleivorans_Hxd3", "Dbaculatum_DSM4028", "Dautotrophicum_HRM2",
               "Carsenatis_LY1", "Pmirabilis_HI4320",
               # P2 independent positives (not a seed source for the trap target they provide):
               "Dpropionicus_DSM2032", "Dacetoxidans_DSM11109", "Dsulfexigens_DSM10523",
               "Tsulfidiphilus_HLEbGr7", "Tarsenitoxydans_3As",
               # P4 independent positives (clean for the trap target they provide; see DOSSIER_P4_AUDIT
               # for per-candidate seed-overlap notes — aae sqr + sazo doxD/tth genus-level overlaps are
               # NOT independent; the targets each genome IS added to reinforce are independent):
               "Bjaponicum_USDA110", "Aaeolicus_VF5", "Atumefaciens_H13-3",
               "Cnecator_H16", "Hydrogenobaculum_Y04AAS1", "Sazoricus_FC6"}
# (Msedula_DSM5348 indep for doxD/doxA but a tth seed; Tviolascens_DSM198 indep for its dsr/soxD-absent
#  roles but a soxB-expanded seed — both excluded from this fully-independent reporting set.)
# Independent FOR THE DISSIMILATORY MODULE (sox/rDsr oxidizers, not dsr-seed sources): Thiobacillus,
# Sulfurimonas. The original dissim "independent positive" framing leaned on Archaeoglobus, which is
# a dsr seed source — see trap_independence.py for the per-target breakdown.

# Sulfate reducers (reductive Dsr → carry dsrD); sulfur oxidizers run dsr in reverse (no dsrD).
# SP8 (2026-05-27): added the 3 new INDEPENDENT reductive SRBs so their dsrD scores "present"/reductive
# (they carry dsrA but, without listing here, dsrD would be left UNSCORED — flagged in vetting).
SULFATE_REDUCERS = {"Dvulgaris_Hildenborough", "Afulgidus_DSM4304",
                    "Doleivorans_Hxd3", "Dbaculatum_DSM4028", "Dautotrophicum_HRM2",
                    # P2: reductive SRBs (DsrD present, vetting-confirmed)
                    "Dpropionicus_DSM2032", "Dacetoxidans_DSM11109", "Dsulfexigens_DSM10523"}
SOX_REVERSE_DSR = {"Tdenitrificans_ATCC25259", "Sdenitrificans_DSM1251",
                   # P2: sulfur oxidizers running reverse-Dsr (DsrD absent, vetting-confirmed)
                   "Tsulfidiphilus_HLEbGr7", "Tviolascens_DSM198"}


def kegg_ko_set(org: str) -> set[str]:
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{org}.ko"
    if not f.exists():
        url = f"https://rest.kegg.jp/link/ko/{org}"
        print(f"  fetching KEGG KO set: {org}", file=sys.stderr)
        with urllib.request.urlopen(url, timeout=60) as r:
            f.write_bytes(r.read())
    kos = set()
    for line in f.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) == 2 and parts[1].startswith("ko:"):
            kos.add(parts[1][3:])
    return kos


def main() -> None:
    cfg = yaml.safe_load(open(TARGETS))
    targets = cfg["targets"]
    ko_by_target = {t["id"]: list(t.get("ko") or []) for t in targets}
    ko_less = [t["id"] for t in targets if not t.get("ko")]   # dsrC, dsrD, soxD, otr

    # Fetch KEGG KO sets for every genome.
    org_kos = {g: kegg_ko_set(org) for g, org in GENOME_ORG.items()}

    rows = [("genome", "target", "expected", "source")]
    for g in GENOME_ORG:
        kos = org_kos[g]
        # KO-anchored targets: present iff KEGG assigns any of the target's KOs.
        for t in targets:
            tid = t["id"]
            if tid in ko_less:
                continue
            present = any(k in kos for k in ko_by_target[tid])
            rows.append((g, tid, "present" if present else "absent", "kegg_v2"))
        # KO-less targets, derived explicitly:
        has_dsrA = "K11180" in kos
        has_soxC = "K17225" in kos
        # dsrC: co-occurs with dsrAB (both reductive and reverse Dsr).
        rows.append((g, "dsrC", "present" if has_dsrA else "absent", "kegg_v2"))
        # dsrD: reductive-direction marker — present in sulfate reducers, absent in SOB
        # running reverse Dsr; only scored where the direction is unambiguous.
        if g in SULFATE_REDUCERS:
            rows.append((g, "dsrD", "present", "kegg_v2"))
        elif g in SOX_REVERSE_DSR or not has_dsrA:
            rows.append((g, "dsrD", "absent", "kegg_v2"))
        # (genomes with dsrA but ambiguous direction → dsrD left unscored)
        # soxD: co-occurs with the SoxCD dehydrogenase module (soxC present).
        rows.append((g, "soxD", "present" if has_soxC else "absent", "kegg_v2"))
        # otr: octaheme tetrathionate reductase — rare; Shewanella oneidensis is the
        # characterized panel positive, the rest are scored absent.
        rows.append((g, "otr", "present" if g == "Soneidensis_MR1" else "absent", "kegg_v2"))

    with open(OUT, "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    n_p = sum(1 for r in rows[1:] if r[2] == "present")
    n_a = sum(1 for r in rows[1:] if r[2] == "absent")
    print(f"[ground_truth] {len(rows)-1} cells across {len(GENOME_ORG)} genomes "
          f"({n_p} present, {n_a} absent; {len(INDEPENDENT)} independent-positive genomes) -> {OUT}")


if __name__ == "__main__":
    main()
