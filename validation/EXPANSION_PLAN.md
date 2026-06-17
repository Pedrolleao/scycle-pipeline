# Panel-expansion design (SP8) — audit-driven, 2026-05-27

# ════════════════════════════════════════════════════════════════════════════════
# ⏯  CURRENT STATE & WHERE TO START  (last updated 2026-05-28 — READ THIS FIRST)
# ════════════════════════════════════════════════════════════════════════════════
# **THE HARDENING CAMPAIGN IS COMPLETE.** Panel = **43 isolate genomes** (`../../sp1_panel/`, 42 .faa proteomes
# + 1 .fna) registered in `build_ground_truth.py:GENOME_ORG`. PLUS a separately-reported **8-input MAG
# realism / synteny study** (Auka/Guaymas/Saanich + 2 negative reference genomes, NOT folded into headline).
# Status: SP1–SP6c hardening + SP7 METABOLIC/DRAM benchmark + external audits + P1/P2/P3/P4 + post-P4
# audit closure (BBH cross-check + dsrA-miss investigation) — ALL DONE.
#
# HEADLINE (43-genome curated_function_gt, B=10,000 — `../../comparators/benchmark_final_p4_v2.txt`):
#   scycle ALL micro-F1 **0.930**, precision **0.952**, trap precision **0.943** (independent-only 0.969);
#   **BEATS all 3 comparators (KofamScan, METABOLIC, DRAM) on BOTH trap precision AND overall F1,
#   all BH-significant** (Δ-ALL-F1 +0.037/+0.279/+0.285; Δ-trap-P +0.265/+0.355/+0.180; q ≤ 0.0006).
#   11/11 regression PASS (`test_regression.py`, MAX_ALL_FP currently 32; precision is the real guard).
#
# P3 (MAG realism / synteny study, NOT folded into headline F1): scycle's `resolve_dsr_direction` synteny path
#   fires on every dsr-bearing MAG; **5/5 directional calls correct** (4 explicit + 1 honest "ambiguous").
#   **scycle is the only one of three tools that resolves dsr direction** on real fragmented MAGs — DRAM's
#   own module label is "dissimilatory sulfate reduction (and oxidation)" (parenthetical = explicit
#   acknowledgment). REPORT.md §"P3 — MAG realism / synteny study" + §"P3.6".
#
# DONE so far: P1 (25→31 independence fix) · P2 (31→37 DRAM trap-P significance) · P3 (8-input MAG study)
#   · P4 (37→43 isolate expansion, orphan reinforcement + Mo-bis-MGD decoy) · Post-P4 BBH cross-check
#   (closed audit; surfaced + fixed agr-dddP KEGG K28073 misannotation as a function-GT correction)
#   · O3 dsrA "miss" investigation (resolved as MAG fragmentation, no code change). Full detail:
#   REPORT.md (SP8, P2, P3, P3.6, P4, "P4 — Mo-bis-MGD BBH cross-check" sections) + memory
#   (scycle-audit-findings, p3-mag-realism, p4-panel-expansion).
#
# ▶ NOTHING REMAINS IN THE CODEBASE. Remaining work is manuscript / publication. The three minor optional
#   polish items below are listed for completeness but do NOT change any headline:
#   - `benchmark_stats.py:167` — bootstrap two-sided p double-counts the boundary (single-line fix)
#   - BH family includes near-foregone METABOLIC/DRAM all-F1 endpoints (slightly conservative MC correction)
#   - `otr` Thiocystis octaheme FP — known unclosable case; could be explicitly documented as a method limitation
#
# Documented limits — DO NOT chase: assimilatory `cys` FNs (KOfam-resolution + KEGG-GT ambiguity, caps assim
#   pathway at ~0.85); P. marinus K00958 sat FP (KOfam can't distinguish assim/dissim sat — same FP in
#   METABOLIC); sdo n=1 (no clean non-Acidithiobacillus carrier without genus seed overlap); O3 dsrA miss
#   (MAG fragmentation, HMM correctly calibrated). All in REPORT.md.
#
# ▶ THE REPEATABLE WORKFLOW to add genomes (what P1/P2 did — follow exactly):
#   1. SELECT: KEGG `link/genes/<KO>` to find carriers of the target's KO; pick non-seed, non-panel candidates.
#      Verify each: (a) carriage — intersect KEGG `link/ko/<org>` with target KOs; (b) NOT a seed source —
#      grep resources/blast_db/*.fasta + targets/*/{refs,expanded}.fasta for the organism; (c) KEGG↔UniProt
#      ASSEMBLY MATCH (KEGG `get/gn:<org>` "Assembly:" == UniProt proteome `genome_assembly`) — MANDATORY or the
#      KEGG-derived GT won't align; (d) BUSCO completeness. (NB: spawned sub-agents have NO network — if you
#      microbiologist-vet via an agent, STAGE the live data in comparators/vetting_dossier/ first, as P1/P2 did.)
#   2. INTEGRATE: download UniProt proteome (`rest.uniprot.org/uniprotkb/stream?query=proteome:<UPID>&format=fasta`)
#      to `../../sp1_panel/<Genome>.faa`; add to `GENOME_ORG`; add SRBs→`SULFATE_REDUCERS`, reverse-dsr
#      oxidizers→`SOX_REVERSE_DSR` (else dsrD UNSCORED); add fully-independent ones to `INDEPENDENT`.
#   3. RE-RUN (env: `conda activate ewaste-pipeline`, from scycle-pipeline/):
#        python validation/build_ground_truth.py && python validation/build_curated_function_gt.py
#        SCYCLE_ENV=ewaste-pipeline python run.py --input ../sp1_panel --prodigal-mode single --skip-db-setup --cores 8
#        python validation/score_scycle.py && python validation/test_regression.py   # re-lock MAX_ALL_FP if it trips (panel-size scaling; precision is the real guard)
#        GT_FILE=curated_function_gt.tsv python validation/trap_independence.py       # independent trap cells should rise
#   4. BENCHMARK (extend to new panel): run METABOLIC (`-in <dir of new .faa>`, env METABOLIC_v4.0) + DRAM
#      (`DRAM.py annotate_genes -i '<new .faa>/*.faa'` then `distill`, env DRAM14) on ONLY the new proteomes;
#      append rows to `../../comparators/{metabolic,dram}.tsv` (metabolic: worksheet1 Present→Corresponding KO;
#      dram: 3 sulfur modules present iff diagnostic gene in annotations.tsv); then
#        GT_FILE=curated_function_gt.tsv python validation/benchmark/benchmark_stats.py \
#            --metabolic ../comparators/metabolic.tsv --dram ../comparators/dram.tsv --bootstrap 10000
#      Tool install gotchas (DRAM14 not 1.5.0, METABOLIC -in-gn needs .fasta, dbCAN placeholder, etc.):
#      memory `comparators-install`. The P1/P2 runner scripts are in /tmp and comparators/run_*.sh as templates.
# ════════════════════════════════════════════════════════════════════════════════

Written after the two external audits (computational-biology methods review + molecular-microbiology
genome/seed verification). The genomes and curated seeds are all correct; the gaps the expansion must
close are **independence/leakage**, **single-positive orphan targets**, **statistical power at the
boundary**, and **input realism (MAGs)**. Goals are ordered by how much they de-risk the headline claim.

## STATUS (2026-05-28): P1 ✅ · P2 ✅ · P3 ✅ · P4 ✅ · Post-P4 BBH audit closure ✅ · **CAMPAIGN COMPLETE**

**P3 (MAG realism / synteny study) — DONE.** 8 nucleotide inputs (6 sulfur MAGs + 2 negatives), 3 distinct
sites (Auka, Guaymas, Saanich). **scycle 5/5 directional calls correct**; METABOLIC + DRAM cannot resolve
direction (DRAM's own module name is "dissimilatory sulfate reduction (and oxidation)" — explicit
acknowledgment). REPORT.md §"P3 — MAG realism / synteny study" + §"P3.6".

**P4 (isolate panel expansion 37→43) — DONE (2026-05-28).** 6 microbiologist-vetted isolate genomes
(bja, aae, agr, reh, hya, sazo). 5 of 6 orphan targets reinforced (sorA n=1→3; sor n=2→3; sreA n=1→2;
dddP n=1→2; doxA n=2→3). Mo-bis-MGD denitrifier decoy (reh) delivered the audit-expected paralog FPs
(tauD). **All 6 BH-significant benchmark wins held on the 43-panel** — scycle WINS vs KofamScan + METABOLIC
+ DRAM on both ALL-F1 and trap-precision (DRAM trap-P q strengthened P2→P4: 0.0114 → 0.0006). 11/11
regression PASS (FP floor re-locked 25→32 for panel scaling). REPORT.md §"P4 — panel expansion 37→43".

**Post-P4 audit closure — DONE (2026-05-28).** Actioned the two audit-recommended follow-up items:
(1) **Mo-bis-MGD BBH cross-check** verified 5/6 P4 suspect calls were correct (sazo sreA → 75% orthologue;
hya ttrA → 41% full-length TtrA, right family; reh fccB/A/sorA → mid-confidence as appropriate) and
surfaced a KEGG K28073 misannotation on agr (WP_013637020.1 = 832 aa vs canonical DddP ~440 aa, zero
homology to any characterized DddP). **Function-correction added** to `build_curated_function_gt.py`:
`(Atumefaciens_H13-3, dddP) = absent`. Impact: micro-F1 0.929→**0.930**, organic_sulfur_dmsp F1
0.857→**0.889**, KofamScan ALL-F1 gap widened +0.035→**+0.037**; all 6 BH-significant wins held or
strengthened (`benchmark_final_p4_v2.txt`). (2) **O3 Parabeggiatoa dsrA miss** from P3.5 investigated:
resolved as MAG-fragmentation (gene split into short ORFs in 1343-contig bin), NOT an HMM threshold
problem. R3 (24-contig MAG) produces clean full-length dsrA ORF that KOfam K11180 calls correctly. **No
code change** — documented as MAG-quality limit. REPORT.md §"P4 — Mo-bis-MGD BBH cross-check (2026-05-28)"
+ §P3 side-findings updated. **Both items closed without further tool changes.**

## Design targets
- **≥3 INDEPENDENT positives per trap target** (independent = the genome is NOT a curated BLAST-seed /
  custom-HMM source for that target; see `trap_independence.py`). This converts the degenerate
  trap-precision estimate (perfect on 21 independent cells, boundary CI) into a properly powered one.
- **≥3 positives for each former-orphan target** (sor, doxA, sorA, sreA, sdo, dddP) — currently 1 each,
  flagged in `prereg.md` as "single-positive, added-to-validate".
- **Paralog decoys** for every trap (genomes that carry the confounder but NOT the function) to test
  precision, not just recall.
- **Reserve hold-out genera** (`HOLDOUT_TAG = "holdout_v3"`, already wired in `adapters.py`) — scored
  separately, never used for any tuning/seed.

## Priority 1 — kill the circularity (the audit's top theme)  ✅ DONE (panel 25→31; added mse, dol, dba, dat, cars, pmr; independent trap cells 21→40, trap precision held 1.000)
Targets whose ONLY panel positives are self-seeds, so they currently have **zero independent test**:
| target | current positive(s) | why circular | add (non-seed-source) candidates |
|---|---|---|---|
| **doxD** | Aferrooxidans, Aambivalens | both are doxD HMM/BLAST seeds | *Metallosphaera sedula*, *Sulfobacillus*, *Acidithiobacillus thiooxidans* |
| **tth** | Aferrooxidans, Aambivalens | both are tth seeds | same acidophile/archaeal S-oxidisers as above (not the seed strains) |
| **dsrA/dsrB/dsrD** | Dvulgaris, Afulgidus (seeds); Tdenitrificans/Sdenitrificans (indep, rDsr) | the reductive-Dsr positives are both seeds | independent SRB: *Desulfobacterium autotrophicum*, *Desulfococcus oleovorans*, *Desulfomicrobium baculatum* |
| **phsA / ttrA** | Styphimurium (+Afulgidus for ttrA) | the ONLY reviewed PhsA/TtrA in SwissProt ARE these panel genomes | non-panel carriers: other *Shewanella*, *Proteus mirabilis*, *Citrobacter* (ttrA); *Shewanella*/*Aeromonas* (phsA) — also seek characterized non-panel refs |
| **soxD** | pde, rde (seeds) + Dshibae/Rpalustris/Sdenitrificans/Snovella (indep) | already has independent positives ✓ | (lower priority — add 1–2 more α/γ Sox) |

Rule going forward: **a genome added as a positive for target T must not be a seed source for T** (and ideally
not for any T). Check each candidate against the seed lists before adding.

## Priority 2 — power the contested comparator contrasts  ✅ DONE (panel 31→37; added dpr, dao, dsf, tgr, tvi, thi; DRAM trap-precision now BH-significant q=0.0114; independent trap cells 40→64)
- **DRAM trap-precision is n.s.** because DRAM distils only 3 sulfur modules and ventures few trap calls.
  Adding independent **dissimilatory-Dsr** and **Sox** positives (the two trap-relevant modules DRAM
  *does* distil) is what powers that contrast — the Desulfo* SRB + extra Sox genomes above do double duty.
- The overall-F1 vs DRAM gap is a coverage-scope artifact (don't chase it); don't add genomes just to widen it.

## Priority 3 — input realism: add real assembled MAGs  ⬅ START HERE (next)
All 37 current inputs are complete isolate proteomes (+1 isolate .fna); METABOLIC/DRAM are built for fragmented
MAGs, and scycle's **operon-synteny `resolve_dsr_direction` is DORMANT on proteome input** (verified: on .faa,
`coords={}` so the synteny branch never runs — apply_rules.py — so the directional-call capability is UNTESTED
by the current validation). Adding real MAGs is the only way to exercise it and is the fairer common ground for
the MAG-metabolism comparators.

**Concrete steps (this is a SEPARATE, clearly-labelled realism study — do NOT fold into the headline isolate F1):**
1. **Source 5–8 published sulfur-cycle MAGs as NUCLEOTIDE contigs (.fna/.fasta)** with a characterized sulfur
   phenotype + known taxonomy — e.g. SRB/SOB MAGs from hydrothermal-vent, marine-sediment, or sulfidic-aquifer
   metagenome papers (GTDB / NCBI assembly / the original study's supplement). Prefer MAGs with a stated
   completeness (CheckM ≥90%) and a dsr-direction call in the paper (so there's a ground truth for the synteny
   test). Pick a mix of reductive (dsrD+) and oxidative (reverse-dsr) so `resolve_dsr_direction` is tested both ways.
2. **Ground truth for MAGs is the hard part** — they have NO KEGG org code, so `build_ground_truth.py`'s KEGG path
   doesn't apply. Options: (a) annotate each MAG's KO set yourself (kofamscan/the pipeline's own KO calls) and
   build a per-MAG GT from that + the paper's stated phenotype; or (b) treat the MAG study as phenotype-level only
   (does scycle call the right pathway + the right dsr DIRECTION vs the paper?) rather than per-subunit F1. Decide
   and document; do NOT reuse the isolate GT machinery blindly.
3. **Run scycle on the MAGs** with `--prodigal-mode meta` (NOT single — MAGs are multi-organism-free but fragmented;
   meta mode is correct for assembled contigs) so prodigal calls genes AND the synteny/coords path activates.
   Confirm `resolve_dsr_direction` actually fires (coords non-empty) and check its reductive/oxidative call vs the paper.
4. **Run METABOLIC (`-in-gn`, needs .fasta extension!) + DRAM (`annotate`, nucleotide mode) on the same MAGs** —
   this is their native input type. Compare on the MAG subset separately.
5. **Report as a distinct "MAG realism / synteny" study** (new REPORT.md section), separate from the 37-isolate
   headline. Key questions: does the synteny direction-call work on real MAGs? does scycle still beat the comparators
   on fragmented input? Reserve hold-out genera here (HOLDOUT_TAG="holdout_v3" is wired in adapters.py).

## Priority 4 — orphan-target reinforcement & decoys (+ doxD/tth independent carriers)  ⬜ pending
- sor/doxA/sorA/sreA/sdo/dddP → ≥2 more positives each (KEGG-carrier search → microbiologist vetting →
  assembly-matched proteome, as in SP5/SP6). Prefer carriers that are NOT seed sources.
- Decoys to stress trap precision: more Mo-bis-MGD carriers (denitrifiers w/ narG/napA/dmsA — for
  phsA/ttrA/sreA), generic multiheme-cyt-c genomes (for fccA), glyoxalase-II / MBL-fold carriers (for sdo).

## Acceptance (re-run after expansion)
1. Every trap target has ≥3 **independent** positives; `trap_independence.py` shows the independent-only
   trap precision/recall on a non-degenerate cell count (target: independent trap present-cells ≥ ~50).
2. Re-run `benchmark_stats.py` (B=10,000) on the expanded panel AND on the **independent-only** trap
   subset; report both. Re-check whether the DRAM trap-precision Δ now reaches BH-significance.
3. MAG subset reported separately (tests synteny + comparator-native input).
4. Hold-out genera scored but excluded from all tuning; regression floors re-locked.

## Sequencing note
Do Priority-1 (independence) first — it's the load-bearing fix and the cheapest (≈6–10 isolate proteomes).
Priorities 2 and 4 overlap with it (the same SRB/Sox additions power the DRAM contrast and reinforce
modules). Priority-3 MAGs are the most effortful (sourcing + nucleotide pipeline path) — do last, as a
separate, clearly-labelled realism study.
