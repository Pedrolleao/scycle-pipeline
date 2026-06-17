# Benchmark pre-registration — scycle-pipeline vs external S-cycle tools

Design adapted from the validated ncycle-pipeline benchmark pre-registration. The
endpoints, trap-target set, vocabulary, and statistics were fixed (in `adapters.py` /
`benchmark_stats.py`) **before** running the comparison; the ground-truth corrections were
made independently of either tool's output (see `build_curated_function_gt.py`). No tuning
of either side after seeing comparison results; no metric/endpoint shopping.

**CHANGELOG / honesty note (2026-05-27 audit).** This prereg was first written against a
**22-genome** panel; the endpoints/trap-set/vocabulary/statistics below were frozen then. The
panel later grew to **25 genomes**: the SP5 additions (*Acidianus ambivalens*, *Starkeya
novella*) and the SP6b addition (*Roseobacter denitrificans*) were made to provide positives for
otherwise-untestable orphan targets (sor/doxA/sorA/sreA/sdo and **dddP — rde is the only dddP
positive on the panel**). This is outcome-directed sample augmentation and is disclosed as such:
the orphan-target results (sor, doxA, sorA, sreA, sdo, dddP) should be read as "single-positive,
added-to-validate", NOT as part of the frozen-design headline. The endpoints and trap set were
NOT changed when the panel grew.

**Panel reached n=43 at SP8 (2026-05-27); reconciled here 2026-06-14.** After the n=25 note above, the
**SP8 expansion** (documented in `build_ground_truth.py` KEGG_ORG map + `EXPANSION_PLAN.md`) added
~18 genomes to reach **43**, specifically **INDEPENDENT positives** — genomes that are NOT a curated
seed / custom-HMM source for ANY target — added to de-circularize the homology-trap claim (e.g.
*Thioalkalivibrio sulfidiphilus* and *Thiocystis violascens* as independent oxidative/reverse-dsrAB
positives; *Proteus mirabilis* for independent phsA+ttrA; independent doxD/doxA). **20 of the 43 are
INDEPENDENT** (the `INDEPENDENT` set in `build_ground_truth.py`), which is exactly what powers the
per-(genome,target) trap-independence metric in `../trap_independence.py`. The frozen confirmatory
DESIGN (endpoints, trap-target set, vocabulary, statistics) was **not** changed by SP8 — only positives
were added, which strengthens (not relaxes) the trap claim's independence. **All "n=22"/"n=25" figures
below are historical;** the final B8 benchmark run (2026-06-14) scored **n=43** (`benchmark_results.tsv`,
`COMPARISON_REPORT.md`). Trap-claim independence is audited in `../trap_independence.py`.

## Hypotheses
- **Primary (confirmatory):** on the homology-trap targets, scycle-pipeline has higher
  **precision** than each comparator. Expected large effect (the pipeline gates the
  shared/generic-KO failure mode) → the panel is well powered here.
- **Secondary:** scycle-pipeline has higher **overall micro-F1** than each comparator.
  Small expected effect → modestly powered at n=25; reported, not relied on.
- **Exploratory:** per-step F1 (report CIs; not part of the corrected family).

## Comparators (default settings, no asymmetric tuning)
| tool | approach | resolves subunits? |
|---|---|---|
| **scycle-pipeline** (this tool) | KO + custom HMM + BLAST gating (+ synteny on nt input) | yes (61 targets) |
| **raw KofamScan** | KO-only, stock `ko_list` thresholds, shared KOs → all targets | yes (data already available) |
| **METABOLIC** | HMM-based N/S/C cycle step detection | partial |
| **DRAM** | MAG metabolism distillation (KEGG/Pfam) | coarse (steps/modules) |
| **SCycDB** | curated S-cycle gene-family DB + DIAMOND best-hit | yes (gene families) — *added 2026-06-12, see amendment* |

Each tool's exact version + command line is recorded in `adapters.py` when run. Default thresholds only.

## Ground truth — curated FUNCTION (not KEGG-KO)
The KEGG-KO GT (`ground_truth.tsv`) makes a raw-KofamScan comparison near-circular. The
benchmark scores against `curated_function_gt.tsv`: the KEGG-KO GT with biology-verified
function-level corrections on the homology-trap cells where KO-presence ≠ functional role.
Each correction is justified by a verified protein function (UniProt), subunit/operon context,
or characterized biology of the reference strain — derived **independently of either tool**.
The KEGG-KO GT is reported only as a contrast to make the circularity visible.

## Inputs
Same 25-genome panel proteomes for all tools. (No hold-out genera are reserved in SP1 yet;
`HOLDOUT_TAG` is wired for when they are.) scycle's operon-synteny step is dormant on proteome
input — this benchmarks the KO+HMM+gating core, the fair common ground.

## Comparison vocabulary (two pre-specified resolutions)
1. **Subunit** — the 61 scycle targets. Coarse tools are mapped *down* (a step call expands to
   all member subunits; credit/blame symmetric). The trap claim lives here. **Trap set** =
   the targets where the pipeline's machinery works beyond raw KO: gated traps
   (dsrA/dsrB/soxD/fccA/sdo/phsA/ttrA/sreA/dmdA), KO-less Pfam+synteny / BLAST calls
   (dsrC/dsrD/soxD/otr), and the custom-HMM / KEGG-sparse marker doxD.
2. **Step** — 12 S-cycle transformation steps (`STEP_MEMBERS`/`STEP_DIAG` in `adapters.py`).
   scycle is mapped *up* (step present iff its diagnostic marker present). Fairest to coarse tools.

Ground truth is collapsed to each resolution by the same rules for all tools.

## Scoring & statistics
- **Scorer:** a cell counts only if it is in the ground truth. Positive class = "present".
- **Unit of independence = genome.** All CIs use a **genome-level cluster bootstrap**
  (resample the 25 genomes with replacement, B=10,000).
- **95% CIs** (percentile) on precision / recall / F1 per tool, per subset, per resolution.
- **Primary significance:** paired **bootstrap of the difference** Δ = metric(scycle) −
  metric(comparator) on the same resampled genomes; win iff the 95% CI of Δ excludes 0.
- **Secondary significance:** **exact McNemar** on cell-level discordant pairs (ignores clustering).
- **Multiple testing:** **Benjamini–Hochberg FDR** across {comparators} × {trap-precision
  (subunit), overall-F1}. Exploratory metrics not corrected.
- **Effect sizes** (ΔF1, Δprecision with CIs) reported alongside every p-value.

## Decision rules
- A comparator is "beaten on the trap" iff Δ trap-precision CI excludes 0 (BH q<0.05) AND
  the point estimate favours scycle.
- Overall-F1 differences are reported with CIs; if a CI includes 0 we state the panel is
  underpowered for that contrast rather than claiming parity or a win.

## Known limitations
- n=25 genomes: well powered for the trap (large effect), modestly powered for the small
  overall gap. Panel expansion is the lever to power the secondary endpoint.
- Trap **recall** is not a confirmatory endpoint: the pipeline's gates are tuned for precision
  and still over-gate a few divergent traps (fccA divergent cyt-c subunit; lone-subunit phsA).
  An fccA←fccB co-occurrence rule (analogous to dsrC←dsrA) would recover most of it; deferred.
- Vocabulary mapping is a bias source; it is fixed here, applied identically to all tools.
- KEGG-based tools may have seen some panel genomes in training references — a shared
  limitation that does not affect the trap-discrimination claim.

## Amendments (2026-06-07) — reporting harmonization with ncycle-pipeline
Recorded as a transparent post-hoc amendment; the **confirmatory family is unchanged** (trap
precision + overall-F1, BH-corrected) and the whole-panel headline is unchanged. Descriptive only.
- **Three reporting frames are now reported consistently across both pipelines:** (A) whole-panel
  micro-F1 (the headline, unchanged); (B) per-(genome,target) independent-subset trap precision
  (`validation/trap_independence.py`); (C) an explicit train/hold-out split + leave-one-clade-out
  CV (`validation/logo_cv.py`). scycle already reported A + B + LOGO-CV; Frame C is added for
  parity — the 20 truly-independent panel genomes (no seed/HMM source for any target; the
  `INDEPENDENT` set in `build_ground_truth.py`) are scored as a held-out set in
  `score_scycle.py`. This is a descriptive cut: it does NOT retag the GT file, change the
  whole-panel headline, or alter the pre-registered benchmark (which stays whole-panel).

## Amendment (2026-06-12) — SCycDB added (domain-DB comparator; mirrors ncycle's NCycDB)
SCycDB (Yu et al. 2020, github.com/qichao1984/SCycDB) was installed and run on the panel
proteomes; provenance + the fixed family→target map are in `adapters.py` (`load_scycdb`,
`SCYC_MAP`). The search is faithful to `SCycDB_FunctionProfiler.PL` protein default
(`diamond blastp -k 1 -e 1e-4` vs `SCycDB_2020Mar`). SCycDB adds its two pre-specified rows
({trap-precision, ALL-F1}) to the same BH-corrected family, exactly as NCycDB does for ncycle;
**no endpoints, decision rules, genome set, or vocabulary rules changed.** Comparator set is now
**{raw KofamScan, METABOLIC v4.0, DRAM v1.4.6, SCycDB}** — the same four-comparator shape as
ncycle (which uses NCycDB as its domain DB), so the two sister benchmarks are now structurally
identical. scycle beats SCycDB on both endpoints (ALL-F1 Δ≈0.37, trap-precision Δ≈0.60; q<0.001).

## Amendment (2026-06-14, B8) — battery harmonization + raw-KofamScan parser fix
Two harmonization actions across both sister benchmarks; **the confirmatory family, decision
rules, genome set, GT, and vocabulary are unchanged** — these are reproducibility/methods fixes.
- **Comparator-set parity confirmed:** scycle = {KofamScan, METABOLIC, DRAM, SCycDB} and ncycle =
  {KofamScan, METABOLIC, DRAM, NCycDB} — identical shape, domain-DB included on both. The earlier
  ncycle note about a divergent comparator set is now resolved.
- **raw-KofamScan baseline parser fix:** the pipeline's HMMER output is now `hmmsearch --domtblout`
  (KO in column 4, full-sequence score in column 8); the benchmark's `load_kofam` had keyed on the
  legacy `hmmscan --tblout` column 1. The adapter now detects either layout. This restores the raw
  KofamScan baseline (it was silently scoring all-absent after the pipeline's hmmsearch switch); it
  does **not** touch scycle's own calls (read from the matrix) and does not alter any endpoint.
- **Final re-run:** B=10,000 genome cluster bootstrap, seed 1234, GT=`curated_function_gt.tsv`,
  scored n=43 panel genomes → `benchmark_results.tsv` + `COMPARISON_REPORT.md`.

## Amendment (2026-06-16) — secondary TN-inclusive metrics (specificity + MCC)
Reported as a transparent **post-hoc secondary** analysis; the pre-registered confirmatory family
(trap-precision + ALL-F1, BH-FDR) is **unchanged** and its numbers are bit-identical on re-run.
Added two metrics that incorporate true negatives (correct absences), which precision/recall/F1 do
not: **specificity** = TN/(TN+FP), and **MCC** (Matthews correlation coefficient) over all four
confusion cells — robust to the present:absent class imbalance. Computed on the same scored GT cells
with the same B=10,000 cluster bootstrap; the MCC contrasts (ALL + trap × 4 comparators) get their
**own** BH-FDR correction, kept separate from the primary family. All 8 favour scycle (q<0.0001).
`metric_of` in `benchmark_stats.py` gains `specificity`/`mcc`; engine TN was already tallied.
