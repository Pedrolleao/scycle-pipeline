# GTDB-scale cross-tool benchmark — 500-genome report (sulfur)

**Status:** complete (499 GTDB r232 representative genomes; 1 enriched genome dropped for an empty
Prodigal proteome). Same locked design as the nitrogen sister tool's run
(`Nitrogen_Cycle/ncycle-pipeline/comparators/gtdb500/GTDB500_REPORT.md`): a **concordance** study —
arbitrary GTDB genomes have no curated ground truth, so we measure where tools **agree/disagree**,
not accuracy. Comparators: **scycle, raw KofamScan, METABOLIC v4.0, SCycDB** (DRAM excluded at scale
by prereg). METABOLIC was added to this GTDB-scale run by **reusing the 380 backbone genomes' calls
from ncycle's GTDB-500 METABOLIC run** (same METABOLIC v4.0; its KO output is metabolism-wide, so the
sulfur KOs are already present) **plus running the 119 S-enriched genomes in ≤12-genome batches** to
stay under METABOLIC's exit-144 scale-failure threshold — giving full 4-vs-4 comparator parity with
the nitrogen sister tool. The accuracy verdict remains the curated-panel benchmark
(`../../validation/benchmark/README.md` + `benchmark_results.tsv`); this run shows the trap
divergence it documents is **pervasive across GTDB**, and adds the **dsrAB reductive-vs-oxidative
direction** split as the sulfur headline.

---

## 1. Panel
499 GTDB r232 species-representatives (`select_gtdb_scyc.py`; `selection.tsv`):
- **119 S-cycle-clade-enriched** across the direction/trap-relevant clades — sulfate reducers
  (*Desulfovibrio, Desulfobacter, Desulfobulbus, Desulfotomaculum, Desulfomicrobium,
  Thermodesulfovibrio, Thermodesulfobacterium*), reverse-Dsr / Sox sulfur oxidizers (*Thiobacillus,
  Chromatiales* PSB, *Beggiatoa/Thiothrix,* SUP05/*Thioglobus, Chlorobi* GSB, *Thioalkalivibrio*),
  DMSP cyclers, and archaea.
- **380 cross-phylum backbone** (random species-reps spanning broad GTDB breadth) — REUSED verbatim
  from the nitrogen tool's gtdb500 backbone, so the two sister tools are benchmarked on the same
  cross-phylum set.
- Fetched via NCBI `datasets`, Prodigal-called uniformly — all four tools see identical proteins.

The dsrAB direction locus is well-populated: scycle present-calls dsrA 85, dsrB 87, dsrD 26,
aprA 90, qmoA 58, soxB 68, sqr 113.

---

## 2. Wall-clock — measured on this 32-core box
- **scycle full pipeline** (incl. hmmscan over the KO HMM DB → matrix): **~30 min** for 499
  proteomes; parallelizes cleanly, *not* the GTDB-scale bottleneck.
- **SCycDB** (DIAMOND blastp vs SCycDB_2020Mar, best-hit): comparable order, single pass.
- **METABOLIC v4.0**: exit-144s at startup on ≥40 GTDB proteomes in a single pass (before building
  `total.faa`), so the **119 S-enriched genomes were run in 10 sequential ≤12-genome batches**
  (~20 min each, ~3.3 h total) — under the failure threshold — and the **380 backbone genomes were
  reused from ncycle's GTDB-500 METABOLIC run** (same v4.0; metabolism-wide KO output, so its sulfur
  KOs were already computed). Hours-scale as on nitrogen, but no longer a blocker.

---

## 3. Concordance findings (full tables in `CONCORDANCE.md` / `CONCORDANCE.tsv`)

### 3.1 Pairwise agreement
Raw present/absent agreement is high but **inflated by shared-absent cells**; the **Jaccard
positive-call agreement** is the honest metric.

| tool pair | Jaccard (ALL) | Jaccard (**trap**) |
|---|---|---|
| scycle vs KofamScan | 0.225 | 0.140 |
| scycle vs METABOLIC | 0.784 | 0.331 |
| scycle vs SCycDB | 0.215 | 0.203 |
| KofamScan vs METABOLIC | 0.238 | 0.327 |
| KofamScan vs SCycDB | 0.401 | 0.351 |
| METABOLIC vs SCycDB | 0.199 | 0.169 |

**Headline:** scycle and **METABOLIC** — the two threshold-gated tools (METABOLIC applies KOfam
adaptive thresholds) — are the most concordant overall (Jaccard 0.784) yet split at the
**homology/direction traps** (0.331), exactly where scycle's BLAST gating and dsrAB-direction
resolution act and METABOLIC's plain KO presence cannot. **Raw KofamScan over-calls broadly** (no
thresholds): it agrees with scycle on only 0.225 of positive calls overall — e.g. dsrA 183/499 vs
scycle's gated 85, tauB 492 vs 18, cysN 483 vs 36. **SCycDB is the other outlier** (Jaccard ~0.20):
its loose DIAMOND best-hit over-calls even harder — cysJ 483/499 (scycle 8), ssuD 461 (6), cysC 498
(53), sat 481 (198), dsrA 332 (85). The two ungated callers (raw KofamScan, SCycDB best-hit) over-
attribute at scale; the threshold-gated tools agree, and scycle alone additionally resolves
direction. This mirrors the nitrogen tool (ncycle most concordant with METABOLIC at 0.617; raw
KofamScan and NCycDB the over-callers).

### 3.2 dsrAB reductive↔oxidative direction split (sulfur headline)
dsrA/dsrB are the **same genes (KO K11180/K11181)** in dissimilatory sulfate reducers (reductive
Dsr, SO₃²⁻→H₂S) and in reverse-Dsr sulfur oxidizers (oxidative rDSR, H₂S/S⁰→SO₃²⁻). Gene presence
alone cannot separate the two metabolic directions; the genomic companions do (dsrD ⇒ reductive;
Sox/sqr with dsrD absent ⇒ oxidative). scycle resolves direction (`apply_rules.resolve_dsr_direction`);
raw KofamScan and SCycDB report dsrA/dsrB presence with **no direction**.

- scycle calls **dsrAB in 89 genomes**, resolving **48 reductive** (sulfate reducers),
  **40 oxidative** (reverse-Dsr sulfur oxidizers), **1 ambiguous**.
- The split tracks textbook biology: the **reductive** calls are *Desulfovibrio, Desulfobacter,
  Desulfobulbus, Desulfotomaculum, Desulfomicrobium, Thermodesulfovibrio, Thermodesulfobacterium*;
  the **oxidative** calls are *Thiobacillus, Chromatiales* PSB, *Beggiatoa/Thiothrix,* SUP05/
  *Thioglobus, Chlorobi* GSB, *Thioalkalivibrio* — the canonical sulfur-oxidizer clades.
- **All 40** oxidative (reverse-Dsr) genomes have ≥1 direction-blind comparator (raw KofamScan,
  METABOLIC, or SCycDB) reporting dsrA/dsrB present — i.e. **every one would be mis-attributed as a
  sulfate reducer** by a tool that reads dsrAB presence as dissimilatory sulfate reduction.
  METABOLIC, despite matching scycle's gated dsrA/dsrB *presence* (85/87), has no reductive-vs-
  oxidative KO and so falls into the trap exactly like the others. This is the dsrAB-direction
  analogue of ncycle's nxrA→narG mis-routing, now observed at GTDB scale.

### 3.3 scycle-distinctive loci (differs from every representing comparator)
Largest: **sdo 250, phsA 46, fccA 21, dsrC 20, dsrD 7** — dominated by the BLAST-gated / KO-less
trap targets the KO-only and best-hit tools cannot separate. Concretely:
- **sdo** (sulfur dioxygenase, K17725): scycle is distinctive here by **withholding** an over-call —
  raw KofamScan fires the KO in 444/499 genomes and METABOLIC in 254, but scycle's BLAST gate keeps
  only 4, so scycle differs from both KO tools in ~250 genomes (SCycDB has no sdo family → gap).
- **phsA** (thiosulfate reductase, shared K08352 with polysulfide reductase psrA): raw KofamScan
  over-fires the shared KO to 346 genomes and METABOLIC to 57; scycle's BLAST gate keeps only 4.
  SCycDB resolves the split the other way (137, classifying many as the unmapped psrA → SCycDB FN).
- **fccA** (sulfide oxidation, K17230 cyt-c): scycle gates the over-calls down (17 vs KofamScan 121).
- **dsrC / dsrD / soxD / otr**: Pfam-only / KO-less targets KofamScan and METABOLIC cannot represent
  via KO (shown as `-` / 0); scycle calls them via Pfam+BLAST.

### 3.4 Comparator coverage gaps (documented, not divergence)
KofamScan/METABOLIC cannot represent the KO-less targets dsrC, dsrD, soxD, otr (no KO → `-`, not a
disagreement). SCycDB has no family for **sdo** (sulfur dioxygenase) or **tth** (thiosulfate
dehydrogenase) → those are SCycDB coverage gaps (predicted absent = FN), and SCycDB splits
phsA/psrA on the shared K08352 (a scycle phsA may classify as the unmapped psrA → SCycDB FN).

---

## 4. Honest caveats
- **Concordance ≠ accuracy.** Without ground truth this measures *divergence*, not who is correct.
  For the trap loci the curated-panel benchmark (bootstrap-CI, CI-lower-bound gates) already
  established scycle as correct; this run shows the divergence is pervasive at scale and that the
  dsrAB direction is recovered along expected clade lines.
- The dsrAB direction call here is **presence-based** (dsrD + qmo/sox companions), the mode that
  applies to pre-called proteomes; on nucleotide/MAG input scycle additionally uses dsrD↔dsrAB
  operon synteny. The 1 "ambiguous" genome carries dsrAB without a resolved companion signal.
- For widespread biosynthetic targets (cysK/cysC/cysH/sat) the higher SCycDB counts may be partly
  biologically real rather than pure over-calls; we do not adjudicate those here.
- Enrichment means the panel is **not** a random GTDB draw (by design, to populate the dsrAB and
  trap loci); the counts are per-genome call concordance, not population frequencies.

---

## 5. Reproduction
```bash
# selection → selection.tsv / enriched_accs.txt (380 backbone reused from nitrogen gtdb500)
python3 comparators/gtdb500_s/select_gtdb_scyc.py
# datasets download → gtdb_dl/ ; prodigal → proteomes/
#   scycle + kofam: scycle.py / snakemake over the 499 proteomes (--keep-going) → scycle_results/
python3 comparators/build_scycdb_tsv.py --proteome-dir comparators/gtdb500_s/proteomes \
        --out comparators/gtdb500_s/scycdb.tsv --outdir comparators/gtdb500_s/scycdb_out
# METABOLIC: 119 S-enriched genomes in 10 sequential ≤12-genome batches (exit-144 above ~40 in one
#   pass); 380 backbone reused from ../Nitrogen_Cycle/ncycle-pipeline/comparators/gtdb500/metabolic.tsv
bash ../comparators/metabolic_gtdb500/run_batches.sh        # → metabolic_gtdb500/kegg_all/*.result.txt
bash ../comparators/metabolic_gtdb500/merge_metabolic.sh    # 380 reused + 119 fresh → metabolic.tsv (499)
python3 validation/benchmark/concordance.py --results comparators/gtdb500_s/results \
        --scycdb comparators/gtdb500_s/scycdb.tsv \
        --metabolic comparators/gtdb500_s/metabolic.tsv \
        --panel-map comparators/gtdb500_s/selection.tsv \
        --out comparators/gtdb500_s/CONCORDANCE.md
```

---

## 6. Outcome
The harness ran end-to-end at 499 genomes and the concordance signal is clear, on-message, and a
faithful sulfur counterpart to the nitrogen gtdb500 run: **scycle is most concordant with the
threshold-gated METABOLIC (Jaccard 0.784 overall) and diverges at the traps (0.331); the ungated
callers — raw KofamScan and SCycDB best-hit — over-call broadly (Jaccard ~0.20–0.23); and the dsrAB
reductive-vs-oxidative direction — invisible to all three comparators — is resolved by scycle into
48 sulfate reducers vs 40 reverse-Dsr sulfur oxidizers along textbook clade lines, with all 40
oxidizers mis-attributable as reducers by a direction-blind tool.** With METABOLIC added via the
reuse-plus-batched run, this run now has **full 4-vs-4 comparator parity** with the nitrogen sister
tool. No further work is required to close this benchmark line.
