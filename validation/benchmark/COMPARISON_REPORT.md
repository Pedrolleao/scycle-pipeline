# scycle-pipeline — Cross-Tool Comparison Report

Pre-registered benchmark (`prereg.md`) of scycle-pipeline against external sulfur-cycle annotation
tools, scored against the independent **curated-function** ground truth (`curated_function_gt.tsv`),
n=43 panel genomes. Genome-level cluster bootstrap B=10,000, seed 1234. Sister report to
ncycle-pipeline's `COMPARISON_REPORT.md` — identical battery, structure, and statistics.

## 1. Headline

| subunit metric | scycle | raw KofamScan | METABOLIC | DRAM | SCycDB |
|---|---|---|---|---|---|
| ALL micro-F1 | **0.930** | 0.550 | 0.651 | 0.645 | 0.563 |
| ALL precision | **0.952** | 0.399 | 0.904 | 0.719 | 0.404 |
| ALL recall | 0.909 | 0.884 | 0.508 | 0.585 | **0.933** |
| **homology-trap precision** | **0.973** | 0.211 | 0.618 | 0.794 | 0.378 |
| homology-trap recall | 0.901 | 0.617 | 0.420 | 0.617 | 0.914 |
| non-trap micro-F1 | **0.929** | 0.595 | 0.676 | 0.638 | 0.568 |
| ALL specificity (TN rate) | **0.986** | 0.592 | 0.984 | 0.930 | 0.578 |
| **ALL MCC** (TP/FP/FN/TN) | **0.910** | 0.404 | 0.616 | 0.556 | 0.435 |

The two TN-inclusive rows (added 2026-06-16, §3.1) make the false-positive control explicit:
**specificity** (true-negative rate) exposes the over-callers — raw KofamScan 0.592 and SCycDB 0.578
mis-call ~40% of truly-absent proteins as present, vs scycle 0.986 — and **MCC** (all four cells, robust
to the 614:2009 present:absent imbalance) summarises overall behaviour at scycle 0.910 vs 0.40–0.62.
METABOLIC's specificity (0.984) is near scycle's: it is a *conservative under-caller*, not an
over-caller, but its low recall still caps MCC at 0.616.

scycle leads on overall F1 and decisively on **homology-trap precision** (0.973 vs 0.21–0.79). Note
scycle is **not** 1.000 here (unlike ncycle on its curated GT): the sulfur curated-function GT is an
independent, harder reference, and scycle's precision-tuned gates still over-gate a few divergent
traps (e.g. fccA, lone-subunit phsA), costing some trap recall — disclosed in the prereg.

| tool | trap failure signature | mechanism |
|---|---|---|
| raw KofamScan | **massively over-calls** (P 0.21, R 0.62) | shared/generic KOs (K17230 cyt-c → fccA; K08352 phsA/psrA → phsA) map to *all* their targets |
| METABOLIC | **under-calls** (P 0.62, **R 0.42**) | stricter bundled-KOfam m-cutoff + KO-less Dsr subunits it cannot carry |
| DRAM | **under-covers** (R 0.62; only 3 sulfur modules distilled) | distillate resolves dsr / asr / Sox modules only → many traps unrepresented |
| SCycDB | **over-calls / mis-routes** (P 0.38, R 0.91) | DIAMOND best-hit over-attributes broadly (the domain-DB best-hit failure mode) |

## 2. Tools compared

| tool | approach | resolution | status |
|---|---|---|---|
| **scycle-pipeline** (this tool) | KOfam HMM + clade-specific custom HMMs + curated BLAST gates (+ operon synteny on nucleotide input) | subunit (61 targets) | — |
| **raw KofamScan** | KO-only, stock `ko_list` thresholds; shared KOs map to *all* their targets; no gating | subunit (by KO) | ✅ |
| **METABOLIC** v4.0 (Zhou et al. 2022) | HMM-based function profiling (`METABOLIC-G.pl`, `-kofam-db full`) → per-genome KO presence | subunit (by KO) | ✅ |
| **DRAM** v1.4.6 (Shaffer et al. 2020) | MAG metabolism distillation (KEGG modules) | coarse (modules → steps) | ✅ |
| **SCycDB** (Yu et al. 2020) | curated S-cycle gene-family DB (207 families) + DIAMOND best-hit | subunit (gene families) | ✅ |

## 3. Pre-registered confirmatory result (subunit, BH-FDR corrected)

All 8 contrasts (4 comparators × {ALL-F1, trap-precision}) favour scycle after Benjamini–Hochberg FDR.

| comparator | metric | Δ (scycle − cmp) | 95% CI | bootstrap p | verdict |
|---|---|---|---|---|---|
| raw KofamScan | ALL F1 | +0.380 | [0.346, 0.420] | <0.0001 | **scycle WINS** |
| METABOLIC | ALL F1 | +0.279 | [0.248, 0.309] | <0.0001 | **scycle WINS** |
| DRAM | ALL F1 | +0.285 | [0.230, 0.342] | <0.0001 | **scycle WINS** |
| SCycDB | ALL F1 | +0.367 | [0.326, 0.411] | <0.0001 | **scycle WINS** |
| raw KofamScan | trap precision | +0.762 | [0.693, 0.833] | <0.0001 | **scycle WINS** |
| METABOLIC | trap precision | +0.355 | [0.241, 0.495] | <0.0001 | **scycle WINS** |
| DRAM | trap precision | +0.180 | [0.061, 0.323] | 0.0012 | **scycle WINS** |
| SCycDB | trap precision | +0.596 | [0.522, 0.680] | <0.0001 | **scycle WINS** |

### 3.1 Secondary (post-hoc) MCC family — BH-FDR corrected

TN-inclusive confirmation, scored on the same cells with the same B=10,000 cluster bootstrap, BH-FDR
over its own 8-contrast family (kept **separate** from the frozen pre-registered family above so the
confirmatory result is unchanged). All 8 favour scycle.

| comparator | metric | Δ (scycle − cmp) | 95% CI | BH q | verdict |
|---|---|---|---|---|---|
| raw KofamScan | ALL MCC | +0.506 | [0.470, 0.541] | <0.0001 | **scycle WINS** |
| METABOLIC | ALL MCC | +0.294 | [0.259, 0.325] | <0.0001 | **scycle WINS** |
| DRAM | ALL MCC | +0.354 | [0.284, 0.424] | <0.0001 | **scycle WINS** |
| SCycDB | ALL MCC | +0.475 | [0.433, 0.520] | <0.0001 | **scycle WINS** |
| raw KofamScan | trap MCC | +0.765 | [0.690, 0.845] | <0.0001 | **scycle WINS** |
| METABOLIC | trap MCC | +0.482 | [0.393, 0.598] | <0.0001 | **scycle WINS** |
| DRAM | trap MCC | +0.270 | [0.162, 0.408] | <0.0001 | **scycle WINS** |
| SCycDB | trap MCC | +0.441 | [0.369, 0.527] | <0.0001 | **scycle WINS** |

## 4. Per-step micro-F1 (step resolution, fairest to coarse tools)

| step | scycle | KofamScan | METABOLIC | DRAM | SCycDB |
|---|---|---|---|---|---|
| ALL | **0.942** | 0.547 | 0.774 | 0.529 | 0.504 |

Even at step resolution (where coarse tools are mapped up and given maximum credit), scycle leads.

---

**Provenance.** Engine `benchmark_stats.py`; loaders + tool versions/commands in `adapters.py`;
data `benchmark_results.tsv`. Re-run: `make benchmark` (or
`GT_FILE=curated_function_gt.tsv python benchmark_stats.py --bootstrap 10000 --seed 1234 --metabolic
../../../comparators/metabolic.tsv --dram ../../../comparators/dram.tsv --scycdb scycdb.tsv`).
This report reflects the 2026-06-14 (B8) harmonized re-run; see `prereg.md` amendments.
