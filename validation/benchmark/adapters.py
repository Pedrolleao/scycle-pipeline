#!/usr/bin/env python3
"""
adapters.py — load each tool's calls into ONE shared vocabulary, plus ground truth.

Every loader returns subunit-level predictions: {(genome, target_id): bool_present},
restricted to the training panel (hold-outs excluded). The benchmark engine collapses to
the `step` resolution itself (STEP_DIAG below), so the collapse rule is applied identically
to every tool and to the ground truth.

Working now: load_truth, load_scycle (our pipeline), load_kofam (raw KO baseline).
Scaffolded: load_metabolic / load_dram — parse a normalized TSV from that tool's native
output and map identifiers to our targets; return {} when given no input path.

GROUND TRUTH: defaults to the curated-FUNCTION GT (`curated_function_gt.tsv`), the
non-circular reference for the "beats raw KofamScan?" claim. Set env GT_FILE=ground_truth.tsv
to score against the KEGG-KO GT instead (shows the circularity the curated GT removes).
"""
from __future__ import annotations

import csv
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GT = ROOT / "validation" / os.environ.get("GT_FILE", "curated_function_gt.tsv")
MATRIX = ROOT / "scycle_results" / "scycle_matrix.tsv"
TARGETS = ROOT / "config" / "targets.yaml"
KO_LIST = ROOT / "resources" / ".cache" / "ko_list"
RESULTS = ROOT / "scycle_results"

HOLDOUT_TAG = "holdout_v3"   # none in the SP1 panel yet; reserved

# Homology-trap targets — the confirmatory-claim subset (subunit resolution). These are the
# targets where the pipeline's machinery does work beyond raw KO: requires_blast_for_confirmation
# gates (dsrA/dsrB/soxD/fccA/sdo/phsA/ttrA/sreA/dmdA), KO-less Pfam+synteny / BLAST calls
# (dsrC/dsrD/soxD/otr), and the custom-HMM / KEGG-sparse marker doxD. Raw KO either over-fires a
# shared/generic KO here (K17230 cyt-c → fccA; K08352 phsA/psrA → phsA) or cannot call at all.
TRAP = {"dsrA", "dsrB", "dsrC", "dsrD", "soxD", "fccA", "sdo",
        "phsA", "ttrA", "sreA", "dmdA", "otr", "doxD"}

# ── step resolution ──────────────────────────────────────────────────────────
# STEP_MEMBERS: every subunit belonging to a step (maps a coarse tool's step call DOWN to
# subunits). STEP_DIAG: the diagnostic marker(s) whose presence defines the step (collapses
# subunit calls + truth UP to steps).
STEP_MEMBERS: dict[str, list[str]] = {
    "dissim_sulfate_reduction": ["sat", "aprA", "aprB", "qmoA", "qmoB", "qmoC",
                                  "dsrA", "dsrB", "dsrC", "dsrD", "dsrM", "dsrK"],
    "sox_thiosulfate_oxidation": ["soxA", "soxB", "soxC", "soxD", "soxX", "soxY", "soxZ"],
    "sulfide_oxidation":         ["sqr", "fccA", "fccB", "sdo"],
    "sulfite_oxidation":         ["soeA", "soeB", "soeC", "sorA"],
    "assim_sulfate_reduction":   ["cysN", "cysD", "cysNC", "cysC", "cysH",
                                  "cysJ", "cysI", "sir", "cysK", "cysM"],
    "thiosulfate_reduction":     ["phsA", "phsB", "phsC"],
    "tetrathionate_reduction":   ["ttrA", "ttrB", "ttrC"],
    "tetrathionate_oxidation":   ["tsdA", "doxA", "doxD", "tth"],
    "polysulfide_transformation": ["sreA", "sseA", "otr"],
    "dmsp_demethylation":        ["dmdA"],
    "dmsp_cleavage":             ["dddP", "mddA", "mtoX"],
    "organosulfonate":           ["tauD", "tauA", "tauB", "tauC", "ssuD", "ssuE"],
}
STEP_DIAG: dict[str, list[str]] = {
    "dissim_sulfate_reduction":  ["dsrA", "aprA"],
    "sox_thiosulfate_oxidation": ["soxB", "soxC"],
    "sulfide_oxidation":         ["sqr", "fccB"],
    "sulfite_oxidation":         ["soeA", "sorA"],
    "assim_sulfate_reduction":   ["cysH", "cysI"],
    "thiosulfate_reduction":     ["phsA"],
    "tetrathionate_reduction":   ["ttrA"],
    "tetrathionate_oxidation":   ["tsdA", "doxD", "tth"],
    "polysulfide_transformation": ["sreA", "otr"],
    "dmsp_demethylation":        ["dmdA"],
    "dmsp_cleavage":             ["dddP", "mddA"],
    "organosulfonate":           ["ssuD", "tauD"],
}
TARGET_STEP = {t: step for step, members in STEP_MEMBERS.items() for t in members}


# ── shared helpers ───────────────────────────────────────────────────────────

def _cfg_targets() -> list[dict]:
    import yaml
    return yaml.safe_load(open(TARGETS))["targets"]


def target_kos() -> dict[str, list[str]]:
    return {t["id"]: (t.get("ko") or []) for t in _cfg_targets()}


def all_target_ids() -> list[str]:
    return [t["id"] for t in _cfg_targets()]


def load_truth() -> tuple[dict[tuple[str, str], str], set[str]]:
    """Return ({(genome,target): 'present'|'absent'}, training_genomes), hold-outs dropped."""
    truth, genomes = {}, set()
    with open(GT) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("source") == HOLDOUT_TAG:
                continue
            truth[(r["genome"], r["target"])] = r["expected"]
            genomes.add(r["genome"])
    return truth, genomes


# ── tool loaders (subunit-level predictions) ─────────────────────────────────

def load_scycle(genomes: set[str]) -> dict[tuple[str, str], bool]:
    """Our pipeline: matrix status codes 1/2 = present."""
    pred: dict[tuple[str, str], bool] = {}
    with open(MATRIX) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["sample"]
            if g not in genomes:
                continue
            for col, val in r.items():
                if col.startswith("target__"):
                    try:
                        pred[(g, col[len("target__"):])] = int(val) in (1, 2)
                    except (TypeError, ValueError):
                        pred[(g, col[len("target__"):])] = False
    return pred


def load_kofam(genomes: set[str]) -> dict[tuple[str, str], bool]:
    """Raw KofamScan baseline: target present iff ANY of its KOs clears the stock
    ko_list threshold (no gating, no custom HMMs; shared KOs → all their targets)."""
    ko_thr: dict[str, float] = {}
    with open(KO_LIST) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[1] not in ("", "-"):
                try:
                    ko_thr[p[0]] = float(p[1])
                except ValueError:
                    pass
    tkos = target_kos()
    pred: dict[tuple[str, str], bool] = {}
    for g in genomes:
        tbl = RESULTS / g / "hmm" / f"{g}.hmmscan.tsv"
        present_kos = set()
        if tbl.exists():
            with open(tbl) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    c = line.split()
                    if len(c) < 8:
                        continue
                    # KO name + full-sequence score column depend on the HMMER layout:
                    # hmmsearch --domtblout (current pipeline) → query=col3, score=col7;
                    # legacy hmmscan --tblout → target=col0, score=col5. Detect by which
                    # column holds a known KO so the raw-KO baseline survives either format.
                    if c[3] in ko_thr:
                        ko, score = c[3], c[7]
                    elif c[0] in ko_thr:
                        ko, score = c[0], c[5]
                    else:
                        continue
                    try:
                        if float(score) >= ko_thr[ko]:
                            present_kos.add(ko)
                    except ValueError:
                        pass
        for tid, kos in tkos.items():
            pred[(g, tid)] = any(k in present_kos for k in kos)
    return pred


# ── external comparator adapters ──
# RUN PROVENANCE (2026-05-27, all on the same 25 panel proteomes, default settings):
#   METABOLIC v4.0 (github AnantharamanLab@9723633) — `METABOLIC-G.pl -in <proteomes> -t 12
#     -kofam-db full`. metabolic.tsv = worksheet1 (HMM-based function presence): for each genome
#     column == "Present", emit (genome, that row's 'Corresponding KO'). KO-less METABOLIC HMMs
#     (e.g. dsrD, dsrMKJOP) cannot be carried by this genome,ko adapter (same limit as raw kofam).
#   DRAM v1.4.6 (env DRAM14; KOfam + distillation forms, Pfam/UniRef/dbCAN intentionally omitted —
#     irrelevant to sulfur) — `DRAM.py annotate_genes -i '<proteomes>/*.faa'` then `DRAM.py distill`.
#     dram.tsv = DRAM's distillate 'Sulfur' modules; a module is present in a genome iff DRAM's
#     metabolism_summary detected that module's diagnostic gene (STEP_DIAG). DRAM distills only the
#     3 modules in DRAM_STEP_MAP below.

def _read_norm_tsv(path: Path | None) -> list[dict]:
    if not path or not Path(path).exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_metabolic(path: Path | None, genomes: set[str]) -> dict[tuple[str, str], bool]:
    """METABOLIC (Zhou et al. 2022) reports KO presence per genome; export a normalized TSV
    with columns genome, ko. A KO present → all its scycle targets present (same KO→target map
    as KofamScan). Returns {} if no path."""
    rows = _read_norm_tsv(path)
    if not rows:
        return {}
    tkos = target_kos()
    present_kos: dict[str, set[str]] = {}
    for r in rows:
        present_kos.setdefault(r["genome"], set()).add(r["ko"])
    return {(g, tid): any(k in present_kos.get(g, set()) for k in kos)
            for g in genomes for tid, kos in tkos.items()}


# DRAM distillate function name → scycle STEP. The engine expands a present step to its
# STEP_MEMBERS at subunit resolution.
# Aligned to DRAM 1.4.6's ACTUAL `genome_summary_form` 'Sulfur'-header module labels (verified
# 2026-05-27 against DRAM_data/forms/genome_summary_form). DRAM's default distillate resolves ONLY
# these three sulfur modules. It has NO distillate module for sulfide oxidation (sqr/fcc/sdo),
# sulfite oxidation (soe/sor), tetrathionate, polysulfide, DMSP, or organosulfonate — so DRAM emits
# nothing for those steps (a genuine DRAM-coverage limitation, scored as absent, not a mapping gap).
DRAM_STEP_MAP: dict[str, str] = {
    "dissimilatory sulfate reduction, sulfate => h2s": "dissim_sulfate_reduction",
    "assimilatory sulfate reduction, sulfate => h2s": "assim_sulfate_reduction",
    "thiosulfate oxidation by sox complex, thiosulfate => sulfate": "sox_thiosulfate_oxidation",
}


def load_dram(path: Path | None, genomes: set[str]) -> dict[tuple[str, str], bool]:
    """DRAM (Shaffer et al. 2020): export a normalized TSV with columns genome, function,
    present. A present step is expanded to all its STEP_MEMBERS subunits. Returns {} if no path."""
    rows = _read_norm_tsv(path)
    if not rows:
        return {}
    step_present: dict[str, set[str]] = {}
    for r in rows:
        step = DRAM_STEP_MAP.get(r["function"].strip().lower())
        if step and str(r.get("present", "")).lower() in ("1", "true", "yes", "present"):
            step_present.setdefault(r["genome"], set()).add(step)
    pred: dict[tuple[str, str], bool] = {}
    for g in genomes:
        present_subunits = {t for s in step_present.get(g, set()) for t in STEP_MEMBERS[s]}
        for tid in all_target_ids():
            pred[(g, tid)] = tid in present_subunits
    return pred


# ── SCycDB domain-database comparator (mirror of ncycle's NCycDB) ──
# RUN PROVENANCE (2026-06-12, all on the same 43 panel proteomes, default settings):
#   SCycDB (Yu et al. 2020, github.com/qichao1984/SCycDB) — SCycDB_2020Mar.faa (911,805 seqs;
#   579,056 family-mapped) + id2gene.2020Mar.map (207 gene families). DIAMOND, faithful to
#   SCycDB_FunctionProfiler.PL protein default: `diamond makedb --in SCycDB_2020Mar.faa` then per
#   proteome `diamond blastp -k 1 -e 1e-4`. Best-hit subject → family via id2gene.2020Mar.map;
#   scycdb.tsv = (sample, family, count). Built by ../../comparators/build_scycdb_tsv.py
#   (DB + per-genome hits under comparators/SCyc/ + comparators/scycdb_out/).

# SCycDB gene-family name → scycle target id. COMPLETE over the 59 of 61 scycle targets that
# SCycDB covers; families with no scycle target are intentionally omitted (SCycDB predicts them,
# but they are not scycle targets — e.g. asrABC, dmsABC, psrABC, hdr*, met*, sgp*, tus*). Two
# scycle targets are SCycDB COVERAGE GAPS (no SCycDB family → SCycDB predicts absent = FN):
#   sdo (sulfur dioxygenase, K17725) and tth (thiosulfate dehydrogenase, K27925).
# NB: SCycDB resolves FINER than scycle in places — it splits phsA (thiosulfate reductase) from
# psrA (polysulfide reductase), which share K08352; a protein scycle calls phsA may be classified
# by SCycDB as psrA (unmapped) → that surfaces as a SCycDB FN, the honest reflection of SCycDB's
# call. Vocabulary is a fixed, reviewable table applied identically to all genomes (prereg).
SCYC_MAP: dict[str, str] = {
    # dissimilatory sulfate reduction
    "sat": "sat", "aprA": "aprA", "aprB": "aprB",
    "dsrA": "dsrA", "dsrB": "dsrB", "dsrC": "dsrC", "dsrD": "dsrD", "dsrM": "dsrM", "dsrK": "dsrK",
    "qmoA": "qmoA", "qmoB": "qmoB", "qmoC": "qmoC",
    # sulfur oxidation (SCycDB has no sdo family → sdo is a coverage gap)
    "soxA": "soxA", "soxB": "soxB", "soxC": "soxC", "soxD": "soxD",
    "soxX": "soxX", "soxY": "soxY", "soxZ": "soxZ",
    "sqr": "sqr", "fccA": "fccA", "fccB": "fccB",
    "soeA": "soeA", "soeB": "soeB", "soeC": "soeC",
    "sorA": "sorA", "tsdA": "tsdA", "sor": "sor",
    # assimilatory sulfate reduction (cysN_cysC = bifunctional fusion → scycle cysNC)
    "cysN": "cysN", "cysD": "cysD", "cysN_cysC": "cysNC", "cysC": "cysC", "cysH": "cysH",
    "cysJ": "cysJ", "cysI": "cysI", "sir": "sir", "cysK": "cysK", "cysM": "cysM",
    # thiosulfate / polysulfide (no tth family → tth is a coverage gap)
    "phsA": "phsA", "phsB": "phsB", "phsC": "phsC",
    "ttrA": "ttrA", "ttrB": "ttrB", "ttrC": "ttrC",
    "doxA": "doxA", "doxD": "doxD", "otr": "otr", "sseA": "sseA", "sreA": "sreA",
    # organic sulfur / DMSP
    "dmdA": "dmdA", "dddP": "dddP", "mddA": "mddA", "mtoX": "mtoX",
    # sulfonate / taurine
    "tauD": "tauD", "tauA": "tauA", "tauB": "tauB", "tauC": "tauC", "ssuD": "ssuD", "ssuE": "ssuE",
}


def load_scycdb(path: Path | None, genomes: set[str]) -> dict[tuple[str, str], bool]:
    """SCycDB (Yu et al. 2020) + DIAMOND. Reads the normalized TSV (sample, family,
    count); family→target via SCYC_MAP; present iff count>0. Returns {} if no path."""
    rows = _read_norm_tsv(path)
    if not rows:
        return {}
    hit: dict[str, set[str]] = {}
    for r in rows:
        tid = SCYC_MAP.get(r["family"])
        if tid and float(r.get("count", 1) or 0) > 0:
            hit.setdefault(r.get("genome", r.get("sample")), set()).add(tid)
    return {(g, tid): tid in hit.get(g, set())
            for g in genomes for tid in all_target_ids()}


# Registry: name → (loader, needs_path). The engine calls each; empty result = skipped.
LOADERS = {
    "scycle":   (load_scycle,   False),
    "kofam":    (load_kofam,    False),
    "metabolic": (load_metabolic, True),
    "dram":     (load_dram,     True),
    "scycdb":   (load_scycdb,   True),
}
