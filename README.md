# scycle-pipeline

Maps MAGs / isolate proteomes to their participation in the **sulfur cycle**.
Modeled on `Nitrogen_Cycle/ncycle-pipeline` (and its `Holomicrobiome-ewaste`
ancestor); detection is **KO-primary** (KOfam HMMs + adaptive per-KO thresholds)
with custom clade HMMs and DIAMOND-BLAST gating for the homology traps. Covers
**61 targets** (across 6 sulfur-cycle process modules), **11 obligatory complexes**,
and **9 process-completeness synergies**.

**Status: reference layer built, runs end-to-end.** Engine cloned + generalized from
the validated ncycle-pipeline; sulfur biology authored fresh (`config/targets.yaml`,
[`../Info-sulfur.md`](../Info-sulfur.md)); curated BLAST seeds + trained custom HMMs
for `soxB/C/D/X`, `sqr`, `sdo`, `tth`, `doxD`, `dsrA`, `dsrB` **lifted from the
validated ewaste-pipeline**. Smoke-validated on a 5-genome panel: *D. vulgaris* →
complete sulfate reduction (1.0, dsrAB **reductive**), *P. denitrificans* → complete
Sox oxidation (1.0), *A. ferrooxidans* → sqr/sdo/doxD/tth sulfur oxidation,
*E. coli* → complete assimilatory reduction, *S. pneumoniae* → negative. The accuracy
validation campaign (reference panel + ground truth + regression gate) is the deferred
follow-on phase — see [`../Info-sulfur.md`](../Info-sulfur.md) "Hardening backlog".

## Run

```bash
cd scycle-pipeline
python run.py --input <dir-of-.faa-or-.fna> --cores 8
# First run creates the standalone `scycle-pipeline` conda env from envs/scycle.yaml.
# To reuse an existing compatible env instead (e.g. the ewaste one):
#   SCYCLE_ENV=ewaste-pipeline python run.py --input <dir> --skip-db-setup
```

`run.py` auto-detects protein (`.faa`) vs nucleotide (`.fna`, → Prodigal) input,
builds the databases on first run, then dispatches Snakemake.

## How it works

1. **Gene calls** — proteomes used directly; nucleotide assemblies → Prodigal.
2. **HMM scan** — `hmmscan` against `resources/hmm/scycle_targets.hmm`, a
   concatenation of the **KOfam profile HMM for every KO in `config/targets.yaml`**
   + the **custom clade HMMs** (`soxB`, `soxD`, `soxX`, `tth`, `doxD` — lifted from
   the ewaste curation, calibrated trusted cutoffs) + a Pfam fallback for KO-less
   targets (`dsrC`/PF04358, `dsrD`/PF08679, `soxD`/PF00034). Thresholds (KOfam
   `ko_list`, custom TCs, per-KO overrides) live in `resources/hmm/tc_cutoffs.tsv`.
3. **BLAST gating** — `diamond blastp` against **curated, clade-spread UniProt seeds**
   (`resources/blast_db/`), per-target `blast_identity_min`, tagged `>{target_id}||{acc}`.
4. **Calls** — `apply_rules.py` integrates evidence per target. Signature
   precedence **custom HMM > KO > Pfam**; targets with
   `requires_blast_for_confirmation` (the homology traps: `dsrAB`, `soxC`, `soxD`,
   `sqr`, `sdo`, `phsA`, `ttrA`, `sreA`) need a clade-specific BLAST hit or are
   `disqualified`. The **dsrAB direction** (reductive sulfate reduction vs
   reverse/oxidative rDSR) is resolved by `resolve_dsr_direction` (dsrD/qmo presence,
   refined by dsrD↔dsrAB operon synteny on nucleotide/MAG input).
5. **Reports** — per-sample `calls/scycle_calls.tsv`, `complex_completeness.tsv`,
   `synergy_completeness.tsv`, `report/gap_analysis.txt`; cross-sample
   `multisample_matrix.tsv` + heatmap + per-pathway/complex/synergy/process figures.

## Pathways & targets

Six process modules / 61 markers (see [`../Info-sulfur.md`](../Info-sulfur.md) for the
full per-gene atlas, KO/Pfam anchors, and trap rationale):

1. **Dissimilatory sulfate reduction** (SO₄²⁻→H₂S) — sat, aprAB, dsrAB(CD), dsrMK, qmoABC
2. **Sulfur oxidation** (H₂S/S⁰/S₂O₃²⁻→SO₄²⁻) — Sox (soxABCDXYZ), sqr, fccAB, sdo, soeABC, sorA, tsdA, sor
3. **Assimilatory sulfate reduction** (→cysteine) — cysND(C), cysC, cysH, cysJI, sir, cysKM
4. **Thiosulfate / tetrathionate / polysulfide** — phsABC, ttrABC, doxAD, tth, otr, sseA, sreA
5. **Organic sulfur (DMSP/DMS)** — dmdA, dddP, mddA, mtoX
6. **Sulfonate / taurine** — tauD, tauABC, ssuD, ssuE

## Build / test

```bash
python workflow/scripts/build_hmm_db.py     # KOfam KO profiles + custom HMMs + Pfam fallbacks
python workflow/scripts/build_blast_db.py    # curated UniProt seeds → DIAMOND DBs
```

The KOfam cache (`resources/.cache/`) is symlinked to the ncycle-pipeline download
(~1.5 GB profiles.tar.gz, pinned release 2026-05-24) to avoid a re-fetch.
