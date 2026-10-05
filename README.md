# scycle-pipeline

Maps MAGs / isolate proteomes to their participation in the **sulfur cycle**.
Modeled on `Nitrogen_Cycle/ncycle-pipeline` (and its `Holomicrobiome-ewaste`
ancestor); detection is **KO-primary** (KOfam HMMs + adaptive per-KO thresholds)
with custom clade HMMs and DIAMOND-BLAST gating for the homology traps. Covers
**61 targets** (across 6 sulfur-cycle process modules), **11 obligatory complexes**,
and **9 process-completeness synergies**.

**Status: reference layer built, runs end-to-end.** Engine cloned + generalized from
the validated ncycle-pipeline; sulfur biology authored fresh (`config/targets.yaml`,
[`Info-sulfur.md`](Info-sulfur.md)); curated BLAST seeds + trained custom HMMs
for `soxB/C/D/X`, `sqr`, `sdo`, `tth`, `doxD`, `dsrA`, `dsrB` **lifted from the
validated ewaste-pipeline**. Smoke-validated on a 5-genome panel: *D. vulgaris* →
complete sulfate reduction (1.0, dsrAB **reductive**), *P. denitrificans* → complete
Sox oxidation (1.0), *A. ferrooxidans* → sqr/sdo/doxD/tth sulfur oxidation,
*E. coli* → complete assimilatory reduction, *S. pneumoniae* → negative. The accuracy
validation campaign (reference panel + ground truth + regression gate) is the deferred
follow-on phase — see [`Info-sulfur.md`](Info-sulfur.md) "Hardening backlog".

## Install

Developed and validated on Linux (x86-64). Needs `git`, conda or mamba and `curl`;
about 3 GB of disk for the environment. Network access is needed for the install and
for the five genomes of the smoke test (NCBI); the databases are built from files in
the repository and the pipeline itself runs offline.

```bash
git clone https://github.com/Pedrolleao/scycle-pipeline.git
cd scycle-pipeline
conda env create -f envs/scycle.yaml       # env `cycle-pipeline`; or: mamba env create …
conda activate cycle-pipeline
make test_protein     # fetch 5 genomes, build the databases, run end to end
make regression       # run the 43-genome reference panel (sp1_panel/) and check the accuracy floors
```

`make regression` ending in `OK: all … checks passed` means the install reproduces the
validated calls. The exact environment of the validation (linux-64) is `envs/ncycle.lock.yml` of the
nitrogen sister repository, for when the open version ranges of `envs/scycle.yaml`
resolve to something that behaves differently.

## Run

```bash
python run.py --input <dir-of-.faa-or-.fna> --cores 8
# Outside the conda env, run.py re-runs itself inside `cycle-pipeline` (and creates it
# from envs/scycle.yaml if it does not exist). To use another env with the same
# dependencies:
#   SCYCLE_ENV=<env-name> python run.py --input <dir>
```

`run.py` auto-detects protein (`.faa`) vs nucleotide (`.fna`, → Prodigal) input,
builds the databases on first run, then dispatches Snakemake. It asks whether
nucleotide input is isolate genomes or metagenome assemblies unless
`--prodigal-mode single|meta` is given, and writes the samples it found into the
`samples:` block of `config/config.yaml` — so `git status` shows that file as modified
after a run.

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
   `synergy_completeness.tsv`, `report/gap_analysis.txt`, `report/scycle_map.*`;
   cross-sample `multisample_matrix.tsv`, figures, and an interactive `report.html`
   (see **Outputs** below).

## Outputs

All paths are under `paths.results_dir` (`results/` by default). Figures are written
as SVG (vector) and PNG (300 DPI).

**Per sample — `<sample>/`**

| file | what it is |
|---|---|
| `calls/scycle_calls.tsv` | one row per target: status, evidence source, protein, HMM hit + E-value, BLAST reference + identity. For nucleotide input the columns `contig`, `start`, `end`, `strand` locate the called gene (empty for a pre-called proteome). One protein is reported per target; `other_copies` lists any further proteins that reach the same status (paralog copies) |
| `calls/complex_completeness.tsv`, `calls/synergy_completeness.tsv` | completeness of the 11 complexes and 9 process modules |
| `calls/scycle_loci.tsv` | *(nucleotide input only)* called genes grouped into loci: genes on one contig with at most 5 other genes between them. `copy` says whether a gene is the one reported in `scycle_calls.tsv` or an additional copy |
| `report/gap_analysis.txt` | plain-text summary |
| `report/scycle_map.svg/.png` | the genome's calls drawn on the sulfur cycle: each reaction arrow is solid (a complete route found), dashed (partial) or grey (absent), with the genes behind it. The sat / aprAB / dsrAB arrows point SO₄²⁻ → H₂S unless the dsrAB direction call is `oxidative` (reverse Dsr), in which case they are drawn the other way |
| `report/loci.svg/.png` | *(nucleotide input only)* gene-arrow maps of every locus, grouped by pathway, to a common bp scale. Genes outside the target set are blank; a bar marks a contig end (where an operon may run off the assembly) |

**Across samples**

| file | what it is |
|---|---|
| `multisample_matrix.tsv` | genomes × (targets, complexes, modules) |
| `multisample_heatmap.svg/.png` | overview dot grid; genomes ordered by gene-content similarity |
| `figures/pathway_<pathway>.svg/.png` | one dot grid per pathway |
| `figures/complexes.svg/.png`, `figures/synergies.svg/.png` | complex / process-module completeness |
| `figures/scycle_maps.svg/.png` | every genome's S-cycle map side by side (up to 48 genomes) |
| `report.html` | self-contained interactive report (no network needed): the gene grid and the complex / module grid with hover evidence, row search / ordering, and a per-genome panel with the S-cycle map, locus maps and the full calls table. Light and dark themes. Every figure in it (gene grid, complex / module grid, cycle map, each locus map) has a **Save PNG (300 dpi)** button: it downloads that figure as currently shown — row filter and order, hidden pathways, selected genome, light or dark theme — with its title and legend, rendered at 300 dpi (a grid too large for a browser canvas is saved at the highest resolution that fits, and says so). The page follows the group's *Simple Terminal* design system (`design/Simple`): JetBrains Mono, hairline `[ bracketed ]` frames, its dark palette or its Light variant according to the system theme, with a LIGHT / DARK selector in the top-right corner to pin either. The font is inlined from `workflow/scripts/fonts/` (SIL OFL 1.1, licence alongside), so the report looks the same offline and the PNG export uses it too; pathway colours stay the validated palette of the static figures. |

**Reading the glyphs** (same in every figure and in `report.html`): solid disc =
confirmed; half-filled = domain-only (HMM signature, no BLAST support); ring with a
cross = disqualified (failed the homology-trap gate); faint ring = absent. In the
complex / module grids: solid = complete, ring with `n/N` = partial, faint ring with
a cross = ruled out. Colour always means pathway; the palette lives in
`workflow/scripts/_domain.py`.

**dsrAB direction and the two dsr modules.** `complete_sulfate_reduction` and
`reverse_dsr_sulfur_oxidation` need the same genes, so they are told apart by the
dsrAB direction call (`reductive` / `oxidative`, tagged on the dsrA / dsrB
`evidence_source`). When a direction is called, the module that needs the opposite
direction is scored `absent` and listed under `forbids_violated` in
`synergy_completeness.tsv` (shown as "ruled out" in the figures); with no call, or an
ambiguous one, both modules are scored on gene presence alone.

The figure and report scripts are shared, byte-identical, with the nitrogen sister
pipeline; only `workflow/scripts/_domain.py` (palette, labels, file names) and
`workflow/scripts/_cycle_model.py` (the cycle diagram) are sulfur-specific.

## Pathways & targets

Six process modules / 61 markers (see [`Info-sulfur.md`](Info-sulfur.md) for the
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

**The database inputs are pinned in the repository.** The 58 KOfam profiles and
their thresholds are built from `resources/kofam_pinned/` (release of 2026-05-24, the
one the tool was validated on), the Pfam fallback profiles from `resources/pfam_pinned/`
and the BLAST seeds from `resources/seeds_pinned/`. genome.jp serves KOfam from a
rolling URL and keeps no old releases — the release of 2026-09-29 has a different
threshold for 49 of the 58 KOs — and UniProt entries are revised and deleted, so
databases built from fresh downloads are not the validated ones. `build_hmm_db.py
--upstream` and `build_blast_db.py --upstream` download the current data anyway;
re-validate (`make regression`) before trusting the result.

**What a clone does not contain.** Pipeline results and the third-party tools and
databases of the comparator benchmark (METABOLIC, DRAM, SCycDB, GTDB proteomes). The scripts under `comparators/` and
the study configs `config/config_*.yaml` are the record of how the validation was run:
they carry paths of the machine it ran on and need editing to be re-run elsewhere.
Their outputs — the tables under `validation/` and `comparators/` — are committed.
