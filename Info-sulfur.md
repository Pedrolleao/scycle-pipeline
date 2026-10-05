# Info-sulfur.md — the sulfur-cycle atlas behind `scycle-pipeline`

Narrative rationale for every target, complex, and synergy in
`config/targets.yaml`. Companion to the nitrogen-cycle
`Info-nitrogen.md`; the detection machinery is identical (KO-primary KOfam HMMs +
adaptive thresholds, custom clade HMMs for the hard families, DIAMOND-BLAST gating).
All KO numbers were verified against the cached KOfam `ko_list` (release 2026-05-24).

The sulfur cycle is a redox wheel spanning sulfate (SO₄²⁻, +6) → sulfite (SO₃²⁻, +4)
→ elemental sulfur / thiosulfate (S⁰ / S₂O₃²⁻, ~0) → sulfide (H₂S, −2), plus the
organic-sulfur shunt (DMSP/DMS, sulfonates). The pipeline resolves a genome's role
across **6 process modules / 61 marker genes**.

---

## 1. Dissimilatory sulfate reduction — SO₄²⁻ → H₂S (energy metabolism)

The defining metabolism of sulfate-reducing bacteria/archaea (SRB/SRA;
*Desulfovibrio*, *Archaeoglobus*). Sulfate is activated to APS, reduced to sulfite,
then to sulfide by the Dsr system.

| gene | KO | role |
|---|---|---|
| sat | K00958 | sulfate adenylyltransferase (SO₄²⁻ → APS) — **shared with assimilation** |
| aprA / aprB | K00394 / K00395 | adenylylsulfate (APS) reductase α / β (APS → sulfite) |
| dsrA / dsrB | K11180 / K11181 | dissimilatory sulfite reductase α / β (siroheme; SO₃²⁻ → H₂S) |
| dsrC | — (PF04358) | DsrC sulfur-carrier; **no KOfam KO** |
| dsrD | — (PF08679) | DsrD — **reductive-direction marker** (see trap below) |
| dsrM / dsrK | K27187 / K27188 | DsrMKJOP membrane complex (recycles DsrC trisulfide) |
| qmoA / qmoB / qmoC | K16885 / K16886 / K16887 | Qmo — couples APS reductase to menaquinone |

⚠️ **K11179 is NOT dsrC** — the KOfam definition is "tRNA 2-thiouridine synthesizing
protein E". dsrC and dsrD therefore have no KO and are anchored on their Pfams.

## 2. Sulfur oxidation — H₂S / S⁰ / S₂O₃²⁻ → SO₄²⁻

Sulfur-oxidizing bacteria/archaea (SOB; *Paracoccus*, *Acidithiobacillus*,
*Beggiatoa*, green/purple sulfur bacteria). Two backbones: the periplasmic **Sox**
thiosulfate-oxidizing system and the **reverse-Dsr** route.

| gene | KO | role |
|---|---|---|
| soxA / soxX | K17222 / K17223 | Sox core c-type cytochromes (custom HMM for soxX) |
| soxB | K17224 | sulfate thiohydrolase (custom HMM; clade-calibrated) |
| soxC / soxD | K17225 / — | sulfane dehydrogenase (soxD custom HMM, no KO) |
| soxY / soxZ | K17226 / K17227 | sulfur carrier / chelator (obligate Sox-core subunits) |
| sqr | K17218 | sulfide:quinone oxidoreductase (H₂S → S⁰) — ⚠️ BLAST-gated |
| fccA / fccB | K17230 / K17229 | flavocytochrome c sulfide dehydrogenase (alt. to sqr) |
| sdo | K17725 | sulfur (persulfide) dioxygenase (S⁰ → sulfite) — ⚠️ BLAST-gated |
| soeA / soeB / soeC | K21307 / K21308 / K21309 | membrane sulfite:quinone dehydrogenase |
| sorA | K05301 | periplasmic sulfite:cytochrome dehydrogenase |
| tsdA | K19713 | thiosulfate dehydrogenase (S₂O₃²⁻ → tetrathionate) |
| sor | K16952 | sulfur oxygenase/reductase (thermoacidophile S⁰ disproportionation) |

**Branched vs complete Sox:** `soxB` + Sox-core but **no soxCD** is the S⁰-storing
branched Sox of green/purple sulfur bacteria; full `soxABCDXYZ` is the
Kelly–Friedrich complete oxidation with no S⁰ intermediate.

## 3. Assimilatory sulfate reduction — SO₄²⁻ → H₂S → cysteine (biosynthesis)

Near-universal route to biomass sulfur (the building of cysteine). Energetically
distinct from dissimilatory reduction — a separate enzyme set.

| gene | KO | role |
|---|---|---|
| cysN / cysD | K00956 / K00957 | assimilatory ATP sulfurylase (cysNC K00955 = fused) |
| cysC | K00860 | adenylylsulfate kinase (APS → PAPS) |
| cysH | K00390 | (P)APS reductase (PAPS → sulfite) |
| cysJ / cysI | K00380 / K00381 | NADPH sulfite reductase flavoprotein / hemoprotein |
| sir | K00392 | ferredoxin sulfite reductase (cyanobacteria/plants; alt. to cysJI) |
| cysK / cysM | K01738 / K12339 | cysteine synthase A / B (sulfide → cysteine) |

⚠️ cysI/sir share the siroheme NIR_SIR fold (PF01077+PF03460) with dsrAB and nitrite
reductases — the **KO tier separates them** (K00381/K00392 are assimilatory-specific).

## 4. Thiosulfate / tetrathionate / polysulfide

Respiratory and disproportionative middle-redox sulfur metabolism (*Salmonella*
tetrathionate respiration, acidophile tetrathionate cycling).

| gene | KO | role |
|---|---|---|
| phsA / phsB / phsC | K08352 / K08353 / K08354 | thiosulfate/polysulfide reductase (⚠️ phsA BLAST-gated) |
| ttrA / ttrB / ttrC | K08357 / K08358 / K08359 | tetrathionate reductase (⚠️ ttrA BLAST-gated) |
| doxA / doxD | K16936 / K16937 | thiosulfate:quinone oxidoreductase (doxD custom HMM) |
| tth | K27925 | tetrathionate hydrolase (custom HMM; acidophile periplasm) |
| otr | — | octaheme tetrathionate reductase — **Tier-3 BLAST-only**, no KO/Pfam |
| sseA | K01011 | rhodanese / thiosulfate sulfurtransferase (⚠️ broad superfamily) |
| sreA | K17219 | sulfur reductase Mo subunit (⚠️ BLAST-gated; S⁰ → H₂S respiration) |

## 5. Organic sulfur — DMSP / DMS / methanethiol

The marine climate-active branch: dimethylsulfoniopropionate (DMSP) is split between
demethylation (sulfur retained) and cleavage (DMS released to the atmosphere).

| gene | KO | role |
|---|---|---|
| dmdA | K17486 | DMSP demethylase — demethylation branch (suppresses DMS flux) |
| dddP | K28072 / K28073 | DMSP lyase (ddd family) — cleavage branch → DMS |
| mddA | K21310 | methanethiol S-methyltransferase (MeSH → DMS) |
| mtoX | K17285 | methanethiol oxidase (MeSH → formaldehyde + H₂S) |

## 6. Sulfonate / taurine (organosulfonate assimilation)

Sulfur scavenging from organosulfonates under sulfate starvation.

| gene | KO | role |
|---|---|---|
| tauD | K03119 | taurine dioxygenase (α-KG-dependent; releases sulfite) |
| tauA / tauB / tauC | K15551 / K10831 / K15552 | taurine ABC transporter |
| ssuD | K04091 | alkanesulfonate monooxygenase (⚠️ luciferase-like, broad) |
| ssuE | K00299 | FMN reductase partner of ssuD |

---

## The central trap — dsrAB reductive vs reverse/oxidative (rDSR)

`dsrA`/`dsrB` are the **same genes** in sulfate reducers (reductive Dsr, SO₃²⁻→H₂S)
and in sulfur oxidizers (reverse/oxidative rDSR, H₂S/S⁰→SO₃²⁻). No sequence method
separates the two directions — exactly the structure of the N-tool's `nxrA`/`narG`
trap. The **genomic companions** decide:

- **reductive** — `dsrD` present (and usually `qmoABC` / `dsrMK`); Sox machinery absent.
  `dsrD` sits in the canonical *dsrABD* operon, so on nucleotide/MAG input a `dsrD` ORF
  syntenic with `dsrA`/`dsrB` is a strong reductive call.
- **oxidative** — `dsrD` absent, Sox core / `sqr` present.

`apply_rules.resolve_dsr_direction` makes this call — presence-based on pre-called
proteomes, synteny-refined when Prodigal gene coordinates are available — and tags the
`dsrA`/`dsrB` evidence as `…|dsr_reductive` / `…|dsr_oxidative`. The
`complete_sulfate_reduction` vs `reverse_dsr_sulfur_oxidation` synergies encode the
same call at the process level. *Validated on the smoke panel: D. vulgaris (SRB) →
reductive (dsrD + qmoABC present, complete_sulfate_reduction = 1.0).*

### Other shared signatures
- **sat, aprAB** are shared across dissimilatory / assimilatory / oxidative sulfur
  metabolism — direction is inferred at the complex/synergy layer, not per gene.
- **soxC / sorA** share PF03473 (MOSC molybdopterin) — the KOs separate them.
- **sqr / fccB** share PF07992 (Pyr_redox_2 FAD) with `gor`/`lpdA` — BLAST-gated to
  canonical SQR clade seeds (lifted from the validated ewaste curation).
- **sdo** is in the broad metallo-β-lactamase fold (PF00753) — BLAST-gated.

---

## Obligatory complexes & process synergies

**Complexes (all subunits required):** aps_reductase (aprAB), dsr_reductase (dsrAB),
sox_core (soxABXYZ), sox_dehydrogenase (soxCD), qmo_complex (qmoABC), soe_sulfite_dh
(soeABC), assim_atp_sulfurylase (cysND), assim_sulfite_reductase (cysJI),
tetrathionate_reductase (ttrABC), thiosulfate_reductase (phsABC), taurine_transporter
(tauABC).

**Synergies (process completeness):** complete_sulfate_reduction (sat+aprAB+dsrAB),
reverse_dsr_sulfur_oxidation (dsrAB+aprA+sat, dsrD-absent), complete_sox_oxidation
(soxABXYZ), sulfide_oxidation (sqr | fccAB), assimilatory_complete
(cysN+cysH+cysI+cysK), thiosulfate_disproportionation (phsA+soxB),
dmsp_demethylation (dmdA), dmsp_cleavage (dddP), organosulfonate_scavenging
(tauD+ssuD).

---

## Provenance & reuse

The curated `blast_refs_uniprot` for `soxB/C/D/X`, `sqr`, `sdo`, `tth`, `doxD`,
`dsrA`, `dsrB` and the trained custom HMMs in `targets/{soxB,soxD,soxX,doxD,tth}/`
are **lifted verbatim from the validated `Holomicrobiome-ewaste/ewaste-pipeline`**,
which hardened these exact families (with documented Pfam corrections — e.g. dsrA
PF02661→PF01077, sqr PF02754→PF07992, sdo PF03867→PF00753). Empty `blast_refs_uniprot`
with a `ref_query:` hint mark targets whose seeds are to be curated in the hardening
phase — **never guessed** (the project's hardest historical bug class was fabricated
UniProt accessions).

## Hardening backlog (deferred validation phase)

The reference layer is built and runs end-to-end; the following are P2/P3-style tuning
items for the validation campaign (assemble a ~30-genome SRB/SOB/assimilator/decoy
panel + ground truth + a `make regression` floor, mirroring `ncycle-pipeline/ROADMAP.md`):
- **`sqr` / `soxD` over-disqualification** — the BLAST gates reject genuine `sqr` in
  *Desulfovibrio*/*Paracoccus* and `soxD` in some SOB (identity floor too tight for
  divergent clades). Curate per-clade seeds / relax `blast_identity_min` per target.
- **Seed the empty-`[]` targets** (phs/ttr/cys/tau/ddd/qmo BLAST refs) so the
  `narrow-no-IPR` cross-DB noise resolves to real calls.
- **Train custom HMMs** for `otr` (Tier-3 octaheme, currently BLAST-only) and the
  multiheme `tsdA`/`fccA` if BLAST proves insufficient.
- **`dddP` granularity** — KOfam folds the ddd lyase families under K28072/K28073;
  per-lyase resolution (dddL/P/Q/W/D) would need custom HMMs.
