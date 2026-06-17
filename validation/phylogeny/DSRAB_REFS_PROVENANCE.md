# DsrAB type-labeled reference set — provenance & curation rationale

**Purpose.** An ORTHOGONAL, phylogeny-anchored ground truth to call the *direction*
of a dissimilatory (bi)sulfite reductase (DsrAB) — **reductive** (sulfate/sulfite
reducers, SO3²⁻→H2S) vs **oxidative / reverse (rDSR)** (sulfur oxidizers,
H2S/S⁰→SO3²⁻) — from sequence ancestry ALONE. This owes nothing to scycle's
companion-gene logic (`apply_rules.resolve_dsr_direction`: dsrD/qmo/sox/sqr
proximity), so a best-hit / phylogenetic-placement call against this set is a
*non-circular* check on the pipeline's direction verdict.

**Files**
- `dsrab_refs.faa` — 46 curated reference proteins, headers `>{acc}|{gene}|{type}|{organism}`.
- `dsrab_refs.tsv` — `accession  gene  type  subtype  organism  source  notes`.
- `dsrab_refs.dmnd` — DIAMOND v2.2.0 database for the lightweight best-hit tier.
- Build scripts (kept for reproducibility): `_query_dsrab.py` (UniProt search by
  organism × subunit × length window, prefers reviewed), `_materialize_dsrab.py`.

## Composition (verified-fetched 2026-06-13)

| gene | reductive | oxidative | total |
|------|-----------|-----------|-------|
| dsrA | 16        | 8         | 24    |
| dsrB | 13        | 9         | 22    |
| **total** | **29** | **17** | **46** |

Subtypes: reductive_bacterial 24 · reductive_archaeal 5 · oxidative_bacterial 17.

Every sequence was fetched from the UniProt REST API
(`https://rest.uniprot.org/uniprotkb/{acc}.fasta` / `.json`), filtered to the
expected subunit length (dsrA 360–520 aa, dsrB 320–470 aa) and protein/gene-name
match, and (where available) restricted to Swiss-Prot *reviewed* entries — all 46
selected are reviewed. No accession or sequence was invented.

### Clade coverage (Müller et al. 2015 framework, see below)
- **Reductive bacterial-type:** *Nitratidesulfovibrio (Desulfovibrio) vulgaris*
  Hildenborough, *D. desulfuricans*, *Megalodesulfovibrio gigas*,
  *Desulfomicrobium baculatum*, *Desulfobacter vibrioformis*, *Desulfobulbus
  rhabdoformis*, *Desulfosudis (Desulfococcus) oleivorans*, *Desulforamulus
  (Desulfotomaculum) reducens*, *Desulfitobacterium hafniense*, *Desulfosporosinus
  nitroreducens*, *Thermodesulfovibrio yellowstonii* (Nitrospirota, deep-branching),
  *Thermodesulfobacterium commune*, *Desulfarculus baarsii*, *Desulfotignum
  phosphitoxidans* — spans Desulfobacterota (former Deltaproteobacteria),
  Firmicutes, Nitrospirota, Thermodesulfobacteriota.
- **Reductive archaeal-type:** *Archaeoglobus fulgidus*, *A. veneficus*,
  *A. sulfaticallidus* (Euryarchaeota).
- **Oxidative / reverse (rDSR) bacterial-type:** *Allochromatium vinosum* (purple
  sulfur bacterium, the biochemically characterized rDSR), *Chlorobaculum tepidum*,
  *C. thiosulfatiphilum*, *Chlorobium phaeobacteroides* (green sulfur bacteria),
  *Thioalkalivibrio nitratireducens*, *Thiocapsa marina/imhoffii*, *Magnetococcus
  marinus*, *Thiothrix lacustris*, *Halorhodospira neutriphila*.

## Type-clade definitions (the field standard)

**Reference framework — Müller AL, Kjeldsen KU, Rattei T, Pester M, Loy A (2015).
"Phylogenetic and environmental diversity of DsrAB-type dissimilatory (bi)sulfite
reductases." *ISME J* 9:1152–1165.** doi:10.1038/ismej.2014.208.

DsrA and DsrB arose from an ancient gene duplication and co-evolve; a concatenated
DsrAB tree resolves into:
1. **Reductive bacterial-type** — the bulk of sulfate/sulfite-reducing bacteria
   (Desulfobacterota, sulfate-reducing Firmicutes, Nitrospirota, etc.). Catalyzes
   SO3²⁻ → H2S in vivo.
2. **Reductive archaeal-type** — *Archaeoglobus* and relatives; sister-grouping
   reflects the archaeal split, still reductive direction.
3. **Oxidative bacterial-type (rDSR)** — sulfur-oxidizing bacteria (purple/green
   sulfur bacteria, many Chromatiales, some Beta/Alpha/Gammaproteobacteria) that run
   the enzyme in REVERSE (H2S/S⁰ → SO3²⁻). This clade is monophyletic and
   sequence-distinguishable from the reductive clades.
4. Additional **uncultured/environmental "family-level" lineages** (Müller's many
   reductive-type clades of unknown physiology) — NOT represented here by isolates;
   obtain the full DB (below) to place these.

Direction (reductive vs oxidative) maps onto the deep tree topology, which is why a
phylogenetic placement / best-hit against type-labeled references is a defensible,
companion-gene-independent direction call.

## Circularity note (scycle presence gate)

scycle's dsrA/dsrB PRESENCE BLAST gate (`config/targets.yaml`) uses UniProt refs
`P45574, Q59109, O33998, S0FXQ9, Q9F4A3` (dsrA) and `P45575, Q59110, D3RSN2,
S0G5M9, Q9F4A2` (dsrB). **9 of these are included here and are flagged
`[scycle-presence-ref]` in `dsrab_refs.tsv`** (D3RSN2 was not re-fetched). They are
correctly TYPE-labeled (note O33998/O33999 *Allochromatium vinosum* = **oxidative**;
the rest reductive). Critically, the direction reference is expanded FAR beyond the
presence set (46 vs ≤10), spanning every type clade, so a direction call from this
set is independent of the genes scycle uses merely to detect dsrAB *presence*. The
presence gate only asks "is this dsrAB at all?"; direction is decided by clade
membership across the full 46.

## KNOWN CAVEATS (what a reviewer will raise)

1. **HGT of dsrAB.** Müller 2015 explicitly documents lateral transfer of dsrAB,
   producing incongruence between the DsrAB tree and the organismal (16S/genome)
   tree. A sequence can therefore carry a clade label that mismatches the host's
   genome-inferred lineage. Direction (the reductive↔oxidative split) is more
   conserved than fine taxonomy, but HGT means clade placement ≠ taxonomy.
2. **Reductive-bacterial-type dsrAB in sulfur OXIDIZERS.** Some sulfur oxidizers
   (and the reverse) carry dsrAB that branches in the "wrong" physiological clade.
   The most cited cases: *Desulfurivibrio alkaliphilus* (sulfur disproportionation /
   oxidation while carrying reductive-type dsrAB) and certain SUP05/Thioglobus and
   "GSO" lineages. Direction inferred from clade can thus contradict the organism's
   actual in-situ flux — the single biggest defensibility caveat. The orthogonal
   set REDUCES (does not eliminate) circularity; it should be reported as
   "sequence-clade direction" not "physiological direction."
3. **Paraphyly / unresolved environmental clades.** The reductive-type radiation is
   broad and paraphyletic with respect to the oxidative clade in some treatments;
   many environmental family-level lineages have no cultured representative and
   uncertain direction. Isolate-only references (this set) cannot place those — use
   the full Müller DB + placement (below).
4. **Divergent / periplasmic and "non-canonical" Dsr.** Some Dsr-like and
   reductase-type sequences (e.g., certain Acidobacteria, candidate phyla) are
   divergent enough that best-hit identity drops and placement is low-confidence.
5. **dsrAB ≠ assimilatory sulfite reductase.** The assimilatory enzymes (asrAB,
   cysIJ/sir) are a different family; the length/name filter here excludes them, but
   a raw BLAST against a query proteome can pull asr/cys hits — gate on the dsrAB
   length window and clade bitscore, not bare KO/Pfam.
6. **DsrA vs DsrB.** They are paralogous subunits; mixing them in one tree without
   separating by subunit creates two clusters. For direction calls, place dsrA
   against dsrA refs and dsrB against dsrB refs separately (the headers/TSV carry
   the `gene` field for this).

## Obtaining the FULL Müller 2015 DsrAB database (rigorous EPA-ng / gappa tier)

The 46 isolates here anchor every type clade for a fast best-hit / sanity tier. For
the publication-grade phylogenetic-placement tier (EPA-ng + gappa onto a fixed
reference tree), use the **complete curated DsrAB reference package** from Müller
et al. 2015:

- **Primary archive (authors' supplement, Loy lab):**
  https://doi.org/10.1038/ismej.2014.208 — Supplementary Information of the ISME J
  paper contains the curated DsrAB alignment, the reference tree, and the
  family-level clade assignments (reductive bacterial-type, reductive
  archaeal-type, oxidative bacterial-type, and the numbered environmental lineages).
- **Maintained copy / updates:** the Loy lab (Division of Microbial Ecology,
  University of Vienna) distributes the DsrAB database and a curated ARB database;
  see https://www.microbial-ecology.net/ and the lab software/resources page. The
  package is also redistributed within **GraftM** dsrA/dsrB gpkgs and in the
  **DiSCo** (Dissimilatory Sulfur Cycle) classifier
  (https://github.com/ericHester/disco), which ships the Müller-derived reference
  set with reductive/oxidative labels and is the most turnkey way to consume it.
- **Recommended placement workflow:** (1) align query DsrA/DsrB to the Müller
  reference alignment with `hmmalign`/`mafft --add`; (2) `epa-ng --tree
  <ref.tree> --ref-msa <ref.aln> --query <query.aln>`; (3) `gappa examine assign`
  / `gappa examine graft` using the Müller clade taxonomy to read off the
  reductive/oxidative/archaeal label and a placement support value.

Cite Müller et al. 2015 (ISME J 9:1152) for the type-clade framework, plus
Loy/Pester for the curated DB and (if used) Hester et al. for DiSCo.

## Lightweight best-hit tier (this set)

```
diamond blastp --db dsrab_refs.dmnd --query <proteome.faa> \
  --outfmt 6 qseqid sseqid pident length bitscore stitle \
  --max-target-seqs 5 --evalue 1e-20 -k 5
```
Read the `type` token from the top hit's `sseqid` (`{acc}|{gene}|{type}|{org}`);
require a clade-consistent margin (top reductive vs top oxidative bitscore gap) and
ideally agreement across the top-N before calling direction. Treat ties / small
margins as ambiguous (the rDSR↔reductive split can be <60% identity in divergent
lineages — see caveat 2/4).
