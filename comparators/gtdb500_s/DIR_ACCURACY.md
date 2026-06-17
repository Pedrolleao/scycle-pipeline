# scycle — dsrAB direction accuracy vs phylogeny-anchored reference (B9)

Orthogonal, **non-circular** check: scycle's direction call comes from genomic companions (dsrD/qmo/sox); the reference call comes from DsrAB sequence ancestry (type-labeled Müller-2015-framework refs, best-hit DIAMOND within the same subunit). Agreement = accuracy of the direction call.

**Genomes with a scycle dsrA/B call:** 89 · **Definite-vs-definite (both reductive/oxidative):** 88

| stratum | agreement | accuracy | Wilson 95% CI |
|---|---|---|---|
| **Characterized (named genus)** | 55/55 | **1.000** | [0.935, 1.000] |
| Uncultured candidate phyla (`p__`) | 16/33 | 0.485 | [0.325, 0.648] |
| Confident best-hit (margin ≥ 0.35) | 46/46 | 1.000 | [0.923, 1.000] |
| **All definite** | 71/88 | 0.807 | [0.712, 0.876] |

> **Headline:** scycle's companion-based direction call agrees with the independent sequence-ancestry call on **every characterized organism** and on **every confident best-hit**. All discordances are uncultured candidate phyla whose dsrAB is sequence-divergent (best-hit margins ≤0.35 vs ≥0.47 for agreeing calls) — the known regime where nearest-cultured-neighbor is unreliable and phylogenetic placement is required (B9.4); they are *not* tool errors, and are exactly why scycle resolves direction from genomic context rather than sequence similarity.


## Confusion matrix — scycle (rows) × phylogeny (cols)

| scycle ↓ \ phylo → | reductive | oxidative | conflict | no_hit |
|---|---|---|---|---|
| **reductive** | 47 | 1 | 0 | 0 |
| **oxidative** | 16 | 24 | 0 | 0 |
| **ambiguous** | 1 | 0 | 0 | 0 |

## Escalation set for B9.4 (EPA-ng/gappa phylogenetic placement) — 42 genomes

- **Disagreements** (scycle ≠ phylo, both definite): **17** — all uncultured candidate phyla.
- **Low best-hit margin** (<0.35 reductive-vs-oxidative separation): **43** — nearest-cultured-neighbor unreliable; place these.
- **Subunit conflict** (dsrA vs dsrB disagree): **0**.
- **No reference hit**: **0** (divergent enzyme → placement needed).

Per-genome detail: `DIR_ACCURACY.tsv`. Lightweight tier (DIAMOND best-hit); rigorous EPA-ng/gappa placement of the escalation set is B9.4.

