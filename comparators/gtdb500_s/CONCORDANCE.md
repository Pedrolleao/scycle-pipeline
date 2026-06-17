# scycle-pipeline — GTDB-scale Cross-Tool CONCORDANCE

**Genomes:** 499 GTDB-representative (119 S-cycle-clade-enriched + 380 cross-phylum backbone) · **Tools:** scycle, raw KofamScan, METABOLIC, SCycDB · **No ground truth** — agreement/divergence only. (METABOLIC included via the 380-reused + 119-batched GTDB-500 run (see WORKPLAN); DRAM excluded at scale per prereg.)

> Concordance, not accuracy: arbitrary GTDB genomes have no curated truth. High agreement on unambiguous markers + systematic divergence at the homology/direction traps is the expected signature; the trap loci are where scycle's gating and dsrAB-direction resolution act.


## 1. Pairwise agreement

Over targets both tools can represent. **Raw** = identical present/absent call (inflated by shared-absent cells). **Jaccard** = positive-call agreement |both present| / |either present| (the honest concordance metric). Reported for ALL targets and the homology-trap subset.

| tool pair | raw (ALL) | raw (trap) | Jaccard (ALL) | Jaccard (trap) |
|---|---|---|---|---|
| scycle vs kofam | 0.692 | 0.665 | 0.225 | 0.140 |
| scycle vs metabolic | 0.977 | 0.924 | 0.784 | 0.331 |
| scycle vs scycdb | 0.658 | 0.777 | 0.215 | 0.203 |
| kofam vs metabolic | 0.697 | 0.737 | 0.238 | 0.327 |
| kofam vs scycdb | 0.642 | 0.680 | 0.401 | 0.351 |
| metabolic vs scycdb | 0.653 | 0.770 | 0.199 | 0.169 |

## 2. dsrAB reductive↔oxidative DIRECTION trap (headline)

dsrA/dsrB are the SAME genes (KO K11180/K11181) in dissimilatory sulfate reducers (REDUCTIVE Dsr, SO₃²⁻→H₂S) and in reverse-Dsr sulfur oxidizers (OXIDATIVE rDSR, H₂S/S⁰→SO₃²⁻). Gene presence alone cannot separate the two metabolic directions — the genomic companions do (dsrD ⇒ reductive; Sox/sqr with dsrD absent ⇒ oxidative). scycle resolves direction (`apply_rules.resolve_dsr_direction`); the direction-blind comparators (raw KofamScan, METABOLIC, SCycDB) report dsrA/dsrB presence with **no direction**.

- scycle calls **dsrAB** in **89** genomes, resolving direction as **48 reductive** (sulfate reducers), **40 oxidative** (reverse-Dsr sulfur oxidizers), **1 ambiguous**.
- Of the **40** oxidative (reverse-Dsr) genomes, **40** have ≥1 direction-blind comparator (raw KofamScan, METABOLIC, SCycDB) reporting dsrA/dsrB present — i.e. they would be **mis-attributed as sulfate reducers** by a tool that reads dsrAB presence as dissimilatory sulfate reduction. This is the dsrAB-direction analogue of ncycle's nxrA→narG mis-routing, now observed at GTDB scale.


## 3. Targets where scycle is DISTINCTIVE (differs from every comparator)

| target | n genomes scycle-distinct | trap? |
|---|---|---|
| sdo | 250 | trap |
| phsA | 46 | trap |
| fccA | 21 | trap |
| dsrC | 20 | trap |
| dsrD | 7 | trap |
| sat | 4 |  |
| ttrA | 3 | trap |

## 4. Per-target present-call counts (by tool; '-' = tool cannot represent)

Sorted by divergence (max−min present-count among representing tools).

| target | scycle | kofam | metabolic | scycdb | divergence |
|---|---|---|---|---|---|
| tauB | 18 | 492 | 15 | 469 | 477 |
| cysJ | 8 | 19 | 8 | 483 | 475 |
| ssuD | 6 | 75 | 6 | 461 | 455 |
| cysN | 36 | 483 | 36 | 445 | 447 |
| cysC | 53 | 258 | 53 | 498 | 445 |
| cysNC | 16 | 216 | 16 | 460 | 444 |
| sdo | 4 | 444 | 254 | - | 440 |
| cysI | 17 | 238 | 17 | 415 | 398 |
| sir | 5 | 95 | 5 | 395 | 390 |
| sqr | 113 | 461 | 72 | 129 | 389 |
| dddP | 11 | 0 | 2 | 383 | 383 |
| cysM | 111 | 348 | 111 | 476 | 365 |
| tauD | 2 | 3 | 2 | 357 | 355 |
| fccB | 41 | 352 | 0 | 153 | 352 |
| tauA | 17 | 233 | 17 | 367 | 350 |
| tsdA | 26 | 176 | 26 | 375 | 349 |
| aprA | 90 | 419 | 71 | 259 | 348 |
| cysD | 80 | 136 | 79 | 423 | 344 |
| phsA | 4 | 346 | 57 | 137 | 342 |
| ttrA | 2 | 245 | 5 | 344 | 342 |
| aprB | 98 | 411 | 76 | 93 | 335 |
| soxB | 68 | 266 | 30 | 338 | 308 |
| phsB | 0 | 306 | 0 | 34 | 306 |
| soeA | 49 | 345 | 43 | 80 | 302 |
| ttrB | 63 | 352 | 63 | 303 | 289 |
| ssuE | 8 | 173 | 8 | 294 | 286 |
| sat | 198 | 235 | 202 | 481 | 283 |
| dmdA | 15 | 50 | 15 | 294 | 279 |
| tauC | 19 | 276 | 19 | 91 | 257 |
| cysK | 242 | 320 | 241 | 490 | 249 |
| dsrA | 85 | 183 | 85 | 332 | 247 |
| soeB | 39 | 234 | 39 | 20 | 214 |
| ttrC | 4 | 199 | 4 | 99 | 195 |
| soxX | 111 | 203 | 107 | 298 | 191 |
| mddA | 40 | 223 | 40 | 90 | 183 |
| cysH | 134 | 193 | 134 | 314 | 180 |
| soxC | 29 | 48 | 22 | 195 | 173 |
| doxD | 34 | 160 | 28 | 11 | 149 |
| qmoA | 58 | 183 | 58 | 102 | 125 |
| sseA | 127 | 246 | 122 | 150 | 124 |
| dsrK | 87 | 203 | 87 | 110 | 116 |
| dsrB | 87 | 200 | 87 | 176 | 113 |
| soxA | 57 | 86 | 58 | 168 | 111 |
| qmoC | 60 | 166 | 60 | 64 | 106 |
| fccA | 17 | 121 | 43 | 53 | 104 |
| soxY | 69 | 74 | 59 | 158 | 99 |
| dsrC | 80 | - | 0 | 79 | 80 |
| tth | 5 | 85 | 5 | - | 80 |
| phsC | 0 | 15 | 0 | 76 | 76 |
| otr | 25 | - | 0 | 73 | 73 |
| soxD | 21 | - | 0 | 72 | 72 |
| qmoB | 67 | 124 | 67 | 90 | 57 |
| sorA | 17 | 31 | 17 | 69 | 52 |
| dsrM | 93 | 142 | 93 | 108 | 49 |
| soxZ | 72 | 75 | 72 | 120 | 48 |
| soeC | 39 | 83 | 39 | 40 | 44 |
| dsrD | 26 | - | 0 | 38 | 38 |
| mtoX | 22 | 24 | 22 | 3 | 21 |
| doxA | 0 | 10 | 0 | 0 | 10 |
| sor | 7 | 0 | 7 | 8 | 8 |
| sreA | 0 | 0 | 0 | 6 | 6 |
