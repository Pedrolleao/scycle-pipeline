#!/usr/bin/env bash
# Merge the 380 reused (from N's METABOLIC run) + 119 freshly-run enriched genomes
# into the final gtdb500_s/metabolic.tsv (columns: genome, ko). Run AFTER run_batches.sh.
set -eu
BASE=/home/dmin/Grants
WORK=$BASE/Sulfur_Cycle/comparators/metabolic_gtdb500
OUT=$BASE/Sulfur_Cycle/scycle-pipeline/comparators/gtdb500_s/metabolic.tsv
NORM=$BASE/Nitrogen_Cycle/ncycle-pipeline/comparators/build_metabolic_tsv.py

# 1. normalize the 119 freshly-run genomes (kegg_all → genome,ko)
python3 "$NORM" --kegg-dir "$WORK/kegg_all" --out "$WORK/fresh_rows_h.tsv"
tail -n +2 "$WORK/fresh_rows_h.tsv" > "$WORK/fresh_rows.tsv"   # strip header

# 2. concatenate: single header + reuse (380) + fresh (119)
{ printf 'genome\tko\n'; cat "$WORK/reuse_rows.tsv" "$WORK/fresh_rows.tsv"; } > "$OUT"

# 3. validate
echo "=== merged metabolic.tsv ==="
echo "rows (excl header): $(($(wc -l < "$OUT") - 1))"
echo "distinct genomes  : $(tail -n +2 "$OUT" | cut -f1 | sort -u | wc -l)  (expect 499)"
echo "reuse genomes      : $(cut -f1 "$WORK/reuse_rows.tsv" | sort -u | wc -l)"
echo "fresh genomes      : $(cut -f1 "$WORK/fresh_rows.tsv" | sort -u | wc -l)"
# coverage vs the S run
comm -23 "$WORK/s_genomes.txt" <(tail -n +2 "$OUT" | cut -f1 | sort -u) > "$WORK/_uncovered.txt"
echo "S genomes with NO metabolic data: $(wc -l < "$WORK/_uncovered.txt")  (expect 0)"
[ -s "$WORK/_uncovered.txt" ] && { echo "UNCOVERED:"; cat "$WORK/_uncovered.txt"; }
echo "dsrA(K11180)=$(grep -c K11180 "$OUT")  dsrB(K11181)=$(grep -c K11181 "$OUT")  (full 499-genome panel)"
