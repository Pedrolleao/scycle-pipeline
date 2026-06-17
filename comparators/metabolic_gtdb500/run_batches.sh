#!/usr/bin/env bash
# Run METABOLIC-G.pl on the 119 enriched GTDB-500 genomes, batched sequentially.
# Batches of 12 stay under the exit-144 scale-failure threshold (40+ genomes).
# NEVER run two METABOLIC instances concurrently (shared cwd in the repo).
set -u
source "$(conda info --base)/etc/profile.d/conda.sh"; set +u
conda activate METABOLIC_v4.0

BASE=/home/dmin/Grants/Sulfur_Cycle/comparators
WORK=$BASE/metabolic_gtdb500
REPO=$BASE/METABOLIC
LOG=$WORK/run_batches.log
: > "$LOG"

echo "[$(date '+%F %T')] START — $(ls -d "$WORK"/batches/b?? | wc -l) batches" | tee -a "$LOG"
cd "$REPO" || { echo "no METABOLIC repo at $REPO" | tee -a "$LOG"; exit 2; }

ok=0; fail=0
for bd in "$WORK"/batches/b??; do
  b=$(basename "$bd")
  out="$WORK/out/$b"
  rm -rf "$out"
  rm -f "$bd"/total.faa            # defensive: never let a leftover total.faa poison a rerun
  ng=$(ls "$bd"/*.faa 2>/dev/null | wc -l)
  echo "[$(date '+%F %T')] $b START ($ng proteomes)" | tee -a "$LOG"
  perl METABOLIC-G.pl -in "$bd" -o "$out" -t 12 -kofam-db full \
       >>"$WORK/$b.metabolic.log" 2>&1
  rc=$?
  kdir="$out/KEGG_identifier_result"
  nres=$(ls "$kdir"/*.result.txt 2>/dev/null | wc -l)
  if [ "$nres" -gt 0 ]; then
    cp "$kdir"/*.result.txt "$WORK/kegg_all/" 2>/dev/null
    echo "[$(date '+%F %T')] $b DONE rc=$rc results=$nres (expected $ng)" | tee -a "$LOG"
    ok=$((ok+1))
  else
    echo "[$(date '+%F %T')] $b FAILED rc=$rc results=0 — see $b.metabolic.log" | tee -a "$LOG"
    fail=$((fail+1))
  fi
done

tot=$(ls "$WORK/kegg_all"/*.result.txt 2>/dev/null | wc -l)
echo "[$(date '+%F %T')] ALL DONE — batches ok=$ok fail=$fail ; total result.txt collected=$tot / 119" | tee -a "$LOG"
echo "[ALLDONE]" | tee -a "$LOG"
