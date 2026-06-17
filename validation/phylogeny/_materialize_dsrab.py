#!/usr/bin/env python3
import json, os
from collections import Counter
OUTDIR = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(OUTDIR, "_dsrab_seqs.json")))

# scycle presence-gate refs (circularity overlap)
SCYCLE_PRESENCE = {"P45574","Q59109","O33998","S0FXQ9","Q9F4A3",  # dsrA gate
                   "P45575","Q59110","D3RSN2","S0G5M9","Q9F4A2"}  # dsrB gate

faa, tsv = [], ["accession\tgene\ttype\tsubtype\torganism\tsource\tnotes"]
kept = []
for acc, gene, typ, subtype, org, src, note, seq in rows:
    seq = seq.replace("*","").strip()
    if len(seq) < 100:
        continue
    overlap = " [scycle-presence-ref]" if acc in SCYCLE_PRESENCE else ""
    org_clean = org.split(" (strain")[0]
    faa.append(f">{acc}|{gene}|{typ}|{org_clean.replace(' ','_')}")
    for i in range(0, len(seq), 60):
        faa.append(seq[i:i+60])
    tsv.append(f"{acc}\t{gene}\t{typ}\t{subtype}\t{org_clean}\t{src}\t{note}{overlap}")
    kept.append((gene, typ, subtype, acc in SCYCLE_PRESENCE))

open(os.path.join(OUTDIR,"dsrab_refs.faa"),"w").write("\n".join(faa)+"\n")
open(os.path.join(OUTDIR,"dsrab_refs.tsv"),"w").write("\n".join(tsv)+"\n")
print("total:", len(kept))
print("by gene x type:", dict(Counter((g,t) for g,t,_,_ in kept)))
print("by subtype:", dict(Counter(s for _,_,s,_ in kept)))
print("scycle-overlap:", sum(1 for *_,o in kept if o))

if __name__ == "__main__":
    pass
