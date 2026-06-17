#!/usr/bin/env python3
"""Build per-subunit combined ref+query FASTAs with SAFE leaf names + a label map.
Leaf names: R### (reference), Q### (query). Map JSON records type/genome/gene/etc."""
import csv, json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
GT = HERE.parents[1] / "comparators" / "gtdb500_s"
REF_FAA = HERE / "dsrab_refs.faa"
REF_TSV = HERE / "dsrab_refs.tsv"
QRY_FAA = GT / "_dsr_query.faa"
SEL = {r["ncbi_acc"].replace(".", "_"): r for r in csv.DictReader(open(GT / "selection.tsv"), delimiter="\t")}
# scycle direction per genome (from the best-hit scorer's TSV)
SCY = {r["genome"]: r["scycle_dir"] for r in csv.DictReader(open(GT / "DIR_ACCURACY.tsv"), delimiter="\t")}


def read_faa(p):
    seqs, name, buf = {}, None, []
    for line in open(p):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


ref_seqs = read_faa(REF_FAA)
ref_meta = {r["accession"]: r for r in csv.DictReader(open(REF_TSV), delimiter="\t")}
qry_seqs = read_faa(QRY_FAA)

label = {}
for gene in ("dsrA", "dsrB"):
    out = HERE / "trees" / f"{gene}_in.faa"
    out.parent.mkdir(exist_ok=True)
    ri = qi = 0
    with open(out, "w") as fh:
        for name, seq in ref_seqs.items():
            acc, g, typ, org = name.split("|")
            if g != gene:
                continue
            lid = f"R{ri:03d}"; ri += 1
            label[lid] = {"kind": "ref", "type": typ, "gene": gene, "acc": acc, "org": org}
            fh.write(f">{lid}\n{seq}\n")
        for name, seq in qry_seqs.items():
            genome, g, pid = name.split("__")
            if g != gene:
                continue
            lid = f"Q{qi:03d}_{gene}"; qi += 1
            cl = SEL.get(genome, {}).get("clade", "")
            label[lid] = {"kind": "query", "gene": gene, "genome": genome,
                          "clade": cl, "characterized": (not cl.startswith("p__")),
                          "scycle_dir": SCY.get(genome, "?")}
            fh.write(f">{lid}\n{seq}\n")
    print(f"{gene}: {ri} refs + {qi} queries -> {out}", file=sys.stderr)

json.dump(label, open(HERE / "trees" / "label_map.json", "w"), indent=0)
print("wrote label_map.json", file=sys.stderr)
